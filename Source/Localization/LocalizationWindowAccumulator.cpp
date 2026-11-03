#include "LocalizationWindowAccumulator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>

namespace ProbeViewer::Localization
{
namespace
{
double bytesToMiB (size_t bytes)
{
    return static_cast<double> (bytes)
        / (1024.0 * 1024.0);
}
}

void LocalizationWindowAccumulator::releaseState (
    State& state)
{
    std::vector<float>().swap (state.samples);
    state = {};
}

void LocalizationWindowAccumulator::prepareState (
    State& state,
    const StreamSpec& spec,
    double seconds,
    double fallbackSampleRateHz)
{
    releaseState (state);
    state.spec = spec;
    state.reportedSampleRate = spec.sampleRate;

    if (! state.required())
        return;

    if (state.spec.sampleRate <= 0.0
        && fallbackSampleRateHz > 0.0)
    {
        state.spec.sampleRate =
            fallbackSampleRateHz;

        state.sampleRateFallbackUsed = true;

        DBG (
            "Localization accumulator: stream "
            + state.spec.key
            + " reported sampleRate="
            + juce::String (
                state.reportedSampleRate,
                6)
            + " Hz; using "
            + juce::String (
                state.spec.sampleRate,
                3)
            + " Hz capture fallback.");
    }

    if (state.spec.channels <= 0
        || state.spec.sampleRate <= 0.0
        || seconds <= 0.0)
    {
        DBG (
            "Localization accumulator ignored invalid stream "
            + state.spec.key
            + " (reported fs="
            + juce::String (
                state.reportedSampleRate,
                6)
            + ", effective fs="
            + juce::String (
                state.spec.sampleRate,
                6)
            + ")");
        releaseState (state);
        return;
    }

    state.targetSamples =
        juce::jmax (
            1,
            static_cast<int> (
                std::round (
                    state.spec.sampleRate * seconds)));

    const size_t channelCount =
        static_cast<size_t> (state.spec.channels);
    const size_t sampleCount =
        static_cast<size_t> (state.targetSamples);

    if (sampleCount > 0
        && channelCount
            > std::numeric_limits<size_t>::max() / sampleCount)
    {
        DBG (
            "Localization accumulator allocation overflow for "
            + state.spec.key);
        releaseState (state);
        return;
    }

    const size_t totalValues =
        channelCount * sampleCount;

    try
    {
        state.samples.resize (totalValues);
    }
    catch (const std::bad_alloc&)
    {
        DBG (
            "Localization accumulator could not allocate "
            + state.spec.key + " window");
        releaseState (state);
        return;
    }

    DBG (
        "Localization temporal window "
        + juce::String (currentWindow.load() + 1)
        + "/"
        + juce::String (totalWindows.load())
        + " "
        + state.spec.key
        + ": "
        + juce::String (state.spec.channels)
        + " channels x "
        + juce::String (state.targetSamples)
        + " samples = "
        + juce::String (bytesToMiB (state.byteCount()), 1)
        + " MiB");
}

void LocalizationWindowAccumulator::prepareCurrentWindow()
{
    prepareState (
        apState,
        apTemplate,
        apSeconds);

    // NP1 LF is canonically 2500 Hz. This makes the requested zero-rate
    // fallback effective before allocation/capture; previously a zero-rate LF
    // stream was discarded before FeatureExtractor could apply its fallback.
    prepareState (
        lfState,
        lfTemplate,
        lfSeconds,
        2500.0);

    ready.store (false);
    armed.store (
        apState.required() || lfState.required());
}

void LocalizationWindowAccumulator::arm (
    ProbeIdentity probe,
    StreamSpec ap,
    StreamSpec lf,
    double apWindowSeconds,
    double lfWindowSeconds,
    int temporalWindowCount)
{
    std::lock_guard<std::mutex> lock (mutex);

    probeIdentity = std::move (probe);
    apTemplate = std::move (ap);
    lfTemplate = std::move (lf);
    apSeconds = apWindowSeconds;
    lfSeconds = lfWindowSeconds;

    currentWindow.store (0);
    totalWindows.store (
        juce::jmax (1, temporalWindowCount));
    batchActive.store (true);

    prepareCurrentWindow();
}

bool LocalizationWindowAccumulator::continueBatchIfNeeded()
{
    std::lock_guard<std::mutex> lock (mutex);

    if (! batchActive.load()
        || armed.load()
        || ready.load())
    {
        return false;
    }

    if (currentWindow.load() >= totalWindows.load())
    {
        batchActive.store (false);
        return false;
    }

    prepareCurrentWindow();

    if (! armed.load())
        batchActive.store (false);

    return armed.load();
}

void LocalizationWindowAccumulator::reset()
{
    std::lock_guard<std::mutex> lock (mutex);

    armed.store (false);
    ready.store (false);
    batchActive.store (false);
    currentWindow.store (0);
    totalWindows.store (1);

    releaseState (apState);
    releaseState (lfState);

    apTemplate = {};
    lfTemplate = {};
    probeIdentity = {};
}

void LocalizationWindowAccumulator::append (
    State& state,
    const juce::AudioBuffer<float>& source,
    const std::vector<int>& indices,
    int numSamples)
{
    if (state.complete())
        return;

    const int count =
        juce::jmin (
            numSamples,
            state.targetSamples - state.writtenSamples);

    if (count <= 0)
        return;

    for (int channel = 0;
         channel < state.spec.channels
             && channel < static_cast<int> (indices.size());
         ++channel)
    {
        const int sourceChannel =
            indices[static_cast<size_t> (channel)];

        if (sourceChannel < 0
            || sourceChannel >= source.getNumChannels())
        {
            continue;
        }

        std::copy_n (
            source.getReadPointer (sourceChannel),
            count,
            state.samples.begin()
                + static_cast<ptrdiff_t> (
                    channel * state.targetSamples
                    + state.writtenSamples));
    }

    state.writtenSamples += count;
}

void LocalizationWindowAccumulator::addBlock (
    const juce::String& key,
    const juce::AudioBuffer<float>& source,
    const std::vector<int>& indices,
    int numSamples)
{
    if (! armed.load() || ready.load())
        return;

    std::lock_guard<std::mutex> lock (mutex);

    if (key == apState.spec.key)
    {
        append (
            apState,
            source,
            indices,
            numSamples);
    }
    else if (key == lfState.spec.key)
    {
        append (
            lfState,
            source,
            indices,
            numSamples);
    }

    ready.store (
        apState.complete() && lfState.complete());
}

LocalizationInput
LocalizationWindowAccumulator::makeInput()
{
    std::lock_guard<std::mutex> lock (mutex);

    LocalizationInput output;

    if (! ready.load())
        return output;

    output.probe = probeIdentity;
    output.probeName = probeIdentity.displayName;
    output.temporalWindowIndex = currentWindow.load();
    output.temporalWindowCount = totalWindows.load();

    auto moveState =
        [] (State& state)
        {
            BandWindow band;
            band.streamKey = state.spec.key;
            band.sampleRate = state.spec.sampleRate;
            band.reportedSampleRate =
                state.reportedSampleRate;
            band.numChannels = state.spec.channels;
            band.numSamples = state.targetSamples;
            band.electrodeIndices = state.spec.electrodeIndices;
            band.channelMajorSamples = std::move (state.samples);
            return band;
        };

    output.ap = moveState (apState);
    output.lf = moveState (lfState);

    armed.store (false);
    ready.store (false);

    releaseState (apState);
    releaseState (lfState);

    const int next =
        currentWindow.fetch_add (1) + 1;

    if (next >= totalWindows.load())
        batchActive.store (false);

    return output;
}

float LocalizationWindowAccumulator::progress()
    const noexcept
{
    std::lock_guard<std::mutex> lock (mutex);

    const int total =
        juce::jmax (1, totalWindows.load());

    float within = 0.0f;
    int count = 0;

    for (const State* state :
         { &apState, &lfState })
    {
        if (! state->required())
            continue;

        within +=
            state->targetSamples > 0
                ? juce::jlimit (
                    0.0f,
                    1.0f,
                    static_cast<float> (state->writtenSamples)
                        / static_cast<float> (state->targetSamples))
                : 0.0f;
        ++count;
    }

    if (count > 0)
        within /= static_cast<float> (count);

    return juce::jlimit (
        0.0f,
        1.0f,
        (static_cast<float> (currentWindow.load()) + within)
            / static_cast<float> (total));
}
}
