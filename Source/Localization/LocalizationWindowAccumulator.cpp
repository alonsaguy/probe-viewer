#include "LocalizationWindowAccumulator.hpp"
#include <algorithm>

namespace ProbeViewer::Localization
{
void LocalizationWindowAccumulator::arm (juce::String key, double rate, int channels, double seconds)
{
    std::lock_guard<std::mutex> lock (mutex);
    streamKey = std::move (key);
    sampleRate = rate;
    numChannels = channels;
    targetSamples = juce::jmax (1, static_cast<int> (std::round (rate * seconds)));
    writtenSamples = 0;
    samples.assign (static_cast<size_t> (numChannels) * static_cast<size_t> (targetSamples), 0.0f);
    ready.store (false);
    armed.store (true);
}

void LocalizationWindowAccumulator::reset()
{
    std::lock_guard<std::mutex> lock (mutex);
    armed.store (false); ready.store (false); samples.clear(); writtenSamples = 0;
}

void LocalizationWindowAccumulator::addBlock (const juce::AudioBuffer<float>& source,
                                                const std::vector<int>& globalIndices,
                                                int numSamples)
{
    if (! armed.load() || ready.load()) return;
    std::lock_guard<std::mutex> lock (mutex);
    const int n = juce::jmin (numSamples, targetSamples - writtenSamples);
    if (n <= 0) return;
    for (int ch = 0; ch < numChannels && ch < static_cast<int> (globalIndices.size()); ++ch)
    {
        const int src = globalIndices[static_cast<size_t> (ch)];
        if (src < 0 || src >= source.getNumChannels()) continue;
        std::copy_n (source.getReadPointer (src), n,
                     samples.begin() + static_cast<ptrdiff_t> (ch * targetSamples + writtenSamples));
    }
    writtenSamples += n;
    if (writtenSamples >= targetSamples) ready.store (true);
}

LocalizationInput LocalizationWindowAccumulator::makeInput()
{
    std::lock_guard<std::mutex> lock (mutex);
    LocalizationInput out;
    if (! ready.load()) return out;
    out.streamKey = streamKey; out.sampleRate = sampleRate; out.numChannels = numChannels;
    out.numSamples = targetSamples; out.channelMajorSamples = std::move (samples);
    armed.store (false); ready.store (false); writtenSamples = 0;
    return out;
}

float LocalizationWindowAccumulator::progress() const noexcept
{
    std::lock_guard<std::mutex> lock (mutex);
    return targetSamples > 0 ? juce::jlimit (0.0f, 1.0f, static_cast<float> (writtenSamples) / targetSamples) : 0.0f;
}
}
