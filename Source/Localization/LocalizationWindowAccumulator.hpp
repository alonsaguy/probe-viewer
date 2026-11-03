#pragma once
#include "FeatureExtractor.hpp"

#include <atomic>
#include <mutex>
#include <limits>

namespace ProbeViewer::Localization
{
struct StreamSpec
{
    juce::String key;
    double sampleRate = 0.0;
    int channels = 0;
    std::vector<int> electrodeIndices;
};

class LocalizationWindowAccumulator
{
public:
    void arm (
        ProbeIdentity probe,
        StreamSpec ap,
        StreamSpec lf,
        double apWindowSeconds,
        double lfWindowSeconds,
        int temporalWindowCount = 3);

    // Called by ProbeViewerNode after the worker has finished processing the
    // previous window. This allocates the next raw AP/LF pair only then, which
    // prevents multiple ~GB AP windows from coexisting in RAM.
    bool continueBatchIfNeeded();

    void reset();

    bool isArmed() const noexcept
    {
        return armed.load();
    }

    bool isWindowReady() const noexcept
    {
        return ready.load();
    }

    bool hasPendingBatchWindow() const noexcept
    {
        return batchActive.load()
            && ! armed.load()
            && ! ready.load();
    }

    int currentWindowNumber() const noexcept
    {
        return currentWindow.load() + 1;
    }

    int totalWindowCount() const noexcept
    {
        return totalWindows.load();
    }

    void addBlock (
        const juce::String& streamKey,
        const juce::AudioBuffer<float>& source,
        const std::vector<int>& globalChannelIndices,
        int numSamples);

    LocalizationInput makeInput();

    float progress() const noexcept;

private:
    struct State
    {
        StreamSpec spec;
        double reportedSampleRate =
            std::numeric_limits<double>::quiet_NaN();
        bool sampleRateFallbackUsed = false;
        int targetSamples = 0;
        int writtenSamples = 0;
        std::vector<float> samples;

        bool required() const noexcept
        {
            return spec.key.isNotEmpty();
        }

        bool complete() const noexcept
        {
            return ! required()
                || writtenSamples >= targetSamples;
        }

        size_t byteCount() const noexcept
        {
            return samples.size() * sizeof (float);
        }
    };

    void prepareState (
        State&,
        const StreamSpec&,
        double seconds,
        double fallbackSampleRateHz = 0.0);

    void prepareCurrentWindow();

    void append (
        State&,
        const juce::AudioBuffer<float>&,
        const std::vector<int>&,
        int);

    static void releaseState (State&);

    std::atomic<bool> armed { false };
    std::atomic<bool> ready { false };
    std::atomic<bool> batchActive { false };
    std::atomic<int> currentWindow { 0 };
    std::atomic<int> totalWindows { 1 };

    ProbeIdentity probeIdentity;
    StreamSpec apTemplate;
    StreamSpec lfTemplate;
    double apSeconds = 30.0;
    double lfSeconds = 60.0;

    State apState;
    State lfState;

    mutable std::mutex mutex;
};
}
