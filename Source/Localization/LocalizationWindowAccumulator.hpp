#pragma once
#include "FeatureExtractor.hpp"
#include <mutex>

namespace ProbeViewer::Localization
{
class LocalizationWindowAccumulator
{
public:
    void arm (juce::String streamKey, double sampleRate, int numChannels, double windowSeconds);
    void reset();
    bool isArmed() const noexcept { return armed.load(); }
    bool isWindowReady() const noexcept { return ready.load(); }
    void addBlock (const juce::AudioBuffer<float>& source, const std::vector<int>& globalChannelIndices, int numSamples);
    LocalizationInput makeInput();
    float progress() const noexcept;

private:
    std::atomic<bool> armed { false }, ready { false };
    juce::String streamKey;
    double sampleRate = 0.0;
    int numChannels = 0, targetSamples = 0, writtenSamples = 0;
    std::vector<float> samples;
    mutable std::mutex mutex;
};
}
