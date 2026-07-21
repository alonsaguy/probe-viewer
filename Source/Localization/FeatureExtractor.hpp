#pragma once
#include "LocalizationConfig.hpp"
#include <vector>

namespace ProbeViewer::Localization
{
struct LocalizationInput
{
    juce::String streamKey;
    double sampleRate = 0.0;
    int numChannels = 0;
    int numSamples = 0;
    std::vector<float> channelMajorSamples;
    std::vector<float> probeViewerFeatures;
    int featuresPerChannel = 0;
};

struct FeatureMatrix
{
    int rows = 0;
    int cols = 0;
    std::vector<float> values;
};

class FeatureExtractor
{
public:
    FeatureMatrix extract (const LocalizationInput&, FeatureSource) const;
private:
    FeatureMatrix computeRawRms (const LocalizationInput&) const;
};
}
