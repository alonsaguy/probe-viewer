#pragma once

#include <JuceHeader.h>
#include <memory>
#include <vector>

namespace ProbeViewer::Localization
{
struct CcfPoint
{
    float apUm = 0.0f;
    float mlUm = 0.0f;
    float dvUm = 0.0f;
};

struct ChannelLocalization
{
    int channelIndex = -1;
    float depthUm = 0.0f;
    int regionId = 0;
    juce::String acronym { "OUT" };
    juce::Colour colour { juce::Colours::darkgrey };
    float confidence = 0.0f;
    CcfPoint ccf;
};

struct LocalizationResult
{
    bool valid = false;
    juce::String status;
    juce::String modelName;
    juce::String streamKey;
    double createdTimeSeconds = 0.0;
    CcfPoint entry;
    CcfPoint tip;
    std::vector<ChannelLocalization> channels;
};

using LocalizationResultPtr = std::shared_ptr<const LocalizationResult>;
}
