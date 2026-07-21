#pragma once

#include <JuceHeader.h>

namespace ProbeViewer::Localization
{
enum class FeatureSource
{
    probeViewerFeatures = 0,
    rawWindow = 1
};

struct PlannedTrajectory
{
    bool enabled = false;
    juce::Point<float> entryApMlUm { 0.0f, 0.0f };
    float entryDvUm = 0.0f;
    float polarDeg = 0.0f;
    float azimuthDeg = 0.0f;
};

struct LocalizationConfig
{
    bool enabled = false;
    double windowSeconds = 100.0;
    double updateIntervalSeconds = 0.0; // 0 = one shot per button press
    FeatureSource featureSource = FeatureSource::probeViewerFeatures;
    juce::String modelName { "random-ccf-demo" };
    PlannedTrajectory prior;
    uint32_t randomSeed = 0x51A7u;
};
}
