#pragma once
#include "LocalizationTypes.hpp"
#include <JuceHeader.h>
#include <array>
#include <memory>
#include <vector>

namespace ProbeViewer::Localization
{
struct ChannelLocalization
{
    int channelIndex = -1;      // row order: 0 is the top of the probe
    int electrodeIndex = -1;    // Open Ephys neuropixels.electrode_index when available
    float depthUm = 0.0f;       // distance from the tip
    int regionId = 0;           // Allen id
    juce::String acronym { "OUT" };
    juce::Colour colour { juce::Colours::darkgrey };
    float confidence = 0.0f;    // 1 for channels with a signal, 0 otherwise
    CcfPoint ccf;
};

struct SliceVoxel
{
    float uUm = 0.0f;
    float vUm = 0.0f;
    int regionId = 0;
    juce::Colour colour { juce::Colours::black };
};

// A plane of the atlas: voxel colours at (u, v), and the basis that maps xyz (µm) onto it.
struct AtlasSlice
{
    bool valid = false;
    std::array<float, 3> originXyzUm {};
    std::array<float, 3> basisUXyz {};
    std::array<float, 3> basisVXyz {};
    std::array<float, 3> normalXyz {};
    float voxelSizeUm = 200.0f;
    float minUUm = 0.0f;
    float maxUUm = 0.0f;
    float minVUm = 0.0f;
    float maxVUm = 0.0f;
    std::vector<SliceVoxel> voxels;
};

struct LocalizationResult
{
    bool valid = false;
    juce::String status;
    juce::String modelName;
    ProbeIdentity probe;

    juce::String apStreamKey;
    juce::String lfStreamKey;
    juce::String streamKey;

    juce::String initializationSource; // what the plan came from, shown in the slice title
    juce::String initializationFilePath;
    double createdTimeSeconds = 0.0;
    double objectiveScore = 0.0;
    int localizationFeatureCount = 0;
    int snippetCount = 0;              // snippets used by this update

    CcfPoint entry;
    CcfPoint tip;
    std::vector<ChannelLocalization> channels;

    std::vector<std::array<float, 3>> initialTrajectoryXyzUm;   // the plan (µm)
    std::vector<std::array<float, 3>> predictedTrajectoryXyzUm; // the estimate (µm)

    // Tilted coronal slice, and the tilted sagittal one drawn beside it.
    AtlasSlice atlasSlice;
    AtlasSlice secondarySlice;
};

using LocalizationResultPtr = std::shared_ptr<const LocalizationResult>;
}
