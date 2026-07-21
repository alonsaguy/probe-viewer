#include "FeatureExtractor.hpp"
#include <cmath>

namespace ProbeViewer::Localization
{
FeatureMatrix FeatureExtractor::extract (const LocalizationInput& input, FeatureSource source) const
{
    if (source == FeatureSource::probeViewerFeatures
        && input.featuresPerChannel > 0
        && static_cast<int> (input.probeViewerFeatures.size()) == input.numChannels * input.featuresPerChannel)
        return { input.numChannels, input.featuresPerChannel, input.probeViewerFeatures };

    return computeRawRms (input);
}

FeatureMatrix FeatureExtractor::computeRawRms (const LocalizationInput& input) const
{
    FeatureMatrix out { input.numChannels, 1, std::vector<float> (static_cast<size_t> (input.numChannels), 0.0f) };
    for (int ch = 0; ch < input.numChannels; ++ch)
    {
        double sumSq = 0.0;
        const auto offset = static_cast<size_t> (ch) * static_cast<size_t> (input.numSamples);
        for (int i = 0; i < input.numSamples; ++i)
        {
            const double x = input.channelMajorSamples[offset + static_cast<size_t> (i)];
            sumSq += x * x;
        }
        out.values[static_cast<size_t> (ch)] = input.numSamples > 0 ? static_cast<float> (std::sqrt (sumSq / input.numSamples)) : 0.0f;
    }
    return out;
}
}
