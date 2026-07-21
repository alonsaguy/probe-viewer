#include "LocalizationEngine.hpp"
#include <random>

namespace ProbeViewer::Localization
{
LocalizationResult LocalizationEngine::run (const LocalizationInput& input, const LocalizationConfig& config) const
{
    LocalizationResult result;
    result.modelName = config.modelName; result.streamKey = input.streamKey;
    result.createdTimeSeconds = juce::Time::getMillisecondCounterHiRes() / 1000.0;
    if (input.numChannels <= 0) { result.status = "Empty localization input"; return result; }

    const auto features = extractor.extract (input, config.featureSource);
    std::mt19937 rng (config.randomSeed ^ static_cast<uint32_t> (juce::Time::getMillisecondCounter()));
    result.entry = config.prior.enabled ? CcfPoint { config.prior.entryApMlUm.x, config.prior.entryApMlUm.y, config.prior.entryDvUm }
                                        : atlas.sampleRandomPoint (rng);
    result.tip = result.entry;
    result.tip.dvUm = juce::jlimit (0.0f, 8000.0f, result.entry.dvUm + 3840.0f);

    result.channels.reserve (static_cast<size_t> (input.numChannels));
    for (int ch = 0; ch < input.numChannels; ++ch)
    {
        const float t = input.numChannels > 1 ? static_cast<float> (ch) / (input.numChannels - 1) : 0.0f;
        CcfPoint p { result.entry.apUm + t * (result.tip.apUm - result.entry.apUm),
                     result.entry.mlUm + t * (result.tip.mlUm - result.entry.mlUm),
                     result.entry.dvUm + t * (result.tip.dvUm - result.entry.dvUm) };
        const auto region = atlas.regionAt (p);
        float confidence = 0.5f;
        if (features.rows == input.numChannels && features.cols > 0)
            confidence = juce::jlimit (0.05f, 0.95f, 0.5f + 0.1f * std::tanh (features.values[static_cast<size_t> (ch * features.cols)]));
        result.channels.push_back ({ ch, t * 3840.0f, region.id, region.acronym, region.colour, confidence, p });
    }
    result.valid = true; result.status = "Random CCF demo localization complete";
    return result;
}
}
