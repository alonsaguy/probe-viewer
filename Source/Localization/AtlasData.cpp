#include "AtlasData.hpp"

namespace ProbeViewer::Localization
{
CcfPoint AtlasData::sampleRandomPoint (std::mt19937& rng) const
{
    std::uniform_real_distribution<float> ap (0.0f, 13200.0f);
    std::uniform_real_distribution<float> ml (0.0f, 11400.0f);
    std::uniform_real_distribution<float> dv (0.0f, 8000.0f);
    return { ap (rng), ml (rng), dv (rng) };
}

bool AtlasData::contains (const CcfPoint& p) const
{
    return p.apUm >= 0 && p.apUm <= 13200 && p.mlUm >= 0 && p.mlUm <= 11400 && p.dvUm >= 0 && p.dvUm <= 8000;
}

AtlasData::Region AtlasData::regionAt (const CcfPoint& p) const
{
    static const Region regions[] = {
        { 1, "CTX", juce::Colour::fromRGB (85, 170, 255) },
        { 2, "TH",  juce::Colour::fromRGB (255, 120, 170) },
        { 3, "STR", juce::Colour::fromRGB (120, 220, 120) },
        { 4, "MB",  juce::Colour::fromRGB (255, 190, 80) },
        { 5, "HY",  juce::Colour::fromRGB (200, 120, 240) },
        { 6, "CB",  juce::Colour::fromRGB (120, 210, 210) }
    };
    const int idx = juce::jlimit (0, 5, static_cast<int> ((p.dvUm / 8000.0f) * 6.0f));
    return regions[idx];
}
}
