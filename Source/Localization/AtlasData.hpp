#pragma once
#include "LocalizationResult.hpp"
#include <random>

namespace ProbeViewer::Localization
{
class AtlasData
{
public:
    struct Region { int id; juce::String acronym; juce::Colour colour; };
    CcfPoint sampleRandomPoint (std::mt19937& rng) const;
    Region regionAt (const CcfPoint& p) const;
    bool contains (const CcfPoint& p) const;
};
}
