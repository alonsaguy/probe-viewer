#pragma once
#include "AtlasData.hpp"
#include "FeatureExtractor.hpp"
#include "LocalizationConfig.hpp"
#include "LocalizationModel.hpp"

namespace ProbeViewer::Localization
{
class LocalizationEngine
{
public:
    LocalizationResult run (const LocalizationInput&, const LocalizationConfig&) const;
private:
    AtlasData atlas;
    FeatureExtractor extractor;
    LocalizationModel model;
};
}
