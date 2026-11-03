#pragma once

// Real-time ephys-only localization, as ephysatlas.realtime.online.OnlinePipeline: each 5 s
// AP + LF snippet gives one set of channel features (FeatureExtractor, which reproduces
// ephysatlas.realtime.snippet_features); an update aggregates every snippet so far
// (aggregate_snippets) and runs the localizer on the result.

#include "../FeatureExtractor.hpp"
#include "EphysOnlyLocalizer.hpp"
#include "FeatureAggregator.hpp"
#include "PlanFile.hpp"

#include <mutex>
#include <string>
#include <vector>

namespace ProbeViewer::Localization::Realtime
{
struct EngineUpdate
{
    int nSnippets = 0;
    std::vector<double> features; // [C x F] row order, the bundle's features
    std::vector<int> labels;      // [C] aggregated bad-channel labels, source channel order
    LocalizationOutput result;
    double aggregateSeconds = 0.0, localizeSeconds = 0.0;
};

class RealtimeEngine
{
public:
    // The bundle must outlive the engine.
    explicit RealtimeEngine (const LocalizerBundle& bundle);

    void setPlan (std::vector<V3> plannedM);
    bool hasPlan() const;

    // Computes and stores the features of one snippet (consumes the input buffers). Safe to
    // call from a worker thread while another thread runs update().
    bool addSnippet (LocalizationInput& input, std::string& error, double* seconds = nullptr);

    // Aggregates the stored snippets and localizes; seed < 0 uses the bundle's seed.
    EngineUpdate update (long long seed = -1) const;

    int snippetCount() const;
    void reset();

private:
    const LocalizerBundle& bundle;
    EphysOnlyLocalizer localizer;
    LocalizationConfig config;
    FeatureExtractor extractor;
    std::vector<V3> plan;
    std::vector<SnippetFeatureSet> snippets;
    int nChannels = 0;
    mutable std::mutex mutex;
};
} // namespace ProbeViewer::Localization::Realtime
