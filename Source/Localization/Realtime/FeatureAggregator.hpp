#pragma once

// Snippet aggregation of ephysatlas.realtime.aggregate_snippets: per channel median over the
// snippets (ephysatlas.aggregation.aggregate_raw_features), mode of the channel labels, then
// ephysatlas.features.denoise_dataframe (labelled channels set to NaN, denoise_shank with
// fac 0.1) and ephysatlas.data.replace_nan (remaining NaN -> column median).

#include <vector>

namespace ProbeViewer::Localization::Realtime
{
struct SnippetFeatureSet
{
    std::vector<double> values; // [nc x nf], source channel order
    std::vector<int> labels;    // [nc] bad-channel labels; empty = all good
};

// denoise_shank: interpolate the feature onto the grid of unique (x, y) positions
// (scipy griddata linear, nearest outside the hull), TV-Chambolle denoise with
// weight = median(|image|) * fac, and read the result back at the channels. Returns the input
// unchanged when fewer than 5 values are finite.
std::vector<double> denoiseShank (const std::vector<double>& feature, const std::vector<double>& x,
                                  const std::vector<double>& y, double fac);

// isLogAmplitude[f]: feature f is 20 log10 of an amplitude (the rms features), which the
// release aggregates as the median amplitude before taking the log. Returns [nc x nf] in source
// channel order; labelsOut receives the aggregated labels.
std::vector<double> aggregateSnippetFeatures (const std::vector<SnippetFeatureSet>& snippets, int nc, int nf,
                                              const std::vector<char>& isLogAmplitude,
                                              const std::vector<double>& x, const std::vector<double>& y,
                                              double fac = 0.1, std::vector<int>* labelsOut = nullptr);
} // namespace ProbeViewer::Localization::Realtime
