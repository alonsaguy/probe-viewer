#include "RealtimeEngine.hpp"

#include "IblLfFeatures.hpp"

#include <chrono>
#include <limits>

namespace ProbeViewer::Localization::Realtime
{
namespace
{
double secondsSince (std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
}
} // namespace

RealtimeEngine::RealtimeEngine (const LocalizerBundle& b) : bundle (b), localizer (b)
{
}

void RealtimeEngine::setPlan (std::vector<V3> plannedM)
{
    const std::lock_guard<std::mutex> lock (mutex);
    plan = std::move (plannedM);
}

bool RealtimeEngine::hasPlan() const
{
    const std::lock_guard<std::mutex> lock (mutex);
    return ! plan.empty();
}

bool RealtimeEngine::addSnippet (LocalizationInput& input, std::string& error, double* seconds)
{
    const auto t0 = std::chrono::steady_clock::now();
    const int channels = input.ap.numChannels;
    if (channels <= 0 || input.lf.numChannels != channels)
    {
        error = "the AP and LF snippets need the same, non-zero, channel count";
        return false;
    }
    const FeatureMatrix fm = extractor.extract (input, config);
    SnippetFeatureSet s;
    const int nf = bundle.F;
    s.values.assign (static_cast<size_t> (channels) * static_cast<size_t> (nf), std::numeric_limits<double>::quiet_NaN());
    for (int f = 0; f < nf; ++f)
    {
        int col = -1;
        for (int c = 0; c < fm.cols; ++c)
            if (fm.names[static_cast<size_t> (c)].toStdString() == bundle.featureNames[static_cast<size_t> (f)])
                col = c;
        if (col < 0)
        {
            error = "the feature extractor does not compute " + bundle.featureNames[static_cast<size_t> (f)];
            return false;
        }
        for (int r = 0; r < fm.rows; ++r)
        {
            const int channel = config.reverseChannelOrder ? channels - 1 - r : r;
            s.values[static_cast<size_t> (channel) * static_cast<size_t> (nf) + static_cast<size_t> (f)] = fm.at (r, col);
        }
    }
    s.labels = fm.channelLabels;
    {
        const std::lock_guard<std::mutex> lock (mutex);
        if (! snippets.empty() && nChannels != channels)
        {
            error = "the channel count changed between snippets";
            return false;
        }
        nChannels = channels;
        snippets.push_back (std::move (s));
    }
    if (seconds)
        *seconds = secondsSince (t0);
    return true;
}

EngineUpdate RealtimeEngine::update (long long seed) const
{
    EngineUpdate u;
    std::vector<SnippetFeatureSet> used;
    std::vector<V3> planned;
    int channels = 0;
    {
        const std::lock_guard<std::mutex> lock (mutex);
        used = snippets;
        planned = plan;
        channels = nChannels;
    }
    u.nSnippets = static_cast<int> (used.size());
    if (used.empty())
    {
        u.result.error = "no snippet yet";
        return u;
    }
    if (static_cast<int> (planned.size()) != channels)
    {
        u.result.error = "the plan has " + std::to_string (planned.size()) + " channels, the recording " + std::to_string (channels);
        return u;
    }
    const int nf = bundle.F;
    auto t0 = std::chrono::steady_clock::now();
    std::vector<char> isLogAmplitude;
    for (const auto& name : bundle.featureNames)
        isLogAmplitude.push_back (name.rfind ("rms_", 0) == 0);
    const auto geometry = Np1Geometry::make (channels);
    const auto agg = aggregateSnippetFeatures (used, channels, nf, isLogAmplitude, geometry.x, geometry.y, 0.1, &u.labels);
    // Row order for the localizer: row r is channel C - 1 - r.
    u.features.assign (agg.size(), 0.0);
    for (int r = 0; r < channels; ++r)
        for (int f = 0; f < nf; ++f)
            u.features[static_cast<size_t> (r) * static_cast<size_t> (nf) + static_cast<size_t> (f)] =
                agg[static_cast<size_t> (channels - 1 - r) * static_cast<size_t> (nf) + static_cast<size_t> (f)];
    u.aggregateSeconds = secondsSince (t0);
    t0 = std::chrono::steady_clock::now();
    u.result = localizer.localize (u.features, channels, planned, seed < 0 ? bundle.seed : seed);
    u.localizeSeconds = secondsSince (t0);
    return u;
}

int RealtimeEngine::snippetCount() const
{
    const std::lock_guard<std::mutex> lock (mutex);
    return static_cast<int> (snippets.size());
}

void RealtimeEngine::reset()
{
    const std::lock_guard<std::mutex> lock (mutex);
    snippets.clear();
    nChannels = 0;
}
} // namespace ProbeViewer::Localization::Realtime
