#include "LocalizationWorker.hpp"

#include <algorithm>
#include <cmath>

namespace ProbeViewer::Localization
{
namespace
{
using Realtime::V3;

V3 sub (const V3& a, const V3& b) { return { a[0] - b[0], a[1] - b[1], a[2] - b[2] }; }
double dot (const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
V3 normalized (const V3& a)
{
    const double n = std::sqrt (dot (a, a));
    return n > 0.0 ? V3 { a[0] / n, a[1] / n, a[2] / n } : a;
}

bool finite (const V3& p)
{
    return std::isfinite (p[0]) && std::isfinite (p[1]) && std::isfinite (p[2]);
}

// Bundle directory: the configured one, else $PROBE_VIEWER_LOCALIZER_BUNDLE, else the newest
// ~/ephys-atlas/results/realtime/bundles/*_lf_ap_no_detection.
juce::File findBundleDirectory (const juce::String& configured)
{
    if (configured.isNotEmpty())
        return juce::File (configured);
    const auto env = juce::SystemStats::getEnvironmentVariable ("PROBE_VIEWER_LOCALIZER_BUNDLE", {});
    if (env.isNotEmpty())
        return juce::File (env);
    const auto root = juce::File::getSpecialLocation (juce::File::userHomeDirectory)
                          .getChildFile ("ephys-atlas/results/realtime/bundles");
    juce::File newest;
    for (const auto& dir : root.findChildFiles (juce::File::findDirectories, false, "*_lf_ap_no_detection"))
        if (newest == juce::File() || dir.getLastModificationTime() > newest.getLastModificationTime())
            newest = dir;
    return newest;
}

// A plane through the probe that contains the probe axis and one atlas axis (0: ML, a tilted
// coronal slice; 1: AP, a tilted sagittal slice), sampled from the bundle's region labels, as
// iblatlas BrainAtlas.tilted_slice.
AtlasSlice tiltedSlice (const Realtime::LocalizerBundle& bundle, const V3& topUm, const V3& tipUm, int axis, bool beryl)
{
    AtlasSlice slice;
    if (! finite (topUm) || ! finite (tipUm))
        return slice;
    const V3 up = normalized (sub (topUm, tipUm));
    const V3 u = axis == 0 ? V3 { 1.0, 0.0, 0.0 } : V3 { 0.0, 1.0, 0.0 };
    const double along = dot (up, u);
    V3 v = normalized (V3 { up[0] - along * u[0], up[1] - along * u[1], up[2] - along * u[2] });
    if (v[2] < 0.0)
        v = { -v[0], -v[1], -v[2] };
    const V3 origin { 0.5 * (topUm[0] + tipUm[0]), 0.5 * (topUm[1] + tipUm[1]), 0.5 * (topUm[2] + tipUm[2]) };
    const double halfLength = 0.5 * std::sqrt (dot (sub (topUm, tipUm), sub (topUm, tipUm)));
    const float step = 50.0f;
    slice.voxelSizeUm = step;
    slice.minUUm = -3000.0f;
    slice.maxUUm = 3000.0f;
    slice.minVUm = static_cast<float> (-halfLength - 1500.0);
    slice.maxVUm = static_cast<float> (halfLength + 1500.0);
    const V3 normal { u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0] };
    for (int k = 0; k < 3; ++k)
    {
        slice.originXyzUm[static_cast<size_t> (k)] = static_cast<float> (origin[static_cast<size_t> (k)]);
        slice.basisUXyz[static_cast<size_t> (k)] = static_cast<float> (u[static_cast<size_t> (k)]);
        slice.basisVXyz[static_cast<size_t> (k)] = static_cast<float> (v[static_cast<size_t> (k)]);
        slice.normalXyz[static_cast<size_t> (k)] = static_cast<float> (normal[static_cast<size_t> (k)]);
    }
    for (float vv = slice.minVUm; vv <= slice.maxVUm; vv += step)
        for (float uu = slice.minUUm; uu <= slice.maxUUm; uu += step)
        {
            // The bundle's labels cover the left hemisphere (the localizer works there); the atlas
            // is symmetric, so both hemispheres are drawn from it.
            const V3 pM { -std::abs (origin[0] + uu * u[0] + vv * v[0]) * 1e-6, (origin[1] + uu * u[1] + vv * v[1]) * 1e-6,
                          (origin[2] + uu * u[2] + vv * v[2]) * 1e-6 };
            const int row = bundle.regionAt (pM, beryl);
            SliceVoxel voxel;
            voxel.uUm = uu;
            voxel.vUm = vv;
            if (row >= 0 && bundle.regions[static_cast<size_t> (row)].id != 0)
            {
                const auto& r = bundle.regions[static_cast<size_t> (row)];
                voxel.regionId = static_cast<int> (r.id);
                voxel.colour = juce::Colour (r.rgb[0], r.rgb[1], r.rgb[2]);
            }
            slice.voxels.push_back (voxel);
        }
    slice.valid = true;
    return slice;
}

CcfPoint ccfPoint (const V3& xyzUm)
{
    CcfPoint p;
    p.mlUm = static_cast<float> (xyzUm[0]);
    p.apUm = static_cast<float> (xyzUm[1]);
    p.dvUm = static_cast<float> (xyzUm[2]);
    return p;
}
} // namespace

LocalizationWorker::LocalizationWorker() : thread ([this] { loop(); }) {}

LocalizationWorker::~LocalizationWorker()
{
    {
        std::lock_guard<std::mutex> lock (mutex);
        quitting = true;
    }
    cv.notify_one();
    if (thread.joinable())
        thread.join();
}

bool LocalizationWorker::canAcceptJob() const noexcept
{
    std::lock_guard<std::mutex> lock (mutex);
    return ! busy && ! hasJob;
}

bool LocalizationWorker::submitWindow (LocalizationInput input, LocalizationConfig config)
{
    std::lock_guard<std::mutex> lock (mutex);
    if (busy || hasJob)
        return false;
    pendingInput = std::move (input);
    pendingConfig = std::move (config);
    hasJob = true;
    status = "Queued";
    cv.notify_one();
    return true;
}

LocalizationResultPtr LocalizationWorker::getLatestResult() const
{
    std::lock_guard<std::mutex> lock (mutex);
    return latest;
}

juce::String LocalizationWorker::getStatus() const
{
    std::lock_guard<std::mutex> lock (mutex);
    return status;
}

float LocalizationWorker::getProgress() const
{
    return progress.load();
}

int LocalizationWorker::getSnippetCount() const
{
    return snippetCount.load();
}

void LocalizationWorker::reset()
{
    std::lock_guard<std::mutex> lock (mutex);
    latest.reset();
    hasJob = false;
    resetPending = true;
    snippetCount.store (0);
    status = busy ? "Running" : "Idle";
}

bool LocalizationWorker::prepareRealtime (const LocalizationConfig& config, int nChannels, juce::String& error)
{
    const auto dir = findBundleDirectory (config.realtime.bundleDirectory);
    if (! dir.isDirectory())
    {
        error = "no localizer bundle: export one with ibleatools scripts/export_localizer_bundle.py";
        return false;
    }
    if (bundle == nullptr || dir.getFullPathName() != bundleDirectory)
    {
        auto b = std::make_unique<Realtime::LocalizerBundle>();
        std::string err;
        if (! b->load (dir.getFullPathName().toStdString(), err))
        {
            error = "cannot load the localizer bundle " + dir.getFullPathName() + ": " + juce::String (err);
            return false;
        }
        realtime.reset();
        bundle = std::move (b);
        realtime = std::make_unique<Realtime::RealtimeEngine> (*bundle);
        bundleDirectory = dir.getFullPathName();
        planFile.clear();
        snippetCount.store (0);
    }
    if (config.initialTrajectoryFile != planFile || static_cast<int> (plan.size()) != nChannels)
    {
        if (config.initialTrajectoryFile.isEmpty())
        {
            error = "load the planned trajectory (Plan...: ephysatlas plan JSON, Pinpoint export or xyz trajectory)";
            return false;
        }
        std::string err, source;
        std::vector<Realtime::V3> p;
        if (! Realtime::loadPlanFile (config.initialTrajectoryFile.toStdString(), nChannels, p, source, err))
        {
            error = juce::String (err);
            return false;
        }
        plan = std::move (p);
        planFile = config.initialTrajectoryFile;
        planSource = juce::String (source);
        realtime->setPlan (plan);
    }
    return true;
}

LocalizationResult LocalizationWorker::runRealtime (LocalizationInput& input, const LocalizationConfig& config)
{
    LocalizationResult result;
    result.modelName = "ephys-only (LF + AP, real time)";
    result.probe = input.probe;
    result.apStreamKey = input.ap.streamKey;
    result.lfStreamKey = input.lf.streamKey;
    result.streamKey = input.ap.streamKey;
    result.createdTimeSeconds = juce::Time::getMillisecondCounterHiRes() * 1e-3;

    juce::String error;
    if (! prepareRealtime (config, input.ap.numChannels, error))
    {
        result.status = error;
        return result;
    }
    const auto electrodes = input.ap.electrodeIndices;
    progress.store (0.1f);
    std::string err;
    if (! realtime->addSnippet (input, err))
    {
        result.status = "snippet rejected: " + juce::String (err);
        return result;
    }
    snippetCount.store (realtime->snippetCount());
    progress.store (0.8f);
    const auto update = realtime->update();
    progress.store (1.0f);
    const auto& out = update.result;
    result.snippetCount = update.nSnippets;
    if (! out.ok)
    {
        result.status = "localization failed: " + juce::String (out.error);
        return result;
    }

    const bool beryl = config.realtime.berylRegions;
    const int channels = static_cast<int> (out.channelXyzM.size());
    std::vector<V3> xyzUm (static_cast<size_t> (channels));
    for (int r = 0; r < channels; ++r)
        for (int k = 0; k < 3; ++k)
            xyzUm[static_cast<size_t> (r)][static_cast<size_t> (k)] = out.channelXyzM[static_cast<size_t> (r)][static_cast<size_t> (k)] * 1e6;
    for (int r = 0; r < channels; ++r)
    {
        // Row r (top of the probe first) is source channel C - 1 - r.
        const int source = channels - 1 - r;
        ChannelLocalization ch;
        ch.channelIndex = r;
        ch.electrodeIndex = source < static_cast<int> (electrodes.size()) ? electrodes[static_cast<size_t> (source)] : source;
        ch.depthUm = static_cast<float> (20.0 * (source / 2) + 20.0);
        const int row = (beryl ? out.berylRow : out.cosmosRow)[static_cast<size_t> (r)];
        if (row >= 0)
        {
            const auto& region = bundle->regions[static_cast<size_t> (row)];
            ch.regionId = static_cast<int> (region.id);
            ch.acronym = juce::String (region.acronym);
            ch.colour = juce::Colour (region.rgb[0], region.rgb[1], region.rgb[2]);
        }
        ch.confidence = out.valid[static_cast<size_t> (r)] ? 1.0f : 0.0f;
        ch.ccf = ccfPoint (xyzUm[static_cast<size_t> (r)]);
        result.channels.push_back (ch);
    }
    for (const auto& p : out.traceXyzM)
        result.predictedTrajectoryXyzUm.push_back ({ static_cast<float> (p[0] * 1e6), static_cast<float> (p[1] * 1e6), static_cast<float> (p[2] * 1e6) });
    for (size_t r = 0; r < plan.size(); r += 8)
        result.initialTrajectoryXyzUm.push_back ({ static_cast<float> (plan[r][0] * 1e6), static_cast<float> (plan[r][1] * 1e6), static_cast<float> (plan[r][2] * 1e6) });
    if (! plan.empty())
        result.initialTrajectoryXyzUm.push_back ({ static_cast<float> (plan.back()[0] * 1e6), static_cast<float> (plan.back()[1] * 1e6), static_cast<float> (plan.back()[2] * 1e6) });
    if (channels > 0)
    {
        result.entry = ccfPoint (xyzUm.front());
        result.tip = ccfPoint (xyzUm.back());
        result.atlasSlice = tiltedSlice (*bundle, xyzUm.front(), xyzUm.back(), 0, beryl);
        result.secondarySlice = tiltedSlice (*bundle, xyzUm.front(), xyzUm.back(), 1, beryl);
    }
    result.objectiveScore = out.score;
    result.localizationFeatureCount = bundle->F;
    result.initializationSource = juce::File (planFile).getFileName() + " (" + planSource + ") | " + juce::String (update.nSnippets)
                                + " snippets | score " + juce::String (out.score, 3);
    result.initializationFilePath = planFile;
    result.valid = true;
    result.status = "Snippet " + juce::String (update.nSnippets) + " localized (score " + juce::String (out.score, 3) + ")";
    return result;
}

void LocalizationWorker::loop()
{
    for (;;)
    {
        LocalizationInput input;
        LocalizationConfig config;
        {
            std::unique_lock<std::mutex> lock (mutex);
            cv.wait (lock, [this] { return quitting || hasJob; });
            if (quitting)
                return;
            input = std::move (pendingInput);
            config = pendingConfig;
            hasJob = false;
            busy = true;
            status = "Running";
            if (resetPending)
            {
                if (realtime != nullptr)
                    realtime->reset();
                resetPending = false;
            }
        }
        progress.store (0.0f);
        auto result = std::make_shared<LocalizationResult> (runRealtime (input, config));
        {
            std::lock_guard<std::mutex> lock (mutex);
            busy = false;
            status = result->status;
            // A failed update keeps the last good result on screen.
            if (result->valid || latest == nullptr)
                latest = std::move (result);
        }
    }
}
}
