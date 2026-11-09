#pragma once
#include <JuceHeader.h>

namespace ProbeViewer::Localization
{
// Real-time ephys-only localization (Realtime::RealtimeEngine): the localizer of ibleatools
// (ephysatlas.alignment.ephys_only, LF + AP features without spike detection) running on a
// bundle exported by ephysatlas.realtime.bundle.
struct RealtimeLocalizationConfig
{
    // Empty: $PROBE_VIEWER_LOCALIZER_BUNDLE, else the newest
    // ~/ephys-atlas/results/realtime/bundles/*_lf_ap_no_detection.
    juce::String bundleDirectory;

    // One AP + LF snippet per window, as ephysatlas.realtime (5 s); snippets are taken back to
    // back whenever the worker is free, and every update uses all snippets so far.
    double snippetSeconds = 5.0;
    int maxSnippets = 100000;

    // Region level shown on the channels and the slices: Beryl, else Cosmos.
    bool berylRegions = true;
};

struct LocalizationConfig
{
    // Feature rows run from the top of the probe (row r is channel C - 1 - r), as the localizer
    // and the IBL feature tables expect.
    bool reverseChannelOrder = true;

    // Open Ephys continuous data are in µV.
    double inputScaleToVolts = 1.0e-6;

    // Used when an LF stream reports no sample rate (File Reader can).
    double lfFallbackSampleRateHz = 2500.0;

    // cor_ratio on every AP sample (1) or a deterministic subsample.
    int corRatioSampleStride = 1;

    // The planned trajectory: an ephysatlas-plan/1 JSON file, a Pinpoint v5 experiment, a
    // Pinpoint text export or an xyz trajectory (see Realtime/PlanFile.hpp), and the probe to
    // read from a Pinpoint v5 experiment with several probes.
    juce::String initialTrajectoryFile;
    juce::String initialTrajectoryProbe;

    RealtimeLocalizationConfig realtime;
};

} // namespace ProbeViewer::Localization
