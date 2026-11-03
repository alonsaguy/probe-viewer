#pragma once

// Port of specparam 2.0 SpectralModel.fit (fixed aperiodic mode, gaussian peaks, MAE error) as
// used by ephysatlas.features.get_psd_decay_features: peak_width_limits=[10, 15], max_n_peaks=4,
// min_peak_height=0.1, the other settings at their specparam defaults.
//
// The aperiodic fits are linear in their parameters, so they are solved exactly by least squares
// (scipy's unbounded curve_fit converges to the same solution). The bounded gaussian fit uses a
// projected Levenberg-Marquardt run to convergence in place of scipy's trust-region-reflective
// solver; both minimise the same problem from the same guess.

#include <array>
#include <vector>

namespace ProbeViewer::Localization::Realtime
{
struct SpecParamSettings
{
    double peakWidthLo = 10.0, peakWidthHi = 15.0; // two-sided bandwidth limits (Hz)
    int maxNPeaks = 4;
    double minPeakHeight = 0.1;
    double peakThreshold = 2.0;          // in units of the flattened spectrum's std
    double apPercentileThresh = 0.025;   // percentile (not fraction) for the robust aperiodic fit
    double cfBound = 1.5;
    double bwStdEdge = 1.0;
    double gaussOverlapThresh = 0.75;
};

struct SpecParamResult
{
    bool ok = false;
    double offset = 0.0, exponent = 0.0; // log10 power = offset - exponent * log10(f)
    double error = 0.0;                  // mean absolute error in log10 power
    double rSquared = 0.0;
    std::vector<std::array<double, 3>> gaussians; // [centre, height, std], sorted by centre
};

// freqs (Hz, evenly spaced) and power (linear) of one spectrum; fitted over [fLo, fHi] inclusive,
// skipping 0 Hz as specparam does.
SpecParamResult fitSpecParam (const std::vector<double>& freqs, const std::vector<double>& power,
                              double fLo, double fHi, const SpecParamSettings& settings = {});
} // namespace ProbeViewer::Localization::Realtime
