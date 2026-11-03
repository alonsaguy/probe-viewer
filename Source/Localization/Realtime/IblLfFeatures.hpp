#pragma once

// LF-band features of the IBL feature pipeline that work across channels or need spectral
// models, reproduced from ephysatlas.features: the BPR (bipolar referencing) family of csd() and the
// aperiodic / periodic / residual features of get_psd_decay_features().

#include <array>
#include <vector>

namespace ProbeViewer::Localization::Realtime
{
// neuropixel.trace_header(version=1) for the first nc channels.
struct Np1Geometry
{
    std::vector<double> x, y; // µm
    std::vector<int> col, row;
    static Np1Geometry make (int nc);
};

// ephysatlas.features.BANDS order: delta, theta, alpha, beta, gamma, lfp.
constexpr std::array<std::array<double, 2>, 6> kBands { { { 0.0, 4.0 }, { 4.0, 10.0 }, { 8.0, 12.0 }, { 15.0, 30.0 }, { 30.0, 90.0 }, { 0.0, 90.0 } } };

struct BprFeatures
{
    std::vector<double> rmsDb;                  // [nc] rms_lf_bpr: 20 log10 rms of the BPR trace
    std::array<std::vector<double>, 6> bandDb;  // [band][nc] psd_<band>_bpr (dB)
};

// csd(data, fs, geometry, decimate=10)[*_bpr] of destriped LF (volts, [nc x ns]): FIR decimation
// by 10, Cadzow denoising (cadzow_lfp_kwargs, see Cadzow.hpp), per-column first spatial
// difference scaled by sigma / dx (current_source_density n=1), then rms and periodogram band
// powers at fs / 10.
BprFeatures bprFeatures (const std::vector<double>& lf, int nc, int ns, double fs, const Np1Geometry& geometry, int threads);

struct DecayFeatures
{
    double offset, exponent, fitError, rSquared, nPeaks; // NaN when the spectrum is all zero
    std::array<double, 6> residualDb;                    // psd_residual_<band>, kBands order
};

// get_psd_decay_features for one channel: specparam on the Welch PSD (nperseg 2048) over
// [0, 90] Hz, then band powers of the periodogram divided by the fitted aperiodic curve
// (0 Hz excluded). periodFreqs / period: scipy.signal.periodogram of the same samples.
DecayFeatures psdDecayFeatures (const double* x, int ns, double fs,
                                const std::vector<double>& periodFreqs, const std::vector<double>& period);
} // namespace ProbeViewer::Localization::Realtime
