#pragma once

// Port of ibldsp.cadzow.cadzow_denoiser (F-X Cadzow rank reduction over channel windows) with
// the LFP settings of ephysatlas.features.cadzow_lfp_kwargs for NP1, i.e.
// ibldsp.cadzow.LFP_KWARGS_NP1 (ibl-neuropixel >= 1.14): the channels placed on the full
// lateral x depth grid with the empty positions filled, Gavish-Donoho singular-value shrinkage
// (scale 2.25) capped at rank 5, every frequency bin, nswx=64, ovx=32, ppca_k=2.0, npad=0.
// Only the npad = 0, ovx <= nswx / 2 path is implemented.

#include <vector>

namespace ProbeViewer::Localization::Realtime
{
struct CadzowParams
{
    int rank = 5;
    int niter = 1;
    double fmax = 0.0;         // Hz; <= 0 processes every bin (fmax=None)
    int nswx = 64;
    int ovx = 32;
    double gapThreshold = 0.0; // hard rank cut only: > 0 enables the adaptive rank
    double ppcaK = 2.0;        // <= 0 disables the outlier correction
    double shrinkage = 2.25;   // Gavish-Donoho noise scale; <= 0 uses the hard rank cut
    bool fillGrid = true;      // run on the full grid of unique (x, y) positions
    int threads = 8;
};

// wav: [nc x ns] row-major (channel-major) samples; x, y: channel positions (µm). Returns the
// denoised samples, rounded through float32 like the Python implementation.
std::vector<double> cadzowDenoise (const std::vector<double>& wav, int nc, int ns,
                                   const std::vector<double>& x, const std::vector<double>& y,
                                   double fs, const CadzowParams& params);

// Median of the Marchenko-Pastur distribution of aspect ratio beta (unit variance).
double marchenkoPasturMedian (double beta);
} // namespace ProbeViewer::Localization::Realtime
