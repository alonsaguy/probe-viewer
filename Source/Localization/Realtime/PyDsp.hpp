#pragma once

// Double-precision signal processing that reproduces the scipy / ibldsp functions used by the
// IBL feature pipeline (ephysatlas.features), so that C++ features match the Python ones:
// FFT (any length with small factors), firwin, scipy.signal.decimate (FIR, zero phase),
// periodogram, welch, and ibldsp's raised-cosine band power.

#include <complex>
#include <vector>

namespace ProbeViewer::Localization::Realtime
{
using cplx = std::complex<double>;

// Mixed-radix FFT (Cooley-Tukey on factors 2, 3, 4, 5, generic DFT for other factors).
class Fft
{
public:
    explicit Fft (int n);
    int size() const { return n; }
    // In-place forward (sign -1) or inverse without 1/n scaling (sign +1).
    void transform (std::vector<cplx>& data, bool inverse) const;
    // Real input of length n -> n/2 + 1 complex bins (numpy.fft.rfft).
    std::vector<cplx> rfft (const double* x) const;
    // n/2 + 1 bins -> real signal of length n (numpy.fft.irfft, with 1/n scaling).
    std::vector<double> irfft (const std::vector<cplx>& bins) const;

private:
    int n;
    std::vector<int> factors;
    std::vector<cplx> twiddles; // exp(-2 pi i k / n)
    void recurse (const cplx* in, cplx* out, int stride, int m, size_t factorIdx, bool inverse, std::vector<cplx>& scratch) const;
};

// scipy.signal.firwin(numtaps, cutoff, window='hamming') for a low-pass, cutoff relative to Nyquist.
std::vector<double> firwinHamming (int numtaps, double cutoff);

// scipy.signal.decimate(x, q, ftype='fir', zero_phase=True): ceil(n / q) samples.
std::vector<double> decimateFir (const double* x, int n, int q);

// scipy.signal.periodogram(x, fs): boxcar window, constant detrend, one-sided density.
void periodogram (const double* x, int n, double fs, std::vector<double>& freqs, std::vector<double>& psd);

// scipy.signal.welch(x, fs, nperseg): periodic Hann, 50 % overlap, constant detrend, density.
void welch (const double* x, int n, double fs, int nperseg, std::vector<double>& freqs, std::vector<double>& psd);

// ephysatlas.features._get_power_in_band: raised-cosine weights centred on the band (half-width =
// band width), 10 log10 of the weighted mean of psd over all given bins.
double powerInBandDb (const std::vector<double>& freqs, const std::vector<double>& psd, double lo, double hi, size_t firstBin = 0);

// Smallest 2,3,5-smooth integer >= n (scipy.fft.next_fast_len(n, real=True)).
int nextFastLen (int n);
} // namespace ProbeViewer::Localization::Realtime
