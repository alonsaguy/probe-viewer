#include "PyDsp.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ProbeViewer::Localization::Realtime
{
namespace
{
constexpr double kPi = 3.14159265358979323846;
}

// ---------------------------------------------------------------------------------------------
// FFT
// ---------------------------------------------------------------------------------------------

Fft::Fft (int size) : n (size)
{
    if (n <= 0)
        throw std::runtime_error ("FFT size must be positive");
    int m = n;
    for (int p : { 4, 2, 3, 5 })
        while (m % p == 0)
        {
            factors.push_back (p);
            m /= p;
        }
    for (int p = 7; m > 1; p += 2)
        while (m % p == 0)
        {
            factors.push_back (p);
            m /= p;
        }
    for (int p : factors)
        maxFactor = std::max (maxFactor, p);
    twiddles.resize (static_cast<size_t> (n));
    for (int k = 0; k < n; ++k)
        twiddles[static_cast<size_t> (k)] = std::polar (1.0, -2.0 * kPi * k / n);
}

void Fft::prepare (Work& work) const
{
    // One stage's p values per factor (the stages of one recursion path are live together).
    const size_t scratch = static_cast<size_t> (maxFactor) * (factors.size() + 1);
    if (work.in.size() != static_cast<size_t> (n))
        work.in.resize (static_cast<size_t> (n));
    if (work.out.size() != static_cast<size_t> (n))
        work.out.resize (static_cast<size_t> (n));
    if (work.scratch.size() < scratch)
        work.scratch.resize (scratch);
}

void Fft::recurse (const cplx* in, cplx* out, int stride, int m, size_t f, bool inverse, cplx* scratch) const
{
    if (m == 1)
    {
        out[0] = in[0];
        return;
    }
    const int p = factors[f];
    const int mm = m / p;
    for (int q = 0; q < p; ++q)
        recurse (in + static_cast<ptrdiff_t> (q) * stride, out + static_cast<ptrdiff_t> (q) * mm, stride * p, mm, f + 1, inverse, scratch);
    const int twStep = n / m; // twiddle index step for this stage
    cplx* t = scratch + f * static_cast<size_t> (maxFactor);
    for (int k = 0; k < mm; ++k)
    {
        for (int q = 0; q < p; ++q)
        {
            cplx w = twiddles[static_cast<size_t> ((static_cast<long long> (q) * k * twStep) % n)];
            if (inverse)
                w = std::conj (w);
            t[q] = out[q * mm + k] * w;
        }
        for (int r = 0; r < p; ++r)
        {
            cplx s = 0.0;
            for (int q = 0; q < p; ++q)
            {
                cplx w = twiddles[static_cast<size_t> ((static_cast<long long> (q) * r * (n / p)) % n)];
                if (inverse)
                    w = std::conj (w);
                s += t[q] * w;
            }
            out[r * mm + k] = s;
        }
    }
}

void Fft::transform (std::vector<cplx>& data, bool inverse) const
{
    if (static_cast<int> (data.size()) != n)
        throw std::runtime_error ("FFT input has the wrong length");
    std::vector<cplx> out (static_cast<size_t> (n));
    std::vector<cplx> scratch (static_cast<size_t> (maxFactor) * (factors.size() + 1));
    recurse (data.data(), out.data(), 1, n, 0, inverse, scratch.data());
    data.swap (out);
}

void Fft::rfft (const double* x, Work& work) const
{
    prepare (work);
    for (int i = 0; i < n; ++i)
        work.in[static_cast<size_t> (i)] = cplx (x[i], 0.0);
    recurse (work.in.data(), work.out.data(), 1, n, 0, false, work.scratch.data());
}

std::vector<cplx> Fft::rfft (const double* x) const
{
    std::vector<cplx> d (static_cast<size_t> (n));
    for (int i = 0; i < n; ++i)
        d[static_cast<size_t> (i)] = cplx (x[i], 0.0);
    transform (d, false);
    d.resize (static_cast<size_t> (n / 2 + 1));
    return d;
}

std::vector<double> Fft::irfft (const std::vector<cplx>& bins) const
{
    std::vector<cplx> d (static_cast<size_t> (n));
    const int half = n / 2;
    for (int k = 0; k <= half && k < static_cast<int> (bins.size()); ++k)
        d[static_cast<size_t> (k)] = bins[static_cast<size_t> (k)];
    d[0] = cplx (d[0].real(), 0.0);
    if (n % 2 == 0)
        d[static_cast<size_t> (half)] = cplx (d[static_cast<size_t> (half)].real(), 0.0);
    for (int k = 1; k < (n + 1) / 2; ++k)
        d[static_cast<size_t> (n - k)] = std::conj (d[static_cast<size_t> (k)]);
    transform (d, true);
    std::vector<double> out (static_cast<size_t> (n));
    for (int i = 0; i < n; ++i)
        out[static_cast<size_t> (i)] = d[static_cast<size_t> (i)].real() / n;
    return out;
}

// ---------------------------------------------------------------------------------------------
// Filters and spectra
// ---------------------------------------------------------------------------------------------

std::vector<double> firwinHamming (int numtaps, double cutoff)
{
    std::vector<double> h (static_cast<size_t> (numtaps));
    const double alpha = 0.5 * (numtaps - 1);
    double sum = 0.0;
    for (int i = 0; i < numtaps; ++i)
    {
        const double m = i - alpha;
        const double x = cutoff * m;
        const double sinc = x == 0.0 ? 1.0 : std::sin (kPi * x) / (kPi * x);
        const double win = 0.54 - 0.46 * std::cos (2.0 * kPi * i / (numtaps - 1));
        h[static_cast<size_t> (i)] = cutoff * sinc * win;
        sum += h[static_cast<size_t> (i)];
    }
    for (auto& v : h)
        v /= sum;
    return h;
}

std::vector<double> decimateFir (const double* x, int n, int q)
{
    const int halfLen = 10 * q;
    const std::vector<double> b = firwinHamming (2 * halfLen + 1, 1.0 / q);
    const int nOut = n / q + (n % q ? 1 : 0);
    std::vector<double> y (static_cast<size_t> (nOut), 0.0);
    for (int k = 0; k < nOut; ++k)
    {
        const int centre = q * k + halfLen;
        double s = 0.0;
        const int iLo = std::max (0, centre - (n - 1));
        const int iHi = std::min (2 * halfLen, centre);
        for (int i = iLo; i <= iHi; ++i)
            s += b[static_cast<size_t> (i)] * x[centre - i];
        y[static_cast<size_t> (k)] = s;
    }
    return y;
}

void periodogram (const double* x, int n, double fs, std::vector<double>& freqs, std::vector<double>& psd)
{
    double mean = 0.0;
    for (int i = 0; i < n; ++i)
        mean += x[i];
    mean /= n;
    std::vector<double> d (x, x + n);
    for (auto& v : d)
        v -= mean;
    const Fft fft (n);
    const auto X = fft.rfft (d.data());
    const int nb = n / 2 + 1;
    freqs.resize (static_cast<size_t> (nb));
    psd.resize (static_cast<size_t> (nb));
    for (int k = 0; k < nb; ++k)
    {
        freqs[static_cast<size_t> (k)] = k * fs / n;
        double p = std::norm (X[static_cast<size_t> (k)]) / (fs * n);
        const bool nyquist = (n % 2 == 0) && k == nb - 1;
        if (k > 0 && ! nyquist)
            p *= 2.0;
        psd[static_cast<size_t> (k)] = p;
    }
}

void welch (const double* x, int n, double fs, int nperseg, std::vector<double>& freqs, std::vector<double>& psd)
{
    nperseg = std::min (nperseg, n);
    const int step = nperseg - nperseg / 2;
    const int nseg = (n - nperseg) / step + 1;
    std::vector<double> win (static_cast<size_t> (nperseg));
    double scale = 0.0;
    for (int i = 0; i < nperseg; ++i)
    {
        win[static_cast<size_t> (i)] = 0.5 - 0.5 * std::cos (2.0 * kPi * i / nperseg);
        scale += win[static_cast<size_t> (i)] * win[static_cast<size_t> (i)];
    }
    const Fft fft (nperseg);
    const int nb = nperseg / 2 + 1;
    freqs.resize (static_cast<size_t> (nb));
    psd.assign (static_cast<size_t> (nb), 0.0);
    std::vector<double> seg (static_cast<size_t> (nperseg));
    Fft::Work work;
    for (int s = 0; s < nseg; ++s)
    {
        const double* p = x + static_cast<ptrdiff_t> (s) * step;
        double mean = 0.0;
        for (int i = 0; i < nperseg; ++i)
            mean += p[i];
        mean /= nperseg;
        for (int i = 0; i < nperseg; ++i)
            seg[static_cast<size_t> (i)] = (p[i] - mean) * win[static_cast<size_t> (i)];
        fft.rfft (seg.data(), work);
        for (int k = 0; k < nb; ++k)
            psd[static_cast<size_t> (k)] += std::norm (work.out[static_cast<size_t> (k)]);
    }
    for (int k = 0; k < nb; ++k)
    {
        freqs[static_cast<size_t> (k)] = k * fs / nperseg;
        double v = psd[static_cast<size_t> (k)] / nseg / (fs * scale);
        const bool nyquist = (nperseg % 2 == 0) && k == nb - 1;
        if (k > 0 && ! nyquist)
            v *= 2.0;
        psd[static_cast<size_t> (k)] = v;
    }
}

double powerInBandDb (const std::vector<double>& freqs, const std::vector<double>& psd, double lo, double hi, size_t firstBin)
{
    const double width = hi - lo, centre = 0.5 * (lo + hi);
    double num = 0.0, den = 0.0;
    for (size_t k = firstBin; k < freqs.size(); ++k)
    {
        const double x = -std::abs (freqs[k] - centre);
        double w;
        if (x <= -width)
            w = 0.0;
        else if (x >= 0.0)
            w = 1.0;
        else
            w = (1.0 - std::cos ((x + width) / width * kPi)) / 2.0;
        num += psd[k] * w;
        den += w;
    }
    return 10.0 * std::log10 (num / den);
}

int nextFastLen (int n)
{
    for (int m = std::max (1, n);; ++m)
    {
        int k = m;
        for (int p : { 2, 3, 5 })
            while (k % p == 0)
                k /= p;
        if (k == 1)
            return m;
    }
}
} // namespace ProbeViewer::Localization::Realtime
