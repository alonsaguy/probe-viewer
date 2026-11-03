#include "IblLfFeatures.hpp"

#include "Cadzow.hpp"
#include "PyDsp.hpp"
#include "SpecParam.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace ProbeViewer::Localization::Realtime
{
Np1Geometry Np1Geometry::make (int nc)
{
    static const double xs[4] = { 43.0, 11.0, 59.0, 27.0 };
    static const int cols[4] = { 2, 0, 3, 1 };
    Np1Geometry g;
    for (int c = 0; c < nc; ++c)
    {
        g.x.push_back (xs[c % 4]);
        g.row.push_back (c / 2);
        g.y.push_back (20.0 + 20.0 * (c / 2));
        g.col.push_back (cols[c % 4]);
    }
    return g;
}

namespace
{
double medianOf (std::vector<double> v)
{
    const size_t n = v.size(), m = n / 2;
    std::nth_element (v.begin(), v.begin() + static_cast<long> (m), v.end());
    const double hi = v[m];
    if (n % 2 == 1)
        return hi;
    return 0.5 * (hi + *std::max_element (v.begin(), v.begin() + static_cast<long> (m)));
}

// ibldsp.voltage.current_source_density(lfp, h, n=1, method='diff', sigma=1/3, scale=True).
std::vector<double> firstSpatialDifference (const std::vector<double>& lfp, int nc, int ns, const Np1Geometry& g)
{
    std::vector<double> csd (static_cast<size_t> (nc) * static_cast<size_t> (ns), std::numeric_limits<double>::quiet_NaN());
    const double sigma = 1.0 / 3.0;
    std::vector<int> colsUnique (g.col.begin(), g.col.end());
    std::sort (colsUnique.begin(), colsUnique.end());
    colsUnique.erase (std::unique (colsUnique.begin(), colsUnique.end()), colsUnique.end());
    for (int col : colsUnique)
    {
        std::vector<int> itr;
        for (int c = 0; c < nc; ++c)
            if (g.col[static_cast<size_t> (c)] == col)
                itr.push_back (c);
        std::stable_sort (itr.begin(), itr.end(), [&g] (int a, int b) { return g.row[static_cast<size_t> (a)] < g.row[static_cast<size_t> (b)]; });
        if (itr.size() < 3)
            continue;
        std::vector<double> d;
        for (size_t i = 0; i + 1 < itr.size(); ++i)
        {
            const double a = std::hypot (g.x[static_cast<size_t> (itr[i])], g.y[static_cast<size_t> (itr[i])]) / 1e6;
            const double b = std::hypot (g.x[static_cast<size_t> (itr[i + 1])], g.y[static_cast<size_t> (itr[i + 1])]) / 1e6;
            d.push_back (b - a);
        }
        const double dx = medianOf (d);
        for (size_t i = 0; i + 1 < itr.size(); ++i)
        {
            const double* lo = &lfp[static_cast<size_t> (itr[i]) * static_cast<size_t> (ns)];
            const double* hi = &lfp[static_cast<size_t> (itr[i + 1]) * static_cast<size_t> (ns)];
            double* out = &csd[static_cast<size_t> (itr[i]) * static_cast<size_t> (ns)];
            for (int t = 0; t < ns; ++t)
                out[t] = (hi[t] - lo[t]) / dx * sigma;
        }
        // csd[itr[0]] = csd[itr[1]]; csd[itr[-1]] = csd[itr[-2]]
        std::copy_n (&csd[static_cast<size_t> (itr[1]) * static_cast<size_t> (ns)], ns, &csd[static_cast<size_t> (itr[0]) * static_cast<size_t> (ns)]);
        std::copy_n (&csd[static_cast<size_t> (itr[itr.size() - 2]) * static_cast<size_t> (ns)], ns, &csd[static_cast<size_t> (itr.back()) * static_cast<size_t> (ns)]);
    }
    return csd;
}
} // namespace

BprFeatures bprFeatures (const std::vector<double>& lf, int nc, int ns, double fs, const Np1Geometry& g, int threads)
{
    constexpr int q = 10;
    // FIR decimation of every channel.
    const int nsr = ns / q + (ns % q ? 1 : 0);
    std::vector<double> rs (static_cast<size_t> (nc) * static_cast<size_t> (nsr));
    for (int c = 0; c < nc; ++c)
    {
        const auto y = decimateFir (&lf[static_cast<size_t> (c) * static_cast<size_t> (ns)], ns, q);
        std::copy (y.begin(), y.end(), rs.begin() + static_cast<long> (c) * nsr);
    }
    const double fsr = fs / q;
    CadzowParams params;
    params.threads = threads;
    const std::vector<double> den = cadzowDenoise (rs, nc, nsr, g.x, g.y, fsr, params);
    const std::vector<double> csd = firstSpatialDifference (den, nc, nsr, g);

    BprFeatures out;
    out.rmsDb.assign (static_cast<size_t> (nc), std::numeric_limits<double>::quiet_NaN());
    for (auto& b : out.bandDb)
        b.assign (static_cast<size_t> (nc), std::numeric_limits<double>::quiet_NaN());
    std::vector<double> freqs, psd;
    for (int c = 0; c < nc; ++c)
    {
        const double* x = &csd[static_cast<size_t> (c) * static_cast<size_t> (nsr)];
        if (! std::isfinite (x[0]))
            continue;
        double s2 = 0.0;
        for (int t = 0; t < nsr; ++t)
            s2 += x[t] * x[t];
        out.rmsDb[static_cast<size_t> (c)] = 20.0 * std::log10 (std::sqrt (s2 / nsr));
        periodogram (x, nsr, fsr, freqs, psd);
        for (size_t b = 0; b < kBands.size(); ++b)
            out.bandDb[b][static_cast<size_t> (c)] = powerInBandDb (freqs, psd, kBands[b][0], kBands[b][1]);
    }
    return out;
}

DecayFeatures psdDecayFeatures (const double* x, int ns, double fs,
                                const std::vector<double>& periodFreqs, const std::vector<double>& period)
{
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    DecayFeatures out { nan, nan, nan, nan, nan, { nan, nan, nan, nan, nan, nan } };
    std::vector<double> freqs, psd;
    welch (x, ns, fs, 2048, freqs, psd);
    if (std::accumulate (psd.begin(), psd.end(), 0.0) == 0.0)
        return out;
    const auto fit = fitSpecParam (freqs, psd, kBands[5][0], kBands[5][1]);
    if (! fit.ok)
        return out;
    out.offset = fit.offset;
    out.exponent = fit.exponent;
    out.fitError = fit.error;
    out.rSquared = fit.rSquared;
    out.nPeaks = static_cast<double> (fit.gaussians.size());
    // residual = period / 10^(offset - exponent log10 f), with fscale[0] == 0 dropped.
    std::vector<double> residual (period.size(), 0.0);
    for (size_t k = 1; k < period.size(); ++k)
        residual[k] = period[k] / std::pow (10.0, fit.offset - fit.exponent * std::log10 (periodFreqs[k]));
    for (size_t b = 0; b < kBands.size(); ++b)
        out.residualDb[b] = powerInBandDb (periodFreqs, residual, kBands[b][0], kBands[b][1], 1);
    return out;
}
} // namespace ProbeViewer::Localization::Realtime
