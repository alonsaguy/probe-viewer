#include "SpecParam.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>

namespace ProbeViewer::Localization::Realtime
{
namespace
{
struct Line
{
    double offset, exponent;
};

// Least-squares fit of y = offset - exponent * log10(f) on the selected points.
std::optional<Line> fitAperiodic (const std::vector<double>& logf, const std::vector<double>& y, const std::vector<char>* mask = nullptr)
{
    double n = 0.0, sx = 0.0, sy = 0.0;
    for (size_t i = 0; i < y.size(); ++i)
        if (! mask || (*mask)[i])
        {
            n += 1.0;
            sx += logf[i];
            sy += y[i];
        }
    if (n < 2.0)
        return std::nullopt;
    const double mx = sx / n, my = sy / n;
    double sxx = 0.0, sxy = 0.0;
    for (size_t i = 0; i < y.size(); ++i)
        if (! mask || (*mask)[i])
        {
            sxx += (logf[i] - mx) * (logf[i] - mx);
            sxy += (logf[i] - mx) * (y[i] - my);
        }
    if (! (sxx > 0.0))
        return std::nullopt;
    const double slope = sxy / sxx;
    return Line { my - slope * mx, -slope };
}

// numpy.percentile with linear interpolation.
double percentile (std::vector<double> v, double q)
{
    std::sort (v.begin(), v.end());
    const double pos = q / 100.0 * static_cast<double> (v.size() - 1);
    const size_t lo = static_cast<size_t> (std::floor (pos));
    const size_t hi = std::min (lo + 1, v.size() - 1);
    return v[lo] + (pos - static_cast<double> (lo)) * (v[hi] - v[lo]);
}

double populationStd (const std::vector<double>& v)
{
    const double mean = std::accumulate (v.begin(), v.end(), 0.0) / static_cast<double> (v.size());
    double s = 0.0;
    for (double x : v)
        s += (x - mean) * (x - mean);
    return std::sqrt (s / static_cast<double> (v.size()));
}

double gaussian (double f, double ctr, double hgt, double wid)
{
    return hgt * std::exp (-(f - ctr) * (f - ctr) / (2.0 * wid * wid));
}

std::vector<double> peakCurve (const std::vector<double>& freqs, const std::vector<double>& p)
{
    std::vector<double> y (freqs.size(), 0.0);
    for (size_t g = 0; g + 2 < p.size(); g += 3)
        for (size_t i = 0; i < freqs.size(); ++i)
            y[i] += gaussian (freqs[i], p[g], p[g + 1], p[g + 2]);
    return y;
}

// Solves A x = b for symmetric positive definite A (n x n, row-major); false if not SPD.
bool choleskySolve (std::vector<double> a, std::vector<double>& b, size_t n)
{
    for (size_t j = 0; j < n; ++j)
    {
        double d = a[j * n + j];
        for (size_t k = 0; k < j; ++k)
            d -= a[j * n + k] * a[j * n + k];
        if (! (d > 0.0))
            return false;
        d = std::sqrt (d);
        a[j * n + j] = d;
        for (size_t i = j + 1; i < n; ++i)
        {
            double s = a[i * n + j];
            for (size_t k = 0; k < j; ++k)
                s -= a[i * n + k] * a[j * n + k];
            a[i * n + j] = s / d;
        }
    }
    for (size_t i = 0; i < n; ++i)
    {
        double s = b[i];
        for (size_t k = 0; k < i; ++k)
            s -= a[i * n + k] * b[k];
        b[i] = s / a[i * n + i];
    }
    for (size_t i = n; i-- > 0;)
    {
        double s = b[i];
        for (size_t k = i + 1; k < n; ++k)
            s -= a[k * n + i] * b[k];
        b[i] = s / a[i * n + i];
    }
    return true;
}

// Minimises 0.5 ||sum of gaussians(p) - y||^2 subject to lb <= p <= ub by Levenberg-Marquardt on
// the variables not held at a bound by the gradient, projecting each step onto the box.
std::vector<double> fitGaussians (const std::vector<double>& freqs, const std::vector<double>& y, std::vector<double> p,
                                  const std::vector<double>& lb, const std::vector<double>& ub)
{
    const size_t n = p.size(), m = freqs.size();
    auto cost = [&] (const std::vector<double>& q, std::vector<double>& r)
    {
        r = peakCurve (freqs, q);
        double c = 0.0;
        for (size_t i = 0; i < m; ++i)
        {
            r[i] -= y[i];
            c += r[i] * r[i];
        }
        return 0.5 * c;
    };
    for (size_t k = 0; k < n; ++k)
        p[k] = std::clamp (p[k], lb[k], ub[k]);
    std::vector<double> r, rNew, jac (m * n), grad (n), hess (n * n), pNew (n);
    double f = cost (p, r);
    double lambda = -1.0;
    for (int iter = 0; iter < 1000; ++iter)
    {
        for (size_t g = 0; g + 2 < n; g += 3)
            for (size_t i = 0; i < m; ++i)
            {
                const double dx = freqs[i] - p[g], w = p[g + 2];
                const double e = std::exp (-dx * dx / (2.0 * w * w));
                jac[i * n + g] = p[g + 1] * e * dx / (w * w);
                jac[i * n + g + 1] = e;
                jac[i * n + g + 2] = p[g + 1] * e * dx * dx / (w * w * w);
            }
        for (size_t a = 0; a < n; ++a)
        {
            double s = 0.0;
            for (size_t i = 0; i < m; ++i)
                s += jac[i * n + a] * r[i];
            grad[a] = s;
            for (size_t b = a; b < n; ++b)
            {
                double h = 0.0;
                for (size_t i = 0; i < m; ++i)
                    h += jac[i * n + a] * jac[i * n + b];
                hess[a * n + b] = hess[b * n + a] = h;
            }
        }
        std::vector<size_t> freeVars;
        double projGrad = 0.0;
        for (size_t k = 0; k < n; ++k)
        {
            const bool heldLow = p[k] <= lb[k] && grad[k] > 0.0;
            const bool heldHigh = p[k] >= ub[k] && grad[k] < 0.0;
            if (! heldLow && ! heldHigh)
            {
                freeVars.push_back (k);
                projGrad = std::max (projGrad, std::abs (grad[k]));
            }
        }
        if (freeVars.empty() || projGrad <= 1e-14 * (1.0 + f))
            break;
        if (lambda < 0.0)
        {
            double dmax = 0.0;
            for (size_t k = 0; k < n; ++k)
                dmax = std::max (dmax, hess[k * n + k]);
            lambda = 1e-3 * dmax;
        }
        const size_t nf = freeVars.size();
        bool accepted = false;
        double reduction = 0.0, stepNorm = 0.0;
        for (int inner = 0; inner < 60 && ! accepted; ++inner)
        {
            std::vector<double> a (nf * nf), step (nf);
            for (size_t i = 0; i < nf; ++i)
            {
                for (size_t j = 0; j < nf; ++j)
                    a[i * nf + j] = hess[freeVars[i] * n + freeVars[j]];
                a[i * nf + i] += lambda * std::max (hess[freeVars[i] * n + freeVars[i]], 1e-12);
                step[i] = -grad[freeVars[i]];
            }
            if (! choleskySolve (a, step, nf))
            {
                lambda = std::max (lambda * 10.0, 1e-12);
                continue;
            }
            pNew = p;
            for (size_t i = 0; i < nf; ++i)
            {
                const size_t k = freeVars[i];
                pNew[k] = std::clamp (p[k] + step[i], lb[k], ub[k]);
            }
            const double fNew = cost (pNew, rNew);
            if (fNew < f)
            {
                reduction = f - fNew;
                stepNorm = 0.0;
                for (size_t k = 0; k < n; ++k)
                    stepNorm += (pNew[k] - p[k]) * (pNew[k] - p[k]);
                stepNorm = std::sqrt (stepNorm);
                p.swap (pNew);
                r.swap (rNew);
                f = fNew;
                lambda = std::max (lambda / 3.0, 1e-15);
                accepted = true;
            }
            else
                lambda *= 4.0;
        }
        if (! accepted)
            break;
        double pNorm = 0.0;
        for (double v : p)
            pNorm += v * v;
        if (reduction <= 1e-13 * f || stepNorm <= 1e-12 * (1e-12 + std::sqrt (pNorm)))
            break;
    }
    return p;
}
} // namespace

SpecParamResult fitSpecParam (const std::vector<double>& freqsIn, const std::vector<double>& powerIn,
                              double fLo, double fHi, const SpecParamSettings& s)
{
    SpecParamResult out;
    std::vector<double> freqs, power;
    for (size_t i = 0; i < freqsIn.size() && i < powerIn.size(); ++i)
        if (freqsIn[i] >= fLo && freqsIn[i] <= fHi && freqsIn[i] != 0.0)
        {
            freqs.push_back (freqsIn[i]);
            power.push_back (std::log10 (powerIn[i]));
        }
    const size_t m = freqs.size();
    if (m < 3)
        return out;
    for (double v : power)
        if (! std::isfinite (v))
            return out;
    const double rangeLo = freqs.front(), rangeHi = freqs.back(), freqRes = freqs[1] - freqs[0];
    std::vector<double> logf (m);
    for (size_t i = 0; i < m; ++i)
        logf[i] = std::log10 (freqs[i]);
    auto apCurve = [&] (const Line& line)
    {
        std::vector<double> y (m);
        for (size_t i = 0; i < m; ++i)
            y[i] = line.offset - line.exponent * logf[i];
        return y;
    };

    // Robust aperiodic fit: refit on the points at or below a low percentile of the flattened
    // spectrum (negative residuals clipped to zero).
    const auto initial = fitAperiodic (logf, power);
    if (! initial)
        return out;
    std::vector<double> flat (m);
    {
        const auto ap0 = apCurve (*initial);
        for (size_t i = 0; i < m; ++i)
            flat[i] = std::max (power[i] - ap0[i], 0.0);
    }
    const double thresh = percentile (flat, s.apPercentileThresh);
    std::vector<char> mask (m);
    for (size_t i = 0; i < m; ++i)
        mask[i] = flat[i] <= thresh;
    const auto robust = fitAperiodic (logf, power, &mask);
    if (! robust)
        return out;
    const auto apRobust = apCurve (*robust);
    std::vector<double> spectrumFlat (m);
    for (size_t i = 0; i < m; ++i)
        spectrumFlat[i] = power[i] - apRobust[i];

    // Iterative peak search on the flattened spectrum.
    const double stdLo = s.peakWidthLo / 2.0, stdHi = s.peakWidthHi / 2.0;
    std::vector<std::array<double, 3>> guess;
    std::vector<double> flatIter = spectrumFlat;
    while (static_cast<int> (guess.size()) < s.maxNPeaks)
    {
        const size_t maxInd = static_cast<size_t> (std::max_element (flatIter.begin(), flatIter.end()) - flatIter.begin());
        const double maxHeight = flatIter[maxInd];
        if (maxHeight <= s.peakThreshold * populationStd (flatIter))
            break;
        if (! (maxHeight > s.minPeakHeight))
            break;
        const double half = 0.5 * maxHeight;
        std::optional<size_t> shortSide;
        for (size_t v = maxInd; v-- > 1;) // range(max_ind - 1, 0, -1): index 0 is never checked
            if (flatIter[v] <= half)
            {
                shortSide = maxInd - v;
                break;
            }
        for (size_t v = maxInd + 1; v < m; ++v)
            if (flatIter[v] <= half)
            {
                shortSide = shortSide ? std::min (*shortSide, v - maxInd) : v - maxInd;
                break;
            }
        double guessStd = shortSide ? static_cast<double> (*shortSide) * 2.0 * freqRes / (2.0 * std::sqrt (2.0 * std::log (2.0)))
                                    : 0.5 * (s.peakWidthLo + s.peakWidthHi);
        guessStd = std::clamp (guessStd, stdLo, stdHi);
        guess.push_back ({ freqs[maxInd], maxHeight, guessStd });
        for (size_t i = 0; i < m; ++i)
            flatIter[i] -= gaussian (freqs[i], freqs[maxInd], maxHeight, guessStd);
    }

    // Drop guesses too close to the edges, then the lower of each overlapping adjacent pair.
    std::vector<std::array<double, 3>> kept;
    for (const auto& g : guess)
        if (std::abs (g[0] - rangeLo) > g[2] * s.bwStdEdge && std::abs (g[0] - rangeHi) > g[2] * s.bwStdEdge)
            kept.push_back (g);
    std::stable_sort (kept.begin(), kept.end(), [] (const auto& a, const auto& b) { return a[0] < b[0]; });
    std::vector<char> drop (kept.size(), 0);
    for (size_t i = 0; i + 1 < kept.size(); ++i)
        if (kept[i][0] + kept[i][2] * s.gaussOverlapThresh > kept[i + 1][0] - kept[i + 1][2] * s.gaussOverlapThresh)
            drop[kept[i + 1][1] < kept[i][1] ? i + 1 : i] = 1;
    guess.clear();
    for (size_t i = 0; i < kept.size(); ++i)
        if (! drop[i])
            guess.push_back (kept[i]);

    std::vector<double> peakFit (m, 0.0);
    if (! guess.empty())
    {
        std::vector<double> p0, lb, ub;
        for (const auto& g : guess)
        {
            double cfLo = g[0] - 2.0 * s.cfBound * g[2], cfHi = g[0] + 2.0 * s.cfBound * g[2];
            if (! (cfLo > rangeLo))
                cfLo = rangeLo;
            if (! (cfHi < rangeHi))
                cfHi = rangeHi;
            p0.insert (p0.end(), { g[0], g[1], g[2] });
            lb.insert (lb.end(), { cfLo, 0.0, stdLo });
            ub.insert (ub.end(), { cfHi, std::numeric_limits<double>::infinity(), stdHi });
        }
        const auto p = fitGaussians (freqs, spectrumFlat, p0, lb, ub);
        for (size_t g = 0; g + 2 < p.size(); g += 3)
            out.gaussians.push_back ({ p[g], p[g + 1], p[g + 2] });
        std::stable_sort (out.gaussians.begin(), out.gaussians.end(), [] (const auto& a, const auto& b) { return a[0] < b[0]; });
        peakFit = peakCurve (freqs, p);
    }

    // Final aperiodic fit on the peak-removed spectrum, then the goodness of fit of the full model.
    std::vector<double> peakRemoved (m);
    for (size_t i = 0; i < m; ++i)
        peakRemoved[i] = power[i] - peakFit[i];
    const auto finalFit = fitAperiodic (logf, peakRemoved);
    if (! finalFit)
        return out;
    const auto apFinal = apCurve (*finalFit);
    std::vector<double> model (m);
    double mae = 0.0, mp = 0.0, mm = 0.0;
    for (size_t i = 0; i < m; ++i)
    {
        model[i] = peakFit[i] + apFinal[i];
        mae += std::abs (power[i] - model[i]);
        mp += power[i];
        mm += model[i];
    }
    mp /= static_cast<double> (m);
    mm /= static_cast<double> (m);
    double spp = 0.0, smm = 0.0, spm = 0.0;
    for (size_t i = 0; i < m; ++i)
    {
        spp += (power[i] - mp) * (power[i] - mp);
        smm += (model[i] - mm) * (model[i] - mm);
        spm += (power[i] - mp) * (model[i] - mm);
    }
    out.ok = true;
    out.offset = finalFit->offset;
    out.exponent = finalFit->exponent;
    out.error = mae / static_cast<double> (m);
    out.rSquared = spm * spm / (spp * smm);
    return out;
}
} // namespace ProbeViewer::Localization::Realtime
