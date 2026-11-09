#include "BadChannels.hpp"

#include "PyDsp.hpp"
#include "Threads.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <thread>

namespace ProbeViewer::Localization::Realtime
{
namespace
{
template <typename Fn>
void parallelChunks (int n, int threads, Fn&& fn)
{
    threads = std::max (1, std::min (threads, n));
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; ++t)
        pool.emplace_back ([&, t] {
            lowerThreadPriority();
            fn (n * t / threads, n * (t + 1) / threads);
        });
    for (auto& th : pool)
        th.join();
}

// channels_similarity: zero-lag correlation of each (demeaned) channel with the median trace,
// normalised by the trace's own energy. offsets are subtracted before taking the median.
std::vector<double> channelsSimilarity (const float* data, int nc, int ns, const std::vector<double>& offsets, int threads)
{
    std::vector<double> ref (static_cast<size_t> (ns));
    parallelChunks (ns, threads, [&] (int t0, int t1)
    {
        std::vector<double> v (static_cast<size_t> (nc));
        for (int t = t0; t < t1; ++t)
        {
            for (int c = 0; c < nc; ++c)
                v[static_cast<size_t> (c)] = static_cast<double> (data[static_cast<size_t> (c) * static_cast<size_t> (ns) + static_cast<size_t> (t)]) - offsets[static_cast<size_t> (c)];
            const size_t m = v.size() / 2;
            std::nth_element (v.begin(), v.begin() + static_cast<long> (m), v.end());
            double med = v[m];
            if (v.size() % 2 == 0)
                med = 0.5 * (med + *std::max_element (v.begin(), v.begin() + static_cast<long> (m)));
            ref[static_cast<size_t> (t)] = med;
        }
    });
    double refMean = 0.0;
    for (double r : ref)
        refMean += r;
    refMean /= ns;
    double apeak = 0.0;
    for (auto& r : ref)
    {
        r -= refMean;
        apeak += r * r;
    }
    std::vector<double> xcor (static_cast<size_t> (nc));
    parallelChunks (nc, threads, [&] (int c0, int c1)
    {
        for (int c = c0; c < c1; ++c)
        {
            const float* x = data + static_cast<size_t> (c) * static_cast<size_t> (ns);
            double mean = 0.0;
            for (int t = 0; t < ns; ++t)
                mean += x[t];
            mean /= ns;
            double s = 0.0;
            for (int t = 0; t < ns; ++t)
                s += (x[t] - mean) * ref[static_cast<size_t> (t)];
            xcor[static_cast<size_t> (c)] = s / apeak;
        }
    });
    return xcor;
}

// x - medfilt(x padded with ceil(nmed / 2) edge values, nmed).
std::vector<double> detrend (const std::vector<double>& x, int nmed)
{
    const int n = static_cast<int> (x.size()), ntap = (nmed + 1) / 2, half = nmed / 2;
    std::vector<double> xf (static_cast<size_t> (n + 2 * ntap));
    for (int i = 0; i < ntap; ++i)
    {
        xf[static_cast<size_t> (i)] = x.front();
        xf[static_cast<size_t> (n + ntap + i)] = x.back();
    }
    std::copy (x.begin(), x.end(), xf.begin() + ntap);
    std::vector<double> out (x.size()), w (static_cast<size_t> (nmed));
    for (int i = 0; i < n; ++i)
    {
        const int centre = i + ntap;
        for (int k = 0; k < nmed; ++k)
        {
            const int j = centre - half + k;
            w[static_cast<size_t> (k)] = (j >= 0 && j < static_cast<int> (xf.size())) ? xf[static_cast<size_t> (j)] : 0.0;
        }
        std::nth_element (w.begin(), w.begin() + half, w.end());
        out[static_cast<size_t> (i)] = x[static_cast<size_t> (i)] - w[static_cast<size_t> (half)];
    }
    return out;
}
} // namespace

BadChannelRawStats badChannelRawStats (const float* raw, int nc, int ns, double fs, int threads)
{
    BadChannelRawStats out;
    std::vector<double> means (static_cast<size_t> (nc));
    for (int c = 0; c < nc; ++c)
    {
        const float* x = raw + static_cast<size_t> (c) * static_cast<size_t> (ns);
        double s = 0.0;
        for (int t = 0; t < ns; ++t)
            s += x[t];
        means[static_cast<size_t> (c)] = s / ns;
    }
    out.xcor = channelsSimilarity (raw, nc, ns, means, threads);
    out.psdHf.assign (static_cast<size_t> (nc), 0.0);
    parallelChunks (nc, threads, [&] (int c0, int c1)
    {
        std::vector<double> x (static_cast<size_t> (ns)), freqs, psd;
        for (int c = c0; c < c1; ++c)
        {
            const float* r = raw + static_cast<size_t> (c) * static_cast<size_t> (ns);
            for (int t = 0; t < ns; ++t)
                x[static_cast<size_t> (t)] = (r[t] - means[static_cast<size_t> (c)]) * 1e6;
            welch (x.data(), ns, fs, 256, freqs, psd);
            double s = 0.0;
            int n = 0;
            for (size_t k = 0; k < freqs.size(); ++k)
                if (freqs[k] > fs / 2.0 * 0.8)
                {
                    s += psd[k];
                    ++n;
                }
            out.psdHf[static_cast<size_t> (c)] = n > 0 ? s / n : 0.0;
        }
    });
    return out;
}

std::vector<int> badChannelLabels (const BadChannelRawStats& rawStats, const float* hf, int nc, int ns, int threads)
{
    const std::vector<double> xcorf = channelsSimilarity (hf, nc, ns, std::vector<double> (static_cast<size_t> (nc), 0.0), threads);
    const std::vector<double> xcorHf = detrend (rawStats.xcor, 11);
    const std::vector<double> dtf = detrend (xcorf, 11);
    std::vector<int> labels (static_cast<size_t> (nc), 0);
    // Outside the brain: xcor_lf < -0.75 on a run of channels reaching the top of the probe.
    std::vector<int> outside;
    for (int c = 0; c < nc; ++c)
        if (xcorf[static_cast<size_t> (c)] - dtf[static_cast<size_t> (c)] - 1.0 < -0.75)
            outside.push_back (c);
    if (! outside.empty() && outside.back() == nc - 1)
    {
        size_t first = outside.size() - 1;
        while (first > 0 && outside[first - 1] == outside[first] - 1)
            --first;
        for (size_t i = first; i < outside.size(); ++i)
            labels[static_cast<size_t> (outside[i])] = 3;
    }
    for (int c = 0; c < nc; ++c)
        if (-0.5 > xcorHf[static_cast<size_t> (c)])
            labels[static_cast<size_t> (c)] = 1;
    for (int c = 0; c < nc; ++c)
        if (rawStats.psdHf[static_cast<size_t> (c)] > 0.02 || xcorHf[static_cast<size_t> (c)] > 1.0)
            labels[static_cast<size_t> (c)] = 2;
    return labels;
}

void interpolateBadChannels (float* data, int nc, int ns, const std::vector<int>& labels,
                             const std::vector<double>& x, const std::vector<double>& y)
{
    std::vector<int> bad;
    for (int c = 0; c < nc; ++c)
        if (labels[static_cast<size_t> (c)] == 1 || labels[static_cast<size_t> (c)] == 2)
            bad.push_back (c);
    std::vector<double> w (static_cast<size_t> (nc)), acc (static_cast<size_t> (ns));
    for (int i : bad)
    {
        double wsum = 0.0;
        for (int c = 0; c < nc; ++c)
        {
            const double offset = std::abs (std::complex<double> (x[static_cast<size_t> (c)] - x[static_cast<size_t> (i)], y[static_cast<size_t> (c)] - y[static_cast<size_t> (i)]));
            double wc = std::exp (-std::pow (offset / 20.0, 1.3));
            if (labels[static_cast<size_t> (c)] == 1 || labels[static_cast<size_t> (c)] == 2 || wc < 0.005)
                wc = 0.0;
            w[static_cast<size_t> (c)] = wc;
            wsum += wc;
        }
        float* out = data + static_cast<size_t> (i) * static_cast<size_t> (ns);
        if (wsum == 0.0)
        {
            std::fill (out, out + ns, 0.0f);
            continue;
        }
        std::fill (acc.begin(), acc.end(), 0.0);
        for (int c = 0; c < nc; ++c)
        {
            const double wc = w[static_cast<size_t> (c)] / wsum;
            if (! (wc > 0.005)) // the normalised weights are thresholded again, not renormalised
                continue;
            const float* src = data + static_cast<size_t> (c) * static_cast<size_t> (ns);
            for (int t = 0; t < ns; ++t)
                acc[static_cast<size_t> (t)] += wc * src[t];
        }
        for (int t = 0; t < ns; ++t)
            out[t] = static_cast<float> (acc[static_cast<size_t> (t)]);
    }
}
} // namespace ProbeViewer::Localization::Realtime
