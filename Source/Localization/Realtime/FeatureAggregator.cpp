#include "FeatureAggregator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <unordered_map>

namespace ProbeViewer::Localization::Realtime
{
namespace
{
constexpr double kNan = std::numeric_limits<double>::quiet_NaN();

struct Pt
{
    double x, y;
};

// Exact for the small-integer (um) coordinates of probe channels.
double orient (const Pt& a, const Pt& b, const Pt& c)
{
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

// > 0 when d lies inside the circumcircle of the counter-clockwise triangle abc.
double incircle (const Pt& a, const Pt& b, const Pt& c, const Pt& d)
{
    const double adx = a.x - d.x, ady = a.y - d.y, bdx = b.x - d.x, bdy = b.y - d.y, cdx = c.x - d.x, cdy = c.y - d.y;
    return (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy) + (bdx * bdx + bdy * bdy) * (cdx * ady - adx * cdy)
         + (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady);
}

using Tri = std::array<int, 3>;

// Delaunay triangulation (counter-clockwise triangles): a sweep along y builds a triangulation
// of the convex hull, then Lawson flips make it Delaunay. Co-circular quadrilaterals keep the
// sweep's diagonal.
std::vector<Tri> delaunay (const std::vector<Pt>& p)
{
    const int n = static_cast<int> (p.size());
    std::vector<int> order (static_cast<size_t> (n));
    std::iota (order.begin(), order.end(), 0);
    std::sort (order.begin(), order.end(), [&p] (int a, int b) { return p[static_cast<size_t> (a)].y != p[static_cast<size_t> (b)].y ? p[static_cast<size_t> (a)].y < p[static_cast<size_t> (b)].y : p[static_cast<size_t> (a)].x < p[static_cast<size_t> (b)].x; });
    auto P = [&p] (int i) -> const Pt& { return p[static_cast<size_t> (i)]; };
    std::vector<Tri> tris;
    if (n < 3)
        return tris;
    int j = 2;
    while (j < n && orient (P (order[0]), P (order[1]), P (order[static_cast<size_t> (j)])) == 0.0)
        ++j;
    if (j == n)
        return tris; // all collinear
    // Fan from the first off-line point to the collinear prefix (sorted along the line).
    const int apex = order[static_cast<size_t> (j)];
    const bool left = orient (P (order[0]), P (order[static_cast<size_t> (j - 1)]), P (apex)) > 0.0;
    std::vector<int> hull;
    for (int i = 0; i + 1 < j; ++i)
    {
        const int a = order[static_cast<size_t> (i)], b = order[static_cast<size_t> (i + 1)];
        tris.push_back (left ? Tri { a, b, apex } : Tri { b, a, apex });
    }
    if (left)
    {
        for (int i = 0; i < j; ++i)
            hull.push_back (order[static_cast<size_t> (i)]);
        hull.push_back (apex);
    }
    else
    {
        hull.push_back (apex);
        for (int i = j - 1; i >= 0; --i)
            hull.push_back (order[static_cast<size_t> (i)]);
    }
    // Each further point is outside the hull built so far: join it to the visible hull edges.
    for (int k = j + 1; k < n; ++k)
    {
        const int q = order[static_cast<size_t> (k)];
        const int h = static_cast<int> (hull.size());
        std::vector<char> visible (static_cast<size_t> (h));
        for (int i = 0; i < h; ++i)
            visible[static_cast<size_t> (i)] = orient (P (hull[static_cast<size_t> (i)]), P (hull[static_cast<size_t> ((i + 1) % h)]), P (q)) < 0.0;
        int start = -1;
        for (int i = 0; i < h; ++i)
            if (visible[static_cast<size_t> (i)] && ! visible[static_cast<size_t> ((i + h - 1) % h)])
                start = i;
        if (start < 0)
            continue; // collinear with the hull (not expected for a lexicographic sweep)
        int count = 0;
        while (visible[static_cast<size_t> ((start + count) % h)])
        {
            const int a = hull[static_cast<size_t> ((start + count) % h)], b = hull[static_cast<size_t> ((start + count + 1) % h)];
            tris.push_back (Tri { b, a, q });
            ++count;
        }
        // Replace the interior vertices of the visible chain by q.
        std::vector<int> next;
        next.reserve (static_cast<size_t> (h + 1));
        for (int i = 0; i <= h - count; ++i)
            next.push_back (hull[static_cast<size_t> ((start + count + i) % h)]);
        next.push_back (q);
        hull.swap (next);
    }

    // Lawson flips.
    auto key = [] (int u, int v) { return (static_cast<uint64_t> (static_cast<uint32_t> (u)) << 32) | static_cast<uint32_t> (v); };
    std::unordered_map<uint64_t, int> edgeTri;
    edgeTri.reserve (tris.size() * 4);
    std::vector<std::pair<int, int>> stack;
    for (int t = 0; t < static_cast<int> (tris.size()); ++t)
        for (int e = 0; e < 3; ++e)
        {
            const int u = tris[static_cast<size_t> (t)][static_cast<size_t> (e)], v = tris[static_cast<size_t> (t)][static_cast<size_t> ((e + 1) % 3)];
            edgeTri[key (u, v)] = t;
            stack.emplace_back (u, v);
        }
    auto third = [] (const Tri& t, int u, int v)
    {
        for (int w : t)
            if (w != u && w != v)
                return w;
        return -1;
    };
    while (! stack.empty())
    {
        const auto [u, v] = stack.back();
        stack.pop_back();
        const auto it1 = edgeTri.find (key (u, v));
        const auto it2 = edgeTri.find (key (v, u));
        if (it1 == edgeTri.end() || it2 == edgeTri.end())
            continue;
        const int t1 = it1->second, t2 = it2->second;
        const int w = third (tris[static_cast<size_t> (t1)], u, v), z = third (tris[static_cast<size_t> (t2)], u, v);
        if (! (incircle (P (u), P (v), P (w), P (z)) > 0.0))
            continue;
        for (auto e : { key (u, v), key (v, w), key (w, u), key (v, u), key (u, z), key (z, v) })
            edgeTri.erase (e);
        tris[static_cast<size_t> (t1)] = Tri { u, z, w };
        tris[static_cast<size_t> (t2)] = Tri { z, v, w };
        edgeTri[key (u, z)] = t1;
        edgeTri[key (z, w)] = t1;
        edgeTri[key (w, u)] = t1;
        edgeTri[key (z, v)] = t2;
        edgeTri[key (v, w)] = t2;
        edgeTri[key (w, z)] = t2;
        stack.emplace_back (u, z);
        stack.emplace_back (z, v);
        stack.emplace_back (v, w);
        stack.emplace_back (w, u);
    }
    return tris;
}

double medianOf (std::vector<double> v)
{
    const size_t m = v.size() / 2;
    std::nth_element (v.begin(), v.begin() + static_cast<long> (m), v.end());
    const double hi = v[m];
    if (v.size() % 2 == 1)
        return hi;
    return 0.5 * (hi + *std::max_element (v.begin(), v.begin() + static_cast<long> (m)));
}

// skimage.restoration.denoise_tv_chambolle (2-D, eps 2e-4, at most 200 iterations), float64.
std::vector<double> tvChambolle (const std::vector<double>& image, int rows, int cols, double weight)
{
    const size_t size = image.size();
    std::vector<double> p0 (size, 0.0), p1 (size, 0.0), g0 (size, 0.0), g1 (size, 0.0), d (size, 0.0), out = image;
    auto at = [cols] (int r, int c) { return static_cast<size_t> (r) * static_cast<size_t> (cols) + static_cast<size_t> (c); };
    constexpr double tau = 0.25, eps = 2.0e-4;
    double eInit = 0.0, ePrev = 0.0;
    for (int i = 0; i < 200; ++i)
    {
        if (i > 0)
        {
            for (int r = 0; r < rows; ++r)
                for (int c = 0; c < cols; ++c)
                {
                    double v = -p0[at (r, c)] - p1[at (r, c)];
                    if (r > 0)
                        v += p0[at (r - 1, c)];
                    if (c > 0)
                        v += p1[at (r, c - 1)];
                    d[at (r, c)] = v;
                    out[at (r, c)] = image[at (r, c)] + v;
                }
        }
        double e = 0.0;
        for (double v : d)
            e += v * v;
        for (int r = 0; r + 1 < rows; ++r)
            for (int c = 0; c < cols; ++c)
                g0[at (r, c)] = out[at (r + 1, c)] - out[at (r, c)];
        for (int r = 0; r < rows; ++r)
            for (int c = 0; c + 1 < cols; ++c)
                g1[at (r, c)] = out[at (r, c + 1)] - out[at (r, c)];
        for (size_t k = 0; k < size; ++k)
        {
            const double norm = std::sqrt (g0[k] * g0[k] + g1[k] * g1[k]);
            e += weight * norm;
            const double denom = 1.0 + norm * tau / weight;
            p0[k] = (p0[k] - tau * g0[k]) / denom;
            p1[k] = (p1[k] - tau * g1[k]) / denom;
        }
        e /= static_cast<double> (size);
        if (i == 0)
            eInit = ePrev = e;
        else if (std::abs (ePrev - e) < eps * eInit)
            break;
        else
            ePrev = e;
    }
    return out;
}
} // namespace

std::vector<double> denoiseShank (const std::vector<double>& feature, const std::vector<double>& x,
                                  const std::vector<double>& y, double fac)
{
    const size_t nc = feature.size();
    std::vector<int> validIdx;
    for (size_t c = 0; c < nc; ++c)
        if (std::isfinite (feature[c]))
            validIdx.push_back (static_cast<int> (c));
    if (validIdx.size() < 5)
        return feature;
    std::vector<double> xu (x.begin(), x.end()), yu (y.begin(), y.end());
    std::sort (xu.begin(), xu.end());
    xu.erase (std::unique (xu.begin(), xu.end()), xu.end());
    std::sort (yu.begin(), yu.end());
    yu.erase (std::unique (yu.begin(), yu.end()), yu.end());
    const int rows = static_cast<int> (yu.size()), cols = static_cast<int> (xu.size());

    std::vector<Pt> pts;
    for (int c : validIdx)
        pts.push_back ({ x[static_cast<size_t> (c)], y[static_cast<size_t> (c)] });
    const auto tris = delaunay (pts);

    std::vector<double> image (static_cast<size_t> (rows) * static_cast<size_t> (cols), kNan);
    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c)
        {
            const Pt s { xu[static_cast<size_t> (c)], yu[static_cast<size_t> (r)] };
            double value = kNan;
            for (const auto& t : tris)
            {
                const Pt &a = pts[static_cast<size_t> (t[0])], &b = pts[static_cast<size_t> (t[1])], &cc = pts[static_cast<size_t> (t[2])];
                const double o1 = orient (a, b, s), o2 = orient (b, cc, s), o3 = orient (cc, a, s);
                if (o1 < 0.0 || o2 < 0.0 || o3 < 0.0)
                    continue;
                const double area = orient (a, b, cc);
                value = (o2 * feature[static_cast<size_t> (validIdx[static_cast<size_t> (t[0])])]
                         + o3 * feature[static_cast<size_t> (validIdx[static_cast<size_t> (t[1])])]
                         + o1 * feature[static_cast<size_t> (validIdx[static_cast<size_t> (t[2])])]) / area;
                break;
            }
            if (! std::isfinite (value))
            {
                // Outside the hull: nearest valid channel (first in channel order on ties).
                double best = std::numeric_limits<double>::infinity();
                for (size_t k = 0; k < pts.size(); ++k)
                {
                    const double d2 = (pts[k].x - s.x) * (pts[k].x - s.x) + (pts[k].y - s.y) * (pts[k].y - s.y);
                    if (d2 < best)
                    {
                        best = d2;
                        value = feature[static_cast<size_t> (validIdx[k])];
                    }
                }
            }
            image[static_cast<size_t> (r) * static_cast<size_t> (cols) + static_cast<size_t> (c)] = value;
        }

    std::vector<double> absValues (image.size());
    for (size_t k = 0; k < image.size(); ++k)
        absValues[k] = std::abs (image[k]);
    const auto den = tvChambolle (image, rows, cols, medianOf (absValues) * fac);

    std::vector<double> out (nc, kNan);
    for (size_t ch = 0; ch < nc; ++ch)
    {
        const auto ix = std::lower_bound (xu.begin(), xu.end(), x[ch]) - xu.begin();
        const auto iy = std::lower_bound (yu.begin(), yu.end(), y[ch]) - yu.begin();
        out[ch] = den[static_cast<size_t> (iy) * static_cast<size_t> (cols) + static_cast<size_t> (ix)];
    }
    return out;
}

std::vector<double> aggregateSnippetFeatures (const std::vector<SnippetFeatureSet>& snippets, int nc, int nf,
                                              const std::vector<char>& isLogAmplitude,
                                              const std::vector<double>& x, const std::vector<double>& y,
                                              double fac, std::vector<int>* labelsOut)
{
    std::vector<double> agg (static_cast<size_t> (nc) * static_cast<size_t> (nf), kNan);
    std::vector<int> labels (static_cast<size_t> (nc), 0);
    std::vector<double> vals;
    for (int c = 0; c < nc; ++c)
    {
        for (int f = 0; f < nf; ++f)
        {
            const bool logAmp = ! isLogAmplitude.empty() && isLogAmplitude[static_cast<size_t> (f)];
            vals.clear();
            for (const auto& s : snippets)
            {
                const double v = s.values[static_cast<size_t> (c) * static_cast<size_t> (nf) + static_cast<size_t> (f)];
                if (std::isfinite (v))
                    vals.push_back (logAmp ? std::pow (10.0, v / 20.0) : v);
            }
            if (! vals.empty())
            {
                const double m = medianOf (vals);
                agg[static_cast<size_t> (c) * static_cast<size_t> (nf) + static_cast<size_t> (f)] = logAmp ? 20.0 * std::log10 (m) : m;
            }
        }
        // Mode of the labels, the smallest on ties (pandas Series.mode().iloc[0]).
        std::map<int, int> counts;
        for (const auto& s : snippets)
            ++counts[s.labels.empty() ? 0 : s.labels[static_cast<size_t> (c)]];
        int best = 0, bestCount = -1;
        for (const auto& [label, count] : counts)
            if (count > bestCount)
            {
                best = label;
                bestCount = count;
            }
        labels[static_cast<size_t> (c)] = best;
    }
    std::vector<double> column (static_cast<size_t> (nc));
    for (int f = 0; f < nf; ++f)
    {
        for (int c = 0; c < nc; ++c)
            column[static_cast<size_t> (c)] = labels[static_cast<size_t> (c)] != 0 ? kNan : agg[static_cast<size_t> (c) * static_cast<size_t> (nf) + static_cast<size_t> (f)];
        auto den = denoiseShank (column, x, y, fac);
        std::vector<double> finite;
        for (double v : den)
            if (std::isfinite (v))
                finite.push_back (v);
        const double fill = finite.empty() ? kNan : medianOf (finite);
        for (int c = 0; c < nc; ++c)
        {
            const double v = den[static_cast<size_t> (c)];
            agg[static_cast<size_t> (c) * static_cast<size_t> (nf) + static_cast<size_t> (f)] = std::isfinite (v) ? v : fill;
        }
    }
    if (labelsOut)
        *labelsOut = labels;
    return agg;
}
} // namespace ProbeViewer::Localization::Realtime
