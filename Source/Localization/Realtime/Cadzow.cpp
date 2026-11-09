#include "Cadzow.hpp"

#include "PyDsp.hpp"
#include "Threads.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <numeric>
#include <thread>

namespace ProbeViewer::Localization::Realtime
{
namespace
{
constexpr double kPi = 3.14159265358979323846;

// Singular values (descending) and the top-r reconstruction of a complex matrix, by one-sided
// (Hestenes) Jacobi on the columns. A: rows x cols, column-major; m = rows >= n = cols.
struct ColumnSvd
{
    int m = 0, n = 0;
    std::vector<cplx> a;  // A V, columns sigma_j u_j (column-major)
    std::vector<cplx> v;  // V, n x n column-major
    std::vector<double> sigma;
    std::vector<int> order; // columns by decreasing sigma

    void compute (const std::vector<cplx>& A, int rows, int cols)
    {
        m = rows;
        n = cols;
        a = A;
        v.assign (static_cast<size_t> (n * n), cplx (0.0, 0.0));
        for (int j = 0; j < n; ++j)
            v[static_cast<size_t> (j * n + j)] = 1.0;
        for (int sweep = 0; sweep < 40; ++sweep)
        {
            bool rotated = false;
            for (int p = 0; p < n - 1; ++p)
                for (int q = p + 1; q < n; ++q)
                {
                    cplx* ap = &a[static_cast<size_t> (p * m)];
                    cplx* aq = &a[static_cast<size_t> (q * m)];
                    double alpha = 0.0, beta = 0.0;
                    cplx gamma = 0.0;
                    for (int i = 0; i < m; ++i)
                    {
                        alpha += std::norm (ap[i]);
                        beta += std::norm (aq[i]);
                        gamma += std::conj (ap[i]) * aq[i];
                    }
                    const double g = std::abs (gamma);
                    if (g <= 1e-15 * std::sqrt (alpha * beta) || g < 1e-300)
                        continue;
                    rotated = true;
                    const double zeta = (beta - alpha) / (2.0 * g);
                    const double t = (zeta >= 0.0 ? 1.0 : -1.0) / (std::abs (zeta) + std::sqrt (1.0 + zeta * zeta));
                    const double c = 1.0 / std::sqrt (1.0 + t * t), s = c * t;
                    const cplx ph = std::conj (gamma / g);
                    for (int i = 0; i < m; ++i)
                    {
                        const cplx x = ap[i], y = aq[i] * ph;
                        ap[i] = c * x - s * y;
                        aq[i] = s * x + c * y;
                    }
                    cplx* vp = &v[static_cast<size_t> (p * n)];
                    cplx* vq = &v[static_cast<size_t> (q * n)];
                    for (int i = 0; i < n; ++i)
                    {
                        const cplx x = vp[i], y = vq[i] * ph;
                        vp[i] = c * x - s * y;
                        vq[i] = s * x + c * y;
                    }
                }
            if (! rotated)
                break;
        }
        sigma.assign (static_cast<size_t> (n), 0.0);
        for (int j = 0; j < n; ++j)
        {
            double s2 = 0.0;
            for (int i = 0; i < m; ++i)
                s2 += std::norm (a[static_cast<size_t> (j * m + i)]);
            sigma[static_cast<size_t> (j)] = std::sqrt (s2);
        }
        order.resize (static_cast<size_t> (n));
        std::iota (order.begin(), order.end(), 0);
        std::stable_sort (order.begin(), order.end(), [this] (int i, int j) { return sigma[static_cast<size_t> (i)] > sigma[static_cast<size_t> (j)]; });
    }

    // Rank kept by ibldsp.cadzow._apply_rank_threshold.
    int rank (int r, double gapThreshold) const
    {
        if (gapThreshold <= 0.0)
            return std::min (r, n);
        double best = -1.0;
        int argmax = 0;
        for (int i = 0; i + 1 < n; ++i)
        {
            const double ratio = sigma[static_cast<size_t> (order[static_cast<size_t> (i)])] / (sigma[static_cast<size_t> (order[static_cast<size_t> (i + 1)])] + 1e-10);
            if (ratio > best)
            {
                best = ratio;
                argmax = i;
            }
        }
        const int adapt = best >= gapThreshold ? argmax + 1 : r;
        return std::clamp (adapt, 1, r);
    }

    // Multipliers of the singular values (by decreasing order) under ibldsp.cadzow._apply_shrinkage:
    // Gavish-Donoho optimal shrinkage with the noise level median(s) * scale / sqrt(MP median),
    // capped at rank r. shapeRows x shapeCols is the trajectory matrix's shape.
    std::vector<double> shrinkageFactors (int r, double scale, int shapeRows, int shapeCols) const
    {
        constexpr double tiny = 2.2250738585072014e-308;
        const double beta = static_cast<double> (std::min (shapeRows, shapeCols)) / static_cast<double> (std::max (shapeRows, shapeCols));
        std::vector<double> s (static_cast<size_t> (n));
        for (int t = 0; t < n; ++t)
            s[static_cast<size_t> (t)] = sigma[static_cast<size_t> (order[static_cast<size_t> (t)])];
        std::vector<double> sorted = s;
        std::sort (sorted.begin(), sorted.end());
        const size_t mid = sorted.size() / 2;
        const double med = sorted.size() % 2 ? sorted[mid] : 0.5 * (sorted[mid - 1] + sorted[mid]);
        const double noise = std::max (scale * med / std::sqrt (marchenkoPasturMedian (beta)), tiny);
        std::vector<double> factor (static_cast<size_t> (n), 0.0);
        for (int t = 0; t < n && t < r; ++t)
        {
            const double st = s[static_cast<size_t> (t)];
            const double y = st / noise;
            if (! (y > 1.0 + std::sqrt (beta)) || ! (st > 0.0))
                continue;
            const double q = (y * y - beta - 1.0) * (y * y - beta - 1.0) - 4.0 * beta;
            const double eta = std::sqrt (std::max (q, 0.0)) / std::max (y, tiny);
            factor[static_cast<size_t> (t)] = eta * noise / st;
        }
        return factor;
    }

    // Multipliers of a hard rank-k cut.
    std::vector<double> rankFactors (int k) const
    {
        std::vector<double> factor (static_cast<size_t> (n), 0.0);
        for (int t = 0; t < k && t < n; ++t)
            factor[static_cast<size_t> (t)] = 1.0;
        return factor;
    }

    // Element (row, col) of sum_t factor_t a_t v_t^H over the singular triplets by decreasing order.
    cplx at (int row, int col, const std::vector<double>& factor) const
    {
        cplx s = 0.0;
        for (int t = 0; t < n; ++t)
        {
            const double f = factor[static_cast<size_t> (t)];
            if (f == 0.0)
                continue;
            const int j = order[static_cast<size_t> (t)];
            s += f * a[static_cast<size_t> (j * m + row)] * std::conj (v[static_cast<size_t> (j * n + col)]);
        }
        return s;
    }
};

// ibldsp.cadzow.traj_matrix_indices(n)[i][j] = i + (ncols - 1 - j).
struct TrajDims
{
    int nrows, ncols;
    explicit TrajDims (int n) : nrows (n / 2 + 1), ncols ((n + 1) / 2) {}
    int at (int i, int j) const { return i + (ncols - 1 - j); }
};

struct Window
{
    int first = 0, count = 0;
    std::vector<double> gain;
    int rows = 0, cols = 0;
    std::vector<int> entryRow, entryCol, entryChannel; // trajectory entries (row-major order)
    std::vector<int> trcount;                          // entries per channel of the window
};

std::vector<double> sortedUnique (std::vector<double> v)
{
    std::sort (v.begin(), v.end());
    v.erase (std::unique (v.begin(), v.end()), v.end());
    return v;
}

void buildTrajectory (Window& w, const std::vector<double>& x, const std::vector<double>& y)
{
    std::vector<double> xs (x.begin() + w.first, x.begin() + w.first + w.count);
    std::vector<double> ys (y.begin() + w.first, y.begin() + w.first + w.count);
    const auto xu = sortedUnique (xs), yu = sortedUnique (ys);
    std::vector<int> ix (static_cast<size_t> (w.count)), iy (static_cast<size_t> (w.count));
    for (int c = 0; c < w.count; ++c)
    {
        ix[static_cast<size_t> (c)] = static_cast<int> (std::lower_bound (xu.begin(), xu.end(), xs[static_cast<size_t> (c)]) - xu.begin());
        iy[static_cast<size_t> (c)] = static_cast<int> (std::lower_bound (yu.begin(), yu.end(), ys[static_cast<size_t> (c)]) - yu.begin());
    }
    const TrajDims tx (static_cast<int> (xu.size())), ty (static_cast<int> (yu.size()));
    w.rows = tx.nrows * ty.nrows;
    w.cols = tx.ncols * ty.ncols;
    // channel of each (ix, iy) pair: first channel in window order (ismember2d)
    std::vector<int> lookup (xu.size() * yu.size(), -1);
    for (int c = w.count - 1; c >= 0; --c)
        lookup[static_cast<size_t> (ix[static_cast<size_t> (c)]) * yu.size() + static_cast<size_t> (iy[static_cast<size_t> (c)])] = c;
    w.trcount.assign (static_cast<size_t> (w.count), 0);
    for (int R = 0; R < w.rows; ++R)
        for (int C = 0; C < w.cols; ++C)
        {
            const int bi = R / ty.nrows, ri = R % ty.nrows;
            const int bj = C / ty.ncols, cj = C % ty.ncols;
            const int xi = tx.at (bi, bj), yi = ty.at (ri, cj);
            if (xi < 0 || yi < 0 || xi >= static_cast<int> (xu.size()) || yi >= static_cast<int> (yu.size()))
                continue;
            const int ch = lookup[static_cast<size_t> (xi) * yu.size() + static_cast<size_t> (yi)];
            if (ch < 0)
                continue;
            w.entryRow.push_back (R);
            w.entryCol.push_back (C);
            w.entryChannel.push_back (ch);
            ++w.trcount[static_cast<size_t> (ch)];
        }
}

double median (std::vector<double> v)
{
    const size_t n = v.size(), m = n / 2;
    std::nth_element (v.begin(), v.begin() + static_cast<long> (m), v.end());
    const double hi = v[m];
    if (n % 2 == 1)
        return hi;
    return 0.5 * (hi + *std::max_element (v.begin(), v.begin() + static_cast<long> (m)));
}

// One frequency bin of one window: rank reduction (+ outlier correction), averaged back to the
// window's channels (ibldsp.cadzow._process_window for a single bin, niter iterations).
void processBin (const Window& w, std::vector<cplx>& vals /* [count] in/out */, const CadzowParams& p)
{
    const int rows = w.rows, cols = w.cols;
    const bool transpose = rows < cols; // one-sided Jacobi needs rows >= cols
    const int m = transpose ? cols : rows, n = transpose ? rows : cols;
    std::vector<cplx> A (static_cast<size_t> (m * n));
    ColumnSvd svd;
    auto fill = [&] (const std::vector<cplx>& channelValues) {
        std::fill (A.begin(), A.end(), cplx (0.0, 0.0));
        for (size_t e = 0; e < w.entryRow.size(); ++e)
        {
            const int R = w.entryRow[e], C = w.entryCol[e];
            const cplx value = channelValues[static_cast<size_t> (w.entryChannel[e])];
            if (! transpose)
                A[static_cast<size_t> (C * m + R)] = value;
            else
                A[static_cast<size_t> (R * m + C)] = std::conj (value);
        }
    };
    auto reconstruct = [&] (const std::vector<double>& factor, std::vector<cplx>& out) {
        out.assign (static_cast<size_t> (w.count), cplx (0.0, 0.0));
        for (size_t e = 0; e < w.entryRow.size(); ++e)
        {
            const int R = w.entryRow[e], C = w.entryCol[e];
            const cplx value = transpose ? std::conj (svd.at (C, R, factor)) : svd.at (R, C, factor);
            const int ch = w.entryChannel[e];
            out[static_cast<size_t> (ch)] += value / static_cast<double> (w.trcount[static_cast<size_t> (ch)]);
        }
    };
    // Shrinkage, or the hard (optionally adaptive) rank cut.
    auto threshold = [&] () {
        return p.shrinkage > 0.0 ? svd.shrinkageFactors (p.rank, p.shrinkage, rows, cols)
                                 : svd.rankFactors (svd.rank (p.rank, p.gapThreshold));
    };
    for (int iter = 0; iter < p.niter; ++iter)
    {
        fill (vals);
        svd.compute (A, m, n);
        auto k = threshold();
        std::vector<cplx> rec;
        if (p.ppcaK > 0.0)
        {
            reconstruct (k, rec);
            std::vector<double> residual (static_cast<size_t> (w.count));
            for (int c = 0; c < w.count; ++c)
                residual[static_cast<size_t> (c)] = std::abs (vals[static_cast<size_t> (c)] - rec[static_cast<size_t> (c)]);
            const double med = median (residual);
            std::vector<double> dev (residual.size());
            for (size_t c = 0; c < residual.size(); ++c)
                dev[c] = std::abs (residual[c] - med);
            const double mad = median (dev);
            std::vector<cplx> clean = vals;
            for (int c = 0; c < w.count; ++c)
                if (residual[static_cast<size_t> (c)] > med + p.ppcaK * mad)
                    clean[static_cast<size_t> (c)] = rec[static_cast<size_t> (c)];
            fill (clean);
            svd.compute (A, m, n);
            k = threshold();
        }
        reconstruct (k, rec);
        vals = rec;
    }
}
} // namespace

double marchenkoPasturMedian (double beta)
{
    // CDF by the substitution t = lo + (hi - lo) (1 - cos phi) / 2, which removes the square-root
    // edges: dCDF/dphi = ((hi - lo) / 2)^2 sin^2(phi) / (2 pi beta t). Simpson on [0, theta],
    // then bisection for CDF = 1/2 (scipy brentq on quad in ibldsp; same value to ~1e-10).
    const double lo = (1.0 - std::sqrt (beta)) * (1.0 - std::sqrt (beta));
    const double hi = (1.0 + std::sqrt (beta)) * (1.0 + std::sqrt (beta));
    const double half = 0.5 * (hi - lo);
    auto density = [&] (double phi) {
        const double t = lo + half * (1.0 - std::cos (phi));
        const double s = std::sin (phi);
        return half * half * s * s / (2.0 * kPi * beta * t);
    };
    auto cdf = [&] (double theta) {
        constexpr int n = 4000; // even
        const double h = theta / n;
        double sum = density (0.0) + density (theta);
        for (int i = 1; i < n; ++i)
            sum += (i % 2 ? 4.0 : 2.0) * density (i * h);
        return sum * h / 3.0;
    };
    double a = 0.0, b = kPi;
    for (int it = 0; it < 80; ++it)
    {
        const double c = 0.5 * (a + b);
        (cdf (c) < 0.5 ? a : b) = c;
    }
    return lo + half * (1.0 - std::cos (0.5 * (a + b)));
}

std::vector<double> cadzowDenoise (const std::vector<double>& wav, int nc, int ns,
                                   const std::vector<double>& xIn, const std::vector<double>& yIn,
                                   double fs, const CadzowParams& p)
{
    if (p.fillGrid)
    {
        // ibldsp.cadzow._fill_grid: the channels on the grid of unique depth x lateral positions
        // (depth-major), each empty position the mean of its known 4-neighbours, outwards until
        // all are known; then the denoiser on the grid with windows of the same depth extent.
        const auto xu = sortedUnique (std::vector<double> (xIn.begin(), xIn.begin() + nc));
        const auto yu = sortedUnique (std::vector<double> (yIn.begin(), yIn.begin() + nc));
        const int nx = static_cast<int> (xu.size()), ny = static_cast<int> (yu.size()), ng = nx * ny;
        std::vector<int> ireal (static_cast<size_t> (nc));
        std::vector<double> g (static_cast<size_t> (ng) * static_cast<size_t> (ns), 0.0);
        std::vector<char> known (static_cast<size_t> (ng), 0);
        for (int c = 0; c < nc; ++c)
        {
            const int ix = static_cast<int> (std::lower_bound (xu.begin(), xu.end(), xIn[static_cast<size_t> (c)]) - xu.begin());
            const int iy = static_cast<int> (std::lower_bound (yu.begin(), yu.end(), yIn[static_cast<size_t> (c)]) - yu.begin());
            const int k = iy * nx + ix;
            ireal[static_cast<size_t> (c)] = k;
            std::copy_n (&wav[static_cast<size_t> (c) * static_cast<size_t> (ns)], ns, &g[static_cast<size_t> (k) * static_cast<size_t> (ns)]);
            known[static_cast<size_t> (k)] = 1;
        }
        for (bool done = false; ! done;)
        {
            done = true;
            std::vector<char> next = known;
            for (int iy = 0; iy < ny; ++iy)
                for (int ix = 0; ix < nx; ++ix)
                {
                    const int k = iy * nx + ix;
                    if (known[static_cast<size_t> (k)])
                        continue;
                    int neighbours[4], count = 0;
                    if (iy > 0 && known[static_cast<size_t> (k - nx)])
                        neighbours[count++] = k - nx;
                    if (iy + 1 < ny && known[static_cast<size_t> (k + nx)])
                        neighbours[count++] = k + nx;
                    if (ix > 0 && known[static_cast<size_t> (k - 1)])
                        neighbours[count++] = k - 1;
                    if (ix + 1 < nx && known[static_cast<size_t> (k + 1)])
                        neighbours[count++] = k + 1;
                    if (count == 0)
                    {
                        done = false;
                        continue;
                    }
                    double* out = &g[static_cast<size_t> (k) * static_cast<size_t> (ns)];
                    for (int t = 0; t < ns; ++t)
                    {
                        double s = 0.0;
                        for (int q = 0; q < count; ++q)
                            s += g[static_cast<size_t> (neighbours[q]) * static_cast<size_t> (ns) + static_cast<size_t> (t)];
                        out[t] = s / count;
                    }
                    next[static_cast<size_t> (k)] = 1;
                }
            known.swap (next);
        }
        std::vector<double> gx (static_cast<size_t> (ng)), gy (static_cast<size_t> (ng));
        for (int iy = 0; iy < ny; ++iy)
            for (int ix = 0; ix < nx; ++ix)
            {
                gx[static_cast<size_t> (iy * nx + ix)] = xu[static_cast<size_t> (ix)];
                gy[static_cast<size_t> (iy * nx + ix)] = yu[static_cast<size_t> (iy)];
            }
        const int f = ng / nc;
        CadzowParams q = p;
        q.fillGrid = false;
        q.nswx = p.nswx * f;
        q.ovx = p.ovx * f;
        const auto den = cadzowDenoise (g, ng, ns, gx, gy, fs, q);
        std::vector<double> out (static_cast<size_t> (nc) * static_cast<size_t> (ns));
        for (int c = 0; c < nc; ++c)
            std::copy_n (&den[static_cast<size_t> (ireal[static_cast<size_t> (c)]) * static_cast<size_t> (ns)], ns,
                         &out[static_cast<size_t> (c) * static_cast<size_t> (ns)]);
        return out;
    }

    const Fft fft (ns);
    const int nf = ns / 2 + 1;
    // rfft of every channel, plus the phantom row of the npad = 0 path (a copy of channel nc - 2).
    std::vector<std::vector<cplx>> W (static_cast<size_t> (nc + 1));
    for (int c = 0; c < nc; ++c)
        W[static_cast<size_t> (c)] = fft.rfft (&wav[static_cast<size_t> (c) * static_cast<size_t> (ns)]);
    W[static_cast<size_t> (nc)] = W[static_cast<size_t> (nc - 2)];
    std::vector<double> x (xIn.begin(), xIn.begin() + nc), y (yIn.begin(), yIn.begin() + nc);
    x.push_back (xIn[static_cast<size_t> (nc - 2)]);
    y.push_back (yIn[static_cast<size_t> (nc - 2)] + 120.0);
    int imax = p.fmax > 0.0 ? 0 : nf;
    while (imax < nf && imax * fs / ns < p.fmax)
        ++imax;

    // Windows and their partition-of-unity gains (ovx <= nswx / 2).
    const int step = p.nswx - p.ovx;
    const int nwinx = static_cast<int> (std::ceil (static_cast<double> (nc - p.ovx) / step));
    std::vector<double> hanning (static_cast<size_t> (p.ovx));
    const int M = 2 * p.ovx - 1;
    for (int i = 0; i < p.ovx; ++i)
        hanning[static_cast<size_t> (i)] = 0.5 - 0.5 * std::cos (2.0 * kPi * i / (M - 1));
    std::vector<Window> windows;
    for (int wi = 0; wi < nwinx; ++wi)
    {
        Window w;
        w.first = wi * step;
        const int last = w.first + p.nswx;
        w.count = std::min (last, nc + 1) - w.first;
        std::vector<double> gain;
        if (w.first == 0)
        {
            gain.assign (static_cast<size_t> (p.nswx - p.ovx), 1.0);
            gain.insert (gain.end(), hanning.rbegin(), hanning.rend());
        }
        else if (last >= nc)
        {
            gain.assign (hanning.begin(), hanning.end());
            gain.insert (gain.end(), static_cast<size_t> (p.nswx - p.ovx), 1.0);
        }
        else
        {
            gain.assign (hanning.begin(), hanning.end());
            gain.insert (gain.end(), static_cast<size_t> (p.nswx - 2 * p.ovx), 1.0);
            gain.insert (gain.end(), hanning.rbegin(), hanning.rend());
        }
        gain.resize (static_cast<size_t> (w.count));
        w.gain = gain;
        buildTrajectory (w, x, y);
        windows.push_back (std::move (w));
    }

    // Rank reduction of every (window, bin < imax), in parallel.
    std::vector<std::vector<std::vector<cplx>>> out (windows.size()); // [window][bin][channel]
    for (size_t wi = 0; wi < windows.size(); ++wi)
        out[wi].resize (static_cast<size_t> (imax));
    std::atomic<int> next (0);
    const int total = static_cast<int> (windows.size()) * imax;
    auto worker = [&] () {
        for (int task = next++; task < total; task = next++)
        {
            const size_t wi = static_cast<size_t> (task / imax);
            const int k = task % imax;
            const Window& w = windows[wi];
            std::vector<cplx> vals (static_cast<size_t> (w.count));
            for (int c = 0; c < w.count; ++c)
                vals[static_cast<size_t> (c)] = W[static_cast<size_t> (w.first + c)][static_cast<size_t> (k)];
            processBin (w, vals, p);
            out[wi][static_cast<size_t> (k)] = std::move (vals);
        }
    };
    const int nThreads = std::max (1, p.threads);
    std::vector<std::thread> pool;
    for (int t = 1; t < nThreads; ++t)
        pool.emplace_back ([&worker] {
            lowerThreadPriority();
            worker();
        });
    worker();
    for (auto& t : pool)
        t.join();

    // Gain-weighted accumulation (bins >= imax pass through), phantom row dropped, inverse FFT.
    std::vector<std::vector<cplx>> acc (static_cast<size_t> (nc + 1), std::vector<cplx> (static_cast<size_t> (nf), cplx (0.0, 0.0)));
    for (size_t wi = 0; wi < windows.size(); ++wi)
    {
        const Window& w = windows[wi];
        for (int c = 0; c < w.count; ++c)
        {
            auto& row = acc[static_cast<size_t> (w.first + c)];
            const double g = w.gain[static_cast<size_t> (c)];
            for (int k = 0; k < nf; ++k)
            {
                const cplx v = k < imax ? out[wi][static_cast<size_t> (k)][static_cast<size_t> (c)] : W[static_cast<size_t> (w.first + c)][static_cast<size_t> (k)];
                row[static_cast<size_t> (k)] += v * g;
            }
        }
    }
    std::vector<double> result (static_cast<size_t> (nc) * static_cast<size_t> (ns));
    for (int c = 0; c < nc; ++c)
    {
        const auto sig = fft.irfft (acc[static_cast<size_t> (c)]);
        for (int t = 0; t < ns; ++t)
            result[static_cast<size_t> (c) * static_cast<size_t> (ns) + static_cast<size_t> (t)] = static_cast<double> (static_cast<float> (sig[static_cast<size_t> (t)]));
    }
    return result;
}
} // namespace ProbeViewer::Localization::Realtime
