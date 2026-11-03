#include "EphysOnlyLocalizer.hpp"

#include "MiniJson.hpp"
#include "Npy.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>

namespace ProbeViewer::Localization::Realtime
{
namespace
{
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kPi = 3.14159265358979323846;

double deg2rad (double d) { return d * kPi / 180.0; }
double rad2deg (double r) { return r * 180.0 / kPi; }

V3 sub (const V3& a, const V3& b) { return { a[0] - b[0], a[1] - b[1], a[2] - b[2] }; }
V3 add (const V3& a, const V3& b) { return { a[0] + b[0], a[1] + b[1], a[2] + b[2] }; }
V3 mul (const V3& a, double s) { return { a[0] * s, a[1] * s, a[2] * s }; }
double dot (const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
double norm (const V3& a) { return std::sqrt (dot (a, a)); }
V3 cross (const V3& a, const V3& b)
{
    return { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
}
V3 unitVector (const V3& v)
{
    const double n = norm (v);
    if (! std::isfinite (n) || n < 1e-12)
        return { 0.0, 0.0, -1.0 };
    return mul (v, 1.0 / n);
}
bool finite3 (const V3& v) { return std::isfinite (v[0]) && std::isfinite (v[1]) && std::isfinite (v[2]); }
// Finite and not the all-zero placeholder (ephys_only.valid_positions).
bool validPosition (const V3& v) { return finite3 (v) && ! (v[0] == 0.0 && v[1] == 0.0 && v[2] == 0.0); }
V3 mirrorLeft (V3 v)
{
    v[0] = -std::abs (v[0]);
    return v;
}

// Eigenvector of the largest eigenvalue of a symmetric 3x3 matrix (cyclic Jacobi).
V3 principalEigenvector (std::array<std::array<double, 3>, 3> a)
{
    std::array<std::array<double, 3>, 3> v { { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } } };
    for (int sweep = 0; sweep < 50; ++sweep)
    {
        const double off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        if (off < 1e-30)
            break;
        for (int p = 0; p < 2; ++p)
            for (int q = p + 1; q < 3; ++q)
            {
                if (std::abs (a[p][q]) < 1e-300)
                    continue;
                const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                const double t = (theta >= 0 ? 1.0 : -1.0) / (std::abs (theta) + std::sqrt (theta * theta + 1.0));
                const double c = 1.0 / std::sqrt (t * t + 1.0), s = t * c;
                for (int k = 0; k < 3; ++k)
                {
                    const double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; ++k)
                {
                    const double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for (int k = 0; k < 3; ++k)
                {
                    const double vkp = v[k][p], vkq = v[k][q];
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
    }
    int best = 0;
    for (int k = 1; k < 3; ++k)
        if (a[k][k] > a[best][best])
            best = k;
    return { v[0][best], v[1][best], v[2][best] };
}

// First principal direction of points (rows), oriented from the first to the last point.
V3 principalAxis (const std::vector<V3>& pts)
{
    V3 c { 0, 0, 0 };
    for (const auto& p : pts)
        c = add (c, p);
    c = mul (c, 1.0 / static_cast<double> (pts.size()));
    std::array<std::array<double, 3>, 3> cov {};
    for (const auto& p : pts)
    {
        const V3 d = sub (p, c);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                cov[i][j] += d[i] * d[j];
    }
    V3 axis = principalEigenvector (cov);
    if (dot (axis, sub (pts.back(), pts.front())) < 0)
        axis = mul (axis, -1.0);
    return axis;
}

double rhoT (double x, double nu, double loc, double s)
{
    const double z = (x - loc) / s;
    return 0.5 * (nu + 1.0) * std::log1p (z * z / nu);
}

double rhoT2 (double r, double nu, double s)
{
    return 0.5 * (nu + 2.0) * std::log1p (r * r / (nu * s * s));
}

double expSchedule (double start, double end, double alpha)
{
    alpha = std::clamp (alpha, 0.0, 1.0);
    if (start <= 0.0 || end <= 0.0)
        return start * (1.0 - alpha) + end * alpha;
    return start * std::pow (end / start, alpha);
}

std::array<double, 2> clipBend (std::array<double, 2> bend, double maxAngle)
{
    maxAngle = std::max (0.0, maxAngle);
    const double n = std::hypot (bend[0], bend[1]);
    if (n > maxAngle && n > 0.0)
    {
        bend[0] *= maxAngle / n;
        bend[1] *= maxAngle / n;
    }
    return bend;
}

V3 probeAxisFromAngles (double thetaDeg, double phiDeg, bool zNegative)
{
    const double th = deg2rad (thetaDeg), ph = deg2rad (phiDeg);
    const V3 u { std::sin (th) * std::cos (ph), std::sin (th) * std::sin (ph), zNegative ? -std::cos (th) : std::cos (th) };
    return mul (u, 1.0 / norm (u));
}

V3 bendDirection (const V3& direction, const std::array<double, 2>& bendUv, double maxAngle)
{
    const V3 d = unitVector (direction);
    const auto bend = clipBend (bendUv, maxAngle);
    const double angle = std::hypot (bend[0], bend[1]);
    if (angle <= 1e-12)
        return d;
    V3 reference { 0.0, 0.0, 1.0 };
    if (std::abs (dot (d, reference)) > 0.90)
        reference = { 1.0, 0.0, 0.0 };
    const V3 u = unitVector (cross (d, reference));
    const V3 v = unitVector (cross (u, d));
    const V3 tangent = unitVector (add (mul (u, bend[0]), mul (v, bend[1])));
    const double rad = deg2rad (angle);
    return unitVector (add (mul (d, std::cos (rad)), mul (tangent, std::sin (rad))));
}

double median (std::vector<double> v)
{
    if (v.empty())
        return kNaN;
    const size_t n = v.size(), m = n / 2;
    std::nth_element (v.begin(), v.begin() + static_cast<long> (m), v.end());
    const double hi = v[m];
    if (n % 2 == 1)
        return hi;
    return 0.5 * (hi + *std::max_element (v.begin(), v.begin() + static_cast<long> (m)));
}

std::string readText (const std::string& path)
{
    std::ifstream f (path, std::ios::binary);
    if (! f)
        throw std::runtime_error ("cannot open " + path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

std::vector<int> sectionIndices (const std::vector<double>& centers, int k, int n)
{
    const double lo = *std::min_element (centers.begin(), centers.end());
    const double hi = *std::max_element (centers.begin(), centers.end());
    const double span = std::max (hi - lo, 1e-6);
    std::vector<int> out;
    for (size_t i = 0; i < centers.size(); ++i)
    {
        const double s = (centers[i] - lo) / span;
        if (s >= static_cast<double> (k) / n - 1e-9 && s <= static_cast<double> (k + 1) / n + 1e-9)
            out.push_back (static_cast<int> (i));
    }
    return out;
}
} // namespace

// =============================================================================================
// Bundle
// =============================================================================================

bool LocalizerBundle::load (const std::string& directory, std::string& error)
{
    try
    {
        const std::string dir = directory + "/";
        const Json m = Json::parse (readText (dir + "manifest.json"));
        F = static_cast<int> (m["n_features"].num());
        featureNames.clear();
        for (const auto& name : m["feature_names"].items)
            featureNames.push_back (name.str());
        featureSet = m.has ("feature_set") ? m["feature_set"].str() : "";
        modelCommit = m.has ("model_commit") ? m["model_commit"].str() : "";

        const auto& p = m["params"];
        nu = p["nu"].num();
        kappa = p["kappa"].num();
        T = p["T"].num();
        widen = p.has ("plan_prior_widen") ? p["plan_prior_widen"].num() : 1.0;
        const auto& pr = p["prior"];
        const char* axes[3] = { "x", "y", "z" };
        for (int k = 0; k < 3; ++k)
            for (int j = 0; j < 3; ++j)
                prior.centre[k][j] = pr["centre"][axes[k]][j].num();
        prior.axisNu = pr["axis"][0].num();
        prior.axisS = pr["axis"][1].num();
        prior.lengthMedian = pr["length"][0].num();
        prior.lengthSd = pr["length"][1].num();
        prior.bendNu = pr["bend"][0].num();
        prior.bendS = pr["bend"][1].num();

        const auto& sc = m["score_config"];
        score.gridUm = sc["grid_um"].num();
        score.channelSpacingUm = sc["channel_spacing_um"].num();
        score.depthBinUm = sc["depth_bin_um"].num();
        score.minChannelsPerBin = static_cast<int> (sc["min_channels_per_bin"].num());
        score.lambdaOut = sc["lambda_out"].num();
        score.outFreeUm = sc["out_free_um"].num();
        score.outScaleUm = sc["out_scale_um"].num();
        score.planMinSpanUm = sc["plan_min_span_um"].num();
        score.degeneratePriorWiden = sc["degenerate_prior_widen"].num();
        score.nSections = static_cast<int> (sc["n_sections"].num());
        score.minUsablePlanSpanUm = sc["min_usable_plan_span_um"].num();
        score.scalePoolSize = static_cast<int> (sc["scale_pool_size"].num());
        score.scaleK = sc["scale_K"].num();

        const auto& ac = m["anneal_config"];
        anneal.nIterations = static_cast<int> (ac["n_iterations"].num());
        anneal.initialTemperature = ac["initial_temperature"].num();
        anneal.finalTemperature = ac["final_temperature"].num();
        anneal.temperatureDecayIters = static_cast<int> (ac["temperature_decay_iters"].num());
        anneal.initXyzSigmaUm = ac["init_xyz_sigma_um"].num();
        anneal.finalXyzSigmaUm = ac["final_xyz_sigma_um"].num();
        anneal.initAngleSigmaDeg = ac["init_angle_sigma_deg"].num();
        anneal.finalAngleSigmaDeg = ac["final_angle_sigma_deg"].num();
        anneal.initAzimuthSigmaDeg = ac["init_azimuth_sigma_deg"].num();
        anneal.finalAzimuthSigmaDeg = ac["final_azimuth_sigma_deg"].num();
        anneal.initLengthSigmaUm = ac["init_length_sigma_um"].num();
        anneal.finalLengthSigmaUm = ac["final_length_sigma_um"].num();
        anneal.initJointAngleSigmaDeg = ac["init_joint_angle_sigma_deg"].num();
        anneal.finalJointAngleSigmaDeg = ac["final_joint_angle_sigma_deg"].num();
        for (int k = 0; k < 4; ++k)
            anneal.moveProbabilities[static_cast<size_t> (k)] = ac["move_probabilities"][static_cast<size_t> (k)].num();
        anneal.nSections = static_cast<int> (ac["n_sections"].num());
        anneal.maxJointAngleDeg = ac["max_joint_angle_deg"].num();
        anneal.minAngleDeg = ac["min_angle_deg"].num();
        anneal.maxAngleDeg = ac["max_angle_deg"].num();
        anneal.probeLengthUm = ac["probe_length_um"].num();
        anneal.minProbeLengthUm = ac["min_probe_length_um"].num();
        anneal.maxProbeLengthUm = ac["max_probe_length_um"].num();
        anneal.reheatPatience = static_cast<int> (ac["reheat_patience"].num());
        anneal.maxReheats = static_cast<int> (ac["max_reheats"].num());
        anneal.reheatTemperatureFactor = ac["reheat_temperature_factor"].num();
        anneal.reheatMinTemperature = ac["reheat_min_temperature"].num();
        anneal.reheatSigmaFactor = ac["reheat_sigma_factor"].num();
        anneal.reheatDecay = ac["reheat_decay"].num();
        anneal.earlyStopPatience = static_cast<int> (ac["early_stop_patience"].num());
        nRestarts = static_cast<int> (m["n_restarts"].num());
        seed = static_cast<long long> (m["seed"].num());

        auto vec = [&dir, &m] (const std::string& key) {
            return loadNpy (dir + m["files"][key]["file"].str()).as<double>();
        };
        clipLo = vec ("clip_lo");
        clipHi = vec ("clip_hi");
        eMean = vec ("e_mean");
        eStd = vec ("e_std");
        muR = vec ("mu_r");
        W = vec ("W");
        const auto xyz = vec ("grid_xyz_um");
        const auto pred = vec ("grid_pred_std");
        const size_t V = xyz.size() / 3;
        if (pred.size() != V * static_cast<size_t> (F) || W.size() != static_cast<size_t> (F * F))
            throw std::runtime_error ("bundle arrays have inconsistent sizes");

        // Lattice: node index of every candidate voxel (ephys_only.DenseGrid).
        stepUm = m["grid"]["step_um"].num();
        voxelUm.resize (V);
        originUm = { kInf, kInf, kInf };
        for (size_t v = 0; v < V; ++v)
        {
            voxelUm[v] = { xyz[3 * v], xyz[3 * v + 1], xyz[3 * v + 2] };
            for (int k = 0; k < 3; ++k)
                originUm[static_cast<size_t> (k)] = std::min (originUm[static_cast<size_t> (k)], voxelUm[v][static_cast<size_t> (k)]);
        }
        std::vector<std::array<int, 3>> ijk (V);
        shape = { 0, 0, 0 };
        for (size_t v = 0; v < V; ++v)
            for (int k = 0; k < 3; ++k)
            {
                ijk[v][static_cast<size_t> (k)] = static_cast<int> (std::lround ((voxelUm[v][static_cast<size_t> (k)] - originUm[static_cast<size_t> (k)]) / stepUm));
                shape[static_cast<size_t> (k)] = std::max (shape[static_cast<size_t> (k)], ijk[v][static_cast<size_t> (k)] + 1);
            }
        nodeToVoxel.assign (static_cast<size_t> (shape[0]) * static_cast<size_t> (shape[1]) * static_cast<size_t> (shape[2]), -1);
        for (size_t v = 0; v < V; ++v)
            nodeToVoxel[(static_cast<size_t> (ijk[v][0]) * static_cast<size_t> (shape[1]) + static_cast<size_t> (ijk[v][1])) * static_cast<size_t> (shape[2]) + static_cast<size_t> (ijk[v][2])] = static_cast<int> (v);

        // Whitened lattice predictions wg = G W^T.
        wg.assign (V * static_cast<size_t> (F), 0.0);
        for (size_t v = 0; v < V; ++v)
            for (int i = 0; i < F; ++i)
            {
                double s = 0.0;
                for (int j = 0; j < F; ++j)
                    s += pred[v * static_cast<size_t> (F) + static_cast<size_t> (j)] * W[static_cast<size_t> (i * F + j)];
                wg[v * static_cast<size_t> (F) + static_cast<size_t> (i)] = s;
            }

        // Region labels.
        const auto& lab = m["labels"];
        labelStepM = lab["step_um"].num() * 1e-6;
        for (int k = 0; k < 3; ++k)
        {
            labelOriginM[static_cast<size_t> (k)] = lab["origin_m"][static_cast<size_t> (k)].num();
            labelShape[static_cast<size_t> (k)] = static_cast<int> (lab["shape"][static_cast<size_t> (k)].num());
        }
        labelsCosmos = loadNpy (dir + lab["files"]["Cosmos"]["file"].str()).as<uint16_t>();
        labelsBeryl = loadNpy (dir + lab["files"]["Beryl"]["file"].str()).as<uint16_t>();
        const Json table = Json::parse (readText (dir + "regions.json"));
        regions.clear();
        regionRow.clear();
        for (const auto& r : table.items)
        {
            RegionInfo info;
            info.index = static_cast<int> (r["index"].num());
            info.id = static_cast<long long> (r["id"].num());
            info.acronym = r["acronym"].str();
            info.name = r["name"].str();
            for (int k = 0; k < 3; ++k)
                info.rgb[static_cast<size_t> (k)] = static_cast<uint8_t> (r["rgb"][static_cast<size_t> (k)].num());
            regionRow[info.index] = regions.size();
            regions.push_back (info);
        }
        return true;
    }
    catch (const std::exception& e)
    {
        error = e.what();
        return false;
    }
}

int LocalizerBundle::regionAt (const V3& xyzM, bool beryl) const
{
    if (! finite3 (xyzM) || labelsCosmos.empty())
        return -1;
    std::array<size_t, 3> idx {};
    for (int k = 0; k < 3; ++k)
    {
        const long i = std::lround ((xyzM[static_cast<size_t> (k)] - labelOriginM[static_cast<size_t> (k)]) / labelStepM);
        idx[static_cast<size_t> (k)] = static_cast<size_t> (std::clamp<long> (i, 0, labelShape[static_cast<size_t> (k)] - 1));
    }
    const size_t flat = (idx[0] * static_cast<size_t> (labelShape[1]) + idx[1]) * static_cast<size_t> (labelShape[2]) + idx[2];
    const int region = static_cast<int> ((beryl ? labelsBeryl : labelsCosmos)[flat]);
    const auto it = regionRow.find (region);
    return it == regionRow.end() ? -1 : static_cast<int> (it->second);
}

int LocalizerBundle::nearestVoxel (const V3& pUm, double& distUm) const
{
    // Search lattice nodes in growing Chebyshev shells around the point's node until no closer
    // candidate voxel can exist.
    std::array<int, 3> base {};
    for (int k = 0; k < 3; ++k)
        base[static_cast<size_t> (k)] = static_cast<int> (std::floor ((pUm[static_cast<size_t> (k)] - originUm[static_cast<size_t> (k)]) / stepUm));
    const int maxR = std::max ({ shape[0], shape[1], shape[2] }) + std::max ({ std::abs (base[0]), std::abs (base[1]), std::abs (base[2]) }) + 2;
    double best = kInf;
    int bestId = -1;
    for (int r = 0; r <= maxR; ++r)
    {
        if (bestId >= 0 && (r - 1) * stepUm > best)
            break;
        for (int di = -r; di <= r; ++di)
            for (int dj = -r; dj <= r; ++dj)
                for (int dk = -r; dk <= r; ++dk)
                {
                    if (std::max ({ std::abs (di), std::abs (dj), std::abs (dk) }) != r)
                        continue;
                    const int i = base[0] + di, j = base[1] + dj, k = base[2] + dk;
                    if (i < 0 || j < 0 || k < 0 || i >= shape[0] || j >= shape[1] || k >= shape[2])
                        continue;
                    const int v = nodeToVoxel[(static_cast<size_t> (i) * static_cast<size_t> (shape[1]) + static_cast<size_t> (j)) * static_cast<size_t> (shape[2]) + static_cast<size_t> (k)];
                    if (v < 0)
                        continue;
                    const double d = norm (sub (pUm, voxelUm[static_cast<size_t> (v)]));
                    if (d < best)
                    {
                        best = d;
                        bestId = v;
                    }
                }
    }
    distUm = best;
    return bestId;
}

void LocalizerBundle::locate (const V3& pUm, std::array<int, 8>& ids, std::array<double, 8>& w, double& dnnUm) const
{
    const bool finite = finite3 (pUm);
    const V3 q = finite ? pUm : originUm;
    std::array<double, 3> f {};
    std::array<int, 3> b {};
    for (int k = 0; k < 3; ++k)
    {
        const double u = (q[static_cast<size_t> (k)] - originUm[static_cast<size_t> (k)]) / stepUm;
        b[static_cast<size_t> (k)] = static_cast<int> (std::floor (u));
        f[static_cast<size_t> (k)] = u - b[static_cast<size_t> (k)];
    }
    double tot = 0.0;
    int c = 0;
    for (int ci = 0; ci < 2; ++ci)
        for (int cj = 0; cj < 2; ++cj)
            for (int ck = 0; ck < 2; ++ck, ++c)
            {
                const int i = b[0] + ci, j = b[1] + cj, k = b[2] + ck;
                int vid = -1;
                if (i >= 0 && j >= 0 && k >= 0 && i < shape[0] && j < shape[1] && k < shape[2])
                    vid = nodeToVoxel[(static_cast<size_t> (i) * static_cast<size_t> (shape[1]) + static_cast<size_t> (j)) * static_cast<size_t> (shape[2]) + static_cast<size_t> (k)];
                const double wf = (ci ? f[0] : 1.0 - f[0]) * (cj ? f[1] : 1.0 - f[1]) * (ck ? f[2] : 1.0 - f[2]);
                ids[static_cast<size_t> (c)] = vid >= 0 ? vid : 0;
                w[static_cast<size_t> (c)] = vid >= 0 ? wf : 0.0;
                tot += w[static_cast<size_t> (c)];
            }
    double d = kInf;
    const int nearest = nearestVoxel (q, d);
    if (tot <= 1e-12)
    {
        std::fill (w.begin(), w.end(), 0.0);
        ids[0] = nearest;
        w[0] = 1.0;
    }
    else
    {
        for (auto& x : w)
            x /= tot;
    }
    dnnUm = finite ? d : kInf;
}

// =============================================================================================
// Localizer
// =============================================================================================

EphysOnlyLocalizer::Prepared EphysOnlyLocalizer::prepare (const std::vector<double>& features, int nChannels) const
{
    const int F = b.F;
    Prepared prep;
    prep.nChannels = nChannels;
    prep.valid.assign (static_cast<size_t> (nChannels), 0);
    for (int c = 0; c < nChannels; ++c)
    {
        bool finite = true, allZero = true;
        for (int f = 0; f < F; ++f)
        {
            const double x = features[static_cast<size_t> (c * F + f)];
            if (! std::isfinite (x))
                finite = false;
            else if (x != 0.0)
                allZero = false;
        }
        prep.valid[static_cast<size_t> (c)] = finite && ! allZero;
    }
    // Depth bins (ephys_only.DepthBins.from_valid), float32 arithmetic as numpy.
    const float spacing = static_cast<float> (b.score.channelSpacingUm), binUm = static_cast<float> (b.score.depthBinUm);
    float d0 = std::numeric_limits<float>::infinity(), d1 = -std::numeric_limits<float>::infinity();
    for (int c = 0; c < nChannels; ++c)
        if (prep.valid[static_cast<size_t> (c)])
        {
            d0 = std::min (d0, static_cast<float> (c) * spacing);
            d1 = std::max (d1, static_cast<float> (c) * spacing);
        }
    if (! std::isfinite (d0))
        throw std::runtime_error ("no channel carries a signal");
    std::vector<float> edges;
    const auto nEdges = static_cast<size_t> (std::ceil ((static_cast<double> (d1) + binUm - d0) / binUm));
    for (size_t e = 0; e < nEdges; ++e)
        edges.push_back (d0 + static_cast<float> (e) * binUm);
    if (edges.size() < 2)
        edges = { d0, d0 + binUm };
    const int nBins = static_cast<int> (edges.size()) - 1;
    std::vector<std::vector<int>> rows (static_cast<size_t> (nBins));
    for (int c = 0; c < nChannels; ++c)
    {
        const float depth = static_cast<float> (c) * spacing;
        const int id = static_cast<int> (std::upper_bound (edges.begin(), edges.end(), depth) - edges.begin()) - 1;
        const int bin = std::clamp (id, 0, nBins - 1);
        if (prep.valid[static_cast<size_t> (c)])
            rows[static_cast<size_t> (bin)].push_back (c);
    }
    for (int bIdx = 0; bIdx < nBins; ++bIdx)
        if (static_cast<int> (rows[static_cast<size_t> (bIdx)].size()) >= b.score.minChannelsPerBin)
        {
            prep.binRows.push_back (rows[static_cast<size_t> (bIdx)]);
            prep.centersUm.push_back (0.5 * (static_cast<double> (edges[static_cast<size_t> (bIdx)]) + edges[static_cast<size_t> (bIdx) + 1]));
        }
    if (prep.binRows.size() < 2)
        throw std::runtime_error ("fewer than 2 depth bins carry a signal");

    // Clip to the release percentiles, standardise (float32), bin means, whiten.
    const size_t B = prep.binRows.size();
    std::vector<double> x (B * static_cast<size_t> (F), 0.0);
    for (size_t bi = 0; bi < B; ++bi)
        for (int f = 0; f < F; ++f)
        {
            double s = 0.0;
            int n = 0;
            for (int c : prep.binRows[bi])
            {
                const float v = std::clamp (static_cast<float> (features[static_cast<size_t> (c * F + f)]),
                                            static_cast<float> (b.clipLo[static_cast<size_t> (f)]),
                                            static_cast<float> (b.clipHi[static_cast<size_t> (f)]));
                const float z = (v - static_cast<float> (b.eMean[static_cast<size_t> (f)]))
                                / std::max (static_cast<float> (b.eStd[static_cast<size_t> (f)]), 1e-6f);
                s += z;
                ++n;
            }
            x[bi * static_cast<size_t> (F) + static_cast<size_t> (f)] = s / n;
        }
    prep.wx.assign (B * static_cast<size_t> (F), 0.0);
    for (size_t bi = 0; bi < B; ++bi)
        for (int i = 0; i < F; ++i)
        {
            double s = 0.0;
            for (int j = 0; j < F; ++j)
                s += (x[bi * static_cast<size_t> (F) + static_cast<size_t> (j)] - b.muR[static_cast<size_t> (j)]) * b.W[static_cast<size_t> (i * F + j)];
            prep.wx[bi * static_cast<size_t> (F) + static_cast<size_t> (i)] = s;
        }
    return prep;
}

EphysOnlyLocalizer::PriorRef EphysOnlyLocalizer::priorReference (const std::vector<V3>& plannedM, const Prepared& prep) const
{
    PriorRef ref;
    // Planned axis and span (ephys_only.plan_axis_um).
    std::vector<V3> pts;
    for (const auto& p : plannedM)
        if (validPosition (p))
            pts.push_back (mirrorLeft (mul (p, 1e6)));
    if (pts.size() >= 2)
    {
        ref.axisPlan = unitVector (principalAxis (pts));
        ref.spanUm = norm (sub (pts.back(), pts.front()));
    }
    else
    {
        ref.axisPlan = { 0.0, 0.0, -1.0 };
        ref.spanUm = 0.0;
    }
    // Binned planned positions.
    std::vector<V3> P;
    for (const auto& rows : prep.binRows)
    {
        V3 s { 0, 0, 0 };
        int n = 0;
        for (int c : rows)
        {
            const V3 p = plannedM[static_cast<size_t> (c)];
            if (validPosition (p))
            {
                s = add (s, mirrorLeft (mul (p, 1e6)));
                ++n;
            }
        }
        P.push_back (n ? mul (s, 1.0 / n) : V3 { kNaN, kNaN, kNaN });
    }
    std::vector<V3> okp;
    for (const auto& p : P)
        if (finite3 (p))
            okp.push_back (p);
    ref.degenerate = ref.spanUm < b.score.planMinSpanUm || okp.size() < 2;
    if (ref.degenerate)
    {
        std::vector<double> depth;
        for (const auto& rows : prep.binRows)
            depth.push_back (std::accumulate (rows.begin(), rows.end(), 0.0) / rows.size() * b.score.channelSpacingUm);
        ref.LPlan = depth.size() > 1 ? *std::max_element (depth.begin(), depth.end()) - *std::min_element (depth.begin(), depth.end()) : 3770.0;
        V3 c { 0, 0, 0 };
        for (const auto& p : okp)
            c = add (c, p);
        ref.cPlan = okp.empty() ? V3 { 0, 0, 0 } : mul (c, 1.0 / okp.size());
    }
    else
    {
        ref.LPlan = norm (sub (okp.back(), okp.front()));
        V3 c { 0, 0, 0 };
        for (const auto& p : okp)
            c = add (c, p);
        ref.cPlan = mul (c, 1.0 / okp.size());
    }
    ref.LPlan = std::max (ref.LPlan, 1.0);
    ref.widen = (ref.degenerate ? b.score.degeneratePriorWiden : 1.0) * b.widen;
    for (int k = 0; k < b.score.nSections; ++k)
        ref.sections.push_back (sectionIndices (prep.centersUm, k, b.score.nSections));
    return ref;
}

double EphysOnlyLocalizer::priorCost (const std::vector<V3>& tr, const PriorRef& ref) const
{
    const double w = ref.widen;
    V3 mean { 0, 0, 0 };
    for (const auto& p : tr)
        mean = add (mean, p);
    mean = mul (mean, 1.0 / tr.size());
    const V3 dc = sub (mean, ref.cPlan);
    double cost = 0.0;
    for (int k = 0; k < 3; ++k)
    {
        const auto& c = b.prior.centre[static_cast<size_t> (k)];
        cost += rhoT (dc[static_cast<size_t> (k)], c[0], c[1], c[2] * w);
    }
    const V3 chord = sub (tr.back(), tr.front());
    const double chordLen = norm (chord);
    const V3 at = mul (chord, 1.0 / std::max (chordLen, 1e-12));
    const double angle = rad2deg (std::acos (std::clamp (dot (at, ref.axisPlan), -1.0, 1.0)));
    cost += rhoT2 (angle, b.prior.axisNu, b.prior.axisS * w);
    const double ratio = chordLen / ref.LPlan;
    const double z = (ratio - b.prior.lengthMedian) / (b.prior.lengthSd * w);
    cost += 0.5 * z * z;
    // Bends: angles between the principal directions of consecutive sections.
    std::vector<V3> dirs;
    std::vector<bool> has;
    for (const auto& idx : ref.sections)
    {
        if (idx.size() < 2)
        {
            dirs.push_back ({ kNaN, kNaN, kNaN });
            has.push_back (false);
            continue;
        }
        std::vector<V3> pts;
        for (int i : idx)
            pts.push_back (tr[static_cast<size_t> (i)]);
        dirs.push_back (principalAxis (pts));
        has.push_back (true);
    }
    for (size_t k = 0; k + 1 < dirs.size(); ++k)
    {
        double ang = 0.0;
        if (has[k] && has[k + 1])
            ang = rad2deg (std::acos (std::clamp (dot (dirs[k], dirs[k + 1]), -1.0, 1.0)));
        if (! std::isfinite (ang))
            ang = 0.0;
        cost += rhoT2 (ang, b.prior.bendNu, b.prior.bendS);
    }
    return cost;
}

EphysOnlyLocalizer::Components EphysOnlyLocalizer::components (const std::vector<V3>& traceUm, const Prepared& prep, const PriorRef& ref) const
{
    Components c;
    const int F = b.F;
    bool finite = true;
    std::array<int, 8> ids {};
    std::array<double, 8> w {};
    std::vector<double> g (static_cast<size_t> (F));
    for (size_t bi = 0; bi < traceUm.size(); ++bi)
    {
        if (! finite3 (traceUm[bi]))
            finite = false;
        double dnn = kInf;
        b.locate (traceUm[bi], ids, w, dnn);
        std::fill (g.begin(), g.end(), 0.0);
        for (int k = 0; k < 8; ++k)
            if (w[static_cast<size_t> (k)] != 0.0)
            {
                const double* row = &b.wg[static_cast<size_t> (ids[static_cast<size_t> (k)]) * static_cast<size_t> (F)];
                for (int f = 0; f < F; ++f)
                    g[static_cast<size_t> (f)] += w[static_cast<size_t> (k)] * row[f];
            }
        double m = 0.0;
        for (int f = 0; f < F; ++f)
        {
            const double e = prep.wx[bi * static_cast<size_t> (F) + static_cast<size_t> (f)] - g[static_cast<size_t> (f)];
            m += e * e;
        }
        c.likSum += 0.5 * (b.nu + F) * std::log1p (m / (b.nu - 2.0));
        const double out = std::isfinite (dnn) ? std::pow (std::max (0.0, dnn - b.score.outFreeUm) / b.score.outScaleUm, 2.0) : 1e6;
        c.outside += b.score.lambdaOut * out;
    }
    c.likOverT = (c.likSum + c.outside) / b.T;
    c.prior = priorCost (traceUm, ref);
    c.raw = c.likOverT + c.prior;
    if (! finite || ! std::isfinite (c.raw))
        c.raw = kInf;
    return c;
}

std::vector<V3> EphysOnlyLocalizer::tracePoints (TraceParams& params, const std::vector<double>& centersUm, bool zNegative) const
{
    const auto& cfg = b.anneal;
    const int n = std::max (1, cfg.nSections);
    params.bends.resize (static_cast<size_t> (n - 1), { 0.0, 0.0 });
    for (auto& bend : params.bends)
        bend = clipBend (bend, cfg.maxJointAngleDeg);
    std::vector<V3> directions { probeAxisFromAngles (params.theta, params.phi, zNegative) };
    for (int j = 0; j < n - 1; ++j)
        directions.push_back (bendDirection (directions.back(), params.bends[static_cast<size_t> (j)], cfg.maxJointAngleDeg));
    const double length = std::max (params.length, 1e-6);
    const double seg = length / n;
    std::vector<V3> vertices (static_cast<size_t> (n + 1));
    vertices[0] = { -std::abs (params.x), params.y, params.z };
    double remaining = 0.5 * length;
    int section = 0;
    while (remaining > 1e-6 && section < n)
    {
        const double step = std::min (seg, remaining);
        vertices[0] = sub (vertices[0], mul (directions[static_cast<size_t> (section)], step));
        remaining -= step;
        ++section;
    }
    for (int s = 0; s < n; ++s)
        vertices[static_cast<size_t> (s + 1)] = add (vertices[static_cast<size_t> (s)], mul (directions[static_cast<size_t> (s)], seg));
    const double d0 = *std::min_element (centersUm.begin(), centersUm.end());
    const double d1 = *std::max_element (centersUm.begin(), centersUm.end());
    std::vector<V3> points;
    for (double depth : centersUm)
    {
        const double arc = std::clamp ((depth - d0) / std::max (d1 - d0, 1e-6), 0.0, 1.0) * length;
        const int sid = std::min (static_cast<int> (arc / seg), n - 1);
        const double local = arc - sid * seg;
        V3 p = add (vertices[static_cast<size_t> (sid)], mul (directions[static_cast<size_t> (sid)], local));
        p[0] = -std::abs (p[0]);
        points.push_back (p);
    }
    return points;
}

TraceParams EphysOnlyLocalizer::plannedSeed (const std::vector<V3>& plannedM, bool& zNegative) const
{
    const auto& cfg = b.anneal;
    std::vector<V3> pts;
    for (const auto& p : plannedM)
    {
        const V3 q = mirrorLeft (p);
        if (validPosition (q))
            pts.push_back (mul (q, 1e6));
    }
    V3 c { 0, 0, 0 };
    for (const auto& p : pts)
        c = add (c, p);
    c = mul (c, 1.0 / pts.size());
    V3 axis = principalAxis (pts);
    axis = mul (axis, 1.0 / std::max (norm (axis), 1e-12));
    zNegative = ! (axis[2] >= 0);
    TraceParams params;
    params.theta = std::clamp (rad2deg (std::acos (std::clamp (std::abs (axis[2]), 0.0, 1.0))), cfg.minAngleDeg, cfg.maxAngleDeg);
    double phi = std::fmod (rad2deg (std::atan2 (axis[1], axis[0])), 360.0);
    params.phi = phi < 0 ? phi + 360.0 : phi;
    params.length = std::clamp (cfg.probeLengthUm, cfg.minProbeLengthUm, cfg.maxProbeLengthUm);
    params.x = -std::abs (c[0]);
    params.y = c[1];
    params.z = c[2];
    params.bends.assign (static_cast<size_t> (cfg.nSections - 1), { 0.0, 0.0 });
    return params;
}

double EphysOnlyLocalizer::mad (std::vector<double> values)
{
    values.erase (std::remove_if (values.begin(), values.end(), [] (double v) { return ! std::isfinite (v); }), values.end());
    if (values.empty())
        return kNaN;
    const double med = median (values);
    for (auto& v : values)
        v = std::abs (v - med);
    return 1.4826 * median (values);
}

std::vector<std::vector<V3>> EphysOnlyLocalizer::neighbourhoodTraces (const std::vector<V3>& baseUm, int n, uint64_t seed) const
{
    std::mt19937_64 rng (seed);
    std::student_t_distribution<double> t3 (3.0);
    std::normal_distribution<double> normal (0.0, 1.0);
    V3 c { 0, 0, 0 };
    for (const auto& p : baseUm)
        c = add (c, p);
    c = mul (c, 1.0 / baseUm.size());
    const V3 axis = principalAxis (baseUm);
    std::vector<std::vector<V3>> out;
    for (int i = 0; i < n; ++i)
    {
        V3 shift { t3 (rng) * 300.0, t3 (rng) * 500.0, t3 (rng) * 500.0 };
        const double nr = norm (shift);
        if (nr > 2500.0)
            shift = mul (shift, 2500.0 / nr);
        V3 om { normal (rng), normal (rng), normal (rng) };
        om = sub (om, mul (axis, dot (om, axis)));
        om = mul (om, 1.0 / std::max (norm (om), 1e-12));
        const double angle = deg2rad (std::abs (t3 (rng)) * 5.0);
        double stretch = std::clamp (1.0 + 0.05 * normal (rng), 0.85, 1.2);
        // Rodrigues rotation about om by angle.
        std::vector<V3> tr;
        for (const auto& p : baseUm)
        {
            const V3 d = sub (p, c);
            const V3 r = add (add (mul (d, std::cos (angle)), mul (cross (om, d), std::sin (angle))),
                              mul (om, dot (om, d) * (1.0 - std::cos (angle))));
            tr.push_back (mirrorLeft (add (mul (r, stretch), add (c, shift))));
        }
        out.push_back (tr);
    }
    return out;
}

LocalizationOutput EphysOnlyLocalizer::localize (const std::vector<double>& features, int nChannels, const std::vector<V3>& plannedM, long long seed) const
{
    const auto t0 = std::chrono::steady_clock::now();
    LocalizationOutput out;
    try
    {
        // Planned trajectory checks (ephys_only.check_planned_trajectory).
        std::vector<V3> okPlan;
        for (const auto& p : plannedM)
            if (validPosition (p))
                okPlan.push_back (p);
        if (okPlan.size() < 2 || norm (sub (okPlan.back(), okPlan.front())) * 1e6 < b.score.minUsablePlanSpanUm)
            throw std::runtime_error ("planned trajectory unavailable (placeholder or shorter than 100 um)");

        const auto& acfg = b.anneal;
        const Prepared prep = prepare (features, nChannels);
        const PriorRef ref = priorReference (plannedM, prep);

        bool zNegative = true;
        TraceParams seedParams = plannedSeed (plannedM, zNegative);
        const std::vector<V3> seedTrace = tracePoints (seedParams, prep.centersUm, zNegative);

        // Search scale c_p = K / MAD(S) over plan-neighbourhood trajectories.
        std::vector<double> poolRaw;
        for (const auto& tr : neighbourhoodTraces (seedTrace, b.score.scalePoolSize, static_cast<uint64_t> (seed) * 7919ULL + 17ULL))
            poolRaw.push_back (components (tr, prep, ref).raw);
        const double madS = mad (poolRaw);
        double cP = (std::isfinite (madS) && madS > 1e-12) ? b.score.scaleK / madS : 1.0;
        cP = (std::isfinite (cP) && cP > 0) ? std::clamp (cP, 1e-8, 1e8) : 1.0;

        auto scoreOf = [&] (const std::vector<V3>& tr) { return cP * components (tr, prep, ref).raw; };
        const double seedScore = scoreOf (seedTrace);

        const int nJoints = acfg.nSections - 1;
        const double lo = acfg.minAngleDeg, hi = acfg.maxAngleDeg;
        double bestOverall = kInf;
        std::vector<V3> bestTrace = seedTrace;
        TraceParams bestParams = seedParams;
        int bestRestart = 0;
        for (int r = 0; r < b.nRestarts; ++r)
        {
            std::mt19937_64 rng (static_cast<uint64_t> (seed + 1000003LL * r));
            std::uniform_real_distribution<double> uni (0.0, 1.0);
            std::normal_distribution<double> normal (0.0, 1.0);
            TraceParams cur = seedParams, best = seedParams;
            double curScore = seedScore, bestScore = seedScore;
            std::vector<V3> bestTr = seedTrace;
            int lastImprove = 0, lastBestImprove = 0, nReheats = 0, lastReheat = -1;
            double boost = 1.0, reheatStart = kNaN, prevTemp = kNaN;
            const int nIt = acfg.nIterations;
            const int decay = std::max (1, acfg.temperatureDecayIters);
            for (int it = 0; it < nIt; ++it)
            {
                const double baseTemp = expSchedule (acfg.initialTemperature, acfg.finalTemperature, std::min (1.0, static_cast<double> (it) / decay));
                const double a = static_cast<double> (it) / std::max (1, nIt - 1);
                double sigXyz = expSchedule (acfg.initXyzSigmaUm, acfg.finalXyzSigmaUm, a);
                double sigTheta = expSchedule (acfg.initAngleSigmaDeg, acfg.finalAngleSigmaDeg, a);
                double sigPhi = expSchedule (acfg.initAzimuthSigmaDeg, acfg.finalAzimuthSigmaDeg, a);
                double sigLen = expSchedule (acfg.initLengthSigmaUm, acfg.finalLengthSigmaUm, a);
                double sigJoint = expSchedule (acfg.initJointAngleSigmaDeg, acfg.finalJointAngleSigmaDeg, a);
                if ((it - lastImprove) >= acfg.reheatPatience && nReheats < acfg.maxReheats)
                {
                    ++nReheats;
                    lastImprove = lastReheat = it;
                    boost = acfg.reheatSigmaFactor;
                }
                else
                {
                    boost = std::max (1.0, boost * acfg.reheatDecay);
                }
                double temp = baseTemp;
                if (lastReheat >= 0)
                {
                    if (it == lastReheat)
                    {
                        const double previous = std::isfinite (prevTemp) ? prevTemp : baseTemp;
                        reheatStart = std::max (previous * acfg.reheatTemperatureFactor, acfg.reheatMinTemperature);
                    }
                    temp = expSchedule (reheatStart, acfg.finalTemperature, std::min (1.0, static_cast<double> (it - lastReheat) / decay));
                }
                if (boost > 1.0)
                {
                    sigXyz *= boost;
                    sigTheta *= boost;
                    sigPhi *= boost;
                    sigLen *= boost;
                    sigJoint *= boost;
                }
                // Proposal.
                const double u = uni (rng);
                int move = 3;
                double acc = 0.0;
                for (int k = 0; k < 4; ++k)
                {
                    acc += acfg.moveProbabilities[static_cast<size_t> (k)];
                    if (u < acc)
                    {
                        move = k;
                        break;
                    }
                }
                TraceParams prop = cur;
                if (move == 0)
                {
                    prop.x = -std::abs (prop.x + normal (rng) * sigXyz / std::sqrt (3.0));
                    prop.y += normal (rng) * sigXyz / std::sqrt (3.0);
                    prop.z += normal (rng) * sigXyz / std::sqrt (3.0);
                }
                else if (move == 1)
                {
                    if (uni (rng) < 0.5)
                        prop.theta = std::clamp (cur.theta + normal (rng) * sigTheta, lo, hi);
                    else
                        prop.phi = std::fmod (std::fmod (cur.phi + normal (rng) * sigPhi, 360.0) + 360.0, 360.0);
                }
                else if (move == 2)
                {
                    prop.length = std::clamp (cur.length + normal (rng) * sigLen, acfg.minProbeLengthUm, acfg.maxProbeLengthUm);
                }
                else if (nJoints > 0)
                {
                    const int j = static_cast<int> (uni (rng) * nJoints) % nJoints;
                    auto& bend = prop.bends[static_cast<size_t> (j)];
                    bend[0] += normal (rng) * sigJoint;
                    bend[1] += normal (rng) * sigJoint;
                    bend = clipBend (bend, acfg.maxJointAngleDeg);
                }
                prop.x = -std::abs (prop.x);
                prop.theta = std::clamp (prop.theta, lo, hi);
                const std::vector<V3> propTrace = tracePoints (prop, prep.centersUm, zNegative);
                const double s = scoreOf (propTrace);
                const double delta = s - curScore;
                const double pAccept = delta <= 0.0 ? 1.0 : std::exp (-delta / std::max (temp, 1e-12));
                if (delta <= 0.0 || uni (rng) < pAccept)
                {
                    cur = prop;
                    curScore = s;
                    if (s < bestScore)
                    {
                        best = prop;
                        bestScore = s;
                        bestTr = propTrace;
                        lastImprove = lastBestImprove = it;
                    }
                }
                prevTemp = temp;
                if ((it - lastBestImprove) >= acfg.earlyStopPatience)
                    break;
            }
            if (bestScore < bestOverall)
            {
                bestOverall = bestScore;
                bestTrace = bestTr;
                bestParams = best;
                bestRestart = r;
            }
        }
        if (! std::isfinite (bestOverall))
            throw std::runtime_error ("every annealing restart ended with a non-finite score");

        // Channel positions along the trace (ephys_only.channel_positions).
        const auto& d = prep.centersUm;
        const size_t B = d.size();
        std::vector<V3> channelUm (static_cast<size_t> (nChannels));
        for (int c = 0; c < nChannels; ++c)
        {
            const double q = c * b.score.channelSpacingUm;
            V3 p {};
            if (q < d.front())
                p = add (bestTrace[0], mul (sub (bestTrace[1], bestTrace[0]), (q - d.front()) / (d[1] - d[0])));
            else if (q > d.back())
                p = add (bestTrace[B - 1], mul (sub (bestTrace[B - 1], bestTrace[B - 2]), (q - d.back()) / (d[B - 1] - d[B - 2])));
            else
            {
                const size_t k = static_cast<size_t> (std::upper_bound (d.begin(), d.end(), q) - d.begin());
                const size_t k1 = std::min (k, B - 1), k0 = k1 == 0 ? 0 : k1 - 1;
                const double t = d[k1] > d[k0] ? (q - d[k0]) / (d[k1] - d[k0]) : 0.0;
                p = add (bestTrace[k0], mul (sub (bestTrace[k1], bestTrace[k0]), t));
            }
            channelUm[static_cast<size_t> (c)] = p;
        }
        // Hemisphere of the plan.
        std::vector<double> xs;
        for (const auto& p : okPlan)
            xs.push_back (p[0]);
        const bool right = median (xs) > 0;
        auto toMetres = [right] (V3 pUm) {
            V3 m = mul (pUm, 1e-6);
            if (right)
                m[0] = std::abs (m[0]);
            return m;
        };
        for (const auto& p : channelUm)
            out.channelXyzM.push_back (toMetres (p));
        for (const auto& p : bestTrace)
            out.traceXyzM.push_back (toMetres (p));
        out.valid = prep.valid;
        out.binOfChannel.assign (static_cast<size_t> (nChannels), -1);
        for (size_t bi = 0; bi < prep.binRows.size(); ++bi)
            for (int c : prep.binRows[bi])
                out.binOfChannel[static_cast<size_t> (c)] = static_cast<int> (bi);
        out.binCentersUm = prep.centersUm;
        for (const auto& p : out.channelXyzM)
        {
            out.cosmosRow.push_back (b.regionAt (p, false));
            out.berylRow.push_back (b.regionAt (p, true));
        }
        out.score = bestOverall;
        out.rawScore = bestOverall / cP;
        out.cP = cP;
        out.madS = madS;
        out.T = b.T;
        out.widen = ref.widen;
        out.plannedScore = seedScore;
        out.thetaDeg = bestParams.theta;
        out.phiDeg = bestParams.phi;
        out.lengthUm = bestParams.length;
        out.nBins = static_cast<int> (B);
        out.bestRestart = bestRestart;
        out.ok = true;
    }
    catch (const std::exception& e)
    {
        out.ok = false;
        out.error = e.what();
    }
    out.seconds = std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
    return out;
}
} // namespace ProbeViewer::Localization::Realtime
