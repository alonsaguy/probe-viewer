// Self-test of the C++ ephys-only localizer against the Python implementation, using the
// self-test written into a localizer bundle by ephysatlas.realtime.bundle.export_localizer_bundle.
//
// Usage: probe-viewer-localizer-selftest <bundle directory>

#include "../Source/Localization/Realtime/EphysOnlyLocalizer.hpp"
#include "../Source/Localization/Realtime/MiniJson.hpp"
#include "../Source/Localization/Realtime/Npy.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

using namespace ProbeViewer::Localization::Realtime;

namespace
{
int failures = 0;

void check (const char* what, double error, double tolerance)
{
    const bool ok = std::isfinite (error) && error <= tolerance;
    std::printf ("  [%s] %-46s max error %.3g (tolerance %.3g)\n", ok ? "PASS" : "FAIL", what, error, tolerance);
    if (! ok)
        ++failures;
}

std::vector<V3> toV3 (const std::vector<double>& flat)
{
    std::vector<V3> out (flat.size() / 3);
    for (size_t i = 0; i < out.size(); ++i)
        out[i] = { flat[3 * i], flat[3 * i + 1], flat[3 * i + 2] };
    return out;
}

double relErr (double a, double b) { return std::abs (a - b) / std::max (1.0, std::abs (b)); }
} // namespace

int main (int argc, char** argv)
{
    if (argc < 2)
    {
        std::printf ("usage: %s <bundle directory>\n", argv[0]);
        return 2;
    }
    const std::string dir = argv[1];
    LocalizerBundle bundle;
    std::string error;
    if (! bundle.load (dir, error))
    {
        std::printf ("cannot load bundle: %s\n", error.c_str());
        return 2;
    }
    std::printf ("bundle %s: %d features (%s), %zu lattice voxels, T=%g, nu=%g, width=%g\n",
                 dir.c_str(), bundle.F, bundle.featureSet.c_str(), bundle.voxelUm.size(), bundle.T, bundle.nu, bundle.widen);

    std::ifstream mf (dir + "/manifest.json", std::ios::binary);
    std::stringstream ss;
    ss << mf.rdbuf();
    const Json manifest = Json::parse (ss.str());
    if (! manifest.has ("selftest"))
    {
        std::printf ("the bundle has no self-test\n");
        return 2;
    }
    const Json& st = manifest["selftest"];
    const std::string sdir = dir + "/selftest/";
    auto arr = [&] (const std::string& name) { return loadNpy (sdir + st["files"][name]["file"].str()); };

    const NpyArray feat = arr ("features");
    const int C = static_cast<int> (feat.shape[0]);
    const std::vector<double> features = feat.as<double>();
    const std::vector<V3> planned = toV3 (arr ("planned_xyz_m").as<double>());
    EphysOnlyLocalizer loc (bundle);

    // 1. Depth bins and whitened bin means.
    const auto prep = loc.prepare (features, C);
    {
        const auto binOf = arr ("bin_of_channel").as<int>();
        const auto centers = arr ("bin_centers_um").as<double>();
        std::vector<int> mine (static_cast<size_t> (C), -1);
        for (size_t bi = 0; bi < prep.binRows.size(); ++bi)
            for (int c : prep.binRows[bi])
                mine[static_cast<size_t> (c)] = static_cast<int> (bi);
        int mismatches = 0;
        for (int c = 0; c < C; ++c)
            mismatches += mine[static_cast<size_t> (c)] != binOf[static_cast<size_t> (c)];
        check ("bin of each channel (mismatches)", mismatches, 0);
        double e = centers.size() == prep.centersUm.size() ? 0.0 : 1e9;
        for (size_t i = 0; e == 0.0 && i < centers.size(); ++i)
            e = std::max (e, std::abs (centers[i] - prep.centersUm[i]));
        check ("bin centres (um)", e, 1e-3);
        const auto wx = arr ("wx").as<double>();
        double ew = 0.0;
        for (size_t i = 0; i < wx.size(); ++i)
            ew = std::max (ew, std::abs (wx[i] - prep.wx[i]));
        check ("whitened bin means wx", ew, 1e-3);
    }

    // 2. Prior reference.
    const auto ref = loc.priorReference (planned, prep);
    {
        const auto& pr = st["prior_reference"];
        double e = 0.0;
        for (size_t k = 0; k < 3; ++k)
            e = std::max (e, std::abs (pr["c_plan_um"][k].num() - ref.cPlan[k]));
        check ("plan centre (um)", e, 1e-2);
        double ea = 0.0;
        for (size_t k = 0; k < 3; ++k)
            ea = std::max (ea, std::abs (pr["axis_plan"][k].num() - ref.axisPlan[k]));
        check ("plan axis", ea, 1e-6);
        check ("plan length (um)", std::abs (pr["L_plan_um"].num() - ref.LPlan), 1e-2);
        check ("plan degenerate / widen", std::abs ((pr["degenerate"].boolean ? 1 : 0) - (ref.degenerate ? 1 : 0)) + std::abs (pr["widen"].num() - ref.widen), 0);
    }

    // 3. Seed trajectory.
    bool zNegative = true;
    TraceParams seedParams = loc.plannedSeed (planned, zNegative);
    const auto seedTrace = loc.tracePoints (seedParams, prep.centersUm, zNegative);
    {
        const auto expected = toV3 (arr ("seed_trace_um").as<double>());
        double e = 0.0;
        for (size_t i = 0; i < expected.size(); ++i)
            for (size_t k = 0; k < 3; ++k)
                e = std::max (e, std::abs (expected[i][k] - seedTrace[i][k]));
        check ("seed trajectory (um; Python works in float32)", e, 0.5);
        check ("seed z direction", (st["z_sign"].str() == "negative") != zNegative ? 1.0 : 0.0, 0);
    }

    // 4. Trajectory geometry from parameters, and score components of every trajectory.
    const NpyArray tracesArr = arr ("traces_um");
    const auto traces = tracesArr.as<double>();
    const size_t N = tracesArr.shape[0], B = tracesArr.shape[1];
    const int nRigid = static_cast<int> (st["n_rigid"].num());
    {
        const NpyArray paramsArr = arr ("trace_params");
        const auto params = paramsArr.as<double>();
        const size_t P = paramsArr.shape[0], width = paramsArr.shape[1];
        double e = 0.0;
        for (size_t i = 0; i < P; ++i)
        {
            const double* p = &params[i * width];
            TraceParams tp;
            tp.x = p[0];
            tp.y = p[1];
            tp.z = p[2];
            tp.theta = p[3];
            tp.phi = p[4];
            tp.length = p[5];
            for (size_t j = 0; j + 1 < static_cast<size_t> (bundle.anneal.nSections); ++j)
                tp.bends.push_back ({ p[6 + 2 * j], p[7 + 2 * j] });
            const auto tr = loc.tracePoints (tp, prep.centersUm, zNegative);
            for (size_t bi = 0; bi < B; ++bi)
                for (size_t k = 0; k < 3; ++k)
                    e = std::max (e, std::abs (tr[bi][k] - traces[((static_cast<size_t> (nRigid) + i) * B + bi) * 3 + k]));
        }
        check ("trajectories from parameters (um)", e, 0.05);
    }
    {
        const auto lik = arr ("lik_sum").as<double>(), outside = arr ("outside").as<double>();
        const auto prior = arr ("prior").as<double>(), raw = arr ("raw").as<double>();
        double el = 0, eo = 0, ep = 0, er = 0;
        for (size_t i = 0; i < N; ++i)
        {
            std::vector<V3> tr (B);
            for (size_t bi = 0; bi < B; ++bi)
                tr[bi] = { traces[(i * B + bi) * 3], traces[(i * B + bi) * 3 + 1], traces[(i * B + bi) * 3 + 2] };
            const auto c = loc.components (tr, prep, ref);
            el = std::max (el, relErr (c.likSum, lik[i]));
            eo = std::max (eo, relErr (c.outside, outside[i]));
            ep = std::max (ep, relErr (c.prior, prior[i]));
            er = std::max (er, relErr (c.raw, raw[i]));
        }
        check ("likelihood sum (relative)", el, 1e-4);
        check ("outside-lattice penalty (relative)", eo, 1e-4);
        check ("plan prior (relative)", ep, 1e-4);
        check ("score S (relative)", er, 1e-4);
    }

    // 5. Search scale from the same pool of trajectories.
    {
        const NpyArray poolArr = arr ("pool_traces_um");
        const auto pool = poolArr.as<double>();
        std::vector<double> raws;
        for (size_t i = 0; i < poolArr.shape[0]; ++i)
        {
            std::vector<V3> tr (B);
            for (size_t bi = 0; bi < B; ++bi)
                tr[bi] = { pool[(i * B + bi) * 3], pool[(i * B + bi) * 3 + 1], pool[(i * B + bi) * 3 + 2] };
            raws.push_back (loc.components (tr, prep, ref).raw);
        }
        const double madS = EphysOnlyLocalizer::mad (raws);
        check ("MAD of the pool scores (relative)", relErr (madS, st["mad_S"].num()), 1e-4);
        check ("search scale c_p (relative)", relErr (bundle.score.scaleK / madS, st["c_p"].num()), 1e-4);
    }

    // 6. Full search (another random stream than Python: a nearby result is expected).
    {
        const auto out = loc.localize (features, C, planned, bundle.seed);
        if (! out.ok)
        {
            std::printf ("  [FAIL] localize: %s\n", out.error.c_str());
            ++failures;
        }
        else
        {
            const auto py = toV3 (arr ("result_channel_xyz_m").as<double>());
            double d = 0.0;
            for (int c = 0; c < C; ++c)
            {
                const V3 a = out.channelXyzM[static_cast<size_t> (c)], b = py[static_cast<size_t> (c)];
                d += std::sqrt ((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
            }
            std::printf ("  [INFO] full search: %.3f s, score %.4f (Python %.4f), mean distance to the Python result %.0f um\n",
                         out.seconds, out.score, st["result_score"].num(), d / C * 1e6);
        }
    }
    std::printf (failures ? "SELF-TEST FAILED (%d checks)\n" : "SELF-TEST PASSED\n", failures);
    return failures ? 1 : 0;
}
