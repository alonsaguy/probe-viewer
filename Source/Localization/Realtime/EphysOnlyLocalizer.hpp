#pragma once

// C++ port of the IBL ephys-only probe localizer (ibleatools: ephysatlas.alignment.ephys_only,
// method student_t_mahalanobis_plan_posterior), running on a bundle exported by
// ephysatlas.realtime.bundle (lattice predictions, whitening, plan prior, annealing settings,
// region labels). Plain C++17: no JUCE, so it can be tested on its own (LocalizerSelfTest).
//
// Positions: IBL-Allen coordinates (origin bregma; x ML, y AP, z DV). Inputs and outputs in
// metres, row order (top of the probe first); internally µm in the left hemisphere, like the
// Python implementation.

#include <array>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace ProbeViewer::Localization::Realtime
{
using V3 = std::array<double, 3>;

struct ScoreConfig
{
    double gridUm = 200.0;
    double channelSpacingUm = 10.0;
    double depthBinUm = 100.0;
    int minChannelsPerBin = 2;
    double lambdaOut = 20.0;
    double outFreeUm = 300.0;
    double outScaleUm = 500.0;
    double planMinSpanUm = 1000.0;
    double degeneratePriorWiden = 5.0;
    int nSections = 4;
    double minUsablePlanSpanUm = 100.0;
    int scalePoolSize = 64;
    double scaleK = 0.167;
};

struct AnnealConfig
{
    int nIterations = 400;
    double initialTemperature = 0.02;
    double finalTemperature = 1e-4;
    int temperatureDecayIters = 400;
    double initXyzSigmaUm = 500.0, finalXyzSigmaUm = 100.0;
    double initAngleSigmaDeg = 5.0, finalAngleSigmaDeg = 1.0;
    double initAzimuthSigmaDeg = 20.0, finalAzimuthSigmaDeg = 4.0;
    double initLengthSigmaUm = 300.0, finalLengthSigmaUm = 50.0;
    double initJointAngleSigmaDeg = 1.8, finalJointAngleSigmaDeg = 0.1;
    std::array<double, 4> moveProbabilities { 0.30, 0.20, 0.10, 0.40 };
    int nSections = 4;
    double maxJointAngleDeg = 10.0;
    double minAngleDeg = 10.0, maxAngleDeg = 30.0;
    double probeLengthUm = 3840.0, minProbeLengthUm = 3500.0, maxProbeLengthUm = 4500.0;
    int reheatPatience = 80;
    int maxReheats = 1;
    double reheatTemperatureFactor = 1.5;
    double reheatMinTemperature = 0.01;
    double reheatSigmaFactor = 1.75;
    double reheatDecay = 0.97;
    int earlyStopPatience = 100;
};

struct PriorParams
{
    std::array<std::array<double, 3>, 3> centre {}; // per axis x, y, z: (nu, loc, s)
    double axisNu = 5.0, axisS = 4.5;
    double lengthMedian = 1.0, lengthSd = 0.04;
    double bendNu = 4.0, bendS = 4.2;
};

struct RegionInfo
{
    int index = 0; // row of the iblatlas region table
    long long id = 0; // Allen id
    std::string acronym, name;
    std::array<uint8_t, 3> rgb { 200, 200, 200 };
};

// Everything the localizer reads from a bundle directory.
class LocalizerBundle
{
public:
    bool load (const std::string& directory, std::string& error);

    int F = 0;
    std::vector<std::string> featureNames;
    std::vector<double> clipLo, clipHi, eMean, eStd, muR, W; // W row-major [F x F]
    double nu = 15.0, kappa = 1.0, T = 20.0, widen = 1.0;
    PriorParams prior;
    ScoreConfig score;
    AnnealConfig anneal;
    int nRestarts = 10;
    long long seed = 0;
    std::string featureSet, modelCommit;

    // Candidate lattice (µm, left hemisphere).
    double stepUm = 200.0;
    V3 originUm {};
    std::array<int, 3> shape {};
    std::vector<int> nodeToVoxel; // [nx * ny * nz], -1 where no candidate voxel
    std::vector<V3> voxelUm;      // [V]
    std::vector<double> wg;       // [V x F] whitened predictions G W^T

    // Region labels on a regular grid (for display).
    double labelStepM = 50e-6;
    V3 labelOriginM {};
    std::array<int, 3> labelShape {};
    std::vector<uint16_t> labelsCosmos, labelsBeryl;
    std::vector<RegionInfo> regions;
    std::map<int, size_t> regionRow; // region index -> row of `regions`

    // Region table row (or -1) of the region at a position (m).
    int regionAt (const V3& xyzM, bool beryl) const;

    // Trilinear corners and weights of a point (µm); weights sum to 1. dnn: distance to the nearest
    // candidate voxel (inf for non-finite points).
    void locate (const V3& pUm, std::array<int, 8>& ids, std::array<double, 8>& w, double& dnnUm) const;
    // Nearest candidate voxel to a point (µm): id and distance.
    int nearestVoxel (const V3& pUm, double& distUm) const;
};

struct TraceParams
{
    double x = 0, y = 0, z = 0, theta = 15, phi = 0, length = 3840;
    std::vector<std::array<double, 2>> bends; // [nSections - 1] (deg)
};

struct LocalizationOutput
{
    bool ok = false;
    std::string error;
    std::vector<V3> channelXyzM;   // [C] in the plan's hemisphere
    std::vector<V3> traceXyzM;     // [B]
    std::vector<uint8_t> valid;    // [C] channels with a signal
    std::vector<int> binOfChannel; // [C], -1 when in no bin
    std::vector<double> binCentersUm;
    std::vector<int> cosmosRow, berylRow; // [C] rows of LocalizerBundle::regions, -1 outside
    double score = 0, rawScore = 0, cP = 1, madS = 0, T = 0, widen = 1, plannedScore = 0;
    double thetaDeg = 0, phiDeg = 0, lengthUm = 0, seconds = 0;
    int nBins = 0, bestRestart = 0;
};

class EphysOnlyLocalizer
{
public:
    explicit EphysOnlyLocalizer (const LocalizerBundle& bundle) : b (bundle) {}

    // features: [C x F] row-major in the bundle's feature order (feature units, NaN allowed);
    // plannedM: [C] planned positions (m). Both in row order (top first).
    LocalizationOutput localize (const std::vector<double>& features,
                                 int nChannels,
                                 const std::vector<V3>& plannedM,
                                 long long seed) const;

    // --- steps, public for the self-test -------------------------------------------------
    struct Prepared
    {
        std::vector<uint8_t> valid;
        std::vector<std::vector<int>> binRows;
        std::vector<double> centersUm;
        std::vector<double> wx; // [B x F]
        int nChannels = 0;
    };
    Prepared prepare (const std::vector<double>& features, int nChannels) const;

    struct PriorRef
    {
        V3 cPlan {}, axisPlan {};
        double LPlan = 1.0, widen = 1.0, spanUm = 0.0;
        bool degenerate = false;
        std::vector<std::vector<int>> sections;
    };
    PriorRef priorReference (const std::vector<V3>& plannedM, const Prepared& prep) const;

    struct Components
    {
        double likSum = 0, outside = 0, likOverT = 0, prior = 0, raw = 0;
    };
    Components components (const std::vector<V3>& traceUm, const Prepared& prep, const PriorRef& ref) const;

    std::vector<V3> tracePoints (TraceParams& params, const std::vector<double>& centersUm, bool zNegative) const;
    TraceParams plannedSeed (const std::vector<V3>& plannedM, bool& zNegative) const;
    static double mad (std::vector<double> values);
    // Rigid + stretch perturbations of a trajectory (the scale pool), with a C++ random stream.
    std::vector<std::vector<V3>> neighbourhoodTraces (const std::vector<V3>& baseUm, int n, uint64_t seed) const;

private:
    const LocalizerBundle& b;
    double priorCost (const std::vector<V3>& traceUm, const PriorRef& ref) const;
};
} // namespace ProbeViewer::Localization::Realtime
