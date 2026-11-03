#pragma once
#include "LocalizationConfig.hpp"
#include "LocalizationTypes.hpp"

#include <cstdint>
#include <functional>
#include <limits>
#include <vector>

namespace ProbeViewer::Localization
{
struct BandWindow
{
    juce::String streamKey;

    // sampleRate is the rate actually used for capture/feature computation.
    // reportedSampleRate preserves the original stream metadata so debug can
    // distinguish a real 2500-Hz stream from a 2500-Hz fallback.
    double sampleRate = 0.0;
    double reportedSampleRate = std::numeric_limits<double>::quiet_NaN();

    int numChannels = 0;
    int numSamples = 0;
    std::vector<int> electrodeIndices;
    std::vector<float> channelMajorSamples;
};

struct LocalizationInput
{
    BandWindow ap;
    BandWindow lf;
    ProbeIdentity probe;
    juce::String probeName;

    // Temporal robustness batch metadata. Raw windows are processed one at a
    // time so that three long AP/LF windows are never resident in RAM together.
    int temporalWindowIndex = 0; // zero-based
    int temporalWindowCount = 1;
};

struct FeatureMatrix
{
    int rows = 0;
    int cols = 0;
    std::vector<juce::String> names;
    std::vector<float> values;

    // ibldsp bad-channel labels per source channel (0 good, 1 dead, 2 noisy, 3 outside the
    // brain), detected on the AP snippet; empty when no AP was processed.
    std::vector<int> channelLabels;

    float at (int row, int col) const
    {
        return values[static_cast<size_t> (row * cols + col)];
    }
};


using FeatureProgressCallback =
    std::function<void (float, const juce::String&)>;

class FeatureExtractor
{
public:
    // The input is intentionally non-const. Feature extraction owns the
    // captured AP/LF windows for this localization job and preprocesses them
    // in-place to avoid allocating full-size duplicate arrays.
    FeatureMatrix extract (
        LocalizationInput&,
        const LocalizationConfig&,
        FeatureProgressCallback progressCallback = {}) const;

private:
    static double rms (const float*, int);

    static void scaleInPlace (
        std::vector<float>&,
        float scale);

    static void commonMedianReference (
        std::vector<float>&,
        int,
        int,
        const std::vector<int>& labels);

    // IBL-compatible temporal preprocessing. For the standard IBL rates
    // (AP=30 kHz, LF=2.5 kHz), these use the same scipy.signal.butter SOS
    // coefficients and scipy-style odd-reflection sosfiltfilt boundary handling.
    static void iblApButterworthZeroPhase (
        std::vector<float>&,
        int channels,
        int samples,
        double sampleRate);

    static void iblLfButterworthZeroPhase (
        std::vector<float>&,
        int channels,
        int samples,
        double sampleRate);

    // NP1 ADC correction matching neuropixel.adc_shifts(version=1) followed by
    // ibldsp.fourier.fshift(..., axis=1). The original hardware/electrode index
    // is used when available.
    static void applyNp1SampleShift (
        std::vector<float>&,
        int channels,
        int samples,
        const std::vector<int>& electrodeIndices);

    // Retained only as a safe fallback for non-standard sample rates.
    static void temporalHighpassZeroPhase (
        std::vector<float>&,
        int,
        int,
        double,
        double,
        int);

    static void temporalLowpassZeroPhase (
        std::vector<float>&,
        int,
        int,
        double,
        double,
        int);

    static double cosineBandWeight (
        double,
        double,
        double);

    static std::vector<double> oneSidedPeriodogram (
        const float*,
        int,
        double,
        std::vector<double>&);

    static double bandPowerDb (
        const std::vector<double>&,
        const std::vector<double>&,
        double,
        double);
};
}
