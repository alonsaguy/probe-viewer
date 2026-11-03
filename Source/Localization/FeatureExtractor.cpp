#include "FeatureExtractor.hpp"
#include "Realtime/BadChannels.hpp"
#include "Realtime/IblLfFeatures.hpp"

#include <thread>
#include "../kissfft/kiss_fftr.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <thread>

namespace ProbeViewer::Localization
{
namespace
{
constexpr double eps = 1.0e-30;

// Deliberately leave CPU capacity for Open Ephys, the UI and Windows.
// Increase only if profiling shows substantial headroom.
constexpr int maxFeatureThreads = 4;

// CAR is expensive at AP sample rates. Processing contiguous time blocks
// avoids scheduling hundreds of thousands of tiny parallel tasks and lets
// each worker reuse one median scratch buffer for many samples.
constexpr int medianReferenceBlockSamples = 2048;

bool is235Smooth (
    int n)
{
    if (n <= 0)
        return false;

    while ((n % 2) == 0)
        n /= 2;

    while ((n % 3) == 0)
        n /= 3;

    while ((n % 5) == 0)
        n /= 5;

    return n == 1;
}

int efficientFftLengthAtMost (
    int samples)
{
    if (samples <= 1)
        return samples;

    // Real FFT length must be even.
    int n =
        (samples % 2 == 0)
            ? samples
            : samples - 1;

    // Usually this loop only removes a handful of samples.
    while (n > 2
           && ! is235Smooth (n))
    {
        n -= 2;
    }

    return juce::jmax (
        2,
        n);
}

inline size_t indexOf (
    int channel,
    int sample,
    int samples)
{
    return static_cast<size_t> (channel)
         * static_cast<size_t> (samples)
         + static_cast<size_t> (sample);
}

int featureThreadCount (int count)
{
    if (count <= 0)
        return 0;

    const unsigned int hardwareThreads =
        juce::jmax (
            1u,
            std::thread::hardware_concurrency());

    return juce::jmax (
        1,
        juce::jmin (
            count,
            juce::jmin (
                maxFeatureThreads,
                static_cast<int> (hardwareThreads))));
}

template <typename Function>
void parallelFor (
    int count,
    Function&& function)
{
    if (count <= 0)
        return;

    const int threadCount =
        featureThreadCount (count);

    if (threadCount <= 1 || count < 8)
    {
        for (int i = 0; i < count; ++i)
            function (i);
        return;
    }

    std::atomic<int> next { 0 };
    std::vector<std::thread> workers;
    workers.reserve (
        static_cast<size_t> (threadCount));

    for (int worker = 0;
         worker < threadCount;
         ++worker)
    {
        workers.emplace_back (
            [&]()
            {
                for (;;)
                {
                    const int i =
                        next.fetch_add (
                            1,
                            std::memory_order_relaxed);

                    if (i >= count)
                        break;

                    function (i);
                }
            });
    }

    for (auto& worker : workers)
        worker.join();
}

template <typename Function>
void parallelForBlocks (
    int totalCount,
    int blockSize,
    Function&& function)
{
    if (totalCount <= 0)
        return;

    blockSize = juce::jmax (1, blockSize);

    const int blockCount =
        (totalCount + blockSize - 1)
        / blockSize;

    parallelFor (
        blockCount,
        [&] (int block)
        {
            const int begin =
                block * blockSize;

            const int end =
                juce::jmin (
                    totalCount,
                    begin + blockSize);

            function (begin, end);
        });
}

double bytesToMiB (size_t bytes)
{
    return static_cast<double> (bytes)
        / (1024.0 * 1024.0);
}

void logBufferSize (
    const char* name,
    const std::vector<float>& values)
{
    DBG (
        "Feature extraction "
        + juce::String (name)
        + " buffer: "
        + juce::String (
            bytesToMiB (
                values.size()
                * sizeof (float)),
            1)
        + " MiB");
}



struct SosSection
{
    double b0 = 1.0;
    double b1 = 0.0;
    double b2 = 0.0;
    double a1 = 0.0;
    double a2 = 0.0;
};

using SosFilter = std::vector<SosSection>;

SosFilter iblApSos30000()
{
    // scipy.signal.butter(
    //     3, 300, btype="highpass", fs=30000, output="sos")
    return {
        {
            0.9390916523119581,
            -0.9390916523119581,
            0.0,
            -0.9390625058174923,
            0.0
        },
        {
            1.0,
            -2.0,
            1.0,
            -1.9352943868599919,
            0.9391207988064240
        }
    };
}

SosFilter iblLfSos2500()
{
    // scipy.signal.butter(
    //     3, [0.5, 300], btype="bandpass", fs=2500, output="sos")
    return {
        {
            0.02851666410014575,
            0.05703332820029149,
            0.02851666410014575,
            -1.0875417986334095,
            0.4911339591905614
        },
        {
            1.0,
            0.0,
            -1.0,
            -1.4327717887942428,
            0.4334848244064954
        },
        {
            1.0,
            -2.0,
            1.0,
            -1.9987445728814415,
            0.9987461535244078
        }
    };
}

std::array<double, 2> lfilterSteadyStateZi (
    const SosSection& s)
{
    // Direct-form-II-transposed steady state for a unit step, equivalent to
    // scipy.signal.lfilter_zi for a normalized second-order section.
    const double denominator =
        1.0 + s.a1 + s.a2;

    const double gain =
        std::abs (denominator) > 1.0e-18
            ? (s.b0 + s.b1 + s.b2) / denominator
            : 0.0;

    return {
        gain - s.b0,
        s.b2 - s.a2 * gain
    };
}

std::vector<std::array<double, 2>>
sosfiltZi (const SosFilter& sos)
{
    std::vector<std::array<double, 2>> zi;
    zi.reserve (sos.size());

    double scale = 1.0;

    for (const auto& section : sos)
    {
        auto sectionZi =
            lfilterSteadyStateZi (section);

        sectionZi[0] *= scale;
        sectionZi[1] *= scale;
        zi.push_back (sectionZi);

        const double denominator =
            1.0 + section.a1 + section.a2;

        const double gain =
            std::abs (denominator) > 1.0e-18
                ? (section.b0
                   + section.b1
                   + section.b2)
                    / denominator
                : 0.0;

        scale *= gain;
    }

    return zi;
}

void sosFilterInPlace (
    std::vector<double>& signal,
    const SosFilter& sos,
    const std::vector<std::array<double, 2>>& zi,
    double initialValue)
{
    for (size_t sectionIndex = 0;
         sectionIndex < sos.size();
         ++sectionIndex)
    {
        const auto& s =
            sos[sectionIndex];

        double z1 =
            zi[sectionIndex][0]
            * initialValue;

        double z2 =
            zi[sectionIndex][1]
            * initialValue;

        for (auto& x : signal)
        {
            const double y =
                s.b0 * x + z1;

            const double nextZ1 =
                s.b1 * x
                - s.a1 * y
                + z2;

            const double nextZ2 =
                s.b2 * x
                - s.a2 * y;

            x = y;
            z1 = nextZ1;
            z2 = nextZ2;
        }
    }
}

void sosFiltFiltOneChannel (
    float* trace,
    int samples,
    const SosFilter& sos,
    int padLength)
{
    if (trace == nullptr
        || samples <= padLength + 1
        || sos.empty())
    {
        return;
    }

    padLength =
        juce::jlimit (
            0,
            samples - 2,
            padLength);

    std::vector<double> extended (
        static_cast<size_t> (
            samples + 2 * padLength));

    // scipy.signal.sosfiltfilt's default padtype is "odd":
    // left = 2*x[0] - reverse(x[1:pad+1])
    for (int i = 0; i < padLength; ++i)
    {
        extended[static_cast<size_t> (i)] =
            2.0 * static_cast<double> (trace[0])
            - static_cast<double> (
                trace[padLength - i]);
    }

    for (int i = 0; i < samples; ++i)
    {
        extended[
            static_cast<size_t> (
                padLength + i)] =
            trace[i];
    }

    // right = 2*x[-1] - reverse(x[-pad-1:-1])
    for (int i = 0; i < padLength; ++i)
    {
        extended[
            static_cast<size_t> (
                padLength + samples + i)] =
            2.0 * static_cast<double> (
                trace[samples - 1])
            - static_cast<double> (
                trace[samples - 2 - i]);
    }

    const auto zi =
        sosfiltZi (sos);

    sosFilterInPlace (
        extended,
        sos,
        zi,
        extended.front());

    std::reverse (
        extended.begin(),
        extended.end());

    sosFilterInPlace (
        extended,
        sos,
        zi,
        extended.front());

    std::reverse (
        extended.begin(),
        extended.end());

    for (int i = 0; i < samples; ++i)
    {
        trace[i] =
            static_cast<float> (
                extended[
                    static_cast<size_t> (
                        padLength + i)]);
    }
}

float np1SampleShiftFromElectrode (
    int electrodeIndex)
{
    // neuropixel.adc_shifts(version=1):
    // 32 ADCs x 12 channels, 13 conversion cycles.
    //
    // For each group of 24 hardware channels, even/odd channels are separate
    // ADCs and the serial position is floor((channel % 24) / 2).
    if (electrodeIndex < 0)
        return 0.0f;

    const int positionWithinAdc =
        (electrodeIndex % 24) / 2;


    return static_cast<float> (
        positionWithinAdc)
        / 13.0f;
}


float pythonCompatibleCorRatio (
    const std::vector<float>& ap,
    int channel,
    int channels,
    int samples,
    int sampleStride)
{
    // Exact port of ephysatlas.features.xcor_acor_ratio for the returned
    // (unmirrored) channels. Python mirrors 12 channels at each edge, evaluates
    // covariance diagonals -8..+8 INCLUDING the zero/autocovariance diagonal,
    // masks mirrored self-pairs on non-zero diagonals, then divides the
    // nan-mean by the autocovariance.
    constexpr int nMirror = 12;
    constexpr int nDiags = 8;

    if (channel < 0 || channel >= channels || channels <= 1 || samples <= 0)
        return std::numeric_limits<float>::quiet_NaN();

    sampleStride = juce::jmax (1, sampleStride);

    auto mirroredToSource =
        [channels] (int mirroredIndex) -> int
        {
            constexpr int m = 12;
            if (mirroredIndex < m)
                return m - mirroredIndex;

            const int middleEnd = m + channels;
            if (mirroredIndex < middleEnd)
                return mirroredIndex - m;

            // np.r_[v[n_mirror:0:-1], v, v[-2:-n_mirror-2:-1]]
            const int k = mirroredIndex - middleEnd;
            return channels - 2 - k;
        };

    const int mirroredCentre = nMirror + channel;
    const int centreSource = mirroredToSource (mirroredCentre);

    auto dotProduct =
        [&] (int a, int b) -> double
        {
            if (a < 0 || a >= channels || b < 0 || b >= channels)
                return std::numeric_limits<double>::quiet_NaN();

            const float* xa =
                ap.data()
                + static_cast<size_t> (a)
                    * static_cast<size_t> (samples);

            const float* xb =
                ap.data()
                + static_cast<size_t> (b)
                    * static_cast<size_t> (samples);

            double sum = 0.0;
            int count = 0;

            for (int sample = 0; sample < samples; sample += sampleStride)
            {
                const double va = xa[sample];
                const double vb = xb[sample];

                // NumPy matrix multiplication propagates non-finite values.
                if (! std::isfinite (va) || ! std::isfinite (vb))
                    return std::numeric_limits<double>::quiet_NaN();

                sum += va * vb;
                ++count;
            }

            return count > 0
                ? sum
                : std::numeric_limits<double>::quiet_NaN();
        };

    const double autoCov = dotProduct (centreSource, centreSource);
    if (! std::isfinite (autoCov) || std::abs (autoCov) <= 1.0e-30)
        return std::numeric_limits<float>::quiet_NaN();

    double sum = 0.0;
    int count = 0;

    for (int offset = -nDiags; offset <= nDiags; ++offset)
    {
        const int neighbourMirrored = mirroredCentre + offset;
        const int neighbourSource = mirroredToSource (neighbourMirrored);

        // Python explicitly sets repeated mirror/self pairs to NaN for
        // non-zero diagonals.
        if (offset != 0 && neighbourSource == centreSource)
            continue;

        const double covariance = dotProduct (centreSource, neighbourSource);
        if (! std::isfinite (covariance))
            continue;

        sum += covariance;
        ++count;
    }

    if (count <= 0)
        return std::numeric_limits<float>::quiet_NaN();

    return static_cast<float> (
        (sum / static_cast<double> (count)) / autoCov);
}

}

double FeatureExtractor::rms (
    const float* x,
    int n)
{
    if (x == nullptr || n <= 0)
    {
        return std::numeric_limits<double>::
            quiet_NaN();
    }

    double sumSquares = 0.0;

    for (int i = 0; i < n; ++i)
    {
        const double value =
            static_cast<double> (x[i]);

        sumSquares += value * value;
    }

    return std::sqrt (
        sumSquares
        / static_cast<double> (n));
}

void FeatureExtractor::scaleInPlace (
    std::vector<float>& data,
    float scale)
{
    if (data.empty() || scale == 1.0f)
        return;

    // Parallelize by contiguous blocks, not by individual samples.
    constexpr int scaleBlockSize = 1 << 18;

    parallelForBlocks (
        static_cast<int> (data.size()),
        scaleBlockSize,
        [&] (int begin, int end)
        {
            for (int i = begin;
                 i < end;
                 ++i)
            {
                data[static_cast<size_t> (i)]
                    *= scale;
            }
        });
}

void FeatureExtractor::commonMedianReference (
    std::vector<float>& data,
    int channels,
    int samples,
    const std::vector<int>& labels)
{
    // ibldsp car(): median over the channels inside the brain (label != 3), subtracted from
    // those channels only; channels outside the brain are left as they are.
    std::vector<int> inside;
    for (int channel = 0; channel < channels; ++channel)
        if (labels.empty() || labels[static_cast<size_t> (channel)] != 3)
            inside.push_back (channel);
    const int count = static_cast<int> (inside.size());
    if (count <= 1 || samples <= 0)
        return;

    parallelForBlocks (
        samples,
        medianReferenceBlockSamples,
        [&] (int beginSample, int endSample)
        {
            // One scratch allocation per block/worker invocation rather than
            // allocating a scratch vector independently for every time sample.
            std::vector<float> values (static_cast<size_t> (count));
            for (int sample = beginSample; sample < endSample; ++sample)
            {
                for (int k = 0; k < count; ++k)
                    values[static_cast<size_t> (k)] = data[indexOf (inside[static_cast<size_t> (k)], sample, samples)];
                const int middle = count / 2;
                std::nth_element (values.begin(), values.begin() + middle, values.end());
                float median = values[static_cast<size_t> (middle)];
                if ((count % 2) == 0)
                {
                    const float upper = median;
                    const float lower = *std::max_element (values.begin(), values.begin() + middle);
                    median = 0.5f * (lower + upper);
                }
                for (int channel : inside)
                    data[indexOf (channel, sample, samples)] -= median;
            }
        });
}

void FeatureExtractor::temporalHighpassZeroPhase (
    std::vector<float>& data,
    int channels,
    int samples,
    double sampleRate,
    double cutoffHz,
    int passes)
{
    if (samples < 3
        || sampleRate <= 0.0
        || cutoffHz <= 0.0)
    {
        return;
    }

    const double dt =
        1.0 / sampleRate;

    const double rc =
        1.0
        / (2.0
           * juce::MathConstants<double>::pi
           * cutoffHz);

    const float alpha =
        static_cast<float> (
            rc / (rc + dt));

    passes =
        juce::jmax (
            1,
            passes);

    parallelFor (
        channels,
        [&] (int channel)
        {
            float* trace =
                data.data()
                + static_cast<size_t> (channel)
                    * static_cast<size_t> (samples);

            auto filter =
                [&] (bool reversed)
                {
                    const int first =
                        reversed
                            ? samples - 1
                            : 0;

                    const int step =
                        reversed
                            ? -1
                            : 1;

                    const int end =
                        reversed
                            ? -1
                            : samples;

                    float previousInput =
                        trace[first];

                    float previousOutput =
                        0.0f;

                    trace[first] = 0.0f;

                    for (int i = first + step;
                         i != end;
                         i += step)
                    {
                        const float currentInput =
                            trace[i];

                        const float output =
                            alpha
                            * (previousOutput
                               + currentInput
                               - previousInput);

                        trace[i] = output;

                        previousInput =
                            currentInput;

                        previousOutput =
                            output;
                    }
                };

            for (int pass = 0;
                 pass < passes;
                 ++pass)
            {
                filter (false);
                filter (true);
            }
        });
}

void FeatureExtractor::temporalLowpassZeroPhase (
    std::vector<float>& data,
    int channels,
    int samples,
    double sampleRate,
    double cutoffHz,
    int passes)
{
    if (samples < 3
        || sampleRate <= 0.0
        || cutoffHz <= 0.0)
    {
        return;
    }

    const double dt =
        1.0 / sampleRate;

    const double rc =
        1.0
        / (2.0
           * juce::MathConstants<double>::pi
           * cutoffHz);

    const float alpha =
        static_cast<float> (
            dt / (rc + dt));

    passes =
        juce::jmax (
            1,
            passes);

    parallelFor (
        channels,
        [&] (int channel)
        {
            float* trace =
                data.data()
                + static_cast<size_t> (channel)
                    * static_cast<size_t> (samples);

            auto filter =
                [&] (bool reversed)
                {
                    const int first =
                        reversed
                            ? samples - 1
                            : 0;

                    const int step =
                        reversed
                            ? -1
                            : 1;

                    const int end =
                        reversed
                            ? -1
                            : samples;

                    float previousOutput =
                        trace[first];

                    for (int i = first + step;
                         i != end;
                         i += step)
                    {
                        previousOutput =
                            previousOutput
                            + alpha
                                * (trace[i]
                                   - previousOutput);

                        trace[i] =
                            previousOutput;
                    }
                };

            for (int pass = 0;
                 pass < passes;
                 ++pass)
            {
                filter (false);
                filter (true);
            }
        });
}


void FeatureExtractor::iblApButterworthZeroPhase (
    std::vector<float>& data,
    int channels,
    int samples,
    double sampleRate)
{
    if (data.empty()
        || channels <= 0
        || samples <= 0)
    {
        return;
    }

    // Exact SciPy coefficients for the IBL AP rate.
    if (std::abs (sampleRate - 30000.0) > 1.0)
    {
        DBG (
            "IBL AP Butterworth: non-standard sample rate "
            + juce::String (sampleRate, 3)
            + " Hz; using legacy high-pass fallback.");

        temporalHighpassZeroPhase (
            data,
            channels,
            samples,
            sampleRate,
            300.0,
            1);

        return;
    }

    const auto sos =
        iblApSos30000();

    // SciPy default padlen:
    // 3 * (2*n_sections + 1 - min(zero_b2, zero_a2)) = 12.
    constexpr int padLength = 12;

    parallelFor (
        channels,
        [&] (int channel)
        {
            float* trace =
                data.data()
                + static_cast<size_t> (channel)
                    * static_cast<size_t> (samples);

            sosFiltFiltOneChannel (
                trace,
                samples,
                sos,
                padLength);
        });
}

void FeatureExtractor::iblLfButterworthZeroPhase (
    std::vector<float>& data,
    int channels,
    int samples,
    double sampleRate)
{
    if (data.empty()
        || channels <= 0
        || samples <= 0)
    {
        return;
    }

    // Exact SciPy coefficients for the standard IBL LF rate.
    if (std::abs (sampleRate - 2500.0) > 1.0)
    {
        DBG (
            "IBL LF Butterworth: non-standard sample rate "
            + juce::String (sampleRate, 3)
            + " Hz; using legacy 0.5-300 Hz fallback. "
              "This run is not filter-equivalent to IBL.");

        temporalHighpassZeroPhase (
            data,
            channels,
            samples,
            sampleRate,
            0.5,
            1);

        temporalLowpassZeroPhase (
            data,
            channels,
            samples,
            sampleRate,
            300.0,
            1);

        return;
    }

    const auto sos =
        iblLfSos2500();

    // SciPy default padlen for the 3-section band-pass = 21.
    constexpr int padLength = 21;

    parallelFor (
        channels,
        [&] (int channel)
        {
            float* trace =
                data.data()
                + static_cast<size_t> (channel)
                    * static_cast<size_t> (samples);

            sosFiltFiltOneChannel (
                trace,
                samples,
                sos,
                padLength);
        });
}

void FeatureExtractor::applyNp1SampleShift (
    std::vector<float>& data,
    int channels,
    int samples,
    const std::vector<int>& electrodeIndices)
{
    if (data.empty()
        || channels <= 0
        || samples <= 1)
    {
        return;
    }

    parallelFor (
        channels,
        [&] (int channel)
        {
            const int electrodeIndex =
                channel
                    < static_cast<int> (
                        electrodeIndices.size())
                ? electrodeIndices[
                    static_cast<size_t> (
                        channel)]
                : channel;

            const float shift =
                np1SampleShiftFromElectrode (
                    electrodeIndex);

            if (std::abs (shift) < 1.0e-8f)
                return;

            float* trace =
                data.data()
                + static_cast<size_t> (channel)
                    * static_cast<size_t> (samples);

            // kiss_fftr requires an even transform length. A nominal 60-s LF
            // window can contain 150001 samples when fs is ~2500.02 Hz.
            const int nfft =
                efficientFftLengthAtMost (
                    samples);

            if (nfft <= 1)
                return;

            if (nfft != samples)
            {
                DBG (
                    "NP1 sample shift: using fast FFT length "
                    + juce::String (nfft)
                    + " instead of "
                    + juce::String (samples)
                    + " samples.");
            }

            const int bins =
                nfft / 2 + 1;

            std::vector<kiss_fft_scalar> fftInput (
                static_cast<size_t> (nfft));

            std::vector<kiss_fft_cpx> spectrum (
                static_cast<size_t> (bins));

            for (int i = 0; i < nfft; ++i)
            {
                fftInput[static_cast<size_t> (i)] =
                    static_cast<kiss_fft_scalar> (
                        trace[i]);
            }

            kiss_fftr_cfg forward =
                kiss_fftr_alloc (
                    nfft,
                    0,
                    nullptr,
                    nullptr);

            kiss_fftr_cfg inverse =
                kiss_fftr_alloc (
                    nfft,
                    1,
                    nullptr,
                    nullptr);

            if (forward == nullptr
                || inverse == nullptr)
            {
                if (forward != nullptr)
                    std::free (forward);

                if (inverse != nullptr)
                    std::free (inverse);

                return;
            }

            kiss_fftr (
                forward,
                fftInput.data(),
                spectrum.data());

            // Match ibldsp.fourier.fshift:
            // W *= exp(i * angle(rfft(delta[n-1])) * shift)
            // angle(delta shifted by +1) = -2*pi*k/N.
            for (int k = 0; k < bins; ++k)
            {
                const double phase =
                    -2.0
                    * juce::MathConstants<double>::pi
                    * static_cast<double> (k)
                    * static_cast<double> (shift)
                    / static_cast<double> (nfft);

                const double c =
                    std::cos (phase);

                const double s =
                    std::sin (phase);

                const double real =
                    spectrum[
                        static_cast<size_t> (k)].r;

                const double imag =
                    spectrum[
                        static_cast<size_t> (k)].i;

                spectrum[
                    static_cast<size_t> (k)].r =
                    static_cast<kiss_fft_scalar> (
                        real * c - imag * s);

                spectrum[
                    static_cast<size_t> (k)].i =
                    static_cast<kiss_fft_scalar> (
                        real * s + imag * c);
            }

            kiss_fftri (
                inverse,
                spectrum.data(),
                fftInput.data());

            const float normalization =
                1.0f
                / static_cast<float> (nfft);

            for (int i = 0; i < nfft; ++i)
            {
                trace[i] =
                    static_cast<float> (
                        fftInput[
                            static_cast<size_t> (i)])
                    * normalization;
            }

            std::free (forward);
            std::free (inverse);
        });
}

double FeatureExtractor::cosineBandWeight (
    double frequency,
    double low,
    double high)
{
    const double width =
        high - low;

    const double centre =
        0.5 * (low + high);

    const double distance =
        std::abs (
            frequency - centre);

    if (width <= 0.0
        || distance >= width)
    {
        return 0.0;
    }

    return 0.5
        * (1.0
           + std::cos (
               juce::MathConstants<double>::pi
               * distance
               / width));
}

std::vector<double>
FeatureExtractor::oneSidedPeriodogram (
    const float* x,
    int n,
    double sampleRate,
    std::vector<double>& frequencies)
{
    frequencies.clear();

    if (x == nullptr
        || n <= 1
        || sampleRate <= 0.0)
    {
        return {};
    }

    // scipy.signal.periodogram uses the complete input length. Open Ephys
    // localization windows are even at the canonical AP/LF rates. If a
    // non-standard capture supplies an odd length, remove only the last sample
    // because kiss_fftr requires an even transform length.
    const int nfft =
        (n % 2 == 0) ? n : n - 1;

    if (nfft <= 1)
        return {};

    const int bins = nfft / 2 + 1;

    std::vector<kiss_fft_scalar> input (
        static_cast<size_t> (nfft),
        0.0f);

    std::vector<kiss_fft_cpx> output (
        static_cast<size_t> (bins));

    double mean = 0.0;

    for (int i = 0; i < nfft; ++i)
        mean += x[i];

    mean /=
        static_cast<double> (nfft);

    for (int i = 0; i < nfft; ++i)
    {
        input[static_cast<size_t> (i)] =
            static_cast<kiss_fft_scalar> (
                static_cast<double> (x[i])
                - mean);
    }

    kiss_fftr_cfg fftConfig =
        kiss_fftr_alloc (
            nfft,
            0,
            nullptr,
            nullptr);

    if (fftConfig == nullptr)
        return {};

    kiss_fftr (
        fftConfig,
        input.data(),
        output.data());

    std::free (fftConfig);

    frequencies.resize (
        static_cast<size_t> (bins));

    std::vector<double> periodogram (
        static_cast<size_t> (bins),
        0.0);

    const double scale =
        1.0
        / (sampleRate
           * static_cast<double> (nfft));

    for (int bin = 0;
         bin < bins;
         ++bin)
    {
        const double real =
            output[static_cast<size_t> (bin)].r;

        const double imaginary =
            output[static_cast<size_t> (bin)].i;

        double value =
            (real * real
             + imaginary * imaginary)
            * scale;

        const bool isDc =
            bin == 0;

        const bool isNyquist =
            (nfft % 2 == 0)
            && bin == nfft / 2;

        if (! isDc && ! isNyquist)
            value *= 2.0;

        periodogram[static_cast<size_t> (bin)] =
            value;

        frequencies[static_cast<size_t> (bin)] =
            sampleRate
            * static_cast<double> (bin)
            / static_cast<double> (nfft);
    }

    return periodogram;
}

double FeatureExtractor::bandPowerDb (
    const std::vector<double>& frequencies,
    const std::vector<double>& periodogram,
    double low,
    double high)
{
    double weightedPower = 0.0;
    double weightSum = 0.0;

    for (size_t i = 0;
         i < frequencies.size()
             && i < periodogram.size();
         ++i)
    {
        const double weight =
            cosineBandWeight (
                frequencies[i],
                low,
                high);

        weightedPower +=
            periodogram[i] * weight;

        weightSum += weight;
    }

    return 10.0
        * std::log10 (
            std::max (
                weightedPower
                / std::max (
                    weightSum,
                    eps),
                eps));
}

FeatureMatrix FeatureExtractor::extract (
    LocalizationInput& input,
    const LocalizationConfig& config,
    FeatureProgressCallback progressCallback) const
{
    auto report =
        [&] (
            float progress,
            const juce::String& status)
        {
            if (progressCallback)
            {
                progressCallback (
                    juce::jlimit (
                        0.0f,
                        1.0f,
                        progress),
                    status);
            }

            DBG (
                "Feature extraction: "
                + status);
        };

    const int channels =
        juce::jmax (
            input.ap.numChannels,
            input.lf.numChannels);

    FeatureMatrix output;
    output.rows = channels;
    // Exact first-30 order from ephysatlas.spatial_encoder.utils.FEATURE_LIST.
    // alpha_mean/alpha_std/spike_count require spike detection/localization and
    // are intentionally left NaN for now.
    output.names = {
        "rms_lf",
        "psd_lfp",
        "psd_alpha",
        "psd_beta",
        "psd_gamma",
        "psd_delta",
        "psd_theta",
        "psd_lfp_bpr",
        "psd_alpha_bpr",
        "psd_beta_bpr",
        "psd_gamma_bpr",
        "psd_delta_bpr",
        "psd_theta_bpr",
        "rms_lf_bpr",
        "psd_residual_lfp",
        "psd_residual_alpha",
        "psd_residual_beta",
        "psd_residual_gamma",
        "psd_residual_delta",
        "psd_residual_theta",
        "decay_fit_error",
        "decay_fit_r_squared",
        "decay_n_peaks",
        "aperiodic_exponent",
        "aperiodic_offset",
        "cor_ratio",
        "rms_ap",
        "alpha_mean",
        "alpha_std",
        "spike_count"
    };

    output.cols =
        static_cast<int> (
            output.names.size());

    output.values.assign (
        static_cast<size_t> (
            output.rows
            * output.cols),
        std::numeric_limits<float>::
            quiet_NaN());

    report (
        0.00f,
        "Preparing captured buffers");

    // ---------------------------------------------------------------
    // AP
    // ---------------------------------------------------------------
    if (! input.ap.channelMajorSamples.empty())
    {
        auto& ap =
            input.ap.channelMajorSamples;

        logBufferSize ("AP", ap);

        report (
            0.04f,
            "AP: scaling");

        scaleInPlace (
            ap,
            static_cast<float> (
                config.inputScaleToVolts));

        report (
            0.10f,
            "AP: IBL 3rd-order Butterworth high-pass");

        // ibldsp.voltage.detect_bad_channels on the raw snippet; the labels drive the
        // interpolation and CAR of both bands, as in ephysatlas.realtime.snippet_features.
        const int labelThreads = juce::jmax (1, static_cast<int> (std::thread::hardware_concurrency()) - 2);
        const auto badChannelStats = Realtime::badChannelRawStats (
            ap.data(), input.ap.numChannels, input.ap.numSamples, input.ap.sampleRate, labelThreads);

        iblApButterworthZeroPhase (
            ap,
            input.ap.numChannels,
            input.ap.numSamples,
            input.ap.sampleRate);

        output.channelLabels = Realtime::badChannelLabels (
            badChannelStats, ap.data(), input.ap.numChannels, input.ap.numSamples, labelThreads);

        report (
            0.24f,
            "AP: NP1 ADC fractional sample-shift correction");

        applyNp1SampleShift (
            ap,
            input.ap.numChannels,
            input.ap.numSamples,
            input.ap.electrodeIndices);

        report (
            0.34f,
            "AP: IBL median common reference");

        {
            const auto geometry = Realtime::Np1Geometry::make (input.ap.numChannels);
            Realtime::interpolateBadChannels (
                ap.data(), input.ap.numChannels, input.ap.numSamples, output.channelLabels, geometry.x, geometry.y);
        }

        commonMedianReference (
            ap,
            input.ap.numChannels,
            input.ap.numSamples,
            output.channelLabels);

        report (
            0.55f,
            "AP: computing RMS");

        parallelFor (
            input.ap.numChannels,
            [&] (int sourceChannel)
            {
                const int outputChannel =
                    config.reverseChannelOrder
                        ? channels
                            - 1
                            - sourceChannel
                        : sourceChannel;

                if (outputChannel < 0
                    || outputChannel >= channels)
                {
                    return;
                }

                const float* x =
                    ap.data()
                    + static_cast<size_t> (
                        sourceChannel)
                        * static_cast<size_t> (
                            input.ap.numSamples);

                const double channelRms =
                    rms (
                        x,
                        input.ap.numSamples);

                const size_t row =
                    static_cast<size_t> (
                        outputChannel
                        * output.cols);

                // FEATURE_LIST[26] = rms_ap
                output.values[row + 26] =
                    static_cast<float> (
                        20.0
                        * std::log10 (
                            std::max (
                                channelRms,
                                eps)));

                // FEATURE_LIST[25] = cor_ratio.
                // This keeps the Python +/-8-neighbor covariance definition,
                // evaluated on a deterministic temporal subsample for speed.
                output.values[row + 25] =
                    pythonCompatibleCorRatio (
                        ap,
                        sourceChannel,
                        input.ap.numChannels,
                        input.ap.numSamples,
                        config.corRatioSampleStride);
            });

        // AP is no longer needed after rms_ap has been produced. Releasing it
        // before LF preprocessing/FFTs substantially reduces peak resident RAM.
        report (
            0.60f,
            "AP: releasing captured buffer");

        std::vector<float>().swap (
            input.ap.channelMajorSamples);
    }

    // ---------------------------------------------------------------
    // LF
    // ---------------------------------------------------------------
    if (! input.lf.channelMajorSamples.empty())
    {
        auto& lf =
            input.lf.channelMajorSamples;

        // File Reader can occasionally expose a valid LF buffer while
        // reporting sampleRate == 0. RMS does not depend on sample rate, but
        // the Butterworth filter and periodogram do. Use the standard NP1/IBL
        // 2500-Hz rate only for this non-positive metadata failure.
        const double configuredLfFallback =
            config.lfFallbackSampleRateHz > 0.0
                ? config.lfFallbackSampleRateHz
                : 2500.0;

        const double reportedLfSampleRate =
            std::isfinite (
                input.lf.reportedSampleRate)
                ? input.lf.reportedSampleRate
                : input.lf.sampleRate;

        const double lfSampleRate =
            input.lf.sampleRate > 0.0
                ? input.lf.sampleRate
                : configuredLfFallback;

        const bool lfFallbackUsed =
            reportedLfSampleRate <= 0.0
            && std::abs (
                lfSampleRate
                - configuredLfFallback)
                < 1.0e-6;


        if (reportedLfSampleRate <= 0.0)
        {
            DBG (
                "LF reported sample rate is "
                + juce::String (reportedLfSampleRate, 6)
                + " Hz for stream "
                + input.lf.streamKey
                + "; using fallback "
                + juce::String (lfSampleRate, 3)
                + " Hz for filtering and PSD.");

            report (
                0.615f,
                "LF: invalid reported sample rate; using "
                    + juce::String (lfSampleRate, 0)
                    + " Hz fallback");
        }

        logBufferSize ("LF", lf);

        report (
            0.63f,
            "LF: scaling");

        scaleInPlace (
            lf,
            static_cast<float> (
                config.inputScaleToVolts));

        report (
            0.67f,
            "LF: IBL 3rd-order 0.5-300 Hz Butterworth band-pass");

        iblLfButterworthZeroPhase (
            lf,
            input.lf.numChannels,
            input.lf.numSamples,
            lfSampleRate);

        report (
            0.73f,
            "LF: NP1 ADC fractional sample-shift correction");

        applyNp1SampleShift (
            lf,
            input.lf.numChannels,
            input.lf.numSamples,
            input.lf.electrodeIndices);

        // IBL would interpolate channels labelled dead/noisy here, before CAR.
        // Probe Viewer currently receives no IBL-compatible channel_labels in
        // BandWindow, so we deliberately do not invent a different detector.
        // This is the only substantive destripe_lfp preprocessing step that
        // cannot be reproduced from the metadata currently available here.
        report (
            0.78f,
            "LF: IBL median common reference");

        std::vector<int> lfLabels = output.channelLabels;
        lfLabels.resize (static_cast<size_t> (input.lf.numChannels), 0);
        {
            const auto geometry = Realtime::Np1Geometry::make (input.lf.numChannels);
            Realtime::interpolateBadChannels (
                lf.data(), input.lf.numChannels, input.lf.numSamples, lfLabels, geometry.x, geometry.y);
        }

        commonMedianReference (
            lf,
            input.lf.numChannels,
            input.lf.numSamples,
            lfLabels);

        report (
            0.84f,
            "LF: RMS and PSD features");

        // This preserves the original full-window periodogram definition.
        // Parallelism is capped to avoid allocating too many simultaneous FFT
        // workspaces and to keep the machine responsive.
        parallelFor (
            input.lf.numChannels,
            [&] (int sourceChannel)
            {
                const int outputChannel =
                    config.reverseChannelOrder
                        ? channels
                            - 1
                            - sourceChannel
                        : sourceChannel;

                if (outputChannel < 0
                    || outputChannel >= channels)
                {
                    return;
                }

                const size_t row =
                    static_cast<size_t> (
                        outputChannel
                        * output.cols);

                const float* x =
                    lf.data()
                    + static_cast<size_t> (
                        sourceChannel)
                        * static_cast<size_t> (
                            input.lf.numSamples);

                const double channelRms =
                    rms (
                        x,
                        input.lf.numSamples);

                // FEATURE_LIST[0] = rms_lf
                output.values[row + 0] =
                    static_cast<float> (
                        20.0
                        * std::log10 (
                            std::max (
                                channelRms,
                                eps)));

                std::vector<double> frequencies;

                const auto periodogram =
                    oneSidedPeriodogram (
                        x,
                        input.lf.numSamples,
                        lfSampleRate,
                        frequencies);

                // Base LF PSD features in exact FEATURE_LIST order.
                output.values[row + 1] =
                    static_cast<float> (
                        bandPowerDb (
                            frequencies,
                            periodogram,
                            0.0,
                            90.0));

                output.values[row + 2] =
                    static_cast<float> (
                        bandPowerDb (
                            frequencies,
                            periodogram,
                            8.0,
                            12.0));

                output.values[row + 3] =
                    static_cast<float> (
                        bandPowerDb (
                            frequencies,
                            periodogram,
                            15.0,
                            30.0));

                output.values[row + 4] =
                    static_cast<float> (
                        bandPowerDb (
                            frequencies,
                            periodogram,
                            30.0,
                            90.0));

                output.values[row + 5] =
                    static_cast<float> (
                        bandPowerDb (
                            frequencies,
                            periodogram,
                            0.0,
                            4.0));

                output.values[row + 6] =
                    static_cast<float> (
                        bandPowerDb (
                            frequencies,
                            periodogram,
                            4.0,
                            10.0));

                // Spectral-decay / residual features: specparam on the Welch PSD as
                // ephysatlas.features.get_psd_decay_features (Realtime::psdDecayFeatures).
                // Columns 14..19: psd_residual_lfp, alpha, beta, gamma, delta, theta;
                // decay.residualDb order: delta, theta, alpha, beta, gamma, lfp.
                const std::vector<double> xDouble (x, x + input.lf.numSamples);
                const auto decay =
                    Realtime::psdDecayFeatures (
                        xDouble.data(),
                        input.lf.numSamples,
                        lfSampleRate,
                        frequencies,
                        periodogram);
                constexpr int residualBandOfColumn[6] = { 5, 2, 3, 4, 0, 1 };
                for (int k = 0; k < 6; ++k)
                    output.values[row + 14 + static_cast<size_t> (k)] =
                        static_cast<float> (decay.residualDb[static_cast<size_t> (residualBandOfColumn[k])]);
                output.values[row + 20] = static_cast<float> (decay.fitError);
                output.values[row + 21] = static_cast<float> (decay.rSquared);
                output.values[row + 22] = static_cast<float> (decay.nPeaks);
                output.values[row + 23] = static_cast<float> (decay.exponent);
                output.values[row + 24] = static_cast<float> (decay.offset);

                // [27:30] alpha_mean, alpha_std, spike_count stay NaN because
                // the reference pipeline obtains them from Dartsort spike
                // detection/localization and waveform extraction.
            });

        // BPR family exactly as ephysatlas.features.csd: FIR decimation by 10, Cadzow denoising
        // across channels, per-column first spatial difference (sigma / dx), then RMS and
        // periodogram band powers at fs / 10 (Realtime::bprFeatures).
        report (
            0.90f,
            "LF: BPR (decimation, Cadzow denoising, first spatial difference)");
        {
            const int nc = input.lf.numChannels;
            const int ns = input.lf.numSamples;
            std::vector<double> lfDouble (lf.begin(), lf.end());
            const auto geometry = Realtime::Np1Geometry::make (nc);
            const auto csd = Realtime::bprFeatures (
                lfDouble, nc, ns, lfSampleRate, geometry,
                juce::jmax (1, static_cast<int> (std::thread::hardware_concurrency()) - 2));
            // FeatureExtractor columns 7..13: psd_lfp, alpha, beta, gamma, delta, theta (bpr),
            // rms_lf_bpr; csd.bandDb order: delta, theta, alpha, beta, gamma, lfp.
            constexpr int bandOfColumn[6] = { 5, 2, 3, 4, 0, 1 };
            for (int sourceChannel = 0; sourceChannel < nc; ++sourceChannel)
            {
                const int outputChannel = config.reverseChannelOrder ? channels - 1 - sourceChannel : sourceChannel;
                if (outputChannel < 0 || outputChannel >= channels)
                    continue;
                const size_t row = static_cast<size_t> (outputChannel * output.cols);
                for (int k = 0; k < 6; ++k)
                    output.values[row + 7 + static_cast<size_t> (k)] =
                        static_cast<float> (csd.bandDb[static_cast<size_t> (bandOfColumn[k])][static_cast<size_t> (sourceChannel)]);
                output.values[row + 13] = static_cast<float> (csd.rmsDb[static_cast<size_t> (sourceChannel)]);
            }
        }

        report (
            0.98f,
            "LF: releasing captured buffer");

        std::vector<float>().swap (
            input.lf.channelMajorSamples);
    }

    report (
        1.0f,
        "Feature extraction complete");

    return output;
}


}
