#pragma once

// Port of ibldsp.voltage.detect_bad_channels (AP band, default thresholds) and
// ibldsp.voltage.interpolate_bad_channels, as used by ephysatlas.realtime.snippet_features:
// labels are detected on the raw AP snippet and used to destripe both bands.
// Labels: 0 good, 1 dead, 2 noisy, 3 outside the brain.

#include <vector>

namespace ProbeViewer::Localization::Realtime
{
// Statistics of the raw (unfiltered) AP snippet, taken before the destripe high-pass so the
// high-pass can then run in place.
struct BadChannelRawStats
{
    std::vector<double> xcor;  // similarity of each channel with the median reference
    std::vector<double> psdHf; // mean Welch PSD above 80 % of Nyquist (uV^2 / Hz)
};

// raw: [nc x ns] volts, channel-major.
BadChannelRawStats badChannelRawStats (const float* raw, int nc, int ns, double fs, int threads);

// hf: the same snippet after scipy.signal.sosfiltfilt(butter(3, 300 Hz, highpass)).
std::vector<int> badChannelLabels (const BadChannelRawStats& rawStats, const float* hf, int nc, int ns, int threads);

// Replaces dead and noisy channels (labels 1, 2) by the kriging average of good neighbours
// (p = 1.3, 20 um), in place. x, y: channel positions (um).
void interpolateBadChannels (float* data, int nc, int ns, const std::vector<int>& labels,
                             const std::vector<double>& x, const std::vector<double>& y);
} // namespace ProbeViewer::Localization::Realtime
