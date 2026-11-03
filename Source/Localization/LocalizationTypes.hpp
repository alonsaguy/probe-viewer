#pragma once

#include <JuceHeader.h>

namespace ProbeViewer::Localization
{
/** A position in IBL-Allen coordinates (origin bregma), micrometres: x ML, y AP, z DV. */
struct CcfPoint
{
    float apUm = 0.0f;
    float mlUm = 0.0f;
    float dvUm = 0.0f;
};

/** Stable physical-probe identity carried through the localization pipeline.
    The current implementation still localizes one probe at a time, but all
    results are keyed by this identity so that a queued/multi-probe worker can
    be added later without changing the result/UI contracts. */
struct ProbeIdentity
{
    juce::String id;            // Prefer hardware serial number; otherwise a stable device/source key.
    juce::String displayName;   // Open Ephys / device-facing name (e.g. ProbeA).
    juce::String deviceName;
    juce::String serialNumber;
    juce::String manufacturer;
    juce::String deviceIdentifier;
    int displayIndex = -1;      // Order of first appearance among physical probe devices.

    bool isValid() const noexcept { return id.isNotEmpty(); }
};
}
