#pragma once
#include <JuceHeader.h>

namespace ProbeViewer::Localization
{
class LocalizationModel
{
public:
    bool load (const juce::File&) { return false; }
    bool isLoaded() const noexcept { return false; }
    juce::String description() const { return "Not implemented"; }
};
}
