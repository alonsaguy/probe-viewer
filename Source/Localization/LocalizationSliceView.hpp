#pragma once

#include "LocalizationResult.hpp"

#include <JuceHeader.h>
#include <map>

namespace ProbeViewer
{
class LocalizationSliceView : public juce::Component
{
public:
    LocalizationSliceView();
    ~LocalizationSliceView() override;

    void setLocalizationResult (
        Localization::LocalizationResultPtr result);

    void clear();
    void paint (juce::Graphics& g) override;

private:
    static juce::String resultKey (
        const Localization::LocalizationResult& result);

    std::map<
        juce::String,
        Localization::LocalizationResultPtr>
        resultsByProbe;

    Localization::LocalizationResultPtr referenceResult;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (
        LocalizationSliceView)
};
} // namespace ProbeViewer
