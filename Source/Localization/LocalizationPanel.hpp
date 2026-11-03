#pragma once
#include <JuceHeader.h>

#include <memory>

namespace ProbeViewer { class ProbeViewerNode; }
namespace ProbeViewer::Localization
{
// Editor controls of the real-time localization: start (continuous snippets), reset, the planned
// trajectory (ephysatlas-plan/1 JSON) and the localizer bundle directory, and the status.
class LocalizationPanel : public juce::Component, private juce::Timer
{
public:
    explicit LocalizationPanel (ProbeViewerNode&);
    void resized() override;
    void paint (juce::Graphics&) override;

private:
    void timerCallback() override;
    void choosePlan();
    void chooseBundle();

    ProbeViewerNode& node;
    juce::TextButton localize { "Localize" }, reset { "Reset" }, planButton { "Plan..." }, bundleButton { "Bundle..." };
    juce::Label planLabel, bundleLabel, status;
    std::unique_ptr<juce::FileChooser> chooser;
};
}
