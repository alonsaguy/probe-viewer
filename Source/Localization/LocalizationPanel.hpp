#pragma once
#include <JuceHeader.h>

namespace ProbeViewer { class ProbeViewerNode; }
namespace ProbeViewer::Localization
{
class LocalizationPanel : public juce::Component, private juce::Timer
{
public:
    explicit LocalizationPanel (ProbeViewerNode&);
    void resized() override;
private:
    void timerCallback() override;
    ProbeViewerNode& node;
    juce::ToggleButton enable { "Enable" };
    juce::TextButton localize { "Localize" }, reset { "Reset" };
    juce::ComboBox featureSource;
    juce::TextEditor modelName;
    juce::Slider windowSeconds;
    juce::Label status;
};
}
