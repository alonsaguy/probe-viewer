#include "LocalizationPanel.hpp"
#include "../ProbeViewerNode.h"

namespace ProbeViewer::Localization
{
LocalizationPanel::LocalizationPanel (ProbeViewerNode& n)
    : node (n)
{
    enable.setButtonText ("Enable localization");
    localize.setButtonText ("Localize");
    reset.setButtonText ("Reset");

    addAndMakeVisible (enable);
    addAndMakeVisible (localize);
    addAndMakeVisible (reset);
    addAndMakeVisible (featureSource);
    addAndMakeVisible (modelName);
    addAndMakeVisible (windowSeconds);
    addAndMakeVisible (status);

    featureSource.addItem ("Probe Viewer features", 1);
    featureSource.addItem ("Raw-window RMS", 2);
    featureSource.setSelectedId (1, juce::dontSendNotification);

    modelName.setText ("random-ccf-demo", false);

    windowSeconds.setSliderStyle (juce::Slider::LinearBar);
    windowSeconds.setTextBoxStyle (juce::Slider::TextBoxRight, false, 58, 20);
    windowSeconds.setRange (1.0, 300.0, 1.0);
    windowSeconds.setValue (100.0, juce::dontSendNotification);
    windowSeconds.setTextValueSuffix (" s");

    status.setJustificationType (juce::Justification::centredLeft);
    status.setText ("Disabled", juce::dontSendNotification);

    enable.onClick = [this]
    {
        node.setLocalizationEnabled (enable.getToggleState());
    };

    localize.onClick = [this]
    {
        node.setLocalizationModelName (modelName.getText());
        node.setLocalizationFeatureSource (
            featureSource.getSelectedId() == 1
                ? FeatureSource::probeViewerFeatures
                : FeatureSource::rawWindow);
        node.setLocalizationWindowSeconds (windowSeconds.getValue());
        node.requestLocalization();
    };

    reset.onClick = [this]
    {
        node.resetLocalization();
    };

    startTimerHz (5);
}

void LocalizationPanel::resized()
{
    auto area = getLocalBounds().reduced (2);

    auto topRow = area.removeFromTop (22);
    enable.setBounds (topRow.removeFromLeft (132));
    status.setBounds (topRow);

    auto sourceRow = area.removeFromTop (24);
    featureSource.setBounds (sourceRow.removeFromLeft (136));
    windowSeconds.setBounds (sourceRow.reduced (2, 0));

    modelName.setBounds (area.removeFromTop (24).reduced (0, 1));

    auto buttonRow = area.removeFromTop (26);
    localize.setBounds (buttonRow.removeFromLeft (buttonRow.getWidth() / 2).reduced (1));
    reset.setBounds (buttonRow.reduced (1));
}

void LocalizationPanel::timerCallback()
{
    status.setText (node.getLocalizationStatus()
                        + "  "
                        + juce::String (node.getLocalizationProgress() * 100.0f, 0)
                        + "%",
                    juce::dontSendNotification);
}
} // namespace ProbeViewer::Localization
