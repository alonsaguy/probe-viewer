#include "LocalizationPanel.hpp"

#include "../ProbeViewerNode.h"

namespace ProbeViewer::Localization
{
LocalizationPanel::LocalizationPanel (ProbeViewerNode& n) : node (n)
{
    for (auto* b : { &localize, &reset, &planButton, &bundleButton })
        addAndMakeVisible (b);
    for (auto* l : { &planLabel, &bundleLabel, &status })
    {
        l->setFont (juce::FontOptions (12.0f));
        l->setColour (juce::Label::textColourId, juce::Colours::black);
        l->setMinimumHorizontalScale (0.7f);
        addAndMakeVisible (l);
    }
    status.setJustificationType (juce::Justification::topLeft);
    localize.setTooltip ("Localize the probe continuously: a 5 s AP + LF snippet whenever the last one is done");
    reset.setTooltip ("Forget the snippets and the result");
    planButton.setTooltip ("Planned trajectory: an ephysatlas-plan/1 JSON file (ephysatlas.realtime.plan)");
    bundleButton.setTooltip ("Localizer bundle directory (ibleatools scripts/export_localizer_bundle.py)");
    localize.onClick = [this] { node.requestLocalization(); };
    reset.onClick = [this] { node.resetLocalization(); };
    planButton.onClick = [this] { choosePlan(); };
    bundleButton.onClick = [this] { chooseBundle(); };
    startTimerHz (4);
    timerCallback();
}

void LocalizationPanel::paint (juce::Graphics& g)
{
    g.setColour (juce::Colours::black.withAlpha (0.08f));
    g.fillRoundedRectangle (getLocalBounds().toFloat(), 4.0f);
}

void LocalizationPanel::resized()
{
    auto area = getLocalBounds().reduced (4);
    auto row = area.removeFromTop (24);
    localize.setBounds (row.removeFromLeft (80));
    row.removeFromLeft (4);
    reset.setBounds (row.removeFromLeft (60));
    area.removeFromTop (4);
    row = area.removeFromTop (22);
    planButton.setBounds (row.removeFromLeft (70));
    planLabel.setBounds (row.withTrimmedLeft (4));
    area.removeFromTop (2);
    row = area.removeFromTop (22);
    bundleButton.setBounds (row.removeFromLeft (70));
    bundleLabel.setBounds (row.withTrimmedLeft (4));
    area.removeFromTop (2);
    status.setBounds (area);
}

void LocalizationPanel::choosePlan()
{
    chooser = std::make_unique<juce::FileChooser> ("Planned trajectory (ephysatlas-plan/1)", juce::File(), "*.json");
    chooser->launchAsync (juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectFiles,
                          [this] (const juce::FileChooser& fc)
                          {
                              const auto f = fc.getResult();
                              if (f.existsAsFile())
                                  node.setLocalizationPlannedTrajectoryFile (f.getFullPathName());
                          });
}

void LocalizationPanel::chooseBundle()
{
    chooser = std::make_unique<juce::FileChooser> ("Localizer bundle directory", juce::File());
    chooser->launchAsync (juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectDirectories,
                          [this] (const juce::FileChooser& fc)
                          {
                              const auto f = fc.getResult();
                              if (f.isDirectory())
                                  node.setLocalizationBundleDirectory (f.getFullPathName());
                          });
}

void LocalizationPanel::timerCallback()
{
    const auto plan = node.getLocalizationPlannedTrajectoryFile();
    planLabel.setText (plan.isEmpty() ? juce::String ("no plan loaded") : juce::File (plan).getFileName(), juce::dontSendNotification);
    const auto bundle = node.getLocalizationBundleDirectory();
    bundleLabel.setText (bundle.isEmpty() ? juce::String ("default bundle") : juce::File (bundle).getFileName(), juce::dontSendNotification);
    status.setText (node.getLocalizationStatus(), juce::dontSendNotification);
}
}
