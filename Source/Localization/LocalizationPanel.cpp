#include "LocalizationPanel.hpp"

#include "../ProbeViewerNode.h"
#include "Realtime/PlanFile.hpp"

namespace ProbeViewer::Localization
{
namespace
{
// Lower-case letters and digits only: "Probe A" and "ProbeA-AP" share "probea".
juce::String squash (const juce::String& s)
{
    return s.toLowerCase().retainCharacters ("abcdefghijklmnopqrstuvwxyz0123456789");
}
} // namespace

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
    planButton.setTooltip ("Planned trajectory: an ephysatlas-plan/1 JSON file (ephysatlas.realtime.plan), a Pinpoint v5 "
                           "experiment file, a Pinpoint text export or an xyz track");
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
    chooser = std::make_unique<juce::FileChooser> ("Planned trajectory (ephysatlas plan, Pinpoint file or xyz track)", juce::File(), "*");
    chooser->launchAsync (juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectFiles,
                          [this] (const juce::FileChooser& fc)
                          {
                              const auto f = fc.getResult();
                              if (! f.existsAsFile())
                                  return;
                              const auto path = f.getFullPathName();
                              const auto names = Realtime::pinpointProbeNames (path.toStdString());
                              if (names.size() < 2)
                              {
                                  node.setLocalizationPlannedTrajectoryFile (path, names.empty() ? juce::String() : juce::String (names.front()));
                                  return;
                              }
                              // A Pinpoint experiment of several probes: pick one, ticking the displayed stream's.
                              const auto stream = squash (node.getDisplayedStream());
                              juce::PopupMenu menu;
                              menu.addSectionHeader ("Probe of " + f.getFileName());
                              for (size_t i = 0; i < names.size(); ++i)
                              {
                                  const auto name = squash (juce::String (names[i]));
                                  menu.addItem (static_cast<int> (i) + 1, juce::String (names[i]), true, name.isNotEmpty() && stream.contains (name));
                              }
                              juce::Component::SafePointer<LocalizationPanel> safe (this);
                              menu.showMenuAsync (juce::PopupMenu::Options().withTargetComponent (&planButton),
                                                  [safe, path, names] (int result)
                                                  {
                                                      if (safe != nullptr && result > 0)
                                                          safe->node.setLocalizationPlannedTrajectoryFile (
                                                              path, juce::String (names[static_cast<size_t> (result - 1)]));
                                                  });
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
    const auto probe = node.getLocalizationPlannedTrajectoryProbe();
    planLabel.setText (plan.isEmpty() ? juce::String ("no plan loaded")
                                      : juce::File (plan).getFileName() + (probe.isEmpty() ? juce::String() : " (" + probe + ")"),
                       juce::dontSendNotification);
    const auto bundle = node.getLocalizationBundleDirectory();
    bundleLabel.setText (bundle.isEmpty() ? juce::String ("default bundle") : juce::File (bundle).getFileName(), juce::dontSendNotification);
    status.setText (node.getLocalizationStatus(), juce::dontSendNotification);
}
}
