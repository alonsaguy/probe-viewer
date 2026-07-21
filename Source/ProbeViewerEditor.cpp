/*
 ------------------------------------------------------------------
 
 This file is part of the Open Ephys GUI
 Copyright (C) 2013 Open Ephys
 
 ------------------------------------------------------------------
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 */

#include "ProbeViewerEditor.h"

#include "ProbeViewerCanvas.h"
#include "ProbeViewerNode.h"

using namespace ProbeViewer;

ProbeViewerEditor::ProbeViewerEditor (GenericProcessor* parentNode)
    : VisualizerEditor (parentNode, "ProbeViewer"),
      hasNoInputs (true)
{
    probeViewerProcessor = (ProbeViewerNode*) parentNode;

    // Keep the existing stream controls on the left and place localization
    // controls in a new section to their right. This avoids relying on extra
    // vertical editor space, which is limited in the Open Ephys signal chain.
    desiredWidth = 430;

    addSelectedStreamParameterEditor (Parameter::PROCESSOR_SCOPE, "display_stream", 15, 30);
    auto* streamEditor = getParameterEditor ("display_stream");
    streamEditor->setLayout (ParameterEditor::nameOnTop);
    streamEditor->setSize (140, 40);

    streamSampleRateLabel = std::make_unique<Label> ("Stream Sample Rate Label", "Sample Rate:");
    streamSampleRateLabel->setFont (FontOptions ("Inter", "Medium", 14.0f));
    streamSampleRateLabel->setBounds (15, 90, 140, 24);
    addAndMakeVisible (streamSampleRateLabel.get());

    localizationPanel = std::make_unique<Localization::LocalizationPanel> (*probeViewerProcessor);
    localizationPanel->setBounds (170, 8, 250, 138);
    addAndMakeVisible (localizationPanel.get());

    // Poll only on the GUI thread. The worker never directly touches canvas
    // objects or JUCE components.
    startTimerHz (5);
}

ProbeViewerEditor::~ProbeViewerEditor()
{
    stopTimer();
}

Visualizer* ProbeViewerEditor::createNewCanvas()
{
    return new ProbeViewerCanvas (probeViewerProcessor);
}

void ProbeViewerEditor::updateSettings()
{
    if (canvas != nullptr)
    {
        static_cast<ProbeViewerCanvas*> (canvas.get())->updateChannelBrowsers();
    }
}

void ProbeViewerEditor::displayStreamChanged()
{
    if (canvas != nullptr)
    {
        static_cast<ProbeViewerCanvas*> (canvas.get())->updateSettings();
    }

    float rate = probeViewerProcessor->getStreamSampleRate();

    String sampleRateLabelText = "Sample Rate: ";

    if (rate > 0)
        sampleRateLabelText += String (rate);
    else
        sampleRateLabelText += "<NA>";

    streamSampleRateLabel->setText (sampleRateLabelText, dontSendNotification);
}

/** Sets custom depths and regions */
void ProbeViewerEditor::setRegions (String streamKey,
                                    Array<int>& electrodeInds,
                                    Array<String>& regionNames,
                                    Array<Colour>& regionColours)
{
    if (canvas != nullptr)
    {
        static_cast<ProbeViewerCanvas*> (canvas.get())->setRegions (streamKey,
                                                                    electrodeInds,
                                                                    regionNames,
                                                                    regionColours);
    }
}

void ProbeViewerEditor::timerCallback()
{
    auto result = probeViewerProcessor->getLatestLocalizationResult();

    if (result == nullptr
        || ! result->valid
        || result == displayedLocalizationResult)
    {
        return;
    }

    displayedLocalizationResult = result;

    Array<int> electrodeInds;
    Array<String> regionNames;
    Array<Colour> regionColours;

    electrodeInds.ensureStorageAllocated (static_cast<int> (result->channels.size()));
    regionNames.ensureStorageAllocated (static_cast<int> (result->channels.size()));
    regionColours.ensureStorageAllocated (static_cast<int> (result->channels.size()));

    for (const auto& channel : result->channels)
    {
        electrodeInds.add (channel.channelIndex);
        regionNames.add (channel.acronym);
        regionColours.add (channel.colour);
    }

    setRegions (result->streamKey,
                electrodeInds,
                regionNames,
                regionColours);
}
