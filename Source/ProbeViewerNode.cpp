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

#include "ProbeViewerNode.h"

#include "ProbeViewerCanvas.h"
#include "ProbeViewerEditor.h"

#include <cmath>
#include <limits>
#include <vector>

using namespace ProbeViewer;

namespace
{
std::vector<int> getElectrodeIndices (const DataStream* stream)
{
    std::vector<int> result;
    if (stream == nullptr)
        return result;

    const auto channels = stream->getContinuousChannels();
    result.reserve (static_cast<size_t> (channels.size()));

    for (int i = 0; i < channels.size(); ++i)
    {
        uint16 electrodeIndex = static_cast<uint16> (i);
        const auto* channel = channels[i];
        const int metadataIndex = channel->findMetadata (
            MetadataDescriptor::MetadataType::UINT16,
            1,
            "neuropixels.electrode_index");
        if (metadataIndex >= 0)
        {
            const MetadataValue* metadataValue = channel->getMetadataValue (metadataIndex);
            if (metadataValue != nullptr)
                metadataValue->getValue (electrodeIndex);
        }
        result.push_back (static_cast<int> (electrodeIndex));
    }
    return result;
}

String deviceIdentityKey (const DataStream* stream)
{
    if (stream == nullptr || ! stream->hasDevice() || stream->device == nullptr)
        return {};

    const auto* device = stream->device;
    if (device->serial_number.isNotEmpty() && device->serial_number != "0000000")
        return "serial:" + device->serial_number;
    if (device->getIdentifier().isNotEmpty())
        return "identifier:" + device->getIdentifier();
    if (device->getName().isNotEmpty())
        return "device:" + device->getName() + "|" + device->manufacturer;
    return {};
}

int probeDisplayIndexFromName (String name)
{
    name = name.trim();
    if (name.length() == 6 && name.startsWithIgnoreCase ("Probe"))
    {
        const auto c = name.getLastCharacter();
        if (c >= 'A' && c <= 'Z')
            return static_cast<int> (c - 'A');
        if (c >= 'a' && c <= 'z')
            return static_cast<int> (c - 'a');
    }
    return -1;
}

Localization::ProbeIdentity makeProbeIdentity (
    const DataStream* selected,
    const Array<const DataStream*>& allStreams)
{
    Localization::ProbeIdentity probe;
    if (selected == nullptr)
        return probe;

    if (selected->hasDevice() && selected->device != nullptr)
    {
        const auto* device = selected->device;
        probe.deviceName = device->getName();
        probe.serialNumber = device->serial_number;
        probe.manufacturer = device->manufacturer;
        probe.deviceIdentifier = device->getIdentifier();
        probe.id = deviceIdentityKey (selected);
        probe.displayName = probe.deviceName.isNotEmpty()
            ? probe.deviceName
            : selected->getName();
        probe.displayIndex = probeDisplayIndexFromName (probe.displayName);

        if (probe.displayIndex < 0)
        {
            StringArray uniqueDevices;
            for (auto* stream : allStreams)
            {
                const auto key = deviceIdentityKey (stream);
                if (key.isNotEmpty())
                    uniqueDevices.addIfNotAlreadyThere (key);
            }
            probe.displayIndex = uniqueDevices.indexOf (probe.id);
        }
        return probe;
    }

    // Metadata-poor files/readers do not expose a physical DeviceInfo. Keep a
    // deterministic per-selected-stream identity now; stream pairing below is
    // still guarded against ambiguity. This field makes result/UI contracts
    // future-multi-probe ready without pretending that hardware identity exists.
    probe.id = "stream:" + selected->getKey();
    probe.displayName = selected->getName().isNotEmpty()
        ? selected->getName()
        : "Probe";
    probe.displayIndex = 0;
    return probe;
}

bool sameElectrodeLayout (const DataStream* a, const DataStream* b)
{
    if (a == nullptr || b == nullptr || a->getChannelCount() != b->getChannelCount())
        return false;
    return getElectrodeIndices (a) == getElectrodeIndices (b);
}

Localization::StreamSpec makeStreamSpec (const DataStream* stream)
{
    if (stream == nullptr)
        return {};

    Localization::StreamSpec spec;
    spec.key = stream->getKey();
    spec.sampleRate = stream->getSampleRate();
    spec.electrodeIndices = getElectrodeIndices (stream);

    // Only the leading ephys channels: a sync word saved as a continuous channel (SpikeGLX
    // SY, "SYNC0" in File Reader recordings, typed as an electrode there) is not localized.
    const auto channels = stream->getContinuousChannels();
    int ephys = 0;
    while (ephys < channels.size()
           && channels[ephys]->getChannelType() == ContinuousChannel::Type::ELECTRODE
           && ! channels[ephys]->getName().startsWithIgnoreCase ("SYNC"))
        ++ephys;
    spec.channels = ephys;
    if (static_cast<int> (spec.electrodeIndices.size()) > ephys)
        spec.electrodeIndices.resize (static_cast<size_t> (ephys));
    return spec;
}

bool discoverBandPair (
    const DataStream* selected,
    const Array<const DataStream*>& streams,
    Localization::StreamSpec& ap,
    Localization::StreamSpec& lf,
    String& failureReason)
{
    ap = {};
    lf = {};
    failureReason.clear();
    if (selected == nullptr)
    {
        failureReason = "selected stream is unavailable";
        return false;
    }

    Array<const DataStream*> candidates;
    const auto selectedDeviceKey = deviceIdentityKey (selected);

    if (selectedDeviceKey.isNotEmpty())
    {
        // Primary path: DeviceInfo identifies the physical probe. No stream-name
        // parsing is involved. AP/LF are distinguished only inside that device.
        for (auto* stream : streams)
        {
            if (deviceIdentityKey (stream) == selectedDeviceKey
                && sameElectrodeLayout (selected, stream)
                && stream->getSampleRate() > 0.0f)
            {
                candidates.add (stream);
            }
        }
    }
    else
    {
        // File Reader / synthetic fallback: use the exact Neuropixels electrode
        // layout. If more than two such streams exist, pairing would be unsafe
        // for a future multi-probe recording, so refuse instead of guessing.
        for (auto* stream : streams)
        {
            if (sameElectrodeLayout (selected, stream)
                && stream->getSampleRate() > 0.0f)
            {
                candidates.add (stream);
            }
        }

        if (candidates.size() > 2)
        {
            failureReason =
                "streams have no DeviceInfo and electrode metadata matches more than one AP/LF pair";
            return false;
        }
    }

    if (candidates.size() < 2)
    {
        failureReason =
            selectedDeviceKey.isNotEmpty()
                ? "physical probe does not expose two matching-rate streams"
                : "could not find an unambiguous second stream with the same electrode layout";
        return false;
    }

    const DataStream* highest = nullptr;
    const DataStream* lowest = nullptr;
    float highestRate = -1.0f;
    float lowestRate = std::numeric_limits<float>::max();

    for (auto* stream : candidates)
    {
        const float rate = stream->getSampleRate();
        if (rate > highestRate)
        {
            highestRate = rate;
            highest = stream;
        }
        if (rate < lowestRate)
        {
            lowestRate = rate;
            lowest = stream;
        }
    }

    if (highest == nullptr || lowest == nullptr || highest == lowest
        || std::abs (highestRate - lowestRate) < 1.0f)
    {
        failureReason = "matching streams do not contain distinct AP/LF sample rates";
        return false;
    }

    ap = makeStreamSpec (highest);
    lf = makeStreamSpec (lowest);
    return true;
}
}


ProbeViewerNode::ProbeViewerNode()
    : GenericProcessor ("Probe Viewer")
{
    streamToDraw = String();
    localizationWorker = std::make_unique<Localization::LocalizationWorker>();
    localizationAccumulator = std::make_unique<Localization::LocalizationWindowAccumulator>();
}

ProbeViewerNode::~ProbeViewerNode()
{
    localizationAccumulator.reset();
    localizationWorker.reset();
}

void ProbeViewerNode::registerParameters()
{
    addIntParameter (Parameter::PROCESSOR_SCOPE, "trigger_line", "Trigger Line", "The TTL trigger line", -1, -1, 15);
    addSelectedStreamParameter (Parameter::PROCESSOR_SCOPE, "display_stream", "Display Stream", "The stream to display", {}, 0, true, false);
}

AudioProcessorEditor* ProbeViewerNode::createEditor()
{
    editor = std::make_unique<ProbeViewerEditor> (this);
    return editor.get();
}

void ProbeViewerNode::process (AudioBuffer<float>& buffer)
{
    for (auto stream : getDataStreams())
    {
        String streamKey = stream->getKey();
        uint16 streamId = stream->getStreamId();
        CircularBuffer* streamBuffer = dataBufferMap[streamKey];
        for (int localId = 0; localId < stream->getChannelCount(); localId++)
        {
            int globalId = stream->getContinuousChannels()[localId]->getGlobalIndex();
            streamBuffer->addData (buffer, localId, globalId, getNumSamplesInBlock (streamId), getFirstSampleNumberForBlock (streamId));
        }

        // Localization: copy the block into the snippet being captured, if any.
        if (localizationAccumulator->isArmed())
        {
            std::vector<int> globalChannelIndices;
            globalChannelIndices.reserve (static_cast<size_t> (stream->getChannelCount()));
            for (int localId = 0; localId < stream->getChannelCount(); localId++)
                globalChannelIndices.push_back (stream->getContinuousChannels()[localId]->getGlobalIndex());
            localizationAccumulator->addBlock (streamKey, buffer, globalChannelIndices, getNumSamplesInBlock (streamId));
        }
    }

    // A complete snippet goes to the localization worker (which runs on its own thread).
    if (localizationAccumulator->isWindowReady() && localizationWorker->canAcceptJob())
        localizationWorker->submitWindow (localizationAccumulator->makeInput(), localizationConfig);

    checkForEvents();
}

void ProbeViewerNode::handleTTLEvent (TTLEventPtr event)
{
    const int eventState = event->getState() ? 1 : 0;
    const int eventLine = event->getLine();
    const int64 sampleNumber = event->getSampleNumber();
    const uint16 streamId = event->getChannelInfo()->getStreamId();

    if (eventState && eventLine == (int) getParameter ("trigger_line")->getValue())
    {
        dataBufferMap[getDataStream (streamId)->getKey()]->setTrigger (sampleNumber);
    }
}

void ProbeViewerNode::updateSettings()
{
    for (auto buffer : dataBuffers)
    {
        buffer->prepareToUpdate();
    }

    for (auto stream : getDataStreams())
    {
        String streamKey = stream->getKey();

        if (dataBufferMap.count (streamKey) == 0)
        {
            dataBuffers.add (new CircularBuffer (streamKey, stream->getSampleRate(), bufferLengthSeconds));
            dataBufferMap[streamKey] = dataBuffers.getLast();
        }
        else
        {
            dataBufferMap[streamKey]->sampleRate = stream->getSampleRate();
        }

        dataBufferMap[streamKey]->updateChannelInfo (stream->getContinuousChannels());
    }

    Array<CircularBuffer*> toDelete;

    for (auto dataBuffer : dataBuffers)
    {
        if (dataBuffer->isNeeded)
        {
            dataBuffer->update();
        }
        else
        {
            dataBufferMap.erase (dataBuffer->key);
            toDelete.add (dataBuffer);
        }
    }

    for (auto dataBuffer : toDelete)
    {
        dataBuffers.removeObject (dataBuffer, true);
    }
}

void ProbeViewerNode::parameterValueChanged (Parameter* param)
{
    if (param->getName() == "display_stream")
    {
        String streamKey = param->getValueAsString();
        if (auto stream = getDataStream (streamKey))
            setDisplayedStream (streamKey);
        else
            setDisplayedStream (String());
    }
}

bool ProbeViewerNode::startAcquisition()
{
    ((ProbeViewerEditor*) getEditor())->enable();
    return true;
}

bool ProbeViewerNode::stopAcquisition()
{
    ((ProbeViewerEditor*) getEditor())->disable();
    return true;
}

void ProbeViewerNode::setDisplayedStream (const String& streamKey)
{
    if (streamToDraw != streamKey)
        localizationAccumulator->reset();
    streamToDraw = streamKey;
    ((ProbeViewerEditor*) getEditor())->displayStreamChanged();
}

String ProbeViewerNode::getDisplayedStream() const
{
    return streamToDraw;
}

float ProbeViewerNode::getStreamSampleRate()
{
    if (streamToDraw.isNotEmpty())
        return getDataStream (streamToDraw)->getSampleRate();
    else
        return 0.0f;
}

int ProbeViewerNode::getNumStreamChannels()
{
    if (streamToDraw.isNotEmpty())
        return getDataStream (streamToDraw)->getChannelCount();
    else
        return 0;
}

CircularBuffer* ProbeViewerNode::getCircularBufferPtr()
{
    if (streamToDraw.isNotEmpty() && dataBufferMap.count (streamToDraw) > 0)
        return dataBufferMap[streamToDraw];
    else
        return nullptr;
}

String ProbeViewerNode::handleConfigMessage (const String& msg)
{
    // message format
    // "<probe_name>;<start_index_1>-<end_index_1>,<region_ID_1>,<hex_color_1>;<start_index_2>-<end_index_2>,...

    // Example:
    // ProbeA;0-69,PT,FF909F;70-97,PVT,FF909F;98-161,-,000000;162-173,-,000000,174-185,SF,90CBED;...

    //LOGD(msg);
    LOGD ("Probe Viewer ", getNodeId(), " received message of length ", msg.length());

    int firstSemicolon = msg.indexOf (";");
    String probeName;
    String electrodeInfo;

    if (firstSemicolon > -1)
    {
        probeName = msg.substring (0, firstSemicolon);
        electrodeInfo = msg.substring (firstSemicolon + 1);
        LOGD ("Probe name: ", probeName);
    }
    else
    {
        LOGD ("No probe name detected.");
        return "No probe name detected.";
    }

    Array<String> streamKeys;

    for (auto stream : dataStreams)
    {
        String streamName = stream->getName();

        if (streamName.endsWith ("-AP"))
        {
            streamName = streamName.substring (0, streamName.length() - 3);
        }
        else if (streamName.endsWith ("-LFP"))
        {
            streamName = streamName.substring (0, streamName.length() - 4);
        }

        if (probeName.equalsIgnoreCase (streamName))
            streamKeys.add (stream->getKey());
    }

    LOGD ("Number of matching streams: ", streamKeys.size());

    if (streamKeys.size() == 0)
    {
        return "No matching stream detected.";
    }

    Array<int> electrodeInds;
    Array<String> regionNames;
    Array<Colour> regionColours;

    StringArray tokens = StringArray::fromTokens (electrodeInfo, ";", "");

    for (auto token : tokens)
    {
        //LOGC(token);

        StringArray rangeInfo = StringArray::fromTokens (token, ",", "");

        int firstElectrode = -1;
        int lastElectrode = -1;
        String regionName = "";
        Colour regionColour = Colours::black;

        int hyphen = rangeInfo[0].indexOf ("-");

        if (hyphen > -1)
        {
            firstElectrode = rangeInfo[0].substring (0, hyphen).getIntValue();
            lastElectrode = rangeInfo[0].substring (hyphen + 1).getIntValue();
            //LOGD("First electrode: ", firstElectrode);
            //LOGD("Last electrode: ", lastElectrode);
        }

        if (rangeInfo.size() > 1)
        {
            regionName = rangeInfo[1];

            if (regionName.startsWith ("SSp"))
                regionName = "SSp";
            //LOGD("Region name: ", regionName);
        }

        if (rangeInfo.size() > 2)
        {
            regionColour = Colour::fromString ("#FF" + rangeInfo[2].toUpperCase());
            //LOGD("Original color: ", rangeInfo[2]);
            //LOGD("Region colour: ", regionColour.toString());
        }

        if (firstElectrode > -1 && lastElectrode > -1)
        {
            //LOGD("Adding new range.");

            for (int i = firstElectrode; i < lastElectrode + 1; i++)
            {
                electrodeInds.add (i);
                regionNames.add (regionName);
                regionColours.add (regionColour);
            }
        }
    }

    if (electrodeInds.size() > 0)
    {
        ProbeViewerEditor* ed = (ProbeViewerEditor*) getEditor();

        for (auto streamKey : streamKeys)
            ed->setRegions (streamKey, electrodeInds, regionNames, regionColours);
    }

    return "Success";
}

void ProbeViewerNode::setLocalizationBundleDirectory (const String& directory)
{
    localizationConfig.realtime.bundleDirectory = directory.trim();
}

String ProbeViewerNode::getLocalizationBundleDirectory() const
{
    return localizationConfig.realtime.bundleDirectory;
}

void ProbeViewerNode::setLocalizationPlannedTrajectoryFile (const String& filename)
{
    localizationConfig.initialTrajectoryFile = filename.trim();
}

String ProbeViewerNode::getLocalizationPlannedTrajectoryFile() const
{
    return localizationConfig.initialTrajectoryFile;
}

void ProbeViewerNode::requestLocalization()
{
    if (streamToDraw.isEmpty() || localizationWorker == nullptr
        || ! localizationWorker->canAcceptJob())
    {
        LOGD ("Probe Viewer localization request ignored: no stream or worker busy.");
        return;
    }

    auto* selected = getDataStream (streamToDraw);
    if (selected == nullptr)
    {
        LOGD ("Probe Viewer localization request ignored: selected stream is unavailable.");
        return;
    }

    const auto streams = getDataStreams();
    Localization::StreamSpec ap;
    Localization::StreamSpec lf;
    String failureReason;
    if (! discoverBandPair (selected, streams, ap, lf, failureReason))
    {
        LOGD ("Probe Viewer localization could not identify AP/LF streams: ", failureReason);
        for (auto* stream : streams)
        {
            const auto deviceName = stream->hasDevice() && stream->device != nullptr
                ? stream->device->getName()
                : String ("<no device>");
            LOGD ("  stream=", stream->getName(),
                  " key=", stream->getKey(),
                  " device=", deviceName,
                  " rate=", stream->getSampleRate(),
                  " channels=", stream->getChannelCount());
        }
        return;
    }

    const auto probe = makeProbeIdentity (selected, streams);

    const int pairedChannels = juce::jmin (ap.channels, lf.channels);
    if (pairedChannels <= 0)
    {
        LOGD ("Probe Viewer localization rejected streams with no shared channels.");
        return;
    }

    ap.channels = pairedChannels;
    lf.channels = pairedChannels;
    if (static_cast<int> (ap.electrodeIndices.size()) > pairedChannels)
        ap.electrodeIndices.resize (static_cast<size_t> (pairedChannels));
    if (static_cast<int> (lf.electrodeIndices.size()) > pairedChannels)
        lf.electrodeIndices.resize (static_cast<size_t> (pairedChannels));

    // One snippet per window, taken back to back for as long as acquisition runs.
    const double apSeconds = localizationConfig.realtime.snippetSeconds;
    const double lfSeconds = localizationConfig.realtime.snippetSeconds;
    const int windowCount = localizationConfig.realtime.maxSnippets;

    const double estimatedBytes = sizeof (float)
        * (apSeconds * ap.sampleRate * ap.channels
           + lfSeconds * lf.sampleRate * lf.channels);
    constexpr double maxInMemoryBytes = 2.0 * 1024.0 * 1024.0 * 1024.0;
    if (estimatedBytes > maxInMemoryBytes)
    {
        LOGD ("Probe Viewer localization window would require ",
              estimatedBytes / (1024.0 * 1024.0),
              " MiB. Reduce the window or use chunked accumulation.");
        return;
    }

    localizationAccumulator->arm (
        probe,
        ap,
        lf,
        apSeconds,
        lfSeconds,
        windowCount);

    LOGD ("Probe Viewer localization armed. Probe=", probe.displayName,
          " id=", probe.id,
          "; AP key=", ap.key, " rate=", ap.sampleRate,
          "; LF key=", lf.key, " rate=", lf.sampleRate,
          "; channels=", pairedChannels,
          "; plan=", localizationConfig.initialTrajectoryFile.isEmpty()
              ? String ("none yet")
              : localizationConfig.initialTrajectoryFile);
}

void ProbeViewerNode::serviceLocalization()
{
    // Allocate the next raw window only after the worker has finished processing the previous
    // one, which keeps peak RAM bounded to one AP/LF pair. Called from the message thread (the
    // editor's timer), not from process(): a 5 s AP window is ~230 MB, too much to allocate on
    // the audio thread every snippet.
    if (localizationAccumulator != nullptr
        && localizationWorker != nullptr
        && ! localizationAccumulator->isArmed()
        && localizationAccumulator->hasPendingBatchWindow()
        && localizationWorker->canAcceptJob())
    {
        localizationAccumulator->continueBatchIfNeeded();
    }
}

void ProbeViewerNode::resetLocalization()
{
    if (localizationAccumulator != nullptr)
        localizationAccumulator->reset();
    if (localizationWorker != nullptr)
        localizationWorker->reset();
}

Localization::LocalizationResultPtr ProbeViewerNode::getLatestLocalizationResult() const
{
    if (localizationWorker != nullptr)
        return localizationWorker->getLatestResult();
    return {};
}

String ProbeViewerNode::getLocalizationStatus() const
{
    if (localizationWorker == nullptr)
        return "Unavailable";
    String text = localizationWorker->getStatus();
    if (localizationAccumulator != nullptr && localizationAccumulator->isArmed())
        text += "\nCollecting snippet " + String (localizationWorker->getSnippetCount() + 1);
    else if (localizationAccumulator != nullptr && localizationAccumulator->isWindowReady())
        text += "\nSnippet " + String (localizationWorker->getSnippetCount() + 1) + " waiting";
    return text;
}

const float ProbeViewerNode::bufferLengthSeconds = 10.0f;
