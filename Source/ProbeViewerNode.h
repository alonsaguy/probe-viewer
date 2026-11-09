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

#ifndef __PROBEVIEWERNODE_H__
#define __PROBEVIEWERNODE_H__

#include "ProcessorHeaders.h"
#include "Utilities/CircularBuffer.hpp"

#include "Localization/LocalizationConfig.hpp"
#include "Localization/LocalizationResult.hpp"
#include "Localization/LocalizationWindowAccumulator.hpp"
#include "Localization/LocalizationWorker.hpp"

#include <memory>

namespace ProbeViewer
{

class ProbeViewerNode : public GenericProcessor
{
public:
    /** Constructor */
    ProbeViewerNode();

    /** Destructor */
    virtual ~ProbeViewerNode() override;

    void registerParameters() override;

    /** Creates the editor UI */
    AudioProcessorEditor* createEditor() override;

    /** Pushes samples to the data buffer*/
    void process (AudioBuffer<float>& buffer) override;

    /** Handles incoming events used as triggers */
    void handleTTLEvent (TTLEventPtr event);

    /** Updates settings */
    void updateSettings() override;

    void parameterValueChanged (Parameter*) override;

    /** Updates the displayed stream, then calls updateSettings() */
    void setDisplayedStream (const String& streamKey);

    /** Gets the displayed stream key*/
    String getDisplayedStream() const;

    /** Returns the sample rate of the currently selected stream*/
    float getStreamSampleRate();

    /** Returns the number of samples in the currently selected stream*/
    int getNumStreamChannels();

    /** Enables the editor */
    bool startAcquisition() override;

    /** Disables the editor*/
    bool stopAcquisition() override;

    /** Returns the buffer for the currently selected stream*/
    CircularBuffer* getCircularBufferPtr();

    /** Responds to config messages with region info */
    String handleConfigMessage (const String& msg) override;

    /** Real-time probe localization (Localization/): the localizer bundle directory (empty: the
        default one) and the planned trajectory file (ephysatlas plan JSON, Pinpoint file or xyz
        trajectory), with the probe to read from a Pinpoint v5 experiment of several probes */
    void setLocalizationBundleDirectory (const String& directory);
    String getLocalizationBundleDirectory() const;
    void setLocalizationPlannedTrajectoryFile (const String& filename, const String& probe = {});
    String getLocalizationPlannedTrajectoryFile() const;
    String getLocalizationPlannedTrajectoryProbe() const;

    /** Starts localizing the displayed probe: AP + LF snippets are captured back to back from now
        on, and each one updates the estimate */
    void requestLocalization();

    /** Starts capturing the next snippet once the worker is free (called on the message thread) */
    void serviceLocalization();

    /** Forgets the snippets and the estimate */
    void resetLocalization();

    /** The latest localization (null before the first one) */
    Localization::LocalizationResultPtr getLatestLocalizationResult() const;

    /** What the localization is doing, for the editor */
    String getLocalizationStatus() const;

private:
    static const float bufferLengthSeconds;

    OwnedArray<CircularBuffer> dataBuffers;

    std::map<String, CircularBuffer*> dataBufferMap;

    String streamToDraw;
    int lastChannelInStream;

    Array<bool> channelsToDraw;

    Localization::LocalizationConfig localizationConfig;
    std::unique_ptr<Localization::LocalizationWorker> localizationWorker;
    std::unique_ptr<Localization::LocalizationWindowAccumulator> localizationAccumulator;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (ProbeViewerNode);
};

} // namespace ProbeViewer

#endif /* __PROBEVIEWERNODE_H__ */
