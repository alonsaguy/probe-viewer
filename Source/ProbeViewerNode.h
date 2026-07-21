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

#include <atomic>
#include <memory>
#include <mutex>

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

    /** Pushes samples to the data buffer and, when armed, to the
        localization accumulator. No feature extraction or localization
        is performed on the acquisition thread. */
    void process (AudioBuffer<float>& buffer) override;

    /** Handles incoming events used as triggers */
    void handleTTLEvent (TTLEventPtr event);

    /** Updates settings */
    void updateSettings() override;

    void parameterValueChanged (Parameter*) override;

    /** Updates the displayed stream, then calls updateSettings() */
    void setDisplayedStream (const String& streamKey);

    /** Gets the displayed stream key */
    String getDisplayedStream() const;

    /** Returns the sample rate of the currently selected stream */
    float getStreamSampleRate();

    /** Returns the number of channels in the currently selected stream */
    int getNumStreamChannels();

    /** Enables the editor */
    bool startAcquisition() override;

    /** Disables the editor */
    bool stopAcquisition() override;

    /** Returns the buffer for the currently selected stream */
    CircularBuffer* getCircularBufferPtr();

    /** Responds to config messages with region info */
    String handleConfigMessage (const String& msg) override;

    // -----------------------------------------------------------------
    // Localization control API. These methods are called by the editor's
    // LocalizationPanel. Localization itself runs on LocalizationWorker.
    // -----------------------------------------------------------------
    void setLocalizationEnabled (bool enabled);
    void setLocalizationWindowSeconds (double seconds);
    void setLocalizationFeatureSource (Localization::FeatureSource source);
    void setLocalizationModelName (const String& modelName);

    /** Arms one localization window. Samples begin accumulating only after
        this method is called, normally by pressing the Localize button. */
    void requestLocalization();

    /** Clears the current window, queued result, and displayed-result state. */
    void resetLocalization();

    /** Thread-safe snapshot of the latest completed result. */
    Localization::LocalizationResultPtr getLatestLocalizationResult() const;

    /** Short status string for the localization panel. */
    String getLocalizationStatus() const;

    /** Fraction of the requested window collected, in [0, 1]. */
    float getLocalizationProgress() const;

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
