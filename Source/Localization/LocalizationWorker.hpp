#pragma once
#include "LocalizationConfig.hpp"
#include "LocalizationResult.hpp"
#include "Realtime/RealtimeEngine.hpp"

#include <atomic>
#include <condition_variable>
#include <memory>
#include <thread>

namespace ProbeViewer::Localization
{
// Runs the real-time ephys-only localizer on a background thread, one captured window at a time:
// each window is one AP + LF snippet, whose features are added before the probe is localized
// from every snippet so far.
class LocalizationWorker
{
public:
    LocalizationWorker();
    ~LocalizationWorker();
    bool canAcceptJob() const noexcept;
    bool submitWindow (LocalizationInput, LocalizationConfig);
    LocalizationResultPtr getLatestResult() const;
    juce::String getStatus() const;
    float getProgress() const;
    int getSnippetCount() const;
    // Forgets the result and the real-time snippets (new recording, probe or plan).
    void reset();

private:
    void loop();
    LocalizationResult runRealtime (LocalizationInput&, const LocalizationConfig&);
    bool prepareRealtime (const LocalizationConfig&, int nChannels, juce::String& error);

    mutable std::mutex mutex;
    std::condition_variable cv;
    std::thread thread;
    bool quitting = false, hasJob = false, busy = false, resetPending = false;
    LocalizationInput pendingInput;
    LocalizationConfig pendingConfig;
    LocalizationResultPtr latest;
    juce::String status { "Idle" };
    std::atomic<float> progress { 0.0f };
    std::atomic<int> snippetCount { 0 };

    // Real-time localizer state (worker thread only).
    std::unique_ptr<Realtime::LocalizerBundle> bundle;
    std::unique_ptr<Realtime::RealtimeEngine> realtime;
    juce::String bundleDirectory, planFile, planSource;
    std::vector<Realtime::V3> plan;
};
}
