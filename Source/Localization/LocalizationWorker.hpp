#pragma once
#include "LocalizationEngine.hpp"
#include <condition_variable>
#include <thread>

namespace ProbeViewer::Localization
{
class LocalizationWorker
{
public:
    LocalizationWorker();
    ~LocalizationWorker();
    bool canAcceptJob() const noexcept;
    bool submitWindow (LocalizationInput, LocalizationConfig);
    LocalizationResultPtr getLatestResult() const;
    juce::String getStatus() const;
    void reset();
private:
    void loop();
    mutable std::mutex mutex;
    std::condition_variable cv;
    std::thread thread;
    bool quitting = false, hasJob = false, busy = false;
    LocalizationInput pendingInput;
    LocalizationConfig pendingConfig;
    LocalizationResultPtr latest;
    juce::String status { "Idle" };
    LocalizationEngine engine;
};
}
