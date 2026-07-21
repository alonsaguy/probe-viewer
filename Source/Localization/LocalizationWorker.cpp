#include "LocalizationWorker.hpp"

namespace ProbeViewer::Localization
{
LocalizationWorker::LocalizationWorker() : thread ([this] { loop(); }) {}
LocalizationWorker::~LocalizationWorker()
{
    { std::lock_guard<std::mutex> lock (mutex); quitting = true; }
    cv.notify_one();
    if (thread.joinable()) thread.join();
}

bool LocalizationWorker::canAcceptJob() const noexcept
{
    std::lock_guard<std::mutex> lock (mutex); return ! busy && ! hasJob;
}

bool LocalizationWorker::submitWindow (LocalizationInput input, LocalizationConfig config)
{
    std::lock_guard<std::mutex> lock (mutex);
    if (busy || hasJob) return false;
    pendingInput = std::move (input); pendingConfig = std::move (config);
    hasJob = true; status = "Queued"; cv.notify_one(); return true;
}

LocalizationResultPtr LocalizationWorker::getLatestResult() const
{ std::lock_guard<std::mutex> lock (mutex); return latest; }
juce::String LocalizationWorker::getStatus() const
{ std::lock_guard<std::mutex> lock (mutex); return status; }
void LocalizationWorker::reset()
{ std::lock_guard<std::mutex> lock (mutex); latest.reset(); hasJob = false; status = busy ? "Running" : "Idle"; }

void LocalizationWorker::loop()
{
    for (;;)
    {
        LocalizationInput input; LocalizationConfig config;
        {
            std::unique_lock<std::mutex> lock (mutex);
            cv.wait (lock, [this] { return quitting || hasJob; });
            if (quitting) return;
            input = std::move (pendingInput); config = pendingConfig; hasJob = false; busy = true; status = "Running";
        }
        auto result = std::make_shared<LocalizationResult> (engine.run (input, config));
        {
            std::lock_guard<std::mutex> lock (mutex);
            latest = std::move (result); busy = false; status = latest->status;
        }
    }
}
}
