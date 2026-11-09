#pragma once

// CPU budget of the localization, which shares the machine with acquisition and the UI.

namespace ProbeViewer::Localization::Realtime
{
// Threads of each parallel step: $PROBE_VIEWER_LOCALIZER_THREADS when set (1 or more), else half
// the logical processors, at most 8.
int localizerThreads();

// Puts the calling thread below normal OS priority, so that the localization only takes the CPU
// time that acquisition, the UI and other programs leave idle.
void lowerThreadPriority();
} // namespace ProbeViewer::Localization::Realtime
