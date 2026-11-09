#include "Threads.hpp"

#include <algorithm>
#include <cstdlib>
#include <thread>

#if defined(_WIN32)
 #ifndef WIN32_LEAN_AND_MEAN
  #define WIN32_LEAN_AND_MEAN
 #endif
 #ifndef NOMINMAX
  #define NOMINMAX
 #endif
 #include <windows.h>
#elif defined(__APPLE__)
 #include <pthread.h>
#else
 #include <sys/resource.h>
 #include <sys/syscall.h>
 #include <unistd.h>
#endif

namespace ProbeViewer::Localization::Realtime
{
int localizerThreads()
{
    static const int threads = []
    {
        if (const char* env = std::getenv ("PROBE_VIEWER_LOCALIZER_THREADS"))
        {
            const int n = std::atoi (env);
            if (n >= 1)
                return n;
        }
        const int logical = static_cast<int> (std::thread::hardware_concurrency());
        return std::clamp (logical / 2, 1, 8);
    }();
    return threads;
}

void lowerThreadPriority()
{
#if defined(_WIN32)
    SetThreadPriority (GetCurrentThread(), THREAD_PRIORITY_BELOW_NORMAL);
#elif defined(__APPLE__)
    pthread_set_qos_class_self_np (QOS_CLASS_UTILITY, 0);
#else
    // Per-thread nice value on Linux; threads created by this one inherit it.
    setpriority (PRIO_PROCESS, static_cast<id_t> (syscall (SYS_gettid)), 10);
#endif
}
} // namespace ProbeViewer::Localization::Realtime
