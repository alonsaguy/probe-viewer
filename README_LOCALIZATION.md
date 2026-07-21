# Probe Viewer localization scaffold

This scaffold extends the existing Probe Viewer processor. It does not create a second Open Ephys processor.

## Important integration note

The supplied `ProbeViewerNode.cpp` contains a placeholder `handleConfigMessage()` because the repository's current function is long and unrelated to localization. Copy the original repository implementation of `handleConfigMessage()` back into that function after replacing the file. Everything else can be copied directly.

## Behavior

1. Enable localization in the editor.
2. Choose the displayed stream.
3. Press **Localize**.
4. `process()` only copies the selected stream into a preallocated 100-second accumulator.
5. Once full, ownership of the window moves to `LocalizationWorker`.
6. The worker extracts cached Probe Viewer features when supplied; otherwise it falls back to raw per-channel RMS.
7. `LocalizationEngine` returns a random demo trajectory in approximate CCF bounds.
8. The editor timer publishes predicted channel regions through the existing `ProbeViewerCanvas::setRegions()` mechanism.

## Probe Viewer feature reuse

The current Probe Viewer computes display features inside visualization classes but does not expose a stable processor-side feature snapshot. Therefore `LocalizationInput::probeViewerFeatures` is the integration point. Until the display feature cache is exposed, selecting “Probe Viewer features” safely falls back to raw-window RMS. A later update should add a thread-safe snapshot method to the heatmap feature owner and fill:

- `LocalizationInput::probeViewerFeatures`
- `LocalizationInput::featuresPerChannel`

Do not read GUI-owned arrays directly from `process()`.

## Build

No CMake change is needed because the repository uses recursive globbing for `.cpp`, `.h`, and `.hpp` files under `Source/`.

Regenerate CMake so Visual Studio sees the new files:

```powershell
cd Build
cmake -G "Visual Studio 17 2022" -A x64 ..
```

## Standalone engine test

`LocalizationEngine` has no dependency on `GenericProcessor`, `Visualizer`, or the Open Ephys callback. A small console target can instantiate `LocalizationInput`, call `LocalizationEngine::run`, and print the returned CCF points. JUCE remains used for strings/colors; these can later be replaced by standard types if a completely JUCE-free core is preferred.
