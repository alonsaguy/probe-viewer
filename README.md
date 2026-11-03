# Probe Viewer

![probe-viewer-editor](https://open-ephys.github.io/gui-docs/_images/probeviewer-01.png)

Visualizes the signal from a high-channel-count linear probe as a heatmap (**time** x **channels**). Can be configured to display RMS signal, spike rate, or power in different frequency bands.

## Installation

This plugin can be added via the Open Ephys GUI Plugin Installer. To access the Plugin Installer, press **ctrl-P** or **⌘P** from inside the GUI. Once the installer is loaded, browse to the "Probe Viewer" plugin and click "Install."

## Usage

### Plugin Editor

`Displayed Stream` controls the currently selected stream to be viewed. The selected stream can be changed while acquisition is active.

### Visualizer Window

![probe-viewer-visualizer](https://open-ephys.github.io/gui-docs/_images/probeviewer-02.png)

#### Channel Browser

Used to select a region of the probe to view.

#### Display Options

  * **Render Mode** - RMS amplitude, Frequency Band Power, or Spike Rate.
  * **Colour Scheme** - Inferno, Plasma, Magma, Viridis, Jet 

## Real-time probe localization (experimental)

The plugin can estimate, while you record, which brain region each channel of a Neuropixels 1.0 probe is in. It uses the LFP and AP bands only (no spike sorting or spike detection) and starts from the planned trajectory of the insertion. Every 5 s snippet of AP + LFP data updates the estimate; the regions appear next to the channels in the Channel Browser, and the **Localization** tab of the visualizer shows the estimated track (solid) and the plan (dashed) on tilted coronal and sagittal atlas slices.

The method is the ephys-only localizer of the [IBL ephys atlas](https://github.com/int-brain-lab/ibleatools) (`ephysatlas.alignment.ephys_only`), restricted to its 27 LFP and AP features: the plugin computes the same features as the IBL Python pipeline (bad-channel detection, destriping, Cadzow-denoised bipolar LFP, aperiodic fits), aggregates the snippets as the atlas does, and searches the atlas for the trajectory whose predicted features best explain them, under a prior on how far real insertions end up from their plan. The predictions come from the IBL channel model [`int-brain-lab/ea-encoder-channel`](https://huggingface.co/int-brain-lab/ea-encoder-channel) (release 2026_W45).

**Requirements**

- A Neuropixels 1.0 probe (384 channels, bank 0) with its AP (30 kHz) and LFP (2.5 kHz) streams, from the Neuropixels plugin or a recording played back with the File Reader.
- The localizer bundle (~40 MB, the model's predictions on the atlas): download `probe-viewer-localizer-2026_W45.zip` from the [releases page](https://github.com/alonsaguy/probe-viewer/releases) and unzip it.
- The planned trajectory of the insertion, in one of these formats:
  - an `ephysatlas-plan/1` JSON file, made with `ephysatlas.realtime.plan` in ibleatools (`plan_from_trajectory` for IBL-convention coordinates: entry point, depth, theta, phi);
  - a Pinpoint export (entry AP, ML, DV relative to bregma, insertion depth, theta, phi; µm or mm);
  - a text file with at least two `x y z` points of the track (IBL coordinates: ML, AP, DV from bregma; µm or mm), the deepest one being the tip.

**Usage**

1. In the plugin editor, click **Bundle...** and choose the unzipped bundle folder (or set the `PROBE_VIEWER_LOCALIZER_BUNDLE` environment variable to it).
2. Click **Plan...** and choose the planned trajectory file.
3. Select the probe's AP or LFP stream as the displayed stream, start acquisition and click **Localize**. The status line shows the snippets as they are collected and localized; **Reset** forgets them (e.g. after moving the probe).

Processing a snippet takes 7-10 s on a 22-core PC; snippets are taken back to back whenever the previous one is done, so the estimate updates every 12-15 s. It is about as good after the first snippet as after five minutes.

**Accuracy** on 196 held-out IBL insertions (features of the IBL release, compared with the histology-based alignment; medians over insertions):

| | channel distance to histology | Cosmos region accuracy | Beryl region accuracy |
|---|---|---|---|
| real-time localizer (27 features) | 536 µm | 0.74 | 0.35 |
| planned trajectory | 754 µm | 0.57 | 0.21 |

The estimate is closer to the histology than the plan for 76% of the insertions. On held-out test recordings the plugin's C++ features match those of the IBL Python pipeline (correlation ≥ 0.995 for every feature, the same bad channels).

**Limitations**

- Neuropixels 1.0 only, one probe at a time.
- The search starts from the plan: plans off by up to ~0.25 mm or 10° cost nothing, but a plan 1 mm or more off is only partly corrected.
- For Pinpoint and xyz files the channels are placed along the planned track as given, with no Needles-frame correction; for IBL (Needles-convention) coordinates make the plan JSON with `plan_from_trajectory`.
- Validated on recordings played back from disk; it has not yet been tested during a live experiment.

A self-test checks the C++ localizer against the Python one on the example stored in the bundle: configure with `-DPROBE_VIEWER_BUILD_SELFTEST=ON` and run `probe-viewer-localizer-selftest <bundle folder>`.

## Building from source

First, follow the instructions on [this page](https://open-ephys.github.io/gui-docs/Developer-Guide/Compiling-the-GUI.html) to build the Open Ephys GUI.

**Important:** This plugin is intended for use with the latest version of the GUI (0.6.0 and higher). The GUI should be compiled from the [`main`](https://github.com/open-ephys/plugin-gui/tree/main) branch, rather than the former `master` branch.

Then, clone this repository into a directory at the same level as the `plugin-GUI`, e.g.:
 
```
Code
├── plugin-GUI
│   ├── Build
│   ├── Source
│   └── ...
├── OEPlugins
│   └── probe-viewer
│       ├── Build
│       ├── Source
│       └── ...
```

### Windows

**Requirements:** [Visual Studio](https://visualstudio.microsoft.com/) and [CMake](https://cmake.org/install/)

From the `Build` directory, enter:

```bash
cmake -G "Visual Studio 17 2022" -A x64 ..
```

Next, launch Visual Studio and open the `OE_PLUGIN_probe-viewer.sln` file that was just created. Select the appropriate configuration (Debug/Release) and build the solution.

Selecting the `INSTALL` project and manually building it will copy the `.dll` and any other required files into the GUI's `plugins` directory. The next time you launch the GUI from Visual Studio, the Probe Viewer plugin should be available.


### Linux

**Requirements:** [CMake](https://cmake.org/install/)

From the `Build` directory, enter:

```bash
cmake -G "Unix Makefiles" ..
cd Debug
make -j
make install
```

This will build the plugin and copy the `.so` file into the GUI's `plugins` directory. The next time you launch the compiled version of the GUI, the Probe Viewer plugin should be available.


### macOS

**Requirements:** [Xcode](https://developer.apple.com/xcode/) and [CMake](https://cmake.org/install/)

From the `Build` directory, enter:

```bash
cmake -G "Xcode" ..
```

Next, launch Xcode and open the `probe-viewer.xcodeproj` file that now lives in the “Build” directory.

Running the `ALL_BUILD` scheme will compile the plugin; running the `INSTALL` scheme will install the `.bundle` file to `/Users/<username>/Library/Application Support/open-ephys/plugins-api`. The Probe Viewer plugin should be available the next time you launch the GUI from Xcode.