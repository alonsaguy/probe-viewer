#pragma once

// The planned trajectory the real-time localizer starts from, as planned channel positions.
// Positions in metres, IBL-Allen coordinates (origin bregma; x ML, y AP, z DV), row order: top
// of the probe first, row r is channel C - 1 - r.

#include <array>
#include <string>
#include <vector>

namespace ProbeViewer::Localization::Realtime
{
using PlanPoint = std::array<double, 3>;

// An ephysatlas-plan/1 file (ephysatlas.realtime.plan.save_plan): the planned channel positions.
bool loadPlanJson (const std::string& path, std::vector<PlanPoint>& plannedM, std::string& error);

// Any supported plan file, with nChannels channels of a Neuropixels 1.0 probe:
// - ephysatlas-plan/1 JSON, used as is (its channel count must be nChannels);
// - a Pinpoint export: entry AP, ML, DV (from bregma), insertion depth, theta (pitch from
//   horizontal) and phi (yaw, 0 = anterior), in µm or mm; the tip is depth along that axis;
// - a trajectory text file: at least two "x y z" lines (µm or mm, IBL x ML, y AP, z DV), the
//   deepest point taken as the tip.
// For the last two the channels are placed along the track from the tip at the NP1 axial
// positions (20 µm rows), with no Needles-frame correction: for an IBL (Needles convention)
// trajectory make the JSON with ephysatlas.realtime.plan.plan_from_trajectory instead.
// source receives a description of what was read.
bool loadPlanFile (const std::string& path, int nChannels, std::vector<PlanPoint>& plannedM,
                   std::string& source, std::string& error);
} // namespace ProbeViewer::Localization::Realtime
