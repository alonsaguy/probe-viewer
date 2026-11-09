#include "PlanFile.hpp"

#include "MiniJson.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <regex>
#include <sstream>

namespace ProbeViewer::Localization::Realtime
{
namespace
{
constexpr double kPi = 3.14159265358979323846;

std::string readText (const std::string& path)
{
    std::ifstream f (path, std::ios::binary);
    if (! f)
        throw std::runtime_error ("cannot open " + path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

std::string lower (std::string s)
{
    std::transform (s.begin(), s.end(), s.begin(), [] (unsigned char c) { return static_cast<char> (std::tolower (c)); });
    return s;
}

// The first number labelled by one of the labels ("ap: 1.2", "\"AP\" = -2100", "ap (um), 300").
bool labeledNumber (const std::string& text, const std::vector<std::string>& labels, double& value)
{
    for (const auto& label : labels)
    {
        const std::regex pattern ("(?:^|[^a-z0-9_])" + label
                                  + "[\"']?\\s*(?:\\([^)]*\\))?\\s*(?:[:=]|,)?\\s*(-?[0-9]+(?:\\.[0-9]+)?)");
        std::smatch match;
        if (std::regex_search (text, match, pattern) && match.size() >= 2)
        {
            value = std::stod (match[1].str());
            return true;
        }
    }
    return false;
}

// Pinpoint/Urchin convention: phi (yaw) = 0 and theta (pitch) = 0 point anterior, flat in the
// horizontal plane; increasing pitch tilts toward the brain. x ML, y AP, z DV (up positive).
PlanPoint pinpointAxis (double thetaDeg, double phiDeg)
{
    const double theta = thetaDeg * kPi / 180.0, phi = phiDeg * kPi / 180.0;
    return { std::sin (phi) * std::cos (theta), std::cos (phi) * std::cos (theta), -std::sin (theta) };
}

// Entry and tip (µm) of a Pinpoint export, or false.
bool parsePinpoint (const std::string& text, PlanPoint& entry, PlanPoint& tip)
{
    double ap = 0, ml = 0, dv = 0, depth = 0, theta = 0, phi = 0;
    if (! (labeledNumber (text, { "ap", "anterior[-_ ]?posterior" }, ap)
           && labeledNumber (text, { "ml", "medial[-_ ]?lateral" }, ml)
           && labeledNumber (text, { "dv", "dorsal[-_ ]?ventral" }, dv)
           && labeledNumber (text, { "depth", "insertion[-_ ]?depth" }, depth)
           && labeledNumber (text, { "theta", "elevation" }, theta)
           && labeledNumber (text, { "phi", "azimuth" }, phi)))
        return false;
    const bool mm = text.find ("mm") != std::string::npos;
    const bool um = text.find ("um") != std::string::npos || text.find ("micron") != std::string::npos
                 || text.find ("\xc2\xb5m") != std::string::npos;
    const double largest = std::max ({ std::abs (ap), std::abs (ml), std::abs (dv), std::abs (depth) });
    const double scale = mm || (! um && largest < 50.0) ? 1000.0 : 1.0;
    entry = { ml * scale, ap * scale, dv * scale };
    const auto axis = pinpointAxis (theta, phi);
    for (int k = 0; k < 3; ++k)
        tip[static_cast<size_t> (k)] = entry[static_cast<size_t> (k)] + axis[static_cast<size_t> (k)] * depth * scale;
    return true;
}

// Points of an "x y z" per line trajectory file (µm, or mm when every |value| < 50).
std::vector<PlanPoint> parseXyz (const std::string& text)
{
    std::vector<PlanPoint> points;
    std::istringstream lines (text);
    std::string line;
    while (std::getline (lines, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::replace_if (line.begin(), line.end(), [] (char c) { return c == ',' || c == ';' || c == '\t'; }, ' ');
        std::istringstream tokens (line);
        PlanPoint p {};
        if (tokens >> p[0] >> p[1] >> p[2])
            points.push_back (p);
    }
    double largest = 0.0;
    for (const auto& p : points)
        largest = std::max ({ largest, std::abs (p[0]), std::abs (p[1]), std::abs (p[2]) });
    if (! points.empty() && largest < 50.0)
        for (auto& p : points)
            for (auto& v : p)
                v *= 1000.0;
    return points;
}

// A Pinpoint v5 experiment file (AllenNeuralDynamics/pinpoint).
bool isPinpointV5 (const Json& doc)
{
    return doc.type == Json::Type::Object && doc.has ("probes") && doc.has ("probeInterfaceProbes") && doc.has ("atlas");
}

// Planned channel positions (m, IBL-Allen, top first) of one probe of a Pinpoint v5 experiment.
// Pinpoint keeps the tip in atlas ASR mm (AP, DV, ML from the atlas origin) and the rotation as
// [roll, yaw, pitch] (rad; Babylon RotationYawPitchRoll), under which probe-local +z, up the
// shank from the tip, is (cos yaw cos pitch, -sin pitch, sin yaw cos pitch) in ASR. A contact
// sits its probeinterface y above the bottom of the probe contour, the tip. The allen_mouse atlas
// is the CCF, read from the IBL bregma; the qiu2018_mouse (MRI) atlas is read from its own bregma
// and scaled to the CCF as iblatlas' MRITorontoAtlas (ML 0.952, AP 1.031, DV 0.885).
std::vector<PlanPoint> pinpointV5Channels (const Json& doc, const Json& probe, int nChannels)
{
    const std::string atlas = doc["atlas"]["name"].str();
    std::vector<double> bregmaMm; // ASR
    PlanPoint scale { 1.0, 1.0, 1.0 }; // ML, AP, DV
    if (atlas.rfind ("allen_mouse", 0) == 0)
        bregmaMm = { 5.400, 0.332, 5.739 };
    else if (atlas.rfind ("qiu2018_mouse", 0) == 0)
    {
        const Json& manifest = doc["atlas"]["manifest"];
        bregmaMm = (manifest.has ("bregma") && manifest["bregma"].type == Json::Type::Array ? manifest["bregma"] : doc["referenceCoordinate"]).numbers();
        scale = { 0.952, 1.031, 0.885 };
    }
    else
        throw std::runtime_error ("its atlas is " + atlas + ": plan in the allen_mouse or qiu2018_mouse atlas");
    if (bregmaMm.size() != 3)
        throw std::runtime_error ("no bregma for the " + atlas + " atlas");

    const Json& definition = doc["probeInterfaceProbes"][probe["probeInterfaceIdentifier"].str()];
    if (definition.has ("annotations") && definition["annotations"].has ("model_name"))
    {
        const std::string model = definition["annotations"]["model_name"].str();
        if (model.rfind ("NP1", 0) != 0)
            throw std::runtime_error (probe["name"].str() + " is a " + model + " probe, the localizer needs a Neuropixels 1.0");
    }
    const std::string units = definition.has ("si_units") ? definition["si_units"].str() : "um";
    const double toMm = units == "m" ? 1000.0 : units == "mm" ? 1.0 : 1e-3;
    double tipY = std::numeric_limits<double>::infinity();
    for (const auto& p : definition["probe_planar_contour"].items)
        tipY = std::min (tipY, p[1].num());
    const auto& contacts = definition["contact_positions"].items;
    if (static_cast<int> (contacts.size()) < nChannels || ! std::isfinite (tipY))
        throw std::runtime_error (probe["name"].str() + " has " + std::to_string (contacts.size()) + " contacts, the recording "
                                  + std::to_string (nChannels) + " channels");

    const auto tip = probe["tipPosition"].numbers();
    const auto rotation = probe["rotation"].numbers();
    if (tip.size() != 3 || rotation.size() != 3)
        throw std::runtime_error ("bad tipPosition or rotation for " + probe["name"].str());
    const double yaw = rotation[1], pitch = rotation[2];
    const PlanPoint up { std::cos (yaw) * std::cos (pitch), -std::sin (pitch), std::sin (yaw) * std::cos (pitch) };

    std::vector<PlanPoint> rows (static_cast<size_t> (nChannels));
    for (int r = 0; r < nChannels; ++r)
    {
        const int channel = nChannels - 1 - r;
        const double height = (contacts[static_cast<size_t> (channel)][1].num() - tipY) * toMm;
        PlanPoint d {};
        for (size_t k = 0; k < 3; ++k)
            d[k] = tip[k] + height * up[k] - bregmaMm[k];
        rows[static_cast<size_t> (r)] = { d[2] / scale[0] * 1e-3, -d[0] / scale[1] * 1e-3, -d[1] / scale[2] * 1e-3 };
    }
    return rows;
}

// Channel c of a Neuropixels 1.0 probe sits 20 (c / 2) + 20 µm up the track from the tip.
std::vector<PlanPoint> np1Channels (const PlanPoint& entryUm, const PlanPoint& tipUm, int nChannels)
{
    PlanPoint u { entryUm[0] - tipUm[0], entryUm[1] - tipUm[1], entryUm[2] - tipUm[2] };
    const double norm = std::sqrt (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    for (auto& v : u)
        v /= norm;
    std::vector<PlanPoint> rows (static_cast<size_t> (nChannels));
    for (int r = 0; r < nChannels; ++r)
    {
        const int channel = nChannels - 1 - r;
        const double axial = 20.0 * (channel / 2) + 20.0;
        for (int k = 0; k < 3; ++k)
            rows[static_cast<size_t> (r)][static_cast<size_t> (k)] = (tipUm[static_cast<size_t> (k)] + axial * u[static_cast<size_t> (k)]) * 1e-6;
    }
    return rows;
}
} // namespace

bool loadPlanJson (const std::string& path, std::vector<PlanPoint>& plannedM, std::string& error)
{
    try
    {
        const Json doc = Json::parse (readText (path));
        if (! doc.has ("format") || doc["format"].str() != "ephysatlas-plan/1")
            throw std::runtime_error (path + " is not an ephysatlas-plan/1 file");
        plannedM.clear();
        for (const auto& row : doc["planned_xyz_m"].items)
            plannedM.push_back ({ row[0].num(), row[1].num(), row[2].num() });
        if (plannedM.empty())
            throw std::runtime_error (path + " has no planned positions");
        return true;
    }
    catch (const std::exception& e)
    {
        error = e.what();
        return false;
    }
}

std::vector<std::string> pinpointProbeNames (const std::string& path)
{
    std::vector<std::string> names;
    try
    {
        const std::string text = readText (path);
        if (text.find ("\"probeInterfaceProbes\"") == std::string::npos)
            return names;
        const Json doc = Json::parse (text);
        if (isPinpointV5 (doc))
            for (const auto& probe : doc["probes"].items)
                names.push_back (probe["name"].str());
    }
    catch (const std::exception&)
    {
        names.clear();
    }
    return names;
}

bool loadPlanFile (const std::string& path, int nChannels, std::vector<PlanPoint>& plannedM,
                   std::string& source, std::string& error, const std::string& probeName)
{
    std::string text;
    try
    {
        text = readText (path);
    }
    catch (const std::exception& e)
    {
        error = e.what();
        return false;
    }
    if (text.find ("ephysatlas-plan/1") != std::string::npos)
    {
        if (! loadPlanJson (path, plannedM, error))
            return false;
        if (static_cast<int> (plannedM.size()) != nChannels)
        {
            error = "the plan has " + std::to_string (plannedM.size()) + " channels, the recording " + std::to_string (nChannels);
            return false;
        }
        source = "ephysatlas plan";
        return true;
    }
    if (text.find ("\"probeInterfaceProbes\"") != std::string::npos)
    {
        try
        {
            const Json doc = Json::parse (text);
            if (isPinpointV5 (doc))
            {
                const auto& probes = doc["probes"].items;
                const Json* chosen = nullptr;
                std::string names;
                for (const auto& probe : probes)
                {
                    const std::string& name = probe["name"].str();
                    names += (names.empty() ? "" : ", ") + name;
                    if (probeName.empty() ? probes.size() == 1 : name == probeName)
                        chosen = &probe;
                }
                if (chosen == nullptr)
                {
                    error = probes.empty() ? path + " has no probes"
                          : probeName.empty() ? path + " has " + std::to_string (probes.size()) + " probes (" + names + "): choose one with Plan..."
                                              : path + " has no probe " + probeName + " (" + names + ")";
                    return false;
                }
                plannedM = pinpointV5Channels (doc, *chosen, nChannels);
                source = "Pinpoint " + (doc.has ("version") ? doc["version"].str() : std::string ("v5")) + ", " + (*chosen)["name"].str();
                return true;
            }
        }
        catch (const std::exception& e)
        {
            error = path + " (Pinpoint v5): " + e.what();
            return false;
        }
    }
    const std::string lowered = lower (text);
    PlanPoint entry {}, tip {};
    if (parsePinpoint (lowered, entry, tip))
        source = "Pinpoint export";
    else
    {
        const auto points = parseXyz (text);
        if (points.size() < 2)
        {
            error = path + " is not an ephysatlas plan, a Pinpoint export or an xyz trajectory";
            return false;
        }
        const auto deepest = std::min_element (points.begin(), points.end(), [] (const PlanPoint& a, const PlanPoint& b) { return a[2] < b[2]; });
        const auto highest = std::max_element (points.begin(), points.end(), [] (const PlanPoint& a, const PlanPoint& b) { return a[2] < b[2]; });
        tip = *deepest;
        entry = *highest;
        source = "xyz trajectory";
    }
    const double dx = entry[0] - tip[0], dy = entry[1] - tip[1], dz = entry[2] - tip[2];
    if (std::sqrt (dx * dx + dy * dy + dz * dz) < 100.0)
    {
        error = "the planned track in " + path + " is shorter than 100 um";
        return false;
    }
    plannedM = np1Channels (entry, tip, nChannels);
    return true;
}
} // namespace ProbeViewer::Localization::Realtime
