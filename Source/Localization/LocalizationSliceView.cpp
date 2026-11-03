#include "LocalizationSliceView.hpp"

#include <algorithm>
#include <cmath>

namespace ProbeViewer
{
namespace
{
using P = std::array<float, 3>;
using AtlasSliceRef = Localization::AtlasSlice;

float dot (const P& a, const P& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

P sub (const P& a, const P& b)
{
    return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
}
}

LocalizationSliceView::LocalizationSliceView() = default;

LocalizationSliceView::~LocalizationSliceView() = default;

juce::String LocalizationSliceView::resultKey (
    const Localization::LocalizationResult& result)
{
    if (result.probe.id.isNotEmpty())
        return result.probe.id;
    if (result.probe.displayName.isNotEmpty())
        return result.probe.displayName;
    if (result.streamKey.isNotEmpty())
        return result.streamKey;
    return "probe";
}

void LocalizationSliceView::setLocalizationResult (
    Localization::LocalizationResultPtr result)
{
    if (result == nullptr || ! result->valid)
        return;

    resultsByProbe[resultKey (*result)] = result;
    referenceResult = result;
    repaint();
}

void LocalizationSliceView::clear()
{
    resultsByProbe.clear();
    referenceResult.reset();
    repaint();
}

void LocalizationSliceView::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colour::fromRGB (18, 18, 18));

    if (referenceResult == nullptr)
    {
        g.setColour (juce::Colours::lightgrey);
        g.setFont (16.0f);
        g.drawFittedText (
            "Run localization to display the atlas slice.",
            getLocalBounds().reduced (20),
            juce::Justification::centred,
            2);
        return;
    }

    const auto& atlasSlice = referenceResult->atlasSlice;

    if (! atlasSlice.valid || atlasSlice.voxels.empty())
    {
        g.setColour (juce::Colours::orange);
        g.setFont (16.0f);
        g.drawFittedText (
            "Localization completed, but the atlas slice is empty.",
            getLocalBounds().reduced (20),
            juce::Justification::centred,
            2);
        return;
    }

    juce::String title ("Localization slice");
    if (referenceResult->probe.displayName.isNotEmpty())
    {
        title += " - ";
        title += referenceResult->probe.displayName;
    }
    title += "  |  init: ";
    title += referenceResult->initializationSource;

    g.setColour (juce::Colours::white);
    g.setFont (juce::FontOptions (16.0f).withStyle ("Bold"));
    g.drawText (title, 14, 7, getWidth() - 28, 22,
        juce::Justification::centredLeft);

    auto available = getLocalBounds().reduced (14).toFloat();
    available.removeFromTop (28.0f);
    available.removeFromBottom (22.0f);

    // One slice in a 2:1 box, or two slices side by side (tilted coronal and sagittal).
    auto drawSlice = [&] (const AtlasSliceRef& slice, juce::Rectangle<float> box, const juce::String& label)
    {
        const float worldWidthUm = slice.maxUUm - slice.minUUm;
        const float worldHeightUm = slice.maxVUm - slice.minVUm;
        const float pixelsPerUm = std::min (
            box.getWidth() / std::max (worldWidthUm, 1.0f),
            box.getHeight() / std::max (worldHeightUm, 1.0f));
        const float usedWidth = worldWidthUm * pixelsPerUm;
        const float usedHeight = worldHeightUm * pixelsPerUm;
        juce::Rectangle<float> plotBounds (usedWidth, usedHeight);
        plotBounds.setCentre (box.getCentre());
        const float x0 = plotBounds.getX();
        const float y0 = plotBounds.getY();

        auto mapPoint = [&] (float uUm, float vUm)
        {
            return juce::Point<float> {
                x0 + (uUm - slice.minUUm) * pixelsPerUm,
                y0 + (slice.maxVUm - vUm) * pixelsPerUm
            };
        };

        const float voxelPixels = std::max (1.0f, slice.voxelSizeUm * pixelsPerUm);
        for (const auto& voxel : slice.voxels)
        {
            const auto p = mapPoint (voxel.uUm, voxel.vUm);
            g.setColour (voxel.regionId > 0 ? voxel.colour : juce::Colours::black);
            g.fillRect (
                p.x - 0.5f * voxelPixels,
                p.y - 0.5f * voxelPixels,
                voxelPixels + 0.5f,
                voxelPixels + 0.5f);
        }

        auto drawTrajectory = [&] (
            const std::vector<P>& xyz,
            juce::Colour colour,
            bool dashed,
            float thickness)
        {
            if (xyz.size() < 2)
                return;
            juce::Path path;
            bool started = false;
            for (const auto& point : xyz)
            {
                const auto relative = sub (point, slice.originXyzUm);
                const float u = dot (relative, slice.basisUXyz);
                const float v = dot (relative, slice.basisVXyz);
                const auto p = mapPoint (u, v);
                if (! started)
                {
                    path.startNewSubPath (p);
                    started = true;
                }
                else
                {
                    path.lineTo (p);
                }
            }
            g.setColour (colour);
            if (dashed)
            {
                const float dashLengths[] = { 7.0f, 5.0f };
                juce::Path dashedPath;
                juce::PathStrokeType (thickness).createDashedStroke (
                    dashedPath, path, dashLengths, 2);
                g.fillPath (dashedPath);
            }
            else
            {
                g.strokePath (path, juce::PathStrokeType (thickness));
            }
        };

        int colourIndex = 0;
        for (const auto& item : resultsByProbe)
        {
            const auto& result = item.second;
            if (result == nullptr)
                continue;
            const auto predictionColour = juce::Colour::fromHSV (
                std::fmod (0.08f + 0.23f * static_cast<float> (colourIndex), 1.0f),
                0.85f, 1.0f, 1.0f);
            drawTrajectory (result->initialTrajectoryXyzUm, juce::Colours::white.withAlpha (0.75f), true, 2.0f);
            drawTrajectory (result->predictedTrajectoryXyzUm, predictionColour, false, 3.0f);
            ++colourIndex;
        }

        g.setColour (juce::Colours::white.withAlpha (0.85f));
        g.drawRect (plotBounds, 1.0f);
        if (label.isNotEmpty())
        {
            g.setFont (12.0f);
            g.drawText (label, plotBounds.reduced (4.0f).removeFromTop (16.0f), juce::Justification::centredLeft);
        }
        return plotBounds;
    };

    const auto& secondary = referenceResult->secondarySlice;
    juce::Rectangle<float> legendRow;
    if (secondary.valid && ! secondary.voxels.empty())
    {
        auto left = available.removeFromLeft (available.getWidth() / 2.0f).reduced (4.0f, 0.0f);
        auto right = available.reduced (4.0f, 0.0f);
        const auto a = drawSlice (atlasSlice, left, "coronal (tilted along the probe), ML horizontal");
        const auto b = drawSlice (secondary, right, "sagittal (tilted along the probe), AP horizontal");
        legendRow = a.getUnion (b);
    }
    else
    {
        const float plotHeight = std::min (available.getHeight(), available.getWidth() / 2.0f);
        juce::Rectangle<float> box (2.0f * plotHeight, plotHeight);
        box.setCentre (available.getCentre());
        legendRow = drawSlice (atlasSlice, box, {});
    }

    g.setColour (juce::Colours::white.withAlpha (0.85f));
    g.setFont (12.0f);
    g.drawText (
        "dashed: initialization   solid: prediction",
        static_cast<int> (legendRow.getX()),
        static_cast<int> (legendRow.getBottom() + 4.0f),
        static_cast<int> (legendRow.getWidth()),
        18,
        juce::Justification::centredLeft);
}
}
