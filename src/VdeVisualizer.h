// VDE (vcg-direct) visualizer.  Owns the per-collapse snapshot renderer and
// the live-callback factory; main_cli holds the instance for the vcg-direct
// branch.

#pragma once

#if defined(QMAT_WITH_POLYSCOPE) && defined(QMAT_WITH_VCGLIB)

#include "MatVisualizerCommon.h"
#include "VcgDirectSimplifier.h"   // for VcgDirectSnapshot, LiveUpdateCallback

class SlabMesh;

class VdeVisualizer : public MatVisualizer {
public:
    // Runs SetupSimplificationViewer, then clears sm.on_collapse_cb and
    // sets vs_.vcg_direct_active so the shared ImGui handler branches right.
    void Setup(SlabMesh& sm) override;

    LiveUpdateCallback MakeLiveCallback();

    // Draws snap into "MAT Faces"/"MAT Edges"/"MAT Verts" + struct overlays.
    void Render(const VcgDirectSnapshot& snap);
};

#endif  // QMAT_WITH_POLYSCOPE && QMAT_WITH_VCGLIB
