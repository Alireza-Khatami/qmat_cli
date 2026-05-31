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
    // Initialises polyscope, registers the input mesh + the live MAT
    // (overwritten on the first collapse), installs the VDE ImGui panel,
    // and clears sm.on_collapse_cb (vcg-direct doesn't drive the slab cb).
    void Setup(SlabMesh& sm) override;

    LiveUpdateCallback MakeLiveCallback();

    // Draws snap into "MAT Faces"/"MAT Edges"/"MAT Verts" + struct overlays.
    void Render(const VcgDirectSnapshot& snap);

    // Latest snapshot, cached so the ImGui panel's pick handler can resolve
    // a clicked "MAT Verts" idx to ancestor positions.
    VcgDirectSnapshot latest_snap;
};

#endif  // QMAT_WITH_POLYSCOPE && QMAT_WITH_VCGLIB
