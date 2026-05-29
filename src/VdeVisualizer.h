// VDE (vcg-direct) visualizer.  Owns the per-collapse snapshot renderer and
// the live-callback factory; main_cli holds the instance for the vcg-direct
// branch.

#pragma once

#if defined(QMAT_WITH_POLYSCOPE) && defined(QMAT_WITH_VCGLIB)

#include "MatVisualizerCommon.h"
#include "VcgDirectSimplifier.h"   // for VcgDirectSnapshot, LiveUpdateCallback

class SlabMesh;

class VdeVisualizer {
public:
    // Reads vcg_direct_active = true, sets sm.on_collapse_cb = nullptr.
    // SetupSimplificationViewer (still in main_cli) must have run first.
    void Bind(SlabMesh& sm);

    LiveUpdateCallback MakeLiveCallback();

    // Draws snap into "MAT Faces"/"MAT Edges"/"MAT Verts" + struct overlays.
    void Render(const VcgDirectSnapshot& snap);

    ViewerState&       State()       { return vs_; }
    const ViewerState& State() const { return vs_; }

private:
    ViewerState vs_;
};

#endif  // QMAT_WITH_POLYSCOPE && QMAT_WITH_VCGLIB
