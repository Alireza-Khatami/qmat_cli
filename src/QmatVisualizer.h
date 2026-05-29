// QMAT visualizer.  Owns the QMAT-specific helpers, ImGui panel install,
// and the per-collapse sm.on_collapse_cb wiring.

#pragma once

#ifdef QMAT_WITH_POLYSCOPE

#include "MatVisualizerCommon.h"

class SlabMesh;

class QmatVisualizer : public MatVisualizer {
public:
    void Setup(SlabMesh& sm) override;

    // Called from main() after Simplify() returns: clears the per-collapse
    // callback, rebuilds the live MAT structures + rejection viz + struct
    // overlays, and hides the per-collapse Collapsed Edge highlight.
    void RenderFinal(SlabMesh& sm);
};

// Shared scene + ImGui panel install.  Both QmatVisualizer::Setup and
// VdeVisualizer::Setup call this.  Does NOT touch sm.on_collapse_cb.
void InstallSharedScene(SlabMesh& sm, ViewerState& vs);

#endif  // QMAT_WITH_POLYSCOPE
