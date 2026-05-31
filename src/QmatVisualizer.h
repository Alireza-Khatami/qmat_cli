// QMAT visualizer.  Owns the QMAT-specific helpers, ImGui panel install,
// and the per-collapse sm.on_collapse_cb wiring.

#pragma once

#include <string>

class SlabMesh;

// ── QMAT-only CLI-side exporter ──────────────────────────────────────────
// Post-simplification visualize_info JSON: positions + connectivity +
// cluster_type + topo_type + struct_ids + rejection reasons + ancestry +
// id->name/colour legends.  Must be called BEFORE SlabMesh::Export() (that
// calls AdjustStorage(), which compacts storage and invalidates the
// edge_last_rejection map).  Defined in QmatVisualizer.cpp's unguarded
// section so it links in CLI-only builds too.
void ExportSimpVisualizeInfo(const SlabMesh& sm, const std::string& path);

#ifdef QMAT_WITH_POLYSCOPE

#include "MatVisualizerCommon.h"

class QmatVisualizer : public MatVisualizer {
public:
    void Setup(SlabMesh& sm) override;

    // Called from main() after Simplify() returns: clears the per-collapse
    // callback, rebuilds the live MAT structures + rejection viz + struct
    // overlays, and hides the per-collapse Collapsed Edge highlight.
    void RenderFinal(SlabMesh& sm);
};

#endif  // QMAT_WITH_POLYSCOPE
