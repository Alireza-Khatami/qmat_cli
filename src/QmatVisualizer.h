// QMAT visualizer.  Wraps the slab-side SetupSimplificationViewer (initial
// register, ImGui panel, picking, per-collapse callback).  Phase 3a skeleton:
// the helpers it depends on still live in main_cli.cpp; later sub-phases
// will move them into the .cpp file.

#pragma once

#ifdef QMAT_WITH_POLYSCOPE

#include "MatVisualizerCommon.h"

class SlabMesh;

class QmatVisualizer : public MatVisualizer {
public:
    void Setup(SlabMesh& sm) override;
};

#endif  // QMAT_WITH_POLYSCOPE
