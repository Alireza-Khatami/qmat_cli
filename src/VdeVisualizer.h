// VDE (vcg-direct) visualizer.  Owns the per-collapse snapshot renderer and
// the live-callback factory; main_cli holds the instance for the vcg-direct
// branch.

#pragma once

#ifdef QMAT_WITH_VCGLIB

#include <string>

#include "VcgDirectSimplifier.h"   // for VcgDirectSnapshot, LiveUpdateCallback

// ── VDE-only CLI-side exporters ──────────────────────────────────────────
// All take a const VcgDirectSnapshot& and reconstruct the corresponding
// QMAT output file from snapshot data alone (no SlabMesh access needed).
// Defined in VdeVisualizer.cpp's unguarded section so they link in
// non-polyscope builds too.

// OFF mesh: snapshot.vertices + snapshot.faces.
void ExportSnapshotAsOff(const VcgDirectSnapshot& snap, const std::string& path);

// Post-simplification visualize_info JSON, same schema as
// ExportSimpVisualizeInfo but sourced from the snapshot.
void ExportSnapshotVisualizeInfo(const VcgDirectSnapshot& snap, const std::string& path);

// .mat_typed: per-vertex "v x y z r MS_<type_name>" + "e a b" + "f v0 v1 v2".
void ExportSnapshotMatTyped(const VcgDirectSnapshot& snap, const std::string& path);

// .ma: "n_v n_e n_f\n" + "v x y z r" + "e a b" + "f v0 v1 v2".  The user
// passes a prefix; this routine appends "___v_N___e_N___f_N.ma" to match
// SlabMesh::Export's filename convention.
void ExportSnapshotMa(const VcgDirectSnapshot& snap, const std::string& path_prefix);

// _cluster.ply: positions + per-vertex RGB by cluster_type + faces + edges.
void ExportSnapshotClusterPLY(const VcgDirectSnapshot& snap, const std::string& path);

// _rejection_skeleton.ply: cylinders per active edge, coloured by
// edge_last_rejection.  radius_frac = cylinder radius as a fraction of the
// snapshot's bbox diagonal (matches SlabMesh::ExportSkeletonPLY default).
void ExportSnapshotRejectionSkeleton(const VcgDirectSnapshot& snap,
                                     const std::string& path,
                                     double radius_frac = 0.001);

#endif  // QMAT_WITH_VCGLIB

#if defined(QMAT_WITH_POLYSCOPE) && defined(QMAT_WITH_VCGLIB)

#include "MatVisualizerCommon.h"

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
