#ifndef _QEM_REJECTION_VIZ_H
#define _QEM_REJECTION_VIZ_H

// ─────────────────────────────────────────────────────────────────────────────
// QemRejectionViz — Polyscope overlay for the VCG-path QEM collapse rejections.
//
// This is a small, self-contained module that the *existing* QMAT Polyscope
// viewer (main_cli.cpp) calls into.  It registers its structures into the SAME
// Polyscope scene/window as everything else (it is NOT a second viewer): an
// extra curve network "MAT QEM Rejection Edges" coloured by QemRejectionReason,
// plus per-reason highlight overlays (sliver triangle, flipped face, the
// non-manifold one-ring).  main_cli just calls UpdateEdgeColors / HandlePick /
// DrawPanel at the right spots, keeping the heavy code out of main_cli.cpp.
//
// Compiles to nothing unless both ONLY_USE_QEM_CONDITION_CHECKS (the VCG path)
// and QMAT_WITH_POLYSCOPE (the viewer) are defined.
// ─────────────────────────────────────────────────────────────────────────────

#if defined(ONLY_USE_QEM_CONDITION_CHECKS) && defined(QMAT_WITH_POLYSCOPE)

#include <vector>
#include <cstddef>

class SlabMesh;
namespace polyscope { class Structure; }   // forward decl — avoids pulling polyscope here

namespace qemviz {

// Pick-mapping + selection state for the "MAT QEM Rejection Edges" overlay.
// Plain data, no polyscope/imgui dependency, so it can be embedded directly in
// main_cli's ViewerState.
struct State {
	std::vector<unsigned> eid_order;       // pick edge slot → slab edge ID
	std::size_t           node_count = 0;  // pick offset (curve-network nodes precede edges)
	int                   selected_eid = -1;
	bool                  show_spheres = false;
};

// Register / refresh the QEM rejection curve network (each active edge coloured
// by its last QemRejectionReason; white = never rejected).  Preserves the
// enabled-state across refreshes.  Rebuilds st.eid_order / st.node_count.
void UpdateEdgeColors(const SlabMesh& sm, State& st);

// Highlight the geometry that explains slab edge `eid`'s rejection
// (sliver triangle / flipped face before+after / non-manifold one-ring).
void ShowPrimitives(const SlabMesh& sm, unsigned eid, bool show_spheres);

// Hide every QEM rejection overlay structure.
void ClearPrimitives();

// If `struct_ptr` is the QEM rejection curve network, consume the pick (update
// the selection and show its primitives) and return true; otherwise return
// false so the caller can keep handling the pick.  Takes the polyscope
// Structure* from pick::getSelection() (forward-declared, so this header stays
// free of polyscope headers).
bool HandlePick(const SlabMesh& sm, State& st, polyscope::Structure* struct_ptr, std::size_t local_idx);

// Draw the "QEM Rejection" section of the ImGui panel for the current selection
// (reason name, the free-text explanation, and any numeric metrics).
void DrawPanel(const SlabMesh& sm, State& st);

} // namespace qemviz

#endif  // ONLY_USE_QEM_CONDITION_CHECKS && QMAT_WITH_POLYSCOPE
#endif  // _QEM_REJECTION_VIZ_H
