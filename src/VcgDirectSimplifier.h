// VcgDirectSimplifier — drive vcg's REAL TriEdgeCollapseQuadric on the medial
// mesh, bypassing QMAT's port (VcgQuadricSimplifier).  The point is to compare
// against / fall back to the canonical MeshLab implementation without copying
// any code: vcglib is pulled via FetchContent and we use exactly the recipe
// from vcglib/apps/tridecimator/tridecimator.cpp.
//
// Only enabled when CMake configures with -DQMAT_WITH_VCGLIB=ON (default).
// Selected at runtime via main_cli's --simplifier vcg-direct flag.
//
// IMPORTANT: this file MUST NOT include any vcg headers (those drag Eigen
// through with macros that conflict with our Wm4 / CGAL world).  All vcg types
// stay hidden inside VcgDirectSimplifier.cpp.

#pragma once

#include <array>
#include <functional>
#include <vector>

class SlabMesh;

// Fires immediately after each collapse executes, with the current surviving
// mesh state (compact 0..N-1 indices, same layout as VcgDirectResult).
using LiveUpdateCallback = std::function<void(
	const std::vector<std::array<double, 3>>& vertices,
	const std::vector<std::array<int,    3>>& faces)>;

struct VcgDirectParams
{
	// One of these stops the run.  TargetFaceNum wins when both are set.
	// -1 means "unset"; if both are -1, target = FN()/2 (half the faces).
	int    TargetFaceNum         = -1;
	int    TargetVertexNum       = -1;

	// MeshLab QEM defaults — see meshlabplugins/filter_meshing/meshfilter.cpp
	// for the GUI defaults and tri_edge_collapse_quadric.h for the internal
	// struct defaults (which DIFFER from the GUI: see project_meshlab_qem_defaults).
	double QualityThr               = 0.3;   // sliver veto threshold
	double BoundaryQuadricWeight    = 1.0;   // MeshLab UI slider; we multiply by 0.5
	                                         // before assigning to vcg (matches
	                                         // meshfilter.cpp:14 effective weight).
	bool   OptimalPlacement         = true;  // "Optimal position" checkbox
	bool   PreserveTopology         = true;  // "Preserve topology" checkbox
	                                         // (vcg LinkConditions gate)
	bool   NormalCheck              = true;  // "Preserve normal" checkbox
};

struct VcgDirectResult
{
	// Simplified mesh in compact 0..N-1 index order — directly registrable in
	// polyscope via registerSurfaceMesh(vertices, faces).  Vertices are in the
	// same world scale used to build the vcg mesh (sphere.center * bb_diagonal
	// _length), so polyscope coordinates line up with "MAT Faces".
	std::vector<std::array<double, 3>> vertices;
	std::vector<std::array<int,    3>> faces;

	int initial_vertex_count = 0;
	int initial_face_count   = 0;
	int final_vertex_count   = 0;
	int final_face_count     = 0;
	int collapses_performed  = 0;
};

// Build a vcg::TriMesh from `sm` (active vertices/faces only, scaled the same
// way as _mat_initial.off), run vcg's TriEdgeCollapseQuadric driven by `params`,
// and write the result into `out`.  Returns false if the slab mesh has no
// triangular faces (nothing to do).
//
// If `live_callback` is non-empty it fires once per executed collapse with the
// current surviving mesh state (compact 0..N-1 indices).  Default is no-op.
bool RunVcgDirectSimplify(const SlabMesh& sm,
                          const VcgDirectParams& params,
                          VcgDirectResult& out,
                          LiveUpdateCallback live_callback = {});
