// VcgDirectSimplifier — see VcgDirectSimplifier.h.
//
// Recipe taken verbatim from vcglib/apps/tridecimator/tridecimator.cpp:
//   - minimal vcg::TriMesh<MyVertex, MyFace> with VFAdj + per-vertex Quadric
//   - vcg::tri::TriEdgeCollapseQuadric<MyMesh, BasicVertexPair, MyEdgeCollapse,
//                                    QInfoStandard<MyVertex>>
//   - vcg::LocalOptimization<MyMesh>::Init / DoOptimization / Finalize
//
// We deliberately keep all vcg includes INSIDE this .cpp so the rest of QMAT
// (which uses Wm4 / CGAL) doesn't see them.

#include "VcgDirectSimplifier.h"

#include "SlabMesh.h"
#include "Mesh.h"      // for pmesh->bb_diagonal_length

#include <algorithm>
#include <iostream>

// ── vcglib ──────────────────────────────────────────────────────────────────
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/update/normal.h>

// ─────────────────────────────────────────────────────────────────────────────
// Minimal vcg::TriMesh definition — same as vcglib/apps/tridecimator.
// Vertex stores the Quadric inline; QInfoStandard wraps Qd() so TECQ can find it.
// ─────────────────────────────────────────────────────────────────────────────
namespace {

class VDVertex;
class VDEdge;
class VDFace;

struct VDUsedTypes : public vcg::UsedTypes<
	vcg::Use<VDVertex>::AsVertexType,
	vcg::Use<VDEdge  >::AsEdgeType,
	vcg::Use<VDFace  >::AsFaceType> {};

class VDVertex : public vcg::Vertex<VDUsedTypes,
	vcg::vertex::VFAdj,
	vcg::vertex::Coord3f,
	vcg::vertex::Normal3f,
	vcg::vertex::Mark,
	vcg::vertex::Qualityf,
	vcg::vertex::BitFlags>
{
public:
	vcg::math::Quadric<double>& Qd() { return q; }
private:
	vcg::math::Quadric<double> q;
};

class VDEdge : public vcg::Edge<VDUsedTypes> {};

class VDFace : public vcg::Face<VDUsedTypes,
	vcg::face::VFAdj,
	vcg::face::VertexRef,
	vcg::face::Normal3f,
	vcg::face::BitFlags> {};

class VDMesh : public vcg::tri::TriMesh<std::vector<VDVertex>, std::vector<VDFace>> {};

typedef vcg::tri::BasicVertexPair<VDVertex> VDVertexPair;

// Per-collapse live-update plumbing.  vcg constructs VDEdgeCollapse instances
// internally and there's no per-instance state we control, so the callback and
// scratch buffers live as TU-level statics, set/cleared in RunVcgDirectSimplify.
LiveUpdateCallback                 g_live_cb;
std::vector<std::array<double, 3>> g_live_verts;
std::vector<std::array<int,    3>> g_live_faces;

void ExtractResult(VDMesh& mesh,
                   std::vector<std::array<double, 3>>& verts,
                   std::vector<std::array<int,    3>>& faces);

class VDEdgeCollapse : public vcg::tri::TriEdgeCollapseQuadric<
	VDMesh, VDVertexPair, VDEdgeCollapse, vcg::tri::QInfoStandard<VDVertex>>
{
public:
	typedef vcg::tri::TriEdgeCollapseQuadric<
		VDMesh, VDVertexPair, VDEdgeCollapse, vcg::tri::QInfoStandard<VDVertex>> TECQ;
	inline VDEdgeCollapse(const VDVertexPair& p, int i, vcg::BaseParameterClass* pp) : TECQ(p, i, pp) {}

	// LocalModification::Execute is virtual — override to fire the live hook
	// after the collapse has actually been applied to the mesh.
	void Execute(VDMesh& m, vcg::BaseParameterClass* pp)
	{
		TECQ::Execute(m, pp);
		if (g_live_cb)
		{
			ExtractResult(m, g_live_verts, g_live_faces);
			g_live_cb(g_live_verts, g_live_faces);
		}
	}
};

// Build VDMesh from the active vertices/faces of `sm`.  Coordinates are scaled
// by bb_diagonal_length so they line up with what _mat_initial.off / the
// existing visualizer use.  Triangle vertex order follows std::set iteration —
// the same order used by ExportMatAsOff in main_cli.cpp.
void BuildFromSlab(const SlabMesh& sm, VDMesh& mesh, int& nv_active, int& nf_tri)
{
	nv_active = 0;
	nf_tri    = 0;

	const double scale = sm.pmesh ? sm.pmesh->bb_diagonal_length : 1.0;

	// Slab vertex idx -> compact vcg vertex idx (-1 if vertex is deleted).
	std::vector<int> slabToVcg(sm.vertices.size(), -1);
	for (size_t i = 0; i < sm.vertices.size(); ++i)
		if (sm.vertices[i].first) ++nv_active;

	for (size_t i = 0; i < sm.faces.size(); ++i)
		if (sm.faces[i].first && sm.faces[i].second &&
		    sm.faces[i].second->vertices_.size() == 3)
			++nf_tri;

	if (nv_active == 0 || nf_tri == 0) return;

	auto vi = vcg::tri::Allocator<VDMesh>::AddVertices(mesh, nv_active);
	int written = 0;
	for (size_t i = 0; i < sm.vertices.size(); ++i)
	{
		if (!sm.vertices[i].first) continue;
		slabToVcg[i] = written++;
		const auto& c = sm.vertices[i].second->sphere.center;
		(*vi).P() = VDMesh::CoordType(
			(float)(c[0] * scale),
			(float)(c[1] * scale),
			(float)(c[2] * scale));
		++vi;
	}

	auto fi = vcg::tri::Allocator<VDMesh>::AddFaces(mesh, nf_tri);
	for (size_t i = 0; i < sm.faces.size(); ++i)
	{
		if (!sm.faces[i].first || !sm.faces[i].second) continue;
		const auto& vs = sm.faces[i].second->vertices_;
		if (vs.size() != 3) continue;
		auto it = vs.begin();
		int a = slabToVcg[*it++];
		int b = slabToVcg[*it++];
		int c = slabToVcg[*it];
		if (a < 0 || b < 0 || c < 0) continue;
		(*fi).V(0) = &mesh.vert[a];
		(*fi).V(1) = &mesh.vert[b];
		(*fi).V(2) = &mesh.vert[c];
		++fi;
	}

	// Topology bookkeeping vcg needs before LocalOptimization::Init.
	vcg::tri::UpdateBounding<VDMesh>::Box(mesh);
	vcg::tri::UpdateTopology<VDMesh>::VertexFace(mesh);
	vcg::tri::UpdateFlags<VDMesh>::FaceBorderFromVF(mesh);
	vcg::tri::UpdateNormal<VDMesh>::PerFaceNormalized(mesh);
}

// Pull surviving vertices/faces out of `mesh` (IsD() == false), remapping vcg
// indices to a compact 0..N-1 range.  Clears the output buffers first so the
// helper is safe to call repeatedly (e.g. from the per-collapse hook).
void ExtractResult(VDMesh& mesh,
                   std::vector<std::array<double, 3>>& verts,
                   std::vector<std::array<int,    3>>& faces)
{
	verts.clear();
	faces.clear();

	std::vector<int> vcgToOut(mesh.vert.size(), -1);
	int next = 0;
	verts.reserve(mesh.VN());
	for (size_t i = 0; i < mesh.vert.size(); ++i)
	{
		if (mesh.vert[i].IsD()) continue;
		vcgToOut[i] = next++;
		const auto& p = mesh.vert[i].P();
		verts.push_back({ (double)p[0], (double)p[1], (double)p[2] });
	}

	faces.reserve(mesh.FN());
	for (size_t i = 0; i < mesh.face.size(); ++i)
	{
		if (mesh.face[i].IsD()) continue;
		const size_t i0 = vcg::tri::Index(mesh, *mesh.face[i].V(0));
		const size_t i1 = vcg::tri::Index(mesh, *mesh.face[i].V(1));
		const size_t i2 = vcg::tri::Index(mesh, *mesh.face[i].V(2));
		faces.push_back({ vcgToOut[i0], vcgToOut[i1], vcgToOut[i2] });
	}
}

} // anonymous namespace

bool RunVcgDirectSimplify(const SlabMesh& sm,
                          const VcgDirectParams& params,
                          VcgDirectResult& out,
                          LiveUpdateCallback live_callback)
{
	VDMesh mesh;
	int nv_active = 0, nf_tri = 0;
	BuildFromSlab(sm, mesh, nv_active, nf_tri);

	out.initial_vertex_count = nv_active;
	out.initial_face_count   = nf_tri;

	if (mesh.FN() == 0)
	{
		std::cerr << "[VcgDirectSimplify] no triangular faces; nothing to do\n";
		return false;
	}

	// Map our params onto vcg's TriEdgeCollapseQuadricParameter.  Defaults
	// for unspecified knobs come from the struct itself, so anything we don't
	// touch matches vcg's built-in defaults.
	vcg::tri::TriEdgeCollapseQuadricParameter qparams;
	qparams.QualityThr             = params.QualityThr;
	qparams.BoundaryQuadricWeight  = params.BoundaryQuadricWeight * 0.5; // UI -> internal
	qparams.OptimalPlacement       = params.OptimalPlacement;
	qparams.PreserveTopology       = params.PreserveTopology;
	qparams.NormalCheck            = params.NormalCheck;
	if (qparams.NormalCheck)
		qparams.NormalThrRad = vcg::math::ToRad(45.0);   // MeshLab default (M_PI/4)

	// Resolve target face count.  Half the current face count is a reasonable
	// fallback when the caller passes neither a face nor a vertex target.
	int target_faces = params.TargetFaceNum;
	if (target_faces < 0 && params.TargetVertexNum > 0)
		target_faces = std::max(1, params.TargetVertexNum * 2);  // V≈F/2 for closed meshes
	if (target_faces < 0)
		target_faces = std::max(1, mesh.FN() / 2);
	target_faces = std::min(target_faces, mesh.FN());

	std::cerr << "[VcgDirectSimplify] V=" << mesh.VN() << " F=" << mesh.FN()
	          << " -> target faces " << target_faces
	          << " | OptimalPlacement=" << qparams.OptimalPlacement
	          << " PreserveTopology="   << qparams.PreserveTopology
	          << " NormalCheck="        << qparams.NormalCheck
	          << " QualityThr="         << qparams.QualityThr
	          << " BoundaryQuadricWeight=" << qparams.BoundaryQuadricWeight
	          << "\n";

	// Install the per-collapse hook for VDEdgeCollapse::Execute to call.
	g_live_cb = std::move(live_callback);

	vcg::LocalOptimization<VDMesh> DeciSession(mesh, &qparams);
	DeciSession.Init<VDEdgeCollapse>();
	std::cerr << "[VcgDirectSimplify] initial heap size " << DeciSession.h.size() << "\n";

	DeciSession.SetTargetSimplices(target_faces);
	DeciSession.SetTimeBudget(0.5f);

	while (DeciSession.DoOptimization() && mesh.FN() > target_faces) { /* iterate */ }

	DeciSession.Finalize<VDEdgeCollapse>();
	out.collapses_performed = DeciSession.nPerformedOps;

	// Clear the hook so we don't dangle into the next call.
	g_live_cb = nullptr;

	ExtractResult(mesh, out.vertices, out.faces);
	out.final_vertex_count = (int)out.vertices.size();
	out.final_face_count   = (int)out.faces.size();

	std::cerr << "[VcgDirectSimplify] done V=" << out.final_vertex_count
	          << " F=" << out.final_face_count
	          << " collapses=" << out.collapses_performed << "\n";
	return true;
}
