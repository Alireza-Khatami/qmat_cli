#ifndef _SLABMESH_H
#define _SLABMESH_H

#include "PrimMesh.h"
#include "MatCollapseHistory.h"
#include <unordered_map>
#include <array>
#include <vector>
#include <optional>
#ifdef QMAT_WITH_POLYSCOPE
#  include <functional>
#endif

class SlabPrim
{
public:
	Wm4::Matrix4d slab_A;
	Wm4::Vector4d slab_b;
	double slab_c;

	Wm4::Matrix4d add_A;
	Wm4::Vector4d add_b;
	double add_c;

	// hyperbolic weight
	double hyperbolic_weight;

	SlabPrim() : slab_c(0.0), add_c(0.0), hyperbolic_weight(0.0){}
};

class SlabVertex : public PrimVertex, public SlabPrim
{
public:
	bool is_pole;
	bool is_non_manifold;
	bool is_disk;
	bool is_boundary;

	// Transferred from NonManifoldMesh_Vertex via _mat_topo.txt sidecar.
	// ── topology flags ────────────────────────────────────────────────────────
	bool is_spike;  // all bplist boundary points form a mesh-edge clique
	bool topo_is_sheet;         // all incident MAT edges are 2-manifold
	bool topo_is_seam;          // >= 1 incident MAT edge shared by > 2 faces
	bool topo_is_junction;      // > 2 seam edges converge here
	bool topo_is_boundary;      // >= 1 incident MAT edge with exactly 1 face

	// ── boundary point information (needed for simplification) ────────────────
	// Indices of the input mesh boundary points (from the Delaunay tetrahedron)
	// that gave rise to this MAT vertex. Named nmn_bplist to distinguish from
	// PrimVertex::bplist which is used for Hausdorff tracking.
	std::set<unsigned> nmn_bplist;

	// Connected components of nmn_bplist restricted to input mesh edges.
	// Each entry is one cluster (set of bp indices reachable from each other
	// via mesh edges).  Populated by SlabMesh::ClusterNMNBplist().
	std::vector<std::set<unsigned>> nmn_bplist_clusters;

	// Cluster-count category derived from nmn_bplist_clusters.size().
	// T0 and T5 should never occur — they flag a logic bug if seen.
	enum class ClusterType : uint8_t {
		// ── Voronoi/Delaunay path (computed from bp cluster count) ───────────
		T0           = 0,  // 0 clusters — impossible (safety sentinel)
		T1_spike     = 1,  // 1 cluster + exactly 4 bpoints — true spike (Delaunay tetrahedron)
		T2           = 2,  // 2 components
		T3           = 3,  // 3 components
		T4           = 4,  // exactly 4 components
		T5           = 5,  // >4 clusters — impossible (safety sentinel)
		T1_non_spike = 6,  // 1 cluster + >4 bpoints — sheet/flat boundary, not a spike
		// ── MatStruct path (loaded from external .ma file) ──────────────────
		// Maps from: enum MedialType { MUNKNOWN=-1, SHEET=0, SEAM=1, BOUNDARY=2, JUNCTION=3 }
		MS_Unknown           = 7,  // not classified
		MS_Sheet             = 8,  // SHEET only (interior 2-manifold)
		MS_Seam              = 9,  // SEAM only (internal feature / crease)
		MS_Boundary          = 10, // BOUNDARY only (not in any sheet/seam/junction struct)
		MS_Junction          = 11, // JUNCTION only
		// Compound: base type + boundary modifier
		MS_Sheet_Boundary    = 12, // SHEET + BOUNDARY
		MS_Seam_Boundary     = 13, // SEAM  + BOUNDARY
		MS_Junction_Boundary = 14, // JUNCTION + BOUNDARY
	};
	ClusterType nmn_cluster_type = ClusterType::T0;

	// Struct IDs assigned by LoadMatstructMA (empty = not part of any named struct).
	// A vertex can belong to multiple structure blocks (e.g. a junction vertex shared
	// by two different named structures). Inherited by child vertices after collapse.
	std::set<int> struct_ids;

	// Maximum face-count observed across all incident MAT edges during
	// DetermineTopology().  Useful for debugging: shows the "strongest" edge
	// type incident to this vertex (1=boundary, 2=sheet, >2=seam).
	signed nf = 0;

	// Set by MarkSharpFeatureVertices() before simplification.
	// True if this vertex sits at a turning angle exceeding the feature-angle
	// threshold along its seam or boundary chain.  Such vertices are excluded
	// from all collapses to preserve the original sharp feature geometry.
	bool sharpNotContractable = false;

	// Lineage: the set of original (pre-simplification) MAT vertex IDs that
	// collapsed into this vertex.  Seeded with {index} when the vertex is first
	// created (LoadInputNMM / LoadMatstructMA / InsertVertex).  On every
	// MergeVertices(v1, v2, v_new): original_ancestors(v_new) = union of
	// original_ancestors(v1) and original_ancestors(v2).  After simplification,
	// the union across all surviving vertices partitions the set of original
	// vertex IDs (every original ID appears in exactly one survivor's set).
	std::set<unsigned> original_ancestors;

	SlabVertex() : is_spike(false),
	               topo_is_sheet(false), topo_is_seam(false),
	               topo_is_junction(false), topo_is_boundary(false),
	               nmn_cluster_type(ClusterType::T0),
	               nf(0),
	               sharpNotContractable(false) {}

};

class SlabEdge : public PrimEdge, public SlabPrim
{
public:
    // Struct IDs assigned by LoadMatstructMA (empty = not part of any named struct).
    // An edge can belong to multiple structure blocks. Inherited (union) on collapse.
    std::set<int> struct_ids;

    // Topological classification of this MAT edge, derived from the imported
    // structure data during LoadMatstructMA — NOT from face-count geometry,
    // which can drift after collapses.  Source rules:
    //   - listed in a SEAM block      → Seam
    //   - listed in a BOUNDARY block  → Boundary
    //   - listed in BOTH              → Seam_Boundary
    //   - listed in NEITHER, but every incident face is in a SHEET struct → Sheet
    //   - no incident faces           → Orphan
    //   - otherwise                   → Unknown
    // Inherited across collapses by flag-union (see CombineTopoType): Seam-ness
    // and Boundary-ness OR together; Sheet survives only when neither stronger
    // flag is raised.
    enum class TopoType : uint8_t {
        Unknown       = 0,
        Sheet         = 1,
        Seam          = 2,
        Boundary      = 3,
        Seam_Boundary = 4,
        Orphan        = 5,
    };
    TopoType topo_type = TopoType::Unknown;

    // Union two edge topo types under flag semantics: seam-ness and boundary-
    // ness OR; sheet survives only if neither stronger flag is set.  Used both
    // during LoadMatstructMA (when the same edge appears in multiple blocks)
    // and during collapse (when two parent edges fold into one new edge).
    static TopoType CombineTopoType(TopoType a, TopoType b)
    {
        using TT = TopoType;
        auto is_seam = [](TT t) { return t == TT::Seam     || t == TT::Seam_Boundary; };
        auto is_bnd  = [](TT t) { return t == TT::Boundary || t == TT::Seam_Boundary; };
        auto is_sht  = [](TT t) { return t == TT::Sheet; };
        const bool s = is_seam(a) || is_seam(b);
        const bool b_flag = is_bnd(a) || is_bnd(b);
        const bool h = is_sht(a) || is_sht(b);
        if (s && b_flag) return TT::Seam_Boundary;
        if (s)           return TT::Seam;
        if (b_flag)      return TT::Boundary;
        if (h)           return TT::Sheet;
        if (a == TT::Orphan && b == TT::Orphan) return TT::Orphan;
        return TT::Unknown;
    }

    static const char* TopoTypeName(TopoType t)
    {
        switch (t) {
            case TopoType::Unknown:       return "Unknown";
            case TopoType::Sheet:         return "Sheet";
            case TopoType::Seam:          return "Seam";
            case TopoType::Boundary:      return "Boundary";
            case TopoType::Seam_Boundary: return "Seam_Boundary";
            case TopoType::Orphan:        return "Orphan";
        }
        return "Unknown";
    }
};

class SlabFace : public PrimFace, public SlabPrim
{
public:
    // Struct ID assigned by LoadMatstructMA (-1 = not part of any named struct).
    int struct_id = -1;
};

typedef std::pair<bool, SlabVertex*> Bool_SlabVertexPointer;
typedef std::pair<bool, SlabEdge*> Bool_SlabEdgePointer;
typedef std::pair<bool, SlabFace*> Bool_SlabFacePointer;

class SlabMesh : public PrimMesh
{
public:
	std::vector<Bool_SlabVertexPointer> vertices;
	std::vector<Bool_SlabEdgePointer> edges;
	std::vector<Bool_SlabFacePointer> faces;

public:
	// 1. preserve method one
	// 2. preserve method two
	// 3. preserve method three
	int preserve_boundary_method;

	// 0. no hyperbolic weight
	// 1. add hyperbolic distance to the related edges
	// 2. add hyperbolic area to the related face
	// 3. add ratio of hyperbolic and Euclid to the related edges
	int hyperbolic_weight_type;

	// 1.compute the boundary vertices only
	// 2.compute the boundary vertices and its related vertices
	int boundary_compute_scale;

	// the influence factor to control the ratio between hyperbolic and Euclid distance
	double k;

	// whether clear the error before initialization
	bool clear_error;

	bool preserve_saved_vertices;

	bool compute_hausdorff;

	bool prevent_inversion;

	// When true, CanMerge allows collapsing an edge where exactly one endpoint
	// is a steep tetrahedron vertex (bypasses the sheet check for that vertex).
	// After the merge the result vertex is always marked non-steep.
	bool allow_steep_collapse;

	// Turning-angle threshold (degrees) used by:
	//   MarkSharpFeatureVertices() — pre-marks sharp chain vertices
	//   WouldExceedCurvatureThreshold() — runtime check during collapse
	// Set from --feature-angle CLI option before Simplify() is called.
	double feature_angle_threshold = 30.0;

	bool initial_boundary_preserve;

	double m_min[3];
	double m_max[3];

	unsigned simplified_inside_edges;
	unsigned simplified_boundary_edges;

	double bound_weight;

	// Voronoi neighbor graph of the input boundary points.
	// voronoi_neighbors[bp_id] = the set of boundary point IDs whose Voronoi
	// cells share a face with bp_id (i.e. connected by a Delaunay edge).
	// Sized to the total number of boundary points; empty set = isolated point.
	// Populated by ThreeDimensionalShape::LoadInputNMM (load path) and
	// ThreeDimensionalShape::ComputeInputNMM (compute path).
	std::vector<std::set<unsigned>> voronoi_neighbors;

	// Feature edge sets on the input surface mesh.
	// Populated by ThreeDimensionalShape::ComputeFeatureEdges() which is
	// called from both ComputeInputNMM and LoadInputNMM.
	// Each entry is a sorted (min_id, max_id) pair of input mesh vertex ids.
	std::set<std::array<int,2>> sharp_edges;
	std::set<std::array<int,2>> concave_edges;
	std::set<int> feature_corners;

	// Snapshot of pre-simplification vertex positions, keyed by original vertex id
	// (index in the vector == the id stored inside SlabVertex::original_ancestors).
	// Populated at load time by LoadInputNMM / LoadMatstructMA.  Needed because
	// DeleteVertex frees vertices[vid].second on collapse; this snapshot is the
	// only reliable way to look up an ancestor's position after simplification.
	std::vector<std::array<double,3>> original_positions;

	// Output prefix used for exported files (set from main_cli after loading).
	// e.g. "bear/bear" → files written as "bear/bear_post_spike.off" etc.
	std::string export_prefix;
	std::string current_phase; // set before each Simplify phase; used to name per-phase rejection logs

	// Full collapse history: merge tree + keyframe snapshots.
	// Populated during Simplify() / MinCostEdgeCollapse().
	MatCollapseHistory history;

public:
	void AdjustStorage();

public:
	bool ValidVertex(unsigned vid);
	bool Edge(unsigned vid0, unsigned vid1, unsigned & eid);
	bool Face(const std::set<unsigned> & vset, unsigned & fid);
	void UpdateCentroid(unsigned fid);
	void ComputeFacesCentroid();
	void UpdateNormal(unsigned fid);
	void ComputeFacesNormal();
	void UpdateVertexNormal(unsigned vid);
	void ComputeVerticesNormal();
	void GetNeighborVertices(unsigned vid, std::set<unsigned> & neighborvertices);
	void GetLinkedEdges(unsigned eid, std::set<unsigned> & neighboredges);
	void GetAdjacentFaces(unsigned fid, std::set<unsigned> & neighborfaces);
	bool Contractible(unsigned vid_src, unsigned vid_tgt);
	bool Contractible(unsigned vid_src1, unsigned vid_src2, const Vector3d &v_tgt,
	                  std::array<std::array<std::array<double,3>,3>,2>* out_flipped_face = nullptr);
	bool MergeVertices(unsigned vid_src1, unsigned vid_src2, unsigned &vid_tgt);

	unsigned VertexIncidentEdgeCount(unsigned vid);
	unsigned VertexIncidentFaceCount(unsigned vid);
	unsigned EdgeIncidentFaceCount(unsigned eid);


public:
	void DeleteFace(unsigned fid);
	void DeleteEdge(unsigned eid);
	void DeleteVertex(unsigned vid);

	void InsertVertex(SlabVertex *vertex, unsigned &vid);
	void InsertEdge(unsigned vid0, unsigned vid1, unsigned & eid);
	void InsertFace(std::set<unsigned> vset);

	void ComputeEdgeCone(unsigned eid);
	void ComputeEdgesCone();

	void ComputeVertexProperty(unsigned vid);
	void ComputeVerticesProperty();
	void ComputeFaceSimpleTriangles(unsigned fid);
	void ComputeFacesSimpleTriangles();

public:
	void initBoundaryCollapseQueue();
	void initCollapseQueue();

	// Populate spike_collapse_queue with all edges that have at least one T1 endpoint.
	// Call after ClusterNMNBplist() + DetermineTopology().
	void initSpikeCollapseQueue();

	// Write the current active MAT vertices and faces to an OFF file.
	void ExportOff(const std::string& path) const;
	void Simplify(int threshold);
	void SimplifyBoudary(int threshold);
	bool MinCostBoundaryEdgeCollapse(unsigned & eid);
	// Context passed to CanMerge/MinCostEdgeCollapse to select the merge policy.
	enum class CollapseContext {
		Spike,     // T1 endpoint — bypass boundary/topology blocks, force collapse
		Boundary,  // boundary edge collapse queue
		Main       // regular simplification queue
	};
	bool MinCostEdgeCollapse(unsigned& eid, CollapseContext ctx = CollapseContext::Main);
	void EvaluateEdgeCollapseCost(unsigned eid);
	void EvaluateEdgeHausdorffCost(unsigned eid);
	void ReEvaluateEdgeHausdorffCost(unsigned eid);

public: 
	void DistinguishVertexType();
	unsigned GetSavedPointNumber();
	unsigned GetConnectPointNumber();
	void InsertSavedPoint(unsigned vid);
	double NearestPoint(Vector3d point, unsigned vid);

public:
	void PreservBoundaryMethodOne();
	void PreservBoundaryMethodTwo();
	void PreservBoundaryMethodThree();
	void PreservBoundaryMethodFour();
	void PreservBoundaryMethodFive();

	void GetEnvelopeSet(const Vector4d & lamder, const set<unsigned> & neighbor_v, const set< std::set<unsigned> > & adj_faces, vector<Sphere> & sph_vec, vector<Cone> & con_vec, vector<SimpleTriangle> & st_vec);
	double EvaluateVertexDistanceErrorEnvelope(Vector4d & lamder, set<unsigned> & neighbor_vertices, set< set<unsigned> > & neighbor_faces, set<unsigned> & bplist);

	double GetHyperbolicLength(unsigned eid);
	double GetRatioHyperbolicEuclid(unsigned eid);

	void ExportSimplifyResult();
	void Export(string fname);

public:
	void clear();
	void RecomputerVertexType();
	void computebb();

	void CleanIsolatedVertices();
	void InitialTopologyProperty(unsigned vid);
	void InitialTopologyProperty();

	// Recomputes topo_is_sheet/seam/junction/boundary and is_spike
	// on every active SlabVertex by inspecting the face-valence of each incident
	// edge. Safe to call after simplification to refresh stale conservative flags.
	// Disabled: depended on removed topo_type field; kept commented for reference.
	// void DetermineTopology();

	// Recomputes topology flags and topo_type for a single vertex by inspecting
	// only its incident edges.  Call after MergeVertices() on the new vid_tgt.
	// Disabled: depended on removed topo_type field; kept commented for reference.
	// void RecomputeVertexTopology(unsigned vid);

	// Computes nmn_bplist_clusters for every active vertex using union-find
	// on input mesh edges.  Call once after LoadInputNMM.
	void ClusterNMNBplist();

	bool CanMerge(unsigned vid1, unsigned vid2, CollapseContext ctx) const;

	// Reasons an edge collapse can be rejected (used for logging).
	enum class RejectionReason : uint8_t {
		StaleEdge,                      // edge has been deleted since insertion
		InvalidVertex,                  // one or both endpoints deleted
		DifferentTopoType,              // CanMerge condition 1
		DifferentClusterType,           // CanMerge condition 2
		BplistNotNeighbors,             // CanMerge condition 3
		NoPmesh,                        // pmesh is null
		TopoNotContractable,            // edge flagged topo_contractable=false
		InversionWouldOccur,            // no valid sphere position avoids inversion
		// Non-manifold sub-reasons (one per test in WouldCreateNonManifold):
		NonManifold_BoundaryEdgePair,   // test 1/2: both side-edges of triangle are boundary
		NonManifold_SharedThirdVert,    // test 3:   2-face edge shares same third vertex (vl==vr)
		NonManifold_BoundaryVertEdge,   // test 4:   both endpoints boundary but shared edge is not
		NonManifold_LinkCondition,      // test 5:   one-ring intersection beyond safe third verts
		WouldCreateFoldOver,            // geometric edge crossing (fold-over / polyline self-intersection)
		SharpNotContractable,           // vertex marked as sharp feature by MarkSharpFeatureVertices()
		WouldExceedCurvatureThreshold,  // post-collapse turning angle would exceed feature_angle_threshold
		struct_ids_sets_different,        // struct edge: edge->struct_ids != both endpoints' struct_ids
		BoundaryHole,                     // edge is part of a triangular boundary hole
	};

	// Returns the RGB colour (0-255 per channel) that represents a rejection
	// reason in both ExportSkeletonPLY and the Polyscope live viewer.
	// Single source of truth — never duplicate this table.

	

	


	static std::array<uint8_t,3> RejectionReasonColorU8(RejectionReason rr)
	{
		using RR = RejectionReason;
		switch (rr) {
			// ── Less critical — greys ─────────────────────────────────────────
			case RR::StaleEdge:            return { 200, 200, 200 }; // light grey
			case RR::InvalidVertex:        return { 120, 120, 120 }; // dark grey
			case RR::DifferentTopoType:    return { 180, 180, 180 }; // mid grey
			case RR::DifferentClusterType: return { 255, 140,   0 }; // ORANGE
			case RR::BplistNotNeighbors:   return { 130, 130, 130 }; // mid grey
			case RR::NoPmesh:              return { 100, 100, 100 }; // dark grey
			case RR::InversionWouldOccur:  return {   0,   0,   0 }; // BLACK
			// ── Important — pure rainbow, easy to distinguish ─────────────────
			case RR::TopoNotContractable:           return { 255,   0,   0 }; // RED
			case RR::NonManifold_BoundaryEdgePair:  return { 255, 215,   0 }; // GOLD
			case RR::NonManifold_SharedThirdVert:   return { 255, 255,   0 }; // YELLOW
			case RR::NonManifold_BoundaryVertEdge:  return {   0, 255,   0 }; // GREEN
			case RR::NonManifold_LinkCondition:     return {   0, 255, 255 }; // CYAN
			case RR::WouldCreateFoldOver:           return {   0,   0, 255 }; // BLUE
			case RR::SharpNotContractable:          return { 148,   0, 211 }; // VIOLET
			case RR::WouldExceedCurvatureThreshold:      return { 255,   0, 255 }; // MAGENTA
			case RR::struct_ids_sets_different:          return { 165,  42,  42 }; // BROWN
			case RR::BoundaryHole:                       return { 255, 105, 180 }; // HOT PINK
			// ── Default — white = never attempted ────────────────────────────
			default:                                     return { 255, 255, 255 };
		}
	}
	
//   ┌─────────────────────────────────────────────┬──────────────────────────────┐
//   │                   Reason                    │           Color              │
//   ├─────────────────────────────────────────────┼──────────────────────────────┤
//   │ InversionWouldOccur                         │ ⚫ BLACK 0,0,0               │
//   ├─────────────────────────────────────────────┼──────────────────────────────┤
//   │ TopoNotContractable                         │ 🔴 RED 255,0,0               │
//   ├─────────────────────────────────────────────┼──────────────────────────────┤
//   │ DifferentClusterType                        │ 🟠 ORANGE 255,140,0          │
//   ├─────────────────────────────────────────────┼──────────────────────────────┤
//   │ NonManifold_BoundaryEdgePair                │ 🟡 GOLD 255,215,0            │
//   ├─────────────────────────────────────────────┼──────────────────────────────┤
//   │ NonManifold_SharedThirdVert                 │ 🟡 YELLOW 255,255,0          │
//   ├─────────────────────────────────────────────┼──────────────────────────────┤
//   │ NonManifold_BoundaryVertEdge                │ 🟢 GREEN 0,255,0             │
//   ├─────────────────────────────────────────────┼──────────────────────────────┤
//   │ NonManifold_LinkCondition                   │ 🩵 CYAN 0,255,255            │
//   ├─────────────────────────────────────────────┼──────────────────────────────┤
//   │ WouldCreateFoldOver                         │ 🔵 BLUE 0,0,255              │
//   ├─────────────────────────────────────────────┼──────────────────────────────┤
//   │ SharpNotContractable                        │ 🟣 VIOLET 148,0,211          │
//   ├─────────────────────────────────────────────┼──────────────────────────────┤
//   │ WouldExceedCurvatureThreshold               │ 🟣 MAGENTA 255,0,255         │
//   ├─────────────────────────────────────────────┼──────────────────────────────┤
//   │ StaleEdge                                   │ ⚪ Light grey 200,200,200    │
//   ├─────────────────────────────────────────────┼──────────────────────────────┤
//   │ DifferentTopoType                           │ 🔘 Mid grey 180,180,180      │
//   ├─────────────────────────────────────────────┼──────────────────────────────┤
//   │ BplistNotNeighbors                          │ 🔘 Mid grey 130,130,130      │
//   ├─────────────────────────────────────────────┼──────────────────────────────┤
//   │ InvalidVertex                               │ ⚫ Dark grey 120,120,120     │
//   ├─────────────────────────────────────────────┼──────────────────────────────┤
//   │ NoPmesh                                     │ ⚫ Dark grey 100,100,100     │
//   ├─────────────────────────────────────────────┼──────────────────────────────┤
//   │ Never attempted                             │ ⬜ WHITE 255,255,255         │
//   └─────────────────────────────────────────────┴──────────────────────────────┘

	// Slab-mesh primitives that directly caused a collapse rejection.
	// Stored alongside edge_last_rejection so the Polyscope viewer can highlight
	// the offending geometry when the user picks a rejected edge.
	struct ReasonPrimitives {
		std::vector<unsigned>             vertices;  // slab vertex IDs
		std::vector<std::array<unsigned,2>> edges;   // [v0, v1] slab vertex ID pairs
		std::vector<std::array<unsigned,3>> faces;   // [v0, v1, v2] slab vertex ID triples
		// World-space position the collapse would have targeted. nullopt = not set.
		std::optional<std::array<double,3>> targ_ver;
		// World-space positions of the triangle whose normal flipped (InversionWouldOccur).
		// [0] = triangle BEFORE collapse (original positions), [1] = triangle AFTER collapse.
		std::optional<std::array<std::array<std::array<double,3>,3>,2>> flipped_face;
	};

	// Returns true if collapsing the edge (vid0, vid1) would produce a
	// non-manifold configuration — i.e. the collapse should be rejected.
	// Translated from the PMP is_collapse_ok() link-condition check.
	bool WouldCreateNonManifold(unsigned vid0, unsigned vid1,
	                            RejectionReason*  out_reason = nullptr,
	                            ReasonPrimitives* out_prims  = nullptr) const;

	// Returns true if collapsing edge (vid0,vid1) to position v_tgt would cause
	// any resulting new edge to geometrically cross an existing edge (fold-over).
	// Uses a signed-volume straddling test — works for skew 3D segments.
	bool WouldCreateFoldOver(unsigned vid0, unsigned vid1,
	                         const Wm4::Vector3d& v_tgt,
	                         ReasonPrimitives* out_prims = nullptr) const;

	// Returns true if the edge (vid0,vid1) should be rejected because it is NOT
	// on a straight portion of a boundary/seam chain.
	// For each endpoint, finds its one other same-type chain neighbour and checks
	// whether the turning angle at that endpoint exceeds feature_angle_threshold.
	// Reject if either endpoint is bent, has no other same-type neighbour on the
	// far side, or has multiple (junction-like — inherently sharp).
	bool WouldExceedCurvatureThreshold(unsigned vid0, unsigned vid1,
	                                   ReasonPrimitives* out_prims = nullptr) const;
	
	// Pre-simplification pass: marks feature vertices as sharpNotContractable.
	//   MS_Junction: always marked (branch points are inherently sharp).
	//   MS_Boundary / MS_Seam: marked if >= 2 same-type neighbours AND turning
	//     angle along the chain exceeds angle_deg_threshold degrees.
	// Call once after DetermineTopology(), before Simplify().
	void MarkSharpFeatureVertices(double angle_deg_threshold = 30.0);
	void ComputeStructureCollapsibility();

	bool CanMerge(unsigned vid1, unsigned vid2,
	              RejectionReason*  out_reason = nullptr,
	              ReasonPrimitives* out_prims  = nullptr) const;

	// Appends one rejection record to {export_prefix}_rejection_log.txt.
	// Also stores prims in edge_reason_primitives[eid].
	// queue_name should be "spike", "boundary", or "edge".
	void LogCollapseRejection(const char* queue_name,
	                          unsigned eid, unsigned v1, unsigned v2,
	                          double cost, RejectionReason reason,
	                          ReasonPrimitives prims = {}) const;

	// Per-edge last rejection reason recorded during Simplify().
	// Mutable so LogCollapseRejection (const) can update it.
	mutable std::unordered_map<unsigned, RejectionReason>    edge_last_rejection;
	// Per-edge primitives that caused the last rejection.
	mutable std::unordered_map<unsigned, ReasonPrimitives>   edge_reason_primitives;

	// Export remaining active edges as cylinders to a PLY file (ASCII, per-vertex RGB).
	// Each cylinder is coloured by the last rejection reason recorded for that edge.
	// radius controls the cylinder radius (default 0.001 world-space units).
	void ExportSkeletonPLY(const std::string& path, double radius = 0.001) const;

	// Export the current MAT as a .mat_typed file: same format as .ma (header +
	// v/e/f lines with denormalized coords) but each vertex line has the
	// ClusterType name appended as a trailing token, e.g. "v x y z r MS_Seam".
	// Indices are compacted on-the-fly (no AdjustStorage side-effects).
	void ExportTypedMA(const std::string& path) const;

	// Export the current MAT as a PLY with per-vertex RGB colours matching the
	// cluster-type colour palette used in the Polyscope visualiser.
	// Includes vertex, face, and edge elements so the full MAT topology is visible.
	void ExportClusterPLY(const std::string& path) const;

#ifdef QMAT_WITH_POLYSCOPE
public:
	// Called just before each accepted edge collapse with the pre-merge sphere
	// data.  Set from main_cli to drive live Polyscope visualization.
	// Signature: (v1_id, v1_pos, v1_r, v2_id, v2_pos, v2_r, result_sphere)
	std::function<void(unsigned, const Wm4::Vector3d&, double,
	                   unsigned, const Wm4::Vector3d&, double,
	                   const Sphere&)> on_collapse_cb;
#endif
};

#endif