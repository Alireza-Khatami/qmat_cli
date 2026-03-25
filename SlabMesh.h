#ifndef _SLABMESH_H
#define _SLABMESH_H

#include "PrimMesh.h"
#include "MatCollapseHistory.h"
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
		T0           = 0,  // 0 clusters — impossible (safety sentinel)
		T1_spike     = 1,  // 1 cluster + exactly 4 bpoints — true spike (Delaunay tetrahedron)
		T2           = 2,  // 2 components
		T3           = 3,  // 3 components
		T4           = 4,  // exactly 4 components
		T5           = 5,  // >4 clusters — impossible (safety sentinel)
		T1_non_spike = 6   // 1 cluster + >4 bpoints — sheet/flat boundary, not a spike
	};
	ClusterType nmn_cluster_type = ClusterType::T0;

	// Topology type derived from the four topo_is_* flags set by DetermineTopology().
	// Priority: junction > seam > boundary > sheet.
	enum class TopoType : uint8_t {
		Unknown  = 0,  // no incident active edges (isolated vertex)
		Sheet    = 1,  // only 2-manifold edges incident
		Boundary = 2,  // >= 1 boundary edge (nf==1), no seam edge
		Seam     = 3,  // >= 1 seam edge (nf>2), no boundary edge
		Junction = 4,  // >= 3 seam edges (nf>2) incident
	};
	TopoType topo_type = TopoType::Unknown;

	// Maximum face-count observed across all incident MAT edges during
	// DetermineTopology().  Useful for debugging: shows the "strongest" edge
	// type incident to this vertex (1=boundary, 2=sheet, >2=seam).
	signed nf = 0;

	SlabVertex() : is_spike(false),
	               topo_is_sheet(false), topo_is_seam(false),
	               topo_is_junction(false), topo_is_boundary(false),
	               nmn_cluster_type(ClusterType::T0),
	               topo_type(TopoType::Unknown),
	               nf(0) {}

};

class SlabEdge : public PrimEdge, public SlabPrim
{
};

class SlabFace : public PrimFace, public SlabPrim
{
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

	// Output prefix used for exported files (set from main_cli after loading).
	// e.g. "bear/bear" → files written as "bear/bear_post_spike.off" etc.
	std::string export_prefix;

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
	bool Contractible(unsigned vid_src1, unsigned vid_src2, Vector3d &v_tgt);
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
	void LabelVertices();   // assign topo_label on every active vertex

	// Recomputes topo_is_sheet/seam/junction/boundary and is_spike
	// on every active SlabVertex by inspecting the face-valence of each incident
	// edge. Safe to call after simplification to refresh stale conservative flags.
	void DetermineTopology();

	// Recomputes topology flags and topo_type for a single vertex by inspecting
	// only its incident edges.  Call after MergeVertices() on the new vid_tgt.
	void RecomputeVertexTopology(unsigned vid);

	// Computes nmn_bplist_clusters for every active vertex using union-find
	// on input mesh edges.  Call once after LoadInputNMM.
	void ClusterNMNBplist();

	bool CanMerge(unsigned vid1, unsigned vid2) const;
	bool CanMerge(unsigned vid1, unsigned vid2, CollapseContext ctx) const;

	// Returns true if collapsing the edge (vid0, vid1) would produce a
	// non-manifold configuration — i.e. the collapse should be rejected.
	// Translated from the PMP is_collapse_ok() link-condition check.
	bool WouldCreateNonManifold(unsigned vid0, unsigned vid1) const;

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