#ifndef _THREE_DIMENSIONAL_SHAPE_H_
#define _THREE_DIMENSIONAL_SHAPE_H_

#include "Mesh.h"
#include "NonManifoldMesh/nonmanifoldmesh.h"
#include "SlabMesh.h"

class ThreeDimensionalShape
{
public:
	ThreeDimensionalShape() : slab_initial(false) {}

	void ComputeInputNMM();
	
	// load the user simplified ma
	void LoadInputNMM(std::string fname);

#ifdef USE_MATSTRUCT_INITIALIZATION
	// Load initial MAT from an externally-generated typed .ma file.
	// Format: header "nv ne nf", then vertices "v x y z r T", edges "e v0 v1",
	// faces "f v0 v1 v2".  T is a MedialType integer (-1=unknown, 0=sheet,
	// 1=seam, 2=boundary, 3=junction) and is mapped to SlabVertex::ClusterType MS_*
	// values.  No sidecar bplist or Voronoi-neighbor files are loaded.
	void LoadMatstructMA(std::string fname);
#endif

	long LoadSlabMesh();

	// initial the matrix of each face and vertex for slab mesh
	void InitialSlabMesh();	

	double NearestPoint(Vector3d point, unsigned vid);

	// compute the Hausdorff distance
	void ComputeHausdorffDistance();	

	void PruningSlabMesh();

	// Classify every active MAT vertex topologically by inspecting its incident
	// edges and the number of faces each edge belongs to.
	// Sets topo_is_sheet / topo_is_seam / topo_is_junction / topo_is_boundary
	// on each NonManifoldMesh_Vertex.
	void DetermineTopology();

	// Detect sharp/concave edges and corner vertices on the input surface mesh
	// using dihedral-angle thresholds (MATFP method).
	// Results stored in slab_mesh.sharp_edges / concave_edges / feature_corners.
	void ComputeFeatureEdges();

	// Cluster input mesh faces by surface orientation using face-adjacency BFS.
	// Two-pass BFS on the face-edge-adjacency graph:
	//   Pass 1 (strict, angle_threshold_deg=25°): separates flat face patches.
	//   Pass 2 (loose, edge_merge_deg=60°): merges remaining unlabeled faces

public:
	Mesh input;		// the mesh of the input shape

	unsigned num_vor_v, num_vor_e, num_vor_f;

	NonManifoldMesh input_nmm;

	SlabMesh slab_mesh;

	bool slab_initial;

	// Voronoi neighbor graph of the boundary points.
	// voronoi_neighbors[bp_id] = set of boundary point IDs whose Voronoi cells
	// share a face with bp_id's cell (equivalently: connected by a Delaunay edge).
	// Populated by ComputeInputNMM() and loaded by LoadInputNMM().
	std::vector<std::set<unsigned>> voronoi_neighbors;

};
#endif // _THREE_DIMENSIONAL_SHAPE_H_