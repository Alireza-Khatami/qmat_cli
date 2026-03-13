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

	// Cluster input mesh faces by surface orientation using face-adjacency BFS.
	// Two-pass BFS on the face-edge-adjacency graph:
	//   Pass 1 (strict, angle_threshold_deg=25°): separates flat face patches.
	//   Pass 2 (loose, edge_merge_deg=60°): merges remaining unlabeled faces
	//           (edge/corner faces) without touching already-labeled patches.
	// Neighbors of a face = all faces sharing an edge with it.
	// Requires compute_normals() to have been called first.
	// Returns vector-of-vectors of vertex indices, sorted largest-first.
	//   Vertices shared between clusters are assigned to the larger cluster.
	// Writes <meshname>_boundary_clusters.txt for external visualisation.
	std::vector<std::vector<unsigned>> ClusterBoundaryPoints(
		double angle_threshold_deg = 25.0,
		double edge_merge_deg      = 60.0);

public:
	Mesh input;		// the mesh of the input shape

	unsigned num_vor_v, num_vor_e, num_vor_f;

	NonManifoldMesh input_nmm;

	SlabMesh slab_mesh;

	bool slab_initial;

	// ── boundary clustering results (populated by ClusterBoundaryPoints) ──────
	// Each inner vector holds the input mesh vertex indices belonging to that cluster.
	// Sorted largest-first.
	std::vector<std::vector<unsigned>> boundary_clusters;

	// Maps each input mesh vertex index → its cluster id in boundary_clusters.
	// -1 if the vertex was not assigned to any cluster.
	std::vector<int> vertex_cluster_id;

	// Voronoi neighbor graph of the boundary points.
	// voronoi_neighbors[bp_id] = set of boundary point IDs whose Voronoi cells
	// share a face with bp_id's cell (equivalently: connected by a Delaunay edge).
	// Populated by ComputeInputNMM() and loaded by LoadInputNMM().
	std::vector<std::set<unsigned>> voronoi_neighbors;

};
#endif // _THREE_DIMENSIONAL_SHAPE_H_