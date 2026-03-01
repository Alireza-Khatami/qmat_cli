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

	// Cluster input mesh vertices (boundary sample points) by position + normal
	// using CGAL::cluster_point_set (region-growing on a k-NN graph).
	// Requires CGAL >= 5.3 and that compute_normals() has been called first.
	//   k_neighbors : size of the neighbourhood graph (default 12)
	//   smoothness  : balance position vs. normal (0 = position only; default 0.5)
	// Returns vector-of-vectors: result[cluster_id] = { vertex indices into pVertexList }
	// Also writes <outputPrefix>_boundary_clusters.txt for external visualisation.
	std::vector<std::vector<unsigned>> ClusterBoundaryPoints(
		int    k_neighbors = 12);

public:
	Mesh input;		// the mesh of the input shape

	unsigned num_vor_v, num_vor_e, num_vor_f;

	NonManifoldMesh input_nmm;

	SlabMesh slab_mesh;

	bool slab_initial;

};
#endif // _THREE_DIMENSIONAL_SHAPE_H_