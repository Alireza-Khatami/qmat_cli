#pragma once

#include "Mesh.h"
#include <set>
#include <utility>

// Detect sharp edges, concave edges, and corner vertices on a surface mesh.
//
// Parameters:
//   mesh            - the input surface mesh (MPMesh / CGAL Polyhedron_3)
//   thres_concave   - dot-product threshold for classifying an edge as concave
//                     (c_B-c_A).normalized . n_A > threshold => concave
//   angle_sharp_deg - dihedral angle threshold in degrees; convex edges whose
//                     dihedral angle exceeds this value are classified as sharp
//   sharp_edges     - output: sorted (min,max) vertex-id pairs of sharp edges
//   concave_edges   - output: sorted (min,max) vertex-id pairs of concave edges
//   corners         - output: vertex ids of corner vertices (>= 3 sharp edges)
void find_feature_edges(
    const MPMesh& mesh,
    double thres_concave,
    double angle_sharp_deg,
    std::set<std::pair<int,int>>& sharp_edges,
    std::set<std::pair<int,int>>& concave_edges,
    std::set<int>& corners);
