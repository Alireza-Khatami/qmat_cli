---
  Your two options

  Option A — Extract just the algorithm (recommended)

  The actual math in find_feature_edges only needs:
  - Triangle normals (dot products, cross products)
  - Dihedral angle between adjacent face normals
  - No Geogram, no libigl, no MATFP internals needed

  You already have CGAL with halfedge access on your input mesh. The algorithm can be rewritten in ~100 lines using your
  existing types. Much cleaner.

  Option B — Pull in the full MATFP dependency

  1. Clone MATFP as a git submodule
  2. Install Geogram and Libigl via vcpkg
  3. Link all MATFP internal targets in your CMakeLists
  4. Adapt Args struct calls to your own config

  This brings in a large transitive dependency tree for one file.

  ---
  My recommendation: go with Option A. I can rewrite find_feature_edges using your CGAL mesh directly — no new
  dependencies. Want me to do that?

✻ Cooked for 50s



● Both files are written. Here's a summary of what was implemented:

  sharp_feature_detection.h — declares a single free function:
  void find_feature_edges(const MPMesh&, double thres_concave, double angle_sharp_deg,
      std::set<std::pair<int,int>>& sharp_edges,
      std::set<std::pair<int,int>>& concave_edges,
      std::set<int>& corners);

  sharp_feature_detection.cpp — rewritten using MPMesh halfedge structure directly, no Geogram/Libigl/MATFP dependencies:

  - Edge iteration: halfedges_begin()/end(), each edge processed once via &*h > &*h->opposite() guard
  - Border skip: open-boundary halfedges are skipped
  - Face normal: cross product of edge vectors from h->facet()->halfedge()
  - Face centroid: average of 3 vertex positions
  - Concave test: same dot-product criterion as original MATFP — (c_B−c_A).norm·n_A > thres_concave OR (c_A−c_B).norm·n_B
  > thres_concave
  - Sharp test: acos(n_A·n_B) > angle_sharp_rad (only runs if not concave)
  - Corners: vertices with >2 sharp-edge neighbors
  - Output edges identified by sorted (min_id, max_id) pairs using MPVertex::id
s



● Yes, find_feature_edges detects both — it outputs two separate sets:

  - sharp_edges — convex dihedral angle exceeds thres_convex (default 30°)
  - concave_edges — dot product test (c_B - c_A).normalized · n_A > thres_concave (default 0.18)