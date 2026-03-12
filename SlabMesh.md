 ⎿  Read 187 lines

∴ Thinking…

  Now I have a clear picture. Here's the plan:

  Vertex label enum (add to SlabMesh.h or PrimMesh.h):
  enum class VertexTopoLabel : uint8_t {
      REGULAR = 0,
      BOUNDARY = 1,
      NON_MANIFOLD_EDGE = 2,
      NON_MANIFOLD_CORNER = 3   // on both a boundary edge AND a non-manifold edge
  };

  Steps:
  1. Add VertexTopoLabel topo_label field to PrimVertex (in PrimMesh.h)

  2. Create a ComputeTopoLabels() method in SlabMesh that derives these labels from the existing fake_boundary_vertex
  and non_manifold_vertex flags, then invoke it from main_cli.cpp before running Simplify
  3. In SlabMesh::Simplify, validate that both vertices in an edge share the same topological label before attempting
  any collapse operation

  Wait, I'm noticing that DistinguishVertexType() already populates those vertex classification flags, so I can leverage
   that existing logic rather than recomputing from scratch.

  Since LoadSlabMesh() already calls DistinguishVertexType() by line 2723, the vertex flags are ready to use by the time
   we reach line 314 in main_cli.cpp. The real work is just converting those flags into a unified label and then
  filtering collapses based on label compatibility. The cleanest path forward is to add a VertexTopoLabel enum and a
  topo_label field to track vertex classification, then implement ComputeTopoLabels() to populate these labels from the
  existing flags—this keeps the logic encapsulated within SlabMesh and requires minimal changes elsewhere.

  The label mapping follows the same logic as the Python version: NM_CORNER combines fake boundary and non-manifold
  status, NM_EDGE is non-manifold without fake boundary, BOUNDARY is fake boundary without non-manifold, and REGULAR is
  neither. Once labels are computed, I'll add the compatibility check into Simplify's collapse loops...

  Adding label compatibility checks... I see the user wants me to add vertex labeling in main_cli.cpp before the
  Simplify call, mirroring the logic from the Python script's check_non_manifold function to categorize vertices as
  boundary, non-manifold edge, or non-manifold corner vertices. Then I need to modify the SlabMesh::Simplify function to
   only allow edge collapses when both vertices share the same label.

  To implement this in C++, I'll add a VertexTopoLabel enum and a topo_label field to PrimVertex, then create a
  LabelVertices() method in SlabMesh that assigns labels based on the existing vertex flags, call it before Simplify,
  and update the collapse condition to check label compatibility.

  Now I need to verify that DistinguishVertexType is being called in LoadSlabMesh to ensure the flags are properly
  initialized.



  ● All three edits are done. Here's a summary of what was added:

  ClusterBoundaryPoints() — now stores results into the shape-level fields:
  - boundary_clusters ← the per-cluster vertex lists
  - vertex_cluster_id ← per-vertex cluster id lookup (used by sidecar writer)

  ComputeInputNMM() — inserts just before Export():
  1. input.compute_normals() — required by face-adjacency BFS
  2. ClusterBoundaryPoints(25.0, 60.0) — classifies faces + tags is_steep_tetrahedron
  3. DetermineTopology() — sets the 4 topo flags

  Then immediately after Export(), writes <meshname>_mat_topo.txt with one line per active vertex:
  <idx> <steep> <sheet> <seam> <junction> <topo_bnd>  <bp_count> <bp0>...  <cl_count> <cl0>...

  LoadInputNMM() — reads the sidecar at the very end (after slab_mesh is fully built) and applies all fields — is_steep_tetrahedron, four topo flags,
  bplist, bplist_clusters — directly to each SlabVertex.




 ---
  1. InsertSavedPoint (line 928–937) — missing field copies

  It already copies bplist, sphere, slab_A/b/c, etc., but skips all new fields:
  // currently copies bplist ✓ but NOT these:
  // bplist_clusters, is_steep_tetrahedron, topo_is_sheet/seam/junction/boundary
  Fix: add copies for those fields alongside line 930.

  ---
  2. MinCostBoundaryEdgeCollapse (line 1070) — redundant but harmless

  vertices[vid_tgt].second->bplist = temp_bplist;  // line 1070
  MergeVertices now already unions bplist and bplist_clusters into vid_tgt. This line overwrites bplist with the same
  data — redundant but correct. bplist_clusters is already handled. No change needed.

  ---
  3. MinCostEdgeCollapse (line 1381) — bp reassignment missing cluster propagation

  When compute_hausdorff is on, each boundary point in temp_bplist gets re-homed to its new closest MAT vertex
  (min_index):
  vertices[min_index].second->bplist.insert(temp_ind);  // line 1381
  // ← bplist_clusters[temp_ind] is NOT propagated to min_index
  After MergeVertices, vid_tgt has the full merged bplist_clusters. When a bp moves to min_index, its cluster id needs
  to follow it. Fix: after line 1381, also write the cluster id into min_index->bplist_clusters.

  ---
  Summary table:

  ┌───────────────────────────────────────┬───────────────────────────────────────┬─────────────────────────┐
  │               Location                │                 Issue                 │         Action          │
  ├───────────────────────────────────────┼───────────────────────────────────────┼─────────────────────────┤
  │ InsertSavedPoint                      │ Missing copies of 6 new fields        │ Add copies              │
  ├───────────────────────────────────────┼───────────────────────────────────────┼─────────────────────────┤
  │ MinCostBoundaryEdgeCollapse line 1070 │ bplist overwrite now redundant        │ No change needed        │
  ├───────────────────────────────────────┼───────────────────────────────────────┼─────────────────────────┤
  │ MinCostEdgeCollapse line 1381         │ bp re-homing doesn't carry cluster id │ Add cluster propagation │
  └───────────────────────────────────────┴───────────────────────────────────────┴─────────────────────────┘

  Want me to implement these two fixes (InsertSavedPoint + MinCostEdgeCollapse line 1381)?


             File            │          Location          │                         Change                         │
  ├───────────────────────────┼────────────────────────────┼────────────────────────────────────────────────────────┤
  │ SlabMesh.h                │ SlabVertex field           │ bplist → nmn_bplist                                    │
  ├───────────────────────────┼────────────────────────────┼────────────────────────────────────────────────────────┤
  │ SlabMesh.cpp              │ MergeVertices              │ svt->bplist / sv1->bplist / sv2->bplist → nmn_bplist   │
  ├───────────────────────────┼────────────────────────────┼────────────────────────────────────────────────────────┤
  │ SlabMesh.cpp              │ InsertSavedPoint           │ Added nmn_bplist copy (original bplist copy left       │
  │                           │                            │ untouched)                                             │
  ├───────────────────────────┼────────────────────────────┼────────────────────────────────────────────────────────┤
  │ ThreeDimensionalShape.cpp │ LoadInputNMM sidecar       │ sv->bplist = bplist → sv->nmn_bplist = bplist   



  
  Our new nmn_bplist field — correctly renamed everywhere:
  - SlabMesh.h:43 — declaration ✓
  - SlabMesh.cpp:347-348 — MergeVertices union ✓
  - SlabMesh.cpp:931 — InsertSavedPoint copy ✓
  - ThreeDimensionalShape.cpp:495 — sidecar loader ✓

  Original PrimVertex::bplist — all untouched (correctly left as bplist):
  - SlabMesh.cpp:930, 1059, 1061, 1077, 1323, 1325, 1388, 1777, 1779, 1934, 1936 — all Hausdorff tracking code ✓
  - ThreeDimensionalShape.cpp:768 — ComputeHausdorffDistance ✓
  - GLWidget.cpp:2094-2095 — visualization, reads from input_nmm vertices ✓

  Unrelated bplist fields — untouched:
  - NonManifoldMesh_Vertex::bplist (nonmanifoldmesh.h:28, .cpp:870-871) ✓
  - pole_bplist on Mesh.h and nonmanifoldmesh.h ✓
  - EvaluateVertexDistanceErrorEnvelope parameter named bplist (just a local param name) ✓

  Everything is correct. No accidental touches to original code.