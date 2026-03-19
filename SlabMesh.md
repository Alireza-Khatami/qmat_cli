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




  ##################################################



  ● Done. The new CanMerge enforces exactly two conditions:

  1. Pure sheet check — both vertices must have topo_is_sheet=true and all three of topo_is_seam, topo_is_junction,
  topo_is_boundary must be false. A single non-sheet flag on either vertex blocks the collapse.
  2. Voronoi neighbor check — iterates every bp1 in v1->nmn_bplist, looks up its Voronoi neighbors, and checks if any
  bp2 from v2->nmn_bplist appears in that set. Returns true on the first hit, false if no pair is found.




  ● Now I have everything I need. Here's the plan:

  1. SlabMesh.h — add CollapseRecord struct + collapse_history + record_collapses
  2. SlabMesh.cpp — record each collapse in MinCostEdgeCollapse and MinCostBoundaryEdgeCollapse
  3. main_cli.cpp — add --visualize flag, snapshot initial MAT, then launch a Polyscope step-through viewer after
  simplification




  ┌───────┬─────────────────────┬──────────────────────────────────────────────────────┐
  │ Level │        Where        │                    What it checks                    │
  ├───────┼─────────────────────┼──────────────────────────────────────────────────────┤
  │ 1     │ Simplify loop       │ Same topo_label                                      │
  ├───────┼─────────────────────┼──────────────────────────────────────────────────────┤
  │ 2a    │ CanMerge            │ Both are pure sheet (no seam/junction/boundary flag) │
  ├───────┼─────────────────────┼──────────────────────────────────────────────────────┤
  │ 2b    │ CanMerge            │ Voronoi adjacency between bplists                    │
  ├───────┼─────────────────────┼──────────────────────────────────────────────────────┤
  │ 3     │ MinCostEdgeCollapse │ topo_contractable edge flag                          │
  ├───────┼─────────────────────┼──────────────────────────────────────────────────────┤
  │ 4     │ MinCostEdgeCollapse │ No geometric inversion (Contractible)                │
  └───────┴─────────────────────┴──────────────────────────────────────────────────────┘



● Done. Summary of changes:

  - ComputeBpClusters (new static helper) — union-find over a bp set using mesh edge connectivity; returns the connected
   components
  - ClusterNMNBplist — now just calls ComputeBpClusters per vertex (same result, no logic change)
  - MergeVertices — after unioning nmn_bplist, calls ComputeBpClusters on the merged set so only surface patches with a
  cross-edge get joined; disconnected patches stay separate
  - CanMerge — condition 2 now uses the halfedge circulator directly instead of voronoi_neighbors; allows collapse only
  if at least one bp from v1's region is mesh-adjacent to at least one bp from v2's region



    What it does:
  - Iterates over every bp in v1's full bplist
  - For each bp, walks all its mesh neighbours via halfedge circulator
  - If any neighbour is in v2's full bplist → returns true (collapse allowed)
  - If none found → returns false (collapse rejected)

  One thing to note: it currently checks the flat bplist, not the clusters explicitly. Checking the flat bplist is
  equivalent — since clusters are just subsets of the bplist, finding any cross-edge between the two bplists is the same
   as finding a cross-edge between any cluster of v1 and any cluster of v2. The result is identical.



#can merge conditions :


  Only two combinations are allowed:
  - Both vertices are T2 (each has exactly 2 disconnected surface patches)
  - Exactly one vertex is T1 (one connected patch) and the other is anything except T1

  Everything else (T2+T3, T3+T3, T1+T1, etc.) → reject

  ---
  Condition 2 — boundary exception

  If either vertex has topo_is_boundary = true:
  - The merge is only allowed if at least one of the two is T1
  - If neither is T1 → reject

  (So a boundary T2 cannot merge with another T2; it needs a T1 partner)

  ---
  Condition 3 — surface adjacency

  At least one boundary point (nmn_bplist) from v1 must share a mesh edge on the input surface with at least one boundary
  point from v2.

  Checked via CGAL halfedge circulator: walk all mesh neighbours of each bp in v1's bplist, if any neighbour is in v2's
  bplist → allow. If no cross-edge found → reject

  ---
  All conditions must pass for the collapse to be allowed.




● Yes — the same filtering applies to all of them. Currently:

  - Sheet: all incident edges have exactly 2 faces
  - Seam: at least 1 incident edge has more than 2 faces
  - Junction: more than 2 seam edges converge
  - Boundary: at least 1 incident edge has exactly 1 face

  With T1-face filtering, the rule becomes: before classifying any edge, filter out faces containing T1 vertices, then
  apply the same thresholds to the remaining face count:

  - An edge whose T1-filtered face count is 2 → normal 2-manifold edge (contributes to sheet)
  - An edge whose T1-filtered face count > 2 → seam edge
  - An edge whose T1-filtered face count = 1 → boundary edge
  - An edge whose T1-filtered face count = 0 → ignored entirely (pure spike edge, doesn't contribute to any
  classification)

  And the vertex-level flags follow from those filtered edge classifications exactly as before.

  Shall I implement this? I'll need to look at where DetermineTopology currently sets these flags to make the targeted
  change.