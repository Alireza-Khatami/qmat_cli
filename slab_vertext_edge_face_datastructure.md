---
  SlabVertex

  Inherits from PrimVertex (connectivity + QEM data) and SlabPrim (quadric matrices).

  Connectivity (from PrimVertex):
  - edges_ — set<unsigned>: IDs of all incident MAT edges
  - faces_ — set<unsigned>: IDs of all incident MAT faces

  Geometry:
  - sphere — the inscribed sphere (center + radius) representing this MAT vertex
  - normal, index, collaspe_weight, mean_square_error

  QEM quadric (from SlabPrim):
  - slab_A, slab_b, slab_c — accumulated quadric for collapse cost
  - add_A, add_b, add_c — additional quadric contribution

  Boundary point information:
  - nmn_bplist — set<unsigned>: all input surface mesh vertex IDs (boundary points) whose Voronoi cell this MAT vertex
  represents
  - nmn_bplist_clusters — vector<set<unsigned>>: the connected components of nmn_bplist restricted to surface mesh
  edges; each entry is one cluster
  - nmn_cluster_type — ClusterType enum (T0–T5, T1_spike, T1_non_spike): derived from nmn_bplist_clusters.size()

  Topology flags (set by DetermineTopology / RecomputeVertexTopology):
  - topo_type — TopoType enum: Unknown / Sheet / Boundary / Seam / Junction
  - topo_is_sheet, topo_is_boundary, topo_is_seam, topo_is_junction — raw boolean flags before priority resolution
  - nf — signed int: max face-count across all incident edges (1=boundary, 2=sheet, >2=seam, -1=Unknown)
  - is_spike — true if all bplist points form a mesh-edge clique (Delaunay tetrahedron spike)

  Other flags:
  - boundary_vertex, fake_boundary_vertex, saved_vertex, non_manifold_vertex
  - is_pole, is_non_manifold, is_disk, is_boundary

  ---
  SlabEdge

  Inherits from PrimEdge and SlabPrim. Has no extra fields of its own.

  Connectivity (from PrimEdge):
  - vertices_ — pair<unsigned, unsigned>: the two endpoint vertex IDs
  - faces_ — set<unsigned>: IDs of all incident MAT faces

  Geometry:
  - sphere, cone — the interpolated sphere/cone geometry along the edge
  - collapse_cost, qem_error — cost used by the QEM priority queue

  Flags:
  - boundary_edge, fake_boundary_edge, non_manifold_edge
  - topo_contractable — pre-computed flag; false = this edge is marked unsafe to collapse (set before the main phase,
  checked in MinCostEdgeCollapse alongside CanMerge)

  ---
  SlabFace

  Inherits from PrimFace and SlabPrim. Has no extra fields of its own.

  Connectivity (from PrimFace):
  - vertices_ — set<unsigned>: the 3 vertex IDs of this triangle
  - edges_ — set<unsigned>: the 3 edge IDs of this triangle

  Geometry:
  - centroid — sphere at the face centroid
  - normal — face normal vector
  - st[2], valid_st — simple triangle geometry for inversion checks

  ---
  Storage in SlabMesh

  All three are stored as tombstone vectors — a pair<bool, T*> per element:

  vertices[i].first   // true = active,  false = deleted (pointer is dangling)
  vertices[i].second  // pointer to SlabVertex (only safe to dereference if .first == true)

  Same pattern for edges and faces. Deleted elements are never erased from the vector — their index is permanently their
   ID.



---
  How Edge Collapse Changes the Geometry

  When MergeVertices(vid_src1, vid_src2, vid_tgt) is called it works in five stages:

  ---
  Stage 1 — Allocate the new vertex

  InsertVertex appends a fresh SlabVertex to the end of the vertices vector and assigns it a new ID (vid_tgt). Its
  sphere, quadric, bplist, clusters, and topology are then filled in from the two sources (the merge logic described in
  previous conversations).

  ---
  Stage 2 — Plan which faces and edges survive (pre-deletion)

  Before deleting anything, the code pre-builds two lists by scanning the old connectivity:

  Faces to recreate (tri_vec):
  For every face incident to vid_src1 or vid_src2, if that face does not already contain vid_tgt (i.e. it is not one of
  the two triangles that collapse), the face is kept but with the source vertex ID replaced by vid_tgt:

  old face {vid_src1, A, B}  →  new face {vid_tgt, A, B}
  old face {vid_src2, C, D}  →  new face {vid_tgt, C, D}

  Faces that do contain vid_tgt at this point are the two triangles sharing the collapsing edge — they get implicitly
  discarded (not added to tri_vec).

  Edges to recreate (edge_vec):
  Same logic for edges: every edge of vid_src1 or vid_src2 that does not already touch vid_tgt is kept, with its
  endpoint remapped:

  old edge (vid_src1, X)  →  new edge (vid_tgt, X)
  old edge (vid_src2, Y)  →  new edge (vid_tgt, Y)

  The shared edge (vid_src1, vid_src2) itself touches both sources so it has vid_tgt as one endpoint after remapping and
   gets excluded.

  ---
  Stage 3 — Delete both source vertices

  DeleteVertex(vid_src1) and DeleteVertex(vid_src2) are called. Each one:

  1. Collects all incident edges into a local set, then calls DeleteEdge on each
  2. DeleteEdge removes the edge from both endpoint vertices' edges_ sets, then calls DeleteFace on all incident faces
  3. DeleteFace removes the face from all three vertices' faces_ sets and from all three edges' faces_ sets, then
  deletes the SlabFace object and sets .first = false
  4. After all edges (and their faces) are gone, deletes the SlabVertex object and sets .first = false

  So after this stage everything connected to vid_src1 and vid_src2 — all their edges and all their triangles — is
  completely gone from the mesh.

  ---
  Stage 4 — Recreate the surviving faces

  InsertFace(vset) is called for each entry in tri_vec. For each face it:

  1. Checks if the face already exists (deduplication) — skips if so
  2. Calls InsertEdge for each of the 3 pairs of vertices — which either finds the existing edge (if two neighbours were
   already connected) or creates a new one
  3. Creates a new SlabFace, registers it in each vertex's faces_ and each edge's faces_
  4. Recomputes the face's centroid, normal, and SimpleTriangle geometry

  ---
  Stage 5 — Recreate the surviving edges

  InsertEdge(va, vb) is called for each entry in edge_vec. Since InsertFace already called InsertEdge for all face
  edges, most of these are no-ops (the Edge(vid0,vid1,eid) check at the top returns early if the edge already exists).
  Only bare edges (edges not part of any face — e.g. boundary wire edges) actually get inserted here. When a new edge is
   created it also calls ComputeEdgeCone to set its slab geometry.

  ---
  What does NOT get updated automatically

  - Quadric matrices (slab_A, slab_b, slab_c) on the new edges and faces — these are recomputed explicitly in
  MinCostEdgeCollapse via EvaluateEdgeCollapseCost after the merge
  - Topology flags (topo_type, nf) — RecomputeVertexTopology(vid_tgt) is called explicitly right after the merge
  - Collapse costs of neighbouring edges — the priority queue is updated after the merge by re-evaluating all edges
  incident to vid_tgt