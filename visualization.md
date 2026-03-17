
  How it works:
  - Add --visualize to any --simplify run:
  qmat_cli model.obj --simplify 500 --visualize
  - A Polyscope window opens before simplification starts
  - Every accepted collapse fires on_collapse_cb which rebuilds the live MAT arrays and calls frameTick() — so you see the mesh shrinking in real time
  - The red edge / red·orange·green points show exactly which edge is being collapsed and where the result lands
  - Right panel controls:
    - Pause — halts between collapses (spins on frameTick() until unpaused)
    - Step — advances exactly one collapse when paused
    - Update every N slider — throttle display updates (e.g. every 10 collapses) for speed
  - Left panel layers: Initial MAT (grey, hidden by default), live MAT Faces/Edges, Collapsed Edge, v1/v2/result points — all toggleable
  - When simplification finishes, the final MAT is shown and you can inspect/rotate freely until you close the window

---

## Selecting a MAT vertex and viewing its nmn_bplist

### What is nmn_bplist?
Each MAT vertex was produced by a set of input mesh boundary points (samples from the
Delaunay triangulation).  `nmn_bplist` stores the indices of those boundary points.
Visualising it shows *which part of the input mesh surface* a given MAT vertex "belongs to".

### Data structures involved
- `SlabVertex::nmn_bplist`  — `std::set<unsigned>`, indices into `pmesh->pVertexList`
- `pmesh->pVertexList[i]->point()` — the 3-D position of boundary point i on the input mesh
- `ViewerState::idx_to_vid`  — `std::vector<unsigned>` that maps the Polyscope point-cloud
  index (0, 1, 2, …) back to the original slab-mesh vertex id.  Rebuilt every time
  `UpdateMatStructures` is called so it always reflects the current simplified state.
- `ViewerState::selected_vid` — the currently highlighted vertex id (-1 = none)

### How idx_to_vid is built
`BuildMatArrays(sm)` iterates `sm.vertices` in order.  For every active entry
(`.first == true`) it appends the vertex position to the output array **and** pushes
the raw slab-mesh index `i` onto `idx_to_vid`.  The two vectors therefore stay in
lock-step: `verts[k]` is the position of the vertex whose id is `idx_to_vid[k]`.
`UpdateMatStructures` copies this vector into `ViewerState` so the ImGui callback
can use it.

### Click → bplist routine (frame by frame)
1. Every frame Polyscope calls `polyscope::state::userCallback` (the ImGui panel).
2. Inside the callback, `polyscope::pick::haveSelection()` is checked.
3. If true, `polyscope::pick::getSelection()` returns `{Structure*, local_idx}` —
   the structure the user clicked and the index of the clicked point within it.
4. We check whether the clicked structure is the `"MAT Verts"` point cloud.
5. If yes, `local_idx` is used as an index into `vs.idx_to_vid` to get the original
   slab-mesh vertex id: `unsigned vid = vs.idx_to_vid[local_idx]`.
6. `BplistPositions(sm, vid)` is called — it iterates `sm.vertices[vid].second->nmn_bplist`
   and for each boundary-point index `bp` reads
   `sm.pmesh->pVertexList[bp]->point()` to get the 3-D position on the input mesh.
7. Those positions are registered as the `"BPList selected"` point cloud (cyan dots).
8. The panel displays the vertex id and `nmn_bplist` size.
9. "Clear selection" sets `selected_vid = -1` and hides the cyan cloud.

### Per-collapse bplist display (automatic, no click needed)
At every visualised collapse step the callback also shows the bplists of the two
vertices being merged:
- `raw_v1` and `raw_v2` are passed directly from `MinCostEdgeCollapse` /
  `MinCostBoundaryEdgeCollapse` through `on_collapse_cb` — no position scan needed.
- `BplistPositions(sm, raw_v1)` → registered as `"BPList v1"` (blue dots)
- `BplistPositions(sm, raw_v2)` → registered as `"BPList v2"` (orange dots)
- After the collapse the two sets are unified into the result vertex, which is
  why the two coloured clouds visually "merge" into one as simplification proceeds.




s is a flat array where indices are never reused — when a vertex gets collapsed it is just marked inactive (.first = false), but its slot stays.
   So the array might look like:

  index:  0      1      2      3      4      5
  .first: true   false  true   true   false  true

  Polyscope's point cloud needs a contiguous array — no gaps. So BuildMatArrays builds one by skipping the inactive slots:

  iteration:  i=0 (active) → verts[0], idx_to_vid[0] = 0
              i=1 (dead)   → skip
              i=2 (active) → verts[1], idx_to_vid[1] = 2
              i=3 (active) → verts[2], idx_to_vid[2] = 3
              i=4 (dead)   → skip
              i=5 (active) → verts[3], idx_to_vid[3] = 5

  Now Polyscope shows 4 points (indices 0–3). When the user clicks point 1 in the viewer, Polyscope reports local_idx = 1. Without idx_to_vid you'd look up
  sm.vertices[1] — which is dead. With it you do idx_to_vid[1] = 2 → sm.vertices[2] — the correct live vertex.

  In short: idx_to_vid is the bridge between Polyscope's compact 0-based indexing and the slab mesh's sparse array with holes.



  he RegisterInputMesh call is there specifically to answer your earlier question — you said the nmn_bplist points "seem to be way off the original mesh
  surface" and asked if pmesh vertices are actually the real mesh vertices.

  By showing the input mesh in the same Polyscope window alongside the cyan bplist dots, you can immediately see whether they land on the surface or float
  in space. It's a visual sanity check.

  If the points DO land on the surface → the code is correct, pmesh->pVertexList indices are right.
  If they DON'T → there's a bug in how the bp IDs map to vertex positions and we need to investigate further.

  Without the mesh visible there's no reference to judge against. 


  ves in normalized space (0..1 roughly), while pmesh->pVertexList positions are in the original world space.

  That's why:
  - The input mesh appears large
  - The bplist dots appear "way off" — they're in world space but plotted against a normalized MAT

  The fix: divide the bplist positions by bb_diagonal_length in BplistPositions, and do the same for the input mesh vertices in RegisterInputMesh:
