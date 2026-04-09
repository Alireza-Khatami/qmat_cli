# MAT Edge Collapse Rejection Checks

All checks are listed in the **exact order they execute** in the code.
There are three collapse paths: **spike**, **boundary**, and **main (general)**.
Each path runs its own priority queue and calls a dedicated collapse function.

---

## Simplification Phases

```
Phase 0 — Spike collapse     (spike_collapse_queue  → MinCostEdgeCollapse, CollapseContext::Spike)
Phase 1 — Main collapse      (topo_collapse_queue   → MinCostEdgeCollapse, CollapseContext::Main)
          Boundary sub-pass  (boundary_edge_collapses_queue → MinCostBoundaryEdgeCollapse)
```

---

## Phase 0 — Spike Collapse Queue Loop

These two checks happen **in the queue loop itself**, before `MinCostEdgeCollapse` is called.

### 1. `StaleEdge`
- **Where:** `src/SlabMesh.cpp:2289–2290` — spike queue loop (`Simplify`)
- **Condition:** `!edges[eid].first` — the edge was already deleted (collapsed away by a prior operation and invalidated in the queue)
- **Primitives captured:** none (edge no longer exists)
- **Action:** `continue` — skip this queue entry

### 2. `InvalidVertex`
- **Where:** `src/SlabMesh.cpp:2293–2300` — spike queue loop (`Simplify`)
- **Condition:** `!ValidVertex(v1) || !ValidVertex(v2)` — one or both endpoints have been deleted
- **Primitives captured:** whichever of `v1`, `v2` still exist in the vertex table
- **Action:** `continue`

---

## `MinCostEdgeCollapse` — Spike Context (`CollapseContext::Spike`)
**Function:** `src/SlabMesh.cpp:1355`

Only the non-manifold check is run for spikes; cluster/topo guards are bypassed intentionally.

### 3. `NonManifold_*` (spike path — via `WouldCreateNonManifold` directly)
- **Where:** `src/SlabMesh.cpp:1373–1376`
- See [WouldCreateNonManifold sub-checks](#wouldcreatenonmanifold-sub-checks) below.

---

## `MinCostEdgeCollapse` — Main / Boundary Context (`CollapseContext::Main`)
**Function:** `src/SlabMesh.cpp:1355`

### 4. `CanMerge` gate (runs checks 4a–4e in order)
- **Called at:** `src/SlabMesh.cpp:1370–1371` (main path) / `src/SlabMesh.cpp:1134–1135` (boundary path)
- **CanMerge defined at:** `src/SlabMesh.cpp:4309`

#### 4a. `SharpNotContractable` — Junction vertex
- **Where:** `src/SlabMesh.cpp:4319–4334` (`CanMerge`, condition 0, first guard)
- **Condition:** either endpoint has cluster type `MS_Junction` or `MS_Junction_Boundary`
- **Rationale:** junction vertices are branch points of the seam/boundary skeleton; collapsing them would destroy the 1-skeleton topology
- **Primitives captured:** the junction endpoint(s)

#### 4b. `SharpNotContractable` — Pre-marked sharp vertex
- **Where:** `src/SlabMesh.cpp:4340–4354` (`CanMerge`, condition 0, second guard)
- **Condition:** both endpoints share a boundary-class type (`MS_Boundary`, `MS_Seam`, `MS_Sheet_Boundary`, `MS_Seam_Boundary`) **and** at least one has `sharpNotContractable = true` (set by `MarkSharpFeatureVertices` in the pre-pass)
- **Rationale:** sharp corners on feature curves must be preserved
- **Primitives captured:** the endpoint(s) with `sharpNotContractable = true`

#### 4c. `DifferentClusterType`
- **Where:** `src/SlabMesh.cpp:4358–4363` (`CanMerge`, condition 2)
- **Condition:** `v1->nmn_cluster_type != v2->nmn_cluster_type`
- **Rationale:** collapsing vertices of different semantic types (e.g. a sheet vertex into a seam vertex) would corrupt the MAT type labelling
- **Primitives captured:** both endpoints

#### 4d–4g. `WouldCreateNonManifold` (four sub-checks)
- **Called at:** `src/SlabMesh.cpp:4366`
- Delegated from `CanMerge` to `WouldCreateNonManifold`. See below.

---

## `WouldCreateNonManifold` Sub-checks
**Function:** `src/SlabMesh.cpp:3824`

Called from `CanMerge` (main path, `src/SlabMesh.cpp:4366`) or directly (spike path, `src/SlabMesh.cpp:1375`). Translated from the PMP `is_collapse_ok()` link-condition test.

### 5. `NonManifold_BoundaryEdgePair`
- **Where:** `src/SlabMesh.cpp:3890–3903` (PMP test 1 & 2)
- **Condition:** for any third vertex `ft` of a face incident to edge `(v0,v1)`: **both** the opposite edges `(v1,ft)` and `(ft,v0)` are boundary edges (each incident to exactly one face)
- **Rationale:** collapsing would merge two boundary chains at a single point, creating a non-manifold pinch
- **Primitives captured:** `vertices={ft}`, `edges={(v1,ft),(ft,v0)}`, `faces={(v0,v1,ft)}`

### 6. `NonManifold_SharedThirdVert`
- **Where:** `src/SlabMesh.cpp:3910–3920` (PMP test 3)
- **Condition:** the collapsed edge has exactly 2 incident faces **and** both share the same third vertex (manifold case: `vl == vr`)
- **Rationale:** collapse would reduce a closed triangular "cap" to a degenerate edge — a topological sphere pinch
- **Primitives captured:** `vertices={ft}` (the shared third vertex), `faces={(v0,v1,ft),(v0,v1,ft)}`

### 7. `NonManifold_BoundaryVertEdge`
- **Where:** `src/SlabMesh.cpp:3922–3931` (PMP test 4)
- **Condition:** both `v0` and `v1` are boundary vertices (each incident to a boundary edge) **but** the shared edge `(v0,v1)` itself is not a boundary edge
- **Rationale:** collapse would connect two separate boundary chains through an interior path, creating a non-manifold vertex
- **Primitives captured:** `vertices={v0,v1}`, `edges={(v0,v1)}`

### 8. `NonManifold_LinkCondition`
- **Where:** `src/SlabMesh.cpp:3962–3979` (PMP test 5)
- **Condition:** the one-rings of `v0` and `v1` share a common neighbour `nbr` that is **not** a face-third-vertex of the shared edge (not an exempt vertex)
- **Rationale:** the link condition — if a non-exempt shared neighbour exists, collapse would create a duplicate edge `(result,nbr)`, violating the manifold property
- **Primitives captured:** `vertices={nbr}`, `edges={(v0,nbr),(v1,nbr)}`

> **Note:** if the shared edge has 0 incident faces (e.g. a seam edge from a `.ma` file with no MAT triangles), the link condition is trivially satisfied and this check is skipped (`src/SlabMesh.cpp:3942`).

---

## Back in `MinCostEdgeCollapse` — After `CanMerge`

### 9. `TopoNotContractable`
- **Where:** `src/SlabMesh.cpp:1386–1392`
- **Condition:** `!edges[eid].second->topo_contractable` — a pre-computed flag; set to `false` on edges that are statically known to be unsafe before the simplification loop begins
- **How the flag is set — two cases:**
  1. **Degree-1 isolated endpoint** (`src/SlabMesh.cpp:2203–2205` and `3494–3496`): if the other endpoint of the edge has exactly 1 edge and 0 incident faces, it is a dangling vertex attached by only this edge with no triangles. Collapsing would orphan it or destroy the last connection.
  2. **Triangular boundary hole** (`src/SlabMesh.cpp:3471–3479`): if three edges form a closed triangular loop where all three are boundary edges (≤ 1 incident face each) and the loop encloses a hole, all three edges are marked non-contractable. Collapsing any one of them would collapse the hole boundary and corrupt the topology.
- **Note:** unlike `WouldCreateNonManifold` (checked at collapse time), this flag is computed statically during topology initialisation — it is not re-evaluated per collapse attempt.
- **Primitives captured:** `vertices={v1,v2}`, `edges={(v1,v2)}`

### 10. `InversionWouldOccur`
- **Where:** `src/SlabMesh.cpp:1445–1449` (main path) / `src/SlabMesh.cpp:1147–1150` (boundary path)
- **Condition:** `Contractible(v1, v2, sphere.center)` returns false **and** none of the three candidate positions (at `v1`, at `v2`, at midpoint) are contractible either — i.e. every candidate collapse position would invert adjacent slab faces
- **Rationale:** the QEM-optimal target position sits inside an inverted region; no safe fallback exists
- **Primitives captured:** `vertices={v1,v2}`, `edges={(v1,v2)}`
- **Action:** re-queues the edge with a large penalty cost (`+1e9`), does **not** return false immediately — the edge re-enters the queue
- **Note (boundary path):** simpler check at `src/SlabMesh.cpp:1145` — no candidate fallbacks, immediate rejection

### 11. `WouldCreateFoldOver`
- **Where:** `src/SlabMesh.cpp:1521–1525` (main path only — not run on boundary path)
- **Function defined at:** `src/SlabMesh.cpp:4079`
- **Condition:** `WouldCreateFoldOver(v1, v2, sphere.center)` — any of the new edges produced by the collapse (from the merged vertex to each neighbour) would geometrically cross an existing non-adjacent edge
- **Rationale:** catches boundary-loop self-intersections that topology checks miss (they are blind to 1-D loop geometry)
- **Detection:** minimum-distance test between pairs of 3-D segments; flags crossing only when both closest points are strictly interior (not at shared endpoints)
- **Primitives captured:** `vertices={X,Y,Z}` where `X` is the new neighbour whose edge would cross existing edge `(Y,Z)`, `edges={(Y,Z)}`

### 12. `WouldExceedCurvatureThreshold`
- **Where:** `src/SlabMesh.cpp:1535–1539` (main path only — not run on boundary path)
- **Function defined at:** `src/SlabMesh.cpp:4160`
- **Applies to:** endpoints whose cluster type is `MS_Boundary`, `MS_Seam`, `MS_Seam_Boundary`, or `MS_Sheet_Boundary`
- **Condition:** for either endpoint `v`, the turning angle at `v` between its far same-type chain neighbour and the other collapsing vertex exceeds `feature_angle_threshold`
  - Also rejects if an endpoint has **≥ 2** same-type neighbours (junction-like)
  - Skips the angle check if an endpoint is a chain end (0 same-type far neighbours)
- **Rationale:** prevents collapsing through a sharp bend in a seam/boundary curve, which would smooth away a geometric feature
- **Primitives captured:** `vertices={v, far}` (the bent endpoint and its far chain neighbour), `edges={(v,far),(v,partner)}`

---

## `MinCostBoundaryEdgeCollapse` — Boundary Path
**Function:** `src/SlabMesh.cpp:1124`

Used for edges in the boundary sub-queue. Runs a subset of the main checks.

| Step | Check | Location | Notes |
|------|-------|----------|-------|
| 1 | `CanMerge` (4a–4g) | `src/SlabMesh.cpp:1134` | identical to main path |
| 2 | `InversionWouldOccur` | `src/SlabMesh.cpp:1145–1150` | simpler: just `Contractible(v1,v2,sphere.center)` — no candidate fallbacks, immediate rejection |

> `WouldCreateFoldOver` and `WouldExceedCurvatureThreshold` are **not** run on the boundary path.

---

## Pre-pass: `MarkSharpFeatureVertices`
**Function:** `src/SlabMesh.cpp:4236`

Runs **once before** `Simplify()`. Sets `sharpNotContractable = true` on vertices that check 4b would always reject, so they are cheaply filtered without re-running the angle computation every collapse.

| Type | Rule |
|------|------|
| `MS_Junction`, `MS_Junction_Boundary` | Always marked (branch points by definition) |
| `MS_Boundary`, `MS_Seam`, `MS_Seam_Boundary`, `MS_Sheet_Boundary` | Marked if **≥ 2** same-type neighbours **and** the minimum turning angle across all neighbour pairs exceeds `feature_angle_threshold` |

---

## Rejection Color Key (Polyscope / `ExportSkeletonPLY`)

| Reason | Color |
|--------|-------|
| `StaleEdge` | Light grey `(200,200,200)` |
| `InvalidVertex` | Dark grey `(120,120,120)` |
| `DifferentTopoType` | Mid grey `(180,180,180)` |
| `DifferentClusterType` | **ORANGE** `(255,140,0)` |
| `BplistNotNeighbors` | Mid grey `(130,130,130)` |
| `NoPmesh` | Dark grey `(100,100,100)` |
| `InversionWouldOccur` | Dull yellow `(180,180,0)` |
| `TopoNotContractable` | **RED** `(255,0,0)` |
| `NonManifold_BoundaryEdgePair` | **GOLD** `(255,215,0)` |
| `NonManifold_SharedThirdVert` | **YELLOW** `(255,255,0)` |
| `NonManifold_BoundaryVertEdge` | **GREEN** `(0,255,0)` |
| `NonManifold_LinkCondition` | **CYAN** `(0,255,255)` |
| `WouldCreateFoldOver` | **BLUE** `(0,0,255)` |
| `SharpNotContractable` | **VIOLET** `(148,0,211)` |
| `WouldExceedCurvatureThreshold` | **MAGENTA** `(255,0,255)` |
| Never attempted | White `(255,255,255)` |





 ┌──────────────────────────────────────┬────────────────────────────────────┐
  │                Check                 │              Location              │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ StaleEdge                            │ SlabMesh.cpp:2289–2290             │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ InvalidVertex                        │ SlabMesh.cpp:2293–2300             │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ Spike WouldCreateNonManifold         │ SlabMesh.cpp:1373–1376             │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ CanMerge call (main)                 │ SlabMesh.cpp:1370–1371             │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ CanMerge call (boundary)             │ SlabMesh.cpp:1134–1135             │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ CanMerge definition                  │ SlabMesh.cpp:4309                  │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ 4a SharpNotContractable (junction)   │ SlabMesh.cpp:4319–4334             │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ 4b SharpNotContractable (pre-marked) │ SlabMesh.cpp:4340–4354             │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ 4c DifferentClusterType              │ SlabMesh.cpp:4358–4363             │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ WouldCreateNonManifold definition    │ SlabMesh.cpp:3824                  │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ NonManifold_BoundaryEdgePair         │ SlabMesh.cpp:3890–3903             │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ NonManifold_SharedThirdVert          │ SlabMesh.cpp:3910–3920             │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ NonManifold_BoundaryVertEdge         │ SlabMesh.cpp:3922–3931             │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ NonManifold_LinkCondition            │ SlabMesh.cpp:3962–3979             │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ TopoNotContractable                  │ SlabMesh.cpp:1386–1392             │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ InversionWouldOccur (main)           │ SlabMesh.cpp:1445–1449             │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ InversionWouldOccur (boundary)       │ SlabMesh.cpp:1147–1150             │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ WouldCreateFoldOver                  │ SlabMesh.cpp:1521–1525 (def: 4079) │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ WouldExceedCurvatureThreshold        │ SlabMesh.cpp:1535–1539 (def: 4160) │
  ├──────────────────────────────────────┼────────────────────────────────────┤
  │ MarkSharpFeatureVertices             │ SlabMesh.cpp:4236                  │
  └──────────────────────────────────────┴────────────────────────────────-----




   Crash / Undefined Behavior (do NOT comment out)

  ┌─────────────────────────────────┬───────────────────────────────────────────────────────────────────────────────┐
  │              Check              │                              Why it would crash                               │
  ├─────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────┤
  │ StaleEdge (SlabMesh.cpp:2289)   │ Proceeds to dereference edges[eid].second on a deleted edge — UB/crash        │
  ├─────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────┤
  │ InvalidVertex                   │ Proceeds to access vertices[v1/v2].second on deleted vertices — UB/crash      │
  │ (SlabMesh.cpp:2293)             │                                                                               │
  ├─────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────┤
  │ WouldCreateNonManifold          │ MergeVertices internally assumes valid manifold adjacency; corrupted topology │
  │ (SlabMesh.cpp:3824)             │  can crash in its edge/face update loops                                      │
  ├─────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────┤
  │ TopoNotContractable             │ Flag was set because the collapse is structurally invalid at the mesh level;  │
  │ (SlabMesh.cpp:1386)             │ MergeVertices can corrupt the half-edge/adjacency structure                   │
  └─────────────────────────────────┴───────────────────────────────────────────────────────────────────────────────┘

  ---
  Wrong Results, No Crash

  ┌───────────────────────────────┬──────────────────────────────────────────────────────────────────────────────────┐
  │             Check             │                                 What goes wrong                                  │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────┤
  │ SharpNotContractable          │ Sharp feature corners and junction branch points get smoothed away — incorrect   │
  │                               │ MAT topology                                                                     │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────┤
  │ DifferentClusterType          │ Vertices of different semantic types get merged — cluster type labels become     │
  │                               │ meaningless                                                                      │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────┤
  │ InversionWouldOccur           │ Slab faces get inverted — incorrect geometry, broken Contractible assumptions    │
  │                               │ downstream                                                                       │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────┤
  │ WouldCreateFoldOver           │ Boundary loops self-intersect — incorrect 1-skeleton                             │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────┤
  │ WouldExceedCurvatureThreshold │ Sharp bends in seam/boundary chains collapse away — loss of geometric features   │
  └───────────────────────────────┴──────────────────────────────────────────────────────────────────────────────────┘

  ---
  Short answer: commenting out StaleEdge, InvalidVertex, WouldCreateNonManifold, or TopoNotContractable will likely
  crash or corrupt the mesh data structure. The rest will give wrong results but won't break the program's control flow.

