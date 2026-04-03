# MAT Simplification: Fold-Over Bug Analysis

## Problem Description

During MAT (Medial Axis Transform) simplification of the frame-shaped surface (the space between two squares), sequential edge collapses of sheet vertices produce geometrically crossing edges after several steps, even though each individual collapse appears locally valid.

### Visual Description

**Before (left state):**
- Sheet vertices (MS_Sheet, red circles) form a loop around the outer square boundary
- Junction vertices (MS_Junction, pink circles) sit at the corners of the inner square
- Edges (blue lines) connect adjacent sheet vertices, forming a correct non-self-intersecting octagonal ring

**After several collapses (right state):**
- Only 4 sheet vertices remain, roughly at the cardinal points (top, left, right, bottom) of the outer square
- The edges now geometrically cross each other — the connectivity forms a "twisted diamond" that self-intersects
- The junction vertices are still at the 4 inner square corners
- The edges connecting sheet vertices to junction vertices cross in the middle

---

## What Is Happening

The collapses are always between **adjacent vertices** (vertices connected by an edge). There are no non-adjacent collapses. The self-intersection develops gradually over many sequential steps.

**The accumulation mechanism:**
- Each collapse merges two adjacent sheet vertices into one
- The surviving vertex inherits all edges from both source vertices
- After several collapses, a single vertex has accumulated connections spanning a wider neighborhood
- Eventually, the surviving vertex has edges to junction vertices on geometrically opposite sides of the shape
- The straight-line connections between the sheet vertex positions and junction vertex positions now cross each other

Each individual collapse passes all local validity checks. The global geometric invalidity only emerges after multiple steps — this is a **cumulative effect** that no single-step check detects.

---

## Name of the Phenomenon

This is called a **"fold-over"** (also: "foldover" or "mesh fold-over during edge collapse") — a well-known problem in geometry processing.

**Related formal names:**

| Name | Description |
|------|-------------|
| **Fold-over prevention** | Standard term in mesh simplification literature |
| **Flip prevention during edge collapse** | Alternate common name |
| **Topological disk condition** | Formal requirement that the one-ring of a collapsed edge forms a valid disk; the geometric version is violated here even when the topological version passes |
| **Link condition** (Dey et al. 1999) | The graph-theoretic version; `WouldCreateNonManifold` implements this but only catches topological violations, not geometric ones |
| **Geometric validity / inversion-free collapse** | Broader requirement from Garland & Heckbert's QEM paper (1997); flagged as a concern separate from topology |
| **Orientation consistency** | Requirement that all face normals remain consistent after collapse — what `Contractible` tries to check |
| **Curve fold-over / polyline self-intersection** | The 1D boundary loop equivalent; the specific phenomenon happening here |

**Closest classical reference:** Hoppe et al. (1993) *"Mesh Optimization"* introduced the **"no fold-over condition"**: for every edge collapse, the new vertex position must not cause any triangle in the one-ring to reverse orientation, or any edge to cross another edge. For the 1D boundary loop, the equivalent condition is that the loop must remain non-self-intersecting after collapse.

---

## Root Cause in the Code

### 1. `prevent_inversion = false` (main_cli.cpp:1632)

```cpp
shape.slab_mesh.prevent_inversion = false;
```

This flag disables the entire inversion check block inside `MinCostEdgeCollapse` and `MinCostBoundaryEdgeCollapse`. The `if (prevent_inversion == true)` guard (SlabMesh.cpp:1378) never executes. As a result, **no geometric inversion check runs at all** during the sheet/boundary/seam collapse phases in the CLI path.

### 2. `Contractible` checks the wrong thing for this geometry

Even if `prevent_inversion` were `true`, `Contractible` would not catch this problem.

**What `Contractible` actually does** (SlabMesh.cpp:1051-1118):

For every MAT triangle face touching `vid_src1` (skipping faces that also contain `vid_src2`):
1. Compute the triangle normal **before** collapse using the 3 sphere centers → `pnorm`
2. Compute the triangle normal **after** collapse by replacing `vid_src1`'s center with the target position → `anorm`
3. Return `false` if `pnorm.Dot(anorm) < 0` — i.e. the triangle **flipped its normal direction**

**Why this is blind to the fold-over:**

The MAT of a frame-shaped surface (the region between two squares) is essentially **1-dimensional — a loop of edges with no triangular faces**. The sheet vertices form a boundary ring with no MAT triangles between them. `Contractible` only examines faces (triangles). If there are no faces between the sheet vertices, there is nothing to check — the function returns `true` immediately.

`Contractible` was designed for a **surface MAT** (2-manifold mesh with triangular faces). For a **boundary loop of sheet vertices** (1-manifold / curve), it is completely blind to the problem. A fold-over on a curve is a **2D edge crossing**, not a 3D triangle normal flip, and `Contractible` has no awareness of edge-edge intersections.

### 3. `WouldCreateNonManifold` checks topology, not geometry

`WouldCreateNonManifold` (SlabMesh.cpp:3880-4005) implements the **link condition** — it checks graph-theoretic manifold properties:
- Test 1/2: Both side-edges of a triangle are boundary edges
- Test 3: Two incident faces share the same third vertex
- Test 4: Both endpoints are boundary vertices but the shared edge is not
- Test 5: One-ring intersection (shared neighbour beyond safe third verts)

None of these tests check the geometric positions of the vertices. An edge collapse that creates a geometrically crossing edge passes all five tests if the graph connectivity remains a valid manifold.

### 4. `CanMerge` only checks cluster type

With the current configuration, `CanMerge` (SlabMesh.cpp:4004) only has one active condition:

```cpp
// Condition 2: same cluster type — ONLY ACTIVE CONDITION
if (v1->nmn_cluster_type != v2->nmn_cluster_type)
{ if (out_reason) *out_reason = RejectionReason::DifferentClusterType; return false; }
```

All other conditions (topo type, bplist neighbours, link condition) are commented out. Two sheet vertices that are adjacent will always pass this check.

---

## Summary of What Each Check Does and Does Not Catch

| Check | What it catches | Does it catch fold-over? |
|-------|----------------|--------------------------|
| `CanMerge` | Different cluster/topo type | No |
| `WouldCreateNonManifold` | Graph topology violations (link condition) | No — purely topological |
| `Contractible` | Triangle normal flips (only on faces) | No — blind to 1D loop geometry |
| `prevent_inversion` flag | Gates the `Contractible` call | Irrelevant here since `Contractible` is wrong check anyway |

---

## What a Correct Fix Would Look Like

For a **1D boundary loop** (sheet vertex ring), the correct fold-over check is:

> After collapsing edge (A, B) into vertex C at position P, check that **no pair of edges in the updated loop crosses geometrically**.

More specifically, the local version: for each edge (C, X) that results from the collapse, verify that it does not intersect any other edge (Y, Z) in the one-ring of C.

This requires a **2D or 3D edge-edge intersection test** on the sphere centers, not a face normal flip test. The projected positions of the sphere centers onto the original surface (or in 3D space) must form a non-self-intersecting curve after every collapse.

This is analogous to the **"no fold-over condition"** from Hoppe et al. (1993) applied to a polyline (1-manifold) rather than a triangle mesh (2-manifold).

---

## Files and Line Numbers for Reference

| Location | Line | Relevance |
|----------|------|-----------|
| `main_cli.cpp` | 1632 | `prevent_inversion = false` — disables all inversion checks |
| `SlabMesh.cpp` | 294–337 | `Contractible(vid_src, vid_tgt)` — first version, normal flip check |
| `SlabMesh.cpp` | 1051–1118 | `Contractible(vid_src1, vid_src2, v_tgt)` — second version with target position |
| `SlabMesh.cpp` | 1378 | `if (prevent_inversion == true)` gate in `MinCostEdgeCollapse` |
| `SlabMesh.cpp` | 1138 | `if (prevent_inversion == true)` gate in `MinCostBoundaryEdgeCollapse` |
| `SlabMesh.cpp` | 3880–4005 | `WouldCreateNonManifold` — topology only, no geometry |
| `SlabMesh.cpp` | 4004–4044 | `CanMerge` — only active condition is same cluster type |
| `SlabMesh.h` | 149 | `bool prevent_inversion` member declaration |
| `SlabMesh.h` | 311 | `WouldCreateNonManifold` declaration |
