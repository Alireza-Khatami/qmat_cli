## Condition Checks — Current State

---

### Active

| Check | Function | Location |
|-------|----------|----------|
| `DifferentClusterType` | `CanMerge` | `src/SlabMesh.cpp:4344` |
| `TopoNotContractable` | `MinCostEdgeCollapse` | `src/SlabMesh.cpp:1376` |
| `InversionWouldOccur` | `MinCostEdgeCollapse` (main path) | `src/SlabMesh.cpp:1385` |
| `InversionWouldOccur` | `MinCostBoundaryEdgeCollapse` (boundary path) | `src/SlabMesh.cpp:1143` |

---

### Not Active

| Check | Function | Status |
|-------|----------|--------|
| `StaleEdge` | spike queue loop | deleted (spike phase removed) |
| `InvalidVertex` | spike queue loop | deleted (spike phase removed) |
| `SharpNotContractable` (junction vertex) | `CanMerge` | commented out (`src/SlabMesh.cpp:4301`) |
| `SharpNotContractable` (pre-marked sharp vertex) | `CanMerge` | commented out (`src/SlabMesh.cpp:4323`) |
| `WouldCreateNonManifold` | `CanMerge` | commented out (`src/SlabMesh.cpp:4351`) |
| `WouldCreateFoldOver` | `MinCostEdgeCollapse` | commented out (`src/SlabMesh.cpp:1520`) |
| `WouldExceedCurvatureThreshold` | `MinCostEdgeCollapse` | commented out (`src/SlabMesh.cpp:1526`) |
