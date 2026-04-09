## Condition Checks — Current State

---

### Active

| Check | Where | Location |
|-------|-------|----------|
| `DifferentClusterType` | inside `CanMerge` | `src/SlabMesh.cpp:4377` |
| `TopoNotContractable` | `MinCostEdgeCollapse` | `src/SlabMesh.cpp:1376` |
| `InversionWouldOccur` | `MinCostEdgeCollapse` | `src/SlabMesh.cpp:1385` |
| `InversionWouldOccur` | `MinCostBoundaryEdgeCollapse` | `src/SlabMesh.cpp:1143` |

---

### Not Active

| Check | Where | Status |
|-------|-------|--------|
| `StaleEdge` | spike queue loop | deleted (spike phase removed) |
| `InvalidVertex` | spike queue loop | deleted (spike phase removed) |
| `SharpNotContractable` (junction vertex) | inside `CanMerge` | commented out |
| `SharpNotContractable` (pre-marked sharp vertex) | inside `CanMerge` | commented out |
| `WouldCreateNonManifold` | inside `CanMerge` | commented out |
| `WouldExceedCurvatureThreshold` | `MinCostEdgeCollapse` | commented out |
