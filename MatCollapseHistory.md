
  MatCollapseHistory.h/.cpp

  ┌────────────────────────────────┬──────────────────────────────────────────────────────────────────────────────┐
  │              Part              │                                 What it does                                 │
  ├────────────────────────────────┼──────────────────────────────────────────────────────────────────────────────┤
  │ CollapseRecord                 │ One entry per collapse: step, vid_src, vid_tgt, bplist_after                 │
  ├────────────────────────────────┼──────────────────────────────────────────────────────────────────────────────┤
  │ MeshSnapshot                   │ Full lightweight mesh state at a keyframe step                               │
  ├────────────────────────────────┼──────────────────────────────────────────────────────────────────────────────┤
  │ Record(step, src, tgt, bplist) │ Called by Simplify at every collapse                                         │
  ├────────────────────────────────┼──────────────────────────────────────────────────────────────────────────────┤
  │ TakeKeyframe(step, sm)         │ Called every keyframe_interval collapses                                     │
  ├────────────────────────────────┼──────────────────────────────────────────────────────────────────────────────┤
  │ GetAncestors(vid)              │ Returns all original vertices whose chain leads to vid (DFS on merge forest) │
  ├────────────────────────────────┼──────────────────────────────────────────────────────────────────────────────┤
  │ GetLineage(vid)                │ Returns the sub-sequence of collapses relevant to vid, in order              │
  ├────────────────────────────────┼──────────────────────────────────────────────────────────────────────────────┤
  │ GetKeyframeBefore(step)        │ Nearest snapshot ≤ step for global scrubbing                                 │
  └────────────────────────────────┴──────────────────────────────────────────────────────────────────────────────┘

  Next steps when you're ready:
  1. Add a MatCollapseHistory history field to SlabMesh
  2. Call history.Record(...) + history.TakeKeyframe(...) inside Simplify
  3. Wire up the Polyscope UI (slider, prev/next, ancestry point cloud)




  
● All done. Here's what was implemented:

  MatCollapseHistory.h/.cpp
  - CollapseRecord now stores bplist_src1, bplist_src2 (before merge) + bplist_after (after merge)
  - Record() takes all three bplist snapshots and registers both sources in the merge forest

  SlabMesh.h
  - #include "MatCollapseHistory.h" + MatCollapseHistory history field added

  SlabMesh.cpp — both MinCostBoundaryEdgeCollapse and MinCostEdgeCollapse
  - Capture bplist_src1/src2 from v1/v2 before MergeVertices
  - Call history.Record(...) with bplist_after from vid_tgt after the merge
  - Auto-keyframe every keyframe_interval collapses

  main_cli.cpp — History panel in the ImGui UI

  ┌───────────────────────────┬────────────────────────────────────────────────────────────────────────────┐
  │          Control          │                                  Function                                  │
  ├───────────────────────────┼────────────────────────────────────────────────────────────────────────────┤
  │ Step slider               │ Scrubs through keyframe snapshots, shown as "History Snapshot" point cloud │
  ├───────────────────────────┼────────────────────────────────────────────────────────────────────────────┤
  │ Live button               │ Returns to live mesh view                                                  │
  ├───────────────────────────┼────────────────────────────────────────────────────────────────────────────┤
  │ Show ancestry cloud       │ Green point cloud of all vertices ever merged into selected vertex         │
  ├───────────────────────────┼────────────────────────────────────────────────────────────────────────────┤
  │ Merge step slider + < >   │ Steps through the lineage of selected vertex                               │
  ├───────────────────────────┼────────────────────────────────────────────────────────────────────────────┤
  │ Show bplists at this step │ Blue = src1 bplist, orange = src2 bplist, cyan = merged result    
