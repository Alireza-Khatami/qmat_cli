  // set in openmeshfile() / importMA()
  input_nmm.pmesh  = &input;   // both point to the same input mesh
  slab_mesh.pmesh  = &input;
  input_nmm.domain = input.domain;

  Key observations:

  1. input_nmm is built once from the DT and exported. The commented-out lines in LoadInputNMM() show that input_nmm was
   originally meant to be reloaded from file too, but that was abandoned — only slab_mesh gets populated from the .ma
  file now.
  2. slab_mesh is the working copy used for simplification (Simplify(), DistinguishVertexType(), etc.). input_nmm is
  essentially discarded after export.
  3. BoundaryPoints in input_nmm and the boundary points implicitly used by slab_mesh both come from input.pVertexList —
   same source, populated separately.

  Bottom line: slab_mesh IS input_nmm's geometry, just reloaded from the .ma file into a different data structure that
  supports edge-collapse simplification. After LoadInputNMM(), input_nmm is effectively unused.




 When two MAT vertices v1 and v2 are merged, their nmn_bplist sets get unioned. But instead of just taking the raw
  union, you want to merge only the clusters that are actually connected to each other across the two vertices — meaning
   a cluster from v1 and a cluster from v2 get joined into one only if there exists at least one mesh edge between a
  point in the v1-cluster and a point in the v2-cluster.

  Clusters that have no such cross-edge stay as separate clusters in the merged result.

  So for example:
  - v1 has clusters: {A, B}, {C}
  - v2 has clusters: {D}, {E, F}
  - If mesh edge A–D exists → {A, B} and {D} merge → {A, B, D}
  - If no edge connects {C} to {E, F} → they stay separate: {C} and {E, F}
  - Result: {A, B, D}, {C}, {E, F}

  And the CanMerge criterion you're building toward is likely: t