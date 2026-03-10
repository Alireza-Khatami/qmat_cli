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
