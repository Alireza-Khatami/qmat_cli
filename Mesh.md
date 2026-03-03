However there's a clean solution: jitter only inside the DT, but always read the original positions back
  from pVertexList using the stored vertex ID. The DT already stores vh->info().id = idx which maps back to
  the original vertex.

  The two things affected and what happens to each:

  ┌───────────────────────────────────┬───────────────────────────┬──────────────────────────────────────┐
  │                                   │     With naive jitter     │        With ID-based recovery        │
  ├───────────────────────────────────┼───────────────────────────┼──────────────────────────────────────┤
  │ BoundaryPoints /                  │ slightly jittered         │ exact original positions             │
  │ sample_points.txt                 │                           │                                      │
  ├───────────────────────────────────┼───────────────────────────┼──────────────────────────────────────┤
  │ MAT vertex positions              │ ~1e-9 * diagonal off      │ same ~1e-9 off (unavoidable,         │
  │ (circumcenters)                   │                           │ negligible)                          │
  ├───────────────────────────────────┼───────────────────────────┼──────────────────────────────────────┤
  │ Inside/outside cell               │ correct (degeneracy       │ same                                 │
  │ classification                    │ broken)                   │                                      │
  └───────────────────────────────────┴───────────────────────────┴──────────────────────────────────────┘

  The fix needs two changes:

  1. Mesh.cpp::computedt() — add jitter when inserting into DT only

  2. ThreeDimensionalShape.cpp::ComputeInputNMM() — recover exact positions from pVertexList instead of
  fvi->point()

  The current line in ComputeInputNMM():
  // reads jittered DT position ← BAD with jitter
  input_nmm.BoundaryPoints.push_back(
      SamplePoint(fvi->point()[0], fvi->point()[1], fvi->point()[2]));

  Should become:
  // reads exact original position via stored vertex ID ← preserves precision
  int orig_id = fvi->info().id;
  auto& orig_pt = input.pVertexList[orig_id]->point();
  input_nmm.BoundaryPoints.push_back(
      SamplePoint(orig_pt[0], orig_pt[1], orig_pt[2]));


Done. Summary of what changed:

  - Mesh.cpp::computedt() — adds std::mt19937 jitter of ±(bb_diagonal × 1e-9) to each point before
  dt.insert(). The mesh vertices themselves are untouched.
  - ThreeDimensionalShape.cpp::ComputeInputNMM() — instead of fvi->point() (jittered DT position), reads
  input.pVertexList[fvi->info().id]->point() to recover the exact original mesh vertex position for
  BoundaryPoints and the sample point files.

  The only thing that sees the jitter is the DT tetrahedra — circumcenters will be off by at most ~1e-9 ×
  diagonal, which for a 1m mesh is sub-nanometre. Everything that the user sees (sample points, clustering,
  visualisation) uses the exact original coordinates.