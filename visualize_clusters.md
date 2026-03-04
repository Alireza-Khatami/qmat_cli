● Done. Here's a summary of what changed and why:

  New functions added:
  - load_vertex_samples() — inlined from read_qmat_output.py (no trimesh dependency), returns only what's needed
  - load_ma_edges() — lightweight .ma loader, only parses v and e lines
  - build_mat_cluster_sets() — builds sample_to_cluster dict, then for each MAT vertex computes frozenset{cluster_id(s) for s in samples(v)}
  - find_components() — union-find on same-cluster edges (strict set equality as you defined)

  Visualization:
  - MAT curve network (MAT Edges): same-cluster edges colored by connected component, cut edges dark grey (0.18, 0.18, 0.18)
  - MAT vertex cloud (MAT Vertices): colored by component (matching edge colors)
  - Boundary clusters (Cluster N): still registered but hidden by default — you can toggle them on in Polyscope's GUI to visually compare surface clusters
  with MAT components

  Side panel shows: component count, same vs cut edge counts, and each component's cluster-set signature (e.g. [0, 2] meaning those MAT vertices sit between
   surface clusters 0 and 2)

  Entry point now takes 3 args or falls back to the hardcoded constants at the top of the file.

✻ Brewed for 3m 51s