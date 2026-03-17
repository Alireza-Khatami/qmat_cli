
  Removed entirely:
  - load_clusters(), build_mat_cluster_sets(), find_components(), export_cluster_ply(), _palette()
  - CLUSTER_FILE constant and CLI arg
  - All CSet point clouds (one-per-cluster-set grouping)
  - Cluster-colored edges
  - Boundary cluster point-clouds display
  - All cluster stats from the ImGui panel

  Updated:
  - load_vertex_samples() now also parses the s line positions — so boundary point 3D positions come directly from
  _vertex_samples.txt instead of needing the cluster file
  - CLI args are now <vertex_samples.txt> <file.ma> [mesh_file]

  Kept / improved:
  - load_mat_topo(), load_voronoi_neighbors(), can_merge() unchanged
  - On click: green = selected vertex, red = its boundary points, gold = all mergeable candidates, orange = combined
  boundary points of all mergeable candidates
  - ImGui panel shows topo flags, nmn_bplist count, mergeable count, sidecar load status

✻ Crunched for 2m 29s