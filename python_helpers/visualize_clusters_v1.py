"""
visualize_clusters.py  –  Visualise MAT clustering induced by boundary clusters.

Usage
-----
    python visualize_clusters.py <boundary_clusters.txt> <vertex_samples.txt> <file.ma> [mesh_file]

    # or edit the constants below and run without arguments.

Algorithm
---------
    1. Load boundary clusters  →  sample_id → cluster_id  mapping.
    2. Load vertex-samples     →  for each MAT vertex, the set of associated sample IDs.
    3. For each MAT vertex compute its *cluster-set*:
           cluster_set(v) = frozenset{ cluster_id(s)  for s in samples(v) }
    4. A MAT edge (v1, v2) is a *same-cluster edge* iff
           cluster_set(v1) == cluster_set(v2)   (strict set equality)
    5. Run Union-Find on same-cluster edges  →  connected components.
    6. Visualise with Polyscope:
           • MAT curve network  – edges coloured by component; cut edges dark grey.
           • MAT vertex cloud   – coloured by component.
           • Boundary clusters  – one point-cloud per cluster (hidden by default).
           • Input mesh         – semi-transparent surface (if provided).
           • Vertex selection   – click index to highlight vertex + its boundary samples.

File formats expected
---------------------
    boundary_clusters.txt  :  # comment
                               <num_clusters>  <num_points>
                               <cluster_id>  <vertex_index>  x y z  nx ny nz  …

    vertex_samples.txt     :  # comment
                               <total_count>
                               v <vertex_idx>  cx cy cz  radius  <num_samples>
                               s <sample_id>  px py pz
                               …

    *.ma                   :  <n_verts> <n_edges> <n_faces>
                               v  x y z radius
                               e  v1 v2
                               f  v1 v2 v3

    mesh file              :  any format trimesh supports (.obj/.off/.ply/…) — optional 4th argument
"""

import sys
import os
import numpy as np


# ── hard-coded fallbacks (edit if running without CLI arguments) ──────────────
CLUSTER_FILE        = r"cube_subdiv_500_boundary_clusters.txt"
VERTEX_SAMPLES_FILE = r"cube_subdiv_500_vertex_samples.txt"
MA_FILE             = r"cube_subdiv_500.ma"
MESH_FILE           = r""   # optional surface mesh (.obj / .off / .ply)


# ─────────────────────────────────────────────────────────────────────────────
# I/O – MAT files
# ─────────────────────────────────────────────────────────────────────────────

def load_clusters(filepath):
    """
    Parse a *_boundary_clusters.txt file.

    Returns
    -------
    clusters : dict  { cluster_id (int) -> dict with keys
                        'positions'  (N,3) float64,
                        'normals'    (N,3) float64,
                        'indices'    (N,)  int32  }
    meta     : dict  { 'num_clusters', 'num_points', 'filepath' }
    """
    clusters = {}
    num_clusters = num_points = 0
    header_read = False

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if not header_read:
                num_clusters = int(parts[0])
                num_points   = int(parts[1])
                header_read  = True
                continue
            # data row: cluster_id  vertex_index  x  y  z  nx  ny  nz
            cid = int(parts[0])
            vi  = int(parts[1])
            x, y, z     = float(parts[2]), float(parts[3]), float(parts[4])
            nx, ny, nz  = float(parts[5]), float(parts[6]), float(parts[7])
            if cid not in clusters:
                clusters[cid] = {'positions': [], 'normals': [], 'indices': []}
            clusters[cid]['positions'].append([x, y, z])
            clusters[cid]['normals'].append([nx, ny, nz])
            clusters[cid]['indices'].append(vi)

    for cid, d in clusters.items():
        d['positions'] = np.array(d['positions'], dtype=np.float64)
        d['normals']   = np.array(d['normals'],   dtype=np.float64)
        d['indices']   = np.array(d['indices'],   dtype=np.int32)

    meta = dict(num_clusters=num_clusters, num_points=num_points,
                filepath=filepath)
    return clusters, meta


def load_vertex_samples(filepath):
    """
    Parse a *_vertex_samples.txt file.

    Returns
    -------
    vertex_indices  : (M,)   int32
    vertex_centers  : (M, 3) float64
    vertex_bplists  : (M,)   object  – per-vertex int32 array of sample IDs
    """
    v_indices  = []
    v_centers  = []
    v_bplists  = []
    cur_ids    = None

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) == 1:          # count line
                continue
            if parts[0] == 'v':
                if cur_ids is not None:
                    v_bplists.append(np.array(cur_ids, dtype=np.int32))
                v_indices.append(int(parts[1]))
                v_centers.append([float(parts[2]), float(parts[3]), float(parts[4])])
                cur_ids = []
            elif parts[0] == 's':
                cur_ids.append(int(parts[1]))

    if cur_ids is not None:
        v_bplists.append(np.array(cur_ids, dtype=np.int32))

    vertex_indices = np.array(v_indices, dtype=np.int32)
    vertex_centers = np.array(v_centers, dtype=np.float64)
    bplists_arr    = np.empty(len(v_bplists), dtype=object)
    for i, bl in enumerate(v_bplists):
        bplists_arr[i] = bl

    return vertex_indices, vertex_centers, bplists_arr


def load_ma_edges(path):
    """
    Load vertex positions and edges from a .ma file.

    Returns
    -------
    coords : (V, 3) float64
    edges  : (E, 2) int32
    """
    coords = []
    edges  = []
    with open(path) as f:
        lines = f.readlines()
    for line in lines[1:]:
        parts = line.split()
        if not parts:
            continue
        if parts[0] == 'v':
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        elif parts[0] == 'e':
            edges.append([int(parts[1]), int(parts[2])])
    return (np.array(coords, dtype=np.float64),
            np.array(edges,  dtype=np.int32) if edges else np.empty((0, 2), dtype=np.int32))


# ─────────────────────────────────────────────────────────────────────────────
# I/O – surface mesh
# ─────────────────────────────────────────────────────────────────────────────

def load_mesh(path):
    """
    Load a surface mesh file via trimesh (triangulates automatically).

    Returns
    -------
    verts : (V, 3) float64
    faces : (F, 3) int32
    """
    import trimesh
    mesh = trimesh.load(path, force='mesh')
    return mesh.vertices.astype(np.float64), mesh.faces.astype(np.int32)


# ─────────────────────────────────────────────────────────────────────────────
# Cluster-set logic
# ─────────────────────────────────────────────────────────────────────────────

def build_mat_cluster_sets(clusters, vertex_bplists):
    """
    For each MAT vertex compute the frozenset of cluster IDs of its samples.

    A sample that does not appear in any cluster gets the sentinel ID -1.

    Returns
    -------
    sample_to_cluster : dict  { sample_id -> cluster_id }
    mat_sets          : list  of frozensets, one per MAT vertex
    """
    sample_to_cluster = {}
    for cid, d in clusters.items():
        for vi in d['indices']:
            sample_to_cluster[int(vi)] = cid

    mat_sets = []
    for bplist in vertex_bplists:
        cset = frozenset(sample_to_cluster.get(int(s), -1) for s in bplist)
        mat_sets.append(cset)

    return sample_to_cluster, mat_sets


def find_components(n, edges, mat_sets):
    """
    Union-Find on same-cluster edges (strict cluster-set equality).

    Returns
    -------
    components : (n,)   int32  – component ID per MAT vertex
    edge_same  : (E,)   bool   – True when the edge connects same-cluster vertices
    n_comps    : int           – total number of components
    """
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        pa, pb = find(a), find(b)
        if pa != pb:
            parent[pa] = pb

    edge_same = []
    for v1, v2 in edges:
        same = (mat_sets[v1] == mat_sets[v2])
        edge_same.append(same)
        if same:
            union(v1, v2)

    root_to_comp = {}
    comp_id = 0
    components = np.zeros(n, dtype=np.int32)
    for i in range(n):
        r = find(i)
        if r not in root_to_comp:
            root_to_comp[r] = comp_id
            comp_id += 1
        components[i] = root_to_comp[r]

    return components, np.array(edge_same, dtype=bool), comp_id


# ─────────────────────────────────────────────────────────────────────────────
# Colour palette  (cycles for large counts)
# ─────────────────────────────────────────────────────────────────────────────

def _palette(n):
    """Return n visually distinct RGB colours as (n, 3) float32."""
    base = np.array([
        [0.93, 0.17, 0.17],  # red
        [0.17, 0.53, 0.93],  # blue
        [0.17, 0.82, 0.17],  # green
        [0.93, 0.73, 0.10],  # yellow
        [0.63, 0.17, 0.93],  # purple
        [0.10, 0.82, 0.82],  # cyan
        [0.93, 0.47, 0.10],  # orange
        [0.90, 0.17, 0.60],  # pink
        [0.50, 0.82, 0.17],  # lime
        [0.17, 0.40, 0.40],  # teal
    ], dtype=np.float32)
    if n <= len(base):
        return base[:n]
    rng   = np.random.default_rng(42)
    tiled = np.tile(base, (n // len(base) + 1, 1))[:n]
    jitter = rng.uniform(-0.12, 0.12, tiled.shape).astype(np.float32)
    return np.clip(tiled + jitter, 0.0, 1.0)


# ─────────────────────────────────────────────────────────────────────────────
# Polyscope visualisation
# ─────────────────────────────────────────────────────────────────────────────

def visualize(clusters, meta,
              vertex_centers, vertex_bplists,
              ma_edges, mat_sets,
              components, edge_same, n_comps,
              mesh_coords=None, mesh_faces=None):
    """
    Interactive Polyscope viewer.

    Structures
    ----------
    MAT Edges        – curve network; same-cluster edges coloured by component,
                       cut edges in dark grey.
    MAT Vertices     – point cloud coloured by component.
    Cluster N (n=…)  – boundary point-cloud per surface cluster (hidden by default).
    Input Mesh       – semi-transparent surface mesh (if provided).
    Selected Vertex  – highlighted when a vertex index is chosen in the panel.
    Responsible Pts  – boundary samples that define the selected MAT sphere.

    Side panel
    ----------
    Shows edge/component statistics, the cluster-set signature of each
    connected component, and an interactive vertex selector.
    """
    try:
        import polyscope as ps
        import polyscope.imgui as psim
    except ImportError:
        print("Polyscope not installed.  Run:  pip install polyscope")
        return

    ps.init()
    ps.set_program_name("MAT Cluster Viewer")
    ps.set_ground_plane_mode("none")
    ps.set_background_color((0.12, 0.12, 0.16))

    n_palette = max(n_comps, meta['num_clusters'], 1)
    colours   = _palette(n_palette)

    COLOR_CUT = np.array([0.18, 0.18, 0.18], dtype=np.float32)  # dark grey

    # ── build sample_id → 3-D position lookup ─────────────────────────────────
    sample_positions = {}
    for cid, d in clusters.items():
        for vi, pos in zip(d['indices'], d['positions']):
            sample_positions[int(vi)] = pos

    # ── input surface mesh ────────────────────────────────────────────────────
    if mesh_coords is not None and mesh_faces is not None and len(mesh_faces) > 0:
        sm = ps.register_surface_mesh("Input Mesh", mesh_coords, mesh_faces)
        sm.set_color((0.55, 0.70, 0.85))
        sm.set_transparency(0.35)
        sm.set_smooth_shade(True)

    # ── MAT curve network ─────────────────────────────────────────────────────
    if len(ma_edges) > 0:
        edge_colors = np.array(
            [colours[components[v1] % len(colours)] if same else COLOR_CUT
             for (v1, v2), same in zip(ma_edges, edge_same)],
            dtype=np.float32)

        cn = ps.register_curve_network("MAT Edges", vertex_centers, ma_edges)
        cn.add_color_quantity("component", edge_colors,
                              defined_on='edges', enabled=True)
        cn.set_radius(0.0008, relative=True)

    # ── MAT vertex point cloud ────────────────────────────────────────────────
    vertex_colors = np.array(
        [colours[c % len(colours)] for c in components], dtype=np.float32)
    pc_mat = ps.register_point_cloud("MAT Vertices", vertex_centers)
    pc_mat.add_color_quantity("component", vertex_colors, enabled=True)
    pc_mat.set_radius(0.0012, relative=True)

    # ── Boundary cluster point-clouds (hidden by default) ────────────────────
    for cid, d in sorted(clusters.items()):
        pts = d['positions']
        col = tuple(colours[cid % len(colours)].tolist())
        pc  = ps.register_point_cloud(f"Cluster {cid}  (n={len(pts)})", pts)
        pc.set_color(col)
        pc.set_radius(0.0015, relative=True)
        pc.set_enabled(False)

    # ── pre-compute component info for the panel ──────────────────────────────
    comp_sizes    = {}
    comp_cset_rep = {}   # one representative cluster-set per component
    for i, c in enumerate(components):
        comp_sizes[c]    = comp_sizes.get(c, 0) + 1
        comp_cset_rep[c] = mat_sets[i]

    fname   = os.path.basename(meta['filepath'])
    n_same  = int(edge_same.sum())
    n_cut   = len(edge_same) - n_same
    n_verts = len(vertex_centers)

    # ── selection state (mutable via single-element lists for closure) ─────────
    state = {
        "sel_changed": [False],
        "sel_idx":     [-1],     # -1 = nothing selected yet
        "last_pick":   [None],   # (structure_name, elem_idx) of previous pick
    }

    def _refresh_selection(idx):
        """Register/update the two highlight clouds for the selected vertex."""
        cx, cy, cz = vertex_centers[idx]
        pc_sel = ps.register_point_cloud(
            "Selected MAT Vertex", np.array([[cx, cy, cz]]))
        pc_sel.set_color((0.0, 1.0, 0.53))
        pc_sel.set_radius(0.0025, relative=True)

        # collect positions of all boundary samples for this vertex
        bp_ids = vertex_bplists[idx]
        bp_pos = np.array(
            [sample_positions[int(sid)]
             for sid in bp_ids if int(sid) in sample_positions],
            dtype=np.float64)

        if len(bp_pos) > 0:
            pc_bp = ps.register_point_cloud("Responsible Boundary Pts", bp_pos)
            pc_bp.set_color((1.0, 0.13, 0.27))
            pc_bp.set_radius(0.003, relative=True)
        else:
            try:
                ps.remove_point_cloud("Responsible Boundary Pts")
            except Exception:
                pass

    # ── ImGui callback ────────────────────────────────────────────────────────
    def callback():
        psim.PushItemWidth(340)
        psim.TextUnformatted("── MAT Cluster Viewer ────────────────────────")
        psim.Separator()
        psim.TextUnformatted(f"Cluster file     : {fname}")
        psim.TextUnformatted(f"Surface clusters : {meta['num_clusters']}")
        psim.TextUnformatted(f"MAT vertices     : {n_verts}")
        psim.TextUnformatted(f"MAT edges        : {len(ma_edges)}")
        psim.Separator()
        psim.TextUnformatted(f"Same-cluster edges : {n_same}")
        psim.TextUnformatted(f"Cut edges          : {n_cut}")
        psim.TextUnformatted(f"MAT components     : {n_comps}")
        psim.Separator()
        psim.TextUnformatted("Components  (verts | cluster-set signature):")
        for cid_key in sorted(comp_sizes.keys()):
            sig = sorted(comp_cset_rep[cid_key])
            psim.TextUnformatted(
                f"  [{cid_key:3d}]  {comp_sizes[cid_key]:5d} verts  |  {sig}")

        # ── pick detection (click on MAT Vertices in the viewport) ───────────
        psim.Separator()
        psim.TextUnformatted("── Vertex Inspector ──────────────────────────")
        psim.TextUnformatted("Left-click (no drag) a point in 'MAT Vertices'.")

        try:
            if ps.have_selection():
                pick = ps.get_selection()
                if pick is not None and pick.is_hit:
                    pick_name = pick.structure_name
                    pick_idx  = int(pick.local_index)
                    pick_key  = (pick_name, pick_idx)
                    if pick_key != state["last_pick"][0]:
                        state["last_pick"][0] = pick_key
                        if pick_name == "MAT Vertices" and 0 <= pick_idx < n_verts:
                            state["sel_idx"][0]     = pick_idx
                            state["sel_changed"][0] = True
        except Exception:
            pass

        # refresh highlight when needed
        if state["sel_changed"][0]:
            state["sel_changed"][0] = False
            _refresh_selection(state["sel_idx"][0])

        # show info for the current selection
        idx = state["sel_idx"][0]
        if idx < 0:
            psim.TextUnformatted("  (no vertex selected)")
        else:
            cx, cy, cz = vertex_centers[idx]
            psim.TextUnformatted(f"  Index      : {idx}")
            psim.TextUnformatted(f"  Position   : ({cx:.4f}, {cy:.4f}, {cz:.4f})")
            psim.TextUnformatted(f"  Component  : {components[idx]}")
            psim.TextUnformatted(f"  Cluster-set: {sorted(mat_sets[idx])}")
            bp_ids = vertex_bplists[idx]
            psim.TextUnformatted(f"  Boundary pts (responsible): {len(bp_ids)}")

        psim.Separator()
        psim.TextUnformatted("Tip: enable 'Cluster N' clouds to compare")
        psim.TextUnformatted("     surface clusters with MAT components.")
        if mesh_coords is not None:
            psim.TextUnformatted("Tip: toggle 'Input Mesh' transparency in")
            psim.TextUnformatted("     its Polyscope panel entry.")
        psim.PopItemWidth()

    ps.set_user_callback(callback)
    ps.show()


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    cluster_file        = sys.argv[1] if len(sys.argv) > 1 else CLUSTER_FILE
    vertex_samples_file = sys.argv[2] if len(sys.argv) > 2 else VERTEX_SAMPLES_FILE
    ma_file             = sys.argv[3] if len(sys.argv) > 3 else MA_FILE
    mesh_file           = sys.argv[4] if len(sys.argv) > 4 else MESH_FILE

    if not cluster_file or not vertex_samples_file or not ma_file:
        print("Usage:  python visualize_clusters.py "
              "<boundary_clusters.txt> <vertex_samples.txt> <file.ma> [mesh_file]")
        sys.exit(1)

    for f in [cluster_file, vertex_samples_file, ma_file]:
        if not os.path.exists(f):
            print(f"File not found: {f}")
            sys.exit(1)

    # ── load MAT files ────────────────────────────────────────────────────────
    print(f"Loading clusters     : {cluster_file}")
    clusters, meta = load_clusters(cluster_file)
    print(f"  {meta['num_clusters']} clusters, {meta['num_points']} boundary points")

    print(f"Loading vertex samples: {vertex_samples_file}")
    vertex_indices, vertex_centers, vertex_bplists = \
        load_vertex_samples(vertex_samples_file)
    print(f"  {len(vertex_centers)} MAT vertices")

    print(f"Loading MA file      : {ma_file}")
    _, ma_edges = load_ma_edges(ma_file)
    print(f"  {len(ma_edges)} MAT edges")

    # ── load optional surface mesh ─────────────────────────────────────────────
    mesh_coords = mesh_faces = None
    if mesh_file:
        if not os.path.exists(mesh_file):
            print(f"Mesh file not found (skipping): {mesh_file}")
        else:
            print(f"Loading mesh         : {mesh_file}")
            try:
                mesh_coords, mesh_faces = load_mesh(mesh_file)
                print(f"  {len(mesh_coords)} vertices, {len(mesh_faces)} faces")
            except Exception as exc:
                print(f"  Failed to load mesh: {exc}")

    # ── cluster-set assignment + connected components ─────────────────────────
    sample_to_cluster, mat_sets = build_mat_cluster_sets(clusters, vertex_bplists)

    n_unclustered = sum(1 for s in mat_sets if -1 in s)
    if n_unclustered:
        print(f"  Warning: {n_unclustered} MAT vertices have samples "
              f"not found in any cluster (mapped to sentinel -1)")

    components, edge_same, n_comps = find_components(
        len(vertex_centers), ma_edges, mat_sets)
    n_same = int(edge_same.sum())
    n_cut  = len(edge_same) - n_same
    print(f"  {n_comps} MAT components  "
          f"({n_same} same-cluster edges, {n_cut} cut edges)")

    # ── visualise ─────────────────────────────────────────────────────────────
    visualize(clusters, meta,
              vertex_centers, vertex_bplists,
              ma_edges, mat_sets,
              components, edge_same, n_comps,
              mesh_coords=mesh_coords,
              mesh_faces=mesh_faces)
