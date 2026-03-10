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
    7. Export button in side panel:
           • <stem>_mat.obj        – MAT vertices (x y z), edges (l), faces (f).
           • <stem>_clusters.ply   – MAT vertices coloured by cluster-set (per-vertex RGB).

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
import colorsys
import numpy as np


# ── hard-coded fallbacks (edit if running without CLI arguments) ──────────────
CLUSTER_FILE        = rf"cube_subdiv_500_cluster_face/cube_subdiv_500_boundary_clusters.txt"
VERTEX_SAMPLES_FILE = rf"cube_subdiv_500_cluster_face/cube_subdiv_500_vertex_samples.txt"
MA_FILE             = rf"cube_subdiv_500_cluster_face/cube_subdiv_500.ma"
MESH_FILE           = rf""   # optional surface mesh (.obj / .off / .ply)


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


def load_ma(path):
    """
    Load vertex positions, edges, and faces from a .ma file.

    Returns
    -------
    coords : (V, 3) float64
    edges  : (E, 2) int32
    faces  : (F, 3) int32
    """
    coords = []
    edges  = []
    faces  = []
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
        elif parts[0] == 'f':
            faces.append([int(parts[1]), int(parts[2]), int(parts[3])])
    return (
        np.array(coords, dtype=np.float64),
        np.array(edges,  dtype=np.int32) if edges else np.empty((0, 2), dtype=np.int32),
        np.array(faces,  dtype=np.int32) if faces else np.empty((0, 3), dtype=np.int32),
    )


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
# Colour palette  – N distinct HSV colours, no duplicates
# ─────────────────────────────────────────────────────────────────────────────

def _palette(n):
    """Return n visually distinct RGB colours as (n, 3) float32 using HSV."""
    if n == 0:
        return np.empty((0, 3), dtype=np.float32)
    colours = []
    for i in range(n):
        hue = i / n                              # evenly spread around the wheel
        sat = 0.85 if i % 2 == 0 else 0.60      # alternate saturation
        val = 0.92 if (i // 2) % 2 == 0 else 0.72  # alternate brightness
        r, g, b = colorsys.hsv_to_rgb(hue, sat, val)
        colours.append([r, g, b])
    return np.array(colours, dtype=np.float32)


# ─────────────────────────────────────────────────────────────────────────────
# Export helpers
# ─────────────────────────────────────────────────────────────────────────────

def export_mat_obj(path, coords, edges, faces):
    """
    Write MAT geometry to a Wavefront OBJ file.

    Vertices carry only x y z (radius dropped).
    Edges are written as OBJ line elements (l).
    Faces are written as OBJ face elements (f).
    All indices are 1-based as required by the format.
    """
    with open(path, 'w') as f:
        f.write("# MAT exported by visualize_clusters_v2.py\n")
        f.write(f"# {len(coords)} vertices  {len(edges)} edges  {len(faces)} faces\n\n")
        for x, y, z in coords:
            f.write(f"v {x:.10g} {y:.10g} {z:.10g}\n")
        if len(edges):
            f.write("\n")
            for v1, v2 in edges:
                f.write(f"l {v1 + 1} {v2 + 1}\n")
        if len(faces):
            f.write("\n")
            for v1, v2, v3 in faces:
                f.write(f"f {v1 + 1} {v2 + 1} {v3 + 1}\n")
    print(f"Exported MAT OBJ  → {path}")


def export_cluster_ply(path, vertex_centers, mat_sets, cset_to_color_idx, colours):
    """
    Write MAT vertices as an ASCII PLY point cloud with per-vertex RGB colour.

    Each vertex is coloured by its cluster-set group.
    MeshLab will display this as a coloured point cloud.
    """
    n = len(vertex_centers)
    with open(path, 'w') as f:
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write("comment MAT cluster colours exported by visualize_clusters_v2.py\n")
        f.write(f"element vertex {n}\n")
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write("property uchar red\n")
        f.write("property uchar green\n")
        f.write("property uchar blue\n")
        f.write("end_header\n")
        for pos, cset in zip(vertex_centers, mat_sets):
            cidx = cset_to_color_idx[cset] % len(colours)
            col  = colours[cidx]
            r = int(np.clip(col[0] * 255, 0, 255))
            g = int(np.clip(col[1] * 255, 0, 255))
            b = int(np.clip(col[2] * 255, 0, 255))
            f.write(f"{pos[0]:.10g} {pos[1]:.10g} {pos[2]:.10g} {r} {g} {b}\n")
    print(f"Exported cluster PLY → {path}")


# ─────────────────────────────────────────────────────────────────────────────
# Polyscope visualisation
# ─────────────────────────────────────────────────────────────────────────────

def visualize(clusters, meta,
              vertex_centers, vertex_bplists,
              ma_coords, ma_edges, ma_faces, mat_sets,
              mesh_coords=None, mesh_faces=None,
              ma_file=""):
    """
    Interactive Polyscope viewer.

    Structures
    ----------
    CSet [a,b,…]     – one point cloud per unique cluster-set frozenset,
                       each with its own colour.  Edge connectivity is ignored
                       for grouping: two vertices belong to the same cloud iff
                       their cluster-set is identical.
    MAT Edges        – curve network; same-cluster-set edges keep the group
                       colour, cross-group edges are dark grey.
    Cluster N (n=…)  – boundary point-cloud per surface cluster (hidden by default).
    Input Mesh       – semi-transparent surface mesh (if provided).
    Selected Vertex  – highlighted on click; shows responsible boundary pts.

    Side panel
    ----------
    Shows group / edge statistics and an interactive vertex inspector.
    Export button writes <stem>_mat.obj and <stem>_clusters.ply.
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

    COLOR_CUT = np.array([0.18, 0.18, 0.18], dtype=np.float32)

    # ── build sample_id → 3-D position lookup ─────────────────────────────────
    sample_positions = {}
    for cid, d in clusters.items():
        for vi, pos in zip(d['indices'], d['positions']):
            sample_positions[int(vi)] = pos

    # ── group MAT vertices by cluster-set (no edge connectivity used) ─────────
    cset_to_verts = {}
    for i, cset in enumerate(mat_sets):
        cset_to_verts.setdefault(cset, []).append(i)

    # deterministic order: sort each frozenset, then sort the list of tuples
    cset_list = sorted(cset_to_verts.keys(), key=lambda s: tuple(sorted(s)))
    cset_to_color_idx = {cset: i for i, cset in enumerate(cset_list)}
    n_groups  = len(cset_list)
    n_verts   = len(vertex_centers)

    # One distinct colour per cluster-set group (at least as many as clusters)
    colours = _palette(max(n_groups, meta['num_clusters'], 1))

    # ── input surface mesh ────────────────────────────────────────────────────
    if mesh_coords is not None and mesh_faces is not None and len(mesh_faces) > 0:
        sm = ps.register_surface_mesh("Input Mesh", mesh_coords, mesh_faces)
        sm.set_color((0.55, 0.70, 0.85))
        sm.set_transparency(0.35)
        sm.set_smooth_shade(True)

    # ── one point cloud per cluster-set group ─────────────────────────────────
    # cset_cloud_info maps cloud_name → global_index_array  (for pick resolution)
    cset_cloud_info = {}
    cset_list_ordered = sorted(cset_list, key=lambda s: len(cset_to_verts[s]), reverse=True)
    for cset in cset_list_ordered:
        global_ids = np.array(cset_to_verts[cset], dtype=np.int32)
        cidx       = cset_to_color_idx[cset]
        col        = tuple(colours[cidx % len(colours)].tolist())
        name       = f"CSet {sorted(cset)} n={len(global_ids)}"
        if len(global_ids) >= 40:
            should_be_enabled = False       
        else:
            should_be_enabled = True   
        pc = ps.register_point_cloud(name, vertex_centers[global_ids], enabled=True)
        pc.set_color(col)
        # if should_be_enabled:
        #     pc.set_radius(0.0050, relative=True)
        # else:
        pc.set_radius(0.0012, relative=True)
        cset_cloud_info[name] = global_ids

    # ── MAT curve network (coloured by cluster-set; cut edges dark grey) ──────
    if len(ma_edges) > 0:
        edge_colors = np.array(
            [colours[cset_to_color_idx[mat_sets[v1]] % len(colours)]
             if mat_sets[v1] == mat_sets[v2] else COLOR_CUT
             for v1, v2 in ma_edges],
            dtype=np.float32)
        cn = ps.register_curve_network("MAT Edges", vertex_centers, ma_edges)
        cn.add_color_quantity("cluster-set", edge_colors,
                              defined_on='edges', enabled=True)
        cn.set_radius(0.0008, relative=True)
        n_same = int(sum(mat_sets[v1] == mat_sets[v2] for v1, v2 in ma_edges))
        n_cut  = len(ma_edges) - n_same
    else:
        n_same = n_cut = 0

    # ── Boundary cluster point-clouds (hidden by default) ────────────────────
    for cid, d in sorted(clusters.items()):
        pts = d['positions']
        col = tuple(colours[cid % len(colours)].tolist())
        pc  = ps.register_point_cloud(f"Cluster {cid}  (n={len(pts)})", pts)
        pc.set_color(col)
        pc.set_radius(0.0015, relative=True)
        pc.set_enabled(False)

    fname = os.path.basename(meta['filepath'])

    # ── export path helpers ───────────────────────────────────────────────────
    if ma_file:
        ma_dir  = os.path.dirname(os.path.abspath(ma_file))
        ma_stem = os.path.splitext(os.path.basename(ma_file))[0]
    else:
        ma_dir  = os.getcwd()
        ma_stem = "mat_export"
    obj_export_path = os.path.join(ma_dir, ma_stem + "_mat.obj")
    ply_export_path = os.path.join(ma_dir, ma_stem + "_clusters.ply")

    # ── selection state ───────────────────────────────────────────────────────
    state = {
        "sel_changed":  [False],
        "sel_idx":      [-1],
        "last_pick":    [None],
        "export_done":  [False],
        "export_msg":   [""],
    }

    def _refresh_selection(idx):
        cx, cy, cz = vertex_centers[idx]
        pc_sel = ps.register_point_cloud(
            "Selected MAT Vertex", np.array([[cx, cy, cz]]))
        pc_sel.set_color((0.0, 1.0, 0.53))
        pc_sel.set_radius(0.0025, relative=True)

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
        psim.TextUnformatted("── MAT Cluster Viewer (by cluster-set) ───────")
        psim.Separator()
        psim.TextUnformatted(f"Cluster file     : {fname}")
        psim.TextUnformatted(f"Surface clusters : {meta['num_clusters']}")
        psim.TextUnformatted(f"MAT vertices     : {n_verts}")
        psim.TextUnformatted(f"MAT edges        : {len(ma_edges)}")
        psim.TextUnformatted(f"MAT faces        : {len(ma_faces)}")
        psim.Separator()
        psim.TextUnformatted(f"Unique cluster-sets : {n_groups}")
        psim.TextUnformatted(f"Same-cset edges     : {n_same}")
        psim.TextUnformatted(f"Cross-cset edges    : {n_cut}")
        psim.Separator()
        psim.TextUnformatted("Groups  (verts | cluster-set):")
        for cset in cset_list:
            cidx = cset_to_color_idx[cset]
            psim.TextUnformatted(
                f"  [{cidx:3d}]  {len(cset_to_verts[cset]):5d} verts  |  {sorted(cset)}")

        # ── pick detection ────────────────────────────────────────────────────
        psim.Separator()
        psim.TextUnformatted("── Vertex Inspector ──────────────────────────")
        psim.TextUnformatted("Left-click (no drag) any CSet cloud to inspect.")

        try:
            if ps.have_selection():
                pick = ps.get_selection()
                if pick is not None and pick.is_hit:
                    pick_name  = pick.structure_name
                    pick_local = int(pick.local_index)
                    pick_key   = (pick_name, pick_local)
                    if pick_key != state["last_pick"][0]:
                        state["last_pick"][0] = pick_key
                        if pick_name in cset_cloud_info:
                            global_ids = cset_cloud_info[pick_name]
                            if 0 <= pick_local < len(global_ids):
                                state["sel_idx"][0]     = int(global_ids[pick_local])
                                state["sel_changed"][0] = True
        except Exception:
            pass

        if state["sel_changed"][0]:
            state["sel_changed"][0] = False
            _refresh_selection(state["sel_idx"][0])

        idx = state["sel_idx"][0]
        if idx < 0:
            psim.TextUnformatted("  (no vertex selected)")
        else:
            cx, cy, cz = vertex_centers[idx]
            psim.TextUnformatted(f"  Index      : {idx}")
            psim.TextUnformatted(f"  Position   : ({cx:.4f}, {cy:.4f}, {cz:.4f})")
            psim.TextUnformatted(f"  Cluster-set: {sorted(mat_sets[idx])}")
            bp_ids = vertex_bplists[idx]
            psim.TextUnformatted(f"  Boundary pts (responsible): {len(bp_ids)}")

        # ── export ────────────────────────────────────────────────────────────
        psim.Separator()
        psim.TextUnformatted("── Export ────────────────────────────────────")
        if psim.Button("Export MAT (.obj) + Clusters (.ply)"):
            try:
                export_mat_obj(obj_export_path, ma_coords, ma_edges, ma_faces)
                export_cluster_ply(ply_export_path, vertex_centers,
                                   mat_sets, cset_to_color_idx, colours)
                state["export_done"][0] = True
                state["export_msg"][0]  = (
                    f"OK  {os.path.basename(obj_export_path)}\n"
                    f"    {os.path.basename(ply_export_path)}"
                )
            except Exception as exc:
                state["export_done"][0] = True
                state["export_msg"][0]  = f"Export failed: {exc}"

        if state["export_done"][0]:
            psim.TextUnformatted(state["export_msg"][0])

        psim.Separator()
        psim.TextUnformatted("Tip: enable 'Cluster N' clouds to compare")
        psim.TextUnformatted("     surface clusters with MAT groups.")
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
    ma_coords, ma_edges, ma_faces = load_ma(ma_file)
    print(f"  {len(ma_coords)} vertices  {len(ma_edges)} edges  {len(ma_faces)} faces")

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

    # ── cluster-set assignment ────────────────────────────────────────────────
    sample_to_cluster, mat_sets = build_mat_cluster_sets(clusters, vertex_bplists)

    n_unclustered = sum(1 for s in mat_sets if -1 in s)
    if n_unclustered:
        print(f"  Warning: {n_unclustered} MAT vertices have samples "
              f"not found in any cluster (mapped to sentinel -1)")

    n_unique = len({frozenset(s) for s in mat_sets})
    print(f"  {n_unique} unique cluster-sets across {len(vertex_centers)} MAT vertices")

    # ── visualise ─────────────────────────────────────────────────────────────
    visualize(clusters, meta,
              vertex_centers, vertex_bplists,
              ma_coords, ma_edges, ma_faces, mat_sets,
              mesh_coords=mesh_coords,
              mesh_faces=mesh_faces,
              ma_file=ma_file)
