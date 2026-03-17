"""
visualize_voronoi_mergeability.py  –  Visualise which MAT vertices can merge.

Usage
-----
    python visualize_clusters_vornoi_neighbours.py <vertex_samples.txt> <file.ma> [mesh_file]

    # or edit the constants below and run without arguments.

Algorithm
---------
    1. Load vertex-samples  →  per-MAT-vertex: center, radius, boundary-point IDs + positions.
    2. Load *_mat_topo.txt  →  per-MAT-vertex: topology flags + nmn_bplist (Delaunay tet BPs).
    3. Load *_voronoi_neighbors.txt  →  Delaunay-edge adjacency between boundary points.
    4. Visualise with Polyscope:
           • MAT Vertices   – single grey point cloud (all vertices, clickable).
           • MAT Edges      – curve network.
           • Input Mesh     – semi-transparent surface (if provided).
           • On click:
               – Selected Vertex          (bright green)
               – Its boundary points      (red)
               – All mergeable candidates (gold)  [CanMerge == True]
               – Their boundary points    (orange, combined)
    5. Export button in side panel:
           • <stem>_mat.obj  – MAT vertices (x y z), edges (l), faces (f).

CanMerge conditions (mirror of C++ SlabMesh::CanMerge)
-------------------------------------------------------
    1. Both vertices must be *pure sheet*:
           topo_is_sheet == True  AND  seam == junction == boundary == False
    2. At least one boundary-point pair (one from each vertex's nmn_bplist)
       must be Voronoi neighbors (share a Delaunay edge).
"""

import sys
import os
import numpy as np


# ── hard-coded fallbacks (edit if running without CLI arguments) ──────────────
VERTEX_SAMPLES_FILE = rf"cube_subdiv/500/cube_subdiv_500_vertex_samples.txt"
MA_FILE             = rf"cube_subdiv/500/cube_subdiv_500.ma"
# VERTEX_SAMPLES_FILE = rf"00002362_470c2865307645c78ae7a0cb_trimesh_005/00002362_470c2865307645c78ae7a0cb_trimesh_005_vertex_samples.txt"
# MA_FILE             = rf"00002362_470c2865307645c78ae7a0cb_trimesh_005/00002362_470c2865307645c78ae7a0cb_trimesh_005.ma"
MESH_FILE           = rf""   # optional surface mesh (.obj / .off / .ply)


# ─────────────────────────────────────────────────────────────────────────────
# I/O
# ─────────────────────────────────────────────────────────────────────────────

def load_vertex_samples(filepath):
    """
    Parse a *_vertex_samples.txt file.

    Format:
        # comments
        <total_count>
        v <vertex_idx>  cx cy cz  radius  <num_samples>
        s <sample_id>  px py pz
        ...

    Returns
    -------
    vertex_indices   : (M,)   int32   – actual MAT vertex index per entry
    vertex_centers   : (M, 3) float64 – MAT vertex sphere centre
    vertex_bplists   : (M,)   object  – per-vertex int32 array of sample IDs
    sample_positions : dict { sample_id (int) -> [x, y, z] float64 }
    """
    v_indices  = []
    v_centers  = []
    v_bplists  = []
    cur_ids    = None
    sample_positions = {}

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
                sid = int(parts[1])
                cur_ids.append(sid)
                sample_positions[sid] = [float(parts[2]), float(parts[3]), float(parts[4])]

    if cur_ids is not None:
        v_bplists.append(np.array(cur_ids, dtype=np.int32))

    vertex_indices = np.array(v_indices, dtype=np.int32)
    vertex_centers = np.array(v_centers, dtype=np.float64)
    bplists_arr    = np.empty(len(v_bplists), dtype=object)
    for i, bl in enumerate(v_bplists):
        bplists_arr[i] = bl

    return vertex_indices, vertex_centers, bplists_arr, sample_positions


def load_mat_topo(path):
    """
    Parse a *_mat_topo.txt sidecar file.

    Format (one data line per active MAT vertex):
        <num_active_vertices>
        <idx> <steep> <sheet> <seam> <junction> <topo_boundary> <bp_count>
              <bp_id0> <cl_id0>  <bp_id1> <cl_id1>  ...

    Returns
    -------
    topo_flags  : dict { mat_vertex_idx (int) -> dict with bool keys
                         'steep','sheet','seam','junction','boundary' }
    nmn_bplists : dict { mat_vertex_idx (int) -> set of bp_id (int) }
    """
    topo_flags  = {}
    nmn_bplists = {}
    if not path or not os.path.exists(path):
        return topo_flags, nmn_bplists

    with open(path) as f:
        count = None
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if count is None:
                count = int(parts[0])
                continue
            idx    = int(parts[0])
            steep  = bool(int(parts[1]))
            sheet  = bool(int(parts[2]))
            seam   = bool(int(parts[3]))
            junc   = bool(int(parts[4]))
            tbound = bool(int(parts[5]))
            bp_cnt = int(parts[6])
            bp_set = set()
            offset = 7
            for _ in range(bp_cnt):
                bp_set.add(int(parts[offset]))
                offset += 2   # skip cluster_id
            topo_flags[idx]  = dict(steep=steep, sheet=sheet, seam=seam,
                                    junction=junc, boundary=tbound)
            nmn_bplists[idx] = bp_set

    return topo_flags, nmn_bplists


def load_voronoi_neighbors(path):
    """
    Parse a *_voronoi_neighbors.txt sidecar file.

    Format:
        # comments
        <num_records>
        <bp_id> <neighbor_count> <nb0> <nb1> ...

    Returns
    -------
    voronoi_neighbors : dict { bp_id (int) -> set of neighbor bp_id (int) }
    """
    vn = {}
    if not path or not os.path.exists(path):
        return vn

    with open(path) as f:
        count = None
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if count is None:
                count = int(parts[0])
                continue
            bp_id    = int(parts[0])
            nb_count = int(parts[1])
            nbrs = set()
            for i in range(nb_count):
                nbrs.add(int(parts[2 + i]))
            vn[bp_id] = nbrs

    return vn


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


def load_mesh(path):
    """Load a surface mesh via trimesh. Returns (verts, faces)."""
    import trimesh
    mesh = trimesh.load(path, force='mesh')
    return mesh.vertices.astype(np.float64), mesh.faces.astype(np.int32)


# ─────────────────────────────────────────────────────────────────────────────
# CanMerge  –  Python mirror of C++ SlabMesh::CanMerge
# ─────────────────────────────────────────────────────────────────────────────

def can_merge(local_idx_1, local_idx_2,
              topo_local, nmn_local, voronoi_neighbors):
    """
    Return True iff the two MAT vertices (by local array index) may be merged.

    Conditions (mirror of C++ SlabMesh::CanMerge):
      1. Both must be pure sheet: sheet==True, seam==junction==boundary==False.
      2. At least one bp pair (one from each nmn_bplist) must be Voronoi neighbors.
    """
    if topo_local is None or nmn_local is None or not voronoi_neighbors:
        return False
    f1 = topo_local[local_idx_1]
    f2 = topo_local[local_idx_2]
    if f1 is None or f2 is None:
        return False
    if not (f1['sheet'] and not f1['seam'] and not f1['junction'] and not f1['boundary']):
        return False
    if not (f2['sheet'] and not f2['seam'] and not f2['junction'] and not f2['boundary']):
        return False
    bps1 = nmn_local[local_idx_1]
    bps2 = nmn_local[local_idx_2]
    for bp1 in bps1:
        nbrs = voronoi_neighbors.get(bp1)
        if nbrs is None:
            continue
        for bp2 in bps2:
            if bp2 in nbrs:
                return True
    return False


# ─────────────────────────────────────────────────────────────────────────────
# Export
# ─────────────────────────────────────────────────────────────────────────────

def export_mat_obj(path, coords, edges, faces):
    """Write MAT geometry to a Wavefront OBJ file (1-based indices)."""
    with open(path, 'w') as f:
        f.write("# MAT exported by visualize_clusters_vornoi_neighbours.py\n")
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
    print(f"Exported MAT OBJ → {path}")


# ─────────────────────────────────────────────────────────────────────────────
# Polyscope visualisation
# ─────────────────────────────────────────────────────────────────────────────

def visualize(vertex_indices, vertex_centers, vertex_bplists, sample_positions,
              ma_coords, ma_edges, ma_faces,
              mesh_coords=None, mesh_faces=None,
              ma_file="",
              topo_local=None, nmn_local=None, voronoi_neighbors=None):
    """
    Interactive Polyscope viewer focused on CanMerge relationships.

    Structures
    ----------
    MAT Vertices          – all MAT vertices as a single grey point cloud (clickable).
    MAT Edges             – curve network.
    Input Mesh            – semi-transparent surface (if provided).
    Selected Vertex       – bright green; click any point on MAT Vertices.
    Selected Boundary Pts – red; boundary pts of the selected vertex.
    Mergeable Candidates  – gold; all vertices that pass CanMerge with selected.
    Candidate Boundary Pts– orange; combined boundary pts of all mergeable candidates.
    """
    try:
        import polyscope as ps
        import polyscope.imgui as psim
    except ImportError:
        print("Polyscope not installed.  Run:  pip install polyscope")
        return

    ps.init()
    ps.set_program_name("MAT Mergeability Viewer")
    ps.set_ground_plane_mode("none")
    ps.set_background_color((0.12, 0.12, 0.16))

    n_verts = len(vertex_centers)

    # ── input surface mesh ────────────────────────────────────────────────────
    if mesh_coords is not None and mesh_faces is not None and len(mesh_faces) > 0:
        sm = ps.register_surface_mesh("Input Mesh", mesh_coords, mesh_faces)
        sm.set_color((0.55, 0.70, 0.85))
        sm.set_transparency(0.35)
        sm.set_smooth_shade(True)

    # ── MAT edges ─────────────────────────────────────────────────────────────
    # if len(ma_edges) > 0:
        # cn = ps.register_curve_network("MAT Edges", vertex_centers, ma_edges)
        # cn.set_color((0.45, 0.65, 0.90))
        # cn.set_radius(0.0008, relative=True)

    # ── MAT mesh ─────────────────────────────────────────────────────────────
    if len(ma_faces) > 0:
        mm = ps.register_surface_mesh("MAT Mesh", ma_coords, ma_faces)
        mm.set_color((0.75, 0.75, 0.75))
        mm.set_transparency(0.25)
        mm.set_smooth_shade(True)


    # ── all MAT vertices (single cloud, clickable) ────────────────────────────
    pc_all = ps.register_point_cloud("MAT Vertices", vertex_centers)
    pc_all.set_color((0.75, 0.75, 0.75))
    pc_all.set_radius(0.0012, relative=True)

    # ── export path ───────────────────────────────────────────────────────────
    if ma_file:
        ma_dir  = os.path.dirname(os.path.abspath(ma_file))
        ma_stem = os.path.splitext(os.path.basename(ma_file))[0]
    else:
        ma_dir  = os.getcwd()
        ma_stem = "mat_export"
    obj_export_path = os.path.join(ma_dir, ma_stem + "_mat.obj")

    # ── selection state ───────────────────────────────────────────────────────
    state = {
        "sel_changed":  [False],
        "sel_idx":      [-1],
        "last_pick":    [None],
        "export_done":  [False],
        "export_msg":   [""],
        "n_mergeable":  [-1],
    }

    CLOUD_MAT      = "MAT Vertices"
    CLOUD_SEL      = "Selected Vertex"
    CLOUD_SEL_BP   = "Selected Boundary Pts"
    CLOUD_MERGE    = "Mergeable Candidates"
    CLOUD_MERGE_BP = "Candidate Boundary Pts"

    def _bplist_positions(local_idx):
        """Return (N,3) array of 3-D positions for the vertex_bplists boundary pts."""
        bp_ids = vertex_bplists[local_idx]
        return np.array(
            [sample_positions[int(sid)]
             for sid in bp_ids if int(sid) in sample_positions],
            dtype=np.float64)

    def _refresh_selection(idx):
        # selected vertex
        pc_sel = ps.register_point_cloud(CLOUD_SEL,
                                         np.array([vertex_centers[idx]]))
        pc_sel.set_color((0.0, 1.0, 0.53))
        pc_sel.set_radius(0.0025, relative=True)

        # boundary points of selected vertex
        bp_pos = _bplist_positions(idx)
        if len(bp_pos) > 0:
            pc_bp = ps.register_point_cloud(CLOUD_SEL_BP, bp_pos)
            pc_bp.set_color((1.0, 0.13, 0.27))
            pc_bp.set_radius(0.003, relative=True)
        else:
            try:
                ps.remove_point_cloud(CLOUD_SEL_BP)
            except Exception:
                pass

        # mergeable candidates + their boundary points
        if topo_local is not None and nmn_local is not None and voronoi_neighbors:
            merge_indices = [
                j for j in range(n_verts)
                if j != idx and can_merge(idx, j, topo_local, nmn_local, voronoi_neighbors)
            ]
            state["n_mergeable"][0] = len(merge_indices)

            if merge_indices:
                merge_pos = vertex_centers[np.array(merge_indices, dtype=np.int32)]
                pc_m = ps.register_point_cloud(CLOUD_MERGE, merge_pos)
                pc_m.set_color((1.0, 0.75, 0.0))
                pc_m.set_radius(0.0018, relative=True)

                # combined boundary points of all mergeable candidates
                all_bp = []
                for j in merge_indices:
                    pos = _bplist_positions(j)
                    if len(pos) > 0:
                        all_bp.append(pos)
                if all_bp:
                    pc_mbp = ps.register_point_cloud(
                        CLOUD_MERGE_BP, np.vstack(all_bp))
                    pc_mbp.set_color((1.0, 0.50, 0.0))
                    pc_mbp.set_radius(0.0022, relative=True)
                else:
                    try:
                        ps.remove_point_cloud(CLOUD_MERGE_BP)
                    except Exception:
                        pass
            else:
                for name in (CLOUD_MERGE, CLOUD_MERGE_BP):
                    try:
                        ps.remove_point_cloud(name)
                    except Exception:
                        pass
        else:
            state["n_mergeable"][0] = -1

    # ── ImGui callback ────────────────────────────────────────────────────────
    def callback():
        psim.PushItemWidth(340)
        psim.TextUnformatted("── MAT Mergeability Viewer ───────────────────")
        psim.Separator()
        psim.TextUnformatted(f"MAT vertices : {n_verts}")
        psim.TextUnformatted(f"MAT edges    : {len(ma_edges)}")
        psim.TextUnformatted(f"MAT faces    : {len(ma_faces)}")
        topo_loaded = topo_local is not None
        vn_loaded   = bool(voronoi_neighbors)
        psim.TextUnformatted(f"Topo sidecar : {'loaded' if topo_loaded else 'NOT FOUND'}")
        psim.TextUnformatted(f"Voronoi nbrs : {'loaded' if vn_loaded else 'NOT FOUND'}")

        # ── pick detection ────────────────────────────────────────────────────
        psim.Separator()
        psim.TextUnformatted("── Vertex Inspector ──────────────────────────")
        psim.TextUnformatted("Left-click any point on 'MAT Vertices' to inspect.")

        try:
            if ps.have_selection():
                pick = ps.get_selection()
                if pick is not None and pick.is_hit:
                    pick_name  = pick.structure_name
                    pick_local = int(pick.local_index)
                    pick_key   = (pick_name, pick_local)
                    if pick_key != state["last_pick"][0]:
                        state["last_pick"][0] = pick_key
                        if pick_name == CLOUD_MAT and 0 <= pick_local < n_verts:
                            state["sel_idx"][0]     = pick_local
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
            psim.TextUnformatted(f"  Index      : {idx}  (MAT id {int(vertex_indices[idx])})")
            psim.TextUnformatted(f"  Position   : ({cx:.4f}, {cy:.4f}, {cz:.4f})")
            psim.TextUnformatted(f"  Boundary pts: {len(vertex_bplists[idx])}  [red]")
            if topo_local is not None and topo_local[idx] is not None:
                fl = topo_local[idx]
                flags_str = " ".join(k for k in ('steep','sheet','seam','junction','boundary')
                                     if fl[k]) or "(none)"
                psim.TextUnformatted(f"  Topo flags : {flags_str}")
                nmn_cnt = len(nmn_local[idx]) if nmn_local is not None else 0
                psim.TextUnformatted(f"  nmn_bplist : {nmn_cnt} pts")
            nm = state["n_mergeable"][0]
            if nm >= 0:
                psim.TextUnformatted(f"  Mergeable  : {nm} candidate(s)  [gold]")
                if nm > 0:
                    psim.TextUnformatted(f"  Cand. BPs  : shown in orange")
            elif not topo_loaded or not vn_loaded:
                psim.TextUnformatted("  Mergeable  : (sidecars not loaded)")

        # ── export ────────────────────────────────────────────────────────────
        psim.Separator()
        psim.TextUnformatted("── Export ────────────────────────────────────")
        if psim.Button("Export MAT (.obj)"):
            try:
                export_mat_obj(obj_export_path, ma_coords, ma_edges, ma_faces)
                state["export_done"][0] = True
                state["export_msg"][0]  = f"OK  {os.path.basename(obj_export_path)}"
            except Exception as exc:
                state["export_done"][0] = True
                state["export_msg"][0]  = f"Export failed: {exc}"

        if state["export_done"][0]:
            psim.TextUnformatted(state["export_msg"][0])

        if mesh_coords is not None:
            psim.Separator()
            psim.TextUnformatted("Tip: toggle 'Input Mesh' transparency in its panel.")
        psim.PopItemWidth()

    ps.set_user_callback(callback)
    ps.show()


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    vertex_samples_file = sys.argv[1] if len(sys.argv) > 1 else VERTEX_SAMPLES_FILE
    ma_file             = sys.argv[2] if len(sys.argv) > 2 else MA_FILE
    mesh_file           = sys.argv[3] if len(sys.argv) > 3 else MESH_FILE

    if not vertex_samples_file or not ma_file:
        print("Usage:  python visualize_clusters_vornoi_neighbours.py "
              "<vertex_samples.txt> <file.ma> [mesh_file]")
        sys.exit(1)

    for f in [vertex_samples_file, ma_file]:
        if not os.path.exists(f):
            print(f"File not found: {f}")
            sys.exit(1)

    # ── load core files ───────────────────────────────────────────────────────
    print(f"Loading vertex samples: {vertex_samples_file}")
    vertex_indices, vertex_centers, vertex_bplists, sample_positions = \
        load_vertex_samples(vertex_samples_file)
    print(f"  {len(vertex_centers)} MAT vertices, {len(sample_positions)} sample positions")

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

    # ── load CanMerge sidecar files ───────────────────────────────────────────
    ma_stem_path = os.path.splitext(os.path.abspath(ma_file))[0]
    topo_path = ma_stem_path + "_mat_topo.txt"
    vn_path   = ma_stem_path + "_voronoi_neighbors.txt"

    print(f"Loading MAT topo sidecar : {topo_path}")
    topo_flags_dict, nmn_bplists_dict = load_mat_topo(topo_path)
    print(f"  {len(topo_flags_dict)} vertices loaded" if topo_flags_dict
          else "  (not found — CanMerge display disabled)")

    print(f"Loading Voronoi neighbors: {vn_path}")
    voronoi_neighbors = load_voronoi_neighbors(vn_path)
    print(f"  {len(voronoi_neighbors)} boundary points with neighbors" if voronoi_neighbors
          else "  (not found — CanMerge display disabled)")

    # Convert to local-index lists (vertex_indices[i] = actual MAT vertex idx)
    topo_local = nmn_local = None
    if topo_flags_dict:
        topo_local = [topo_flags_dict.get(int(vi)) for vi in vertex_indices]
        nmn_local  = [nmn_bplists_dict.get(int(vi), set()) for vi in vertex_indices]

    # ── visualise ─────────────────────────────────────────────────────────────
    visualize(vertex_indices, vertex_centers, vertex_bplists, sample_positions,
              ma_coords, ma_edges, ma_faces,
              mesh_coords=mesh_coords,
              mesh_faces=mesh_faces,
              ma_file=ma_file,
              topo_local=topo_local,
              nmn_local=nmn_local,
              voronoi_neighbors=voronoi_neighbors if voronoi_neighbors else None)
