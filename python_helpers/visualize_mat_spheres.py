"""
visualize_mat_spheres.py  –  Visualise mesh + MAT + sphere envelope.

Scans a root folder (default: cube_subdiv/) for numeric sub-folders.
Each sub-folder is treated as one simplification level.

Usage
-----
    python visualize_mat_spheres.py [root_folder]

    root_folder defaults to  cube_subdiv/  (relative to the script).

What is shown
-------------
    Input Mesh      – semi-transparent surface.
    MAT Mesh        – semi-transparent MAT faces.
    MAT Edges       – curve network (shown when no face mesh exists).
    MAT Vertices    – point cloud at sphere centres.
    Sphere Envelope – union of all MAT spheres as a single triangle mesh.

Side panel
----------
    Dropdown to switch simplification level.  All structures are replaced
    on selection change.  Export buttons write to the source folder.
"""

import sys
import os
import re
import glob
import numpy as np
import trimesh


# ── root folder containing numeric sub-folders ───────────────────────────────
_SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR     = os.path.join(_SCRIPT_DIR, "cube_subdiv")


# ─────────────────────────────────────────────────────────────────────────────
# Folder scanning
# ─────────────────────────────────────────────────────────────────────────────

def scan_levels(root):
    """
    Walk *root* and return a sorted list of Level dicts, one per numeric
    sub-folder.

    Each dict has:
        label    : str   – display label, e.g. "500"
        folder   : str   – absolute path to sub-folder
        ma_file  : str   – path to the simplified .ma  (___v_V___e_E___f_F pattern)
                           falls back to any *.ma if the pattern is not found
        mesh_file: str   – path to .off / .obj mesh file, searched in the
                           sub-folder first, then up to root, then siblings
    """
    levels = []
    if not os.path.isdir(root):
        return levels

    for entry in sorted(os.listdir(root)):
        folder = os.path.join(root, entry)
        if not os.path.isdir(folder):
            continue
        if not re.fullmatch(r'\d+', entry):
            continue   # skip non-numeric folders

        # ── find simplified .ma (prefer the one with ___v_…___e_…___f_… ) ──
        ma_file = ""
        simplified = glob.glob(os.path.join(folder, "*___v_*___e_*___f_*.ma"))
        if simplified:
            ma_file = simplified[0]
        else:
            fallback = glob.glob(os.path.join(folder, "*.ma"))
            if fallback:
                ma_file = fallback[0]

        # ── find mesh (.off / .obj) ──────────────────────────────────────────
        mesh_file = ""
        for ext in ("*.off", "*.obj"):
            hits = glob.glob(os.path.join(folder, ext))
            if hits:
                mesh_file = hits[0]
                break
        if not mesh_file:
            # search sibling folders (e.g. mesh only lives in folder 500)
            for sibling in sorted(os.listdir(root)):
                sibling_path = os.path.join(root, sibling)
                if not os.path.isdir(sibling_path):
                    continue
                for ext in ("*.off", "*.obj"):
                    hits = glob.glob(os.path.join(sibling_path, ext))
                    if hits:
                        mesh_file = hits[0]
                        break
                if mesh_file:
                    break
        if not mesh_file:
            # last resort: root itself
            for ext in ("*.off", "*.obj"):
                hits = glob.glob(os.path.join(root, ext))
                if hits:
                    mesh_file = hits[0]
                    break

        levels.append(dict(label=entry, folder=folder,
                           ma_file=ma_file, mesh_file=mesh_file))

    return levels


# ─────────────────────────────────────────────────────────────────────────────
# I/O
# ─────────────────────────────────────────────────────────────────────────────

def load_ma(path):
    """
    Parse a .ma file.

    Returns
    -------
    coords : (V, 3) float64
    radii  : (V,)   float64
    edges  : (E, 2) int32
    faces  : (F, 3) int32
    """
    coords, radii, edges, faces = [], [], [], []
    with open(path) as f:
        lines = f.readlines()
    for line in lines[1:]:
        parts = line.split()
        if not parts:
            continue
        if parts[0] == 'v':
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            radii.append(float(parts[4]))
        elif parts[0] == 'e':
            edges.append([int(parts[1]), int(parts[2])])
        elif parts[0] == 'f':
            faces.append([int(parts[1]), int(parts[2]), int(parts[3])])
    return (
        np.array(coords, dtype=np.float64),
        np.array(radii,  dtype=np.float64),
        np.array(edges,  dtype=np.int32) if edges else np.empty((0, 2), dtype=np.int32),
        np.array(faces,  dtype=np.int32) if faces else np.empty((0, 3), dtype=np.int32),
    )


def load_mesh(path):
    """Load any trimesh-supported surface mesh. Returns (verts, faces)."""
    mesh = trimesh.load(path, force='mesh')
    return mesh.vertices.astype(np.float64), mesh.faces.astype(np.int32)


# ─────────────────────────────────────────────────────────────────────────────
# MAT vertex topology classification
# ─────────────────────────────────────────────────────────────────────────────

def classify_mat_vertices(n_verts, edges, faces):
    """
    Classify every MAT vertex by its local topological role, derived purely
    from edge face-valence (number of incident faces per edge).

    Definitions (mirror of C++ DetermineTopology):
        sheet    – all incident edges are 2-manifold (face_valence == 2)
        boundary – at least one incident edge has face_valence == 1
        seam     – at least one incident edge has face_valence > 2
        junction – more than 2 seam edges meet at this vertex
        isolated – no incident edges at all

    Priority when multiple flags are set (highest wins for the colour label):
        junction > seam > boundary > sheet > isolated

    Returns
    -------
    labels : list of str, length n_verts
             one of 'junction' | 'seam' | 'boundary' | 'sheet' | 'isolated'
    counts : dict { label -> int }
    """
    from collections import defaultdict

    # edge key → face count
    edge_face_count = defaultdict(int)
    for v1, v2, v3 in faces:
        for a, b in ((v1, v2), (v2, v3), (v1, v3)):
            edge_face_count[(min(a, b), max(a, b))] += 1

    # per-vertex: sets of incident edge keys
    vert_edges = defaultdict(set)
    for v1, v2 in edges:
        key = (min(v1, v2), max(v1, v2))
        vert_edges[v1].add(key)
        vert_edges[v2].add(key)

    labels = []
    for vid in range(n_verts):
        inc = vert_edges[vid]
        if not inc:
            labels.append('isolated')
            continue

        n_seam     = sum(1 for k in inc if edge_face_count[k] > 2)
        n_boundary = sum(1 for k in inc if edge_face_count[k] == 1)

        if n_seam > 2:
            labels.append('junction')
        elif n_seam > 0:
            labels.append('seam')
        elif n_boundary > 0:
            labels.append('boundary')
        else:
            labels.append('sheet')

    counts = {lbl: labels.count(lbl)
              for lbl in ('sheet', 'boundary', 'seam', 'junction', 'isolated')}
    return labels, counts


# ─────────────────────────────────────────────────────────────────────────────
# Sphere envelope
# ─────────────────────────────────────────────────────────────────────────────

def unify_mesh_spheres(centers, radii, max_tris_total=400_000):
    """
    Build a single triangle mesh of one icosphere per MAT vertex.
    Subdivision level chosen adaptively to stay under max_tris_total triangles.
    """
    n = len(centers)
    if n == 0:
        return trimesh.Trimesh()

    for sd in (2, 1, 0):
        if 20 * (4 ** sd) * n <= max_tris_total:
            break

    print(f"  sphere subdivision: {sd}  "
          f"({20 * 4**sd} tris/sphere × {n} = {20 * 4**sd * n} total)")

    base  = trimesh.creation.icosphere(subdivisions=sd)
    parts = []
    for c, r in zip(centers, radii):
        r = float(r)
        if r <= 0.0:
            continue
        s = base.copy()
        s.apply_scale(r)
        s.apply_translation(c.astype(float))
        parts.append(s)

    return trimesh.util.concatenate(parts) if parts else trimesh.Trimesh()


# ─────────────────────────────────────────────────────────────────────────────
# Polyscope helpers
# ─────────────────────────────────────────────────────────────────────────────

# Names of all registered structures so we can remove them on reload
_STRUCT_NAMES = [
    "Input Mesh",
    "MAT Mesh",
    "MAT Edges",
    "MAT Vertices",
    "Sphere Envelope",
]


def _load_level(level):
    """Load and return all data for one level dict. Returns None on error."""
    if not level["ma_file"] or not os.path.exists(level["ma_file"]):
        print(f"  [warn] no .ma file for level {level['label']}")
        return None
    if not level["mesh_file"] or not os.path.exists(level["mesh_file"]):
        print(f"  [warn] no mesh file for level {level['label']}")
        return None

    print(f"Loading mesh : {level['mesh_file']}")
    mesh_verts, mesh_faces = load_mesh(level["mesh_file"])
    print(f"  {len(mesh_verts)} verts, {len(mesh_faces)} faces")

    print(f"Loading MAT  : {level['ma_file']}")
    mat_coords, mat_radii, mat_edges, mat_faces = load_ma(level["ma_file"])
    print(f"  {len(mat_coords)} verts  {len(mat_edges)} edges  {len(mat_faces)} faces")
    if len(mat_radii):
        print(f"  radius range: [{mat_radii.min():.4f}, {mat_radii.max():.4f}]")

    print("Classifying MAT vertex topology ...")
    topo_labels, topo_counts = classify_mat_vertices(
        len(mat_coords), mat_edges, mat_faces)
    for lbl, cnt in topo_counts.items():
        if cnt:
            print(f"  {lbl:10s}: {cnt}")

    print("Building sphere envelope ...")
    sphere_mesh = unify_mesh_spheres(mat_coords, mat_radii)
    print(f"  {len(sphere_mesh.vertices)} verts, {len(sphere_mesh.faces)} tris")

    return dict(
        mesh_verts=mesh_verts, mesh_faces=mesh_faces,
        mat_coords=mat_coords, mat_radii=mat_radii,
        mat_edges=mat_edges,   mat_faces=mat_faces,
        sphere_mesh=sphere_mesh,
        topo_labels=topo_labels, topo_counts=topo_counts,
        ma_stem=os.path.splitext(os.path.basename(level["ma_file"]))[0],
        ma_dir=level["folder"],
    )


def _register_structures(ps, data):
    """Register all Polyscope structures from a loaded data dict."""
    import polyscope as ps_mod   # may already be imported; this is safe
    ps_mod = ps

    sm = ps.register_surface_mesh("Input Mesh",
                                  data["mesh_verts"], data["mesh_faces"])
    sm.set_color((0.55, 0.70, 0.85))
    sm.set_transparency(0.40)
    sm.set_smooth_shade(True)

    if len(data["mat_faces"]) > 0:
        mm = ps.register_surface_mesh("MAT Mesh",
                                      data["mat_coords"], data["mat_faces"])
        mm.set_color((0.90, 0.60, 0.20))
        mm.set_transparency(0.30)
        mm.set_smooth_shade(False)

    if len(data["mat_edges"]) > 0:
        cn = ps.register_curve_network("MAT Edges",
                                       data["mat_coords"], data["mat_edges"])
        cn.set_color((1.0, 0.80, 0.30))
        cn.set_radius(0.0008, relative=True)
        cn.set_enabled(len(data["mat_faces"]) == 0)

    pc = ps.register_point_cloud("MAT Vertices", data["mat_coords"])
    pc.set_color((1.0, 1.0, 0.4))
    pc.set_radius(0.0012, relative=True)
    pc.set_enabled(False)   # topology clouds below are more informative

    # ── topology-classified point clouds ─────────────────────────────────────
    # colour scheme: sheet=blue, boundary=cyan, seam=orange, junction=magenta
    TOPO_STYLE = {
        'sheet':    ((0.30, 0.55, 1.00), 0.0010),
        'boundary': ((0.00, 0.90, 0.85), 0.0016),
        'seam':     ((1.00, 0.50, 0.10), 0.0018),
        'junction': ((0.90, 0.10, 0.90), 0.0022),
        'isolated': ((0.55, 0.55, 0.55), 0.0010),
    }
    labels = data["topo_labels"]
    coords = data["mat_coords"]
    for lbl, (col, rad) in TOPO_STYLE.items():
        idx = [i for i, l in enumerate(labels) if l == lbl]
        if not idx:
            continue
        name = f"Topo: {lbl}  (n={len(idx)})"
        tp = ps.register_point_cloud(name, coords[np.array(idx, dtype=np.int32)])
        tp.set_color(col)
        tp.set_radius(rad, relative=True)
        tp.set_enabled(lbl != 'isolated')   # hide isolated by default

    sph = data["sphere_mesh"]
    if len(sph.vertices) > 0:
        env = ps.register_surface_mesh(
            "Sphere Envelope",
            np.array(sph.vertices, dtype=np.float64),
            np.array(sph.faces,    dtype=np.int32),
            enabled=False
            )
        env.set_color((0.30, 0.85, 0.55))
        env.set_transparency(0.50)
        env.set_smooth_shade(True)


# ─────────────────────────────────────────────────────────────────────────────
# Polyscope visualisation
# ─────────────────────────────────────────────────────────────────────────────

def visualize(levels, initial_idx=0):
    try:
        import polyscope as ps
        import polyscope.imgui as psim
    except ImportError:
        print("Polyscope not installed.  Run:  pip install polyscope")
        return

    ps.init()
    ps.set_program_name("MAT Sphere Viewer")
    ps.set_ground_plane_mode("none")
    ps.set_background_color((0.10, 0.10, 0.14))

    # ── load initial level ────────────────────────────────────────────────────
    state = {
        "cur_idx":      [initial_idx],
        "pending_idx":  [initial_idx],   # set when dropdown changes
        "data":         [None],
        "export_done":  [False],
        "export_msg":   [""],
    }

    def _do_load(idx):
        ps.remove_all_structures()
        state["export_done"][0] = False
        state["export_msg"][0]  = ""
        data = _load_level(levels[idx])
        state["data"][0] = data
        if data:
            _register_structures(ps, data)
        state["cur_idx"][0] = idx

    _do_load(initial_idx)

    level_labels = [lv["label"] for lv in levels]

    # ── ImGui callback ────────────────────────────────────────────────────────
    def callback():
        psim.PushItemWidth(220)
        psim.TextUnformatted("── MAT Sphere Viewer ─────────────────────────")
        psim.Separator()

        # ── level selector ────────────────────────────────────────────────────
        psim.TextUnformatted("Simplification level (sample count):")
        changed, new_idx = psim.Combo("##level", state["cur_idx"][0], level_labels)
        if changed:
            state["pending_idx"][0] = new_idx

        # reload happens outside the combo to avoid mid-frame structure changes
        if state["pending_idx"][0] != state["cur_idx"][0]:
            _do_load(state["pending_idx"][0])

        psim.Separator()

        data = state["data"][0]
        if data is None:
            psim.TextUnformatted("  [no data loaded]")
            psim.PopItemWidth()
            return

        psim.TextUnformatted(f"Folder       : {levels[state['cur_idx'][0]]['folder']}")
        psim.TextUnformatted(f"MA file      : {os.path.basename(levels[state['cur_idx'][0]]['ma_file'])}")
        psim.Separator()
        psim.TextUnformatted(f"Input mesh   : {len(data['mesh_verts'])} verts  {len(data['mesh_faces'])} faces")
        psim.TextUnformatted(f"MAT vertices : {len(data['mat_coords'])}")
        psim.TextUnformatted(f"MAT edges    : {len(data['mat_edges'])}")
        psim.TextUnformatted(f"MAT faces    : {len(data['mat_faces'])}")
        sph = data["sphere_mesh"]
        psim.TextUnformatted(f"Sphere tris  : {len(sph.faces) if len(sph.vertices) > 0 else 0}")
        psim.Separator()
        psim.TextUnformatted("Vertex topology:")
        tc = data["topo_counts"]
        psim.TextUnformatted(f"  sheet    : {tc.get('sheet',    0):6d}  (blue)")
        psim.TextUnformatted(f"  boundary : {tc.get('boundary', 0):6d}  (cyan)")
        psim.TextUnformatted(f"  seam     : {tc.get('seam',     0):6d}  (orange)")
        psim.TextUnformatted(f"  junction : {tc.get('junction', 0):6d}  (magenta)")
        psim.TextUnformatted(f"  isolated : {tc.get('isolated', 0):6d}  (grey)")
        psim.Separator()
        psim.TextUnformatted("Layers (toggle in left panel):")
        psim.TextUnformatted("  Input Mesh      – blue,    transparent")
        psim.TextUnformatted("  MAT Mesh        – orange,  transparent")
        psim.TextUnformatted("  MAT Edges       – yellow curve")
        psim.TextUnformatted("  Topo: sheet     – blue points")
        psim.TextUnformatted("  Topo: boundary  – cyan points")
        psim.TextUnformatted("  Topo: seam      – orange points")
        psim.TextUnformatted("  Topo: junction  – magenta points")
        psim.TextUnformatted("  Sphere Envelope – green,   transparent")
        psim.Separator()

        # ── export buttons ────────────────────────────────────────────────────
        ma_dir  = data["ma_dir"]
        ma_stem = data["ma_stem"]

        if psim.Button("Export MAT (.obj)"):
            path = os.path.join(ma_dir, ma_stem + "_mat.obj")
            try:
                with open(path, 'w') as f:
                    f.write("# MAT exported by visualize_mat_spheres.py\n")
                    f.write(f"# {len(data['mat_coords'])} verts  "
                            f"{len(data['mat_edges'])} edges  "
                            f"{len(data['mat_faces'])} faces\n\n")
                    for x, y, z in data["mat_coords"]:
                        f.write(f"v {x:.10g} {y:.10g} {z:.10g}\n")
                    if len(data["mat_edges"]):
                        f.write("\n")
                        for v1, v2 in data["mat_edges"]:
                            f.write(f"l {v1+1} {v2+1}\n")
                    if len(data["mat_faces"]):
                        f.write("\n")
                        for v1, v2, v3 in data["mat_faces"]:
                            f.write(f"f {v1+1} {v2+1} {v3+1}\n")
                state["export_done"][0] = True
                state["export_msg"][0]  = f"OK  {os.path.basename(path)}"
            except Exception as exc:
                state["export_done"][0] = True
                state["export_msg"][0]  = f"Failed: {exc}"

        if len(data["mat_faces"]) > 0:
            psim.SameLine()
            if psim.Button("Export MAT (.off)"):
                path = os.path.join(ma_dir, ma_stem + "_mat.off")
                try:
                    with open(path, 'w') as f:
                        f.write("OFF\n")
                        f.write(f"{len(data['mat_coords'])} {len(data['mat_faces'])} 0\n")
                        for x, y, z in data["mat_coords"]:
                            f.write(f"{x:.10g} {y:.10g} {z:.10g}\n")
                        for v1, v2, v3 in data["mat_faces"]:
                            f.write(f"3 {v1} {v2} {v3}\n")
                    state["export_done"][0] = True
                    state["export_msg"][0]  = f"OK  {os.path.basename(path)}"
                except Exception as exc:
                    state["export_done"][0] = True
                    state["export_msg"][0]  = f"Failed: {exc}"

        if len(sph.vertices) > 0:
            psim.SameLine()
            if psim.Button("Export Spheres (.obj)"):
                path = os.path.join(ma_dir, ma_stem + "_spheres.obj")
                try:
                    sph.export(path)
                    state["export_done"][0] = True
                    state["export_msg"][0]  = f"OK  {os.path.basename(path)}"
                except Exception as exc:
                    state["export_done"][0] = True
                    state["export_msg"][0]  = f"Failed: {exc}"

        if state["export_done"][0]:
            psim.TextUnformatted(state["export_msg"][0])

        psim.PopItemWidth()

    ps.set_user_callback(callback)
    ps.show()


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) > 1 else ROOT_DIR

    if not os.path.isdir(root):
        print(f"Root folder not found: {root}")
        sys.exit(1)

    print(f"Scanning: {root}")
    levels = scan_levels(root)

    if not levels:
        print("No numeric sub-folders found.")
        sys.exit(1)

    for lv in levels:
        status = "OK" if lv["ma_file"] else "NO .ma"
        mesh_s = os.path.basename(lv["mesh_file"]) if lv["mesh_file"] else "NO MESH"
        print(f"  [{lv['label']}]  {os.path.basename(lv['ma_file']) if lv['ma_file'] else '—'}  |  {mesh_s}  [{status}]")

    visualize(levels, initial_idx=0)
