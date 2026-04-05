"""
qmat_batch_viewer.py  –  Interactive QMAT batch-output viewer

Scans the folder structure produced by batch_qmat.sh:

  <output_root>/
    <mesh_name>/
      <N>/
        <mesh_name>_v_<N>.ma                        raw medial axis
        <mesh_name>_v_<N>___v_V___e_E___f_F.ma      simplified MAT (SlabMesh::Export)
        <mesh_name>_v_<N>_sampledpoints.txt
        <mesh_name>_v_<N>_vertex_samples.txt

Usage:
    python qmat_batch_viewer.py <output_root> [--mesh_folder <input_folder>]

    --mesh_folder  optional path to the original .obj/.off files so the
                   source mesh can be overlaid in the viewer.

Polyscope side-panel UI
------------------------
  Mesh dropdown            – switch between mesh names found in <output_root>
  Simplification dropdown  – switch between available N values for that mesh
  MAT Vertex Inspector     – click any orange MAT vertex to inspect it

Structures shown
-----------------
  Sample Points            light-blue point cloud
  MAT Vertices             point cloud coloured by topology
                             blue   = regular interior
                             green  = boundary
                             red    = non-manifold edge
                             black  = non-manifold corner
  Raw MA Mesh              orange semi-transparent surface (raw .ma faces)
  Simplified MAT Spheres   green semi-transparent sphere mesh (unify_mesh_spheres)
  Simplified MAT Mesh      red-orange semi-transparent surface (simplified faces)
  Simplified MAT Edges     magenta curve network (simplified edges)
  Original Mesh            grey semi-transparent (only when --mesh_folder given)
"""

import os
import sys
import glob
import argparse
from collections import defaultdict

import numpy as np
import trimesh


# ─────────────────────────────────────────────────────────────────────────────
# Topology helpers
# ─────────────────────────────────────────────────────────────────────────────

def check_non_manifold(faces):
    """
    Classify edges and vertices by topology.

    Returns a dict:
        non_manifold_edges          list of (v1,v2) shared by >2 faces
        boundary_edges              list of (v1,v2) shared by exactly 1 face
        non_manifold_corner_vertices sorted list of vertex IDs on both
        edge_face_counts            {edge: int}
        is_manifold, is_closed      bool
    """
    if len(faces) == 0:
        return dict(non_manifold_edges=[], boundary_edges=[],
                    non_manifold_corner_vertices=[], edge_face_counts={},
                    is_manifold=True, is_closed=True)

    edge_to_faces = defaultdict(list)
    for fi, face in enumerate(faces):
        n = len(face)
        for i in range(n):
            e = tuple(sorted([int(face[i]), int(face[(i + 1) % n])]))
            edge_to_faces[e].append(fi)

    nm_edges, bnd_edges, efc = [], [], {}
    for e, fs in edge_to_faces.items():
        c = len(fs);  efc[e] = c
        if   c > 2:  nm_edges.append(e)
        elif c == 1: bnd_edges.append(e)

    bnd_verts  = {v for e in bnd_edges  for v in e}
    nm_e_verts = {v for e in nm_edges   for v in e}
    nm_c_verts = sorted(bnd_verts & nm_e_verts)

    return dict(non_manifold_edges=nm_edges, boundary_edges=bnd_edges,
                non_manifold_corner_vertices=nm_c_verts,
                edge_face_counts=efc,
                is_manifold=len(nm_edges) == 0 and len(nm_c_verts) == 0,
                is_closed=len(bnd_edges) == 0)


# ─────────────────────────────────────────────────────────────────────────────
# File I/O
# ─────────────────────────────────────────────────────────────────────────────

def load_sampled_points(filepath):
    """
    Reads <prefix>_sampledpoints.txt.
    Returns  sample_ids (N,) int,  sample_points (N,3) float.
    """
    ids, coords = [], []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split()
            if len(parts) == 1: continue          # count-only line
            ids.append(int(parts[0]))
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return np.array(ids, dtype=np.int32), np.array(coords, dtype=np.float64)


def load_vertex_samples(filepath):
    """
    Reads <prefix>_vertex_samples.txt.
    Returns  vertex_indices (M,),  vertex_centers (M,3),  vertex_radii (M,),
             vertex_bplists (M,) object,  vertex_bp_coords (M,) object.
    """
    v_idx, v_ctr, v_rad, v_bl, v_bc = [], [], [], [], []
    cur_ids = cur_coords = None

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split()
            if parts[0] == 'v':
                if cur_ids is not None:
                    v_bl.append(np.array(cur_ids,    dtype=np.int32))
                    v_bc.append(np.array(cur_coords, dtype=np.float64))
                v_idx.append(int(parts[1]))
                v_ctr.append([float(parts[2]), float(parts[3]), float(parts[4])])
                v_rad.append(float(parts[5]))
                cur_ids, cur_coords = [], []
            elif parts[0] == 's':
                cur_ids.append(int(parts[1]))
                cur_coords.append([float(parts[2]), float(parts[3]), float(parts[4])])

    if cur_ids is not None:
        v_bl.append(np.array(cur_ids,    dtype=np.int32))
        v_bc.append(np.array(cur_coords, dtype=np.float64))

    vertex_bplists   = np.empty(len(v_bl), dtype=object)
    vertex_bp_coords = np.empty(len(v_bc), dtype=object)
    for i, (bl, bc) in enumerate(zip(v_bl, v_bc)):
        vertex_bplists[i]   = bl
        vertex_bp_coords[i] = bc

    return (np.array(v_idx, dtype=np.int32),
            np.array(v_ctr, dtype=np.float64),
            np.array(v_rad, dtype=np.float64),
            vertex_bplists, vertex_bp_coords)


def load_ma(path):
    """
    Load a .ma file (raw or simplified MAT).

    Format:
        <n_verts> <n_edges> <n_faces>
        v  x  y  z  radius
        e  v1  v2
        f  v1  v2  v3          (quads are fan-triangulated)

    Returns
    -------
    spheres : (V, 4) float32   [x, y, z, radius]
    edges   : (E, 2) int32
    faces   : (F, 3) int32
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    with open(path) as f:
        lines = f.readlines()

    h = lines[0].split()
    n_v, n_e, n_f = int(h[0]), int(h[1]), int(h[2])

    spheres    = np.zeros((n_v, 4), dtype=np.float32)
    edges      = np.zeros((n_e, 2), dtype=np.int32)
    faces_list = []
    vi = ei = 0

    for line in lines[1:]:
        line = line.strip()
        if not line: continue
        p = line.split()
        if p[0] == 'v' and vi < n_v:
            spheres[vi] = [float(p[1]), float(p[2]), float(p[3]), float(p[4])]
            vi += 1
        elif p[0] == 'e' and ei < n_e:
            edges[ei] = [int(p[1]), int(p[2])]
            ei += 1
        elif p[0] == 'f':
            verts = [int(x) for x in p[1:]]
            if len(verts) == 3:
                faces_list.append(verts)
            elif len(verts) >= 4:                       # fan triangulation
                for k in range(1, len(verts) - 1):
                    faces_list.append([verts[0], verts[k], verts[k + 1]])

    faces = (np.array(faces_list, dtype=np.int32)
             if faces_list else np.empty((0, 3), dtype=np.int32))
    return spheres, edges, faces


# ─────────────────────────────────────────────────────────────────────────────
# Folder scanning
# ─────────────────────────────────────────────────────────────────────────────

def scan_output_folder(root_dir):
    """
    Walk <root_dir>/<mesh_name>/<N>/ and collect all available meshes and
    simplification levels.

    Returns
    -------
    mesh_names : sorted list of str
    simp_map   : {mesh_name: sorted list of int}
    """
    mesh_names, simp_map = [], {}
    if not os.path.isdir(root_dir):
        return mesh_names, simp_map

    for entry in sorted(os.scandir(root_dir), key=lambda e: e.name):
        if not entry.is_dir(): continue
        ns = []
        for sub in sorted(os.scandir(entry.path), key=lambda e: e.name):
            if sub.is_dir():
                try:    ns.append(int(sub.name))
                except: pass
        if ns:
            mesh_names.append(entry.name)
            simp_map[entry.name] = sorted(ns)

    return mesh_names, simp_map


def find_simplified_ma(folder, prefix):
    """
    Find the SlabMesh::Export output: <prefix>___v_V___e_E___f_F.ma
    Returns the path string or None.
    """
    hits = sorted(glob.glob(os.path.join(folder, f"{prefix}___v_*___e_*___f_*.ma")))
    return hits[-1] if hits else None


# ─────────────────────────────────────────────────────────────────────────────
# Sphere-mesh construction
# ─────────────────────────────────────────────────────────────────────────────

def unify_mesh_spheres(centers, radii, max_tris_total=400_000):
    """
    Build a single triangle mesh consisting of one icosphere per MAT vertex,
    each scaled to its circumradius.  Subdivision level is chosen adaptively
    so the total triangle count stays below max_tris_total.

    Parameters
    ----------
    centers        : (N, 3) float
    radii          : (N,)   float
    max_tris_total : int    budget for total triangles across all spheres

    Returns
    -------
    trimesh.Trimesh  (may be empty if all radii ≤ 0)
    """
    n = len(centers)
    if n == 0:
        return trimesh.Trimesh()

    # icosphere(sd) has  20 * 4**sd  faces
    for sd in (2, 1, 0):
        if 20 * (4 ** sd) * n <= max_tris_total:
            break

    base = trimesh.creation.icosphere(subdivisions=sd)
    parts = []
    for c, r in zip(centers, radii):
        r = float(r)
        if r <= 0.0: continue
        s = base.copy()
        s.apply_scale(r)
        s.apply_translation(c.astype(float))
        parts.append(s)

    return trimesh.util.concatenate(parts) if parts else trimesh.Trimesh()


# ─────────────────────────────────────────────────────────────────────────────
# Per-selection data loader
# ─────────────────────────────────────────────────────────────────────────────

def load_data_for_selection(root_dir, mesh_name, N, mesh_folder=None):
    """
    Load all available data for (mesh_name, N).

    Returns a dict with keys (each optional – present only when the file exists):
        sample_ids, sample_points
        vertex_indices, vertex_centers, vertex_radii,
          vertex_bplists, vertex_bp_coords
        raw_spheres, raw_edges, raw_faces
        simp_spheres, simp_edges, simp_faces, simp_ma_path
        original_mesh  (trimesh.Trimesh, only when mesh_folder given)
    """
    folder = os.path.join(root_dir, mesh_name, str(N))
    prefix = f"{mesh_name}_v_{N}"
    d = {}

    # raw MA
    raw_path = os.path.join(folder, f"{prefix}.ma")
    if os.path.exists(raw_path):
        try:
            d['raw_spheres'], d['raw_edges'], d['raw_faces'] = load_ma(raw_path)
        except Exception as e:
            print(f"  [warn] raw MA: {e}")

    # simplified MAT
    simp_path = find_simplified_ma(folder, prefix)
    if simp_path:
        try:
            d['simp_spheres'], d['simp_edges'], d['simp_faces'] = load_ma(simp_path)
            d['simp_ma_path'] = simp_path
        except Exception as e:
            print(f"  [warn] simplified MA: {e}")

    # sample points
    sp = os.path.join(folder, f"{prefix}_sampledpoints.txt")
    if os.path.exists(sp):
        try:
            d['sample_ids'], d['sample_points'] = load_sampled_points(sp)
        except Exception as e:
            print(f"  [warn] sample points: {e}")

    # vertex samples
    vs = os.path.join(folder, f"{prefix}_vertex_samples.txt")
    if os.path.exists(vs):
        try:
            (d['vertex_indices'], d['vertex_centers'], d['vertex_radii'],
             d['vertex_bplists'], d['vertex_bp_coords']) = load_vertex_samples(vs)
        except Exception as e:
            print(f"  [warn] vertex samples: {e}")

    # original mesh
    if mesh_folder:
        for ext in ('.off', '.obj', '.OFF', '.OBJ'):
            mp = os.path.join(mesh_folder, mesh_name + ext)
            if os.path.exists(mp):
                try:
                    d['original_mesh'] = trimesh.load(mp, process=False, force='mesh')
                except Exception as e:
                    print(f"  [warn] original mesh: {e}")
                break

    return d


# ─────────────────────────────────────────────────────────────────────────────
# Polyscope visualisation
# ─────────────────────────────────────────────────────────────────────────────

_ALL_STRUCTURES = [
    "Sample Points", "MAT Vertices",
    "Raw MA Mesh",
    "Simplified MAT Spheres", "Simplified MAT Mesh", "Simplified MAT Edges",
    "Simplified MAT Vertices",
    "Original Mesh",
    "Selected MAT Vertex", "Associated Samples",
]

def _remove_all_structures():
    """Remove every known structure, checking existence first to avoid
    polyscope logging C++ 'not registered' exceptions."""
    import polyscope as ps
    _PC  = (ps.has_point_cloud,   ps.remove_point_cloud)
    _SM  = (ps.has_surface_mesh,  ps.remove_surface_mesh)
    _CN  = (ps.has_curve_network, ps.remove_curve_network)
    for name in _ALL_STRUCTURES:
        for has_fn, rm_fn in (_PC, _SM, _CN):
            try:
                if has_fn(name):
                    rm_fn(name)
            except Exception:
                pass


def visualize(root_dir, mesh_names, simp_map, mesh_folder=None):
    try:
        import polyscope as ps
        import polyscope.imgui as psim
    except ImportError:
        print("Polyscope not installed.  Run:  pip install polyscope")
        return

    # colours
    C_SAMPLE    = (0.49, 0.78, 0.89)
    C_RAW_MESH  = (1.00, 0.55, 0.00)
    C_SIMP_SPH  = (0.30, 0.85, 0.40)
    C_SIMP_MESH = (0.90, 0.30, 0.20)
    C_SIMP_EDGE = (0.90, 0.20, 0.80)
    C_ORIG      = (0.80, 0.80, 0.80)
    C_REGULAR   = np.array([0.20, 0.50, 1.00], dtype=np.float32)
    C_BOUNDARY  = np.array([0.10, 0.85, 0.20], dtype=np.float32)
    C_NM_EDGE   = np.array([1.00, 0.00, 0.00], dtype=np.float32)
    C_NM_COR    = np.array([0.00, 0.00, 0.00], dtype=np.float32)

    MESH_NAME = "cube_subdevided_fixed"

    # mutable state (single-element lists so inner functions can mutate them)
    state = {
        "simp_idx":       [0],
        "needs_reload":   [True],
        "reload_frame":   [0],
        "data":           [{}],
        "sel_idx":        [-1],
        "sel_changed":    [False],
        "bnd_verts":      [set()],
        "nm_edge_verts":  [set()],
        "nm_corner_verts":[set()],
    }

    def current_mesh():
        return MESH_NAME

    def current_simps():
        return simp_map[MESH_NAME]

    def current_N():
        simps = current_simps()
        return simps[min(state["simp_idx"][0], len(simps) - 1)]

    # ── reload ────────────────────────────────────────────────────────────────
    def reload():
        mesh_name = current_mesh()
        N         = current_N()
        print(f"\nLoading  {mesh_name}  /  N={N} …")

        _remove_all_structures()
        state["sel_idx"][0]     = -1
        state["sel_changed"][0] = False
        state["bnd_verts"][0]       = set()
        state["nm_edge_verts"][0]   = set()
        state["nm_corner_verts"][0] = set()

        d = load_data_for_selection(root_dir, mesh_name, N, mesh_folder)
        state["data"][0] = d

        # sample points
        if 'sample_points' in d:
            pc = ps.register_point_cloud("Sample Points", d['sample_points'], enabled=False)
            pc.set_color(C_SAMPLE)
            pc.set_radius(0.0015, relative=True)

        # MAT vertices coloured by topology
        if 'vertex_centers' in d:
            vc = d['vertex_centers']
            vi = d['vertex_indices']

            if 'raw_faces' in d and len(d['raw_faces']) > 0:
                topo = check_non_manifold(d['raw_faces'])
                bnd  = {v for e in topo['boundary_edges']    for v in e}
                nm_e = {v for e in topo['non_manifold_edges'] for v in e}
                nm_c = set(topo['non_manifold_corner_vertices'])
            else:
                bnd = nm_e = nm_c = set()

            state["bnd_verts"][0]       = bnd
            state["nm_edge_verts"][0]   = nm_e
            state["nm_corner_verts"][0] = nm_c

            cols = np.tile(C_REGULAR, (len(vc), 1))
            for i, v in enumerate(vi):
                v = int(v)
                if   v in nm_c: cols[i] = C_NM_COR
                elif v in nm_e: cols[i] = C_NM_EDGE
                elif v in bnd:  cols[i] = C_BOUNDARY

            pc = ps.register_point_cloud("MAT Vertices", vc)
            pc.add_color_quantity("topology", cols, enabled=True)
            pc.set_radius(0.001, relative=True)

        # raw MA mesh
        if 'raw_spheres' in d and 'raw_faces' in d and len(d['raw_faces']) > 0:
            m = ps.register_surface_mesh(
                "Raw MA Mesh", d['raw_spheres'][:, :3], d['raw_faces'])
            m.set_color(C_RAW_MESH);  m.set_transparency(0.45)

        # simplified MAT
        if 'simp_spheres' in d:
            sc = d['simp_spheres'][:, :3]
            sr = d['simp_spheres'][:, 3]

            # ── topology colours for the simplified MAT ───────────────────────
            simp_cols = np.tile(C_REGULAR, (len(sc), 1))
            if 'simp_faces' in d and len(d['simp_faces']) > 0:
                stopo = check_non_manifold(d['simp_faces'])
                s_bnd  = {v for e in stopo['boundary_edges']     for v in e}
                s_nm_e = {v for e in stopo['non_manifold_edges']  for v in e}
                s_nm_c = set(stopo['non_manifold_corner_vertices'])
                for vi in range(len(sc)):
                    if   vi in s_nm_c: simp_cols[vi] = C_NM_COR
                    elif vi in s_nm_e: simp_cols[vi] = C_NM_EDGE
                    elif vi in s_bnd:  simp_cols[vi] = C_BOUNDARY

            # unify_mesh_spheres – one icosphere per MAT vertex
            print(f"  Building sphere mesh for {len(sc)} vertices …")
            sph_mesh = unify_mesh_spheres(sc, sr)
            if len(sph_mesh.vertices) > 0:
                sm = ps.register_surface_mesh(
                    "Simplified MAT Spheres",
                    np.array(sph_mesh.vertices, dtype=np.float64),
                    np.array(sph_mesh.faces,    dtype=np.int32))
                sm.set_color(C_SIMP_SPH);  sm.set_transparency(0.35)

            if 'simp_faces' in d and len(d['simp_faces']) > 0:
                m2 = ps.register_surface_mesh(
                    "Simplified MAT Mesh", sc, d['simp_faces'])
                m2.add_color_quantity("topology", simp_cols, defined_on='vertices', enabled=True)
                m2.set_transparency(0.30)

            # point cloud coloured by topology so individual vertices are visible
            pc_simp = ps.register_point_cloud("Simplified MAT Vertices", sc)
            pc_simp.add_color_quantity("topology", simp_cols, enabled=True)
            pc_simp.set_radius(0.002, relative=True)

            if 'simp_edges' in d and len(d['simp_edges']) > 0:
                cn = ps.register_curve_network(
                    "Simplified MAT Edges", sc, d['simp_edges'])
                cn.set_color(C_SIMP_EDGE)
                cn.set_radius(0.0008, relative=True)

        # original mesh
        if 'original_mesh' in d:
            om = d['original_mesh']
            m3 = ps.register_surface_mesh(
                "Original Mesh",
                np.array(om.vertices, dtype=np.float64),
                np.array(om.faces,    dtype=np.int32),
                color=C_ORIG, edge_color=(0.5, 0.5, 0.5))
            m3.set_transparency(0.50)

        print("  Done.")

    # ── per-frame callback ────────────────────────────────────────────────────
    def callback():
        psim.PushItemWidth(280)

        # ── simplification dropdown ────────────────────────────────────────────
        psim.TextUnformatted(f"Mesh: {MESH_NAME}")
        psim.Separator()

        simps     = current_simps()
        simp_strs = [str(s) for s in simps]
        changed_s, new_s_idx = psim.Combo("Simplification (N)", state["simp_idx"][0], simp_strs)
        if changed_s and new_s_idx != state["simp_idx"][0]:
            state["simp_idx"][0]     = new_s_idx
            state["needs_reload"][0] = True
            state["reload_frame"][0] = 2

        # Defer the actual reload by `reload_frame` frames so polyscope can
        # finish the current draw before we clear and re-register structures.
        if state["needs_reload"][0]:
            if state["reload_frame"][0] > 0:
                state["reload_frame"][0] -= 1
            else:
                state["needs_reload"][0] = False
                reload()

        # ── pick handling ─────────────────────────────────────────────────────
        d = state["data"][0]
        sel = ps.get_selection()
        if sel is not None and sel.structure_name == "MAT Vertices":
            orig_idx = int(sel.local_index)
            if orig_idx != state["sel_idx"][0]:
                state["sel_idx"][0]     = orig_idx
                state["sel_changed"][0] = True

        if state["sel_changed"][0]:
            state["sel_changed"][0] = False
            idx = state["sel_idx"][0]
            if 'vertex_centers' in d:
                cx, cy, cz = d['vertex_centers'][idx]
                pc_s = ps.register_point_cloud(
                    "Selected MAT Vertex", np.array([[cx, cy, cz]]))
                pc_s.set_color((0.0, 1.0, 0.53))
                pc_s.set_radius(0.0025, relative=True)

                bp = d['vertex_bp_coords'][idx]
                if len(bp) > 0:
                    pc_bp = ps.register_point_cloud("Associated Samples", bp)
                    pc_bp.set_color((1.0, 0.13, 0.27))
                    pc_bp.set_radius(0.003, relative=True)
                else:
                    try: ps.remove_point_cloud("Associated Samples")
                    except: pass

        # ── info panel ────────────────────────────────────────────────────────
        psim.Separator()
        psim.TextUnformatted("── Simplified MAT Info ──────────")
        if 'simp_ma_path' in d:
            psim.TextUnformatted(f"File : {os.path.basename(d['simp_ma_path'])}")
        if 'simp_spheres' in d:
            psim.TextUnformatted(
                f"V={len(d['simp_spheres'])}  "
                f"E={len(d.get('simp_edges', []))}  "
                f"F={len(d.get('simp_faces', []))}")

        # ── MAT vertex inspector ──────────────────────────────────────────────
        psim.Separator()
        psim.TextUnformatted("── MAT Vertex Inspector ─────────")
        psim.Separator()

        if state["sel_idx"][0] < 0 or 'vertex_centers' not in d:
            psim.TextUnformatted("Click a MAT vertex (orange) to inspect it.")
        else:
            idx    = state["sel_idx"][0]
            vi_int = int(d['vertex_indices'][idx])
            bnd    = state["bnd_verts"][0]
            nm_e   = state["nm_edge_verts"][0]
            nm_c   = state["nm_corner_verts"][0]

            if   vi_int in nm_c: topo = "Non-manifold corner  [black]"
            elif vi_int in nm_e: topo = "Non-manifold edge    [red]"
            elif vi_int in bnd:  topo = "Boundary             [green]"
            else:                topo = "Regular interior     [blue]"

            cx, cy, cz = d['vertex_centers'][idx]
            psim.TextUnformatted(f"Index    : {vi_int}")
            psim.TextUnformatted(f"Topology : {topo}")
            psim.TextUnformatted(f"X : {cx:.6f}")
            psim.TextUnformatted(f"Y : {cy:.6f}")
            psim.TextUnformatted(f"Z : {cz:.6f}")
            psim.TextUnformatted(f"R : {d['vertex_radii'][idx]:.6f}")
            psim.Separator()
            bp_ids = d['vertex_bplists'][idx]
            psim.TextUnformatted(f"Samples : {len(bp_ids)}")
            preview = list(bp_ids[:20])
            psim.TextUnformatted(
                f"IDs : {preview}" + (" …" if len(bp_ids) > 20 else ""))

        psim.PopItemWidth()

    # ── launch ────────────────────────────────────────────────────────────────
    ps.init()
    ps.set_program_name("QMAT Batch Viewer")
    ps.set_ground_plane_mode("none")
    ps.set_background_color((0.12, 0.12, 0.16))
    ps.set_user_callback(callback)
    ps.show()


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # ap = argparse.ArgumentParser(
    #     description="Interactive viewer for batch_qmat.sh output folders.")
    # ap.add_argument("output_root",
    #                 help="Root folder produced by batch_qmat.sh  "
    #                      "(contains <mesh_name>/<N>/ sub-folders)")
    # ap.add_argument("--mesh_folder", default=None,
    #                 help="Folder with original .obj/.off files for overlay "
    #                      "(optional)")
    # args = ap.parse_args()
    output_root = "C:/Users/alirz/Projects/Graphics/QMAT_old working version  exe file/qmat_x64/qmat/batch_qmat_output"
    mesh_folder = "C:/Users/alirz/Projects/Graphics/QMAT_old working version  exe file/qmat_x64/qmat/mesh"
    print(f"Scanning: {output_root}")
    mesh_names, simp_map = scan_output_folder(output_root)

    if not mesh_names:
        print("No mesh folders found.  Make sure batch_qmat.sh has been run first.")
        sys.exit(1)

    print(f"Found {len(mesh_names)} mesh(es):")
    for m in mesh_names:
        print(f"  {m}  →  {simp_map[m]}")
    print()

    visualize(output_root, mesh_names, simp_map, mesh_folder)
