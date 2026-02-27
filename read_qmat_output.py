"""
Reads QMAT output files and loads their contents into numpy arrays.

Usage:
    python read_qmat_output.py <meshname> [--visualize]

    where <meshname> is the base path used by QMAT (e.g. "output/bunny").
    The script expects:
        <meshname>_sampledpoints.txt
        <meshname>_vertex_samples.txt

    Pass --visualize (or -v) to open the interactive Polyscope viewer.
    In the viewer, click any orange MAT vertex to:
      - highlight it in green
      - show its associated sample points in red
      - display its info in the side panel

Variables produced
------------------
From _sampledpoints.txt:
    sample_ids      : (N,)   int    – sample point IDs
    sample_points   : (N, 3) float  – sample point coordinates (x, y, z)

From _vertex_samples.txt:
    vertex_indices  : (M,)   int    – Voronoi vertex indices
    vertex_centers  : (M, 3) float  – circumcenter coordinates (x, y, z)
    vertex_radii    : (M,)   float  – circumradii
    vertex_bplists  : (M,)   object – per-vertex int array of associated sample IDs
    vertex_bp_coords: (M,)   object – per-vertex float array of shape (K_i, 3)
"""

import os
import sys
import argparse
import numpy as np
import trimesh


# ─────────────────────────────────────────────────────────────────────────────
# I/O
# ─────────────────────────────────────────────────────────────────────────────
def check_non_manifold(faces):
    """
    Check for non-manifold edges and vertices in the MA mesh.
    Returns:
        dict containing:
            - 'non_manifold_edges': list of edges (as vertex pairs) shared by more than 2 faces
            - 'boundary_edges': list of edges shared by only 1 face
            - 'non_manifold_corner_vertices': list of vertices with non-manifold topology
            - 'edge_face_counts': dict mapping edges to their face counts
            - 'is_manifold': bool indicating if mesh is manifold
            - 'is_closed': bool indicating if mesh has no boundary
    """
    from collections import defaultdict
    
    if len(faces) == 0:
        return {
            'non_manifold_edges': [],
            'boundary_edges': [],
            'non_manifold_corner_vertices': [],
            'edge_face_counts': {},
            'is_manifold': True,
            'is_closed': True
        }
    # Build edge to face mapping
    edge_to_faces = defaultdict(list)
    for face_idx, face in enumerate(faces):
        num_verts = len(face)
        for i in range(num_verts):
            v1, v2 = face[i], face[(i + 1) % num_verts]
            # Use sorted tuple as edge key to ensure consistency
            edge = tuple(sorted([v1, v2]))
            edge_to_faces[edge].append(face_idx)
    # Identify non-manifold and boundary edges
    non_manifold_edges = []
    boundary_edges = []
    edge_face_counts = {}
    for edge, faces in edge_to_faces.items():
        count = len(faces)
        edge_face_counts[edge] = count
        if count > 2:
            non_manifold_edges.append(edge)
        elif count == 1:
            boundary_edges.append(edge)
    # A vertex is non-manifold if it is part of a non-manifold edge and also a boundary edge
    vertex_to_faces = defaultdict(list)

    boundary_vertices = set()
    for edge in boundary_edges:
        for v in edge:
            boundary_vertices.add(v.item())
    #set of non-manifold edges vertices
    non_manifold_edge_vertices = set()
    for edge in non_manifold_edges:
        for v in edge:
            non_manifold_edge_vertices.add(v.item())
    # non-manifold vertices = vertices on both boundary AND non-manifold edges
    non_manifold_corner_vertices = boundary_vertices.intersection(non_manifold_edge_vertices)
    non_manifold_corner_vertices = sorted(non_manifold_corner_vertices)
    
    # non_manifold_corner_vertices = coords[non_manifold_corner_vertices] if non_manifold_corner_vertices else np.empty((0, 3), dtype=np.float32)

    # Store results as instance attributes for visualization
    non_manifold_edges = non_manifold_edges
    boundary_edges = boundary_edges
    non_manifold_corner_vertices = non_manifold_corner_vertices
    return {
        'non_manifold_edges': non_manifold_edges,
        'boundary_edges': boundary_edges,
        'non_manifold_corner_vertices': non_manifold_corner_vertices,
        'edge_face_counts': edge_face_counts,
        'is_manifold': len(non_manifold_edges) == 0 and len(non_manifold_corner_vertices) == 0,
        'is_closed': len(boundary_edges) == 0
    }

def load_sampled_points(filepath):
    """
    Reads <meshname>_sampledpoints.txt.
    Returns:
        sample_ids    : (N,)   int
        sample_points : (N, 3) float
    """
    ids    = []
    coords = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) == 1:          # count line
                continue
            ids.append(int(parts[0]))
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    return np.array(ids, dtype=np.int32), np.array(coords, dtype=np.float64)

def load_vertex_samples(filepath):
    """
    Reads <meshname>_vertex_samples.txt.
    Returns:
        vertex_indices  : (M,)   int
        vertex_centers  : (M, 3) float
        vertex_radii    : (M,)   float
        vertex_bplists  : (M,)   object – int arrays
        vertex_bp_coords: (M,)   object – float (K_i, 3) arrays
    """
    v_indices   = []
    v_centers   = []
    v_radii     = []
    v_bplists   = []
    v_bp_coords = []

    cur_ids    = None
    cur_coords = None

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()

            if parts[0] == "v":
                if cur_ids is not None:
                    v_bplists.append(np.array(cur_ids,    dtype=np.int32))
                    v_bp_coords.append(np.array(cur_coords, dtype=np.float64))
                v_indices.append(int(parts[1]))
                v_centers.append([float(parts[2]), float(parts[3]), float(parts[4])])
                v_radii.append(float(parts[5]))
                cur_ids    = []
                cur_coords = []

            elif parts[0] == "s":
                cur_ids.append(int(parts[1]))
                cur_coords.append([float(parts[2]), float(parts[3]), float(parts[4])])

    if cur_ids is not None:
        v_bplists.append(np.array(cur_ids,    dtype=np.int32))
        v_bp_coords.append(np.array(cur_coords, dtype=np.float64))

    vertex_indices   = np.array(v_indices, dtype=np.int32)
    vertex_centers   = np.array(v_centers, dtype=np.float64)
    vertex_radii     = np.array(v_radii,   dtype=np.float64)
    vertex_bplists   = np.empty(len(v_bplists),   dtype=object)
    vertex_bp_coords = np.empty(len(v_bp_coords), dtype=object)
    for i, (bl, bc) in enumerate(zip(v_bplists, v_bp_coords)):
        vertex_bplists[i]   = bl
        vertex_bp_coords[i] = bc

    return vertex_indices, vertex_centers, vertex_radii, vertex_bplists, vertex_bp_coords

def load_ma(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load a .ma (medial axis) file.

    .ma format:
    - First line: n_verts n_edges n_faces
    - Vertex lines: v x y z radius
    - Edge lines: e v1 v2
    - Face lines: f v1 v2 v3

    Args:
        path: Path to the .ma file

    Returns:
        spheres: numpy array of shape (n_verts, 4) containing [x, y, z, radius]
        edges: numpy array of shape (n_edges, 2) containing vertex indices
        faces: numpy array of shape (n_faces, 3) containing vertex indices
        meta: QmatMeta with counts derived from file contents
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"File {path} does not exist.")
    if not path.endswith('.ma'):
        raise ValueError(f"File {path} is not a valid .ma file.")

    with open(path, 'r') as f:
        lines = f.readlines()

    # Parse header
    header = lines[0].split()
    n_verts, n_edges, n_faces = int(header[0]), int(header[1]), int(header[2])

    # Pre-allocate arrays
    spheres = np.zeros((n_verts, 4), dtype=np.float32)
    edges = np.zeros((n_edges, 2), dtype=np.int32)
    faces = np.zeros((n_faces, 3), dtype=np.int32)

    sphere_idx = 0
    edge_idx = 0
    face_idx = 0

    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if parts[0] == 'v':
            spheres[sphere_idx] = [float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])]
            sphere_idx += 1
        elif parts[0] == 'e':
            edges[edge_idx] = [int(parts[1]), int(parts[2])]
            edge_idx += 1
        elif parts[0] == 'f':
            faces[face_idx] = [int(parts[1]), int(parts[2]), int(parts[3])]
            face_idx += 1

    return spheres, edges, faces

# ─────────────────────────────────────────────────────────────────────────────
# Polyscope visualisation
# ─────────────────────────────────────────────────────────────────────────────
def visualize(sample_points, vertex_centers, vertex_radii,
              vertex_bplists, vertex_bp_coords, vertex_indices, ma_faces=None,
              topo_info_dict=None, trimesh_mesh=None):
    """
    Interactive Polyscope viewer.

    Point clouds
    ------------
    "Sample Points"        – all input sample points          (light blue)
    "MAT Vertices"         – all Voronoi / MAT vertices       (orange / yellow / magenta)
    "Selected MAT Vertex"  – the clicked MAT vertex           (green star)
    "Associated Samples"   – samples linked to that vertex    (red)

    Vertex colours (topology)
    -------------------------
    Blue     (0.2, 0.5,  1.0) – regular interior vertex
    Green    (0.1, 0.85, 0.2) – boundary vertex (on a boundary edge)
    red     (0.0, 0.85, 0.85) – non-manifold edge vertex (edge shared by >2 faces, not a corner)
    black   (0.0, 0.0,  0.0) – non-manifold corner vertex (on both a boundary and a non-manifold edge)

    How to interact
    ---------------
    Left-click any MAT vertex in the 3D view.
    The side panel will show its index, center, radius, topology type and sample IDs.
    """
    try:
        import polyscope as ps
        import polyscope.imgui as psim
    except ImportError:
        print("Polyscope is not installed.  Run:  pip install polyscope")
        return

    # ── precompute per-vertex topology colours ────────────────────────────────
    COLOR_NORMAL       = np.array([0.2,  0.5,  1.0 ], dtype=np.float32)  # blue
    COLOR_BOUNDARY     = np.array([0.1,  0.85, 0.2 ], dtype=np.float32)  # green
    COLOR_NM_EDGE      = np.array([1.0,  0.0,  0.0 ], dtype=np.float32)  # red   – non-manifold edge vertex
    COLOR_NM_CORNER    = np.array([0.0,  0.0,  0.0 ], dtype=np.float32)  # black  – non-manifold corner vertex

    if topo_info_dict is not None:
        boundary_vertex_set = set()
        for edge in topo_info_dict.get('boundary_edges', []):
            for v in edge:
                boundary_vertex_set.add(int(v))
        # vertices that sit on a non-manifold edge (edge shared by >2 faces)
        nm_edge_vertex_set = set()
        for edge in topo_info_dict.get('non_manifold_edges', []):
            for v in edge:
                nm_edge_vertex_set.add(int(v))
        # corners = boundary ∩ non-manifold-edge  (already a subset of nm_edge_vertex_set)
        nm_corner_vertex_set = set(
            int(v) for v in topo_info_dict.get('non_manifold_corner_vertices', []))
        # union used for the red highlight cloud
        non_manifold_vertex_set = nm_edge_vertex_set | nm_corner_vertex_set
    else:
        boundary_vertex_set     = set()
        nm_edge_vertex_set      = set()
        nm_corner_vertex_set    = set()
        non_manifold_vertex_set = set()

    # one colour row per MAT vertex (indexed the same as vertex_centers)
    vertex_colors = np.tile(COLOR_NORMAL, (len(vertex_centers), 1))
    for i, vi in enumerate(vertex_indices):
        vi_int = int(vi)
        if vi_int in nm_corner_vertex_set:
            vertex_colors[i] = COLOR_NM_CORNER      # black: corner
        elif vi_int in nm_edge_vertex_set:
            vertex_colors[i] = COLOR_NM_EDGE        # red: non-manifold edge (not a corner)
        elif vi_int in boundary_vertex_set:
            vertex_colors[i] = COLOR_BOUNDARY       # green: boundary only

    # separate cloud for non-manifold corner vertices (shown in red, no stride)
    nm_orig_indices = np.array(
        [i for i, vi in enumerate(vertex_indices) if int(vi) in non_manifold_vertex_set],
        dtype=np.int32)
    nm_positions = vertex_centers[nm_orig_indices] if len(nm_orig_indices) > 0 \
        else np.empty((0, 3), dtype=np.float64)

    ps.init()
    ps.set_program_name("QMAT Viewer")
    ps.set_ground_plane_mode("none")
    ps.set_background_color((0.12, 0.12, 0.16))

    max_stride = max(2, min(1000, max(len(sample_points), len(vertex_centers)) // 10))

    # ── mutable state shared with the callback ───────────────────────────────
    state = {
        "sel_idx":        -1,     # index into ORIGINAL vertex_centers arrays
        "sel_changed":    False,
        "stride":         [1],    # list so psim.SliderInt can mutate it
        "stride_changed": True,   # True on first frame to do the initial register
        # indices of currently displayed MAT vertices in the original array
        "mat_display_indices": np.arange(len(vertex_centers), dtype=np.int32),
        "save_path":      ["output_mat.off"],
        "save_status":    [""],
    }

    def _register_clouds(stride):
        """Re-register both point clouds using the given stride."""
        pc_s = ps.register_point_cloud("Sample Points", sample_points)
        pc_s.set_color((0.49, 0.78, 0.89))
        pc_s.set_radius(0.0015, relative=True)

        mat_disp_idx = np.arange(0, len(vertex_centers), stride, dtype=np.int32)
        state["mat_display_indices"] = mat_disp_idx
        pc_m = ps.register_point_cloud("MAT Vertices", vertex_centers[mat_disp_idx])
        pc_m.add_color_quantity("topology", vertex_colors[mat_disp_idx], enabled=True)
        pc_m.set_radius(0.001, relative=True)
        mesh_m = ps.register_surface_mesh("MAT mesh", vertex_centers, ma_faces)
        mesh_m.set_color((1.0, 0.55, 0.0))

        trimesh_mesh_ps = ps.register_surface_mesh("Original Mesh", trimesh_mesh.vertices, trimesh_mesh.faces,\
                                                    color=(0.8, 0.8, 0.8), edge_color=(0.5, 0.5, 0.5),\
                                                        transparency= .5)

        # non-manifold corner vertices — always shown at full resolution in red
        # if len(nm_positions) > 0:
        #     pc_nm = ps.register_point_cloud("Non-Manifold Vertices", nm_positions)
        #     pc_nm.set_color((1.0, 0.1, 0.1))
        #     pc_nm.set_radius(0.0018, relative=True)

    # ── per-frame callback (ImGui + pick handling) ───────────────────────────
    def callback():
        # --- stride slider ---
        psim.PushItemWidth(220)
        psim.Separator()
        psim.TextUnformatted("── Display Density ──")
        changed, new_stride = psim.SliderInt(
            "Show every N-th point", state["stride"][0], v_min=1, v_max=max_stride)
        if changed:
            state["stride"][0]    = new_stride
            state["stride_changed"] = True
            state["sel_idx"]      = -1   # clear selection when density changes

        n_shown_s = max(1, len(sample_points)  // state["stride"][0])
        n_shown_m = max(1, len(vertex_centers) // state["stride"][0])
        psim.TextUnformatted(
            f"Showing ~{n_shown_s} / {len(sample_points)} samples  |  "
            f"~{n_shown_m} / {len(vertex_centers)} MAT verts")

        # apply stride change
        if state["stride_changed"]:
            state["stride_changed"] = False
            _register_clouds(state["stride"][0])

        # --- pick handling ---
        sel = ps.get_selection()
        if sel is not None:
            struct_name = sel.structure_name
            elem_idx    = sel.local_index
            if struct_name == "MAT Vertices":
                # elem_idx is into the subsampled array — map to original index
                orig_idx = int(state["mat_display_indices"][elem_idx])
                if orig_idx != state["sel_idx"]:
                    state["sel_idx"]    = orig_idx
                    state["sel_changed"] = True
            elif struct_name == "Non-Manifold Vertices":
                # elem_idx is directly into nm_orig_indices
                orig_idx = int(nm_orig_indices[elem_idx])
                if orig_idx != state["sel_idx"]:
                    state["sel_idx"]    = orig_idx
                    state["sel_changed"] = True

        if state["sel_changed"]:
            state["sel_changed"] = False
            idx = state["sel_idx"]

            cx, cy, cz = vertex_centers[idx]
            pc_sel = ps.register_point_cloud(
                "Selected MAT Vertex", np.array([[cx, cy, cz]]))
            pc_sel.set_color((0.0, 1.0, 0.53))
            pc_sel.set_radius(0.002, relative=True)

            bp = vertex_bp_coords[idx]
            if len(bp) > 0:
                pc_bp = ps.register_point_cloud("Associated Samples", bp)
                pc_bp.set_color((1.0, 0.13, 0.27))
                pc_bp.set_radius(0.003, relative=True)
            else:
                try:
                    ps.remove_point_cloud("Associated Samples")
                except Exception:
                    pass

        # --- MAT vertex inspector ---
        psim.Separator()
        psim.TextUnformatted("── MAT Vertex Inspector ──")
        psim.Separator()

        if state["sel_idx"] < 0:
            psim.TextUnformatted("Click a MAT vertex to inspect it.")
        else:
            idx = state["sel_idx"]
            cx, cy, cz = vertex_centers[idx]
            vi_int = int(vertex_indices[idx])
            if vi_int in nm_corner_vertex_set:
                topo_label = "Non-manifold corner  [purple]"
            elif vi_int in nm_edge_vertex_set:
                topo_label = "Non-manifold edge    [cyan]"
            elif vi_int in boundary_vertex_set:
                topo_label = "Boundary             [green]"
            else:
                topo_label = "Regular interior     [blue]"
            psim.TextUnformatted(f"MAT vertex index : {vi_int}")
            psim.TextUnformatted(f"Topology         : {topo_label}")
            psim.TextUnformatted(f"Center X : {cx:.6f}")
            psim.TextUnformatted(f"Center Y : {cy:.6f}")
            psim.TextUnformatted(f"Center Z : {cz:.6f}")
            psim.TextUnformatted(f"Radius   : {vertex_radii[idx]:.6f}")
            psim.Separator()
            bp_ids = vertex_bplists[idx]
            psim.TextUnformatted(f"Associated samples : {len(bp_ids)}")
            psim.TextUnformatted(f"Sample IDs : {list(bp_ids)}")

        # --- Export ---
        psim.Separator()
        psim.TextUnformatted("── Export MAT ──")
        psim.Separator()
        changed, new_path = psim.InputText("Output path", state["save_path"][0])
        if changed:
            state["save_path"][0] = new_path
        if psim.Button("Save .off"):
            try:
                mat_mesh = trimesh.Trimesh(
                    vertices=vertex_centers,
                    faces=ma_faces,
                    process=False)
                mat_mesh.export(state["save_path"][0], file_type='off')
                state["save_status"][0] = f"Saved to {state['save_path'][0]}"
            except Exception as e:
                state["save_status"][0] = f"Error: {e}"
        if state["save_status"][0]:
            psim.TextUnformatted(state["save_status"][0])

        psim.PopItemWidth()

    ps.set_user_callback(callback)
    ps.show()

def load_mat_from_file( path):
        if path.endswith('.ma'):
            spheres, edges, faces = load_ma(path)
            coords = np.array(spheres[:, :3], dtype=np.float32)
            radii =  np.array(spheres[:, 3], dtype=np.float32)
            edges =  np.array(edges, dtype=np.int32)
            faces =  np.array(faces, dtype=np.int32)
        elif path.endswith('.off'):
            mesh = trimesh.load(path, file_type='off')
            coords = np.array(mesh.vertices, dtype=np.float32)
            radii = np.zeros(len(mesh.vertices), dtype=np.float32)  # No radius info in OFF
            edges = np.array(mesh.edges_unique, dtype=np.int32) if mesh.edges_unique is not None else np.empty((0, 2), dtype=np.int32)
            faces = np.array(mesh.faces, dtype=np.int32) if mesh.faces is not None else np.empty((0, 3), dtype=np.int32)
        return coords, faces
# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # parser = argparse.ArgumentParser(
    #     description="Load and optionally visualise QMAT output with Polyscope.")
    # parser.add_argument("meshname",
    #                     help="Base path/prefix used by QMAT (e.g. output/bunny)")
    # parser.add_argument("--visualize", "-v", action="store_true",
    #                     help="Open the interactive Polyscope 3D viewer after loading")
    # args = parser.parse_args()

    # ── load ──────────────────────────────────────────────────────────────────\
    meshname = rf"C:\Users\alirz\Projects\Graphics\QMAT_old working version  exe file\qmat_x64\qmat\cube_subdevided_fixed"
    mesh_file_path = rf"C:\Users\alirz\Projects\Graphics\Neural QMAT\Data\all of models\cube_subdevided_fixed.off"
    trimesh_mesh = trimesh.load(mesh_file_path, file_type='off')
    sample_ids, sample_points = load_sampled_points(
        meshname + "_sampledpoints.txt")
    print(f"Loaded {len(sample_ids)} sample points  "
          f"→  sample_ids{sample_ids.shape},  sample_points{sample_points.shape}")

    vertex_indices, vertex_centers, vertex_radii, vertex_bplists, vertex_bp_coords = \
        load_vertex_samples(meshname + "_vertex_samples.txt")
    print(f"Loaded {len(vertex_indices)} MAT vertices  "
          f"→  vertex_centers{vertex_centers.shape},  vertex_radii{vertex_radii.shape}")
    
    ma_vertices, ma_faces = load_mat_from_file(meshname + ".ma")

    # ── quick console summary of first vertex ─────────────────────────────────
    if len(vertex_indices) > 0:
        i = 0
        print(f"\nFirst MAT vertex [{vertex_indices[i]}]:")
        print(f"  center  = {vertex_centers[i]}")
        print(f"  radius  = {vertex_radii[i]:.6f}")
        print(f"  sample IDs    = {vertex_bplists[i]}")
        print(f"  sample coords =\n{vertex_bp_coords[i]}")

    # ── visualize ─────────────────────────────────────────────────────────────
    # if args.visualize:
    topo_info_dict = check_non_manifold(ma_faces)
    visualize(sample_points, vertex_centers, vertex_radii,
                vertex_bplists, vertex_bp_coords,\
                vertex_indices, ma_faces=ma_faces,\
                trimesh_mesh=trimesh_mesh,\
                topo_info_dict=topo_info_dict)
