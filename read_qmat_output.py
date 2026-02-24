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

import sys
import argparse
import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# I/O
# ─────────────────────────────────────────────────────────────────────────────

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


# ─────────────────────────────────────────────────────────────────────────────
# Polyscope visualisation
# ─────────────────────────────────────────────────────────────────────────────

def visualize(sample_points, vertex_centers, vertex_radii,
              vertex_bplists, vertex_bp_coords, vertex_indices):
    """
    Interactive Polyscope viewer.

    Point clouds
    ------------
    "Sample Points"        – all input sample points          (light blue)
    "MAT Vertices"         – all Voronoi / MAT vertices       (orange)
    "Selected MAT Vertex"  – the clicked MAT vertex           (green star)
    "Associated Samples"   – samples linked to that vertex    (red)

    How to interact
    ---------------
    Left-click any orange MAT vertex in the 3D view.
    The side panel will show its index, center, radius and sample IDs.
    """
    try:
        import polyscope as ps
        import polyscope.imgui as psim
    except ImportError:
        print("Polyscope is not installed.  Run:  pip install polyscope")
        return

    ps.init()
    ps.set_program_name("QMAT Viewer")
    ps.set_ground_plane_mode("none")
    ps.set_background_color((0.12, 0.12, 0.16))

    # ── static point clouds ──────────────────────────────────────────────────
    pc_samples = ps.register_point_cloud("Sample Points", sample_points)
    pc_samples.set_color((0.49, 0.78, 0.89))    # light blue
    pc_samples.set_radius(0.0015, relative=True)

    pc_mat = ps.register_point_cloud("MAT Vertices", vertex_centers)
    pc_mat.set_color((1.0, 0.55, 0.0))          # orange
    pc_mat.set_radius(0.003, relative=True)

    # ── mutable state shared with the callback ───────────────────────────────
    state = {
        "sel_idx":      -1,   # index into vertex_centers arrays
        "sel_changed":  False,
    }

    # ── per-frame callback (ImGui + pick handling) ───────────────────────────
    def callback():
        # --- check for a new selection ---
        sel = ps.get_selection()
        if sel is not None:
            struct_name, elem_idx = sel.structure_name , sel.local_index
            if struct_name == "MAT Vertices" and elem_idx != state["sel_idx"]:
                state["sel_idx"]     = elem_idx
                state["sel_changed"] = True

        if state["sel_changed"]:
            state["sel_changed"] = False
            idx = state["sel_idx"]

            # green marker on the selected MAT vertex
            cx, cy, cz = vertex_centers[idx]
            pc_sel = ps.register_point_cloud(
                "Selected MAT Vertex", np.array([[cx, cy, cz]]))
            pc_sel.set_color((0.0, 1.0, 0.53))   # green
            pc_sel.set_radius(0.007, relative=True)

            # red highlight on its associated sample points
            bp = vertex_bp_coords[idx]
            if len(bp) > 0:
                pc_bp = ps.register_point_cloud("Associated Samples", bp)
                pc_bp.set_color((1.0, 0.13, 0.27))  # red
                pc_bp.set_radius(0.005, relative=True)
            else:
                # remove stale cloud if this vertex has no samples
                try:
                    ps.remove_point_cloud("Associated Samples")
                except Exception:
                    pass

        # --- ImGui info panel in the Polyscope side panel ---
        psim.PushItemWidth(220)
        psim.Separator()
        psim.TextUnformatted("── MAT Vertex Inspector ──")
        psim.Separator()

        if state["sel_idx"] < 0:
            psim.TextUnformatted("Click an orange MAT vertex to inspect it.")
        else:
            idx = state["sel_idx"]
            cx, cy, cz = vertex_centers[idx]
            psim.TextUnformatted(f"MAT vertex index : {int(vertex_indices[idx])}")
            psim.TextUnformatted(f"Center X : {cx:.6f}")
            psim.TextUnformatted(f"Center Y : {cy:.6f}")
            psim.TextUnformatted(f"Center Z : {cz:.6f}")
            psim.TextUnformatted(f"Radius   : {vertex_radii[idx]:.6f}")
            psim.Separator()
            bp_ids = vertex_bplists[idx]
            psim.TextUnformatted(f"Associated samples : {len(bp_ids)}")
            psim.TextUnformatted(f"Sample IDs : {list(bp_ids)}")

        psim.PopItemWidth()

    ps.set_user_callback(callback)
    ps.show()


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
    sample_ids, sample_points = load_sampled_points(
        meshname + "_sampledpoints.txt")
    print(f"Loaded {len(sample_ids)} sample points  "
          f"→  sample_ids{sample_ids.shape},  sample_points{sample_points.shape}")

    vertex_indices, vertex_centers, vertex_radii, vertex_bplists, vertex_bp_coords = \
        load_vertex_samples(meshname + "_vertex_samples.txt")
    print(f"Loaded {len(vertex_indices)} MAT vertices  "
          f"→  vertex_centers{vertex_centers.shape},  vertex_radii{vertex_radii.shape}")

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
    visualize(sample_points, vertex_centers, vertex_radii,
                  vertex_bplists, vertex_bp_coords, vertex_indices)
