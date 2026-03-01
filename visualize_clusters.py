"""
visualize_clusters.py  –  Visualise a _boundary_clusters.txt file produced by
                          ThreeDimensionalShape::ClusterBoundaryPoints()

Usage
-----
    python visualize_clusters.py <path_to_*_boundary_clusters.txt>

    # or edit the CLUSTER_FILE constant below and run without arguments.

File format
-----------
    # comment lines (ignored)
    <num_clusters>  <num_points>
    <cluster_id>  <vertex_index>  <x>  <y>  <z>  <nx>  <ny>  <nz>
    ...

Polyscope view
--------------
    • One point-cloud per cluster, each coloured differently.
    • Normals shown as a vector quantity on each cloud.
"""

import sys
import os
import numpy as np

# ── hard-coded fallback (edit if running without a CLI argument) ──────────────
CLUSTER_FILE = ""   # e.g. "batch_qmat_output/cube_subdevided_fixed/50/cube_subdevided_fixed_v_50_boundary_clusters.txt"

# ─────────────────────────────────────────────────────────────────────────────
# I/O
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

    # convert to numpy
    for cid, d in clusters.items():
        d['positions'] = np.array(d['positions'], dtype=np.float64)
        d['normals']   = np.array(d['normals'],   dtype=np.float64)
        d['indices']   = np.array(d['indices'],   dtype=np.int32)

    meta = dict(num_clusters=num_clusters, num_points=num_points,
                filepath=filepath)
    return clusters, meta


# ─────────────────────────────────────────────────────────────────────────────
# Colour palette  (cycles for large cluster counts)
# ─────────────────────────────────────────────────────────────────────────────

def _palette(n):
    """Return n visually distinct RGB colours as (n,3) float32."""
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
    # tile + jitter for extra clusters
    rng  = np.random.default_rng(42)
    tiled = np.tile(base, (n // len(base) + 1, 1))[:n]
    jitter = rng.uniform(-0.12, 0.12, tiled.shape).astype(np.float32)
    return np.clip(tiled + jitter, 0.0, 1.0)


# ─────────────────────────────────────────────────────────────────────────────
# Polyscope visualisation
# ─────────────────────────────────────────────────────────────────────────────

def visualize(clusters, meta):
    try:
        import polyscope as ps
        import polyscope.imgui as psim
    except ImportError:
        print("Polyscope not installed.  Run:  pip install polyscope")
        return

    ps.init()
    ps.set_program_name("Boundary Cluster Viewer")
    ps.set_ground_plane_mode("none")
    ps.set_background_color((0.12, 0.12, 0.16))

    colours = _palette(meta['num_clusters'])

    # Register one point cloud per cluster
    for cid, d in sorted(clusters.items()):
        pts  = d['positions']
        nrms = d['normals']
        col  = tuple(colours[cid % len(colours)].tolist())
        name = f"Cluster {cid}  (n={len(pts)})"

        pc = ps.register_point_cloud(name, pts)
        pc.set_color(col)
        pc.set_radius(0.0018, relative=True)
        pc.add_vector_quantity("normal", nrms, enabled=False,
                               vectortype='ambient', length=0.02,
                               color=col)

    # ── imgui info panel ──────────────────────────────────────────────────────
    fname = os.path.basename(meta['filepath'])

    def callback():
        psim.PushItemWidth(300)
        psim.TextUnformatted("── Boundary Cluster Info ────────────────")
        psim.Separator()
        psim.TextUnformatted(f"File     : {fname}")
        psim.TextUnformatted(f"Clusters : {meta['num_clusters']}")
        psim.TextUnformatted(f"Points   : {meta['num_points']}")
        psim.Separator()
        psim.TextUnformatted("Cluster sizes:")
        for cid, d in sorted(clusters.items()):
            col = colours[cid % len(colours)]
            psim.TextUnformatted(f"  [{cid:3d}]  {len(d['positions']):6d} pts")
        psim.PopItemWidth()

    ps.set_user_callback(callback)
    ps.show()


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else CLUSTER_FILE
    if not path:
        print("Usage:  python visualize_clusters.py <*_boundary_clusters.txt>")
        sys.exit(1)
    if not os.path.exists(path):
        print(f"File not found: {path}")
        sys.exit(1)

    print(f"Loading: {path}")
    clusters, meta = load_clusters(path)
    print(f"  {meta['num_clusters']} clusters, {meta['num_points']} points")
    for cid, d in sorted(clusters.items()):
        print(f"  cluster {cid:3d} : {len(d['positions'])} pts")
    print()

    visualize(clusters, meta)
