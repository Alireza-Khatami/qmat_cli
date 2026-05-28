#!/usr/bin/env python3
"""Compare our full initial-heap costs against MeshLab's, at one or more ScaleFactors.

Reports, per ScaleFactor: median/mean log10(ours/MeshLab), fraction within 2x/10x,
and the overlap of the cheapest-200 edge sets. Shows whether a global ScaleFactor
change reconciles the heaps (it does not — see md_files/qem analysis).

Usage:
    py heap_compare.py [base_prefix]
base_prefix defaults to the 01_00040057... test mesh; the script reads
<base>_mat_initial.off and <base>_mat_initial.obj_initial_heap_state_.json.
"""
import os
import re
import io
import json
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qem_check as qc

DEFAULT_BASE = ("output/01_00040057_f8f78dbd17414efda75bc437_trimesh_000/"
                "01_00040057_f8f78dbd17414efda75bc437_trimesh_000")


def load_meshlab_heap(path):
    """Per-edge min cost from a MeshLab initial_heap_state_.json (sanitizing inf/nan)."""
    raw = open(path).read()
    raw = re.sub(r'(?<=:)-?inf', '1e308', raw)
    raw = re.sub(r'(?<=:)-?nan', '1e308', raw)
    ml = json.load(io.StringIO(raw))
    mlc = {}
    for e in ml["entries"]:
        k = (min(e["v0_id"], e["v1_id"]), max(e["v0_id"], e["v1_id"]))
        if k not in mlc or e["cost"] < mlc[k]:
            mlc[k] = e["cost"]
    return mlc


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    V, F = qc.load_off(base + "_mat_initial.off")
    Qd, vfaces, edgef, SF, diag = qc.build(V, F)
    mlc = load_meshlab_heap(base + "_mat_initial.obj_initial_heap_state_.json")
    edges = sorted(edgef.keys())

    for sf, label in [(SF, f"ours (diag={diag:.3f})"), (3.8e4, "reduced 3.8e4")]:
        oc = {(a, b): qc.cost(V, F, Qd, vfaces, edgef, sf, a, b, verbose=False)
              for (a, b) in edges}
        common = [k for k in oc if k in mlc and 0 < mlc[k] < 1e307 and oc[k] > 0]
        lr = np.array([np.log10(oc[k] / mlc[k]) for k in common])
        ours_cheap = set(k for k, _ in sorted(oc.items(), key=lambda x: x[1])[:200])
        ml_cheap = set(k for k, _ in sorted(mlc.items(), key=lambda x: x[1])[:200])
        print(f"\nSF={sf:.3e} [{label}]")
        print(f"  log10(ours/ml): median={np.median(lr):+.2f} mean={lr.mean():+.2f}")
        print(f"  within 2x: {np.mean(np.abs(lr) < np.log10(2)) * 100:.1f}%   "
              f"within 10x: {np.mean(np.abs(lr) < 1) * 100:.1f}%")
        print(f"  cheapest-200 overlap with MeshLab: {len(ours_cheap & ml_cheap)}/200")


if __name__ == "__main__":
    main()
