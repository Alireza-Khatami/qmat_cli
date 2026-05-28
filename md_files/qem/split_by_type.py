#!/usr/bin/env python3
"""Split the ours-vs-MeshLab cost disagreement by edge type (border / interior /
non-manifold) to localize where our QEM cost diverges.

Finding: borders match well (~92% within 10x); interior median ~1x but wide spread;
non-manifold (>2 faces) edges are ~2400x too cheap on our side (0% within 10x).

Usage:
    py split_by_type.py [base_prefix]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qem_check as qc
from heap_compare import load_meshlab_heap, DEFAULT_BASE


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    V, F = qc.load_off(base + "_mat_initial.off")
    Qd, vfaces, edgef, SF, diag = qc.build(V, F)
    mlc = load_meshlab_heap(base + "_mat_initial.obj_initial_heap_state_.json")

    border, interior, nonman = [], [], []
    for (a, b), fs in edgef.items():
        m = mlc.get((a, b))
        o = qc.cost(V, F, Qd, vfaces, edgef, SF, a, b, verbose=False)
        if m is None or m > 1e307 or m <= 0 or o <= 0:
            continue
        lr = np.log10(o / m)
        (border if len(fs) == 1 else (interior if len(fs) == 2 else nonman)).append(lr)

    for name, arr in [("border(1 face)", border),
                      ("interior(2 face)", interior),
                      ("nonman(>2 face)", nonman)]:
        arr = np.array(arr)
        print(f"{name:18s} n={len(arr):5d}  "
              f"median log10(ours/ml)={np.median(arr):+.2f} ({10 ** np.median(arr):.2g}x)  "
              f"within10x={np.mean(np.abs(arr) < 1) * 100:.0f}%")


if __name__ == "__main__":
    main()
