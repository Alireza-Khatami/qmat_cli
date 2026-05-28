#!/usr/bin/env python3
"""Dissect non-manifold edge costs (>2 incident faces) to find why ours is far
cheaper than MeshLab's. For each sampled non-manifold edge, prints the incident
face count, the per-endpoint face counts, and the full ComputePriority breakdown
(gate, Apply at optimum, QuadErr, newQual, MinCos, cost) alongside MeshLab's cost.

The same non-manifold mesh is fed to MeshLab, so the >2-face edges are identical
on both sides; any cost gap is in the quadric/cost evaluation, not the topology.

Usage:
    py nonman_dissect.py [base_prefix] [num_edges]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qem_check as qc
from heap_compare import load_meshlab_heap, DEFAULT_BASE


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    V, F = qc.load_off(base + "_mat_initial.off")
    Qd, vfaces, edgef, SF, diag = qc.build(V, F)
    mlc = load_meshlab_heap(base + "_mat_initial.obj_initial_heap_state_.json")
    print(f"ScaleFactor={SF:.4e} diag={diag:.4f}\n")

    shown = 0
    for (a, b), fs in edgef.items():
        if len(fs) <= 2:
            continue
        m = mlc.get((a, b))
        if m is None or m > 1e307:
            continue
        o = qc.cost(V, F, Qd, vfaces, edgef, SF, a, b, verbose=False)
        print(f"edge {a}-{b}: incident_faces={len(fs)}  "
              f"vfaces(a)={len(vfaces[a])} vfaces(b)={len(vfaces[b])}")
        print(f"   ours={o:.4e}  meshlab={m:.4e}  ratio(ours/ml)={o / m:.2e}")
        qc.cost(V, F, Qd, vfaces, edgef, SF, a, b, verbose=True)
        shown += 1
        if shown >= nmax:
            break


if __name__ == "__main__":
    main()
