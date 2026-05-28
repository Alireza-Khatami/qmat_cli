#!/usr/bin/env python3
"""Localize the interior (2-face) edge cost overshoot under the odd-count border
rule: split interior edges by whether an endpoint touches an odd-count edge
(border or non-manifold junction). If only junction-adjacent edges overshoot, the
cause is border-quadric accumulation at junctions, not the interior edges per se.

Usage:
    py analyze_interior.py [base_prefix]
"""
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qem_check as qc
from heap_compare import load_meshlab_heap, DEFAULT_BASE
from test_border_rule import build_quadrics


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    V, F = qc.load_off(base + "_mat_initial.off")
    mlc = load_meshlab_heap(base + "_mat_initial.obj_initial_heap_state_.json")
    Qd, vfaces, edgef, SF, diag = build_quadrics(V, F, lambda c: c % 2 == 1)

    # vertices incident to any odd-count edge (border or non-manifold)
    odd_vert = set()
    for (a, b), fs in edgef.items():
        if len(fs) % 2 == 1:
            odd_vert.add(a); odd_vert.add(b)

    pure, junctional = [], []
    worst = []  # (ratio, edge, info)
    for (a, b), fs in edgef.items():
        if len(fs) != 2:
            continue
        m = mlc.get((a, b))
        o = qc.cost(V, F, Qd, vfaces, edgef, SF, a, b, verbose=False)
        if m is None or m > 1e307 or m <= 0 or o <= 0:
            continue
        lr = np.log10(o / m)
        touches = (a in odd_vert) or (b in odd_vert)
        (junctional if touches else pure).append(lr)
        if touches:
            worst.append((o / m, (a, b)))

    for name, arr in [("pure interior", pure), ("junction-adjacent", junctional)]:
        arr = np.array(arr)
        if len(arr) == 0:
            print(f"{name:20s} n=0"); continue
        print(f"{name:20s} n={len(arr):5d}  median ours/ml={10 ** np.median(arr):.3g}x  "
              f"within10x={np.mean(np.abs(arr) < 1) * 100:.0f}%")

    print("\nworst-overshoot junction-adjacent interior edges:")
    for ratio, (a, b) in sorted(worst, reverse=True)[:5]:
        # how many odd edges touch each endpoint, and border-quadric count
        def odd_edges_at(v):
            return [(min(v, w), max(v, w)) for w in
                    set(F[fi][k] for fi in vfaces[v] for k in range(3) if F[fi][k] != v)
                    if len(edgef[(min(v, w), max(v, w))]) % 2 == 1]
        oa, ob = odd_edges_at(a), odd_edges_at(b)
        print(f"  {a}-{b}: ratio={ratio:.2e}  odd-edges@a={len(oa)} odd-edges@b={len(ob)}")


if __name__ == "__main__":
    main()
