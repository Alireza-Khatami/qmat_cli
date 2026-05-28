#!/usr/bin/env python3
"""Quantify the QuadErr<=eps clamp gap between us and MeshLab, per edge type.

For flat/coplanar regions MeshLab's QuadErr (=ScaleFactor*Apply at the collapse
position) falls below QuadricEpsilon and is clamped to eps*edgeLength (a tiny tie-
break ~1e-16). Ours often does not, because our Apply is numerically larger. With
ScaleFactor~1.25e7 and eps=1e-15 the clamp boundary sits at Apply~8e-23 -- a FP
knife-edge. This script counts, per edge type, how many edges each side "clamps"
(cost below 1e-15) and the median Apply(opt) we compute, to size the gap.

Usage:
    py clamp_gap.py [base_prefix]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qem_check as qc
from heap_compare import load_meshlab_heap, DEFAULT_BASE
from test_border_rule import build_quadrics

EPS = qc.EPS
CLAMP_COST = 1e-15  # below this, a cost is effectively the clamped tie-break value


def apply_at_opt(V, F, Qd, vfaces, edgef, SF, a, b):
    """Replicate ComputePosition + qq.Apply(opt) to expose the pre-clamp QuadErr."""
    qq = qc.Q(); qq.iadd(Qd[a]); qq.iadd(Qd[b])
    mid = (V[a] + V[b]) / 2
    if Qd[a].apply(mid) + Qd[b].apply(mid) > 2 * EPS:
        try:
            opt = np.linalg.solve(qq.A, -qq.b / 2)
        except np.linalg.LinAlgError:
            opt, _, _, _ = np.linalg.lstsq(qq.A, -qq.b / 2, rcond=None)
    else:
        opt = mid
    return qq.apply(opt), SF * qq.apply(opt)


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    V, F = qc.load_off(base + "_mat_initial.off")
    mlc = load_meshlab_heap(base + "_mat_initial.obj_initial_heap_state_.json")
    Qd, vfaces, edgef, SF, diag = build_quadrics(V, F, lambda c: c % 2 == 1)
    print(f"ScaleFactor={SF:.4e}  eps={EPS:.0e}  clamp@Apply<={EPS/SF:.2e}\n")

    types = {1: "border(1)", 2: "interior(2)", 3: "nonman(>=3)"}
    rows = {1: [], 2: [], 3: []}
    for (a, b), fs in edgef.items():
        m = mlc.get((a, b))
        if m is None or m > 1e307 or m <= 0:
            continue
        o = qc.cost(V, F, Qd, vfaces, edgef, SF, a, b, verbose=False)
        applyv, quaderr = apply_at_opt(V, F, Qd, vfaces, edgef, SF, a, b)
        rows[min(len(fs), 3)].append((o, m, applyv, quaderr))

    for c, name in types.items():
        r = rows[c]
        if not r:
            continue
        o = np.array([x[0] for x in r]); m = np.array([x[1] for x in r])
        quaderr = np.array([x[3] for x in r])
        ours_clamp = np.mean(o < CLAMP_COST) * 100
        ml_clamp = np.mean(m < CLAMP_COST) * 100
        we_clamp_q = np.mean(quaderr <= EPS) * 100
        print(f"{name:14s} n={len(r):5d}  ML clamps={ml_clamp:5.1f}%  ours clamp={ours_clamp:5.1f}%"
              f"  (our QuadErr<=eps: {we_clamp_q:4.1f}%)  median our QuadErr={np.median(quaderr):.2e}")


if __name__ == "__main__":
    main()
