#!/usr/bin/env python3
"""Test a midpoint-fallback guard for the rank-deficient quadric solve.

For near-flat fans the summed quadric is near-singular; the SOLVE branch then
returns a wild, far-from-edge position with a huge Apply -> huge cost, while
MeshLab keeps such edges clamped (Apply(mid) stays below the gate). This tests
falling back to the midpoint when the solved optimum lands farther than
guard_factor * edge_length from the midpoint, and reports per-type agreement
with MeshLab for several guard factors (inf = no guard = current behavior).

Usage:
    py test_solve_guard.py [base_prefix]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qem_check as qc
from heap_compare import load_meshlab_heap, DEFAULT_BASE
from test_border_rule import build_quadrics

EPS = qc.EPS
QTHR = qc.QTHR
COS_THR = qc.COS_THR


def cost_guarded(V, F, Qd, vfaces, edgef, SF, a, b, guard):
    qq = qc.Q(); qq.iadd(Qd[a]); qq.iadd(Qd[b])
    mid = (V[a] + V[b]) / 2
    if Qd[a].apply(mid) + Qd[b].apply(mid) > 2 * EPS:
        try:
            opt = np.linalg.solve(qq.A, -qq.b / 2)
        except np.linalg.LinAlgError:
            opt, _, _, _ = np.linalg.lstsq(qq.A, -qq.b / 2, rcond=None)
        el = np.linalg.norm(V[a] - V[b])
        if el > 0 and np.linalg.norm(opt - mid) > guard * el:   # wild solve -> midpoint
            opt = mid
    else:
        opt = mid
    newQual = 1e308; mincos = 1e308
    for (va, partner) in ((a, b), (b, a)):
        for fi in vfaces[va]:
            f = F[fi]
            if partner in f:
                continue
            pp = [opt if f[k] == va else V[f[k]] for k in range(3)]
            newQual = min(newQual, qc.qualityface(*pp))
            nb = qc.trinorm(V[f[0]], V[f[1]], V[f[2]])
            mincos = min(mincos, float(qc.trinorm(*pp) @ nb))
    QuadErr = SF * qq.apply(opt)
    if newQual > QTHR:
        newQual = QTHR
    if mincos > COS_THR:
        mincos = COS_THR
    mincos = abs((mincos + 1) / 2)
    QuadErr = max(QuadErr, EPS)
    if QuadErr <= EPS:
        QuadErr *= np.linalg.norm(V[a] - V[b])
    return QuadErr / (newQual * mincos)


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    V, F = qc.load_off(base + "_mat_initial.off")
    mlc = load_meshlab_heap(base + "_mat_initial.obj_initial_heap_state_.json")
    Qd, vfaces, edgef, SF, diag = build_quadrics(V, F, lambda c: c % 2 == 1)

    for guard in [float('inf'), 10.0, 4.0, 2.0, 1.0]:
        buckets = {1: [], 2: [], 3: []}
        oc = {}
        for (a, b), fs in edgef.items():
            m = mlc.get((a, b))
            o = cost_guarded(V, F, Qd, vfaces, edgef, SF, a, b, guard)
            oc[(a, b)] = o
            if m is None or m > 1e307 or m <= 0 or o <= 0:
                continue
            buckets[min(len(fs), 3)].append(np.log10(o / m))
        ours_cheap = set(k for k, _ in sorted(oc.items(), key=lambda x: x[1])[:200])
        ml_cheap = set(k for k, _ in sorted(mlc.items(), key=lambda x: x[1])[:200])
        allr = np.concatenate([np.array(buckets[c]) for c in (1, 2, 3)])
        gtxt = "none" if guard == float('inf') else f"{guard:g}x"
        print(f"\nguard={gtxt:>5}  overall within10x={np.mean(np.abs(allr)<1)*100:.0f}%  "
              f"cheapest200 overlap={len(ours_cheap & ml_cheap)}/200")
        for c, name in [(1, "border"), (2, "interior"), (3, "nonman")]:
            arr = np.array(buckets[c])
            print(f"    {name:9s} median={10**np.median(arr):.3g}x within10x={np.mean(np.abs(arr)<1)*100:.0f}%")


if __name__ == "__main__":
    main()
