#!/usr/bin/env python3
"""Test the vcg FaceBorderFromVF border rule: an edge is a border iff its incident
face count is ODD (not just ==1). Non-manifold edges with 3 faces are therefore
borders and get border quadrics, which we currently omit.

Compares per-edge-type cost agreement with MeshLab under two border rules:
  count==1   (our current C++ rule)
  count%2==1 (faithful vcg rule)

Usage:
    py test_border_rule.py [base_prefix]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qem_check as qc
from heap_compare import load_meshlab_heap, DEFAULT_BASE


def build_quadrics(V, F, border_pred):
    """InitQuadric with a configurable border predicate on incident-face count."""
    diag = np.linalg.norm(V.max(0) - V.min(0))
    SF = 1e8 * (1.0 / diag) ** 6
    from collections import defaultdict
    vfaces, edgef = defaultdict(list), defaultdict(list)
    for fi, f in enumerate(F):
        for k in range(3):
            vfaces[f[k]].append(fi)
            a, b = f[k], f[(k + 1) % 3]
            edgef[(min(a, b), max(a, b))].append(fi)

    Qd = [qc.Q() for _ in range(len(V))]
    for f in F:
        p0, p1, p2 = V[f[0]], V[f[1]], V[f[2]]
        da = np.cross(p1 - p0, p2 - p0); area = np.linalg.norm(da)
        if area < 1e-30:
            continue
        n = da / area; off = n @ p0
        q = qc.Q(); q.byplane(n, off); q.mul(area)
        for j in range(3):
            Qd[f[j]].iadd(q)
        for j in range(3):
            a, b = f[j], f[(j + 1) % 3]
            cnt = len(edgef[(min(a, b), max(a, b))])
            if not border_pred(cnt):
                continue
            e = V[b] - V[a]; el = np.linalg.norm(e)
            if el < 1e-30:
                continue
            bd = np.cross(n, e / el) * qc.BQW; bo = bd @ V[a]
            bq = qc.Q(); bq.byplane(bd, bo)
            Qd[a].iadd(bq); Qd[b].iadd(bq)
    return Qd, vfaces, edgef, SF, diag


def summarize(V, F, border_pred, mlc, label):
    Qd, vfaces, edgef, SF, diag = build_quadrics(V, F, border_pred)
    buckets = {1: [], 2: [], 3: []}
    for (a, b), fs in edgef.items():
        m = mlc.get((a, b))
        o = qc.cost(V, F, Qd, vfaces, edgef, SF, a, b, verbose=False)
        if m is None or m > 1e307 or m <= 0 or o <= 0:
            continue
        c = min(len(fs), 3)
        buckets.setdefault(c, []).append(np.log10(o / m))
    # all-edge cost map + cheapest-set overlap with MeshLab
    oc = {(a, b): qc.cost(V, F, Qd, vfaces, edgef, SF, a, b, verbose=False)
          for (a, b) in edgef}
    allr = np.array([np.log10(oc[k] / mlc[k]) for k in oc
                     if k in mlc and 0 < mlc[k] < 1e307 and oc[k] > 0])
    ours_cheap = set(k for k, _ in sorted(oc.items(), key=lambda x: x[1])[:200])
    ml_cheap = set(k for k, _ in sorted(mlc.items(), key=lambda x: x[1])[:200])
    print(f"\n=== border rule: {label} ===")
    print(f"  OVERALL within10x={np.mean(np.abs(allr) < 1) * 100:.0f}%  "
          f"cheapest-200 overlap={len(ours_cheap & ml_cheap)}/200")
    for c, name in [(1, "border(1)"), (2, "interior(2)"), (3, "nonman(>=3)")]:
        arr = np.array(buckets[c])
        print(f"  {name:14s} n={len(arr):5d}  median ours/ml={10 ** np.median(arr):.3g}x  "
              f"within10x={np.mean(np.abs(arr) < 1) * 100:.0f}%")
    # sample nonman edges
    print("  sample nonman edges (ours vs MeshLab):")
    shown = 0
    for (a, b), fs in edgef.items():
        if len(fs) <= 2:
            continue
        m = mlc.get((a, b))
        if m is None or m > 1e307:
            continue
        o = qc.cost(V, F, Qd, vfaces, edgef, SF, a, b, verbose=False)
        print(f"    {a}-{b} ({len(fs)}f): ours={o:.4e} ml={m:.4e} ratio={o / m:.2e}")
        shown += 1
        if shown >= 4:
            break


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    V, F = qc.load_off(base + "_mat_initial.off")
    mlc = load_meshlab_heap(base + "_mat_initial.obj_initial_heap_state_.json")
    summarize(V, F, lambda c: c == 1, mlc, "count==1 (current C++)")
    summarize(V, F, lambda c: c % 2 == 1, mlc, "count odd (vcg FaceBorderFromVF)")


if __name__ == "__main__":
    main()
