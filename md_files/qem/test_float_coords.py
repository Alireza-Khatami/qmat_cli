#!/usr/bin/env python3
"""Test whether MeshLab's float coordinates (MESHLAB_SCALAR=float) explain the
flat-region clamp gap. MeshLab stores vertex positions as 32-bit float and builds
the plane normal/offset from float positions (quadric accumulator is double).
We use double throughout. For near-flat fans this changes Apply(mid) and thus
which edges trip the gate / clamp.

Compares per-edge-type agreement with MeshLab under three position precisions:
  double        (current)
  float-pos     (positions round-tripped through float32, quadric math double)
This isolates whether float coords reproduce MeshLab's clamp behavior.

Usage:
    py test_float_coords.py [base_prefix]
"""
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qem_check as qc
from heap_compare import load_meshlab_heap, DEFAULT_BASE


def build_quadrics_prec(V, F, border_pred, to_float):
    """InitQuadric with optional float32 quantization of positions (MeshLab-like).

    to_float: if True, every position is read at float32 precision (cast back to
    float64 for the double quadric math), mirroring MESHLAB_SCALAR=float.
    """
    P = V.astype(np.float32).astype(np.float64) if to_float else V
    diag = np.linalg.norm(P.max(0) - P.min(0))
    SF = 1e8 * (1.0 / diag) ** 6
    vfaces, edgef = defaultdict(list), defaultdict(list)
    for fi, f in enumerate(F):
        for k in range(3):
            vfaces[f[k]].append(fi)
            a, b = f[k], f[(k + 1) % 3]
            edgef[(min(a, b), max(a, b))].append(fi)
    Qd = [qc.Q() for _ in range(len(V))]
    for f in F:
        p0, p1, p2 = P[f[0]], P[f[1]], P[f[2]]
        da = np.cross(p1 - p0, p2 - p0); area = np.linalg.norm(da)
        if area < 1e-30:
            continue
        n = da / area; off = n @ p0
        q = qc.Q(); q.byplane(n, off); q.mul(area)
        for j in range(3):
            Qd[f[j]].iadd(q)
        for j in range(3):
            a, b = f[j], f[(j + 1) % 3]
            if not border_pred(len(edgef[(min(a, b), max(a, b))])):
                continue
            e = P[b] - P[a]; el = np.linalg.norm(e)
            if el < 1e-30:
                continue
            bd = np.cross(n, e / el) * qc.BQW; bo = bd @ P[a]
            bq = qc.Q(); bq.byplane(bd, bo)
            Qd[a].iadd(bq); Qd[b].iadd(bq)
    return P, Qd, vfaces, edgef, SF, diag


def summarize(V, F, mlc, to_float, label):
    P, Qd, vfaces, edgef, SF, diag = build_quadrics_prec(V, F, lambda c: c % 2 == 1, to_float)
    buckets = {1: [], 2: [], 3: []}
    oc = {}
    ml_clamp = ours_clamp = ntot = 0
    for (a, b), fs in edgef.items():
        m = mlc.get((a, b))
        o = qc.cost(P, F, Qd, vfaces, edgef, SF, a, b, verbose=False)
        oc[(a, b)] = o
        if m is None or m > 1e307 or m <= 0 or o <= 0:
            continue
        buckets[min(len(fs), 3)].append(np.log10(o / m))
        ntot += 1
        if m < 1e-15:
            ml_clamp += 1
        if o < 1e-15:
            ours_clamp += 1
    ours_cheap = set(k for k, _ in sorted(oc.items(), key=lambda x: x[1])[:200])
    ml_cheap = set(k for k, _ in sorted(mlc.items(), key=lambda x: x[1])[:200])
    allr = np.concatenate([np.array(buckets[c]) for c in (1, 2, 3)])
    print(f"\n=== positions: {label} (diag={diag:.4f}) ===")
    print(f"  overall within10x={np.mean(np.abs(allr)<1)*100:.0f}%  "
          f"cheapest200 overlap={len(ours_cheap & ml_cheap)}/200  "
          f"clamps: ML={ml_clamp} ours={ours_clamp}")
    for c, name in [(1, "border"), (2, "interior"), (3, "nonman")]:
        arr = np.array(buckets[c])
        print(f"    {name:9s} median={10**np.median(arr):.3g}x within10x={np.mean(np.abs(arr)<1)*100:.0f}%")


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    V, F = qc.load_off(base + "_mat_initial.off")
    mlc = load_meshlab_heap(base + "_mat_initial.obj_initial_heap_state_.json")
    summarize(V, F, mlc, False, "double (current)")
    summarize(V, F, mlc, True, "float32 (MeshLab-like)")


if __name__ == "__main__":
    main()
