#!/usr/bin/env python3
"""Replicate vcg's TriEdgeCollapseQuadric::ComputePriority on a .off mesh.

Builds the vcg quadrics (area-weighted face planes + border quadrics) exactly as
vcglib/tri_edge_collapse_quadric.h InitQuadric does, computes the heap cost for a
given edge, and prints the intermediate values. Used to compare our QEM cost
against MeshLab's per-edge cost (see md_files/qem/initial_heap_state_schema.md).

Finding: this faithful re-implementation reproduces OUR C++ costs, not MeshLab's.
MeshLab clamps interior flat-edge QuadErr to QuadricEpsilon (cost ~9e-17); we do
not, because our ScaleFactor (1e8/diag^6, diag=1.414 -> 1.25e7) is ~330x too large
to trigger the clamp. MeshLab behaves as if ScaleFactor ~= 3.8e4 (diag ~= 3.72).

Usage:
    py qem_check.py [mesh.off] [v0 v1 ...]
Defaults to the 01_00040057... test mesh and a few diagnostic edges.
"""
import numpy as np
import sys
from collections import defaultdict

# ── parameters (MeshLab QEM defaults: Preserve Normal ON, Preserve Topology ON) ──
BQW     = 0.5                  # effective BoundaryQuadricWeight (0.5 struct * 1.0 UI)
EPS     = 1e-15                # QuadricEpsilon
QTHR    = 0.3                  # QualityThr
COS_THR = np.cos(np.pi / 4)    # NormalCheck -> NormalThrRad = pi/4


def load_off(path):
    lines = [l.strip() for l in open(path) if l.strip()]
    assert lines[0] == "OFF", "not an OFF file"
    nv, nf, _ = map(int, lines[1].split())
    V = np.array([list(map(float, lines[2 + i].split())) for i in range(nv)])
    F = [list(map(int, lines[2 + nv + i].split()))[1:] for i in range(nf)]
    return V, F


class Q:
    """vcg math::Quadric<double>: f(x) = x^T A x + b.x + c."""
    def __init__(s):
        s.A = np.zeros((3, 3)); s.b = np.zeros(3); s.c = 0.0
    def byplane(s, n, off):
        s.A = np.outer(n, n); s.b = -2 * off * n; s.c = off * off
    def iadd(s, o):
        s.A += o.A; s.b += o.b; s.c += o.c
    def mul(s, w):
        s.A *= w; s.b *= w; s.c *= w
    def apply(s, p):
        return float(p @ s.A @ p + s.b @ p + s.c)


def build(V, F):
    """InitQuadric: per-vertex accumulated quadrics, ScaleFactor, adjacency."""
    diag = np.linalg.norm(V.max(0) - V.min(0))
    SF = 1e8 * (1.0 / diag) ** 6

    vfaces = defaultdict(list)
    edgef = defaultdict(list)
    for fi, f in enumerate(F):
        for k in range(3):
            vfaces[f[k]].append(fi)
            a, b = f[k], f[(k + 1) % 3]
            edgef[(min(a, b), max(a, b))].append(fi)

    def isborder(a, b):
        return len(edgef[(min(a, b), max(a, b))]) == 1

    Qd = [Q() for _ in range(len(V))]
    for f in F:
        p0, p1, p2 = V[f[0]], V[f[1]], V[f[2]]
        da = np.cross(p1 - p0, p2 - p0)
        area = np.linalg.norm(da)
        if area < 1e-30:
            continue
        n = da / area; off = n @ p0
        q = Q(); q.byplane(n, off); q.mul(area)        # UseArea
        for j in range(3):
            Qd[f[j]].iadd(q)
        for j in range(3):                              # border quadrics
            a, b = f[j], f[(j + 1) % 3]
            if not isborder(a, b):
                continue
            e = V[b] - V[a]; el = np.linalg.norm(e)
            if el < 1e-30:
                continue
            bd = np.cross(n, e / el) * BQW; bo = bd @ V[a]
            bq = Q(); bq.byplane(bd, bo)
            Qd[a].iadd(bq); Qd[b].iadd(bq)
    return Qd, vfaces, edgef, SF, diag


def qualityface(p0, p1, p2):
    d10, d20, d12 = p1 - p0, p2 - p0, p1 - p2
    a = np.linalg.norm(np.cross(d10, d20))
    if a == 0:
        return 0.0
    b = max(d10 @ d10, d20 @ d20, d12 @ d12)
    return a / b if b > 0 else 0.0


def trinorm(p0, p1, p2):
    n = np.cross(p1 - p0, p2 - p0); l = np.linalg.norm(n)
    return n / l if l > 1e-30 else np.zeros(3)


def cost(V, F, Qd, vfaces, edgef, SF, v0, v1, verbose=True):
    def isborder(a, b):
        return len(edgef[(min(a, b), max(a, b))]) == 1
    qq = Q(); qq.iadd(Qd[v0]); qq.iadd(Qd[v1])
    mid = (V[v0] + V[v1]) / 2
    # ComputePosition (OptimalPlacement, SVDPlacement=false -> Minimum, used regardless)
    if Qd[v0].apply(mid) + Qd[v1].apply(mid) > 2 * EPS:
        try:
            opt = np.linalg.solve(qq.A, -qq.b / 2)
        except np.linalg.LinAlgError:
            opt, _, _, _ = np.linalg.lstsq(qq.A, -qq.b / 2, rcond=None)
        gate = "SOLVE"
    else:
        opt = mid; gate = "MID"
    # surviving faces: faces of v0 not containing v1, faces of v1 not containing v0
    newQual = 1e308; mincos = 1e308
    for (va, partner) in ((v0, v1), (v1, v0)):
        for fi in vfaces[va]:
            f = F[fi]
            if partner in f:
                continue
            pp = [opt if f[k] == va else V[f[k]] for k in range(3)]
            newQual = min(newQual, qualityface(*pp))
            nb = trinorm(V[f[0]], V[f[1]], V[f[2]])
            na = trinorm(*pp)
            mincos = min(mincos, float(na @ nb))
    QuadErr = SF * qq.apply(opt)
    if newQual > QTHR:
        newQual = QTHR
    if mincos > COS_THR:
        mincos = COS_THR
    mincos = abs((mincos + 1) / 2)
    QuadErr = max(QuadErr, EPS)
    if QuadErr <= EPS:                                  # degenerate-edge tie-break
        QuadErr *= np.linalg.norm(V[v0] - V[v1])
    err = QuadErr / (newQual * mincos)
    if verbose:
        print(f"  edge {v0}-{v1} border={isborder(v0,v1)} gate={gate} "
              f"Apply(opt)={qq.apply(opt):.3e} QuadErr={QuadErr:.3e} "
              f"newQual={newQual:.3e} mincos={mincos:.3e} COST={err:.4e}")
    return err


def main():
    off = sys.argv[1] if len(sys.argv) > 1 else (
        "01_00040057_f8f78dbd17414efda75bc437_trimesh_000_mat_initial.off")
    V, F = load_off(off)
    Qd, vfaces, edgef, SF, diag = build(V, F)
    print(f"mesh={off}")
    print(f"nv={len(V)} nf={len(F)} bboxdiag={diag:.6f} ScaleFactor={SF:.6e} "
          f"BQW={BQW} COS_THR={COS_THR:.4f}")

    if len(sys.argv) > 2:
        ids = list(map(int, sys.argv[2:]))
        for a, b in zip(ids[0::2], ids[1::2]):
            cost(V, F, Qd, vfaces, edgef, SF, a, b)
    else:
        # diagnostic edges with known MeshLab / our values
        for v0, v1, ml, ours in [
            (199, 205, 9.15e-17, 1.28e-12),
            (1293, 1459, 1.29e-08, 7.29e-17),
            (80, 111, 9.27e-17, 5.54e-13),
        ]:
            print(f"edge {v0}-{v1}: MeshLab={ml:.2e}  ours={ours:.2e}")
            cost(V, F, Qd, vfaces, edgef, SF, v0, v1)

        # ScaleFactor reconciliation: MeshLab clamps interior edges; our clamped
        # formula reproduces MeshLab exactly, but our SF is too big to clamp.
        print("\n=== ScaleFactor reconciliation (interior clamped edges) ===")
        for v0, v1, mlcost in [(199, 205, 9.15445085251621e-17),
                               (80, 111, 9.2734775104511604e-17)]:
            qq = Q(); qq.iadd(Qd[v0]); qq.iadd(Qd[v1])
            applymid = qq.apply((V[v0] + V[v1]) / 2)
            el = np.linalg.norm(V[v0] - V[v1])
            mincos = abs((COS_THR + 1) / 2)
            clamped = EPS * el / (QTHR * mincos)
            sf_thr = EPS / applymid if applymid > 0 else float('inf')
            print(f"edge {v0}-{v1}: Apply(mid)={applymid:.3e} our SF*Apply={SF*applymid:.3e} "
                  f"-> clamps if SF<={sf_thr:.3e}")
            print(f"          clamped prediction={clamped:.4e}  MeshLab={mlcost:.4e}  "
                  f"ratio={clamped/mlcost:.3f}")
        print(f"\nour ScaleFactor={SF:.4e} (diag={diag:.4f}); "
              f"diag giving SF=3.8e4 is {(1e8/3.8e4)**(1/6):.3f}")


if __name__ == "__main__":
    main()
