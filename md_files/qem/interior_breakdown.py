#!/usr/bin/env python3
"""Localize the interior-edge 150x overshoot: is it in QuadErr (quadric/optimal
position) or in the cost divisor (newQual*MinCos)?

Compares, per interior (2-face) edge, our heap cost vs MeshLab's, and prints our
full ComputePriority breakdown (Apply(opt), QuadErr, newQual, MinCos, opt vs mid)
recomputed from the .off geometry. Focuses on edges NON-clamped on both sides
(cost between lo and hi) where the fullPivLu/lstsq solve is well-conditioned, so
the Python recompute matches the C++. Also prints per-type cost percentiles from
the two actual heap dumps (no recompute) to show the true distribution.

Usage:
    py interior_breakdown.py [base_prefix]
"""
import os
import re
import sys
import json
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qem_check as qc
from heap_compare import DEFAULT_BASE
from test_border_rule import build_quadrics

EPS, QTHR, COS_THR = qc.EPS, qc.QTHR, qc.COS_THR


def load_heap(path):
    txt = re.sub(r'\b(inf|-inf|nan)\b', '1e308', open(path).read())
    ec = {}
    for e in json.loads(txt)["entries"]:
        k = (min(e["v0_id"], e["v1_id"]), max(e["v0_id"], e["v1_id"]))
        if k not in ec or e["cost"] < ec[k]:
            ec[k] = e["cost"]
    return ec


def min_closest(A, b, pt, qeps=1e-3):
    """vcg MinimumClosestToPoint: truncate small singular values, pin null-space
    toward pt (the midpoint)."""
    be = -b / 2
    U, s, Vt = np.linalg.svd(A)
    sinv = np.array([1.0 / s[0],
                     1.0 / s[1] if s[1] / s[0] > qeps else 0.0,
                     1.0 / s[2] if s[2] / s[0] > qeps else 0.0])
    Apinv = Vt.T @ np.diag(sinv) @ U.T
    return pt + Apinv @ (be - A @ pt)


def cost_at(V, F, Qd, vfaces, edgef, SF, a, b, opt):
    qq = qc.Q(); qq.iadd(Qd[a]); qq.iadd(Qd[b])
    nq, mc = 1e308, 1e308
    for (va, partner) in ((a, b), (b, a)):
        for fi in vfaces[va]:
            f = F[fi]
            if partner in f:
                continue
            pp = [opt if f[k] == va else V[f[k]] for k in range(3)]
            nq = min(nq, qc.qualityface(*pp))
            nb = qc.trinorm(V[f[0]], V[f[1]], V[f[2]])
            mc = min(mc, float(qc.trinorm(*pp) @ nb))
    nq = min(nq, QTHR); mc = min(mc, COS_THR); mc = abs((mc + 1) / 2)
    return max(SF * qq.apply(opt), EPS) / (nq * mc)


def breakdown(V, F, Qd, vfaces, edgef, SF, a, b):
    qq = qc.Q(); qq.iadd(Qd[a]); qq.iadd(Qd[b])
    mid = (V[a] + V[b]) / 2
    if Qd[a].apply(mid) + Qd[b].apply(mid) > 2 * EPS:
        try:
            opt = np.linalg.solve(qq.A, -qq.b / 2)
        except np.linalg.LinAlgError:
            opt, *_ = np.linalg.lstsq(qq.A, -qq.b / 2, rcond=None)
        gate = "SOLVE"
    else:
        opt, gate = mid, "MID"
    nq, mc = 1e308, 1e308
    for (va, partner) in ((a, b), (b, a)):
        for fi in vfaces[va]:
            f = F[fi]
            if partner in f:
                continue
            pp = [opt if f[k] == va else V[f[k]] for k in range(3)]
            nq = min(nq, qc.qualityface(*pp))
            nb = qc.trinorm(V[f[0]], V[f[1]], V[f[2]])
            mc = min(mc, float(qc.trinorm(*pp) @ nb))
    applyv = qq.apply(opt)
    QuadErr = SF * applyv
    nq = min(nq, QTHR); mc = min(mc, COS_THR); mc = abs((mc + 1) / 2)
    QuadErr = max(QuadErr, EPS)
    cost = QuadErr / (nq * mc)
    optdist = np.linalg.norm(opt - mid); el = np.linalg.norm(V[a] - V[b])
    # what the cost would be at the MIDPOINT (i.e. if the solve were skipped)
    nqm, mcm = 1e308, 1e308
    for (va, partner) in ((a, b), (b, a)):
        for fi in vfaces[va]:
            f = F[fi]
            if partner in f:
                continue
            pp = [mid if f[k] == va else V[f[k]] for k in range(3)]
            nqm = min(nqm, qc.qualityface(*pp))
            nb = qc.trinorm(V[f[0]], V[f[1]], V[f[2]])
            mcm = min(mcm, float(qc.trinorm(*pp) @ nb))
    nqm = min(nqm, QTHR); mcm = min(mcm, COS_THR); mcm = abs((mcm + 1) / 2)
    qem_mid = max(SF * qq.apply(mid), EPS)
    cost_mid = qem_mid / (nqm * mcm)
    gatev = Qd[a].apply(mid) + Qd[b].apply(mid)
    return dict(gate=gate, apply=applyv, QuadErr=QuadErr, nq=nq, mc=mc,
                cost=cost, optdist=optdist, el=el, cost_mid=cost_mid, gatev=gatev)


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    ours = load_heap(base + "_initial_heap_state_.json")
    ml = load_heap(base + "_mat_initial.obj_initial_heap_state_.json")
    V, F = qc.load_off(base + "_mat_initial.off")
    Qd, vfaces, edgef, SF, diag = build_quadrics(V, F, lambda c: c % 2 == 1)

    fc = defaultdict(int)
    for f in F:
        for k in range(3):
            fc[(min(f[k], f[(k+1) % 3]), max(f[k], f[(k+1) % 3]))] += 1

    # per-type cost percentiles from the ACTUAL heaps
    print("cost percentiles (actual heaps)   p10        p50        p90")
    for c, name in [(1, "border"), (2, "interior"), (3, "nonman")]:
        oc = np.array([ours[k] for k in ours if fc[k] == c or (c == 3 and fc[k] >= 3)])
        mc = np.array([ml[k] for k in ml if (fc[k] == c or (c == 3 and fc[k] >= 3)) and ml[k] < 1e307])
        op = np.percentile(oc, [10, 50, 90]); mp = np.percentile(mc, [10, 50, 90])
        print(f"  {name:8s} ours {op[0]:.2e} {op[1]:.2e} {op[2]:.2e}   ml {mp[0]:.2e} {mp[1]:.2e} {mp[2]:.2e}")

    # non-clamped interior edges near the median ratio: dump breakdown
    print("\nnon-clamped interior edges (both 1e-9..1e-2), our breakdown:")
    rows = []
    for k in edgef:
        if fc[k] != 2:
            continue
        o, m = ours.get(k), ml.get(k)
        if o is None or m is None or m > 1e307:
            continue
        if not (1e-9 < o < 1e-2 and 1e-9 < m < 1e-2):
            continue
        rows.append((o / m, k, o, m))
    rows.sort()
    mid_rows = rows[len(rows)//2 - 4: len(rows)//2 + 4] if len(rows) >= 8 else rows
    print(f"  ({len(rows)} such edges; showing 8 near median ratio)")

    # the destructive tail: highest-cost interior edges in OUR heap
    print("\nWORST interior edges in our heap (highest cost), our breakdown:")
    worst = sorted(((ours[k], k) for k in edgef
                    if fc[k] == 2 and k in ours and k in ml and ml[k] < 1e307),
                   reverse=True)[:8]
    print(f"  {'edge':>11} {'ours':>10} {'ml':>10} {'cost@mid':>10} {'cost@SVD':>10}")
    for o, k in worst:
        a, b = k
        qq = qc.Q(); qq.iadd(Qd[a]); qq.iadd(Qd[b])
        mid = (V[a] + V[b]) / 2
        svd_opt = min_closest(qq.A, qq.b, mid)
        c_mid = cost_at(V, F, Qd, vfaces, edgef, SF, a, b, mid)
        c_svd = cost_at(V, F, Qd, vfaces, edgef, SF, a, b, svd_opt)
        print(f"  {a:5d}-{b:<5d} {o:10.2e} {ml[k]:10.2e} {c_mid:10.2e} {c_svd:10.2e}")
    print(f"  {'edge':>11} {'ours/ml':>9} {'gate':>5} {'Apply':>10} {'QuadErr':>10} "
          f"{'newQual':>8} {'MinCos':>7} {'optdist/el':>10}")
    for ratio, k, o, m in mid_rows:
        bd = breakdown(V, F, Qd, vfaces, edgef, SF, k[0], k[1])
        rel = bd['optdist'] / bd['el'] if bd['el'] > 0 else 0
        print(f"  {k[0]:5d}-{k[1]:<5d} {ratio:9.2e} {bd['gate']:>5} {bd['apply']:10.2e} "
              f"{bd['QuadErr']:10.2e} {bd['nq']:8.3f} {bd['mc']:7.3f} {rel:10.2e}")


if __name__ == "__main__":
    main()
