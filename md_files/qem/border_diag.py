#!/usr/bin/env python3
"""Per-edge breakdown for border / non-manifold edges, vs MeshLab's stored cost.

For a sample of border (1-face) and non-manifold (>=3-face) edges, prints the full
ComputePriority breakdown under the vcg parity border rule (our C++ rule): the
placement gate, the optimal position and how far it lands from the midpoint
(opt-dist / edge-len), Apply(opt), newQual, MinCos, our cost, and MeshLab's cost.

Goal: localize the border/nonman divergence to PLACEMENT (gate=SOLVE flies) vs
the QUADRIC/divisor (gate=MID, geometry identical to ML => only quadric can move
the cost).

Usage:  py border_diag.py [base_prefix]
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


def breakdown(V, F, Qd, vfaces, edgef, SF, a, b):
    qq = qc.Q(); qq.iadd(Qd[a]); qq.iadd(Qd[b])
    mid = (V[a] + V[b]) / 2
    gatev = Qd[a].apply(mid) + Qd[b].apply(mid)
    if gatev > 2 * EPS:
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
    QuadErr = max(SF * applyv, EPS)
    nqc = min(nq, QTHR); mcc = min(mc, COS_THR); mcc = abs((mcc + 1) / 2)
    cost = QuadErr / (nqc * mcc)
    el = np.linalg.norm(V[a] - V[b])
    optd = np.linalg.norm(opt - mid)
    # also cost at the midpoint (placement-free)
    applymid = qq.apply(mid)
    return dict(gate=gate, opt=opt, optd=optd, el=el, apply=applyv,
                applymid=applymid, QuadErr=QuadErr, nq=nqc, mc=mcc, cost=cost,
                gatev=gatev)


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    ml = load_heap(base + "_mat_initial.obj_initial_heap_state_.json")
    V, F = qc.load_off(base + "_mat_initial.off")
    Qd, vfaces, edgef, SF, diag = build_quadrics(V, F, lambda c: c % 2 == 1)
    print(f"diag={diag:.4f} SF={SF:.4e}")

    fc = {k: len(v) for k, v in edgef.items()}

    def show(title, pred, n=8):
        print(f"\n=== {title} ===")
        print(f"  {'edge':>11} {'nf':>2} {'gate':>5} {'optd/el':>9} {'Apply(opt)':>11} "
              f"{'Apply(mid)':>11} {'newQual':>8} {'MinCos':>7} {'ours':>10} {'ml':>10} {'o/ml':>8}")
        shown = 0
        for k in sorted(edgef):
            if not pred(fc[k]) or k not in ml or ml[k] > 1e307:
                continue
            a, b = k
            d = breakdown(V, F, Qd, vfaces, edgef, SF, a, b)
            r = d['cost'] / ml[k] if ml[k] else 0
            od = d['optd'] / d['el'] if d['el'] else 0
            print(f"  {a:5d}-{b:<5d} {fc[k]:2d} {d['gate']:>5} {od:9.2e} "
                  f"{d['apply']:11.3e} {d['applymid']:11.3e} {d['nq']:8.3f} {d['mc']:7.3f} "
                  f"{d['cost']:10.3e} {ml[k]:10.3e} {r:8.2e}")
            shown += 1
            if shown >= n:
                break

    show("BORDER edges (1 face)", lambda c: c == 1)
    show("NONMAN edges (3 faces)", lambda c: c == 3)


if __name__ == "__main__":
    main()
