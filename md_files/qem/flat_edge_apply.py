#!/usr/bin/env python3
"""Pinpoint the flat-edge clamp divergence using ML's recorded breakdown.

For the cheapest edges (the ones that collapse first and set the order), compare
MeshLab's recorded applyOpt / quadErr (from its enriched heap dump) against OUR
applyOpt recomputed from OUR vertex quadrics at the same survivor position. Shows
whether ML's Apply sits ABOVE the QuadricEpsilon clamp threshold while ours sits
below (=> the clamp fires for us but not ML, flipping the cost by orders of
magnitude and reordering the early collapses).

Usage:  py flat_edge_apply.py [base_prefix]
"""
import os
import re
import sys
import json

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qem_check as qc
from heap_compare import DEFAULT_BASE

ML_HEAP = "_mat_initial.off_initial_heap_state_.json"
ML_VQ   = "_mat_initial.off_vertex_quadrics_.json"
OUR_VQ  = "_vertex_quadrics_.json"
EPS = 1e-15


def loadj(p):
    return json.loads(re.sub(r'\b(inf|-inf|nan)\b', '1e308', open(p).read()))


def qfrom(d):
    q = qc.Q(); q.A = np.array([[d["a"][0], d["a"][1], d["a"][2]],
                                [d["a"][1], d["a"][3], d["a"][4]],
                                [d["a"][2], d["a"][4], d["a"][5]]])
    q.b = np.array(d["b"]); q.c = d["c"]; return q


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    mh = loadj(base + ML_HEAP)
    mv = {v["id"]: v for v in loadj(base + ML_VQ)["verts"]}
    ov = {v["id"]: v for v in loadj(base + OUR_VQ)["verts"]}
    V, F = qc.load_off(base + "_mat_initial.off")
    SF = 1e8 * (1.0 / np.linalg.norm(V.max(0) - V.min(0))) ** 6

    # cheapest ML directed entries that carry a breakdown
    ents = [e for e in mh["entries"] if "applyOpt" in e and 0 < e["cost"] < 1e307]
    ents.sort(key=lambda e: e["cost"])
    print(f"SF={SF:.4e}  clamp: QuadErr<=1e-15  <=>  Apply<= {EPS/SF:.3e}\n")
    print(f"{'edge(v0->v1)':>14} {'ml.applyOpt':>12} {'our.applyOpt':>13} "
          f"{'ml.quadErr':>11} {'our.SF*Apply':>12} {'ml.clmp?':>8} {'our.clmp?':>9}")
    ml_by_dir = {(e["v0_id"], e["v1_id"]): e for e in mh["entries"] if "applyOpt" in e}
    for e in ents[:10]:
        v0, v1 = e["v0_id"], e["v1_id"]
        if v0 not in ov or v1 not in ov:
            continue
        q = qc.Q(); q.iadd(qfrom(ov[v0])); q.iadd(qfrom(ov[v1]))
        our_apply = q.apply(V[v1])
        ml_apply = e["applyOpt"]
        print(f"  {v0:5d}->{v1:<5d} {ml_apply:12.3e} {our_apply:13.3e} "
              f"{e['quadErr']:11.3e} {our_apply*SF:12.3e} "
              f"{str((ml_apply*SF)<=EPS):>8} {str((our_apply*SF)<=EPS):>9}")

    # The edges WE rank cheapest but ML costs ~1e-9 (we clamp, ML doesn't):
    print("\n--- edges where OUR clamp fires but ML's doesn't ---")
    for (v0, v1) in [(633, 634), (633, 631), (1293, 1459), (173, 180)]:
        if v0 not in ov or v1 not in ov:
            continue
        q = qc.Q(); q.iadd(qfrom(ov[v0])); q.iadd(qfrom(ov[v1]))
        our_apply = q.apply(V[v1])
        e = ml_by_dir.get((v0, v1))
        ml_apply = e["applyOpt"] if e else float('nan')
        print(f"  {v0}->{v1}: ml.applyOpt={ml_apply:.3e}  our.applyOpt={our_apply:.3e}  "
              f"(ml clamp={ (ml_apply*SF)<=EPS }, our clamp={ (our_apply*SF)<=EPS })")
        # per-vertex quadric: do 633/634 quadrics themselves match ML?
        for vv in (v0, v1):
            oa = np.array(ov[vv]["a"]); ma = np.array(mv[vv]["a"])
            print(f"     v{vv}: max|a_our-a_ml|={np.max(np.abs(oa-ma)):.3e}  "
                  f"our.pos={np.array(ov[vv]['pos']).round(6)}  off.V={V[vv].round(6)}")


if __name__ == "__main__":
    main()
