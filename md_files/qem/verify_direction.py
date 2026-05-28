#!/usr/bin/env python3
"""Verify the OptimalPlacement=false divergence is the collapse DIRECTION.

1. Quadrics: report ABSOLUTE max diff (not relative) to confirm the earlier
   'divergence' was near-zero float noise.
2. ML heap: does it contain BOTH directions per edge (IsSymmetric=false seeds
   v0->v1 AND v1->v0)?  Compare raw entry count vs unique undirected edges, and
   check for (a,b) and (b,a) pairs with different cost.
3. Costs: ours (one minmax direction) vs ML (min over both directions).

Usage:  py verify_direction.py [base_prefix]
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

ML_VQ   = "_mat_initial.off_vertex_quadrics_.json"
OUR_VQ  = "_vertex_quadrics_.json"
ML_HEAP = "_mat_initial.off_initial_heap_state_.json"
OUR_HEAP= "_initial_heap_state_.json"


def loadj(p):
    return json.loads(re.sub(r'\b(inf|-inf|nan)\b', '1e308', open(p).read()))


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE

    # 1. absolute quadric diff
    ov = {v["id"]: v for v in loadj(base + OUR_VQ)["verts"]}
    mv = {v["id"]: v for v in loadj(base + ML_VQ)["verts"]}
    maxabs = 0.0
    for vid in set(ov) & set(mv):
        a = np.array(ov[vid]["a"] + ov[vid]["b"] + [ov[vid]["c"]])
        b = np.array(mv[vid]["a"] + mv[vid]["b"] + [mv[vid]["c"]])
        maxabs = max(maxabs, np.max(np.abs(a - b)))
    print(f"1. vertex quadrics: max ABSOLUTE diff over all a/b/c = {maxabs:.3e}")
    print(f"   (quadric magnitudes are ~1e-3..1e-2; so this is {'NOISE -> MATCH' if maxabs < 1e-6 else 'REAL'})")

    # 2. ML heap directions
    mh = loadj(base + ML_HEAP)["entries"]
    oh = loadj(base + OUR_HEAP)["entries"]
    print(f"\n2. heap entry counts: ours={len(oh)}  ml={len(mh)}")
    ml_dir = {}
    for e in mh:
        ml_dir[(e["v0_id"], e["v1_id"])] = e["cost"]
    ml_und = set((min(a, b), max(a, b)) for (a, b) in ml_dir)
    print(f"   ml directed entries={len(ml_dir)}  unique undirected={len(ml_und)}")
    both = 0; costdiff = []
    for (a, b) in ml_dir:
        if (b, a) in ml_dir:
            both += 1
            c0, c1 = ml_dir[(a, b)], ml_dir[(b, a)]
            if c0 > 0 and c1 > 0 and max(c0, c1) < 1e307:
                costdiff.append(max(c0, c1) / min(c0, c1))
    print(f"   ml edges present in BOTH directions: {both//2} (of {len(ml_und)})")
    if costdiff:
        cd = np.array(costdiff)
        print(f"   when both dirs present, cost ratio max/min: median={np.median(cd):.2f} p90={np.percentile(cd,90):.2f} max={cd.max():.2f}")
        print("   => direction MATTERS (the two endpoints give different cost)" if np.median(cd) > 1.01
              else "   => direction doesn't change cost")

    # 3. cost agreement ours vs ML (dedup ml to min over directions)
    ml_min = {}
    for (a, b), c in ml_dir.items():
        k = (min(a, b), max(a, b))
        if k not in ml_min or c < ml_min[k]:
            ml_min[k] = c
    our = {(min(e["v0_id"], e["v1_id"]), max(e["v0_id"], e["v1_id"])): e["cost"] for e in oh}
    common = [k for k in our if k in ml_min and 0 < ml_min[k] < 1e307 and our[k] > 0]
    r = np.array([our[k] / ml_min[k] for k in common])
    print(f"\n3. cost ours/ml (ml=min over both dirs), n={len(common)}: "
          f"median={np.median(r):.3f} within2x={np.mean((r>0.5)&(r<2))*100:.0f}% "
          f"exact(==)={np.mean(np.abs(r-1)<1e-6)*100:.0f}%")
    # how often does our minmax direction equal ml's cheaper direction?
    agree = 0; tot = 0
    for k in common:
        a, b = k
        if (a, b) in ml_dir and (b, a) in ml_dir:
            tot += 1
            cheaper_is_ab = ml_dir[(a, b)] <= ml_dir[(b, a)]
            # ours always seeds (min,max)=(a,b) -> collapses to v1=b ; compare cost to ml (a,b)
            if abs(our[k] - ml_dir[(a, b)]) <= 1e-6 * max(our[k], 1e-30):
                agree += 1
    if tot:
        print(f"   our cost matches ml's (a=min -> b=max) directed cost on {agree}/{tot} two-dir edges")


if __name__ == "__main__":
    main()
