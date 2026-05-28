#!/usr/bin/env python3
"""Directed per-edge cost comparison (ours vs MeshLab) now that both sides seed
both directions in OptimalPlacement=false mode. Compares cost(v0->v1) to the
SAME directed (v0->v1) entry in ML (not min-over-directions), which is the real
faithfulness test. Also reports the cheapest-edge order overlap and the final
collapse / vertex counts from the collapse-record logs.

Usage:  py directed_cost.py [base_prefix]
"""
import os
import re
import sys
import json

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from heap_compare import DEFAULT_BASE

OUR_HEAP = "_initial_heap_state_.json"
ML_HEAP  = "_mat_initial.off_initial_heap_state_.json"
OUR_CR   = "_collapse_records.jsonl"
ML_CR    = "_mat_initial.off_collapse_records.jsonl"


def loadj(p):
    return json.loads(re.sub(r'\b(inf|-inf|nan)\b', '1e308', open(p).read()))


def directed(j):
    """(v0,v1) -> min cost for that DIRECTED pair (ML has dup/stale entries)."""
    d = {}
    for e in j["entries"]:
        k = (e["v0_id"], e["v1_id"])
        if k not in d or e["cost"] < d[k]:
            d[k] = e["cost"]
    return d


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    od = directed(loadj(base + OUR_HEAP))
    md = directed(loadj(base + ML_HEAP))
    common = [k for k in od if k in md and 0 < md[k] < 1e307 and 0 < od[k] < 1e307]
    r = np.array([od[k] / md[k] for k in common])
    print(f"DIRECTED cost ours/ml  (n={len(common)} of ours={len(od)}, ml={len(md)})")
    print(f"  median={np.median(r):.5f}  within2x={np.mean((r>.5)&(r<2))*100:.0f}%  "
          f"exact(<1e-6)={np.mean(np.abs(r-1)<1e-6)*100:.0f}%  "
          f"within1pct={np.mean(np.abs(r-1)<0.01)*100:.0f}%")
    bad = sorted(((abs(np.log(od[k]/md[k])), k) for k in common), reverse=True)[:6]
    print("  worst directed mismatches (edge: ours vs ml):")
    for _, k in bad:
        print(f"    {k}: ours={od[k]:.4e} ml={md[k]:.4e} ratio={od[k]/md[k]:.3e}")

    # cheapest-edge overlap (collapse order proxy), directed
    n = 50
    oc = set(k for k, _ in sorted(od.items(), key=lambda x: x[1])[:n])
    mc = set(k for k, _ in sorted(md.items(), key=lambda x: x[1])[:n])
    print(f"\ncheapest-{n} directed overlap = {len(oc & mc)}/{n}")

    # final counts from collapse logs
    def count(p):
        return sum(1 for _ in open(p)) if os.path.exists(p) else -1
    print(f"\ncollapses: ours={count(base + OUR_CR)}  ml={count(base + ML_CR)}")


if __name__ == "__main__":
    main()
