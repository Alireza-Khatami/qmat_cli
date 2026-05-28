#!/usr/bin/env python3
"""Compare our ACTUAL C++ initial heap dump against MeshLab's, directly from the
two JSON files (no Python recompute). Earlier analysis used qem_check.py which
solves the optimum with np.linalg.solve/lstsq -- a different solver than the
Eigen fullPivLu our C++ (and MeshLab) actually use, so those numbers were a
model artifact. This reads the real costs both engines emitted.

Edges are undirected; MeshLab emits both directions + per-face dups, so we keep
the min cost per undirected edge. Edge type (border/interior/nonman) comes from
the .off face incidence.

Usage:
    py compare_actual_heaps.py [base_prefix]
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


def load_heap(path):
    txt = open(path).read()
    txt = re.sub(r'\b(inf|-inf|nan)\b', '1e308', txt)
    d = json.loads(txt)
    edge_cost = {}
    for e in d["entries"]:
        a, b = e["v0_id"], e["v1_id"]
        k = (min(a, b), max(a, b))
        c = e["cost"]
        if k not in edge_cost or c < edge_cost[k]:
            edge_cost[k] = c
    return edge_cost


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    ours = load_heap(base + "_initial_heap_state_.json")
    ml = load_heap(base + "_mat_initial.obj_initial_heap_state_.json")
    V, F = qc.load_off(base + "_mat_initial.off")

    facecnt = defaultdict(int)
    for f in F:
        for k in range(3):
            a, b = f[k], f[(k + 1) % 3]
            facecnt[(min(a, b), max(a, b))] += 1

    common = set(ours) & set(ml)
    print(f"ours edges={len(ours)}  ml edges={len(ml)}  common={len(common)}")
    print(f"only in ours={len(set(ours)-set(ml))}  only in ml={len(set(ml)-set(ours))}\n")

    buckets = {1: [], 2: [], 3: []}
    for k in common:
        o, m = ours[k], ml[k]
        if o <= 0 or m <= 0 or m > 1e307 or o > 1e307:
            continue
        buckets[min(facecnt[k], 3)].append(np.log10(o / m))
    allr = np.concatenate([np.array(buckets[c]) for c in (1, 2, 3) if buckets[c]])
    print(f"overall within10x={np.mean(np.abs(allr)<1)*100:.1f}%  "
          f"within2x={np.mean(np.abs(allr)<np.log10(2))*100:.1f}%  median={10**np.median(allr):.4g}x")
    for c, name in [(1, "border"), (2, "interior"), (3, "nonman")]:
        arr = np.array(buckets[c])
        if not len(arr):
            continue
        print(f"  {name:9s} n={len(arr):5d} median={10**np.median(arr):.4g}x "
              f"within10x={np.mean(np.abs(arr)<1)*100:.0f}% within2x={np.mean(np.abs(arr)<np.log10(2))*100:.0f}%")

    for N in (20, 200):
        oc = set(k for k, _ in sorted(ours.items(), key=lambda x: x[1])[:N])
        mc = set(k for k, _ in sorted(ml.items(), key=lambda x: x[1])[:N])
        print(f"cheapest{N} overlap = {len(oc & mc)}/{N}")


if __name__ == "__main__":
    main()
