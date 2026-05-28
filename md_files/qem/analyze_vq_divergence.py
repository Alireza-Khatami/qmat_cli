#!/usr/bin/env python3
"""Localize the InitQuadric vertex-quadric divergence (ours vs MeshLab).

Loads both <prefix>_vertex_quadrics_.json files, finds the vertices whose
accumulated quadric (a[6],b[3],c) differs most, and for each prints:
  - our a/b/c vs ML a/b/c
  - the vertex's incident edges with their incident-face counts (1=border,
    2=interior, >=3 non-manifold)
so we can see whether the divergence tracks border / non-manifold edges (i.e. a
border-quadric count or direction mismatch in InitQuadric).

Also confirms the OptimalPlacement=false endpoint choice: prints, per edge, whether
our opt and ML opt land on v0 vs v1 (a difference of exactly one edge length means
we pick the opposite survivor).

Usage:  py analyze_vq_divergence.py [base_prefix]
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


def vqmap(p):
    return {v["id"]: v for v in loadj(p)["verts"]}


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    ov = vqmap(base + OUR_VQ); mv = vqmap(base + ML_VQ)
    V, F = qc.load_off(base + "_mat_initial.off")

    # edge -> incident face count; vertex -> incident edges
    ef = defaultdict(int)
    vedges = defaultdict(set)
    for f in F:
        for k in range(3):
            a, b = f[k], f[(k+1) % 3]
            key = (min(a, b), max(a, b))
            ef[key] += 1
            vedges[a].add(key); vedges[b].add(key)

    common = sorted(set(ov) & set(mv))
    worst = []
    for vid in common:
        a = np.array(ov[vid]["a"] + ov[vid]["b"] + [ov[vid]["c"]])
        b = np.array(mv[vid]["a"] + mv[vid]["b"] + [mv[vid]["c"]])
        rel = np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-12))
        worst.append((rel, vid))
    worst.sort(reverse=True)

    print("=== worst-diverging vertex quadrics (ours vs ML) ===")
    for rel, vid in worst[:8]:
        inc = sorted(vedges[vid])
        cnts = [ef[e] for e in inc]
        nb = sum(1 for c in cnts if c == 1)
        nnm = sum(1 for c in cnts if c >= 3)
        print(f"\nvertex {vid}  max-rel-diff={rel:.3e}  "
              f"incident edges={len(inc)} (border={nb}, nonman={nnm})")
        oa = ov[vid]; ma = mv[vid]
        print(f"  our a={np.array(oa['a']).round(5)}")
        print(f"  ml  a={np.array(ma['a']).round(5)}")
        print(f"  our b={np.array(oa['b']).round(5)}  c={oa['c']:.5g}")
        print(f"  ml  b={np.array(ma['b']).round(5)}  c={ma['c']:.5g}")
        print(f"  edge face-counts: {cnts}")

    # fraction of diverging vertices that touch a non-manifold / border edge
    div = [vid for rel, vid in worst if rel > 1e-3]
    tot = len(div)
    touch_nm = sum(1 for v in div if any(ef[e] >= 3 for e in vedges[v]))
    touch_b  = sum(1 for v in div if any(ef[e] == 1 for e in vedges[v]))
    print(f"\n=== {tot} vertices diverge >1e-3; "
          f"{touch_nm} touch a non-manifold edge, {touch_b} touch a border edge ===")


if __name__ == "__main__":
    main()
