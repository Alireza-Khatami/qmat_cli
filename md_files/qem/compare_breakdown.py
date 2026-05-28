#!/usr/bin/env python3
"""Stage-by-stage QEM comparison: ours vs MeshLab, using the enriched dumps.

After rebuilding BOTH MeshLab (with the collapse_logger optPos + breakdown +
vertex-quadric instrumentation) and QMAT, each side writes:
  <prefix>_initial_heap_state_.json   now carries per-edge "opt" (+ ML breakdown)
  <prefix>_vertex_quadrics_.json      per-vertex InitQuadric quadric (a[6],b[3],c)

This script localizes the divergence to one stage:
  1. VERTEX QUADRICS  — if a/b/c differ, InitQuadric itself diverges.
  2. OPTIMAL POSITION — if quadrics match but opt differs, it's placement (the
     rank-deficient border solve flying).
  3. DIVISOR          — if opt matches but cost differs, it's newQual / MinCos.

Usage:  py compare_breakdown.py [base_prefix]
"""
import os
import re
import sys
import json
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from heap_compare import DEFAULT_BASE

OURS_HEAP = "_initial_heap_state_.json"
ML_HEAP   = "_mat_initial.off_initial_heap_state_.json"
OURS_VQ   = "_vertex_quadrics_.json"
ML_VQ     = "_mat_initial.off_vertex_quadrics_.json"


def loadj(path):
    if not os.path.exists(path):
        print(f"  !! missing: {path}")
        return None
    return json.loads(re.sub(r'\b(inf|-inf|nan)\b', '1e308', open(path).read()))


def edge_map(j):
    """key (min,max) -> dict(cost, opt, ...extra ML breakdown if present)."""
    out = {}
    for e in j["entries"]:
        k = (min(e["v0_id"], e["v1_id"]), max(e["v0_id"], e["v1_id"]))
        if k in out and e["cost"] >= out[k]["cost"]:
            continue
        out[k] = e
    return out


def vq_map(j):
    return {v["id"]: v for v in j["verts"]} if j else {}


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE

    # ---- 1. vertex quadrics ----
    print("=== 1. VERTEX QUADRICS (InitQuadric) ===")
    ovq = vq_map(loadj(base + OURS_VQ))
    mvq = vq_map(loadj(base + ML_VQ))
    if ovq and mvq:
        common = sorted(set(ovq) & set(mvq))
        print(f"  ours={len(ovq)} ml={len(mvq)} common={len(common)}")
        worst = []
        for vid in common:
            a = np.array(ovq[vid]["a"] + ovq[vid]["b"] + [ovq[vid]["c"]])
            b = np.array(mvq[vid]["a"] + mvq[vid]["b"] + [mvq[vid]["c"]])
            denom = np.maximum(np.abs(b), 1e-12)
            rel = np.max(np.abs(a - b) / denom)
            worst.append((rel, vid))
        worst.sort(reverse=True)
        rels = np.array([w[0] for w in worst])
        print(f"  max rel-diff over a/b/c: median={np.median(rels):.2e} "
              f"p99={np.percentile(rels,99):.2e} max={rels.max():.2e}")
        print(f"  worst 5 vertices: {[(v, f'{r:.1e}') for r,v in worst[:5]]}")
        print("  => quadrics MATCH" if rels.max() < 1e-4 else
              "  => quadrics DIVERGE (InitQuadric bug)")

    # ---- 2. optimal position ----
    print("\n=== 2. OPTIMAL POSITION (placement) ===")
    oh = loadj(base + OURS_HEAP); mh = loadj(base + ML_HEAP)
    if not (oh and mh):
        return
    oe = edge_map(oh); me = edge_map(mh)
    common = sorted(set(oe) & set(me))
    print(f"  common edges={len(common)}")
    has_opt = "opt" in next(iter(oe.values())) and "opt" in next(iter(me.values()))
    if not has_opt:
        print("  !! 'opt' field absent — rebuild both sides with the new instrumentation.")
        return

    # edge incidence (border/interior/nonman) from the .off
    from qem_check import load_off
    V, F = load_off(base + "_mat_initial.off")
    fc = defaultdict(int)
    for f in F:
        for k in range(3):
            fc[(min(f[k], f[(k+1)%3]), max(f[k], f[(k+1)%3]))] += 1

    def bucket(k):
        c = fc.get(k, 0)
        return "border" if c == 1 else "interior" if c == 2 else "nonman"

    rows = defaultdict(list)
    for k in common:
        o = np.array(oe[k]["opt"]); m = np.array(me[k]["opt"])
        el = np.linalg.norm(V[k[0]] - V[k[1]]) if k[0] < len(V) and k[1] < len(V) else 1.0
        d = np.linalg.norm(o - m) / el if el > 0 else np.linalg.norm(o - m)
        rows[bucket(k)].append((d, k))
    for name in ("border", "interior", "nonman"):
        arr = np.array([r[0] for r in rows[name]]) if rows[name] else np.array([0.])
        print(f"  {name:8s} n={len(rows[name]):5d}  |our_opt - ml_opt|/el: "
              f"median={np.median(arr):.3e} p90={np.percentile(arr,90):.3e} max={arr.max():.3e}")
        worst = sorted(rows[name], reverse=True)[:3]
        print(f"            worst: {[(k, f'{d:.2f}') for d,k in worst]}")

    # ---- 3. divisor (only ML carries the breakdown; compare to its own cost) ----
    print("\n=== 3. ML BREAKDOWN sanity (cost == quadErr/(newQual*minCos)?) ===")
    if "quadErr" in next(iter(me.values())):
        bad = 0
        for k in common[:5] + common[-5:]:
            e = me[k]
            recon = e["quadErr"] / (e["newQual"] * e["minCos"]) if e["newQual"]*e["minCos"] else float('inf')
            r = recon / e["cost"] if e["cost"] else 0
            print(f"  edge {k} cost={e['cost']:.3e} quadErr={e['quadErr']:.3e} "
                  f"newQual={e['newQual']:.3e} minCos={e['minCos']:.3e} recon/cost={r:.3f}")
    else:
        print("  (ML breakdown fields absent — rebuild ML)")


if __name__ == "__main__":
    main()
