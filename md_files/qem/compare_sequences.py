#!/usr/bin/env python3
"""Step through ours vs MeshLab collapse sequences in parallel and find where they
first diverge. For each step compares the collapsed undirected edge {v0,v1}. Reports
how many leading collapses match, and characterizes the first divergence (costs,
whether it's a near-flat / clamp-region edge) to tell noise-cascade from a real bug.

Usage:  py compare_sequences.py [base_prefix]
"""
import os
import re
import sys
import json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from heap_compare import DEFAULT_BASE

OUR_CR = "_collapse_records.jsonl"
ML_CR  = "_mat_initial.off_collapse_records.jsonl"


def seq(path):
    out = []
    for line in open(path):
        r = json.loads(re.sub(r'\b(inf|-inf|nan)\b', '1e308', line))
        e = r["edge"]
        a, b = e["v0"]["id"], e["v1"]["id"]
        out.append((min(a, b), max(a, b), r["cost"], r["idx"]))
    return out


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    o = seq(base + OUR_CR); m = seq(base + ML_CR)
    print(f"collapses: ours={len(o)} ml={len(m)}")

    # leading exact match (same undirected edge in same order)
    lead = 0
    for (oa, ob, oc, _), (ma, mb, mc, _) in zip(o, m):
        if (oa, ob) == (ma, mb):
            lead += 1
        else:
            break
    print(f"leading collapses with identical edge+order: {lead}")

    if lead < len(o) and lead < len(m):
        print("\n--- first divergence ---")
        for label, s in (("OURS", o), ("ML", m)):
            a, b, c, idx = s[lead]
            print(f"  {label:4s} idx={idx} edge=({a},{b}) cost={c:.4e}")
        # context: next few on each side
        print("\n  next 5 OURS:", [((a, b), f'{c:.1e}') for a, b, c, _ in o[lead:lead+5]])
        print("  next 5 ML  :", [((a, b), f'{c:.1e}') for a, b, c, _ in m[lead:lead+5]])
        # is the divergence edge a near-flat (clamp-region) collapse?
        oc = o[lead][2]; mc = m[lead][2]
        print(f"\n  divergence costs: ours={oc:.3e} ml={mc:.3e}  "
              f"(both < 1e-6 => near-flat/clamp tie: {'YES' if max(oc,mc) < 1e-6 else 'NO'})")

    # order-independent: how many undirected edges collapsed by BOTH (set overlap)?
    os_ = set((a, b) for a, b, _, _ in o)
    ms_ = set((a, b) for a, b, _, _ in m)
    print(f"\nundirected collapsed-edge set overlap = {len(os_ & ms_)}  "
          f"(ours-only={len(os_-ms_)}, ml-only={len(ms_-os_)})")


if __name__ == "__main__":
    main()
