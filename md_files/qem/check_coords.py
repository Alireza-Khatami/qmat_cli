#!/usr/bin/env python3
"""Check whether our QEM coordinates (Pos(v) = ToFloat(center*qemScale), dumped in
_vertex_quadrics_.json as 'pos') are bit-identical to the .off coordinates MeshLab
loads (text -> float32). On flat edges the low bits of the position decide whether
the near-zero Apply trips the QuadricEpsilon clamp, so any mismatch matters.

Reports:
  - max |our_pos - off(double)|  and how many .off decimal digits are present
  - whether our_pos == float32(off)  exactly (what ML actually stores)
  - per-vertex detail for the known divergent flat edges (633,634,631,1293,1459)

Usage:  py check_coords.py [base_prefix]
"""
import os
import re
import sys
import json

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from heap_compare import DEFAULT_BASE


def load_off_raw(path):
    """Return {vid: np.array([x,y,z], float64)} and the raw text tokens (to see digits)."""
    lines = [l.strip() for l in open(path) if l.strip()]
    assert lines[0] == "OFF"
    nv, nf, _ = map(int, lines[1].split())
    V = {}
    raw = {}
    for i in range(nv):
        toks = lines[2 + i].split()
        V[i] = np.array(list(map(float, toks[:3])))
        raw[i] = toks[:3]
    return V, raw


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    V, raw = load_off_raw(base + "_mat_initial.off")
    ov = {v["id"]: np.array(v["pos"]) for v in
          json.loads(open(base + "_vertex_quadrics_.json").read())["verts"]}

    common = sorted(set(V) & set(ov))
    dd = np.array([np.max(np.abs(ov[i] - V[i])) for i in common])   # ours vs off (double)
    f32 = np.array([np.max(np.abs(ov[i].astype(np.float32).astype(np.float64)
                                  - V[i].astype(np.float32).astype(np.float64)))
                    for i in common])
    # does our double pos == float32(off) exactly?
    eqf = sum(1 for i in common
              if np.array_equal(ov[i], V[i].astype(np.float32).astype(np.float64)))
    print(f"common verts={len(common)}")
    print(f"|our_pos - off(double)|:  max={dd.max():.3e}  median={np.median(dd):.3e}  nonzero={np.sum(dd>0)}")
    print(f"|float32(our) - float32(off)|: max={f32.max():.3e}")
    print(f"our_pos == float32(off) exactly: {eqf}/{len(common)}")
    print(f"sample .off raw tokens (digit count): {raw[common[0]]}")

    print("\n--- divergent flat-edge vertices ---")
    for vid in (633, 634, 631, 1293, 1459, 173, 180):
        if vid not in ov:
            continue
        o = ov[vid]; v = V[vid]
        print(f"  v{vid}:")
        print(f"     off  raw = {raw[vid]}")
        print(f"     off  dbl = {v!r}")
        print(f"     our  pos = {o!r}")
        print(f"     diff     = {(o - v)}   (== float32(off)? {np.array_equal(o, v.astype(np.float32).astype(np.float64))})")


if __name__ == "__main__":
    main()
