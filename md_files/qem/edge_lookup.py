#!/usr/bin/env python3
"""Cross-check a single edge's cost across three sources:
  1. our C++ heap dump   (<base>_initial_heap_state_.json)
  2. MeshLab heap dump   (<base>_mat_initial.obj_initial_heap_state_.json)
  3. Python recompute    (qem_check.cost on the .off geometry)

This isolates whether a cost gap is (our C++ vs MeshLab) or (Python model vs C++).

Usage:
    py edge_lookup.py v0 v1 [v0 v1 ...]
    py edge_lookup.py --base <prefix> v0 v1 ...
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qem_check as qc
from heap_compare import load_meshlab_heap, DEFAULT_BASE


def main():
    args = sys.argv[1:]
    base = DEFAULT_BASE
    if args and args[0] == "--base":
        base = args[1]
        args = args[2:]
    ids = list(map(int, args))
    pairs = list(zip(ids[0::2], ids[1::2]))

    ours = load_meshlab_heap(base + "_initial_heap_state_.json")          # same parser
    ml = load_meshlab_heap(base + "_mat_initial.obj_initial_heap_state_.json")
    V, F = qc.load_off(base + "_mat_initial.off")
    Qd, vfaces, edgef, SF, diag = qc.build(V, F)

    print(f"{'edge':>12} {'faces':>5} {'C++heap':>13} {'MeshLab':>13} {'python':>13}"
          f" {'py/C++':>9} {'C++/ML':>9}")
    for a, b in pairs:
        k = (min(a, b), max(a, b))
        nf = len(edgef.get(k, []))
        c_cpp = ours.get(k)
        c_ml = ml.get(k)
        c_py = qc.cost(V, F, Qd, vfaces, edgef, SF, a, b, verbose=False)
        def fmt(x): return f"{x:.4e}" if isinstance(x, float) else str(x)
        py_cpp = (c_py / c_cpp) if c_cpp else float('nan')
        cpp_ml = (c_cpp / c_ml) if (c_cpp and c_ml and c_ml < 1e307) else float('nan')
        print(f"{a}-{b:<10} {nf:>5} {fmt(c_cpp):>13} {fmt(c_ml):>13} {fmt(c_py):>13}"
              f" {py_cpp:>9.2e} {cpp_ml:>9.2e}")


if __name__ == "__main__":
    main()
