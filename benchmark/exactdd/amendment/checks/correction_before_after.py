#!/usr/bin/env python3
"""Exact LIMDD (limdd_qisq2) on the random circuits up to 64 qubits, before
and after the correction for more than 64 qubits: the uncorrected pass
(ROOT/ldd/out) against the corrected one (ROOT/amend/out_rnd). Behind the
Setup sentence on the 29 random circuits that the exact LIMDD completes both
before and after the correction, and open item 5.

  correction_before_after.py ROOT      (or ROOT in the environment)
"""
import importlib.util
import os, sys, statistics as S
SP = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("ROOT")
if not SP: sys.exit("usage: correction_before_after.py ROOT")
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("mf", os.path.join(here, "..", "make_figures.py"))
mf = importlib.util.module_from_spec(spec); spec.loader.exec_module(mf)
mf.count_qasm(f"{SP}/ldd/qasm")
ARM = "limdd_qisq2"
def records(src):
    return {c: (s, st) for c, a, s, st, lab in mf.load_source(src)
            if a == ARM and mf.family_of(c) == "random" and (mf.nqubits(c, st) or 0) <= mf.WIDE}
before, after = records(f"{SP}/ldd/out"), records(f"{SP}/amend/out_rnd")
both, same, slower, ratio, unfinished = [], 0, 0, [], []
for c in sorted(before):
    s0, st0 = before[c]
    if s0 != "OK" or st0 is None:
        unfinished.append(f"{c[19:]}({s0}, after: {after[c][0] if c in after else 'not run'})"); continue
    if c not in after or after[c][0] != "OK" or after[c][1] is None: continue
    st1 = after[c][1]; both.append(c)
    nodes0, nodes1 = (st0["final_nodes"], st0["max_nodes"]), (st1["final_nodes"], st1["max_nodes"])
    q = st1["simulation_time"] / st0["simulation_time"]
    same += nodes0 == nodes1; slower += q > 1; ratio.append(q)
    print(f"{c[19:]:8s} final/peak {nodes0[0]:.0f}/{nodes0[1]:.0f} -> {nodes1[0]:.0f}/{nodes1[1]:.0f}"
          f"  time {st0['simulation_time']:.3f} -> {st1['simulation_time']:.3f} s  x{q:.3f}")
print(f"{ARM}, random circuits up to {mf.WIDE} qubits: {len(before)} before the correction, {len(after)} after")
print(f"  not finished before the correction ({len(unfinished)}): {unfinished}")
print(f"  finished in both passes: {len(both)}; same final and peak nodes on {same}; "
      f"corrected slower on {slower}; median time after/before {S.median(ratio) if ratio else float('nan'):.3f}")
