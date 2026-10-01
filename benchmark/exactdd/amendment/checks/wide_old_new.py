#!/usr/bin/env python3
"""LIMDD cells above 64 qubits: the Sep sweep (final.log; gw/out for the other families)
against the corrected runs (ldd/out for the random circuits, amend/out for the others).

  wide_old_new.py ROOT      (or ROOT in the environment)
"""
import os, sys
import importlib.util
SP = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("ROOT")
if not SP: sys.exit("usage: wide_old_new.py ROOT")
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("mf", os.path.join(here, "..", "make_figures.py"))
mf = importlib.util.module_from_spec(spec); spec.loader.exec_module(mf)
for d in (f"{SP}/ldd/qasm", f"{SP}/amend/qasm"): mf.count_qasm(d)
old = {}
for src in (f"{SP}/exactdd/final.log", f"{SP}/gw/out", f"{SP}/gw/out2"):
    for c, a, s, st, lab in mf.load_source(src): old[(c, a)] = (s, st, lab)
new = {}
for src in (f"{SP}/ldd/out", f"{SP}/amend/out_rnd", f"{SP}/amend/out"):
    for c, a, s, st, lab in mf.load_source(src): new[(c, a)] = (s, st, lab)
def fmt(r):
    if r is None: return "missing"
    s, st, lab = r
    if st is None: return f"{s} [{lab}]"
    return (f"{s} {st['final_nodes']:.0f}/{st['max_nodes']:.0f} norm {st['norm']:.6g} "
            f"p {st.get('first_qubit_measurement_prob')} t {st['simulation_time']:.1f} [{lab}]")
changed = {a: [] for a in mf.LIMDD_ARMS}
for c in sorted(set(k[0] for k in list(old) + list(new))):
    if (mf.nqubits(c) or 0) <= 64: continue
    for a in mf.LIMDD_ARMS:
        o, n = old.get((c, a)), new.get((c, a))
        same = (o and n and o[1] and n[1] and o[0] == n[0] == "OK"
                and (o[1]["final_nodes"], o[1]["max_nodes"]) == (n[1]["final_nodes"], n[1]["max_nodes"])
                and abs(o[1]["norm"] - n[1]["norm"]) < 1e-6)
        tag = "same" if same else "CHANGED"
        if not same and o and o[0] == "OK" and n and n[0] == "OK": changed[a].append(c)
        print(f"{c:28s} {mf.nqubits(c):3d} {a:12s} {tag:8s} old {fmt(o)}\n{'':55s} new {fmt(n)}")
for a, v in changed.items(): print(a, "node counts or norm changed (both OK):", len(v), v)
