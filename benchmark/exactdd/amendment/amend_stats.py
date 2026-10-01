#!/usr/bin/env python3
"""Numbers quoted in experiments_limdd_amendment.tex, from the same merged
records and the same correctness rule as make_figures.py.

  amend_stats.py  (same --old/--qasm/NEWDIR arguments as make_figures.py, no OUTDIR)
"""
import importlib.util, os, sys, statistics as S
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("mf", os.path.join(here, "make_figures.py"))
mf = importlib.util.module_from_spec(spec); spec.loader.exec_module(mf)
import argparse
ap = argparse.ArgumentParser()
ap.add_argument("new", nargs="+"); ap.add_argument("--old", action="append", default=[])
ap.add_argument("--qasm", action="append", default=[])
a = ap.parse_args()
for d in a.qasm: mf.count_qasm(d)
recs, rows, _ = mf.load(a.new, a.old, {"qisq2_low", "float_low"})
E, L, EF, LF = "qisq2_low", "limdd_qisq2", "float_low", "limdd_float"
short = lambda c: c[19:] if c.startswith("clifford_T_circuit_") else c
def bad(r, arm, ref): return mf.judge_point(r, arm, ref)[0]
def both(r, x, y): return x in r and y in r
def med(v): return S.median(v) if v else float("nan")
def f(r, arm, k): return float(r[arm][k])

def compare(x, y, label, fams=None):
    pts = [(c, r) for c, r in rows.items() if both(r, x, y) and (fams is None or r["_fam"] in fams)]
    print(f"\n== {label}: {len(pts)} circuits where both finish" + (f" (families {fams})" if fams else ""))
    for k in ("final_nodes", "max_nodes", "simulation_time"):
        rat = [(f(r, x, k) / f(r, y, k), c) for c, r in pts if f(r, x, k) > 0 and f(r, y, k) > 0]
        smaller = sum(1 for q, _ in rat if q > 1); equal = sum(1 for q, _ in rat if q == 1)
        rs = sorted(rat)
        print(f"  {k:16s} y<x on {smaller}/{len(rat)} (equal {equal}); median x/y {med([q for q,_ in rat]):.3g}; "
              f"max {rs[-1][0]:.4g} ({short(rs[-1][1])}), min {rs[0][0]:.3g} ({short(rs[0][1])})")
    return pts

print("circuits:", len(rows))
compare(E, L, "upper row, exact EVDD (x) vs exact LIMDD (y)")
compare(E, L, "upper row, random only", {"random"})
pts = compare(EF, LF, "lower row, float EVDD (x) vs float LIMDD (y)")
compare(EF, LF, "lower row, random only", {"random"})
four = {c: r for c, r in rows.items() if all(k in r for k in (E, L, EF, LF))}
for fams in (None, {"random"}):
    sub = {c: r for c, r in four.items() if fams is None or r["_fam"] in fams}
    print(f"\n== all four finish: {len(sub)} circuits" + (" (random)" if fams else ""))
    for k in ("final_nodes", "max_nodes"):
        ex = [f(r, E, k) / f(r, L, k) for r in sub.values()]
        fl = [f(r, EF, k) / f(r, LF, k) for r in sub.values() if f(r, LF, k) > 0]
        print(f"  {k}: exact median EVDD/LIMDD {med(ex):.3g}, LIMDD smaller {sum(q>1 for q in ex)}/{len(ex)};"
              f" float median {med(fl):.3g}, LIMDD smaller {sum(q>1 for q in fl)}/{len(fl)}, equal {sum(q==1 for q in fl)}")
    # float / exact inside each diagram
    for nm, x, ex_ in (("LIMDD", LF, L), ("EVDD", EF, E)):
        for k in ("final_nodes", "max_nodes", "simulation_time"):
            q = [f(r, x, k) / f(r, ex_, k) for r in sub.values() if f(r, ex_, k) > 0]
            print(f"  {nm} float/exact {k}: median {med(q):.3g}, float larger/slower on {sum(v>1 for v in q)}/{len(q)}, max {max(q):.4g}")

print("\n== wrong float runs (make_figures rule) among circuits where float and an exact arm finish")
for arm, ref in ((LF, L), (EF, E)):
    w = sorted(short(c) for c, r in rows.items() if arm in r and (E in r or L in r) and bad(r, arm, ref))
    n = sum(1 for c, r in rows.items() if arm in r and (E in r or L in r))
    print(f"  {arm}: {len(w)}/{n}: {' '.join(w)}")
print("  lower-row points with both floats wrong / only LIMDD / only EVDD:")
bl, ol, oe = [], [], []
for c, r in pts:
    bL, bE = bad(r, LF, L), bad(r, EF, E)
    (bl if bL and bE else ol if bL else oe if bE else []).append(short(c))
print("   both", bl, "\n   only LIMDD", ol, "\n   only EVDD", oe)

print("\n== random circuits: finished / correct per arm (of 49)")
rnd = sorted(c for c in set(k[0] for k in recs) if mf.family_of(c) == "random")   # sorted: a set's order varies run to run
for arm in (E, L, EF, LF):
    fin = [c for c in rnd if (c, arm) in recs and recs[(c, arm)]["status"] == "OK"]
    if arm in (EF, LF):
        ref = L if arm == LF else E
        okc = [c for c in fin if not bad(rows[c], arm, ref)]
    else: okc = fin
    st = {}
    for c in rnd:
        s = recs[(c, arm)]["status"] if (c, arm) in recs else "missing"
        st[s] = st.get(s, 0) + 1
    print(f"  {arm:12s} finished {len(fin)}, correct {len(okc)}; statuses {st}")
ex_only = [short(c) for c in rnd if (c, L) in recs and recs[(c, L)]["status"] == "OK"
           and not ((c, E) in recs and recs[(c, E)]["status"] == "OK")]
ev_only = [short(c) for c in rnd if (c, E) in recs and recs[(c, E)]["status"] == "OK"
           and not ((c, L) in recs and recs[(c, L)]["status"] == "OK")]
print("  exact LIMDD finishes, exact EVDD not:", sorted(ex_only))
print("  exact EVDD finishes, exact LIMDD not:", sorted(ev_only))
print("  exact LIMDD fails:", sorted((short(c), recs[(c, L)]["status"]) for c in rnd if recs[(c, L)]["status"] != "OK"))
lf_fail = sorted((short(c), recs[(c, LF)]["status"]) for c in rnd if recs[(c, LF)]["status"] != "OK")
print("  float LIMDD fails:", lf_fail)
print("  float LIMDD fails where exact LIMDD finishes:", sum(1 for c in rnd if recs[(c, LF)]["status"] != "OK" and recs[(c, L)]["status"] == "OK"))
