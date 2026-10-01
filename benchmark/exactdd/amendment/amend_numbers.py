#!/usr/bin/env python3
"""Further numbers quoted in experiments_limdd_amendment.tex (revision of
2026-10-01), from the same records and correctness rule as make_figures.py.

  amend_numbers.py --old OLD... --qasm DIR... NEWDIR...   (as regen_amendment.sh calls it)

Prints: per-arm finished/correct counts on the random circuits, the upper-row
and lower-row ratios split by family group, the float/exact inflation inside
each diagram, the failure causes of the new LIMDD runs, and old-vs-new
agreement of the LIMDD arms on circuits up to 64 qubits.
"""
import argparse, glob, importlib.util, os, re, statistics as S
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("mf", os.path.join(here, "make_figures.py"))
mf = importlib.util.module_from_spec(spec); spec.loader.exec_module(mf)

ap = argparse.ArgumentParser()
ap.add_argument("new", nargs="+"); ap.add_argument("--old", action="append", default=[])
ap.add_argument("--qasm", action="append", default=[])
a = ap.parse_args()
for d in a.qasm: mf.count_qasm(d)
recs, rows, _ = mf.load(a.new, a.old, {"qisq2_low", "float_low"})
E, L, EF, LF = "qisq2_low", "limdd_qisq2", "float_low", "limdd_float"
short = lambda c: c[19:] if c.startswith("clifford_T_circuit_") else c
med = lambda v: S.median(v) if v else float("nan")
def val(c, arm, k): return float(rows[c][arm][k])
def ok(c, arm): return c in rows and arm in rows[c]
def bad(c, arm, ref): return mf.judge_point(rows[c], arm, ref)[0]

allc = sorted(set(k[0] for k in recs))
rnd = [c for c in allc if mf.family_of(c) == "random"]
struct = [c for c in allc if mf.family_of(c) in ("grover", "wstate", "hshift")]

print("== random circuits (49): finished / correct")
for arm, ref in ((E, None), (L, None), (EF, E), (LF, L)):
    fin = [c for c in rnd if ok(c, arm)]
    wr = [short(c) for c in fin if ref and bad(c, arm, ref)]
    print(f"  {arm:12s} finished {len(fin)}, wrong {len(wr)} {sorted(wr)}")
print("  exact LIMDD only:", sorted(short(c) for c in rnd if ok(c, L) and not ok(c, E)))
print("  exact EVDD only:", sorted(short(c) for c in rnd if ok(c, E) and not ok(c, L)))
print("  neither exact:", sorted(short(c) for c in rnd if not ok(c, E) and not ok(c, L)))

def ratios(cs, x, y, k):
    return [(val(c, x, k) / val(c, y, k), c) for c in cs if ok(c, x) and ok(c, y)
            and val(c, x, k) > 0 and val(c, y, k) > 0]

def show(label, cs, x, y):
    print(f"\n== {label}: x={x} y={y}")
    for k in ("final_nodes", "max_nodes", "simulation_time"):
        r = sorted(ratios(cs, x, y, k))
        if not r: print("  no data"); continue
        print(f"  {k:16s} n={len(r)} x>y on {sum(q > 1 for q, _ in r)} (eq {sum(q == 1 for q, _ in r)}),"
              f" median x/y {med([q for q, _ in r]):.4g}, min {r[0][0]:.4g} ({short(r[0][1])}),"
              f" max {r[-1][0]:.4g} ({short(r[-1][1])})")

show("upper row, random", rnd, E, L)
show("upper row, structured", struct, E, L)
for fam in ("grover", "wstate", "hshift"):
    show(f"upper row, {fam}", [c for c in struct if mf.family_of(c) == fam], E, L)
show("lower row, random", rnd, EF, LF)
show("lower row, structured", struct, EF, LF)
four = [c for c in rnd if all(ok(c, k) for k in (E, L, EF, LF))]
show(f"all four finish, random ({len(four)})", four, E, L)
show(f"all four finish, random ({len(four)}), float", four, EF, LF)
show("float/exact inside LIMDD, random", rnd, LF, L)
show("float/exact inside EVDD, random", rnd, EF, E)
show("float/exact inside LIMDD, structured", struct, LF, L)
show("float/exact inside EVDD, structured", struct, EF, E)

print("\n== lower-row points (both floats finish), wrong per diagram")
pts = [c for c in allc if ok(c, EF) and ok(c, LF)]
both_w = [short(c) for c in pts if bad(c, LF, L) and bad(c, EF, E)]
only_l = [short(c) for c in pts if bad(c, LF, L) and not bad(c, EF, E)]
only_e = [short(c) for c in pts if bad(c, EF, E) and not bad(c, LF, L)]
print(f"  {len(pts)} points; both wrong {both_w}\n  only LIMDD wrong {only_l}\n  only EVDD wrong {only_e}")
r = sorted(ratios(pts, EF, LF, "final_nodes"))
for q, c in r[-3:] + r[:2]:
    print(f"  final ratio {q:.4g} {short(c)}: EVDD {val(c, EF, 'final_nodes'):.0f} LIMDD {val(c, LF, 'final_nodes'):.0f}"
          f" LIMDD wrong={bad(c, LF, L)} EVDD wrong={bad(c, EF, E)}")

print("\n== failure causes of the new LIMDD runs (status, and the table named in .out)")
for arm in (L, LF):
    for fams, lab in ((rnd, "random"), (struct, "structured")):
        st = {}
        for c in fams:
            rec = recs.get((c, arm))
            if rec is None: s = "missing"
            else: s = rec["status"]
            if s == "TAB_FULL":
                for d in a.new:
                    f = os.path.join(d, f"{c}_{arm}.out")
                    if os.path.exists(f):
                        t = open(f).read()
                        s += "(weights)" if "Amplitude table full" in t else "(nodes)" if "node table" in t.lower() else ""
            if s != "OK": st.setdefault(s, []).append(short(c))
        print(f"  {arm} {lab}: " + "; ".join(f"{k} {len(v)}: {' '.join(sorted(v))}" for k, v in sorted(st.items())))

print("\n== old vs new LIMDD, circuits up to 64 qubits")
old = {}
for src in a.old:
    for circ, arm, status, st, lab in mf.load_source(src):
        old[(circ, arm)] = (status, st)
for arm in (L, LF):
    same, diff, onlynew, onlyold, speed = [], [], [], [], []
    for c in allc:
        if (mf.nqubits(c) or 0) > 64 or mf.family_of(c) is None and not c.startswith(("adder", "ising")): continue
        o = old.get((c, arm)); n = recs.get((c, arm))
        if o is None or n is None: continue
        ofin = o[0] == "OK" and o[1] is not None; nfin = n["status"] == "OK" and n["st"] is not None
        if ofin and nfin:
            ks = ("final_nodes", "max_nodes")
            eq = all(float(o[1][k]) == float(n["st"][k]) for k in ks) and abs(float(o[1]["norm"]) - float(n["st"]["norm"])) < 1e-6
            (same if eq else diff).append(short(c))
            if float(n["st"]["simulation_time"]) > 0:
                speed.append(float(o[1]["simulation_time"]) / float(n["st"]["simulation_time"]))
        elif nfin: onlynew.append(short(c))
        elif ofin: onlyold.append(short(c))
    print(f"  {arm}: both finished {len(same) + len(diff)}, same nodes+norm {len(same)}, differ {len(diff)} {diff}")
    print(f"     finished only new {len(onlynew)} {onlynew}; only old {len(onlyold)} {onlyold}")
    if speed: print(f"     old/new time: median {med(speed):.3g}, min {min(speed):.3g}, max {max(speed):.4g}, n={len(speed)}")
