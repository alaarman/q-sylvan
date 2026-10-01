#!/usr/bin/env python3
"""The 38 float EVDD runs on the random circuits under two rules:
   norm rule (|norm-1| > 1e-3, or sentinel, or 5% off the exact probability) as in make_figures.py,
   Section 6 rule (5% off the exact probability) plus the sentinel check, no norm.
Same for the float LIMDD runs.

  float_evdd_rules.py ROOT      (or ROOT in the environment)
"""
import importlib.util
import os, sys
SP = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("ROOT")
if not SP: sys.exit("usage: float_evdd_rules.py ROOT")
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("mf", os.path.join(here, "..", "make_figures.py"))
mf = importlib.util.module_from_spec(spec); spec.loader.exec_module(mf)
for d in (f"{SP}/ldd/qasm", f"{SP}/amend/qasm"): mf.count_qasm(d)
recs, rows, _ = mf.load([f"{SP}/ldd/out", f"{SP}/amend/out_rnd", f"{SP}/amend/out"],
                        [f"{SP}/exactdd/final.log", f"{SP}/gw/out", f"{SP}/gw/out2"], {"qisq2_low", "float_low"})
short = lambda c: c[19:]
for arm, own in (("float_low", "qisq2_low"), ("limdd_float", "limdd_qisq2")):
    cs = sorted(c for c, r in rows.items() if r["_fam"] == "random" and arm in r)
    norm_w, sec6_w, sent, noref, p_ok_norm_bad = [], [], [], [], []
    for c in cs:
        r = rows[c]; d = r[arm]
        bad, how = mf.judge_point(r, arm, own)
        if bad: norm_w.append(short(c))
        pf = mf.p(d)
        ex = None
        for a in (own,) + tuple(e for e in mf.EXACT if e != own):
            if a in r and mf.p(r[a]) is not None: ex = mf.p(r[a]); break
        if pf is not None and pf >= mf.SENTINEL:
            sec6_w.append(short(c)); sent.append(short(c)); continue
        if ex is None:
            noref.append(f"{short(c)}(norm {d['norm']:.3g}, p {pf})"); continue
        if abs(pf - ex) > 0.05 * max(abs(ex), 1e-12): sec6_w.append(short(c))
        elif bad: p_ok_norm_bad.append(f"{short(c)}(norm {d['norm']:.4g}, p {pf} vs {ex})")
    print(f"{arm}: {len(cs)} finished")
    print(f"  norm rule wrong {len(norm_w)}: {norm_w}")
    print(f"  Section 6 rule + sentinel wrong {len(sec6_w)}: {sec6_w}  (sentinel {sent})")
    print(f"  no exact probability to compare with ({len(noref)}): {noref}")
    print(f"  probability right but norm off ({len(p_ok_norm_bad)}): {p_ok_norm_bad}")
