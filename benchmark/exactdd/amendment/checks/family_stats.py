#!/usr/bin/env python3
"""The per-family counts behind the text of the paper's Section 6 and the LIMDD
section, from a clean sweep alone: for each pair of arms (float against exact,
LIMDD against EVDD) and each family, the circuits both finish, the wrong float
runs (make_figures.judge: a measurement more than 5% off, or a norm more than
1e-3 off 1; 'meas' counts the first kind), the medians of exact/float time,
final and peak nodes, and the runs that did not finish.

  family_stats.py ROOT
ROOT as for regen_amendment.sh: ldd/qasm, amend/qasm, amend/out_clean, and
where present pilot/qasm/<family>/, pilot/out_clean, pilot/reference_probs.csv.
Run with the plotting Python (it loads make_figures.py)."""
import os, sys, statistics as st
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import make_figures as mf, glob
if len(sys.argv) != 2: sys.exit(__doc__)
S = sys.argv[1]
for d in [S + "/ldd/qasm", S + "/amend/qasm"] + glob.glob(S + "/pilot/qasm/*/"): mf.count_qasm(d)
if os.path.exists(S + "/pilot/reference_probs.csv"): mf.read_reference(S + "/pilot/reference_probs.csv")
recs, rows, _ = mf.load([d for d in (S + "/amend/out_clean", S + "/pilot/out_clean") if os.path.isdir(d)], [], set())
mf.add_reference(rows)
mf.UNCHECKED = "grey"
def ok(c, a): return a in rows.get(c, {})
def stat(c, a): r = recs.get((c, a)); return r["status"] if r else "-"
def cmp(fl, ex, label):
    print(f"== {label}: {fl} vs {ex}")
    fams = {}
    for c, r in rows.items(): fams.setdefault(r["_fam"], []).append(c)
    for fam in ["random", "grover", "wstate", "mirror", "iqp", "groveriter", "qwalk", "synth"]:
        cs = sorted(fams.get(fam, []), key=lambda c: mf.nqubits(c) or 0)
        both = [c for c in cs if ok(c, fl) and ok(c, ex)]
        wrong = [c for c in cs if ok(c, fl) and mf.judge(rows[c], fl, ex)[0] in ("meas", "norm")]
        meas = [c for c in cs if ok(c, fl) and mf.judge(rows[c], fl, ex)[0] == "meas"]
        tr = [float(rows[c][ex]["simulation_time"]) / max(float(rows[c][fl]["simulation_time"]), 1e-4) for c in both]
        nr = [float(rows[c][ex]["final_nodes"]) / max(float(rows[c][fl]["final_nodes"]), 1) for c in both]
        pr = [float(rows[c][ex]["max_nodes"]) / max(float(rows[c][fl]["max_nodes"]), 1) for c in both]
        fdnf = sum(1 for c in cs if stat(c, fl) not in ("OK", "-", "SKIPPED_PREDICTED"))
        ednf = sum(1 for c in cs if stat(c, ex) not in ("OK", "-", "SKIPPED_PREDICTED"))
        exfaster = sum(1 for x in tr if x < 1)
        md = lambda v: f"{st.median(v):.2f}" if v else "-"
        print(f"  {fam:10s} n={len(cs):2d} both={len(both):2d} wrong={len(wrong):2d} (meas {len(meas)}) "
              f"exact/float time med {md(tr)} exact faster {exfaster}/{len(both)}; final {md(nr)} peak {md(pr)}; "
              f"dnf float {fdnf} exact {ednf}")
        if fam in ("grover", "wstate"):
            for c in cs[-4:]:
                if ok(c, fl) and ok(c, ex):
                    print(f"      {c:22s} t {float(rows[c][fl]['simulation_time']):8.2f} {float(rows[c][ex]['simulation_time']):8.2f} "
                          f"final {rows[c][fl]['final_nodes']} {rows[c][ex]['final_nodes']} peak {rows[c][fl]['max_nodes']} {rows[c][ex]['max_nodes']} {mf.judge(rows[c], fl, ex)[0]}")
cmp("float_low", "qisq2_low", "Fig 6 low")
cmp("float_l2", "qisq2_low", "Fig 6 L2")
cmp("limdd_float", "limdd_qisq2", "Fig 7 low")
cmp("limdd_float_l2", "limdd_qisq2", "Fig 7 L2")
cmp("qisq2_low", "limdd_qisq2", "Fig 8 top: EVDD exact (as 'float') vs LIMDD exact")
print("== only exact completes / only float completes (LIMDD)")
for fl in ("limdd_float", "limdd_float_l2"):
    oe = [c for c in rows if stat(c, "limdd_qisq2") == "OK" and stat(c, fl) not in ("OK", "-", "SKIPPED_PREDICTED")]
    of = [(c, stat(c, "limdd_qisq2")) for c in rows if stat(c, fl) == "OK" and stat(c, "limdd_qisq2") not in ("OK", "-", "SKIPPED_PREDICTED")]
    ofw = [(c, s, mf.judge(rows[c], fl, "limdd_qisq2")[0]) for c, s in of]
    print(f"  {fl}: only exact finishes {len(oe)}; only float finishes {len(of)}: {ofw}")
print("== Fig 8 bottom: LIMDD float vs EVDD float, low and L2")
for fe, fl in (("float_low", "limdd_float"), ("float_l2", "limdd_float_l2")):
    both = [c for c in rows if ok(c, fe) and ok(c, fl)]
    lf = sum(1 for c in both if float(rows[c][fl]["simulation_time"]) < float(rows[c][fe]["simulation_time"]))
    ld = sum(1 for c in rows if stat(c, fl) not in ("OK", "-", "SKIPPED_PREDICTED") and stat(c, fe) == "OK")
    ed = sum(1 for c in rows if stat(c, fe) not in ("OK", "-", "SKIPPED_PREDICTED") and stat(c, fl) == "OK")
    wl = sum(1 for c in rows if ok(c, fl) and mf.judge(rows[c], fl, "limdd_qisq2")[0] in ("meas", "norm"))
    we = sum(1 for c in rows if ok(c, fe) and mf.judge(rows[c], fe, "qisq2_low")[0] in ("meas", "norm"))
    print(f"  {fl} vs {fe}: both {len(both)}, LIMDD faster {lf}; LIMDD dnf where EVDD ok {ld}; EVDD dnf where LIMDD ok {ed}; wrong LIMDD {wl} EVDD {we}")
print("== exact: LIMDD vs EVDD")
both = [c for c in rows if ok(c, "qisq2_low") and ok(c, "limdd_qisq2")]
print("  both", len(both), "LIMDD fewer peak", sum(1 for c in both if float(rows[c]["limdd_qisq2"]["max_nodes"]) < float(rows[c]["qisq2_low"]["max_nodes"])),
      "LIMDD faster", sum(1 for c in both if float(rows[c]["limdd_qisq2"]["simulation_time"]) < float(rows[c]["qisq2_low"]["simulation_time"])))
