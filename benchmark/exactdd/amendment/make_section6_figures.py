#!/usr/bin/env python3
"""Figures 6 to 8 of the paper (EVDD, float against algebraic, one panel per
family and quantity) remade from the EVDD runs of the 25-26 September sweep,
in the look of the originals, with the runs that did not finish added.

  make_section6_figures.py [--old SRC]... [--qasm DIR]... [NEWDIR]... OUTDIR

Writes OUTDIR/{Grover,w-state,random_circuits}/ with the file names of
figures/final_plots_mm-(in)correct/, so that the folder can stand in for that
one. The sources are those of make_figures.py: --old the September records
(exactdd/final.log, gw/out), NEWDIR the October LIMDD runs, which are read only
to judge a float run whose exact EVDD run did not finish.

What differs from the originals, and has to be said wherever these are used:
the float arm is 'low' normalisation (the September sweep has no 'L2' run), the
Grover and W-state circuits are those of ../gen_grover.py,
../gen_wstate_ancilla.py and ../gen_wstate_clifford.py (12 and 11, against the
99 and 7 points of the originals), and the runs are those of the Apple M1 sweep.
The correctness rule is make_figures.py's: a norm more than 1e-3 off 1, the
runner's error value, or a top-qubit probability more than 5% off the exact one.
The September Grover and W-state records carry no probability, so those float
runs are judged by their norm.

A run that did not finish (a timeout or a full table) sits on a dashed line:
on the right when the float run failed, on top when the algebraic one did, in
the corner when both did; its marker is grey when the float run is the one that
failed, since its correctness is then undefined. A run that was not made is not
drawn.
"""
import argparse, math, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, LogLocator, NullFormatter
import make_figures as mf

FLOAT, EXACT = "float_low", "qisq2_low"
PANELS = [("Grover", "grover"), ("w-state", "wstate"), ("random_circuits", "random")]
QUANT = [("wgt_type_runtime_float_vs_algebraic_nw.pdf", "simulation_time", "runtime (s)"),
         ("wgt_final_nodecount_float_vs_algebraic_nw.pdf", "final_nodes", "final # of nodes"),
         ("wgt_nodecount_float_vs_algebraic_nw.pdf", "max_nodes", "peak # of nodes")]
COL_OK, COL_BAD, COL_NA = "royalblue", "darkorange", "0.55"   # as the originals, plus grey

def panel(rows, fails, fam, field, lab, path):
    pts = {"ok": ([], []), "bad": ([], []), "na": ([], [])}
    nfail = {"float": 0, "algebraic": 0, "both": 0}
    circs = {c for c in rows if rows[c]["_fam"] == fam} | \
            {c for (c, a) in fails if a in (FLOAT, EXACT) and mf.family_of(c) == fam}
    for c in sorted(circs):
        r = rows.get(c, {"_fam": fam})
        fok, eok = FLOAT in r, EXACT in r
        ff, ef = (c, FLOAT) in fails, (c, EXACT) in fails
        if not (fok or ff) or not (eok or ef): continue          # an arm not run
        try:
            x = float(r[FLOAT][field]) if fok else None
            y = float(r[EXACT][field]) if eok else None
        except (KeyError, TypeError, ValueError): continue
        if x is None and y is None: nfail["both"] += 1
        elif x is None: nfail["float"] += 1
        elif y is None: nfail["algebraic"] += 1
        if not fok: k = "na"
        else: k = "bad" if mf.judge_point(r, FLOAT, EXACT)[0] else "ok"
        pts[k][0].append(x); pts[k][1].append(y)
    vals = [v for xs, ys in pts.values() for v in xs + ys if v is not None and v > 0]
    if not vals: print(f"  {path}: no data"); return
    lo, hi = min(vals) * 0.5, max(vals) * 2.0
    anyfail = any(nfail.values())
    gap = 10 ** (0.08 * math.log10(hi / lo))
    F = hi * gap
    top = F * gap if anyfail else hi
    floor = min(vals) * 0.6
    fig, ax = plt.subplots(figsize=(3.765, 2.824))
    ax.plot([lo, hi], [lo, hi], ls="--", color="0.45", lw=1.2, zorder=1)
    if anyfail:
        mf.band(ax, lo, F, gap, top, nfail["both"], 6.5)
    style = {"ok": ("^", COL_OK, "float measurement correct"),
             "bad": ("o", COL_BAD, "float measurement wrong"),
             "na": ("s", COL_NA, "float run did not finish")}
    for k in ("ok", "bad", "na"):
        xs, ys = pts[k]
        if not xs: continue
        m, col, label = style[k]
        ax.scatter([F if v is None else max(v, floor) for v in xs],
                   [F if v is None else max(v, floor) for v in ys],
                   marker=m, color=col, s=18, zorder=3, label=label)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, top); ax.set_ylim(lo, top)
    if anyfail:
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_locator(FixedLocator([t for t in LogLocator().tick_values(lo, hi) if lo <= t <= hi]))
            axis.set_minor_locator(FixedLocator(
                [t for t in LogLocator(subs=range(2, 10)).tick_values(lo, hi) if lo <= t <= hi]))
            axis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(f"{lab} float"); ax.set_ylabel(f"{lab} algebraic")
    edge = 1.0                                 # keep the legend out of the band
    if anyfail: edge = math.log(F / gap ** 0.5 / lo) / math.log(top / lo)
    ax.legend(loc="lower right", bbox_to_anchor=(edge, 0.0), prop={"family": "monospace", "size": 8.5})
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight"); plt.close(fig)
    print(f"  {os.path.relpath(path)}: {sum(len(v[0]) for v in pts.values())} points, "
          f"{len(pts['bad'][0])} wrong, did not finish: float {nfail['float']}, "
          f"algebraic {nfail['algebraic']}, both {nfail['both']}")

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("paths", nargs="+", help="[NEWDIR...] OUTDIR")
    ap.add_argument("--old", action="append", default=[])
    ap.add_argument("--qasm", action="append", default=[])
    a = ap.parse_args()
    out, new = a.paths[-1], a.paths[:-1]
    for d in a.qasm: mf.count_qasm(d)
    recs, rows, _ = mf.load(new, a.old, {FLOAT, EXACT})
    fails = mf.failures(recs)
    for sub, fam in PANELS:
        os.makedirs(os.path.join(out, sub), exist_ok=True)
        for fname, field, lab in QUANT:
            panel(rows, fails, fam, field, lab, os.path.join(out, sub, fname))

if __name__ == "__main__":
    main()
