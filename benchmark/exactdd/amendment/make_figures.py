#!/usr/bin/env python3
"""float vs algebraic, and LIMDD vs EVDD, pooled over benchmark families.

Colour means what it means in the paper's own figures -- whether the float
measurement is correct -- and the family is carried by the marker SHAPE, so
nothing has to be relearned to read these beside Figures 7 to 9.

The adder and the Ising chain are deliberately not plotted. The adder puts
every point on the diagonal (float and algebraic agree at every size), and the
Ising chain saturates the state space by 16 qubits; both are reported in the
text instead of diluting the panels.

Data sources (added on Oct 1 to an earlier version of this script):

  make_figures.py [--old SRC]... [--old-arms LIST] [--qasm DIR]...
                  [--dump CSV] NEWDIR... OUTDIR

  NEWDIR   a run directory of <circ>_<arm>.json and/or <circ>_<arm>.status
           files (run_arm.sh in this directory, or the earlier version of it
           that wrote ldd/out, whose status lines lack the last two columns).
           Its records override every --old record of the same circuit and
           arm, whatever their status: a new TIMEOUT removes an old point
           rather than letting the old one show through.
  --old    a surviving record of the Sep 25-26 sweeps: a final.log-style file
           ("arm circuit status time final peak norm p_top", later lines
           override earlier ones, so the phase 2 and 3 retries win), or a
           directory of .status files in either the run_arm.sh layout (same
           columns, optionally followed by t_count final_width) or the
           layout of the September run_gw.sh, which is not kept here
           ("circuit arm status time final peak norm t_count final_width",
           which has NO top-qubit probability), or a directory of JSON files.
  --old-arms  which arms the old records may supply (default: the EVDD arms
           qisq2_low,float_low; "all" for every arm). Old LIMDD records on
           more than 64 qubits are dropped whatever this says: those runs used
           64-bit qubit masks and simulated a wrong state.

With only positional directories the script behaves as that earlier version.

Correctness rule: a float run is wrong when its norm is off 1 by more than
1e-3, or its top-qubit probability is the error sentinel, or it differs from
the exact one by more than 5%. The exact probability is taken from the panel's
own exact arm, and from the other exact arm when that one did not finish (the
two are exact, so they agree wherever both finished; the script checks this).
When no exact arm finished, or the float record carries no probability (the
run_gw.sh status lines), the run is judged by its norm alone, and the
printout says how many points that applies to.
"""
import argparse, json, os, re, sys, glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT  = None                                    # set in main
VERBOSE = False
ARMS = ["qisq2_low", "limdd_qisq2", "float_low", "limdd_float"]
LIMDD_ARMS = ("limdd_qisq2", "limdd_float")
EXACT = ("limdd_qisq2", "qisq2_low")           # either can serve as the truth
SENTINEL = 1e10
WIDE = 64                                      # one 64-bit Pauli word
NUMERIC = ("simulation_time", "final_nodes", "max_nodes", "norm",
           "first_qubit_measurement_prob", "t_count", "final_width")

# family -> (marker, label); order fixes the legend
FAMILY = [
    ("random",  "o", "random Clifford+$T$"),
    ("grover",  "s", "Grover"),
    ("wstate",  "^", "W state"),
    ("hshift",  "D", "hidden shift"),
]
# AMEND_FAMILIES=random,grover,wstate keeps only those families, in the figures,
# the legend and every number the other scripts take from family_of. Unset, all
# four are kept.
ONLY = [f for f in os.environ.get("AMEND_FAMILIES", "").split(",") if f]
if ONLY:
    unknown = sorted(set(ONLY) - {k for k, _, _ in FAMILY})
    if unknown: sys.exit(f"AMEND_FAMILIES: unknown families {unknown}")
    FAMILY = [f for f in FAMILY if f[0] in ONLY]
MARK = {k: m for k, m, _ in FAMILY}

def family_of(name):
    if name.startswith("clifford_T_circuit"): fam = "random"
    elif name.startswith("grover"):           fam = "grover"
    elif name.startswith(("w-state", "wclif")): fam = "wstate"
    elif name.startswith("hidden-shift"):     fam = "hshift"
    else: return None                          # adder, ising: text only
    return fam if fam in MARK else None

# --------------------------------------------------------------------------
# qubit counts, needed to drop the old LIMDD runs past one Pauli word
# --------------------------------------------------------------------------
QUBITS = {}                                    # circuit -> total qubits, from --qasm

def count_qasm(d):
    for f in glob.glob(os.path.join(d, "*.qasm")):
        n = 0
        for m in re.finditer(r"qreg\s+\w+\s*\[\s*(\d+)\s*\]", open(f).read()):
            n += int(m.group(1))
        if n: QUBITS[os.path.basename(f)[:-5]] = n

def nqubits(circ, st=None):
    if circ in QUBITS: return QUBITS[circ]
    if st and st.get("n_qubits"):
        try: return int(st["n_qubits"])
        except (TypeError, ValueError): pass
    m = re.match(r"clifford_T_circuit_(\d+)_\d+$", circ)
    if m: return int(m.group(1))
    m = re.match(r"grover_n(\d+)_it\d+$", circ)
    if m: return 2 * int(m.group(1)) - 2       # n data + n-2 work qubits (checked on the qasm, n=4..50)
    m = re.search(r"_n(\d+)", circ)            # w-state_K_nM, wclif_K_nM, hidden-shift_nM, adder_nM, ising_nM_s3
    return int(m.group(1)) if m else None

# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------
def split_name(b):
    for a in sorted(ARMS, key=len, reverse=True):
        if b.endswith("_" + a): return b[:-(len(a) + 1)], a
    return None, None

def numeric(st):
    for k in NUMERIC:
        if k in st:
            try: st[k] = float(st[k])
            except (TypeError, ValueError): st[k] = None
    return st

RUNARM_COLS = ["simulation_time", "final_nodes", "max_nodes", "norm",
               "first_qubit_measurement_prob"]
RUNARM_COLS_T = RUNARM_COLS + ["t_count", "final_width"]   # run_arm.sh (Oct 1)
GW_COLS     = ["simulation_time", "final_nodes", "max_nodes", "norm",
               "t_count", "final_width"]

def parse_status(line, where):
    """One status line, in the run_arm.sh or the run_gw.sh layout."""
    t = line.split()
    if len(t) < 3: return None
    if   t[0] in ARMS: arm, circ, cols = t[0], t[1], RUNARM_COLS
    elif t[1] in ARMS: circ, arm, cols = t[0], t[1], GW_COLS
    else: return None
    status = t[2]
    if status != "OK": return circ, arm, status, None
    vals = t[3:]
    if cols is RUNARM_COLS and len(vals) == len(RUNARM_COLS_T): cols = RUNARM_COLS_T
    if len(vals) != len(cols):
        print(f"  WARNING {where}: {len(vals)} fields where {len(cols)} were expected: {line.strip()}")
        return circ, arm, "BADLINE", None
    return circ, arm, status, numeric(dict(zip(cols, vals)))

def read_json(f):
    try: return numeric(dict(json.load(open(f))["statistics"]))
    except Exception: return None

def load_source(path):
    """[(circ, arm, status, stats or None, label)] in the order they override."""
    out = []
    if os.path.isfile(path):
        lab = os.path.basename(path)
        for line in open(path):
            r = parse_status(line, lab)
            if r: out.append(r + (lab,))
        return out
    lab = os.path.basename(os.path.normpath(path))
    par = os.path.basename(os.path.dirname(os.path.normpath(path)))
    lab = f"{par}/{lab}" if par else lab
    seen = {}
    for f in sorted(glob.glob(os.path.join(path, "*.status"))):
        txt = open(f).read().strip()
        r = parse_status(txt, os.path.basename(f)) if txt else None
        if r is None:
            circ, arm = split_name(os.path.basename(f)[:-7])
            if circ is None: continue
            r = (circ, arm, "EMPTY", None)
        seen[(r[0], r[1])] = list(r)
    for f in sorted(glob.glob(os.path.join(path, "*.json"))):
        circ, arm = split_name(os.path.basename(f)[:-5])
        if circ is None: continue
        st = read_json(f)
        k = (circ, arm)
        if k in seen:
            if seen[k][2] == "OK" and st is not None: seen[k][3] = st   # full precision
        elif st is not None:
            seen[k] = [circ, arm, "OK", st]                             # JSON without status
    return [tuple(v) + (lab,) for v in seen.values()]

def load(new_srcs, old_srcs, old_arms):
    recs, dropped = {}, {"arm": 0, "wide": []}
    for src in old_srcs:
        for circ, arm, status, st, lab in load_source(src):
            if arm not in old_arms: dropped["arm"] += 1; continue
            if arm in LIMDD_ARMS and (nqubits(circ, st) or 0) > WIDE:
                dropped["wide"].append(f"{circ}/{arm}"); continue
            recs[(circ, arm)] = dict(status=status, st=st, src=lab, old=True)
    for src in new_srcs:
        for circ, arm, status, st, lab in load_source(src):
            recs[(circ, arm)] = dict(status=status, st=st, src=lab, old=False)
    rows = {}
    for (circ, arm), rec in recs.items():
        fam = family_of(circ)
        if fam is None or rec["status"] != "OK" or rec["st"] is None: continue
        r = rows.setdefault(circ, {"_fam": fam, "_src": {}})
        r[arm] = rec["st"]; r["_src"][arm] = rec["src"]
    return recs, rows, dropped

# --------------------------------------------------------------------------
# correctness
# --------------------------------------------------------------------------
def p(d):
    v = d.get("first_qubit_measurement_prob")
    return None if v is None else float(v)

def norm_of(d):
    try: return float(d.get("norm"))
    except (TypeError, ValueError): return None

def wrong(pf, pa, dfloat=None):
    """The paper's 5% rule on the measurement, AND a check on the norm.

    The 5% rule alone is not enough: on the 30-qubit hidden shift the float
    run returns norm 0 -- the state has been annihilated -- while its
    top-qubit probability coincidentally equals the exact one, so the rule
    calls it correct. A unitary circuit cannot produce an unnormalised state,
    so a norm away from 1 is a wrong answer whatever the measurement says.
    """
    if dfloat is not None:
        nf = norm_of(dfloat)
        if nf is None or nf != nf or abs(nf - 1.0) > 1e-3: return True
    if pf is None or pa is None: return True
    if pf >= SENTINEL or pf != pf: return True
    return abs(pf - pa) > 0.05 * max(abs(pa), 1e-12)

def judge_point(r, judge, ref):
    """(bad, how). how is 'norm!=1' (wrong on the norm alone), 'sentinel'
    (the runner's error value for the probability), 'p' (5% rule against the
    panel's own exact arm), 'p:<arm>' (against the other exact arm, the own
    one did not finish), or 'norm-only' (norm is 1 and there is nothing to
    compare the probability with: no exact arm finished, or the float record
    is a run_gw.sh status line without a probability; counted as correct)."""
    d = r[judge]
    nf = norm_of(d)
    if nf is None or nf != nf or abs(nf - 1.0) > 1e-3: return True, "norm!=1"
    pf = p(d)
    if pf is not None and (pf >= SENTINEL or pf != pf): return True, "sentinel"
    if pf is not None:
        for a in (ref,) + tuple(e for e in EXACT if e != ref):
            if a in r and p(r[a]) is not None:
                return wrong(pf, p(r[a]), d), ("p" if a == ref else "p:" + a)
    return False, "norm-only"

# --------------------------------------------------------------------------
# plotting (unchanged look)
# --------------------------------------------------------------------------
def scatter(rows, xarm, yarm, field, xlabel, ylabel, fname, judge=None, truth=None):
    pts = {}          # (family, bad) -> ([x],[y])
    nfloor = 0
    how = {}          # how each judged point was judged
    perfam = {}       # family -> [points, wrong]
    detail = []       # (circuit, how) of every point drawn red, for -v
    for circ, r in rows.items():
        if xarm not in r or yarm not in r: continue
        try: x, y = float(r[xarm][field]), float(r[yarm][field])
        except (KeyError, TypeError, ValueError): continue
        ref = truth if truth else yarm
        bad, h = judge_point(r, judge, ref) if judge else (False, None)
        if x <= 0 or y <= 0:
            if not judge: continue
            bad = True; nfloor += 1; h = (h or "") + "+floored"
        if h: how[h] = how.get(h, 0) + 1
        if bad: detail.append(f"{circ[19:] if circ.startswith('clifford_T_circuit_') else circ}:{h}")
        pf = perfam.setdefault(r["_fam"], [0, 0]); pf[0] += 1; pf[1] += bad
        pts.setdefault((r["_fam"], bad), ([], []))
        pts[(r["_fam"], bad)][0].append(x); pts[(r["_fam"], bad)][1].append(y)
    allv = [v for (xs, ys) in pts.values() for v in xs + ys if v > 0]
    if not allv: print(f"  {fname}: no data"); return
    floor = min(allv) * 0.35
    lo, hi = min(allv) * 0.5, max(allv) * 2.0

    fig, ax = plt.subplots(figsize=(2.9, 2.6))
    ax.plot([lo, hi], [lo, hi], color="0.6", lw=0.8, zorder=1)
    for (fam, bad), (xs, ys) in sorted(pts.items()):
        xs = [v if v > 0 else floor for v in xs]
        ys = [v if v > 0 else floor for v in ys]
        ax.scatter(xs, ys, s=20, marker=MARK[fam],
                   facecolors="none" if not bad else "tab:red",
                   edgecolors="tab:blue" if not bad else "tab:red",
                   linewidths=0.9, zorder=3 + bad)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel(xlabel, fontsize=7.5, labelpad=1.5)
    ax.set_ylabel(ylabel, fontsize=7.5, labelpad=1.5)
    ax.tick_params(labelsize=6.5, pad=1.5)
    fig.tight_layout(pad=0.25)
    fig.savefig(os.path.join(OUT, fname), bbox_inches="tight"); plt.close(fig)
    n = sum(len(v[0]) for v in pts.values())
    nb = sum(len(v[0]) for k, v in pts.items() if k[1])
    fams = " ".join(f"{k}={v[0]}/{v[1]}w" for k, v in sorted(perfam.items()))
    hw = ", judged by " + " ".join(f"{k}={v}" for k, v in sorted(how.items())) if how else ""
    print(f"  {fname}: {n} points ({nb} wrong{', %d floored' % nfloor if nfloor else ''})"
          f"  [{fams}{hw}]")
    if VERBOSE and detail: print("      red: " + " ".join(sorted(detail)))

def legend(fname):
    """One shared legend, so the panels keep their space."""
    fig = plt.figure(figsize=(6.4, 0.42))
    hs  = [plt.Line2D([], [], ls="", marker=m, mfc="none", mec="0.25", mew=0.9,
                      ms=5, label=lab) for _, m, lab in FAMILY]
    hs += [plt.Line2D([], [], ls="", marker="o", mfc="none", mec="tab:blue",
                      mew=0.9, ms=5, label="measurement correct"),
           plt.Line2D([], [], ls="", marker="o", mfc="tab:red", mec="tab:red",
                      ms=5, label="float measurement wrong")]
    fig.legend(handles=hs, loc="center", ncol=6, fontsize=7, frameon=False,
               handletextpad=0.35, columnspacing=1.1)
    fig.savefig(os.path.join(OUT, fname), bbox_inches="tight"); plt.close(fig)
    print("  wrote", fname)

# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------
def provenance(recs, dropped):
    print("records per arm (status OK / all), by source:")
    for a in ARMS:
        by = {}
        for (c, arm), rec in recs.items():
            if arm != a: continue
            k = ("old " if rec["old"] else "new ") + rec["src"]
            v = by.setdefault(k, [0, 0]); v[1] += 1; v[0] += rec["status"] == "OK"
        print(f"  {a:<12} " + "; ".join(f"{k}: {v[0]}/{v[1]}" for k, v in sorted(by.items())))
    if dropped["arm"]: print(f"  ignored {dropped['arm']} old records of arms outside --old-arms")
    if dropped["wide"]:
        w = sorted(set(dropped["wide"]))
        print(f"  dropped {len(w)} old LIMDD cells on >{WIDE} qubits (64-bit masks, wrong state): "
              + " ".join(w))

def exact_agreement(rows):
    """The two exact arms must agree; say so, or say where they do not."""
    n, bad = 0, []
    for c, r in sorted(rows.items()):
        if all(a in r and p(r[a]) is not None for a in EXACT):
            n += 1
            pa, pb = p(r[EXACT[0]]), p(r[EXACT[1]])
            if abs(pa - pb) > 1e-6 * max(1.0, abs(pb)): bad.append(f"{c} ({pa} vs {pb})")
    print(f"exact arms agree on the top-qubit probability on {n - len(bad)}/{n} circuits"
          + ("" if not bad else ": DISAGREE on " + ", ".join(bad)))

def dump(recs, path):
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["circuit", "family", "qubits", "arm", "status", "source", "old"] + list(NUMERIC))
        for (c, a), rec in sorted(recs.items()):
            st = rec["st"] or {}
            w.writerow([c, family_of(c) or "", nqubits(c, st) or "", a, rec["status"], rec["src"],
                        int(rec["old"])] + ["" if st.get(k) is None else st.get(k) for k in NUMERIC])
    print("  wrote", os.path.basename(path))

def main():
    global OUT, VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("paths", nargs="+", help="NEWDIR... OUTDIR")
    ap.add_argument("--old", action="append", default=[], help="old record: final.log-style file or status/JSON dir")
    ap.add_argument("--old-arms", default="qisq2_low,float_low")
    ap.add_argument("--qasm", action="append", default=[], help="directory of the circuits, for qubit counts")
    ap.add_argument("--dump", help="write every merged record, with its source, to this CSV")
    ap.add_argument("-v", "--verbose", action="store_true", help="list the circuits drawn red in each panel")
    a = ap.parse_args()
    if len(a.paths) < 1: ap.error("need OUTDIR")
    OUT, new = a.paths[-1], a.paths[:-1]
    VERBOSE = a.verbose
    old_arms = set(ARMS) if a.old_arms == "all" else set(a.old_arms.split(","))
    for d in a.qasm: count_qasm(d)
    os.makedirs(OUT, exist_ok=True)
    recs, rows, dropped = load(new, a.old, old_arms)
    provenance(recs, dropped)
    exact_agreement(rows)
    fams = {}
    for r in rows.values(): fams[r["_fam"]] = fams.get(r["_fam"], 0) + 1
    print(f"{len(rows)} circuits: " + ", ".join(f"{k}={v}" for k, v in sorted(fams.items())))
    print("panel: points (wrong, floored)  [family=points/wrong, how the float run was judged]")
    for field, lab in [("simulation_time", "runtime (s)"),
                       ("final_nodes", "final # of nodes"),
                       ("max_nodes", "peak # of nodes")]:
        scatter(rows, "qisq2_low", "limdd_qisq2", field,
                f"{lab} algebraic EVDD", f"{lab} algebraic LIMDD",
                f"limdd_vs_evdd_{field}.pdf")
        scatter(rows, "float_low", "limdd_float", field,
                f"{lab} float EVDD", f"{lab} float LIMDD",
                f"limdd_vs_evdd_float_{field}.pdf",
                judge="limdd_float", truth="limdd_qisq2")
        scatter(rows, "float_low", "qisq2_low", field,
                f"{lab} float", f"{lab} algebraic",
                f"evdd_float_vs_algebraic_{field}.pdf", judge="float_low")
        scatter(rows, "limdd_float", "limdd_qisq2", field,
                f"{lab} float LIMDD", f"{lab} algebraic LIMDD",
                f"limdd_float_vs_algebraic_{field}.pdf", judge="limdd_float")
    legend("panel_legend.pdf")
    if a.dump: dump(recs, a.dump)

if __name__ == "__main__":
    main()
