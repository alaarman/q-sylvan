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
run_gw.sh status lines), only the norm can be checked: if it is 1, the
measurement is not checked and the run is drawn grey, as Figures 6 to 8
(make_section6_figures.py) draw it. Both scripts judge through judge() below.
Since 30 Oct 2026 the panels colour four verdicts apart (kind()): red when
only the measurement is wrong, violet when the norm is wrong and the
measurement right or not checked, orange when both are wrong, grey when the
norm is right and the measurement not checked; in LIMDD against EVDD the fill
is the float LIMDD's verdict and the edge the float EVDD's. A run with 0 nodes
sits in a band below the data, in the earlier look too.
--unchecked right counts such a run as correct instead, as the figures of the
2026-10-01 revision did (regen_amendment.sh passes it for those), and keeps
that revision's legend.
Without it, the float-against-algebraic panels (evdd_ and
limdd_float_vs_algebraic_*) draw a float run whose norm alone is wrong violet,
apart from a wrong measurement, as Figures 6 to 8 do, and a float run that did
not finish hollow grey, since its correctness is unknown; each set of three
gets its own legend, evdd_ and limdd_float_vs_algebraic_legend.pdf.

--rule section6 applies the rule of Section 6 of the paper instead: a float run
is wrong when its top-qubit probability is more than 5% off the exact one, and
its norm is not looked at. The runner's error value (the state is zero) is off
every probability, so it is wrong without an exact run to compare with; any
other run with no exact probability to compare with is not judged and is drawn
grey. --float-evdd picks the float EVDD arm of the lower row and of the EVDD
panels: float_low (the default) or float_l2, the 'L2' normalisation that
Section 6 calls float. --float-limdd likewise picks the float LIMDD arm of the
lower row and of the LIMDD panels: limdd_float (the default, low) or
limdd_float_l2, the LIMDD with every node read at norm 1 (limdd_set_l2), or
limdd_float_rel and limdd_float_l2_rel, the same two under the relative
merging rule (run_arm.sh; the arms of the 19 October rerun).
"""
import argparse, json, math, os, re, sys, glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT  = None                                    # set in main
VERBOSE = False
ARMS = ["qisq2_low", "limdd_qisq2", "float_low", "float_l2", "limdd_float", "limdd_float_l2"]
LIMDD_ARMS = ("limdd_qisq2", "limdd_float", "limdd_float_l2")
EXACT = ("limdd_qisq2", "qisq2_low")           # either can serve as the truth
SENTINEL = 1e10
RULE = "amend"                                 # or "section6" (--rule)
UNCHECKED = "right"                            # or "grey" (--unchecked; its default); the
                                               # numbers scripts keep the revision's "right"
GREY = "0.55"                                  # a float run that was not judged
NORMCOL = "tab:purple"                         # a float run whose norm alone is wrong (split panels)
PAIRCOL = "tab:red"                            # a wrong float run in a panel of two float runs, the
                                               # float LIMDD of LIMDD against EVDD (--pair-wrong-colour)
UNJUDGED = set()                               # panels with a run that was not judged
DRAWN = {}                                     # panel -> (families, verdicts) drawn, for split_legend
WIDE = 64                                      # one 64-bit Pauli word
NUMERIC = ("simulation_time", "final_nodes", "max_nodes", "norm",
           "first_qubit_measurement_prob", "t_count", "final_width")

# family -> (marker, label); order fixes the legend
FAMILY = [
    ("random",  "o", "random Clifford+$T$"),
    ("grover",  "s", "Grover"),
    ("wstate",  "^", "W state"),
    ("hshift",  "D", "hidden shift"),
    # the new families of 6-7 October (pilot/qasm/<family>/ of the data)
    ("mirror",     "v", "mirror"),
    ("iqp",        "p", "IQP"),
    ("groveriter", "h", "Grover, many iterations"),
    ("qwalk",      "X", "quantum walk"),
    ("synth",      "*", "synthesised (SK, gridsynth)"),
]
NEW_FAMILIES = ("mirror", "iqp", "groveriter", "qwalk", "synth")
# AMEND_FAMILIES=random,grover,wstate keeps only those families, in the figures,
# the legend and every number the other scripts take from family_of. Unset, all
# are kept (the new families only appear where their runs are given).
ONLY = [f for f in os.environ.get("AMEND_FAMILIES", "").split(",") if f]
if ONLY:
    unknown = sorted(set(ONLY) - {k for k, _, _ in FAMILY})
    if unknown: sys.exit(f"AMEND_FAMILIES: unknown families {unknown}")
    FAMILY = [f for f in FAMILY if f[0] in ONLY]
MARK = {k: m for k, m, _ in FAMILY}

def family_of(name):
    if name.startswith("clifford_T_circuit"): fam = "random"
    elif name.startswith("grover_iter"):      fam = "groveriter"
    elif name.startswith("grover"):           fam = "grover"
    elif name.startswith(("w-state", "wclif")): fam = "wstate"
    elif name.startswith("hidden-shift"):     fam = "hshift"
    elif name.startswith("mirror"):           fam = "mirror"
    elif name.startswith("iqp"):              fam = "iqp"
    elif name.startswith("qwalk"):            fam = "qwalk"
    elif name.startswith("synth_"):           fam = "synth"
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

REFERENCE = {}                                 # circuit -> state-vector p0 (--reference)

def add_reference(rows):
    """A state-vector probability (reference_probs.py) as the pseudo-arm
    'statevector' of the circuits that have one: judge_section6 takes it only
    when no exact arm has a probability."""
    for c, r in rows.items():
        if c in REFERENCE: r["statevector"] = {"first_qubit_measurement_prob": REFERENCE[c]}

def read_reference(path):
    import csv
    for row in csv.DictReader(open(path)): REFERENCE[row["circuit"]] = float(row["p0"])

def judge_section6(r, judge, ref):
    """Section 6's rule, (bad, how): bad is None when the run is not judged.
    The exact probability comes from the panel's exact arm, else the other
    exact arm, else a state-vector reference (--reference)."""
    pf = p(r[judge])
    if pf is not None and (pf >= SENTINEL or pf != pf): return True, "sentinel"
    if pf is None: return None, "no-p"
    for a in (ref,) + tuple(e for e in EXACT if e != ref) + ("statevector",):
        if a in r and p(r[a]) is not None:
            pa = p(r[a])
            return abs(pf - pa) > 0.05 * max(abs(pa), 1e-12), ("p" if a == ref else "p:" + a)
    return None, "no-exact"

def judge(r, arm, ref):
    """The float run r[arm] by the amendment's rule, as every figure judges it:
    (verdict, how). verdict is 'meas' when its top-qubit probability is the
    runner's error value or more than 5% off the exact one (Section 6's rule,
    judge_section6), 'norm' when its measurement is not wrong but its norm is
    more than 1e-3 off 1, 'unchecked' when its norm is 1 and there is no exact
    probability to compare its own with, and 'right' otherwise. how is
    'norm!=1' whenever the norm is off, else judge_section6's.

    The 5% rule alone is not enough: on the 30-qubit hidden shift the float
    run returns norm 0 -- the state has been annihilated -- while its
    top-qubit probability coincidentally equals the exact one, so the rule
    calls it correct. A unitary circuit cannot produce an unnormalised state,
    so a norm away from 1 is a wrong answer whatever the measurement says."""
    meas, how = judge_section6(r, arm, ref)          # True, False or None
    nf = norm_of(r[arm])
    normbad = nf is None or nf != nf or abs(nf - 1.0) > 1e-3
    if meas: v = "meas"
    elif normbad: v = "norm"
    elif meas is None: v = "unchecked"
    else: v = "right"
    return v, ("norm!=1" if normbad else how)

def judge_point(r, judge_arm, ref):
    """(bad, how) as the panels take it: bad is True for judge()'s 'meas' and
    'norm', and for 'unchecked' None (drawn grey) or, with UNCHECKED 'right',
    False. how is 'norm!=1', 'sentinel' (the runner's error value for the
    probability), 'p' (5% rule against the panel's own exact arm), 'p:<arm>'
    (against the other exact arm, the own one did not finish), or 'norm-only'
    (norm is 1 and there is nothing to compare the probability with: no exact
    arm finished, or the float record is a run_gw.sh status line without a
    probability)."""
    if RULE == "section6": return judge_section6(r, judge_arm, ref)
    v, how = judge(r, judge_arm, ref)
    if v == "unchecked": return (None if UNCHECKED == "grey" else False), "norm-only"
    return v != "right", how

# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------
# A run that did not finish (a timeout or a full table) is drawn on a dashed
# line beyond the data: on the right when the x arm failed, on top when the y
# arm failed, and in the corner when both did. A run that was not made
# (SKIPPED_PREDICTED) is not drawn.
NOT_RUN = ("SKIPPED_PREDICTED", "EMPTY", "BADLINE")

def failures(recs):
    """(circ, arm) -> status, for every plotted family's run that did not finish."""
    return {k: rec["status"] for k, rec in recs.items()
            if rec["status"] != "OK" and rec["status"] not in NOT_RUN
            and family_of(k[0]) is not None}

def kinds():
    """True when the panels colour the four verdicts of kind() apart (the
    paper's figures since 30 Oct 2026); --unchecked right and --rule section6
    keep the earlier look."""
    return UNCHECKED == "grey" and RULE != "section6"

KIND_COL = {"meas": "tab:red", "norm": NORMCOL, "both": "tab:orange", "unjudged": GREY}

def kind(r, arm, ref):
    """judge()'s verdict as the panels colour it: 'meas' (measurement wrong,
    norm right), 'norm' (norm wrong, measurement right or not checked), 'both'
    (measurement and norm wrong), 'unjudged' (norm right, measurement not
    checked: no exact run to compare with) or 'right'."""
    v, how = judge(r, arm, ref)
    if v == "meas" and how == "norm!=1": return "both"
    return "unjudged" if v == "unchecked" else v

def colours(y, x, split=False):
    """y, x: 'wrong', 'right', 'unjudged', or None where the panel judges no
    such float run; a split panel's y can also be 'norm' (only the norm is
    wrong) or 'failed' (the float run did not finish). Fill: the y arm's float
    run (or the only float arm), red when wrong, violet when only its norm is,
    grey when not judged, none when it did not finish. Edge: the x arm's float
    run, where the panel judges one (the lower row), orange when wrong, grey
    when not judged; otherwise the edge goes with the fill, and is grey for a
    float run that did not finish.
    With kinds(), y and x are kind()'s verdicts (or 'failed'), and fill and
    edge take the same colours: red, violet, orange and grey for 'meas',
    'norm', 'both' and 'unjudged'; a right x run has a blue edge, and one
    that did not finish ('failed') a grey one."""
    if kinds():
        face = KIND_COL.get(y, "none")
        if x is not None: edge = GREY if x == "failed" else KIND_COL.get(x, "tab:blue")
        elif face != "none": edge = face
        elif y == "failed": edge = GREY
        else: edge = "tab:blue"
        return face, edge
    face = {"wrong": "tab:red" if split else PAIRCOL, "norm": NORMCOL, "unjudged": GREY}.get(y, "none")
    if x == "wrong": edge = "tab:orange"
    elif x == "unjudged": edge = GREY
    elif y in ("wrong", "norm"): edge = face
    elif y in ("unjudged", "failed") and x is None: edge = GREY
    else: edge = "tab:blue"
    return face, edge

def verdict(bad):
    """judge_point's bad as colours() takes it."""
    return "unjudged" if bad is None else ("wrong" if bad else "right")

def norm_only(r, arm, ref):
    """True when judge() finds the float run's measurement right and only its
    norm off 1. The split panels draw such a run apart, as Figures 6 to 8 did."""
    return RULE != "section6" and judge(r, arm, ref)[0] == "norm"

def band(ax, lo, F, gap, top, nboth, fs):
    """The shaded band beyond the data, its dashed lines and their labels."""
    ax.axvspan(F / gap ** 0.5, top, color="0.94", lw=0, zorder=0)
    ax.axhspan(F / gap ** 0.5, top, color="0.94", lw=0, zorder=0)
    ax.axvline(F, ls=(0, (3, 2)), color="0.45", lw=0.7, zorder=1)
    ax.axhline(F, ls=(0, (3, 2)), color="0.45", lw=0.7, zorder=1)
    ax.text(lo * 1.2, F * gap ** 0.18, "did not finish", fontsize=fs, color="0.35", va="bottom")
    ax.text(F * gap ** 0.18, lo * 1.2, "did not finish", fontsize=fs, color="0.35",
            rotation=90, ha="left", va="bottom")
    if nboth:
        ax.text(F / gap ** 0.35, F * gap ** 0.18, f"{nboth}", fontsize=fs, color="0.35",
                ha="right", va="bottom")

def zero_band(ax, lo, Z, gap, bot, zx, zy, fs):
    """The shaded band below the data of the runs with 0 nodes, on the axes
    that have one: its dotted line at Z and its label."""
    if zx:
        ax.axvspan(bot, Z * gap ** 0.5, color="0.94", lw=0, zorder=0)
        ax.axvline(Z, ls=(0, (1, 1.5)), color="0.45", lw=0.7, zorder=1)
        ax.text(Z / gap ** 0.18, lo * 1.2, "0 nodes", fontsize=fs, color="0.35",
                rotation=90, ha="right", va="bottom")
    if zy:
        ax.axhspan(bot, Z * gap ** 0.5, color="0.94", lw=0, zorder=0)
        ax.axhline(Z, ls=(0, (1, 1.5)), color="0.45", lw=0.7, zorder=1)
        ax.text(lo * 1.2, Z / gap ** 0.18, "0 nodes", fontsize=fs, color="0.35", va="top")

def scatter(rows, xarm, yarm, field, xlabel, ylabel, fname, judge=None, truth=None,
            fails=None, xjudge=None, fams=None, split=False):
    """One panel. fams keeps only those families. split (float against
    algebraic, with --unchecked grey) draws a float run whose norm alone is
    wrong violet, and one that did not finish hollow grey, since its
    correctness is unknown; split_legend() says so."""
    fails = fails or {}
    split = split and UNCHECKED == "grey"      # --unchecked right keeps the revision's look
    pts = {}          # (family, face, edge) -> ([x],[y]); None for a failed arm
    nfloor = 0
    vk = {"y": {}, "x": {}}             # how many points per verdict, for the report
    how = {}          # how each judged point was judged
    perfam = {}       # family -> [points, wrong]
    detail = []       # (circuit, how) of every point drawn red, for -v
    nfail = {"x": 0, "y": 0, "both": 0}
    nj = {"y": 0, "x": 0}                      # float runs not judged (--rule section6)
    short = lambda c: c[19:] if c.startswith("clifford_T_circuit_") else c
    circs = set(rows) | {c for (c, a) in fails if a in (xarm, yarm)}
    for circ in sorted(circs):                 # sorted: a set's order varies run to run
        r = rows.get(circ, {"_fam": family_of(circ)})
        if fams is not None and r["_fam"] not in fams: continue
        xok, yok = xarm in r, yarm in r
        xf, yf = (circ, xarm) in fails, (circ, yarm) in fails
        if not (xok or xf) or not (yok or yf): continue        # an arm not run
        if not xok and not yok: x = y = None; nfail["both"] += 1
        else:
            try:
                x = float(r[xarm][field]) if xok else None
                y = float(r[yarm][field]) if yok else None
            except (KeyError, TypeError, ValueError): continue
            if x is None: nfail["x"] += 1
            if y is None: nfail["y"] += 1
        ref = truth if truth else yarm
        yv = xv = None
        bad, h = False, None
        if judge and judge in r:
            bad, h = judge_point(r, judge, ref); yv = verdict(bad)
            if kinds(): yv = kind(r, judge, ref)
            elif split and bad and norm_only(r, judge, ref): yv = "norm"
        elif split and judge and (circ, judge) in fails: yv = "failed"
        xbad, xh = False, None
        if xjudge and xjudge in r:
            xbad, xh = judge_point(r, xjudge, EXACT[1])
            xv = kind(r, xjudge, EXACT[1]) if kinds() else verdict(xbad)
        elif xjudge and kinds() and (circ, xjudge) in fails: xv = "failed"
        if (x is not None and x <= 0) or (y is not None and y <= 0):
            # 0 nodes: drawn in a band below the data. It keeps its verdict:
            # a LIMDD of 0 nodes is a root edge, which holds a product
            # stabiliser state (mirror-none-plus_*: |+>^n) as well as the zero
            # vector, and judge() tells them apart by the norm. Until 19 Oct
            # 2026 such a point was always drawn wrong; until 30 Oct it was
            # clipped by the axis (and dropped from the exact panels).
            nfloor += 1; h = (h or "") + "+zero"
        if h: how[h] = how.get(h, 0) + 1
        if bad: detail.append(f"{short(circ)}:{h}")
        if xbad: detail.append(f"{short(circ)}:x:{xh}")
        if bad is None: detail.append(f"{short(circ)}:unjudged"); nj["y"] += 1
        if xbad is None: detail.append(f"{short(circ)}:x:unjudged"); nj["x"] += 1
        pf = perfam.setdefault(r["_fam"], [0, 0]); pf[0] += 1; pf[1] += bool(bad)
        if bad is None or xbad is None: UNJUDGED.add(fname)
        face, edge = colours(yv, xv, split)
        drawn = DRAWN.setdefault(fname, (set(), set()))
        drawn[0].add(r["_fam"]); drawn[1].add(yv or "right")
        if xv: drawn[1].add(xv)
        vk["y"][yv] = vk["y"].get(yv, 0) + 1; vk["x"][xv] = vk["x"].get(xv, 0) + 1
        pts.setdefault((r["_fam"], face, edge), ([], []))
        pts[(r["_fam"], face, edge)][0].append(x); pts[(r["_fam"], face, edge)][1].append(y)
    allv = [v for (xs, ys) in pts.values() for v in xs + ys if v is not None and v > 0]
    if not allv: print(f"  {fname}: no data"); return
    lo, hi = min(allv) * 0.5, max(allv) * 2.0
    anyfail = any(nfail.values())
    gap = 10 ** (0.08 * math.log10(hi / lo))  # the band of the runs that did not finish: 8% of an axis
    F = hi * gap                               # its dashed line
    top = F * gap if anyfail else hi
    zx = any(v is not None and v <= 0 for xs, _ in pts.values() for v in xs)
    zy = any(v is not None and v <= 0 for _, ys in pts.values() for v in ys)
    Z = lo / gap                               # the dotted line of the runs with 0 nodes
    bot = Z / gap

    fig, ax = plt.subplots(figsize=(2.9, 2.6))
    ax.plot([lo, hi], [lo, hi], color="0.6", lw=0.8, zorder=1)
    if anyfail:
        band(ax, lo, F, gap, top, nfail["both"], 5.5)
    zero_band(ax, lo, Z, gap, bot, zx, zy, 5.5)
    plain = ("none", "tab:blue", GREY)
    for (fam, face, edge), (xs, ys) in sorted(pts.items()):
        xs = [F if v is None else (v if v > 0 else Z) for v in xs]
        ys = [F if v is None else (v if v > 0 else Z) for v in ys]
        ax.scatter(xs, ys, s=20, marker=MARK[fam], facecolors=face, edgecolors=edge,
                   linewidths=1.1 if xjudge and kinds() else 0.9,
                   zorder=3 + (face != "none") + (edge not in plain and edge != face))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(bot if zx else lo, top); ax.set_ylim(bot if zy else lo, top)
    if anyfail or zx or zy:                    # no tick labels in the bands
        from matplotlib.ticker import FixedLocator, LogLocator, NullFormatter
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_locator(FixedLocator([t for t in LogLocator().tick_values(lo, hi) if lo <= t <= hi]))
            axis.set_minor_locator(FixedLocator(
                [t for t in LogLocator(subs=range(2, 10)).tick_values(lo, hi) if lo <= t <= hi]))
            axis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(xlabel, fontsize=7.5, labelpad=1.5)
    ax.set_ylabel(ylabel, fontsize=7.5, labelpad=1.5)
    ax.tick_params(labelsize=6.5, pad=1.5)
    fig.tight_layout(pad=0.25)
    fig.savefig(os.path.join(OUT, fname), bbox_inches="tight"); plt.close(fig)
    n = sum(len(v[0]) for v in pts.values())
    nb = sum(len(v[0]) for k, v in pts.items() if k[1] == ("tab:red" if split else PAIRCOL))
    nn = sum(len(v[0]) for k, v in pts.items() if k[1] == NORMCOL)
    nfl = sum(len(v[0]) for k, v in pts.items() if k[1] == "none" and k[2] == GREY)
    nxb = sum(len(v[0]) for k, v in pts.items() if k[2] == "tab:orange")
    fams = " ".join(f"{k}={v[0]}/{v[1]}w" for k, v in sorted(perfam.items()))
    hw = ", judged by " + " ".join(f"{k}={v}" for k, v in sorted(how.items())) if how else ""
    fl = (f"; did not finish: x {nfail['x']}, y {nfail['y']}, both {nfail['both']}" if anyfail else "")
    njs = (f", {nj['y']} not judged" if nj["y"] else "") + (f", {nj['x']} x not judged" if nj["x"] else "")
    sp = f" (of which {nn} norm only), {nfl} float did not finish" if split else ""
    if kinds():                                # the verdicts by colour: fill for y, edge for x
        cnt = lambda side, k: vk["y" if side == 1 else "x"].get(k, 0)
        ks = lambda side: ", ".join(f"{k} {cnt(side, k)}" for k in ("meas", "norm", "both", "unjudged"))
        sp = f": {ks(1)}" + (f"; x: {ks(2)}" if xjudge else "") + (f"; {nfl} float did not finish" if split else "")
        nb, nn = sum(cnt(1, k) for k in ("meas", "norm", "both")), 0
        nxb = sum(cnt(2, k) for k in ("meas", "norm", "both"))
    print(f"  {fname}: {n} points ({nb + nn} wrong{sp}{', %d x wrong' % nxb if xjudge else ''}{njs}"
          f"{', %d with 0 nodes' % nfloor if nfloor else ''}{fl})  [{fams}{hw}]")
    if VERBOSE and detail: print("      red/orange: " + " ".join(sorted(detail)))

def legend(fname, grey=False):
    """One shared legend, so the panels keep their space. grey adds the runs
    that were not judged (no exact probability to compare with). The labels
    say what was checked: the measurement alone under --rule section6, the
    measurement and the norm otherwise (worded as in the 2026-10-01 revision
    with --unchecked right)."""
    if kinds(): return pair_legend(fname)
    fig = plt.figure(figsize=(6.4, 0.62))
    hs  = [plt.Line2D([], [], ls="", marker=m, mfc="none", mec="0.25", mew=0.9,
                      ms=5, label=lab) for _, m, lab in FAMILY]
    both = RULE != "section6" and UNCHECKED == "grey"
    xlab = ("float EVDD measurement wrong (lower row)" if RULE == "section6"
            else "float EVDD wrong (lower row)")
    hs += [plt.Line2D([], [], ls="", marker="o", mfc="none", mec="tab:blue",
                      mew=0.9, ms=5, label="float correct" if both else "measurement correct"),
           plt.Line2D([], [], ls="", marker="o", mfc=PAIRCOL, mec=PAIRCOL,
                      ms=5, label="float (LIMDD) wrong" if both else "float (LIMDD) measurement wrong"),
           plt.Line2D([], [], ls="", marker="o", mfc="none", mec="tab:orange",
                      mew=0.9, ms=5, label=xlab)]
    if grey:
        hs += [plt.Line2D([], [], ls="", marker="o", mfc=GREY, mec=GREY, ms=5,
                          label="float measurement not checked" if both
                          else "not judged (no exact run finished)")]
    hs += [plt.Line2D([], [], ls=(0, (3, 2)), color="0.45", lw=0.7,
                      label="did not finish (timeout or full table)")]
    fig.legend(handles=hs, loc="center", ncol=4, fontsize=7, frameon=False,
               handletextpad=0.35, columnspacing=1.1)
    fig.savefig(os.path.join(OUT, fname), bbox_inches="tight"); plt.close(fig)
    print("  wrote", fname)

def split_legend(fname, panels):
    """The legend of a figure of split panels (float against algebraic: Figure
    6 of Section 6 and Figure 9), listing the families and the verdicts that
    its panels draw."""
    fams, verds = set(), set()
    for p in panels:
        f, k = DRAWN.get(p, (set(), set()))
        fams |= f; verds |= k
    dot = lambda fc, ec, lab: plt.Line2D([], [], ls="", marker="o", mfc=fc, mec=ec,
                                         mew=0.9, ms=5, label=lab)
    hs  = [plt.Line2D([], [], ls="", marker=m, mfc="none", mec="0.25", mew=0.9,
                      ms=5, label=lab) for k, m, lab in FAMILY if k in fams]
    hs += [dot("none", "tab:blue", "float correct")]
    if "wrong" in verds: hs += [dot("tab:red", "tab:red", "float measurement wrong")]
    hs += [dot(KIND_COL[k], KIND_COL[k], "float " + KIND_LAB[k]) for k in KIND_LAB if k in verds]
    if "failed" in verds: hs += [dot("none", GREY, "float run did not finish")]
    hs += [plt.Line2D([], [], ls=(0, (3, 2)), color="0.45", lw=0.7,
                      label="did not finish (timeout or full table)")]
    fig = plt.figure(figsize=(6.4, 0.62))
    fig.legend(handles=hs, loc="center", ncol=4, fontsize=7, frameon=False,
               handletextpad=0.35, columnspacing=1.1)
    fig.savefig(os.path.join(OUT, fname), bbox_inches="tight"); plt.close(fig)
    print("  wrote", fname)

KIND_LAB = {"meas": "measurement wrong", "norm": "norm wrong",
            "both": "measurement and norm wrong", "unjudged": "measurement not checked"}

def pair_legend(fname):
    """The legend of LIMDD against EVDD with kinds(): the families, and the
    verdicts that its lower row draws, the float LIMDD's as the fill and the
    float EVDD's as the edge."""
    verds = set()
    for p, (f, k) in DRAWN.items():
        if p.startswith("limdd_vs_evdd_float"): verds |= k
    dot = lambda fc, ec, lab: plt.Line2D([], [], ls="", marker="o", mfc=fc, mec=ec,
                                         mew=1.1, ms=5, label=lab)
    hs  = [plt.Line2D([], [], ls="", marker=m, mfc="none", mec="0.25", mew=0.9,
                      ms=5, label=lab) for _, m, lab in FAMILY]
    hs += [dot("none", "tab:blue", "float correct")]
    hs += [dot(KIND_COL[k], KIND_COL[k], KIND_LAB[k]) for k in KIND_LAB if k in verds]
    hs += [plt.Line2D([], [], ls="", label="bottom row: fill float LIMDD, edge float EVDD"),
           plt.Line2D([], [], ls=(0, (3, 2)), color="0.45", lw=0.7,
                      label="did not finish (timeout or full table)")]
    fig = plt.figure(figsize=(6.4, 0.8))
    fig.legend(handles=hs, loc="center", ncol=4, fontsize=7, frameon=False,
               handletextpad=0.35, columnspacing=1.1)
    fig.savefig(os.path.join(OUT, fname), bbox_inches="tight"); plt.close(fig)
    print("  wrote", fname)

# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------
def provenance(recs, dropped):
    print("records per arm (status OK / all), by source:")
    for a in ARMS:
        if not any(arm == a for (c, arm) in recs): continue   # e.g. float_l2, when not given
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
    global OUT, VERBOSE, RULE, UNCHECKED
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("paths", nargs="+", help="NEWDIR... OUTDIR")
    ap.add_argument("--old", action="append", default=[], help="old record: final.log-style file or status/JSON dir")
    ap.add_argument("--old-arms", default="qisq2_low,float_low")
    ap.add_argument("--qasm", action="append", default=[], help="directory of the circuits, for qubit counts")
    ap.add_argument("--dump", help="write every merged record, with its source, to this CSV")
    ap.add_argument("-v", "--verbose", action="store_true", help="list the circuits drawn red in each panel")
    ap.add_argument("--rule", choices=("amend", "section6"), default="amend",
                    help="correctness rule: the amendment's (norm and 5%%) or Section 6's (5%% alone)")
    ap.add_argument("--float-evdd", choices=("float_low", "float_l2"), default="float_low",
                    help="the float EVDD arm of the lower row and the EVDD panels")
    ap.add_argument("--float-limdd", choices=("limdd_float", "limdd_float_l2", "limdd_float_rel", "limdd_float_l2_rel"),
                    default="limdd_float",
                    help="the float LIMDD arm of the lower row and the LIMDD panels")
    ap.add_argument("--unchecked", choices=("grey", "right"), default="grey",
                    help="a float run with norm 1 and no exact probability to compare with: "
                         "drawn grey, as Figures 6 to 8 draw it, or counted as correct "
                         "(the 2026-10-01 revision)")
    ap.add_argument("--reference", action="append", default=[],
                    help="CSV of state-vector top-qubit probabilities (reference_probs.py), used where "
                         "no exact run has one")
    ap.add_argument("--pair-wrong-colour", default="tab:red",
                    help="fill of a wrong float LIMDD run in the LIMDD-against-EVDD panels, so that red "
                         "can keep one meaning (measurement wrong) in the float-against-algebraic panels "
                         "(the earlier look only: --unchecked right or --rule section6)")
    ap.add_argument("--no-norm-label", action="store_true",
                    help="no '(low)' or '(L2)' in the axis labels (the text says which normalisation)")
    ap.add_argument("--low-label", action="store_true",
                    help="'(low)' after a low float arm in the axis labels, as '(L2)' after an L2 one")
    a = ap.parse_args()
    global PAIRCOL, LIMDD_ARMS
    RULE, UNCHECKED, PAIRCOL = a.rule, a.unchecked, a.pair_wrong_colour
    FE, FL = a.float_evdd, a.float_limdd
    if FL not in ARMS:                         # a relative-rule arm: read its records too
        ARMS.append(FL); LIMDD_ARMS = LIMDD_ARMS + (FL,)
    low = " (low)" if a.low_label else ""     # without it the labels stay as the paper has them
    fl = " (L2)" if FE == "float_l2" else low
    ll = " (L2)" if FL.startswith("limdd_float_l2") else low
    if a.no_norm_label: fl = ll = ""
    if len(a.paths) < 1: ap.error("need OUTDIR")
    OUT, new = a.paths[-1], a.paths[:-1]
    VERBOSE = a.verbose
    old_arms = set(ARMS) if a.old_arms == "all" else set(a.old_arms.split(","))
    for d in a.qasm: count_qasm(d)
    os.makedirs(OUT, exist_ok=True)
    for f in a.reference: read_reference(f)
    recs, rows, dropped = load(new, a.old, old_arms)
    add_reference(rows)
    provenance(recs, dropped)
    exact_agreement(rows)
    fams = {}
    for r in rows.values(): fams[r["_fam"]] = fams.get(r["_fam"], 0) + 1
    print(f"{len(rows)} circuits: " + ", ".join(f"{k}={v}" for k, v in sorted(fams.items())))
    fails = failures(recs)
    print("panel: points (wrong, 0 nodes)  [family=points/wrong, how the float run was judged]")
    for field, lab in [("simulation_time", "runtime (s)"),
                       ("final_nodes", "final # of nodes"),
                       ("max_nodes", "peak # of nodes")]:
        scatter(rows, "qisq2_low", "limdd_qisq2", field,
                f"{lab} algebraic EVDD", f"{lab} algebraic LIMDD",
                f"limdd_vs_evdd_{field}.pdf", fails=fails)
        scatter(rows, FE, FL, field,
                f"{lab} float EVDD{fl}", f"{lab} float LIMDD{ll}",
                f"limdd_vs_evdd_float_{field}.pdf",
                judge=FL, truth="limdd_qisq2", fails=fails, xjudge=FE)
        scatter(rows, FE, "qisq2_low", field,
                f"{lab} float{fl}", f"{lab} algebraic",
                f"evdd_float_vs_algebraic_{field}.pdf", judge=FE, fails=fails, split=True)
        scatter(rows, FL, "limdd_qisq2", field,
                f"{lab} float LIMDD{ll}", f"{lab} algebraic LIMDD",
                f"limdd_float_vs_algebraic_{field}.pdf", judge=FL, fails=fails, split=True)
    legend("panel_legend.pdf", grey=any(f.startswith("limdd_vs_evdd") for f in UNJUDGED))
    if UNCHECKED == "grey":
        for d in ("evdd", "limdd"):
            split_legend(f"{d}_float_vs_algebraic_legend.pdf",
                         [f"{d}_float_vs_algebraic_{f}.pdf" for f in ("simulation_time", "final_nodes", "max_nodes")])
    if a.dump: dump(recs, a.dump)

if __name__ == "__main__":
    main()
