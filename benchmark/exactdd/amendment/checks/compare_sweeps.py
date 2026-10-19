#!/usr/bin/env python3
"""Two sweeps of the same arms against each other, per arm: how many finish,
whether the exact arms keep their node counts (they must, on any machine and
any build of the same code), and the runtime ratio new/old on the circuits
both finish. For a repeat of the paper's sweep: OLD is the stored data, NEW
the repeat (sweep.py's OUT).

  compare_sweeps.py OLD NEW [FAMILIES]
OLD, NEW   roots in the layout of README.md; amend/out_clean and, where
           present, pilot/out_clean are compared; circuits from OLD's qasm dirs
FAMILIES   comma-separated (default: every family), e.g. random,grover,wstate

With the stored data as OLD and the records behind the figures of 5 October
(OLD/ldd/out, amend/out_rnd, amend/out and the October reruns, with the
September EVDD records) as NEW, give --paper5 as OLD's first argument instead:
  compare_sweeps.py --paper5 ROOT   (the check of the clean sweep's node counts)
Run with the plotting Python (it loads make_figures.py)."""
import glob, os, statistics, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import make_figures as mf
PAPER_ARMS = list(mf.ARMS)
mf.ARMS = PAPER_ARMS + ["limdd_qisq2_lazy", "limdd_lazy", "limdd_float_rel", "limdd_float_l2_rel"]   # read these too

a = sys.argv[1:]
if len(a) < 2: sys.exit(__doc__)
def outs(root): return [d for d in (root + "/amend/out_clean", root + "/pilot/out_clean", root + "/amend/out_clean_lazy",
                                     root + "/pilot/out_clean_lazy") if os.path.isdir(d)]
if a[0] == "--paper5":
    S = a[1]; FAMS = {"random", "grover", "wstate"}
    OLDREC = [f"{S}/exactdd/final.log", f"{S}/gw/out", f"{S}/gw/out2"]
    NEW5 = [f"{S}/{d}" for d in ("ldd/out", "amend/out_rnd", "amend/out", "amend/out_l2", "amend/out_low",
                                "amend/out_limdd_l2", "amend/out_limdd_low")]
    old_recs, old_rows, _ = mf.load(NEW5, OLDREC, {"qisq2_low", "float_low"})
    new_recs, new_rows, _ = mf.load([S + "/amend/out_clean"], [], set())
else:
    S = a[0]; FAMS = set(a[2].split(",")) if len(a) > 2 else None
    old_recs, old_rows, _ = mf.load(outs(a[0]), [], set())
    new_recs, new_rows, _ = mf.load(outs(a[1]), [], set())
for d in [f"{S}/ldd/qasm", f"{S}/amend/qasm"] + glob.glob(f"{S}/pilot/qasm/*/"):
    if os.path.isdir(d): mf.count_qasm(d)

def status(recs, c, a):
    r = recs.get((c, a)); return r["status"] if r else "-"

def close(x, y):
    """Equal to the digits a status line prints (6), or one of them unrecorded
    (the September records have no probability)."""
    try: x, y = float(x), float(y)
    except (TypeError, ValueError): return True
    return x == y or abs(x - y) <= 1e-5 * max(abs(x), abs(y))

def st(recs, c, a):
    r = recs.get((c, a)); return r["st"] if r and r["status"] == "OK" and r["st"] else None

EXACT = set(mf.EXACT) | {"limdd_qisq2_lazy"}
circs = sorted({c for (c, a) in new_recs if FAMS is None or mf.family_of(c) in FAMS})
arms = [a for a in mf.ARMS if any(k[1] == a for k in new_recs)]
print(f"{len(circs)} circuits")
bad = 0
for arm in arms:
    ok_old = [c for c in circs if st(old_recs, c, arm)]
    ok_new = [c for c in circs if st(new_recs, c, arm)]
    both = [c for c in ok_old if c in ok_new]
    ratios = [float(st(new_recs, c, arm)["simulation_time"]) / max(float(st(old_recs, c, arm)["simulation_time"]), 1e-3)
              for c in both]
    key = lambda d: (d.get("final_nodes"), d.get("max_nodes"))
    same_nodes = [c for c in both if key(st(new_recs, c, arm)) == key(st(old_recs, c, arm))]
    gained = [c for c in ok_new if c not in ok_old]; lost = [c for c in ok_old if c not in ok_new]
    med = statistics.median(ratios) if ratios else float("nan")
    print(f"{arm:18s} OK old {len(ok_old):3d} new {len(ok_new):3d} | both {len(both):3d}: same final+peak nodes "
          f"{len(same_nodes):3d}, time new/old median {med:.2f} (min {min(ratios, default=0):.2f} max {max(ratios, default=0):.2f})")
    for c in gained: print(f"    now finishes: {c:28s} old {status(old_recs, c, arm):10s} "
                           f"new {float(st(new_recs, c, arm)['simulation_time']):8.1f}s")
    for c in lost:
        print(f"    no longer:    {c:28s} old {float(st(old_recs, c, arm)['simulation_time']):8.1f}s "
              f"new {status(new_recs, c, arm)}")
    for c in both:
        o, n = st(old_recs, c, arm), st(new_recs, c, arm)
        same_p = all(close(o.get(k), n.get(k)) for k in ("first_qubit_measurement_prob", "norm"))
        if key(o) != key(n) or (arm in EXACT and not same_p):
            if arm in EXACT: bad += 1
            print(f"    {'EXACT ' if arm in EXACT else ''}DIFFERS: {c} nodes old {key(o)} new {key(n)}; "
                  f"p0 old {o.get('first_qubit_measurement_prob')} new {n.get('first_qubit_measurement_prob')}; "
                  f"norm old {o.get('norm')} new {n.get('norm')}")
print(f"exact runs that differ: {bad}")
sys.exit(1 if bad else 0)
