#!/usr/bin/env python3
"""
Check arms against a reference arm at sampled basis states.

  amp_check.py <outdir> [--ref qisq2_low] [--tol 1e-9]

<outdir> holds run_arm.sh's <circ>_<arm>.json. The reference arm was run with
QSY_AMP_SAMPLE=K (--amp-sample), which records K basis states drawn from its
final state by |amplitude|^2 together with their amplitudes; the other arms were
run with QSY_AMP_AT=<dir of the reference> (--amp-at) and record their
amplitudes at the same points. A run agrees when every amplitude is within
tol times the largest reference amplitude of the reference's. This is what
checks an arm on the circuits of 25 to 63 qubits, whose state is not decoded,
so that norm and probability read -1 there for the BQD.

Prints one line per (circuit, arm) and a count per arm: agree, differ, no
points (the run did not finish, or ran without --amp-at), no reference.
"""
import argparse, collections, glob, json, math, os

ap = argparse.ArgumentParser()
ap.add_argument("outdir")
ap.add_argument("--ref", default="qisq2_low")
ap.add_argument("--tol", type=float, default=1e-9)
a = ap.parse_args()

# run_arm.sh's arms; a file is <circ>_<arm>.json, and circuit names have underscores
ARMS = ["qisq2_low", "float_low", "float_l2", "limdd_qisq2", "limdd_float", "limdd_float_l2",
        "limdd_heur", "limdd_lazy", "limdd_qisq2_lazy", "bqd_scalar", "bqd_x", "bqd_pauli"]
runs = collections.defaultdict(dict)
for f in glob.glob(os.path.join(a.outdir, "*.json")):
    name = os.path.basename(f)[:-5]
    arm = max((x for x in ARMS if name.endswith("_" + x)), key=len, default=None)
    if arm is None:
        continue
    try:
        d = json.load(open(f))
    except Exception:
        continue
    runs[name[:-len(arm) - 1]][arm] = d

count = collections.Counter()
for circ in sorted(runs):
    ref = runs[circ].get(a.ref)
    pts = ref.get("amp_points") if ref else None
    for arm, d in sorted(runs[circ].items()):
        if arm == a.ref:
            continue
        if not pts:
            count[(arm, "no reference")] += 1
            continue
        if d.get("amp_points") != pts:
            count[(arm, "no points")] += 1
            print("%-30s %-16s no points" % (circ, arm))
            continue
        rv = [complex(*v) for v in ref["amp_values"]]
        dv = [complex(*v) if v[0] is not None and v[1] is not None else complex("nan") for v in d["amp_values"]]
        scale = max(abs(v) for v in rv) or 1.0
        errs = [abs(x - y) / scale for x, y in zip(rv, dv)]
        # max() skips a NaN that is not first, so a null amplitude is checked apart
        err = float("nan") if any(math.isnan(e) for e in errs) else max(errs)
        ok = err <= a.tol
        count[(arm, "agree" if ok else "differ")] += 1
        print("%-30s %-16s %-6s max error %.3g of the largest amplitude (%d points)" % (
            circ, arm, "agree" if ok else "DIFFER", err, len(pts)))

print()
for arm in sorted({k[0] for k in count}):
    print("%-16s agree %3d  differ %3d  no points %3d  no reference %3d" % (
        arm, count[(arm, "agree")], count[(arm, "differ")], count[(arm, "no points")], count[(arm, "no reference")]))
