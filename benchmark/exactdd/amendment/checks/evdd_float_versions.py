#!/usr/bin/env python3
"""Float EVDD (qmdd, low, float): Sep 25-26 sweep (final.log, last line wins)
against the reruns on the code of the LIMDD runs (ROOT/amend/span/out2/qmdd_float_*_w24.json).

  evdd_float_versions.py ROOT      (or ROOT in the environment)
"""
import json, glob, os, re, sys
SP = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("ROOT")
if not SP: sys.exit("usage: evdd_float_versions.py ROOT")
old = {}
for line in open(f"{SP}/exactdd/final.log"):
    t = line.split()
    if len(t) >= 3 and t[0] == "float_low":
        old[t[1]] = t[2:]
ndiff = 0; n = 0
for f in sorted(glob.glob(f"{SP}/amend/span/out2/qmdd_float_*_w24.json")):
    c = re.search(r"qmdd_float_(\d+_\d+)_w24", f).group(1)
    st = json.load(open(f))["statistics"]
    o = old.get("clifford_T_circuit_" + c)
    new = (st["final_nodes"], st["max_nodes"], st["norm"], st["first_qubit_measurement_prob"], st["simulation_time"])
    if o and o[0] == "OK":
        of, op, on, opt = int(o[2]), int(o[3]), float(o[4]), float(o[5])
        same = (of, op) == new[:2] and abs(on - new[2]) < 1e-6
    else:
        same = False
    n += 1; ndiff += not same
    print(f"{c:8s} old {o}  new final {new[0]} peak {new[1]} norm {new[2]} p {new[3]} t {new[4]:.2f}  {'same' if same else 'DIFFER'}")
print(f"{ndiff} of {n} differ")
