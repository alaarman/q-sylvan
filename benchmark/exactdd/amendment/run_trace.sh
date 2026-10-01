#!/bin/bash
# The per-gate trace behind width_vs_tcount.pdf and bits_vs_tcount.pdf.
#   ./run_trace.sh <qubits>_<gates> <qasmdir> <outdir>
# e.g. ./run_trace.sh 20_700 ROOT/ldd/qasm ROOT/amend/trace runs
# <qasmdir>/clifford_T_circuit_20_700.qasm. Exact LIMDD (qisq2, low), one
# worker, the sweep's table sizes and a 20-minute timeout: the flags of the
# 25 September trace runs. -j is added only to keep the run's statistics.
# The ten traces of the 2026-10-01 revision are 20_700 20_750 30_700 30_750
# 40_700 40_800 50_700 60_700 70_700 80_700.
#
# Runners, from the environment, as for run_arm.sh: QSY_W1 (LIMDD_PAULI_WORDS=1)
# up to 64 qubits, QSY_W2 (LIMDD_PAULI_WORDS=2) above.
#
# Writes <outdir>/<id>.csv (the trace make_theorem_figs.py reads), .json, .out
# and .status. An existing .status is printed and kept.
set -u
[ $# -eq 3 ] || { echo "usage: $0 <qubits>_<gates> <qasmdir> <outdir>" >&2; exit 2; }
c=$1; f=$2/clifford_T_circuit_$c.qasm; T=$3; nq=${c%%_*}
[ -f "$f" ] || { echo "missing $f" >&2; exit 2; }
if [ "$nq" -gt 64 ]; then B=${QSY_W2:?set QSY_W2 to run_qasm_on_qmdd of the LIMDD_PAULI_WORDS=2 build}
else                      B=${QSY_W1:?set QSY_W1 to run_qasm_on_qmdd of the LIMDD_PAULI_WORDS=1 build}; fi
mkdir -p "$T"
[ -f "$T/$c.status" ] && { cat "$T/$c.status"; exit 0; }
s=$(date +%s)
timeout 1200 "$B" -c -m --dd=limdd -e qisq2 -s low -w 1 --trace="$T/$c.csv" \
   --node-tab-size=25 --wgt-tab-size=24 -j "$T/$c.json" "$f" > "$T/$c.out" 2>&1; rc=$?
e=$(date +%s)
echo "trace $c rc=$rc wall=$((e-s))s rows=$(wc -l < "$T/$c.csv" | tr -d ' ') build=$(basename "$(dirname "$(dirname "$B")")")" | tee "$T/$c.status"
