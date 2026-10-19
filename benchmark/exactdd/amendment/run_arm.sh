#!/bin/bash
# Run one arm of the LIMDD amendment on one circuit.
#   ./run_arm.sh <arm> <circuit.qasm> <outdir> [timeout_s] [limtab]
# qisq2_low   = EVDD  qisq2 low                          (paper "algebraic")
# float_low   = EVDD  float low --merging=abs -t 1e-14   (the amendment's float EVDD)
# float_l2    = EVDD  float L2  --merging=abs -t 1e-14   (paper "float", Section 6)
# limdd_qisq2 = LIMDD qisq2 low
# limdd_float = LIMDD float low --merging=abs -t 1e-14
# limdd_float_l2 = LIMDD float L2 --merging=abs -t 1e-14 (every node read at norm 1)
# limdd_float_rel, limdd_float_l2_rel = the same two under the relative rule
#               (--merging=hybrid --rel-tol 1e-12 --zero-tol 0), the LIMDD's own
#               default; the absolute rule merges small weights that differ
#               relatively, which the LIMDD's labels cannot absorb
# limdd_heur  = LIMDD float low, no canonical form while simulating (--dd=limdd-heur;
#               the diagram is canonicalised once at the end, inside the timing)
# limdd_lazy  = LIMDD float low, lazy canonical form (--canon=lazy: one node per
#               class, the coset search only where a node may share its class)
# limdd_qisq2_lazy = LIMDD qisq2 low, lazy canonical form
# bqd_scalar, bqd_x, bqd_pauli = BQD qisq2 in that label family, zero rule SM
#               (the default); circuits on more than 63 qubits get status SKIP
# Every arm adds -c -m -w 1 --node-tab-size=25 --wgt-tab-size=24; the timeout
# defaults to 30 minutes. limtab adds --lim-tab-size=<limtab> to the LIMDD arms
# (no run of the 2026-10-01 revision used it).
#
# Runners, from the environment:
#   QSY_W1  run_qasm_on_qmdd of a Release build with LIMDD_PAULI_WORDS=1
#   QSY_W2  the same with LIMDD_PAULI_WORDS=2; the LIMDD arms use it above
#           64 qubits, the EVDD arms never do
#
# Answer checks above the 24 qubits that are decoded, also from the
# environment (see amp_check.py):
#   QSY_AMP_SAMPLE=K  the qmdd and limdd arms sample K basis states of their
#                     final state (--amp-sample); run the reference arm so
#   QSY_AMP_AT=DIR    every arm evaluates its final state at the points of
#                     DIR/<circ>_<QSY_AMP_REF>.json (--amp-at), when it exists;
#                     QSY_AMP_REF defaults to qisq2_low
#
# Derived from the script that wrote ldd/out, with the same runner flags and
# two changes:
#  (1) the qubit count is the sum of the qreg sizes in the file, not a field of
#      the file name (grover_n10_it1, wclif_128_n128, ... have no 4th field, so
#      the original would never pick the 2-word build for them);
#  (2) the status line also carries t_count and final_width.
# Writes <outdir>/<circ>_<arm>.json, .out (tail of the output, the command and
# the exit code) and .status, one line:
#   arm circuit status time final peak norm p_top t_count final_width
# with status OK, TIMEOUT, LIM_FULL, STAB_FULL, TAB_FULL, ERR(rc) or NOJSON.
# An existing .status is printed and kept, so a sweep can be restarted.
set -u
[ $# -ge 3 ] || { echo "usage: $0 <arm> <circuit.qasm> <outdir> [timeout_s] [limtab]" >&2; exit 2; }
arm=$1; C=$2; OUT=$3; TO=${4:-1800}; LIMTAB=${5:-}
unset QSYLVAN_REL_TOL QSYLVAN_ZERO_TOL      # these silently defeat --merging=abs
[ -f "$C" ] || { echo "missing $C" >&2; exit 2; }
base=$(basename "$C" .qasm)
nq=$(grep -oE 'qreg [a-zA-Z_0-9]+\[[0-9]+\]' "$C" | grep -oE '\[[0-9]+\]' | tr -d '[]' | paste -sd+ - | bc)
[ -n "$nq" ] || { echo "no qreg in $C" >&2; exit 2; }
SIZES="-c -m -w 1 --node-tab-size=25 --wgt-tab-size=24"
[ -n "$LIMTAB" ] && case $arm in limdd_*) SIZES="$SIZES --lim-tab-size=$LIMTAB";; esac
A=${QSY_W1:?set QSY_W1 to run_qasm_on_qmdd of the LIMDD_PAULI_WORDS=1 build}
W2() { echo "${QSY_W2:?set QSY_W2 to run_qasm_on_qmdd of the LIMDD_PAULI_WORDS=2 build}"; }
case $arm in
  float_low)   F="--dd=qmdd -e float -s low --merging=abs -t 1e-14" ;;
  float_l2)    F="--dd=qmdd -e float -s l2 --merging=abs -t 1e-14" ;;
  qisq2_low)   F="--dd=qmdd -e qisq2 -s low" ;;
  limdd_float) F="--dd=limdd -e float -s low --merging=abs -t 1e-14"; [ "$nq" -gt 64 ] && A=$(W2) ;;
  limdd_float_l2) F="--dd=limdd -e float -s l2 --merging=abs -t 1e-14"; [ "$nq" -gt 64 ] && A=$(W2) ;;
  limdd_float_rel)    F="--dd=limdd -e float -s low --merging=hybrid --rel-tol 1e-12 --zero-tol 0"; [ "$nq" -gt 64 ] && A=$(W2) ;;
  limdd_float_l2_rel) F="--dd=limdd -e float -s l2 --merging=hybrid --rel-tol 1e-12 --zero-tol 0"; [ "$nq" -gt 64 ] && A=$(W2) ;;
  limdd_heur)  F="--dd=limdd-heur -e float -s low --merging=abs -t 1e-14"; [ "$nq" -gt 64 ] && A=$(W2) ;;
  limdd_lazy)  F="--dd=limdd --canon=lazy -e float -s low --merging=abs -t 1e-14"; [ "$nq" -gt 64 ] && A=$(W2) ;;
  limdd_qisq2_lazy) F="--dd=limdd --canon=lazy -e qisq2 -s low"; [ "$nq" -gt 64 ] && A=$(W2) ;;
  limdd_qisq2) F="--dd=limdd -e qisq2 -s low";              [ "$nq" -gt 64 ] && A=$(W2) ;;
  bqd_scalar|bqd_x|bqd_pauli) F="--dd=bqd --bqd-family=${arm#bqd_} -e qisq2 -s low" ;;
  *) echo "unknown arm $arm" >&2; exit 2 ;;
esac
[ -n "$A" ] || exit 2
mkdir -p "$OUT"
J=$OUT/${base}_${arm}.json
ST=$OUT/${base}_${arm}.status
LOG=$OUT/${base}_${arm}.out
if [ -f "$ST" ]; then cat "$ST"; exit 0; fi
case $arm in bqd_*) if [ "$nq" -gt 63 ]; then
  printf '%-12s %-26s %-9s\n' "$arm" "$base" "SKIP" | tee "$ST"; exit 0; fi ;; esac
AMP=""
if [ -n "${QSY_AMP_SAMPLE:-}" ]; then
  case $arm in bqd_*) ;; *) AMP="--amp-sample=$QSY_AMP_SAMPLE" ;; esac
elif [ -n "${QSY_AMP_AT:-}" ] && [ -f "$QSY_AMP_AT/${base}_${QSY_AMP_REF:-qisq2_low}.json" ]; then
  AMP="--amp-at=$QSY_AMP_AT/${base}_${QSY_AMP_REF:-qisq2_low}.json"
fi
s=$(date +%s)
out=$(timeout "$TO" "$A" $F $SIZES $AMP -j "$J" "$C" 2>&1); rc=$?
e=$(date +%s)
printf '%s\n' "$out" | tail -40 > "$LOG"
echo "cmd: $A ${F} ${SIZES} ${AMP} -j $J $C" >> "$LOG"; echo "rc=$rc wall=$((e-s))s nq=$nq" >> "$LOG"
if   [ $rc -eq 124 ];                     then st=TIMEOUT
elif printf '%s\n' "$out" | grep -qi "LIM table is full";  then st=LIM_FULL
elif printf '%s\n' "$out" | grep -qi "stabiliser table is full"; then st=STAB_FULL
elif printf '%s\n' "$out" | grep -qiE "table full|table is full"; then st=TAB_FULL
elif [ $rc -ne 0 ];                       then st="ERR($rc)"
elif [ ! -s "$J" ];                       then st=NOJSON
else st=OK; fi
g(){ [ -s "$J" ] && python3 -c "import json,sys;d=json.load(open(sys.argv[1]))['statistics'];print(d.get(sys.argv[2],''))" "$J" "$1" 2>/dev/null; }
printf '%-12s %-26s %-9s %10s %10s %10s %12s %s %s %s\n' "$arm" "$base" "$st" \
  "$(g simulation_time)" "$(g final_nodes)" "$(g max_nodes)" "$(g norm)" "$(g first_qubit_measurement_prob)" "$(g t_count)" "$(g final_width)" | tee "$ST"
