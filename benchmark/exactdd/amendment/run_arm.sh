#!/bin/bash
# Run one arm of the LIMDD amendment on one circuit.
#   ./run_arm.sh <arm> <circuit.qasm> <outdir> [timeout_s] [limtab]
# qisq2_low   = EVDD  qisq2 low                          (paper "algebraic")
# float_low   = EVDD  float low --merging=abs -t 1e-14   (the amendment's float EVDD)
# float_l2    = EVDD  float L2  --merging=abs -t 1e-14   (paper "float", Section 6)
# limdd_qisq2 = LIMDD qisq2 low
# limdd_float = LIMDD float low --merging=abs -t 1e-14
# limdd_float_l2 = LIMDD float L2 --merging=abs -t 1e-14 (every node read at norm 1)
# Every arm adds -c -m -w 1 --node-tab-size=25 --wgt-tab-size=24; the timeout
# defaults to 30 minutes. limtab adds --lim-tab-size=<limtab> to the LIMDD arms
# (no run of the 2026-10-01 revision used it).
#
# Runners, from the environment:
#   QSY_W1  run_qasm_on_qmdd of a Release build with LIMDD_PAULI_WORDS=1
#   QSY_W2  the same with LIMDD_PAULI_WORDS=2; the LIMDD arms use it above
#           64 qubits, the EVDD arms never do
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
  limdd_qisq2) F="--dd=limdd -e qisq2 -s low";              [ "$nq" -gt 64 ] && A=$(W2) ;;
  *) echo "unknown arm $arm" >&2; exit 2 ;;
esac
[ -n "$A" ] || exit 2
mkdir -p "$OUT"
J=$OUT/${base}_${arm}.json
ST=$OUT/${base}_${arm}.status
LOG=$OUT/${base}_${arm}.out
if [ -f "$ST" ]; then cat "$ST"; exit 0; fi
s=$(date +%s)
out=$(timeout "$TO" "$A" $F $SIZES -j "$J" "$C" 2>&1); rc=$?
e=$(date +%s)
printf '%s\n' "$out" | tail -40 > "$LOG"
echo "cmd: $A ${F} ${SIZES} -j $J $C" >> "$LOG"; echo "rc=$rc wall=$((e-s))s nq=$nq" >> "$LOG"
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
