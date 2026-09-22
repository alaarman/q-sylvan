#!/bin/bash
# Run one benchmark circuit under one of the four exact-DD configurations.
#   ./run_exactdd.sh <A|B|C|D> <circuit.qasm> <outdir> [binary] [timeout_s]
# A = QMDD float  L2  absolute     (paper "float")
# B = QMDD qisq2  low              (paper "algebraic")
# C = LIMDD float low absolute     (explicitly NOT the hybrid rule)
# D = LIMDD qisq2 low              (weights exact; --merging=abs only labels the JSON)
set -u
CFG=$1; QASM=$2; OUT=$3
BIN=${4:-$(dirname "$0")/../../build-clean/qasm/run_qasm_on_qmdd}
TMO=${5:-1800}
unset QSYLVAN_REL_TOL QSYLVAN_ZERO_TOL      # these silently defeat --merging=abs
mkdir -p "$OUT"
name=$(basename "$QASM" .qasm)
json="$OUT/${name}_${CFG}.json"
case $CFG in
  A) ARGS="--dd=qmdd  -e float -s l2  --merging=abs -t 1e-14 -c -m" ;;
  B) ARGS="--dd=qmdd  -e qisq2 -s low                        -c -m" ;;
  C) ARGS="--dd=limdd -e float -s low --merging=abs -t 1e-14 -c"    ;;
  D) ARGS="--dd=limdd -e qisq2 -s low --merging=abs           -c"   ;;
  *) echo "config must be A|B|C|D" >&2; exit 2 ;;
esac
start=$(date +%s)
timeout "$TMO" "$BIN" $ARGS -j "$json" "$QASM" > "$OUT/${name}_${CFG}.log" 2>&1
rc=$?
echo "$name $CFG rc=$rc wall=$(( $(date +%s) - start ))s $json"
exit $rc
