#!/bin/bash
# Regenerate every figure and number of experiments_limdd_amendment.tex
# (revision of 2026-10-01) from the stored runs. No simulator runs.
#   ./regen_amendment.sh ROOT OUTDIR
#   ROOT    the data directory, laid out as README.md describes
#   OUTDIR  where the PDFs, merged_records.csv and one .out per script go
# Python from $PY, else python3; it needs matplotlib (README.md has the
# versions that wrote the paper's PDFs):
#   PY=/path/to/venv/bin/python ./regen_amendment.sh ROOT OUTDIR
# AMEND_FAMILIES restricts the figures and the numbers to some of the families
# random, grover, wstate and hshift, e.g. the families of Section 6:
#   AMEND_FAMILIES=random,grover,wstate ./regen_amendment.sh ROOT OUTDIR
# Sources, later overriding earlier: LIMDD from ldd/out (random circuits,
# 70/80 qubits on the corrected build, and the six exact timeouts of the
# uncorrected build), amend/out_rnd (random circuits up to 64 qubits on the
# corrected build), amend/out (the other families). EVDD from the surviving
# records of the 25-26 September sweep: exactdd/final.log, gw/out, gw/out2.
# When ROOT has amend/out_l2 and amend/out_low (the float EVDD rerun of
# 4 October, 'L2' and 'low', on the binary of the LIMDD runs), Figures 6 to 8
# are also made from those, into final_plots_mm-(in)correct2 as
# wgt_*_float_l2_vs_algebraic_nw.pdf and wgt_*_float_low_vs_algebraic_nw.pdf,
# and so is an alternative Figure 9 into limdd2/, with the L2 float EVDD and,
# when ROOT has amend/out_limdd_l2 (the LIMDD with every node read at norm 1,
# run_arm.sh limdd_float_l2), the L2 float LIMDD. No rerun is read by any
# other step.
set -u
[ $# -eq 2 ] || { echo "usage: $0 ROOT OUTDIR" >&2; exit 2; }
ROOT=$1; O=$2
PY=${PY:-python3}
export PYTHONDONTWRITEBYTECODE=1           # no __pycache__ in the source tree
S=$(cd "$(dirname "$0")" && pwd)
V=$S/checks
miss=0
for d in ldd/out ldd/qasm amend/out_rnd amend/out amend/qasm amend/trace amend/span/out2 \
         exactdd/final.log gw/out gw/out2; do
  [ -e "$ROOT/$d" ] || { echo "missing $ROOT/$d" >&2; miss=1; }
done
[ $miss -eq 0 ] || exit 1
OLD=(--old "$ROOT/exactdd/final.log" --old "$ROOT/gw/out" --old "$ROOT/gw/out2")
Q=(--qasm "$ROOT/ldd/qasm" --qasm "$ROOT/amend/qasm")
NEW=("$ROOT/ldd/out" "$ROOT/amend/out_rnd" "$ROOT/amend/out")
mkdir -p "$O"
fail=()
step() { local log=$1; shift; "$@" > "$O/$log" 2>&1 || fail+=("$log"); }
step make_figures.out      "$PY" "$S/make_figures.py" -v "${OLD[@]}" "${Q[@]}" --dump "$O/merged_records.csv" "${NEW[@]}" "$O"
step amend_numbers.out     "$PY" "$S/amend_numbers.py" "${OLD[@]}" "${Q[@]}" "${NEW[@]}"
step amend_stats.out       "$PY" "$S/amend_stats.py" "${OLD[@]}" "${Q[@]}" "${NEW[@]}"
step make_theorem_figs.out "$PY" "$S/make_theorem_figs.py" "$ROOT/amend/trace" "$O"
S6=$O/'final_plots_mm-(in)correct2'
step make_section6_figures.out "$PY" "$S/make_section6_figures.py" "${OLD[@]}" "${Q[@]}" "${NEW[@]}" "$S6"
if [ -d "$ROOT/amend/out_l2" ] && [ -d "$ROOT/amend/out_low" ]; then
  step make_section6_figures_l2.out  "$PY" "$S/make_section6_figures.py" --float float_l2 --tag l2 \
       "${OLD[@]}" "${Q[@]}" "${NEW[@]}" "$ROOT/amend/out_l2" "$S6"
  step make_section6_figures_low.out "$PY" "$S/make_section6_figures.py" --float float_low --tag low \
       "${OLD[@]}" "${Q[@]}" "${NEW[@]}" "$ROOT/amend/out_low" "$S6"
  LLOPT=(); LLDIR=()
  if [ -d "$ROOT/amend/out_limdd_l2" ]; then
    LLOPT=(--float-limdd limdd_float_l2); LLDIR=("$ROOT/amend/out_limdd_l2")
  fi
  step make_figures_l2.out "$PY" "$S/make_figures.py" -v --float-evdd float_l2 ${LLOPT[@]+"${LLOPT[@]}"} \
       "${OLD[@]}" "${Q[@]}" "${NEW[@]}" "$ROOT/amend/out_l2" ${LLDIR[@]+"${LLDIR[@]}"} "$O/limdd2"
fi
# checks behind the 2026-10-01 revision
step evdd_float_versions.out "$PY" "$V/evdd_float_versions.py" "$ROOT"
step wide_old_new.out        "$PY" "$V/wide_old_new.py" "$ROOT"
step float_evdd_rules.out    "$PY" "$V/float_evdd_rules.py" "$ROOT"
step correction_before_after.out "$PY" "$V/correction_before_after.py" "$ROOT"
# the 15 PDFs written above, by name, so that stale PDFs in a reused OUTDIR are
# not counted (in the order of a glob in the C locale)
PDFS=()
for p in bits_vs_tcount \
         evdd_float_vs_algebraic_final_nodes evdd_float_vs_algebraic_max_nodes \
         evdd_float_vs_algebraic_simulation_time \
         limdd_float_vs_algebraic_final_nodes limdd_float_vs_algebraic_max_nodes \
         limdd_float_vs_algebraic_simulation_time \
         limdd_vs_evdd_final_nodes limdd_vs_evdd_float_final_nodes \
         limdd_vs_evdd_float_max_nodes limdd_vs_evdd_float_simulation_time \
         limdd_vs_evdd_max_nodes limdd_vs_evdd_simulation_time \
         panel_legend width_vs_tcount; do
  PDFS+=("$O/$p.pdf")
done
step count_markers.out       "$PY" "$V/count_markers.py" "${PDFS[@]}"
ls "$O"
if [ ${#fail[@]} -gt 0 ]; then echo "FAILED (see these files in $O): ${fail[*]}" >&2; exit 1; fi
