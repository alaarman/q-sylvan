#!/bin/bash
# The experiment figures the exact-DD paper typesets from Q-Sylvan runs, made
# from a data root in the paper's file layout; optionally compared with, or
# copied into, a checkout of the paper.
#   ./make_paper_figures.sh ROOT OUTDIR [PAPERDIR [--install]]
#   ROOT      the stored data (REPRODUCE.md) or a repeat of the sweep
#             (sweep.py's OUT)
#   OUTDIR    receives everything regen_amendment.sh writes, and the paper's
#             figures in OUTDIR/figures/, laid out as PAPERDIR/figures/
#   PAPERDIR  a checkout of the paper: every figure of OUTDIR/figures/ is
#             compared with the paper's copy (checks/pdf_same.py, which
#             ignores the creation date) and reported "same" or "DIFF"; the
#             exit status is 1 if any differs
#   --install copies OUTDIR/figures/ over PAPERDIR/figures/ instead
# Python, as for regen_amendment.sh (REPRODUCE.md, "Python environments"):
#   PY         the plotting Python of Figures 6 and 7 (requirements-plot.txt)
#   PY_TCOUNT  the one of the T-count figure (requirements-plot-tcount.txt);
#              default PY. With other versions the figures show the same
#              data but are not byte-identical to the paper's.
# FLOAT_LIMDD  the float LIMDD of the LIMDD panels: abs (the default, and the
#              paper's), the absolute tolerance the float EVDD uses too
#              (run_arm.sh limdd_float, limdd_float_l2), so that the two
#              diagrams differ in nothing else; or rel, the relative merging
#              rule (limdd_float_rel, limdd_float_l2_rel; only where ROOT has
#              those runs, as the rerun of 19 Oct 2026 does). The other one,
#              where ROOT has it, goes to OUTDIR/figures_alt/.
#
# The figures, and where they come from:
#   figures/final_plots_mm-(in)correct2/pooled_*  Figure 6's EVDD panels and
#       legend, low and L2 floats; regen_amendment.sh's paper step
#   figures/limdd2/, figures/limdd2_low/  Figure 6's LIMDD panels and Figure 7
#       (LIMDD against EVDD), L2 and low floats; make_figures.py with the
#       float LIMDD arms FLOAT_LIMDD picks (abs: the paper step's)
#   figures/limdd/{width,bits}_vs_tcount.pdf  the T-count figure, from
#       ROOT/amend/trace_clean (else amend/trace), by PY_TCOUNT; the copies in
#       limdd2/ and limdd2_low/ are replaced by these
# The paper picks low or L2 with \iflowfloats in its macros.tex. The families
# drawn are those of the paper's text (no hidden shift); AMEND_FAMILIES
# overrides. The other PDFs of the paper's figures/limdd/ belong to text after
# an \endinput (the 2026-10-01 revision) and are not made here; nor are the
# figures of Section 6 that compare with SliQSim and Quasimodo (the paper's
# original artifact).
set -u
[ $# -ge 2 ] && [ $# -le 4 ] || { sed -n 2,40p "$0" >&2; exit 2; }
ROOT=$1; O=$2; PAPER=${3:-}; MODE=${4:-compare}
S=$(cd "$(dirname "$0")" && pwd)
PY=${PY:-python3}
PY_TCOUNT=${PY_TCOUNT:-$PY}
export AMEND_FAMILIES=${AMEND_FAMILIES:-random,grover,wstate,mirror,iqp,groveriter,qwalk,synth}
export PYTHONDONTWRITEBYTECODE=1
[ -d "$ROOT/amend/out_clean" ] || { echo "$ROOT/amend/out_clean is missing" >&2; exit 1; }
HAVE_REL=0; ls "$ROOT"/amend/out_clean/*_limdd_float_rel.status > /dev/null 2>&1 && HAVE_REL=1
FLOAT_LIMDD=${FLOAT_LIMDD:-abs}
case $FLOAT_LIMDD in
  rel) [ $HAVE_REL -eq 1 ] || { echo "FLOAT_LIMDD=rel, but $ROOT has no limdd_float_rel runs" >&2; exit 1; } ;;
  abs) ;;
  *) echo "FLOAT_LIMDD is rel or abs" >&2; exit 2 ;;
esac
mkdir -p "$O"
PY=$PY "$S/regen_amendment.sh" "$ROOT" "$O" > "$O/regen.log" 2>&1
rc=$?
[ $rc -eq 0 ] || echo "regen_amendment.sh exited $rc (see $O/regen.log and the .out files in $O)" >&2
TR=$ROOT/amend/trace; [ -d "$ROOT/amend/trace_clean" ] && TR=$ROOT/amend/trace_clean
"$PY_TCOUNT" "$S/make_theorem_figs.py" "$TR" "$O/tcount" > "$O/tcount.out" 2>&1 \
  || { echo "make_theorem_figs.py failed (see $O/tcount.out)" >&2; rc=1; }
F=$O/figures
rm -rf "$F" "$O/figures_alt"; mkdir -p "$F/final_plots_mm-(in)correct2" "$F/limdd2" "$F/limdd2_low" "$F/limdd"
missing=0
put() {  # put SRC DST: copy one figure, counting the ones regen did not make
  if [ -e "$1" ]; then cp "$1" "$2"; else echo "not made: $1" >&2; missing=$((missing + 1)); fi
}
for t in low l2; do
  for q in wgt_type_runtime wgt_final_nodecount wgt_nodecount legend; do
    put "$O/paper/final_plots_mm-(in)correct2/pooled_${q}_float_${t}_vs_algebraic_nw.pdf" "$F/final_plots_mm-(in)correct2/"
  done
done
# the LIMDD panels: under the relative rule made here, under the absolute
# rule taken from the paper step (the same make_figures.py call with the
# absolute-rule arms)
QA=(--qasm "$ROOT/ldd/qasm" --qasm "$ROOT/amend/qasm"); PN=("$ROOT/amend/out_clean")
if [ -d "$ROOT/pilot/out_clean" ]; then
  PN+=("$ROOT/pilot/out_clean")
  for d in "$ROOT"/pilot/qasm/*/; do QA+=(--qasm "$d"); done
  [ -f "$ROOT/pilot/reference_probs.csv" ] && QA+=(--reference "$ROOT/pilot/reference_probs.csv")
fi
if [ $HAVE_REL -eq 1 ]; then
  for v in "limdd2 float_l2 limdd_float_l2_rel" "limdd2_low float_low limdd_float_rel"; do
    read -r d fe fl <<< "$v"
    "$PY" "$S/make_figures.py" -v --no-norm-label --pair-wrong-colour black --float-evdd "$fe" --float-limdd "$fl" \
      "${QA[@]}" "${PN[@]}" "$O/paper_rel/$d" > "$O/paper_rel_$d.out" 2>&1 \
      || { echo "make_figures.py failed for $fl (see $O/paper_rel_$d.out)" >&2; rc=1; }
  done
fi
LIMDD_PANELS="limdd_float_vs_algebraic_simulation_time limdd_float_vs_algebraic_final_nodes
              limdd_float_vs_algebraic_max_nodes limdd_float_vs_algebraic_legend
              limdd_vs_evdd_simulation_time limdd_vs_evdd_final_nodes limdd_vs_evdd_max_nodes
              limdd_vs_evdd_float_simulation_time limdd_vs_evdd_float_final_nodes
              limdd_vs_evdd_float_max_nodes panel_legend"
limdd_panels() {  # limdd_panels SRCROOT DSTROOT: the LIMDD folders, T-count copies included
  for d in limdd2 limdd2_low; do
    mkdir -p "$2/$d"
    for p in $LIMDD_PANELS; do put "$1/$d/$p.pdf" "$2/$d/"; done
    for p in width_vs_tcount bits_vs_tcount; do put "$O/tcount/$p.pdf" "$2/$d/"; done
  done
}
if [ "$FLOAT_LIMDD" = rel ]; then
  limdd_panels "$O/paper_rel" "$F"
  limdd_panels "$O/paper" "$O/figures_alt"           # the absolute rule, for comparison
else
  limdd_panels "$O/paper" "$F"
  [ $HAVE_REL -eq 1 ] && limdd_panels "$O/paper_rel" "$O/figures_alt"
fi
for p in width_vs_tcount bits_vs_tcount; do put "$O/tcount/$p.pdf" "$F/limdd/"; done
[ $missing -eq 0 ] || rc=1
n=$(find "$F" -name '*.pdf' | wc -l | tr -d ' ')
echo "$n figures in $F (float LIMDD: $FLOAT_LIMDD rule)"
[ -d "$O/figures_alt" ] && echo "the other rule's LIMDD panels in $O/figures_alt"
[ -n "$PAPER" ] || exit $rc
if [ "$MODE" = "--install" ]; then
  (cd "$F" && find . -name '*.pdf' | while read -r p; do cp "$p" "$PAPER/figures/$p"; done)
  echo "copied into $PAPER/figures/"
  exit $rc
fi
same=0; diff=0; absent=0
while read -r p; do
  q=$PAPER/figures/${p#./}
  if [ ! -e "$q" ]; then echo "not in the paper: ${p#./}"; absent=$((absent + 1)); continue; fi
  r=$("$PY" "$S/checks/pdf_same.py" "$F/${p#./}" "$q")
  if [ "$r" = same ]; then same=$((same + 1)); else diff=$((diff + 1)); echo "${p#./}: $r"; fi
done < <(cd "$F" && find . -name '*.pdf' | sort)
echo "compared with $PAPER/figures: $same same, $diff differ, $absent not in the paper"
[ $diff -eq 0 ] && [ $absent -eq 0 ] && exit $rc || exit 1
