# Scripts of the LIMDD amendment

The scripts here make every figure and the measured numbers of
`experiments_limdd_amendment.tex`, the LIMDD section added to the exact-DD
paper, as revised on 2026-10-01. That section was written by Claude (AI) and
says so in its header; so were these scripts. Neither has been checked by an
author yet.

**The run data of the 2026-10-01 revision is not in this repository.** It is
kept by the authors and is available from them on request. Every script that
reads it takes the root of a copy of it, called ROOT below, as an argument.
Only `run_arm.sh` and `run_trace.sh` run the simulator, and they are what made
that data. Regenerating the figures and numbers from it takes a few seconds and
no simulator.

| file | what it is |
|---|---|
| `regen_amendment.sh` | the driver: runs the nine Python scripts below on ROOT |
| `make_figures.py` | the scatter panels, `panel_legend.pdf` and `merged_records.csv` |
| `amend_numbers.py`, `amend_stats.py` | the counts, medians and ratios of the prose |
| `make_theorem_figs.py` | the two T-count figures |
| `checks/` | five checks behind particular sentences and the open items |
| `run_arm.sh` | one run of one arm on one circuit (the LIMDD sweep) |
| `run_trace.sh` | one per-gate trace (the T-count figures) |

## What each script makes

All four arms: `qisq2_low` (EVDD, algebraic), `float_low` (EVDD, float),
`limdd_qisq2` (LIMDD, algebraic), `limdd_float` (LIMDD, float). "Random" is the
49 random Clifford+T circuits, "structured" Grover, W state and hidden shift.

| script | output | where it is used in the amendment |
|---|---|---|
| `make_figures.py` | `limdd_vs_evdd_{simulation_time,final_nodes,max_nodes}.pdf` | `fig:limdd-vs-evdd`, top row |
| | `limdd_vs_evdd_float_{...}.pdf` | `fig:limdd-vs-evdd`, bottom row |
| | `panel_legend.pdf` | the legend of `fig:limdd-vs-evdd` |
| | `limdd_float_vs_algebraic_{...}.pdf` | `fig:limdd-float-vs-algebraic` |
| | `evdd_float_vs_algebraic_{...}.pdf` | not typeset; its point counts are in open item 1 |
| | `merged_records.csv` | every merged record with its source. The values quoted for single circuits are rows of it: `tab:limdd-representative`, 4.4 s against 97.4 s, 443264 nodes against 411, the 32068 and 3628 node maxima, the Ising norms 0.84 and 0.68 |
| | `make_figures.out` | points and red points per panel (62; 57 with 5 wrong; 56 with 2 wrong; 60 with 15 wrong), and the agreement of the exact arms on 33 of 33 (Setup) |
| `amend_numbers.py` | `amend_numbers.out` | Results: 33 of 33 smaller, medians 189 and 150, faster on 19 of 33, 43 of 49 against 33, the six exact LIMDD timeouts; structured 29 of 29, medians 1.06 and 2.0; float LIMDD 33 of 49 and its failures (9 node table, 6 coefficient table, 1 unnamed); float against exact (7.0, up to 965, 2.7; EVDD 1.4 and 0.9); the lower row (both floats wrong on 5, EVDD alone on 8 more, the ratios 1.0e4 on 70_900 and 1.6e4 on 40_750). Open item 7: the sweep-to-sweep comparison (52 circuits, median 18, 1.3 to 756; float LIMDD differs on 23 of 49) |
| `amend_stats.py` | `amend_stats.out` | Results: the 28 random circuits on which all four arms finish (medians 193 and 137 exact, 96 and 28 float, fewer final nodes on 26 of 28); wrong float runs per diagram |
| `make_theorem_figs.py` | `width_vs_tcount.pdf`, `bits_vs_tcount.pdf`, `make_theorem_figs.out` | `fig:limdd-tcount`; the .out gives width, bits and slack per trace |
| `checks/evdd_float_versions.py` | `evdd_float_versions.out` | Setup: the float EVDD rerun on the code of the LIMDD runs differs on 6 of 10; the caption of `tab:limdd-representative` (338869, 365296, 0.9993); open item 7 |
| `checks/wide_old_new.py` | `wide_old_new.out` | open item 6: the LIMDD cells above 64 qubits before and after the correction |
| `checks/float_evdd_rules.py` | `float_evdd_rules.out` | Results: float EVDD wrong on 16 of 38, nine on the measurement, seven of them the zero vector, five with the norm off, two without an exact counterpart; open item 8 |
| `checks/correction_before_after.py` | `correction_before_after.out` | Setup and open item 5: the 29 random circuits that the exact LIMDD completes both before and after the correction (same node counts on 29, the corrected code slower on 28, median 1.20), and the six it timed out on before it. The records of `ROOT/ldd/out` up to 64 qubits are the uncorrected pass, which the figures override with `ROOT/amend/out_rnd` |
| `checks/count_markers.py` | `count_markers.out` | open item 1: markers actually drawn in each PDF (62, 57, 60, 56) |

Numbers of the amendment that none of these scripts prints:

- the coefficient ranges (78 to 85 decades, 20 to 26, at most 58 for the EVDD):
  the `[span] looked-up` lines of `ROOT/amend/span/out*/*.out`, from a build
  with a magnitude counter in `cmap_find_or_put`. That patch is kept with the
  data and is not part of q-sylvan;
- the 1406 s after which the exact LIMDD fills the coefficient table on the
  30-qubit Ising chain: the last line of
  `ROOT/amend/out/ising_n30_s3_limdd_qisq2.out`;
- the 11-node, norm-0 float EVDD on the 20-qubit, 750-gate circuit in the
  stored measurements of Section 6: the paper's own results.

## Data layout

ROOT must contain the following. Later sources override earlier ones, record
by record, as listed in `regen_amendment.sh`.

| path under ROOT | contents | made by |
|---|---|---|
| `ldd/qasm/` | the 49 random circuits `clifford_T_circuit_<qubits>_<gates>.qasm` | the paper's artifact, identical by blob hash to git ref `quist/qisq2_tests`, `benchmark/qasm_clifford_T/` (see `../README.md`) |
| `amend/qasm/` | the other 42 circuits: `grover_n*_it1`, `w-state_*_n*`, `wclif_*_n*`, `hidden-shift_n*`, `adder_n*`, `ising_n*_s3` | the generators in `..`, see below |
| `ldd/out/` | LIMDD records of the random circuits, `<circ>_<arm>.status` and `.json`. Used: 70 and 80 qubits (corrected builds) and the six exact timeouts (20 qubits from 800 gates, 30 qubits 900 gates; uncorrected build). The rest is overridden by `amend/out_rnd/` | an earlier version of `run_arm.sh`, whose status lines lack the last two columns |
| `amend/out_rnd/` | LIMDD records of the random circuits up to 64 qubits, corrected build | `run_arm.sh` |
| `amend/out/` | LIMDD records of the five other families, corrected builds | `run_arm.sh` |
| `amend/trace/` | `<qubits>_<gates>.csv`, one per-gate trace each, for 20_700 20_750 30_700 30_750 40_700 40_800 50_700 60_700 70_700 80_700 | `run_trace.sh` |
| `amend/span/out2/` | `qmdd_float_*_w24.json`, the float EVDD rerun on ten random circuits with the code of the LIMDD runs (read by `evdd_float_versions.py`); with `amend/span/out/`, the coefficient ranges | the instrumented build, flags as `run_arm.sh float_low` |
| `exactdd/final.log` | the September records of the random circuits, all four arms, one status line per run (`arm circuit status time final peak norm p_top`), later lines win | the 25-26 September sweep, same EVDD flags as `run_arm.sh` |
| `gw/out/`, `gw/out2/` | the September records of the other families, `<circ>_<arm>.status` with `circuit arm status time final peak norm t_count final_width` (no top-qubit probability) | the 25-26 September sweep |

Only the EVDD arms (`qisq2_low`, `float_low`) are taken from the September
records; their LIMDD records there are superseded by the October runs, and
those above 64 qubits are wrong (64-bit qubit masks). To use EVDD runs made
with `run_arm.sh` instead, add their directory to `NEW` in
`regen_amendment.sh`: a directory there overrides every September record of
the same circuit and arm.

`ROOT/amend/regen_out/` holds the outputs of the regeneration of 2026-10-01,
the reference for the next section. It is not an input.

## Regenerating the figures and numbers

```
PY=/path/to/venv-plot/bin/python ./regen_amendment.sh "$ROOT" "$OUTDIR"
```

OUTDIR receives the 15 PDFs, `merged_records.csv` and one `.out` per script.
The PDFs the paper includes from `figures/limdd/` are the 12 that are not
`evdd_float_vs_algebraic_*`. The script stops if a path of the table above is
missing, and exits non-zero, naming the `.out` file, if a script fails.

`AMEND_FAMILIES` restricts the scatter panels, their legend and the numbers of
`amend_numbers.py` to some of the families `random`, `grover`, `wstate` and
`hshift` (unset: all four). Section 6 of the paper has no hidden shift, so the
figures that sit next to its own are made with

```
AMEND_FAMILIES=random,grover,wstate PY=... ./regen_amendment.sh "$ROOT" "$OUTDIR"
```

The Grover and W-state circuits here are generated by `../gen_grover.py`,
`../gen_wstate_ancilla.py` and `../gen_wstate_clifford.py`. They are not the
Grover and W-state circuits of Section 6, whose panels have 99 and 7 points
against the 12 and 11 circuits here; only the 49 random circuits are the same.

`make_section6_figures.py` remakes Figures 6 to 8 of the paper (EVDD, float
against algebraic, one panel per family and quantity) from the same September
EVDD records, into `OUTDIR/final_plots_mm-(in)correct2/` with the file names of
the paper's `figures/final_plots_mm-(in)correct/`. They differ from the
originals in ways the text has to say: the float arm is 'low' (the September
sweep has no 'L2' run), the Grover and W-state circuits are ours, the runs are
those of the Apple M1, and the correctness rule is that of `make_figures.py`
(the Grover and W-state records carry no probability, so those float runs are
judged by their norm).

In every scatter panel a run that did not finish (a timeout or a full table)
sits on a dashed line in a shaded band beyond the data: on the right when the
x arm failed, on top when the y arm failed, in the corner, with a count, when
both did. A run that was not made (SKIPPED_PREDICTED) is not drawn. In the
lower row of the LIMDD figure the fill marks a wrong float LIMDD and an orange
edge a wrong float EVDD.

### Python versions

The paper's PDFs were written by CPython 3.11.15 (macOS, arm64) with

```
matplotlib==3.11.2 numpy==2.4.6 contourpy==1.3.3 cycler==0.12.1
fonttools==4.65.0 kiwisolver==1.5.1 packaging==26.3 pillow==12.3.0
pyparsing==3.3.3 python-dateutil==2.9.0.post0 six==1.17.0
```

and nothing else installed (FreeType 2.14.3, as the matplotlib wheel ships
it). For example:

```
python3.11 -m venv venv-plot
venv-plot/bin/pip install matplotlib==3.11.2 numpy==2.4.6 contourpy==1.3.3 cycler==0.12.1 \
    fonttools==4.65.0 kiwisolver==1.5.1 packaging==26.3 pillow==12.3.0 \
    pyparsing==3.3.3 python-dateutil==2.9.0.post0 six==1.17.0
```

Other versions give the same numbers, but the PDFs may differ.
`make_figures.py`, `amend_numbers.py`, `amend_stats.py` and the checks
`wide_old_new.py`, `float_evdd_rules.py` and `correction_before_after.py` load
`make_figures.py`, and `make_theorem_figs.py` imports matplotlib itself, so
these need the versions above. `evdd_float_versions.py` and `count_markers.py`
need only the standard library, as does the `python3` that `run_arm.sh` calls.

### What a correct run looks like

On the data of 2026-10-01 with these versions, every PDF is identical byte
for byte to the one in the paper and to the one in `ROOT/amend/regen_out/`,
except its `/CreationDate`. Set `SOURCE_DATE_EPOCH` (e.g. to 0) to fix that
date and make two runs byte-identical.

Every text output is identical to the one in `ROOT/amend/regen_out/` except
in these four files:

- `make_figures.out`: the line that names `merged_records.csv`, which the
  stored file gives with its absolute path and this script by name only.
- `amend_stats.out`: the keys of the `statuses` dicts. They came out in an
  arbitrary order, which `amend_stats.py` now fixes by sorting the circuits.
  The counts are the same.
- `count_markers.out`: every line. The stored file was written by an earlier
  `count_markers.py` that counted only the marker XObjects. The counts of
  open item 1 are those of this version.
- `correction_before_after.out`: not stored; the check is newer than the
  stored outputs. It ends with "finished in both passes: 29; same final and
  peak nodes on 29; corrected slower on 28; median time after/before 1.203".

## Rerunning the simulations

This replaces the data, and takes hours: each cell has a 30-minute limit, and
the October runs went at most two at a time, on an Apple M1 Pro.

### Builds

The LIMDD runs of the 2026-10-01 revision used q-sylvan `8f4c3cf` plus a
correction for circuits on more than 64 qubits (a patch kept with the data).
The correction has since been committed, not byte for byte, as `ab90fa9`
("limdd: gates and limdd_eval reach every qubit a Pauli word holds"), so build
`ab90fa9` or later. Two Release builds of the runner are needed, one per Pauli
word count. On macOS the runner needs argp (`brew install argp-standalone`),
and `qasm/CMakeLists.txt` finds it with `cmake/FindArgp.cmake`. That module
is not in this tree. The first line below takes it from Sylvan's commit
`1570700`, which is in this branch's history. From the top of the clone:

```
git show 1570700:cmake/FindArgp.cmake > cmake/FindArgp.cmake   # macOS only
cmake -S . -B build-w1 -DCMAKE_BUILD_TYPE=Release -DLIMDD_PAULI_WORDS=1
cmake --build build-w1 --target run_qasm_on_qmdd
cmake -S . -B build-w2 -DCMAKE_BUILD_TYPE=Release -DLIMDD_PAULI_WORDS=2
cmake --build build-w2 --target run_qasm_on_qmdd
export QSY_W1=$PWD/build-w1/qasm/run_qasm_on_qmdd
export QSY_W2=$PWD/build-w2/qasm/run_qasm_on_qmdd
```

The exact arms give the same node counts with either code version wherever both
were run. The float arms do not: their results depend on the order of the
arithmetic, and so on the code, for either diagram (Setup; open item 7).

### Runs

`run_arm.sh` and `run_trace.sh` take the runners from `QSY_W1` and `QSY_W2`,
and both need GNU `timeout` (coreutils). `run_arm.sh` also needs `bc`, and
`python3` for the status line. A cell whose `.status` exists is skipped, so an
interrupted sweep can be restarted; delete the `.status` of a cell to run it
again. From this directory, with ROOT set:

```
# LIMDD, random circuits up to 64 qubits (on 2026-10-01 limdd_qisq2 was not
# rerun on the six circuits it had timed out on before the correction)
for f in "$ROOT"/ldd/qasm/clifford_T_circuit_[2-6]0_*.qasm; do
  for arm in limdd_qisq2 limdd_float; do ./run_arm.sh $arm "$f" "$ROOT/amend/out_rnd"; done
done
# LIMDD, 70 and 80 qubits (QSY_W2 is picked by the qubit count)
for f in "$ROOT"/ldd/qasm/clifford_T_circuit_[78]0_*.qasm; do
  for arm in limdd_qisq2 limdd_float; do ./run_arm.sh $arm "$f" "$ROOT/ldd/out"; done
done
# LIMDD, the other families
for f in "$ROOT"/amend/qasm/*.qasm; do
  for arm in limdd_qisq2 limdd_float; do ./run_arm.sh $arm "$f" "$ROOT/amend/out"; done
done
# the traces of fig:limdd-tcount
for c in 20_700 20_750 30_700 30_750 40_700 40_800 50_700 60_700 70_700 80_700; do
  ./run_trace.sh $c "$ROOT/ldd/qasm" "$ROOT/amend/trace"
done
```

`run_arm.sh` also runs the EVDD arms with the flags of the September sweep
(`qisq2_low`, `float_low`). Its fifth argument, `limtab` (after `timeout_s`),
sets `--lim-tab-size` on the LIMDD arms; no run of this revision needed it.

### Circuits

`ROOT/ldd/qasm` is the paper's artifact, see `../README.md`. The 42 circuits
of `ROOT/amend/qasm` come out byte-identical from the generators in `..` with
CPython 3.11.15 and

```
pip install qiskit==2.5.2 numpy==2.4.6 scipy==1.17.1 rustworkx==0.18.1 dill==0.4.1 \
    stevedore==5.9.1 typing_extensions==4.16.0
```

The generators write to the current directory, so run them in an empty one,
QASM below (an absolute path), and compare it with `ROOT/amend/qasm`:

```
G=$(git rev-parse --show-toplevel)/benchmark/exactdd     # inside the clone
mkdir -p "$QASM" && cd "$QASM"
for n in 4 6 8 10 12 14 16 20 24 30 40 50; do python3 "$G/gen_grover.py" $n 1; done
for k in 4 8 16 32 64;     do python3 "$G/gen_wstate_ancilla.py" $k; done
for k in 4 8 16 32 64 128; do python3 "$G/gen_wstate_clifford.py" $k; done
for n in 8 12 16 20 24 30 40; do python3 "$G/gen_realistic.py" hidden-shift $n; done
for b in 4 6 8 10 12 16 20;   do python3 "$G/gen_realistic.py" adder $b; done   # 2b+1 qubits
for n in 10 16 20 30 40;      do python3 "$G/gen_realistic.py" ising $n --steps 3; done
diff -r "$QASM" "$ROOT/amend/qasm"
```
