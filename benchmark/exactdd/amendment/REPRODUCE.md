# Reproducing the exact-DD paper's Q-Sylvan experiments

This file explains how to repeat the runs behind the exact-DD paper's
experiment figures and remake the figures from them. It covers three figures:

- **Figure 6** (`fig:qisq-vs-float_plots`): float against exact weights, EVDD
  and LIMDD;
- **Figure 7** (`fig:limdd-vs-evdd`): LIMDD against EVDD;
- **the T-count figure** (`fig:limdd-tcount`).

It also covers the counts in the text that come from those runs. Everything
here works from a clone of Q-Sylvan plus the data archive, on any machine.

It does not cover the figures of Section 6 that compare with SliQSim and
Quasimodo, which come from the paper's original artifact.

In the commands, `$QSY` is a clone of this repository and `$DATA` the unpacked
data archive. `$PAPER`, a checkout of the paper, is only needed to compare
the remade figures with the paper's copies.

## 0. What you need

**The code.** Check out the tag `exactdd-clean-sweep-2026-10` of this
repository: the commit "limdd: the code of the exact-DD clean sweep of 6-7
October (bin-lazy6)" plus the scripts.

- The paper's runs were made with exactly this code.
- Later commits fix float LIMDD canonicity, which changes the float LIMDD
  results and all runtimes.

```
git clone -b limdd https://github.com/alaarman/q-sylvan.git $QSY
cd $QSY && git checkout exactdd-clean-sweep-2026-10
```

**The data.** `exactdd-data-2026-10.tar.xz`, kept with the paper's sources in
`reproducibility/`, or from the authors. Its SHA-256 is in
`reproducibility/SHA256SUMS` there. Unpacked, it is `$DATA`:

```
tar xJf exactdd-data-2026-10.tar.xz     # -> exactdd-data/
DATA=$PWD/exactdd-data
```

What is in the archive:

- every circuit;
- every run of the clean sweep (`.status`, `.json`, `.out`), with its audit and
  logs;
- the older runs that `regen_amendment.sh` reads for the 2026-10-01 revision;
- the source diff of the build.

`$DATA/DATA_README.md` lists the layout.

**Tools.**

- a C/C++ compiler;
- CMake 3.23 or later;
- GMP, MPFR and MPC (Debian: `libgmp-dev libmpfr-dev libmpc-dev`; macOS:
  `brew install gmp mpfr libmpc argp-standalone`);
- GNU `timeout` (coreutils), `bc` and `python3`.

The plots and the circuit generators need Python 3.11 with pinned packages
(section 8).

## 1. The pipeline

```
circuits (.qasm)               $DATA/ldd/qasm, $DATA/amend/qasm, $DATA/pilot/qasm/<family>/
  -> runner                    run_qasm_on_qmdd, built twice (section 2)
  -> one run                   run_arm.sh ARM CIRCUIT OUTDIR [TIMEOUT]       (section 4)
  -> all runs                  sweep.py DATA OUT                             (section 5)
  -> <circuit>_<arm>.{status,json,out}   OUT/amend/out_clean, OUT/pilot/out_clean
  -> state-vector references   reference_probs.py                            (section 6)
  -> figures                   make_paper_figures.sh ROOT OUTDIR [PAPER]     (section 7)
```

The scripts are in `$QSY/benchmark/exactdd/amendment/`, and the circuit
generators one level up.

## 2. Build

Two Release builds of the runner are needed. They differ in the number of
64-bit words in a Pauli string:

- 1 word: every EVDD run, and the LIMDD runs up to 64 qubits;
- 2 words: the LIMDD runs above 64 qubits.

Always pass `-S` and `-B`. On macOS, argp is not in libc, and the CMake module
that finds it is in this repository's history:

```
cd $QSY
git show 1570700:cmake/FindArgp.cmake > cmake/FindArgp.cmake    # macOS only
cmake -S . -B build-w1 -DCMAKE_BUILD_TYPE=Release -DLIMDD_PAULI_WORDS=1
cmake --build build-w1 --target run_qasm_on_qmdd
cmake -S . -B build-w2 -DCMAKE_BUILD_TYPE=Release -DLIMDD_PAULI_WORDS=2
cmake --build build-w2 --target run_qasm_on_qmdd
export QSY_W1=$QSY/build-w1/qasm/run_qasm_on_qmdd
export QSY_W2=$QSY/build-w2/qasm/run_qasm_on_qmdd
```

The paper's binaries were built with Apple clang 17.0.0 and CMake 4.3.2 on
macOS 26 (arm64). Their SHA-256 sums are in `$DATA/source/SHA256SUMS`. The
binaries themselves are not distributed: a rebuild is not byte-identical, and
need not be.

## 3. Circuits

The archive holds every circuit. To regenerate them and compare:

| set | count | in `$DATA` | made by |
|---|---|---|---|
| random Clifford+T | 49 | `ldd/qasm/` | the paper's artifact. Identical by blob hash to `benchmark/qasm_clifford_T/` of git ref `quist/qisq2_tests` (`git fetch https://github.com/Quist-A/q-sylvan qisq2_tests`) |
| Grover, W state (two constructions), hidden shift, adder, Ising | 42 | `amend/qasm/` | `../gen_grover.py`, `../gen_wstate_ancilla.py`, `../gen_wstate_clifford.py`, `../gen_realistic.py`; commands in `README.md`, section "Circuits" |
| mirror, IQP, Grover with many iterations, quantum walks, synthesised QFT/QAOA/TFIM | 91 | `pilot/qasm/<family>/` | `../gen_new_families.sh` |

```
PYQ=$PWD/venv-qiskit/bin/python $QSY/benchmark/exactdd/gen_new_families.sh /tmp/newfam $DATA/ldd/qasm
diff -r -x '*.json' /tmp/newfam $DATA/pilot/qasm        # no output: all 91 identical
```

Each synthesised circuit has a `.json` sidecar with its generation parameters
(seed, eps, SK depth, T-count, rotation error). Only its `generation_time_s`
changes when regenerated.

The figures draw the random, Grover and W-state circuits and the new families.
Hidden shift, adder and Ising are run but not drawn, because Section 6 has no
such family.

## 4. One run

`run_arm.sh ARM CIRCUIT OUTDIR [TIMEOUT_S]` runs one arm on one circuit and
writes three files:

- `<circuit>_<arm>.json`, the runner's statistics;
- `<circuit>_<arm>.out`, the tail of the output, the command, the exit code and
  the wall time;
- `<circuit>_<arm>.status`, one line:
  `arm circuit status time final_nodes peak_nodes norm p0 t_count final_width`.

The status is `OK`, `TIMEOUT`, `TAB_FULL`, `LIM_FULL`, `STAB_FULL`, `ERR(rc)`
or `NOJSON`. A run whose `.status` exists is not repeated.

```
cd $QSY/benchmark/exactdd/amendment
./run_arm.sh qisq2_low $DATA/ldd/qasm/clifford_T_circuit_20_700.qasm /tmp/myruns 1800
#   qisq2_low clifford_T_circuit_20_700 OK 6.53 83275 124300 1.0 0.5 238 0
cat $DATA/amend/out_clean/clifford_T_circuit_20_700_qisq2_low.status
#   the same nodes, norm and p0; the time is the M1's
```

Every arm adds `-c -m -w 1 --node-tab-size=25 --wgt-tab-size=24`, with a
default timeout of 1800 s. These are the paper's six arms:

| arm | diagram | weights | flags |
|---|---|---|---|
| `qisq2_low` | EVDD | exact (algebraic) | `--dd=qmdd -e qisq2 -s low` |
| `float_low` | EVDD | float, low | `--dd=qmdd -e float -s low --merging=abs -t 1e-14` |
| `float_l2` | EVDD | float, L2 | `--dd=qmdd -e float -s l2 --merging=abs -t 1e-14` |
| `limdd_qisq2` | LIMDD | exact | `--dd=limdd -e qisq2 -s low` |
| `limdd_float` | LIMDD | float, low | `--dd=limdd -e float -s low --merging=abs -t 1e-14` |
| `limdd_float_l2` | LIMDD | float, L2 | `--dd=limdd -e float -s l2 --merging=abs -t 1e-14` |

`run_arm.sh` knows more arms (`limdd_lazy`, `limdd_qisq2_lazy`, and others for
later work); its header lists them. Leave `QSY_AMP_SAMPLE` and `QSY_AMP_AT`
unset: the runner of this tag does not have the options they switch on.

In the JSON:

- `first_qubit_measurement_prob` is the probability that qubit 0 (`q[0]`)
  measures 0;
- `norm` is the norm of the final state.

The T-count figure's ten per-gate traces come from `run_trace.sh`:

```
./run_trace.sh 20_700 $DATA/ldd/qasm /tmp/mytraces
```

## 5. All runs: the sweep

`sweep.py` repeats the paper's sweep: every arm on every circuit, plus the ten
traces. It is the driver that made the paper's runs, kept as
`$DATA/sweep/clean_sweep.py`, with its paths made arguments:

```
cd $QSY/benchmark/exactdd/amendment
./sweep.py $DATA /path/to/OUT --dry      # the plan
./sweep.py $DATA /path/to/OUT > sweep.log
```

It writes `OUT` in the data layout, so `OUT` can be given to the figure
scripts:

- `OUT/amend/out_clean`, `OUT/amend/trace_clean` and `OUT/pilot/out_clean`;
- the lazy canonical form, which is not in the figures, in `OUT/*/out_clean_lazy`;
- links to `$DATA`'s circuits.

**Phases**, by priority:

| phase | what |
|---|---|
| A | the 72 circuits of Figures 6 and 7 among the old benchmarks (49 random, 12 Grover, 11 W state) |
| T | the ten traces |
| B | hidden shift, adder and Ising |
| C | the 91 new circuits |
| D, E | the lazy canonical form |

`--phases`, `--arms` and `--only REGEX` select a part. The new circuits come
in ladders of growing size. After an arm times out or fills a table on one
instance, the larger instances of that ladder are recorded as
`SKIPPED_PREDICTED` and not drawn.

The paper's sweep took about a day on the machine below. Most of that is
timeouts: 30 minutes each, with at most six runs at a time.

**Timing.** Runtimes are only meaningful when nothing else slows the machine.
On 6 October, a first attempt with other load on the machine gave node counts
that matched earlier runs but runtimes up to 12 times longer. So the sweep:

- caps the number of runs at `maxjobs`;
- starts a run only when the simulators' memory, plus the new run's estimate,
  stays within a budget;
- stops and requeues the youngest run if the memory still exceeds a hard
  limit;
- holds new runs while the machine swaps;
- reruns once any run of more than 20 s that spent more than `max_sys` of its
  CPU time in the kernel, or less than `min_cpu` of its wall time on the CPU.

The limits are read from `OUT/sweep/sweep.conf` (JSON) while the sweep runs,
so they can be changed on the fly. The defaults are in the script.

The paper's sweep used:

- budget 9 to 18 GB and hard limit 11 to 15 GB, depending on what else was open;
- `maxjobs` 6, which keeps the runs on the performance cores;
- `max_sys` 0.06 and `min_cpu` 0.85.

Set `maxjobs` to at most the number of performance cores, and the budget to
what the machine can spare.

Memory is measured on macOS (`phys_footprint`, compressed memory included) and
on Linux (resident plus swapped). Elsewhere only `maxjobs` limits the sweep.

**Audit.** `OUT/sweep/audit.csv` has one line per run:

- start and end times;
- wall time and simulation time;
- CPU time, user and system;
- peak footprint and maximum RSS;
- the total memory in use while the run ran, and swap-outs meanwhile.

A run whose CPU time is well below its wall time was slowed by something else.
The paper's sweep has the same file as `$DATA/sweep/clean_audit.csv`.

**The machine of the paper's runs.**

- Apple M1 Pro (8 performance and 2 efficiency cores), 32 GB;
- macOS 26 (Darwin 25.2), on power, kept awake with `caffeinate`;
- each run single-threaded (`-w 1`).

**What a repeat must give.**

- **Exact arms** (`qisq2_low`, `limdd_qisq2` and the lazy `limdd_qisq2_lazy`):
  the same status, node counts, norm and probability on any machine. Only
  timeouts near the limit can change, with the machine's speed.
- **Float arms**: these depend on how the platform rounds:
  - Q-Sylvan's float weight table keys its hash on `long double`, which is
    64-bit on macOS arm64 and wider on x86-64 and on Linux arm64;
  - the C compiler decides whether multiply-adds are fused.

  On another platform, a float run can therefore end with other node counts,
  or another verdict. Section 9 has a check that counts how many do.
- **Runtimes**: those of the machine.

This was checked by repeating 110 runs from a fresh clone: 14 circuits from
every family, every arm, and two traces, one of them above 64 qubits.
`checks/compare_sweeps.py` compared them with the stored runs:

| platform | exact arms | float, low | float, L2 |
|---|---|---|---|
| macOS 26, arm64, Apple clang 17 (the paper's) | identical | identical | identical |
| Debian 12, aarch64 (container), gcc 12.2, glibc 2.36 | identical | identical | other node counts on 5 of 13 EVDD and 2 of 12 LIMDD runs, with the same probability and norm |

On both platforms the 36 figures came out byte-identical, made from the
archive with the pinned Python environments (section 8).

## 6. Reference probabilities

On 9 synthesised and quantum-walk circuits, no exact run (EVDD or LIMDD) gives
a probability. They time out, or in one case the exact coefficients overflow a
double. For these 9 circuits only, the float runs are judged against a
double-precision state-vector simulation.

- Script: `reference_probs.py`. Qiskit parses the circuit, and a sparse state
  vector keeps only the nonzero amplitudes.
- Values: `$DATA/pilot/reference_probs.csv`.
- Validation: `$DATA/pilot/validation_probs.csv`. It covers ten circuits from
  five families on which exact runs exist, and agrees with them to every digit
  printed.

```
venv-qiskit/bin/python reference_probs.py /tmp/ref.csv $DATA/pilot/qasm/synth_qft/synth_qft_n12_e3.qasm ...
```

## 7. Figures

`make_paper_figures.sh ROOT OUTDIR [PAPERDIR]` remakes, from the stored runs
alone, every figure of this pipeline that the paper typesets, in the paper's
layout:

```
cd $QSY/benchmark/exactdd/amendment
PY=$PWD/venv-plot/bin/python PY_TCOUNT=$PWD/venv-tcount/bin/python \
  ./make_paper_figures.sh $DATA /tmp/figs $PAPER
#   36 figures in /tmp/figs/figures
#   compared with $PAPER/figures: 36 same, 0 differ, 0 not in the paper
```

`ROOT` is either `$DATA` or a sweep's `OUT`. `/tmp/figs/figures/` then holds
the paper's figure files:

| paper | file in `figures/` |
|---|---|
| Fig. 6, EVDD panels and legend | `final_plots_mm-(in)correct2/pooled_wgt_{type_runtime,final_nodecount,nodecount}_float_{low,l2}_vs_algebraic_nw.pdf`, `pooled_legend_float_{low,l2}_vs_algebraic_nw.pdf` |
| Fig. 6, LIMDD panels | `limdd2/` (L2 floats) and `limdd2_low/` (low floats): `limdd_float_vs_algebraic_{simulation_time,final_nodes,max_nodes}.pdf` |
| Fig. 7 | `limdd2/`, `limdd2_low/`: `limdd_vs_evdd_*.pdf`, `panel_legend.pdf` |
| T-count figure | `limdd/width_vs_tcount.pdf`, `limdd/bits_vs_tcount.pdf` |

Where the files come from:

- `\iflowfloats` in the paper's `macros.tex` chooses between the low and L2
  versions.
- With `PAPERDIR`, every file is compared with the paper's copy by
  `checks/pdf_same.py`, which ignores the creation date. With `PAPERDIR
  --install`, the files are copied into it instead.
- The script runs `regen_amendment.sh` with the families of the text
  (`AMEND_FAMILIES=random,grover,wstate,mirror,iqp,groveriter,qwalk,synth`).
- The per-panel counts are in `/tmp/figs/paper_*.out`: points, wrong runs, runs
  that did not finish, and how each float run was judged.

`regen_amendment.sh ROOT OUTDIR` makes everything, including the 2026-10-01
revision's figures and numbers when `ROOT` has that revision's runs, as
`$DATA` does. `README.md` describes those steps.

The correctness rule and the colours are defined in the docstrings of
`make_figures.py` and `make_section6_figures.py`. A float run is wrong when:

- its top-qubit probability is more than 5% off the exact one; or
- its norm is more than 1e-3 off 1.

The colours:

- Figure 6: red is a wrong measurement, violet a wrong norm only.
- Figure 7: a black fill marks a wrong float LIMDD, an orange ring a wrong
  float EVDD.

## 8. Python environments

Three environments, each CPython 3.11 with the pinned versions in this
directory:

| file | for | made the paper's files |
|---|---|---|
| `requirements-plot.txt` | Figures 6 and 7, `checks/` | matplotlib 3.7.1, numpy 1.26.4 |
| `requirements-plot-tcount.txt` | the T-count figure (and the 2026-10-01 revision's PDFs) | matplotlib 3.11.2, numpy 2.4.6 |
| `requirements-qiskit.txt` | the generators, `reference_probs.py` | qiskit 2.5.2, pygridsynth 2.0.0 |

```
for e in plot:requirements-plot.txt tcount:requirements-plot-tcount.txt qiskit:requirements-qiskit.txt; do
  python3.11 -m venv venv-${e%%:*} && venv-${e%%:*}/bin/pip install -r ${e#*:}
done
```

With these versions, the figures are byte-identical to the paper's apart from
their creation date. Other versions plot the same data, but the PDF bytes
differ. `run_arm.sh` itself needs only a `python3` with the standard library.

## 9. Checks you can repeat

Run these with the plotting Python. `ROOT` is `$DATA` or a sweep's `OUT`.

- **The exact arms agree.** `make_figures.py` prints `exact arms agree on the
  top-qubit probability on N/N circuits` into `/tmp/figs/paper_limdd2.out`.
  The paper's runs give 145 of 145.
- **The counts in the text.** `checks/family_stats.py ROOT`, also written by
  `make_paper_figures.sh` to `paper_family_stats.out`. Per family and per pair
  of arms (low and L2), it gives:
  - the wrong runs;
  - the faster runs;
  - the runs that did not finish;
  - the circuits where only the exact or only the float run finishes.
- **A repeat against the paper's runs.**
  `checks/compare_sweeps.py $DATA OUT [FAMILIES]`. For each arm it gives:
  - the runs that finish in each sweep;
  - whether the exact arms keep their node counts (they must);
  - the time ratio between the two sweeps.
- **The clean sweep against the runs of 5 October.**
  `checks/compare_sweeps.py --paper5 $DATA`. These are the runs behind the
  figures before the clean sweep. The exact runs have identical node counts on
  every circuit that both finish.
- **Timing quality.** In `sweep/clean_audit.csv` (or `OUT/sweep/audit.csv`),
  compare `cpu_s` with `wall_s`, and `stime_s` with `cpu_s`. Clean runs spend
  1.5 to 3.3% of their CPU time in the kernel. Three of the paper's runs spend
  7.1 to 7.3% (`float_l2` on 40_900, 30_800 and 80_950); they were accepted
  under an early, looser limit and are at most about 5% slow.
- **One figure.** `checks/pdf_same.py A.pdf B.pdf`.
