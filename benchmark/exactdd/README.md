# Random Clifford+T benchmarks for the exact-DD paper

Three sets live here. Exactly one of them is the paper's artifact; the other two
are ours. Keep the distinction — reviewers will ask.

| directory | what it is | provenance |
|---|---|---|
| `orig/qasm/` | **the paper's actual artifact**, 49 files | extracted verbatim from git ref `quist/qisq2_tests`, path `benchmark/qasm_clifford_T/` (ref head `4e7af4c`) |
| `gen_extension/` | **our extension**: the grid points the paper claims but never shipped | generated here by `gen_random_cliffT.py` |
| `gen_truedepth/` | **our reconstruction** of what the paper's prose literally says | generated here by `gen_random_cliffT.py` |

The generator and this README live in the q-sylvan tree; the generated circuit
sets do not — regenerate them with the commands below (the seed makes them
byte-identical).

## What the paper says vs. what the artifact is

`experiments.tex`: *"random circuits ... gates randomly sampled from the
Clifford+T gates with 20 to 100 qubits and circuit depth between 700 and 1000"*
and *"Each circuit is followed by a single-qubit measurement on the top qubit."*

Three mismatches with the surviving artifact:

1. **"depth" is a gate count.** The artifact's `clifford_T_circuit_<n>_<g>.qasm`
   second number is the number of gate lines in the file, confirmed by
   `benchmark/run_benchmarks.sh` on the same ref, which loops
   `num_gates=j*50, j in 14..20` and `num_qubits=i*10, i in 2..8`. The real ASAP
   depth of `clifford_T_circuit_20_700.qasm` is **102**, not 700; of
   `clifford_T_circuit_80_700.qasm` it is **27**.
2. **"20 to 100 qubits" is 20 to 80.** The grid is `{20,30,...,80}`; there are no
   90- or 100-qubit files anywhere on the ref. `gen_extension/` fills that in.
3. **No measurement is in the files.** No `creg`, no `measure` in any of the 49.
   The top-qubit measurement is the runner's `-m`
   (`first_qubit_measurement_prob`), and the paper's own `run_benchmarks.sh`
   does not even pass `-m`.

## The generator

`gen_random_cliffT.py` — dependency-free (stdlib `random` only), no qiskit.

It reproduces the **artifact's** distribution, which is *not* what
`benchmark/generate_qasm_random_clif_T_gate.py` in the q-sylvan tree emits
today. The in-tree script differs on four counts:

- it emits `sx`/`sxdg`, which **the LIMDD path refuses** — so as it stands the
  in-tree script cannot produce a runnable LIMDD benchmark at all;
- it emits `cy`, which appears nowhere in the artifact;
- it has `swap` commented out, though the artifact is ~11% swap;
- its `globals_pipeline.PERCENTAGE_T_GATES` is 85, against the artifact's ~34%.

### Sampling model (measured off all 49 artifact files, 41 650 gates)

Three equiprobable classes; the gate is uniform inside its class.

| class | gates | artifact share |
|---|---|---|
| T | `t` | 33.81% |
| 1-qubit Clifford | `h s sdg x y z` (**no** `sx`/`sxdg`) | 33.25% total, ~5.5% each |
| 2-qubit Clifford | `cx cz swap` (**no** `cy`) | 32.94% total, ~11% each |

So `--t-frac 1/3` with the rest split evenly. Default in the script is exactly
`1/3`; the artifact's 0.3381 is +2.1σ off 1/3 on 41 650 draws, i.e. consistent.

**Two-qubit operand order is random, not ascending.** The artifact draws an
ordered pair `(c,t)`, `c != t`, uniformly: 6 803 of its 13 721 two-qubit gates
have `control < target` and 6 918 have `control > target`. The in-tree
generator instead rejects until `target > control`. The default here is
`--ctrl-order random` to match the artifact; `--ctrl-order ascending`
reproduces the in-tree convention if a control-before-target variable order is
ever wanted.

File shape matches the artifact byte-for-byte in style: `OPENQASM 2.0;` /
`include "qelib1.inc";` / `qreg q[n];`, then one gate per line, **no `creg`, no
`measure`, and no trailing newline** (all 49 artifact files end on the last `;`).

### Depth mode

`--depth D` grows the circuit and stops the moment its ASAP depth first reaches
`D`, so the depth is **exact**, not estimated. Depth is the standard DAG
longest path (`level[q]`, a gate lands at `1 + max` over its operands); the
script re-measures every file it writes from disk and asserts agreement.

The earlier `depth × n / 1.2` estimate for how many gates a target depth needs
is wrong by ~2.7×. Measured gates-per-depth-layer is **≈ 0.30 n**, not
`n/1.2 ≈ 0.83 n`: random greedy packing reaches only ~40% of the perfect-packing
bound (average 4/3 qubits per gate ⇒ 0.75 n gates per layer at best). The
`--depth` mode does not rely on any such estimate.

### Seeds

`seed = 1234 + 1000*qubits + size`, where `size` is the gate count (`--gates`)
or the target depth (`--depth`). Fully reproducible; re-running either command
below regenerates identical files.

## Reproducing the two sets

```
python3 gen_random_cliffT.py --out gen_extension \
        --qubits 90 100 --gates 700 750 800 850 900 950 1000

python3 gen_random_cliffT.py --out gen_truedepth --name clifford_T_depth \
        --qubits 20 40 60 80 --depth 700 1000
```

## Verification against the artifact

Same grid as the artifact (7 widths × 7 gate counts, 41 650 gates), generated
and compared aggregate:

| gate | artifact | ours |
|---|---|---|
| t | 33.81% | 33.25% |
| cz | 11.18% | 11.23% |
| swap | 10.91% | 10.91% |
| cx | 10.85% | 11.22% |
| sdg | 5.67% | 5.53% |
| h | 5.62% | 5.62% |
| y | 5.57% | 5.56% |
| s | 5.50% | 5.61% |
| x | 5.46% | 5.56% |
| z | 5.43% | 5.51% |
| `ctrl<targ` / `ctrl>targ` | 6803 / 6918 | 6838 / 7059 |

Gate sets identical, every share within 0.6 points.

## `gen_extension/` — 14 files, qubits {90,100} × gates {700..1000 step 50}

Named `clifford_T_circuit_<n>_<gates>.qasm`, i.e. the artifact's convention,
because they are the artifact's missing grid points. Measured aggregate over
11 900 gates: t 33.26%, cx 11.19%, cz 10.87%, swap 10.69%, singles 5.1–6.2%.
Their real ASAP depths are 24–40.

These are **not** trivial despite the small gate count: at 100 qubits / 1000
gates the exact QMDD reaches 296 697 final nodes. See the sanity table.

## `gen_truedepth/` — 8 files, qubits {20,40,60,80} × depth {700,1000}

Named `clifford_T_depth_<n>_<depth>.qasm` so they can never be confused with
the artifact's gate-count files. Depth is exact and re-measured from disk:

| file | qubits | gates | measured depth | gates/layer |
|---|---|---|---|---|
| `clifford_T_depth_20_700.qasm` | 20 | 4 911 | 700 | 7.02 |
| `clifford_T_depth_20_1000.qasm` | 20 | 6 967 | 1000 | 6.97 |
| `clifford_T_depth_40_700.qasm` | 40 | 9 076 | 700 | 12.97 |
| `clifford_T_depth_40_1000.qasm` | 40 | 13 637 | 1000 | 13.64 |
| `clifford_T_depth_60_700.qasm` | 60 | 13 692 | 700 | 19.56 |
| `clifford_T_depth_60_1000.qasm` | 60 | 19 758 | 1000 | 19.76 |
| `clifford_T_depth_80_700.qasm` | 80 | 17 383 | 700 | 24.83 |
| `clifford_T_depth_80_1000.qasm` | 80 | 25 349 | 1000 | 25.35 |

Aggregate over 110 773 gates: t 33.39%, cx 11.10%, cz 11.05%, swap 11.25%,
singles 5.45–5.62%.

These are 7–25× more gates than the artifact's files at the same width. Taking
the paper's prose literally is a far harder benchmark than the paper ran.

## `gen_ladder/` — 20-qubit depth ladder (diagnostic, not a benchmark set)

Five files, depth {50,100,150,200,300} at 20 qubits, same generator and seed
rule. Generated to find where each backend gives up; see the caveat below.

## Sanity runs

Flags: `-c -m -w 1 --node-tab-size=25 --wgt-tab-size=23`, `timeout 300`,
one worker. `qisq2` arms are `-e qisq2 -s low`; `float` arms are
`-e float -s low --merging=abs`. Binary is `build-clean` except LIMDD above 64
qubits, which uses `build-w2rel` (`LIMDD_PAULI_WORDS=2`).

### `gen_extension/` — everything completes except one LIMDD run

| circuit | arm | status | time | final nodes | peak | norm | p(top) |
|---|---|---|---|---|---|---|---|
| 90_700 | QMDD qisq2 | OK | 0 s | 4 601 | 4 601 | 1.0 | 0.5 |
| 90_700 | LIMDD qisq2 (w2) | OK | 53 s | **119** | 148 | 1.0 | 0.5 |
| 90_1000 | QMDD qisq2 | OK | 1 s | 14 890 | 19 946 | 1.0 | 0.0 |
| 90_1000 | LIMDD qisq2 (w2) | OK | 73 s | **77** | 90 | 1.0 | 0.0 |
| 100_700 | QMDD qisq2 | OK | 0 s | 4 096 | 4 096 | 1.0 | 0.853553 |
| 100_700 | LIMDD qisq2 (w2) | OK | 16 s | **124** | 124 | 1.0 | 0.5 |
| 100_1000 | QMDD qisq2 | OK | 21 s | 296 697 | 315 129 | 1.0 | 1.0 |
| 100_1000 | LIMDD qisq2 (w2) | TIMEOUT | 300 s | | | | |

These are **not** trivial, contrary to what one might guess from 1000 gates on
100 qubits: the exact QMDD hits 296 697 nodes on `100_1000`. LIMDD is 30–190×
smaller where it finishes, and its one timeout is canonicalisation cost, not a
blow-up. Two caveats worth stating in the paper if these are used: the top-qubit
outcome is deterministic (p = 0.0 or 1.0) on two of the four, and at these widths
each qubit is touched only ~10 times, so the circuits are shallow by construction.

### `gen_truedepth/` — every arm exceeds a 300 s budget

All 12 exact runs (`qisq2`, both QMDD and LIMDD, widths 20/40/60/80, depths
700 and 1000) hit `TIMEOUT` at 300 s. The float pass:

| circuit | arm | status | time | final | peak | norm |
|---|---|---|---|---|---|---|
| depth_20_700 | QMDD float | OK | 0 s | 9 | 6 261 | **0.0** |
| depth_20_700 | LIMDD float | OK | 45 s | 0 | 7 733 | **0.0** |
| depth_40_700 | QMDD float | **TAB_FULL** | 84 s | | | |
| depth_40_700 | LIMDD float | TIMEOUT | 300 s | | | |
| depth_60_700 | QMDD float | TIMEOUT | 300 s | | | |
| depth_60_700 | LIMDD float | TIMEOUT | 300 s | | | |
| depth_80_700 | QMDD float | TIMEOUT | 300 s | | | |
| depth_80_700 | LIMDD float | TIMEOUT | 300 s | | | |
| depth_80_1000 | QMDD float | TIMEOUT | 300 s | | | |
| depth_80_1000 | LIMDD float | TIMEOUT | 300 s | | | |

**`norm = 0.0` is float giving up, not a valid answer.** A unitary circuit
cannot produce a zero state; absolute-tolerance merging has annihilated it, and
the reported `p(top)` is then a division by ~0 (one run printed 1e10).

The 20-qubit ladder in `gen_ladder/` locates the cliff:

| depth | gates | float: final / norm | qisq2: final / norm |
|---|---|---|---|
| 50 | 393 | 548 / 1.0 | 548 / 1.0 |
| 100 | 714 | 854 541 / **1.0218** | 25 589 / 1.0 |
| 150 | 1 083 | 1 / **0.0** | TIMEOUT (120 s) |
| 200 | 1 474 | 7 / **0.0** | TIMEOUT (120 s) |
| 300 | 2 296 | TIMEOUT | TIMEOUT |

At 393 gates float and exact agree node-for-node. By 714 gates — the artifact's
own size — this instance's float run is 33× larger than exact *and* 2.2% off in
norm; by 1 083 gates the float state is gone. Instance variance is large: the
artifact's `clifford_T_circuit_20_700` gave float norm 0.999999, this ladder's
714-gate instance gives 1.0218.

**Conclusion on `gen_truedepth/`:** it is a correct realisation of the paper's
prose and it is far out of reach — 300 s is not enough for any arm, and float is
not merely slow there but wrong. If the true-depth reading is to be used at all,
it needs either much longer budgets or smaller depths (the ladder suggests
exact QMDD at 20 qubits tops out somewhere between depth 100 and 150 within a
2-minute budget). The paper's gate-count reading is the tractable one, which is
presumably why the artifact is what it is.

## IQP circuits and the BQD arm

`gen_iqp.py` writes the IQP fragment of a circuit: a Hadamard on every qubit,
then diagonal gates drawn from the third level of the Clifford hierarchy on at
most three qubits (`t`, `z`, `s`, `sdg`, `cz`, `cs`, `csdg`, `ccz`), and no
closing layer. Every intermediate state is a phase state, which is the class
the binary quotient diagram (`--dd=bqd`, `src/qsylvan_bqd.h`) has an algorithm
for, and the file stops where all three arms can go: the BQD has no Hadamard.

The `cs` and `ccz` gates are what makes the set worth running. With `t` and
`cz` alone the high cofactor at every level is a Pauli times the low one, so
the LIMDD keeps width 1 and n nodes however many T gates there are (measured:
19 nodes at 20 qubits, 29 at 30). The parser and all three arms take the three
gates as one gate each; the LIMDD arm goes through `limdd_cgate` with a
control mask, since they are symmetric and the highest qubit can be the
target. `cp` with an angle is not used because it is not exact on `qisq2`.

The T-count reported by the runner counts `cs` as 3 and `ccz` as 7, their
ancilla-free Clifford+T decompositions, on every arm.

```
python3 gen_iqp.py --out gen_iqp --qubits 20 30 40 50 60 --gates 700 1000
run_qasm_on_qmdd gen_iqp/iqp_20_700.qasm --dd=bqd   -e qisq2 -s low -c -m
run_qasm_on_qmdd gen_iqp/iqp_20_700.qasm --dd=limdd -e qisq2 -s low -c -m
run_qasm_on_qmdd gen_iqp/iqp_20_700.qasm --dd=qmdd  -e qisq2 -s low -c -m
```

The BQD arm reports `norm` and `first_qubit_measurement_prob` as -1 above 24
qubits: it has no summation algorithm, and below that they come from decoding
the whole state, which is how `qasm/test/test_sim_qasm_bqd.py` checks the arm
against the EVDD arm amplitude for amplitude. Its node counts follow the
paper's convention (no level skipping), so they are not the LIMDD's.
