# Hybrid edge-weight merging for Q-Sylvan

Implements the merging rule of *"Numerical Errors in Quantitative System Analysis
With Edge-Weighted Decision Diagrams"*: a **relative** tolerance for non-zero
weights plus a **separate absolute** threshold for collapsing a weight to zero,
replacing the single absolute tolerance used today.

    MakeEdge(w):  if |w| <= zero_tol      -> 0
                  else lookup w' with  |w - w'| <= rel_tol * max(|w|,|w'|)

## Changes

| file | change |
|---|---|
| `src/edge_weight_storage/cmap.{c,h}` | tolerance *mode* (`CMAP_TOL_ABS` default, `CMAP_TOL_HYBRID`); relative equality predicate; log-polar hash key so relatively-close weights share a bucket, mirroring the existing uniform grid |
| `src/sylvan_edge_weights.{c,h}` | `sylvan_edge_weights_set_hybrid_tolerance(rel, zero)`; `..._tolerance_from_env()` reading `QSYLVAN_REL_TOL` / `QSYLVAN_ZERO_TOL` |
| `src/sylvan_evbdd.c` | one call to the env hook in `sylvan_init_evbdd` |
| `examples/alg_run.c` | `--rel-tol=<x> --zero-tol=<x>` |
| `examples/hybrid_demo.c` | new: H-layer benchmark reaching the small-scale-factor regime |

Default behaviour is unchanged: without `--rel-tol` / `QSYLVAN_REL_TOL` the
historical absolute rule runs, bit-for-bit.

## Result 1 — the scale-factor failure (`hybrid-results-hadamard.csv`)

`H` on every qubit of `|0..0>`, twice; the exact result is `|0..0>`.  The state's
root weight is `2^{-n/2}`, an exponentially small *scale factor*.

| n | absolute 1e-14: amp(0..0) | hybrid rel=1e-14, zero=0 | nodes (abs / hyb) |
|---|---|---|---|
| 92 | 1.0000000000000 | 1.0000000000000 | 93 / 93 |
| 94 | **2.0** | 0.99999999999999 | 95 / 95 |
| 100 | **16.0** | 0.99999999999998 | 101 / 101 |
| 120 | **16384.0** | 0.99999999999998 | 121 / 121 |
| 160 | **1.7e10** | 0.99999999999997 | 161 / 161 |

The absolute rule diverges as `sqrt(2)^(n-92)` once `2^{-n/2}` falls below the
tolerance — the Q-Sylvan analogue of mqt-core issue #575.  The hybrid rule stays
at the rounding floor (~`n*eps`) with **identical node counts**.  All four
normalisation strategies fail identically under the absolute rule and are all
fixed by the hybrid rule: this failure is a property of the merging rule, not of
the normalisation.

## Result 2 — zero-collapse is a separate knob

Same instance, n=100, `NORM_MAX`, `rel=1e-14`:

| zero_tol | 0 | 1e-30 | 1e-20 | 1e-16 | 1e-14 |
|---|---|---|---|---|---|
| error | 1.6e-14 | 1.6e-14 | 1.6e-14 | 1.6e-14 | **1.0 (state = 0)** |

`zero_tol` must sit below the smallest legitimate amplitude (`2^{-50} = 8.9e-16`
here).  It is an approximation budget, not a merging tolerance, and conflating
the two is what the single absolute tolerance does.

## Result 3 — Grover (`hybrid-results-grover.csv`), norm / nodes

Under the contractive `NORM_MAX`, hybrid matches absolute exactly (norm 1,
12-14 nodes).  Under the non-contractive `NORM_LOW` / `NORM_MIN` the absolute
rule returns 0.5, 0.0, or 2.5e20; only the pure relative rule (`zero_tol=0`)
recovers norm 1 -- at the cost of the compression that zero-merging provides,
since Grover's diagram depends on it.  Under `NORM_L2` hybrid is often the more
compact (n=8: 67 vs 130 nodes; n=12: 846 vs 1431).

So the two conditions are independent, as the analysis predicts: a contractive
normalisation *and* a relative merging rule are both needed, and neither
substitutes for the other.

## Result 4 -- performance

**Per-lookup cost (`hybrid-results-perlookup.txt`).** The H-layer benchmark gives
*identical* node counts under both rules, so a runtime difference there is the
cost of the rule itself and not of a different amount of DD work.  Comparing the
marginal cost per layer (init is the intercept, so slope isolates per-gate cost):

| n | absolute | hybrid | ratio |
|---|---|---|---|
| 40 | 0.1750 s/layer | 0.1817 s/layer | 1.038 |
| 80 | 0.3650 | 0.3633 | 0.995 |
| 120 | 0.5500 | 0.5567 | 1.012 |
| 160 | 0.7500 | 0.7517 | 1.002 |

Ratios straddle 1.0.  The log-polar hash (`logl`, `atan2l`) and the `hypotl` calls
in the predicate are free in practice: the table lookup is bound by memory access,
not arithmetic.  No squared-magnitude micro-optimisation is warranted.

**End-to-end on MQT Bench (`hybrid-results-mqtbench.csv`).** 198 circuits of the
`MQTBench_2024-10-10-All-Qiskit-IBM-2-10` archive shipped in `benchmark/`, all
families with <= 12 qubits (28 families: ae, dj, ghz, graphstate, grover, qaoa,
qft, qpeexact, qwalk, vqe, wstate, portfolio*, ...), x 2 normalisation strategies
x 3 merging rules, single-threaded, 60 s timeout.  Paired ratios against the
absolute rule:

| strat | rule | time median | time p90 | time total | nodes median | nodes p90 | nodes max |
|---|---|---|---|---|---|---|---|
| max | rel, zero=0 | 1.006 | 1.050 | 1.025 | 1.000 | 1.429 | 9.11 |
| max | rel, zero=1e-16 | 1.009 | 1.046 | 1.008 | 1.000 | 1.091 | 2.12 |
| l2 | rel, zero=0 | 1.011 | 1.050 | 1.035 | 1.000 | 1.167 | 8.15 |
| l2 | rel, zero=1e-16 | 1.011 | 1.048 | 1.014 | 1.000 | 1.100 | 1.82 |

So roughly **1 % slower in the median, 1-3 % on total runtime**, and node counts
are unchanged for the median circuit.  With `zero_tol = 0` a minority of circuits
(31/196 more than 5 % larger; qwalk-noancilla, qpeexact) blow up to as much as
9x, because nothing ever merges with zero; setting `zero_tol = 1e-16` recovers
almost all of that (max 2.12x) at identical accuracy.  **`rel = 1e-14`,
`zero = 1e-16` is the configuration to use.**

Accuracy on this suite is identical for all three rules (worst `|norm - 1|` is
exactly 0 in every case).  That is expected rather than disappointing: at 2-12
qubits the scale factors are `O(1)`, far from the regime where the rules diverge.
The suite is evidence that the hybrid rule *costs* nothing, not that it *gains*
anything; the gains are in Results 1-3.

## Caveats

* `test_qmdd_gates` and `test_qmdd_matrix` fail on this branch -- they **also
  fail on pristine `master`**, identically (`test_qmdd_gates.c:428`,
  `test_pauli_rotation_gates`).  Not caused by, and not fixed by, this change.
* The relative rule shares the existing implementation's grid-boundary
  behaviour: two weights within tolerance can straddle adjacent buckets and fail
  to merge.  Probing neighbouring cells would reduce this for both rules.
* Only the `COMP_HASHMAP` backend is implemented.
* `qft_indep_qiskit_5.qasm` in the shipped MQT Bench archive is a **0-byte file**
  (1 of 217).  Both rules "disagree" on it -- absolute reports norm 0, hybrid
  norm 1 -- but this is a degenerate empty-circuit artefact, not a merging
  result, and it is excluded from the numbers above.
* Two circuits (`grover-noancilla_indep_qiskit_10`, `qwalk-noancilla_indep_qiskit_10`)
  exceeded the 60 s timeout under every rule alike.
* `qasm/CMakeLists.txt` needed two build fixes to compile the QASM front-end with
  current clang: `-Wno-missing-template-arg-list-after-template-kw` for the
  third-party `exprtk.hpp`, and the argp include/link that `examples/` already
  does.  Build-only; no semantic change.
* Build: `cmake -S . -B build -G Ninja -DCMAKE_POLICY_VERSION_MINIMUM=3.5`
  (the in-tree build dirs point at a cmake that no longer exists, and the QASM
  front-end's bundled `exprtk.hpp` does not compile with current clang, so build
  the `qsylvan`, `alg_run` and `hybrid_demo` targets specifically).
