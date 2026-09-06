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

## Caveats

* `test_qmdd_gates` and `test_qmdd_matrix` fail on this branch -- they **also
  fail on pristine `master`**, identically (`test_qmdd_gates.c:428`,
  `test_pauli_rotation_gates`).  Not caused by, and not fixed by, this change.
* The relative rule shares the existing implementation's grid-boundary
  behaviour: two weights within tolerance can straddle adjacent buckets and fail
  to merge.  Probing neighbouring cells would reduce this for both rules.
* Only the `COMP_HASHMAP` backend is implemented.
* Build: `cmake -S . -B build -G Ninja -DCMAKE_POLICY_VERSION_MINIMUM=3.5`
  (the in-tree build dirs point at a cmake that no longer exists, and the QASM
  front-end's bundled `exprtk.hpp` does not compile with current clang, so build
  the `qsylvan`, `alg_run` and `hybrid_demo` targets specifically).
