/*
 * Copyright 2026 System Verification Lab, LIACS, Leiden University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * Operations on LIMDDs: addition, and applying gates.
 *
 * COFACTORS ARE THE WHOLE DIFFICULTY
 *
 * Descending through an EVBDD node is free: the children are right there. A
 * LIMDD edge carries a LIM that acts on every qubit below it, so descending
 * means pushing that LIM through the node first. Its Pauli at this level
 * decides what happens to the two branches -- I and Z leave them in place, X
 * and Y exchange them -- and the rest of the LIM multiplies into both. That
 * is limdd_cofactors, and every operation here is built on it.
 *
 * WHY THE CACHE NEEDS THE LABELS DIVIDED OUT
 *
 * Keying the cache on the two edges as they arrive would almost never hit:
 * A|u> + B|v> and (2A)|u> + (2B)|v> are different keys for what is
 * essentially one computation. Both operations below therefore factor the
 * first operand's label out, cache on what remains, and multiply it back
 * afterwards. Addition also orders its two arguments, since it commutes.
 *
 * CONCURRENCY
 *
 * Lace tasks throughout, splitting on the two branches. The tables they touch
 * are lock-free and the operation cache is Sylvan's own, so nothing here
 * takes a lock.
 *
 * These do NOT protect their intermediate results from collection, so
 * limdd_gc must not run while one is in flight; see qsylvan_limdd_gc.h.
 */

#ifndef QSYLVAN_LIMDD_OPS_H
#define QSYLVAN_LIMDD_OPS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "qsylvan_limdd_canon.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * The two edges at level `var`+1 that `e`, an edge at level `var`, splits into.
 *
 * `e`'s own LIM is pushed through the node it points at: the Pauli at `var`
 * chooses and signs the branches, and what is left of the LIM multiplies into
 * both. `*low` is the |0> component and `*high` the |1> one.
 */
void limdd_cofactors(LIMDD e, uint32_t var, LIMDD *low, LIMDD *high);

/** `w` times `e`. Multiplying by zero gives the zero edge. */
LIMDD limdd_scale(LIMDD e, EVBDD_WGT w);

/**
 * The sum of two edges, both at level `var`.
 *
 * Must be called from a Lace worker; use RUN from outside one.
 */
TASK_DECL_3(LIMDD, limdd_plus, LIMDD, LIMDD, uint32_t);
#define limdd_plus(a, b, var) RUN(limdd_plus, a, b, var)

/**
 * `e` with the 2x2 gate `gateid` applied to qubit `target`.
 *
 * The gate is taken from Q-Sylvan's gate table, so the same identifiers work
 * here as for QMDDs, and the gate's entries are ordinary edge weights -- which
 * is what lets an exact backend stay exact through a circuit.
 */
TASK_DECL_4(LIMDD, limdd_gate, LIMDD, uint32_t, uint32_t, uint32_t);
#define limdd_gate(e, gateid, target, nqubits) \
    RUN(limdd_gate, e, gateid, target, nqubits)

/**
 * `e` with `gateid` applied to `target`, conditioned on every qubit whose bit
 * is set in `control_mask` being |1>.
 *
 * Every control must be ABOVE the target, i.e. control_mask < 2^target. The
 * recursion resolves a control by descending into its |1> branch, so one
 * below the target would still be unresolved when the target's branches are
 * mixed. Q-Sylvan's QMDD has the same restriction and falls back to a
 * matrix-vector product; LIMDD has no matrix representation yet, so callers
 * must reorder (the simulator has --reorder-swaps for this).
 */
TASK_DECL_5(LIMDD, limdd_cgate, LIMDD, uint32_t, uint64_t, uint32_t, uint32_t);
#define limdd_cgate(e, gateid, control_mask, target, nqubits) \
    RUN(limdd_cgate, e, gateid, control_mask, target, nqubits)

/**
 * The probability of measuring |1> on `qubit`, for a normalised state.
 *
 * Computed by summing squared magnitudes over the diagram rather than over
 * basis states, so it costs the size of the diagram and not 2^n.
 */
double limdd_prob_qubit_one(LIMDD e, uint32_t qubit, uint32_t nqubits);

/** The squared norm of `e`, which a normalised state has equal to 1. */
double limdd_norm_squared(LIMDD e, uint32_t nqubits);

/** The all-zero basis state on `nqubits` qubits. */
LIMDD limdd_all_zero_state(uint32_t nqubits);

#ifdef __cplusplus
}
#endif

#endif // QSYLVAN_LIMDD_OPS_H
