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
 * LEVELS ARE CONTEXT
 *
 * An edge does not say what level it is read at, and once a diagram is fully
 * reduced an edge may skip levels -- point from level k at a node of level
 * k' > k, each skipped level holding |0>+|1> under whatever Pauli the label
 * puts there (see qsylvan_limdd_node.h). limdd_cofactors handles that: at a
 * skipped level both cofactors are the edge itself, and the four-way Pauli
 * dispatch it already has does the rest.
 *
 * Two disciplines follow. Operations whose level the operands cannot supply
 * -- addition, the norm, a probability -- carry the level as a parameter and,
 * where they cache, in the cache key. A gate recovers its working level from
 * the edge and the gate's own qubits (the first level that matters: the
 * target's node, or the gate qubit or highest control if the edge skips
 * it), hoists the label's entries above that level out of the recursion,
 * and can therefore keep a level-free key. And a label must act only at or
 * below its edge's level; every argument above rests on that.
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
 *
 * `var` may lie above the target's level, in which case `e` skips it and the
 * two cofactors are the same edge up to the sign or phase that the label's
 * Pauli at `var` gives |0>+|1>. `var` may not lie below the target's level,
 * and the label must act only at or below `var`.
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
 * A controlled gate with the control on either side of the target.
 *
 * limdd_cgate needs the control above the target. When it is not, two
 * identities put it there rather than reaching for a matrix representation
 * LIMDD does not have:
 *
 *   CZ is symmetric, so its two qubits simply exchange roles.
 *   CX reverses under Hadamards on both qubits:
 *   (H (x) H) CNOT_{a->b} (H (x) H) = CNOT_{b->a}.
 *
 * Both are exact and stay inside Clifford+T. Any other gate with the control
 * below the target is refused.
 */
LIMDD limdd_cgate_either(LIMDD e, uint32_t gateid, uint32_t control,
                         uint32_t target, uint32_t nqubits, bool *ok);

/**
 * Exchange two qubits, as three CNOTs, the middle one reversed by the
 * identity above.
 */
LIMDD limdd_swap(LIMDD e, uint32_t a, uint32_t b, uint32_t nqubits);

/**
 * The probability of measuring |1> on `qubit`, for a normalised state.
 *
 * Computed by summing squared magnitudes over the diagram rather than over
 * basis states, so it costs the size of the diagram and not 2^n.
 */
double limdd_prob_qubit_one(LIMDD e, uint32_t qubit, uint32_t nqubits);

/** The squared norm of `e`, which a normalised state has equal to 1. */
double limdd_norm_squared(LIMDD e, uint32_t nqubits);

/**
 * How many nodes the diagram rooted at `e` actually uses.
 *
 * Not the same as limdd_node_table_count, which counts everything the table
 * has ever held; this walks the reachable nodes, so it is comparable with
 * evbdd_countnodes on the QMDD side.
 */
uint64_t limdd_countnodes(LIMDD e);

/** The all-zero basis state on `nqubits` qubits. */
LIMDD limdd_all_zero_state(uint32_t nqubits);

#ifdef __cplusplus
}
#endif

#endif // QSYLVAN_LIMDD_OPS_H
