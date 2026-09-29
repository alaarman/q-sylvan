/*
 * Copyright 2025 Q-Sylvan contributors
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
 * Operations on a BQD: two of the five the paper proves, the two that act
 * on states of full support with scalar labels.
 *
 * The cofactors of a pointwise product are the products of the cofactors,
 * and under the quotient rule the ratio of the products is the product of
 * the ratios with no cross term. With scalar labels on states of full
 * support that is the whole algorithm (prop:prodscalar): the node of fg has
 * low child the node of u.p and high child the node of v.q, the scalars
 * multiply, and memoised on node pairs the work is at most |f||g| entries.
 *
 * A diagonal gate is a product with a phase of full support. For one
 * monomial the gate's operand has the constant one as a child at every
 * level, the levels it skips included, and a product with the constant one
 * returns the other operand, so the recursion is a single path of at most n
 * steps (prop:diag). bqd_apply_diagonal walks that path directly. Z, S, T,
 * CZ, CS, CCZ and their inverses are monomials.
 *
 * Diagrams here are fully reduced, as everywhere in the BQD: a level the
 * function does not depend on is skipped (skip:def:fr). On full support that
 * is one comparison, of a new node's ratio with the constant one
 * (skip:cor:full), so neither operation walks anything to decide it.
 *
 * Together they keep every state a phase state (def:level), the state in
 * the middle of an IQP circuit, which is what thm:size bounds: a Hadamard on
 * every qubit, then diagonal gates. Proved in the paper and NOT here: the
 * character rewrite on support-nested X-BQDs (prop:chi), the product on
 * aligned cosets with translation labels (thm:prodx), and the CNOT with its
 * control decided first as a relabelling (prop:cnot). The Hadamard and the
 * product on any support have no bounded algorithm in the paper; they are in
 * qsylvan_bqd_gates.h, exact and without a bound. A summation query is not
 * anywhere.
 *
 * Both need scalar labels and full support: a zero makes the copy fire and
 * the ratio of the products stop being the product of the ratios, and a
 * translation or a sign label is not a scalar; outside them they return a
 * wrong diagram rather than failing. The translation family has scalar
 * labels wherever the support is full, since every least point is then 0,
 * and its diagram is the scalar family's node for node, so both serve it
 * too. The Pauli family has sign labels there, and a product of two pivot
 * values in [0, pi) need not be in [0, pi), so a node of the product may
 * need the sign repair. So each tests full support of its whole input, a
 * walk memoised on nodes (bqd_has_full_support), and hands a state with a
 * zero to the general product of qsylvan_bqd_gates.h, which is correct on
 * any support and has no bound; bqd_product hands it every state of the
 * Pauli family as well, and bqd_apply_diagonal takes a state of full support
 * in the Pauli family to bqd_xp_diagonal, the same walk with label products
 * and the sign repair at the last qubit of A (skip:alg:xdiag).
 */

#ifndef QSYLVAN_BQD_OPS_H
#define QSYLVAN_BQD_OPS_H

#include "qsylvan_bqd.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * The pointwise product f.g of two diagrams, canonical, as a Lace task
 * memoised on node pairs: prop:prodscalar when both have full support in the
 * scalar or translation family, and bqd_multiply otherwise. The support of
 * both is tested once, here, and the recursion below takes it as given.
 */
TASK_DECL_2(BQD, bqd_product, BQD, BQD);
#define bqd_product(f, g) RUN(bqd_product, f, g)

/**
 * The constant one on the variables from `var` down: the terminal, since it
 * depends on none of them and every level is skipped. The edge of |+>^n up to
 * scale is limdd_one_edge().
 */
LIMDD_TARG bqd_ones(uint32_t var, uint32_t nqubits);

/**
 * The phase state phase^{x_{i_1} ... x_{i_d}} for the monomial on the qubits
 * of `A` (a vector-index mask): one node at each qubit of A and none
 * elsewhere, a path whose low edges are the constant one, and the terminal
 * for a phase of one; in the Pauli family a phase of -1 is a Z at the last
 * qubit of A and no node there. The diagram of a diagonal gate, for
 * bqd_product; bqd_apply_diagonal does not need it on full support.
 */
BQD bqd_monomial(uint64_t A, EVBDD_WGT phase, uint32_t nqubits);

/**
 * Multiply e by phase^{x_A}, the diagonal gate of the monomial on the qubits
 * of `A` (prop:diag). On full support the nodes of e visited are one path,
 * at most one per level down to the last qubit of A, n - ctz(A) at most,
 * whatever |e| is; the walk stops early where e skips a qubit of A
 * (skip:prop:diag). At most one new node per level is made, after the test
 * for full support, which is memoised and so costs only the nodes it has not
 * seen. Without full support it is bqd_multiply with the monomial, and
 * visits nothing. In the Pauli family the walk is bqd_xp_diagonal, which
 * visits the same path. The result is the canonical diagram, the same edge
 * as bqd_product(e, bqd_monomial(...)). From a Lace worker.
 */
BQD bqd_apply_diagonal(BQD e, uint64_t A, EVBDD_WGT phase, uint32_t nqubits);

/** As bqd_apply_diagonal, and set *visits to the number of nodes of e visited. */
BQD bqd_apply_diagonal_counted(BQD e, uint64_t A, EVBDD_WGT phase, uint32_t nqubits,
                               uint32_t *visits);

#ifdef __cplusplus
}
#endif

#endif
