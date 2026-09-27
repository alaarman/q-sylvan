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
 * Operations on a BQD: the ones the paper proves.
 *
 * The cofactors of a pointwise product are the products of the cofactors,
 * and under the quotient rule the ratio of the products is the product of
 * the ratios with no cross term. With scalar labels on states of full
 * support that is the whole algorithm (prop:prodscalar): the node of fg has
 * low child the node of u.p and high child the node of v.q, the scalars
 * multiply, and memoised on node pairs the work is at most |f||g| entries.
 *
 * A diagonal gate is a product with a phase of full support, and one
 * monomial of it is a product with a diagram that is the constant one at
 * every level outside the monomial, so the recursion makes one call per
 * level (prop:diag). Z, S, T, CZ, CS, CCZ and their inverses are monomials.
 *
 * That is exactly the class of phase states, which is what the paper's main
 * size theorem is about, and it is the IQP fragment of a circuit: a layer of
 * Hadamards, then diagonal gates. What is NOT here, because the paper has no
 * algorithm for it: the product with translation labels in general (only the
 * corrected recursion of thm:prodx under its hypotheses), the Hadamard, and
 * any summation query. None of that is a limitation of this port.
 *
 * Everything here needs the scalar family, since the theorems are stated for
 * it, and full support, since a zero would make the copy fire and the ratio
 * of the products stop being the product of the ratios.
 */

#ifndef QSYLVAN_BQD_OPS_H
#define QSYLVAN_BQD_OPS_H

#include "qsylvan_bqd.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * The pointwise product f.g of two full-support diagrams of the scalar
 * family, canonical, as a Lace task memoised on node pairs. `var` is the
 * level both operands are read at: 0 for two root edges.
 */
TASK_DECL_3(BQD, bqd_product, BQD, BQD, uint32_t);
#define bqd_product(f, g, var) RUN(bqd_product, f, g, var)

/** The constant one on the variables from `var` down, an n-node chain. */
LIMDD_TARG bqd_ones(uint32_t var, uint32_t nqubits);

/**
 * The phase state phase^{x_{i_1} ... x_{i_d}} for the monomial on the qubits
 * of `A` (a vector-index mask), which has one node per level: the constant
 * one as a child everywhere outside A. This is the diagram of a diagonal
 * gate, and a product with it costs one call per level.
 */
BQD bqd_monomial(uint64_t A, EVBDD_WGT phase, uint32_t nqubits);

/** Apply the diagonal gate that multiplies by phase^{monomial on A}. */
BQD bqd_apply_diagonal(BQD e, uint64_t A, EVBDD_WGT phase, uint32_t nqubits);

#ifdef __cplusplus
}
#endif

#endif
