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
 * Exponent diagrams, for the BQD's phase multiplications (skip:sec:ratio:phase
 * of the note "Binary Quotient Diagrams with Level Skipping", bqd-skip.tex).
 *
 * A diagonal gate multiplies a function by phi^{pi_A}, and the multiplication
 * off a support that Rebuild needs multiplies it by c^{[x not in supp U]}.
 * Both are beta^eps for a scalar beta and an integer function eps, and
 * PhaseMul recurses on a node and such an exponent together, one level at a
 * time. The exponent is held as a reduced, ordered diagram under the Shannon
 * rule with integer terminals: a Sylvan MTBDD with int64 leaves, whose
 * variable v decides qubit v as a BQD node's does, so that the two recursions
 * meet at the same variables. Its values are taken modulo the order r of beta
 * when beta is a root of unity of order at most 8, the only roots the exact
 * weights have, and are plain integers (r = 0) otherwise. The reduction is
 * what keeps the diagram of beta^eps free of a node whose two exponents differ
 * by a multiple of r, which Exp would otherwise make into a node with ratio 1.
 * On floats the order test can fail for a root of unity, and the exponents
 * are then integers, which is correct and only shares less.
 *
 * Every operation is memoised in the operation cache on its operands and r,
 * under ids of its own, and is a Lace task where it recurses. A point y is a
 * mask of qubits, bit v for qubit v, as the translation and Pauli labels are.
 *
 * Collections. The diagrams live in Sylvan's node table, which Sylvan collects
 * when it fills, in the middle of an operation if need be, and a collection
 * keeps only what a Sylvan root reaches. So every exponent a caller holds
 * across a call that can make a node goes on the mtbdd_refs stack, as in
 * Sylvan's own operations; the operations here protect what they make. A
 * Sylvan collection also clears the operation cache, which is what keeps a
 * memo entry from naming an exponent that is gone. The LIMDD's collection
 * runs between gates only and leaves Sylvan's table alone.
 */

#ifndef QSYLVAN_BQD_EXP_H
#define QSYLVAN_BQD_EXP_H

#include <stdbool.h>
#include <stdint.h>

#include <sylvan_int.h>
#include "sylvan_edge_weights.h"

#ifdef __cplusplus
extern "C" {
#endif

/** An exponent diagram: an MTBDD with int64 leaves, variable v for qubit v. */
typedef MTBDD BQD_EXP;

/** The order of beta as a root of unity, 1 to 8, or 0 if beta^k is not 1 for k <= 8. */
uint32_t bqd_exp_order(EVBDD_WGT beta);

/** beta^v, for any integer v; 1/beta^-v for a negative one. */
EVBDD_WGT bqd_exp_power(EVBDD_WGT beta, int64_t v);

/** The constant v, reduced modulo r when r > 0 into [0, r). */
BQD_EXP bqd_exp_const(int64_t v, uint32_t r);

static inline bool bqd_exp_is_const(BQD_EXP e) { return mtbdd_isleaf(e) != 0; }

/** The value of a constant. */
int64_t bqd_exp_value(BQD_EXP e);

/** The variable of e's root, or UINT32_MAX for a constant, which is below every qubit. */
uint32_t bqd_exp_var(BQD_EXP e);

/** The cofactor x_var = b of e read at var, at or above its root: e itself where it skips var. */
BQD_EXP bqd_exp_cof(BQD_EXP e, uint32_t var, int b);

/** e at the point y, bit v of y for qubit v. */
int64_t bqd_exp_eval(BQD_EXP e, uint64_t y);

/** a - b and a + b, leafwise, modulo r. */
TASK_DECL_3(BQD_EXP, bqd_exp_sub, BQD_EXP, BQD_EXP, uint32_t);
#define bqd_exp_sub(a, b, r) RUN(bqd_exp_sub, a, b, r)
TASK_DECL_3(BQD_EXP, bqd_exp_add, BQD_EXP, BQD_EXP, uint32_t);
#define bqd_exp_add(a, b, r) RUN(bqd_exp_add, a, b, r)

/**
 * Sel(iota, x, v): x where the 0/1 exponent iota is 1, and the constant v,
 * already reduced, where it is 0.
 */
TASK_DECL_3(BQD_EXP, bqd_exp_sel, BQD_EXP, BQD_EXP, int64_t);
#define bqd_exp_sel(iota, x, v) RUN(bqd_exp_sel, iota, x, v)

/** X^t e, e read at y xor t: the two children swapped at every qubit of t. */
TASK_DECL_2(BQD_EXP, bqd_exp_translate, BQD_EXP, uint64_t);
#define bqd_exp_translate(e, t) RUN(bqd_exp_translate, e, t)

/** gamma . prod_{v in vars} x_v, modulo r: one node per qubit of vars, and gamma alone for none. */
BQD_EXP bqd_exp_monomial(uint64_t vars, int64_t gamma, uint32_t r);

/** P_s = 4 (s . x) modulo 8, so that w_8^{P_s} = (-1)^{s . x}: two nodes per qubit of s at most. */
BQD_EXP bqd_exp_parity4(uint64_t s);

/** The exponent node (var; lo, hi), reduced: lo itself where the two agree. */
BQD_EXP bqd_exp_node(uint32_t var, BQD_EXP lo, BQD_EXP hi);

/**
 * For the tests: run a Sylvan collection at every k-th exponent node that the
 * calling worker makes from now on, so that collections happen inside the
 * operations; 0 turns it off, which is the default. Only that worker starts
 * them: two workers that ask sylvan_gc at once can wait on each other.
 */
void bqd_exp_collect_every(uint64_t k);

#ifdef __cplusplus
}
#endif

#endif
