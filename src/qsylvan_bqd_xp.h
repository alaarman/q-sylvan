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
 * The operations of the translation and Pauli families, the X-BQD and the
 * Pauli-BQD, whose labels are c X^t and c Z^s X^t (skip:sec:xp of the note
 * "Binary Quotient Diagrams with Level Skipping", bqd-skip.tex).
 *
 * The scheme is the scalar family's (qsylvan_bqd_gates.h): a pointwise
 * operation acts on each cofactor separately, so it recurses on the cofactors
 * of its arguments and rebuilds its result with Compose. Three things change.
 *
 * Labelled edges. The cofactor of an edge is a label product, and a label
 * product is in general not the canonical edge of the function it denotes, so
 * the operations take LABELLED edges as arguments (skip:def:labelled): any
 * label on a node of a canonical diagram, interned in the LIM table like any
 * other label. Only what they return is canonical. Where the scalar family
 * returns an operand or a cofactor as it is, because it is canonical there,
 * the edge goes through Canon here (skip:prop:xcanon). A scalar multiple of a
 * canonical edge is canonical in every family (skip:lem:xscale), so a scalar
 * still leaves a key as it does in the scalar family.
 *
 * The top. A Z at a level above an edge's node makes its function depend on
 * that level: (Z, terminal) is (1, -1). So a recursion runs at the TOP of its
 * arguments, the higher of the node's level and the highest level with a Z,
 * and not at the level of their nodes. An X above the node acts on nothing,
 * and Norm clears it. An edge skips a level when its top is below it.
 *
 * Compose finds the representative of the node it makes (skip:prop:xcompose):
 * the least point of the high cofactor's support, by a walk down one path
 * (skip:lem:minpoint), and in the Pauli family the sign of the level's pivot,
 * one test of an argument. Two cofactors that are equal skip the level, and
 * in the Pauli family so do two that are opposite, with a Z on the label.
 *
 * Masks here are the LIM's own, bit q for qubit q, and not the vector-index
 * masks of bqd_lim_make, so that a variable is a bit: var v is bit v, and the
 * qubits above it are the bits below. Products and inverses of labels are
 * those of skip:sec:xp:edges, (ZX)(ZX) = -1, and NOT limdd_lim_mul's, which
 * multiplies Pauli matrices and reads x = z = 1 as a Y that squares to I. The
 * two agree in the translation family only.
 *
 * Memo keys. A labelled edge goes into a key as the index of its label's
 * Pauli word beside its node, in one word, and the scalar that the operation
 * lets stay in the key goes in the first word (skip:sec:xp:keys), so forming
 * a key interns nothing. An operation that misses recurses on its operands
 * as they are, scalars and all, since a pointwise operation commutes with a
 * cofactor whatever the scalars, returns what that gives, and stores it
 * divided by the scalar that left the key. It pushes its operands onto the
 * refs stack, as it reads their labels again after a point where a
 * collection can run.
 *
 * Diagonal gates. On full support the translation family's diagram is the
 * scalar family's node for node, and the scalar walk serves it. The Pauli
 * family has a walk of its own, bqd_xp_diagonal (skip:alg:xdiag): the same
 * path, with label products where the scalar walk multiplies scalars, since
 * on full support a Z toggles one pivot sign and no other, and the sign
 * repair and the skip at the levels of A, the repair at the last one only.
 *
 * Exact weights are what these are for. On floats they run, and after a
 * cancellation a rounding error can decide the representative's comparisons
 * with zero and with [0, pi), and Compose's test b = -a, the wrong way: a
 * wrong sign costs sharing, a wrong zero divides by a residue, which can
 * overflow to a non-finite amplitude.
 *
 * Limits (skip:sec:xp:limits). The operations are exact and carry no bound,
 * as the scalar family's Hadamard. They are memoised on pairs of a label and
 * a node, one node can meet as many labels as label products reach it, and
 * every label product interns a LIM. The paper's specialised recursions on
 * coset states, the character rewrite (prop:chi), the product with a relative
 * translation (thm:prodx) and the CNOT as a relabelling (prop:cnot), are not
 * what is implemented here.
 *
 * What that costs, measured on one worker. A gate's cost follows the
 * labelled cofactors it meets, which the diagram does not hold, and not the
 * diagram's size, so fewer nodes are not a faster run. On a 16-qubit
 * Clifford+T circuit of 500 gates the scalar family takes 2.4 s and ends with
 * 383 nodes, the translation family 3.6 s with 151, the Pauli family 2.9 s
 * with 111; on a 20-qubit Clifford circuit of 400 gates 8.5 s with 714
 * nodes, 8.3 s with 129 and 3.3 s with 28. On 20-qubit Clifford+T circuits
 * of 700 gates all three run for minutes. The Pauli family also mints more:
 * about eleven LIMs per node it mints, new Pauli words at a third of the
 * rate of new LIMs, which the runner sizes its tables for, and more weights:
 * a 10-qubit Clifford+T circuit of 149 gates that the other two families
 * run in a weight table of 2^20 needs 2^21 in the Pauli family.
 *
 * Restrict, Project and Gate are the scalar family's at the top of their
 * argument, with Canon on the operand or cofactor the scalar family returns
 * as it is, Mono is built with Compose, and the basis state is the chain of
 * |0...0> under a translation; every other operation of qsylvan_bqd_gates.h
 * is made of these, as in the scalar family.
 *
 * The entry points of qsylvan_bqd_gates.h and qsylvan_bqd_ops.h dispatch here
 * on bqd_family(); the scalar family never comes here. From a Lace worker.
 */

#ifndef QSYLVAN_BQD_XP_H
#define QSYLVAN_BQD_XP_H

#include "qsylvan_bqd.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * The four pointwise operations: the product and the sum a state needs, and
 * the two the quotient rule is made of. XPROD(a, b) is a.b where a is nonzero
 * and b where it is zero, the (.) of def:bqd, and XQUOT(a, b) is a/b where b
 * is nonzero and a where it is zero. The numbers are those of the scalar
 * family's operations.
 */
enum { BQD_OP_MUL = 0, BQD_OP_ADD = 1, BQD_OP_XPROD = 2, BQD_OP_XQUOT = 3 };

/**
 * op(f, g) of two labelled edges, as the canonical edge: Apply of
 * skip:alg:xapply. The operands may be any labelled edges, normal or not; the
 * result is the edge bqd_from_vector builds for the pointwise result.
 */
TASK_DECL_3(BQD, bqd_xp_apply_op, int, BQD, BQD);
#define bqd_xp_apply_op(op, f, g) RUN(bqd_xp_apply_op, op, f, g)

/**
 * The canonical edge, at level `var`, of the function whose cofactors are
 * f_0 = lo and f_1 = hi, two CANONICAL edges read at var + 1: Compose of
 * skip:alg:xcompose.
 */
TASK_DECL_3(BQD, bqd_xp_compose, uint32_t, BQD, BQD);
#define bqd_xp_compose(var, lo, hi) RUN(bqd_xp_compose, var, lo, hi)

/**
 * The cofactor x_var = b of the labelled edge e read at `var`, which must be
 * at or above its top, as a labelled edge read at var + 1: Cof of
 * skip:alg:xpedges. NOT canonical in general; bqd_xp_canon makes it so, which
 * is what bqd_cofactor returns.
 */
TASK_DECL_3(BQD, bqd_xp_cofactor, BQD, uint32_t, int);
#define bqd_xp_cofactor(e, var, b) RUN(bqd_xp_cofactor, e, var, b)

/** The canonical edge of the function a labelled edge denotes: Canon of skip:alg:xcompose. */
TASK_DECL_1(BQD, bqd_xp_canon, BQD);
#define bqd_xp_canon(e) RUN(bqd_xp_canon, e)

/**
 * The canonical edge of [N_0] (.) [N_1], for N_0 and N_1 the children of the
 * node N: Join of skip:alg:xpedges. The high cofactor of N is the high label
 * applied to it, outside the product.
 */
TASK_DECL_1(BQD, bqd_xp_join, LIMDD_TARG);
#define bqd_xp_join(N) RUN(bqd_xp_join, N)

/** c . e: the scalar of the label multiplied by c, canonical when e is. */
BQD bqd_xp_scale(BQD e, EVBDD_WGT c);

/**
 * The top of a labelled edge as a variable: the smaller of its node's
 * variable (limdd_level for the terminal) and the smallest qubit with a Z.
 * The edge may be read at any variable at or above it, and skips every
 * variable above it. The zero edge's top is the terminal's.
 */
uint32_t bqd_xp_top(BQD e);

/**
 * e|_{x_q = b} and e . [x_q = b] of a labelled edge, canonical: Restrict and
 * Project of skip:alg:xrestrict, which bqd_restrict and bqd_project are in
 * these families.
 */
TASK_DECL_3(BQD, bqd_xp_restrict, BQD, uint32_t, int);
#define bqd_xp_restrict(e, q, b) RUN(bqd_xp_restrict, e, q, b)
TASK_DECL_3(BQD, bqd_xp_project, BQD, uint32_t, int);
#define bqd_xp_project(e, q, b) RUN(bqd_xp_project, e, q, b)

/**
 * The gate `gateid` (an index into gates[]) on qubit q of a labelled edge,
 * conditioned on the qubits of `cmask`, bit c for qubit c, all above q, being
 * 1: Gate of skip:alg:xapply, canonical. bqd_gate and bqd_cgate come here in
 * these families for every gate but a monomial phase.
 */
TASK_DECL_4(BQD, bqd_xp_cgate_rec, BQD, uint32_t, uint64_t, uint32_t);
#define bqd_xp_cgate_rec(e, gateid, cmask, q) RUN(bqd_xp_cgate_rec, e, gateid, cmask, q)

/** |x>, x a vector index: X^x on the n nodes of |0...0>, which every x shares. */
BQD bqd_xp_basis_state(uint64_t x, uint32_t nqubits);

/**
 * phase^{x_A} for the qubits of the vector-index mask A: Mono of
 * skip:alg:xrestrict, at most one node per qubit of A, fewer in the Pauli
 * family, where a factor -1 is a Z on a label. From a Lace worker.
 */
BQD bqd_xp_monomial(uint64_t A, EVBDD_WGT phase, uint32_t nqubits);

/**
 * e . phase^{x_A} in the Pauli family, for e of FULL SUPPORT, canonical: the
 * diagonal walk of skip:alg:xdiag, the scalar family's walk (bqd_apply_diagonal)
 * with label products for its scalar products and a sign repair and a skip
 * at the levels of A. One path, n steps at most, no memo; *visits, when not
 * NULL, gets the nodes of e it visits, the path of skip:prop:diag. The
 * translation family needs none: on full support its labels are scalars and
 * the scalar walk serves it. From a Lace worker.
 */
BQD bqd_xp_diagonal(BQD e, uint64_t A, EVBDD_WGT phase, uint32_t nqubits, uint32_t *visits);

#ifdef __cplusplus
}
#endif

#endif
