/*
 * Copyright 2026 Q-Sylvan contributors
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
 * Rule SM, Shannon where misaligned, for the BQD's translation and Pauli
 * families (skip:sec:sm:xp of the note "Binary Quotient Diagrams with Level
 * Skipping", bqd-skip.tex, whose labels are cited as skip:...).
 *
 * The scalar family's rule (qsylvan_bqd_sm.h) with the labels c X^t and
 * c Z^s X^t of qsylvan_bqd_xp.h. A node has the identity on its low edge, to
 * A, which is never zero in these families, since 0 is in the support of
 * every stored node, and its high edge (l_1, R) is read by its tag
 * (skip:def:smxp):
 *
 *     Q (nested)    f_0 = [A],  f_1 = l_1 . ([A] . [R]),  supp R inside supp A
 *     S (Shannon)   f_0 = [A],  f_1 = l_1 . [R]
 *
 * The X part t_1 of l_1 is the least point of supp f_1, and the tag is decided
 * AFTER the high cofactor is moved by it: Q exactly where supp X^{t_1} f_1
 * lies inside supp f_0, which on the stored nodes is supp R inside supp A.
 * The ratio of a Q node is then X^{t_1} f_1 / f_0 on supp f_0 and zero off
 * it, an S node stores the canonical edge of f_1, and neither has a copy. A
 * Pauli node has the tag of the translation family's node of its function,
 * since the two representatives differ by signs and not by supports. The
 * product is the plain one, a level is skipped exactly where f_0 = f_1, and in
 * the Pauli family with a Z where f_1 = -f_0, as under the copy rule.
 *
 * The support of a node is {0} x supp [A] and {1} x (supp [R] ^ t_1) for both
 * tags, the copy rule's shadow, so the least point (minpoint), the indicator,
 * the subset test and the dependence test read the stored edges and do no
 * arithmetic. On every function of full support and every function whose
 * support is an affine subspace the diagram is the copy rule's, node for node
 * (skip:prop:smxpd0); elsewhere a node the copy rule would mix is an S node.
 *
 * One thing does not carry over from the scalar family. Since the tag is
 * decided after the move by the least point, a map that moves that point can
 * change it: a relabelling or a permutation of the lower levels can turn an S
 * node into a Q node, or the reverse (skip:lem:smxptag), and the rebuild of
 * an S node keeps the tag only where the new least point is the image of the
 * old one, and goes through Compose otherwise.
 *
 * Nodes are the scalar SM's: bqd_sm_mk makes every node of such a session,
 * with the tag S and the full flag in the node's flag bits, so bqd_sm_is_s and
 * bqd_sm_full read them, and the full flag answers bqd_has_full_support in
 * O(1). The copy rule's constructors never run in such a session: the public
 * entry points ask bqd_xp_sm() and come here.
 *
 * The scheme is the copy rule's of qsylvan_bqd_xp.h: the operations take
 * labelled edges and return canonical ones, run at the top of their
 * arguments, and are memoised on labelled edges as key words (key_edge) with
 * the scalar out, under operation ids of their own (sylvan_int.h). Unlike the
 * copy rule's, they recurse without that scalar, as the scalar family's do,
 * and put it on the result once, which on exact weights saves a division of
 * two large numbers per miss. What changes is how a node is made and read.
 * Compose moves the high cofactor by its least point, and makes a Q node with
 * the ratio where the shadow says the supports nest, and an S node with the
 * canonical edge of the moved cofactor where they do not, which costs one
 * Canon and no quotient. A Q node read at its level is its VIEW
 * (skip:lem:smxpview), f_0 = F_0 and f_1 = X^t (F_0 . rho) with supp rho
 * inside supp F_0, and an operation that keeps the translation t acts on F_0
 * and rho alone and has NEST make the node, with no high cofactor and no
 * quotient where the least point stays: the product and the ratio of two
 * nested operands with one translation, the sum of two with one ratio, and
 * Canon, which moves a node by a label.
 *
 * What is here, skip:sec:sm:xp:alg: the builder, the decoder and the
 * amplitude query; the shadow walks Ind, Sub with a translation and Dep; the
 * cofactors and Cof1; the constructors Compose, View, Nest and the Shannon
 * rebuild; Canon; Apply with the product, the ratio and the sum; the
 * selections Restrict, Project, X, Pair and Perm, which NEST a Q node's view
 * where its translation leaves their qubits alone; the gate, with the scalar
 * family's one-sided gate (Side), a restriction and a phase on the gate's
 * qubit where the support lies on one side of it, and the pass-through of a
 * ratio that does not depend on that qubit; the phase multiplications of the
 * two families, Exp, MulOff and the Pauli family's diagonal walk on full
 * support, which keep every support and so every tag and least point; the
 * basis states and the monomials. Every public entry point of
 * qsylvan_bqd_gates.h and qsylvan_bqd_ops.h comes here in such a session.
 *
 * Every operation is a Lace task, spawns one of its two independent calls,
 * pushes what it holds across a call onto the refs stack as the copy rule's
 * do, and counts its misses for bqd_cache_fit. From a Lace worker.
 */

#ifndef QSYLVAN_BQD_XP_SM_H
#define QSYLVAN_BQD_XP_SM_H

#include "qsylvan_bqd.h"
#include "qsylvan_bqd_exp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* --- the support shadow ---------------------------------------------------- */

/**
 * The node of the indicator of supp [N], with the X parts of N's high labels
 * and N's tags: 0 for 0, the terminal on full support.
 */
TASK_DECL_1(LIMDD_TARG, bqd_xpsm_ind, LIMDD_TARG);
#define bqd_xpsm_ind(N) RUN(bqd_xpsm_ind, N)

/** Whether supp [X] ^ tau lies inside supp [Y]; X and Y nodes, the terminal, or 0 for zero. */
TASK_DECL_3(int, bqd_xpsm_sub, LIMDD_TARG, LIMDD_TARG, uint64_t);
#define bqd_xpsm_sub(X, Y, tau) RUN(bqd_xpsm_sub, X, Y, tau)

/** Whether a node of variable q, or a Z_q on a high label, is reachable from N, a node or 0. */
TASK_DECL_2(int, bqd_xpsm_dep, LIMDD_TARG, uint32_t);
#define bqd_xpsm_dep(N, q) RUN(bqd_xpsm_dep, N, q)

/**
 * Side: whether supp [N], N a node or the terminal, lies in {x_q = 0}, the
 * only side the support of a stored node can lie on, since it has 0 in it.
 */
TASK_DECL_2(int, bqd_xpsm_side, LIMDD_TARG, uint32_t);
#define bqd_xpsm_side(N, q) RUN(bqd_xpsm_side, N, q)

/* --- cofactors ------------------------------------------------------------- */

/**
 * The canonical edge of [A] . [R] for a Q node N with a nonzero high edge,
 * without the high label, which goes on outside: Cof1. An S node needs none,
 * since its high cofactor is its stored edge.
 */
TASK_DECL_1(BQD, bqd_xpsm_cof1, LIMDD_TARG);

/**
 * The cofactor x_var = b of the labelled edge e read at var, at or above its
 * top, as a labelled edge read at var + 1, not canonical in general:
 * bqd_xp_cofactor with the high side of a node read by its tag.
 */
TASK_DECL_3(BQD, bqd_xpsm_cofactor, BQD, uint32_t, int);
#define bqd_xpsm_cofactor(e, var, b) RUN(bqd_xpsm_cofactor, e, var, b)

/* --- constructors ---------------------------------------------------------- */

/**
 * The canonical edge at level var of the function with the cofactors lo and
 * hi, two canonical edges read at var + 1: Compose of skip:alg:smxpcons.
 */
TASK_DECL_3(BQD, bqd_xpsm_compose, uint32_t, BQD, BQD);
#define bqd_xpsm_compose(var, lo, hi) RUN(bqd_xpsm_compose, var, lo, hi)

/**
 * The view of a nonzero labelled edge F read at level var, at or above its
 * top: the edges *lo and *rho and the translation *t with F_0 = [*lo],
 * F_1 = X^t ([*lo] . [*rho]) and supp *rho inside supp *lo, the ratio carrying
 * none of F's scalar. Defined, and 1, where F's node is below var, as its
 * virtual Q node (N, Ind N), and where it is a Q node at var whose label has
 * no X at var; 0 otherwise. *rho is canonical where its translation is 0.
 */
TASK_DECL_5(int, bqd_xpsm_view, BQD, uint32_t, BQD *, uint64_t *, BQD *);

/**
 * The canonical edge at level var of the function with the cofactors [lo]
 * and X^t ([lo] . [rho]), for a canonical lo and a labelled rho with
 * supp rho inside supp lo, canonical when rho_canon: a Q node with no product
 * and no quotient where the least point of the high cofactor is t, the level
 * skipped where the ratio is the indicator of lo's support, and Compose of the
 * product otherwise.
 */
TASK_DECL_5(BQD, bqd_xpsm_nest, uint32_t, BQD, uint64_t, BQD, int);

/**
 * The node at level var from the canonical low cofactor lo and a labelled
 * high cofactor f1 of the image of an S node under an affine bijection of the
 * levels below var, tt the image of the old node's high translation under its
 * linear part: an S node with one Canon where the least point of f1 is tt,
 * and Compose's general case where it moved.
 */
TASK_DECL_4(BQD, bqd_xpsm_rebuild_s, uint32_t, BQD, BQD, uint64_t);

/** The canonical edge of the function a labelled edge denotes: Canon. */
TASK_DECL_1(BQD, bqd_xpsm_canon, BQD);
#define bqd_xpsm_canon(e) RUN(bqd_xpsm_canon, e)

/* --- pointwise operations -------------------------------------------------- */

/**
 * op(f, g) of two labelled edges, canonical: Apply of skip:alg:smxpapply, op
 * BQD_SM_MUL, BQD_SM_RATIO or BQD_SM_ADD (qsylvan_bqd_sm.h). RATIO(f, g) is
 * g / f on supp f and 0 off it, and is defined only where supp g lies inside
 * supp f, which the caller guarantees.
 */
TASK_DECL_3(BQD, bqd_xpsm_apply, int, BQD, BQD);
#define bqd_xpsm_apply(op, f, g) RUN(bqd_xpsm_apply, op, f, g)

/** rho . [Y != 0], for Y nonzero: rho where its support lies inside Y's already. */
TASK_DECL_2(BQD, bqd_xpsm_restr, BQD, BQD);
#define bqd_xpsm_restr(rho, Y) RUN(bqd_xpsm_restr, rho, Y)

/* --- selections ------------------------------------------------------------- */

/*
 * The signatures are those of bqd_restrict, bqd_project, bqd_x, bqd_pair and
 * bqd_perm (qsylvan_bqd_gates.h), on labelled edges, with canonical results.
 */

/** e|_{x_q = b}, as a function that ignores x_q. */
TASK_DECL_3(BQD, bqd_xpsm_restrict, BQD, uint32_t, int);
#define bqd_xpsm_restrict(e, q, b) RUN(bqd_xpsm_restrict, e, q, b)

/** e . [x_q = b]. */
TASK_DECL_3(BQD, bqd_xpsm_project, BQD, uint32_t, int);
#define bqd_xpsm_project(e, q, b) RUN(bqd_xpsm_project, e, q, b)

/** X on qubit q: Canon of X_q . e. */
TASK_DECL_2(BQD, bqd_xpsm_x, BQD, uint32_t);
#define bqd_xpsm_x(e, q) RUN(bqd_xpsm_x, e, q)

/** [x_b = 0] A|_{x_b = s} + [x_b = 1] B|_{x_b = t}. */
TASK_DECL_5(BQD, bqd_xpsm_pair, BQD, int, BQD, int, uint32_t);
#define bqd_xpsm_pair(A, s, B, t, b) RUN(bqd_xpsm_pair, A, s, B, t, b)

/** The permutation `kind` (BQD_PERM_*) of the qubits qa < qb: Canon(l^pi . PermN(N)). */
TASK_DECL_4(BQD, bqd_xpsm_perm, BQD, uint32_t, uint32_t, uint32_t);
#define bqd_xpsm_perm(e, kind, qa, qb) RUN(bqd_xpsm_perm, e, kind, qa, qb)

/* --- gates ------------------------------------------------------------------ */

/**
 * gates[gid] on qubit q of a labelled edge, conditioned on the controls in
 * cmask (bit c for qubit c, every one above q) that are still pending: the
 * copy rule's Gate at the top of e, with the scalar family's one-sided gate,
 * u_0s e|_{x_q = s} . (u_1s / u_0s)^{x_q} where no control is pending and the
 * support lies in {x_q = s}, and its pass-through of a Q node whose ratio
 * does not depend on x_q and whose translation leaves q and the controls
 * alone. bqd_gate and bqd_cgate come here for every gate but a monomial phase
 * and a flip.
 */
TASK_DECL_4(BQD, bqd_xpsm_gate, BQD, uint32_t, uint64_t, uint32_t);
#define bqd_xpsm_gate(e, gid, cmask, q) RUN(bqd_xpsm_gate, e, gid, cmask, q)

/* --- phases and diagonal gates ---------------------------------------------- */

/** PhaseMulX: [N] . beta^eps in the translation family, r the order of beta. */
TASK_DECL_4(BQD, bqd_xpsm_phase_mul_x, LIMDD_TARG, BQD_EXP, EVBDD_WGT, uint32_t);

/** PhaseMulP: [N] . w_8^eps in the Pauli family, eps modulo 8. */
TASK_DECL_2(BQD, bqd_xpsm_phase_mul_p, LIMDD_TARG, BQD_EXP);

/**
 * e . beta^eps for a labelled edge e, bqd_phase_mul's: the label times
 * PhaseMulX of the translated exponent, or in the Pauli family, for
 * beta = w_8^m, times PhaseMulP of m eps; the Pauli family exits on any other
 * beta. Canonical where e is.
 */
TASK_DECL_3(BQD, bqd_xpsm_phase_mul, BQD, BQD_EXP, EVBDD_WGT);
#define bqd_xpsm_phase_mul(e, eps, beta) RUN(bqd_xpsm_phase_mul, e, eps, beta)

/** beta^eps, of full support, bqd_exp_state's; in the Pauli family beta a power of w_8. */
TASK_DECL_3(BQD, bqd_xpsm_exp_state, BQD_EXP, EVBDD_WGT, uint32_t);

/** MulOffX, the translation family: k on supp [U] and c . k off it; U a node, the terminal, or 0. */
TASK_DECL_3(BQD, bqd_xpsm_mul_off, BQD, EVBDD_WGT, LIMDD_TARG);
#define bqd_xpsm_mul_off(k, c, U) RUN(bqd_xpsm_mul_off, k, c, U)

/**
 * e . phase^{x_A} in the Pauli family for e zero or of FULL SUPPORT, A a
 * vector-index mask: the copy rule's walk (bqd_xp_diagonal), whose diagram on
 * full support is this rule's, with this rule's nodes. One path, no memo;
 * *visits, when not NULL, gets the nodes of e it visits, none for zero. The
 * translation family's is the scalar family's walk, bqd_sm_diag_walk, on full
 * support. From a Lace worker.
 */
BQD bqd_xpsm_diagonal(BQD e, uint64_t A, EVBDD_WGT phase, uint32_t nqubits, uint32_t *visits);

/* --- building, decoding, states -------------------------------------------- */

/** The canonical diagram of a vector of length 2^nqubits, bqd_from_vector's. */
TASK_DECL_2(BQD, bqd_xpsm_from_vector, const EVBDD_WGT *, uint32_t);

/** Decode into a vector of length 2^nqubits, bqd_to_vector's. */
VOID_TASK_DECL_3(bqd_xpsm_to_vector, BQD, uint32_t, EVBDD_WGT *);

/**
 * One amplitude, bqd_eval's: an S node's high cofactor is one path, and a Q
 * node's forks into both children where the low one is not zero.
 */
EVBDD_WGT bqd_xpsm_eval(BQD e, uint32_t nqubits, uint64_t x);

/** |x>, x a vector index: X^x on the chain of |0...0>, Q nodes with no high cofactor. */
BQD bqd_xpsm_basis_state(uint64_t x, uint32_t nqubits);

/**
 * phase^{x_A} in the Pauli family, A a vector-index mask: bqd_xp_monomial
 * with this rule's Compose. The translation family's is the scalar family's,
 * bqd_sm_monomial, whose labels are scalars.
 */
BQD bqd_xpsm_monomial(uint64_t A, EVBDD_WGT phase, uint32_t nqubits);

#ifdef __cplusplus
}
#endif

#endif
