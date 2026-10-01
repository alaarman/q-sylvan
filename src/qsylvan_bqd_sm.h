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
 * Rule SM, Shannon where misaligned, for the BQD's scalar family
 * (skip:sec:sm of the note "Binary Quotient Diagrams with Level Skipping",
 * bqd-skip.tex, whose labels are cited as skip:...).
 *
 * The copy rule of def:bqd stores at every node the low cofactor and the
 * ratio, f_1 = f_0 (.) r, with r a copy of f_1 where f_0 is zero. A node whose
 * low cofactor is zero on part of the high one's support is then a mixture, a
 * ratio on supp f_0 and the high cofactor off it, and an operation that
 * recurses on the stored edges has to multiply off a support to get one from
 * the other. Rule SM has no copy. A node carries a tag that says which of two
 * kinds it is (skip:def:sm), with A, R the targets of its two edges and c_1
 * the high edge's scalar:
 *
 *     Q (nested)    f_0 = [A],  f_1 = [A] . c_1 [R],  supp R inside supp A
 *     S (Shannon)   f_0 = [A],  f_1 = c_1 [R]
 *
 * S exactly where f_0 = 0 (A is then zero and c_1 = 1) or where supp f_1 is
 * not inside supp f_0. The product is the plain one, and R of a Q node is the
 * representative of f_1 / f_0 on supp f_0, zero off it. A level is skipped
 * exactly where f_0 = f_1, and a node below level v read at v is the virtual
 * Q node (N, Ind N), whose ratio is the indicator of its support. The support
 * of a node is {0} x supp [A] and {1} x supp [R] for both tags
 * (skip:lem:smshadow), so the indicator, the subset test and the dependence
 * test walk the stored edges and do no arithmetic. The diagram is canonical
 * (skip:thm:smcanon); on every function of full support and every function
 * whose support is an affine subspace it is the copy rule's, node for node
 * (skip:cor:smd0), and elsewhere a node the copy rule would mix is an S node.
 *
 * Storage. The tag and a fullness flag, set where the node's function has no
 * zero (skip:lem:smshadow (b)), are two of the node's flag bits
 * (qsylvan_limdd_node.h): bit 45 of its low word for S and bit 44 for full.
 * Both are functions of the node's function, so they are part of its identity
 * in the unique table without making a second node for it, and a collection,
 * which keeps every surviving node at its index, keeps them. A node made
 * without them would be a second node for its function, so every node of an
 * SM session is made by bqd_sm_mk, and the copy rule's constructors never run
 * in one: the public entry points ask bqd_sm() first.
 *
 * What is here: the builder, the decoder and the amplitude query
 * (skip:alg:smbuild); the constructors Compose, mk_nested and the Shannon
 * rebuild (skip:alg:smcons); the support walks Ind, Sub, Dep and Side; the
 * cofactors; Apply, the pointwise product, ratio and sum, with the two
 * shortcuts on nested nodes (skip:alg:smapply); the gate with its two
 * shortcuts, the one-sided gate and the ratio that passes through
 * (skip:alg:smgate); the selections Restrict, Project, XQ, Pair and Perm
 * (skip:alg:smsel); and the diagonal gates and phase multiplications, Mono,
 * Diag, its walk on full support, PhaseMul, Exp and MulOff (skip:alg:smdiag).
 * None of them multiplies off a support but MulOff, which nothing here calls.
 * The public entry points of qsylvan_bqd_gates.h and qsylvan_bqd_ops.h come
 * here first in an SM session.
 *
 * Every operation is a Lace task memoised in the operation cache on level-free
 * keys with the scalar taken out, as the copy rule's are, under operation ids
 * of its own (sylvan_int.h), and counts its misses for bqd_cache_fit, except
 * the walk of a diagonal gate on full support, which keeps the copy rule's
 * count of the nodes it visits and has no memo.
 */

#ifndef QSYLVAN_BQD_SM_H
#define QSYLVAN_BQD_SM_H

#include "qsylvan_bqd.h"
#include "qsylvan_bqd_exp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* --- nodes ---------------------------------------------------------------- */

/** The tag S and the fullness flag, as node flags: bits 45 and 44 of the low word. */
#define BQD_SM_FULL (UINT32_C(1) << 12)
#define BQD_SM_S    (UINT32_C(1) << 13)

/** Whether N is an S node; the terminal is not. */
static inline bool
bqd_sm_is_s(LIMDD_TARG N)
{
    return N != LIMDD_TERMINAL && (limdd_node_flags(N) & BQD_SM_S) != 0;
}

/** Whether [N] has no zero: the terminal does, 0 (the zero function) does not. O(1). */
static inline bool
bqd_sm_full(LIMDD_TARG N)
{
    return N == LIMDD_TERMINAL || (N != 0 && (limdd_node_flags(N) & BQD_SM_FULL) != 0);
}

/**
 * The node (var; lo, hi) with tag S when `s` and Q otherwise, and its fullness
 * flag, set exactly for a Q node whose two edges are nonzero and lead to full
 * nodes. The one maker of nodes in an SM session. lo is the identity on its
 * node or zero; a Q node has a nonzero lo, an S node a nonzero hi, and an S
 * node with a zero lo an identity hi. It decides nothing else: that the tag
 * is the right one, and that the level is not one to skip, is the caller's.
 */
LIMDD_TARG bqd_sm_mk(uint32_t var, bool s, BQD lo, BQD hi);

/* --- constructors (skip:alg:smcons) --------------------------------------- */

/**
 * The canonical edge at level var of the function with the cofactors lo and
 * hi, canonical edges read at var + 1: lo itself where the two are equal, an
 * S node where lo is zero or the support of hi is not inside lo's, and a Q
 * node with the ratio otherwise.
 */
TASK_DECL_3(BQD, bqd_sm_compose, uint32_t, BQD, BQD);
#define bqd_sm_compose(var, lo, hi) RUN(bqd_sm_compose, var, lo, hi)

/**
 * The canonical edge at level var of the function with low cofactor lo and
 * high cofactor lo . rho, where rho is a canonical edge that is zero off
 * supp lo: zero where lo is, lo itself where rho is the indicator of its
 * support, and a Q node otherwise. No quotient and no subset test.
 */
TASK_DECL_3(BQD, bqd_sm_mk_nested, uint32_t, BQD, BQD);
#define bqd_sm_mk_nested(var, lo, rho) RUN(bqd_sm_mk_nested, var, lo, rho)

/**
 * The S node at level var from the images lo and hi of the two cofactors of
 * an S node under a map that keeps supports, so that they stay misaligned:
 * the scale of lo, or of hi where lo is zero, comes out onto the edge.
 */
BQD bqd_sm_rebuild_s(uint32_t var, BQD lo, BQD hi);

/* --- the support shadow ---------------------------------------------------- */

/** The node of the indicator of supp [N]: 0 for 0, the terminal on full support. */
TASK_DECL_1(LIMDD_TARG, bqd_sm_ind, LIMDD_TARG);
#define bqd_sm_ind(N) RUN(bqd_sm_ind, N)

/** Whether supp [X] lies inside supp [Y], X and Y nodes or 0 for the zero function. */
TASK_DECL_2(int, bqd_sm_sub, LIMDD_TARG, LIMDD_TARG);
#define bqd_sm_sub(X, Y) RUN(bqd_sm_sub, X, Y)

/** Whether a node of variable q is reachable from N, a node or 0. */
TASK_DECL_2(int, bqd_sm_dep, LIMDD_TARG, uint32_t);
#define bqd_sm_dep(N, q) RUN(bqd_sm_dep, N, q)

/* --- cofactors ------------------------------------------------------------- */

/**
 * The low cofactor and the ratio of [N] read at level var, N not an S node
 * there: its own two edges where N sits at var, and the virtual node
 * (N, Ind N) where N is below var or is the terminal.
 */
VOID_TASK_DECL_4(bqd_sm_nest, LIMDD_TARG, uint32_t, BQD *, BQD *);

/** The high cofactor of node N, as an edge read below it: the stored edge of an S node. */
TASK_DECL_1(BQD, bqd_sm_cof1, LIMDD_TARG);

/** The cofactor x_var = b of e read at var, at or above its node, canonical. */
TASK_DECL_3(BQD, bqd_sm_cof, BQD, uint32_t, int);
#define bqd_sm_cof(e, var, b) RUN(bqd_sm_cof, e, var, b)

/* --- pointwise operations (skip:alg:smapply) ------------------------------- */

/**
 * The three operations of Apply. RATIO(f, g) is g / f on supp f and 0 off it,
 * and is defined only where supp g lies inside supp f, which the caller
 * guarantees.
 */
enum { BQD_SM_MUL = 0, BQD_SM_RATIO = 1, BQD_SM_ADD = 2 };

/** op(f, g) for two canonical edges read at one level; the result is canonical. */
TASK_DECL_3(BQD, bqd_sm_apply, int, BQD, BQD);
#define bqd_sm_apply(op, f, g) RUN(bqd_sm_apply, op, f, g)

/** rho . [Y != 0], for Y nonzero: rho where its support is inside Y's already. */
TASK_DECL_2(BQD, bqd_sm_restr, BQD, BQD);
#define bqd_sm_restr(rho, Y) RUN(bqd_sm_restr, rho, Y)

/* --- selections (skip:alg:smsel) ------------------------------------------- */

/*
 * Each takes every value of its result from one value of its operands, so a
 * Q node above the qubits it acts on maps to mk_nested of the images of its
 * two stored edges (skip:lem:smpair (d)), with no high cofactor and no
 * quotient. An S node maps to Compose of the images of its cofactors where
 * the result can come out nested or redundant (Restrict, Project, Pair), and
 * keeps its tag under a bijection (XQ, Perm). The signatures are those of
 * bqd_restrict, bqd_project, bqd_x, bqd_pair and bqd_perm.
 */

/** e|_{x_q = b}, as a function of the same variables that skips q. */
TASK_DECL_3(BQD, bqd_sm_restrict, BQD, uint32_t, int);
#define bqd_sm_restrict(e, q, b) RUN(bqd_sm_restrict, e, q, b)

/** e . [x_q = b]. */
TASK_DECL_3(BQD, bqd_sm_project, BQD, uint32_t, int);
#define bqd_sm_project(e, q, b) RUN(bqd_sm_project, e, q, b)

/** X on qubit q: the two cofactors at q trade places. */
TASK_DECL_2(BQD, bqd_sm_xq, BQD, uint32_t);
#define bqd_sm_xq(e, q) RUN(bqd_sm_xq, e, q)

/** [x_b = 0] A|_{x_b = s} + [x_b = 1] B|_{x_b = t}. */
TASK_DECL_5(BQD, bqd_sm_pair, BQD, int, BQD, int, uint32_t);
#define bqd_sm_pair(A, s, B, t, b) RUN(bqd_sm_pair, A, s, B, t, b)

/** The permutation `kind` (BQD_PERM_*, qsylvan_bqd_gates.h) of the qubits qa < qb. */
TASK_DECL_4(BQD, bqd_sm_perm, BQD, uint32_t, uint32_t, uint32_t);
#define bqd_sm_perm(e, kind, qa, qb) RUN(bqd_sm_perm, e, kind, qa, qb)

/* --- gates (skip:alg:smgate) ------------------------------------------------ */

/**
 * Side: 0 or 1 where supp [N], N a node or the terminal read at q or above,
 * lies in {x_q = 0} or in {x_q = 1}, and 2 where it lies in neither. The
 * shadow again, with no arithmetic; 2 at once on full support.
 */
TASK_DECL_2(int, bqd_sm_side, LIMDD_TARG, uint32_t);
#define bqd_sm_side(N, q) RUN(bqd_sm_side, N, q)

/**
 * gates[gid] on qubit q, conditioned on the controls in cmask (bit c for
 * qubit c, every one above q) that are still pending: Gate of
 * skip:alg:smgate, bqd_cgate's recursion under SM. Above every pending qubit
 * it has two shortcuts before the recursion on both cofactors: where no
 * control is pending and the support lies on one side of x_q, the gate is a
 * restriction and a phase (skip:lem:smpass); and where a Q node's ratio does
 * not depend on x_q, the gate passes through it and acts on the low cofactor
 * alone.
 */
TASK_DECL_4(BQD, bqd_sm_gate, BQD, uint32_t, uint64_t, uint32_t);
#define bqd_sm_gate(e, gid, cmask, q) RUN(bqd_sm_gate, e, gid, cmask, q)

/* --- diagonal gates and phases (skip:alg:smdiag) ---------------------------- */

/*
 * The masks of a diagonal gate here are variable masks, bit v for qubit v,
 * where the public entry points take vector-index masks, bit n-1-v for qubit
 * v; bqd_sm_var_mask turns one into the other.
 */

static inline uint64_t
bqd_sm_var_mask(uint64_t A, uint32_t nqubits)
{
    uint64_t vars = 0;
    for (uint64_t rest = A; rest != 0; rest &= rest - 1)
        vars |= UINT64_C(1) << (nqubits - 1 - (uint32_t)__builtin_ctzll(rest));
    return vars;
}

/** Mono: phase^{prod of x_v, v in A}, one Q node per variable of A, the terminal for phase 1. */
BQD bqd_sm_monomial(uint64_t A, EVBDD_WGT phase);

/** Diag: e . phase^{prod of x_v, v in A} on any support, memoised. */
TASK_DECL_3(BQD, bqd_sm_diag, BQD, uint64_t, EVBDD_WGT);
#define bqd_sm_diag(e, A, phase) RUN(bqd_sm_diag, e, A, phase)

/**
 * The same on full support by the walk of skip:prop:diag, bqd_apply_diagonal's
 * on the copy rule: one path, no memo, and *visits counts the nodes of e on
 * it. e must have full support.
 */
BQD bqd_sm_diag_walk(BQD e, uint64_t A, EVBDD_WGT phase, uint32_t *visits);

/** PhaseMul: k . beta^eps on any support, r the order of beta (bqd_exp_order). */
TASK_DECL_4(BQD, bqd_sm_phase_mul, BQD, BQD_EXP, EVBDD_WGT, uint32_t);

/** Exp: beta^eps, of full support. */
TASK_DECL_3(BQD, bqd_sm_exp, BQD_EXP, EVBDD_WGT, uint32_t);

/** iota: the indicator of supp [N] as a 0/1 exponent; N a node, the terminal, or 0. */
TASK_DECL_1(BQD_EXP, bqd_sm_iota, LIMDD_TARG);

/** MulOff: k on supp [u] and c . k off it; u a node, the terminal, or 0. */
TASK_DECL_3(BQD, bqd_sm_mul_off, BQD, EVBDD_WGT, LIMDD_TARG);

/* --- building, decoding, states -------------------------------------------- */

/** The canonical SM diagram of a vector of length 2^nqubits, bqd_from_vector's. */
TASK_DECL_2(BQD, bqd_sm_from_vector, const EVBDD_WGT *, uint32_t);

/** Decode into a vector of length 2^nqubits, bqd_to_vector's. */
VOID_TASK_DECL_3(bqd_sm_to_vector, BQD, uint32_t, EVBDD_WGT *);

/**
 * One amplitude, bqd_eval's. An S node's high cofactor is one path; a Q node's
 * forks into both children, as under the copy rule, and stops early where the
 * low one is zero.
 */
EVBDD_WGT bqd_sm_eval(BQD e, uint32_t nqubits, uint64_t x);

/** |x>, x a vector index (qubit q at bit n-1-q): an S node at a 1 and a Q node at a 0. */
BQD bqd_sm_basis_state(uint64_t x, uint32_t nqubits);

#ifdef __cplusplus
}
#endif

#endif
