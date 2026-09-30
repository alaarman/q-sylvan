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
 * Gates on a BQD of the scalar family, on any function: the Hadamard, every
 * other single-qubit gate, and controlled gates, with the signatures of
 * limdd_gate and limdd_cgate so that the two diagrams are called the same way.
 *
 * What the paper proves and what it does not. Its table of operations
 * (tab:ops) has the Hadamard on any function only as a dense fallback, the
 * vector rebuilt, O(n^2 2^n); an algorithm bounded in the sizes of the
 * diagrams is open, and sec:ops:hadamard shows that none is polynomial in its
 * input in general: a level-four phase state with O(n^4) nodes in every order
 * has a Hadamard image of 2^{n - O(sqrt n)} nodes in every order. The
 * pointwise product on arbitrary supports is open in the same sense.
 *
 * So what is here is correct on every function and carries no bound. It is a
 * recursion on the diagram, never on the vector. A node stores the low
 * cofactor f_0 and the ratio r, and the high cofactor is f_1 = f_0 (.) r, the
 * product with the copy where f_0 is zero (def:bqd). Any pointwise operation
 * commutes with taking cofactors, so it is computed on the two cofactors and
 * the result is put back into quotient form: its ratio is the quotient of the
 * two new cofactors, again with the copy, or, where the two are equal, no
 * node at all, since the result does not depend on that variable. Every step
 * is itself such an operation on the cofactors at the higher of its
 * operands' levels, where an operand that skips the level is both of its
 * cofactors, so the recursion ends at the terminal. A gate on qubit q is the
 * same walk down to level q, where the two new cofactors are the gate's
 * linear combination of the old ones: for the Hadamard, (f_0 + f_1)/sqrt2
 * and (f_0 - f_1)/sqrt2. Results are canonical, the same edge
 * bqd_from_vector builds for the same vector, levels skipped where it skips
 * them.
 *
 * The cost is the size of what the recursion meets, the diagrams of the
 * cofactors and of the intermediate sums, memoised on node pairs. Where those
 * stay small, as on the states of structured circuits, it goes far beyond the
 * vector; where they cannot, as on the family above, nothing does.
 *
 * Diagonal gates keep the O(n) walk of prop:diag while the state has full
 * support, and fall back to the product here when it does not, since the walk
 * is only correct on full support.
 *
 * Everything here is for all three label families. In the translation and
 * Pauli families the entry points dispatch to qsylvan_bqd_xp.h, whose
 * recursion is the one above on labelled edges, at the top of its arguments
 * rather than at their nodes, with Canon wherever the scalar family returns
 * an operand or a cofactor as it is. The diagonal gates are the O(n) walk on
 * full support in every family: the scalar family's in the translation
 * family too, whose diagram of a function of full support is the scalar
 * family's, and one with label products and a sign repair in the Pauli
 * family (qsylvan_bqd_xp.h). The operations declared with a RUN macro, and
 * bqd_cgate_either and bqd_swap, which reach the diagram only through such
 * macros, may be called from any thread, as limdd_gate may;
 * bqd_local_matvec, which pushes on the worker's reference stack itself, only
 * from a Lace worker.
 */

#ifndef QSYLVAN_BQD_GATES_H
#define QSYLVAN_BQD_GATES_H

#include "qsylvan_bqd.h"
#include "qsylvan_bqd_exp.h"
#include "qsylvan_bqd_ops.h"

#ifdef __cplusplus
extern "C" {
#endif

/* --- states --------------------------------------------------------------- */

/**
 * The basis state |x>, x a vector index (qubit q at bit n-1-q): n nodes. In
 * the translation and Pauli families X^x on the nodes of |0...0>, which every
 * basis state shares.
 */
BQD bqd_basis_state(uint64_t x, uint32_t nqubits);

/** Whether e is nonzero at every point. Memoised on nodes, so O(|e|) at most. */
bool bqd_has_full_support(BQD e);

/* --- pointwise operations on any support ---------------------------------- */

/**
 * f . g, on any supports, in every family. bqd_product is the O(|f||g|) case
 * of full support in the scalar family.
 */
TASK_DECL_2(BQD, bqd_multiply, BQD, BQD);
#define bqd_multiply(f, g) RUN(bqd_multiply, f, g)

/** f + g, in every family. */
TASK_DECL_2(BQD, bqd_add, BQD, BQD);
#define bqd_add(f, g) RUN(bqd_add, f, g)

/**
 * The canonical edge, at level `var`, of the function whose cofactors are
 * f_0 = lo and f_1 = hi (canonical edges read at level var + 1, which may
 * skip levels below it). Where lo == hi the function does not depend on
 * x_var, and the result is lo itself, which skips var (skip:alg:constructors).
 * In the Pauli family so is lo == -hi, with a Z at var on lo's label
 * (skip:prop:xcompose).
 */
TASK_DECL_3(BQD, bqd_compose, uint32_t, BQD, BQD);
#define bqd_compose(var, lo, hi) RUN(bqd_compose, var, lo, hi)

/**
 * The cofactor x_var = b of the function e denotes read at level `var`, as an
 * edge read at var + 1, canonical. An edge that skips var, whose node is below
 * it or is the terminal, denotes a function that does not depend on x_var,
 * and is its own cofactor on both sides; in the Pauli family unless its label
 * has a Z at var, and the two cofactors are then opposite.
 */
TASK_DECL_3(BQD, bqd_cofactor, BQD, uint32_t, int);
#define bqd_cofactor(e, var, b) RUN(bqd_cofactor, e, var, b)

/** c . e, and -e, in every family. */
BQD bqd_scale(BQD e, EVBDD_WGT c);
BQD bqd_negate(BQD e);

/**
 * e|_{x_q = b}, the cofactor as a function of the same variables that does
 * not depend on x_q, and e . [x_q = b], the other half set to zero. As
 * evbdd_restrict and evbdd_project, and limdd_restrict and limdd_project.
 */
TASK_DECL_3(BQD, bqd_restrict, BQD, uint32_t, int);
#define bqd_restrict(e, q, b) RUN(bqd_restrict, e, q, b)
TASK_DECL_3(BQD, bqd_project, BQD, uint32_t, int);
#define bqd_project(e, q, b) RUN(bqd_project, e, q, b)

/* --- gates without the high cofactor (skip:sec:ratio) --------------------- */

/*
 * A node stores its low cofactor and the ratio of its two cofactors, so an
 * operation that recurses on cofactors makes the high one (Cof1, a product)
 * at every node above the level it acts on, and forms each new ratio by a
 * quotient. The operations below read a node's stored edges instead wherever
 * the gate leaves the pairing of the two cofactors alone: an operation that
 * takes each value of its result from one value of an operand commutes with
 * the pair (low cofactor, ratio) (skip:lem:select). Above the qubits they act
 * on they make no Cof1 and no Apply, and at most one call per node of the
 * diagram. The gate entry points below use them, and so do bqd_restrict,
 * bqd_project and, on a state without full support, bqd_apply_diagonal.
 * Every result is the canonical edge, the one bqd_from_vector builds.
 */

/** The three permutations of two qubits qa < qb, (pi f)(x) = f(pi x). */
enum {
    BQD_PERM_SWAP    = 0,   /* exchange x_qa and x_qb */
    BQD_PERM_CX_DOWN = 1,   /* CX, control qa above target qb: x_qb ^= x_qa */
    BQD_PERM_CX_UP   = 2,   /* CX, control qb below target qa: x_qa ^= x_qb */
};

/**
 * Perm of skip:alg:perm and skip:alg:xratio: the permutation `kind` of the
 * qubits qa < qb. Above qa it maps a node through its two stored edges; at qa
 * it takes the node's two cofactors, and the new ones are Pairs at qb.
 */
TASK_DECL_4(BQD, bqd_perm, BQD, uint32_t, uint32_t, uint32_t);
#define bqd_perm(e, kind, qa, qb) RUN(bqd_perm, e, kind, qa, qb)

/** X on qubit q: XQ of skip:alg:perm, and Canon(X_q . e) in the other families (CanonT). */
TASK_DECL_2(BQD, bqd_x, BQD, uint32_t);
#define bqd_x(e, q) RUN(bqd_x, e, q)

/**
 * Pair of skip:alg:perm: [x_b = 0] A|_{x_b = s} + [x_b = 1] B|_{x_b = t}, the
 * function whose values come from A where x_b is 0 and from B where it is 1.
 */
TASK_DECL_5(BQD, bqd_pair, BQD, int, BQD, int, uint32_t);
#define bqd_pair(A, s, B, t, b) RUN(bqd_pair, A, s, B, t, b)

/** Ind: the node of the indicator of supp [t], the terminal on full support. */
TASK_DECL_1(LIMDD_TARG, bqd_ind, LIMDD_TARG);
#define bqd_ind(t) RUN(bqd_ind, t)

/**
 * MulOff of skip:alg:rebuild: [k] . c^{[x not in supp [u]]}, k on the support
 * of u and c k off it; in the scalar and the translation family.
 */
TASK_DECL_3(BQD, bqd_mul_off, BQD, EVBDD_WGT, LIMDD_TARG);
#define bqd_mul_off(k, c, u) RUN(bqd_mul_off, k, c, u)

/**
 * PhaseMul of skip:alg:phasemul and skip:alg:xphase: [k] . beta^eps on any
 * support, eps an exponent diagram (qsylvan_bqd_exp.h) with its values in
 * the base beta, which in the Pauli family must be a power of w_8. No Cof1
 * and no Apply. eps is protected by the caller.
 */
TASK_DECL_3(BQD, bqd_phase_mul, BQD, BQD_EXP, EVBDD_WGT);
#define bqd_phase_mul(k, eps, beta) RUN(bqd_phase_mul, k, eps, beta)

/** Exp of skip:alg:phasemul: beta^eps, of full support, r the order of beta; not in the Pauli family. */
TASK_DECL_3(BQD, bqd_exp_state, BQD_EXP, EVBDD_WGT, uint32_t);
#define bqd_exp_state(eps, beta, r) RUN(bqd_exp_state, eps, beta, r)

/**
 * DiagR: e . phase^{x_A} on any support by PhaseMul, A a vector-index mask as
 * bqd_apply_diagonal's, which calls it off full support. In the Pauli family a
 * phase that is not a power of w_8 takes the product with the monomial.
 */
TASK_DECL_4(BQD, bqd_diag_any, BQD, uint64_t, EVBDD_WGT, uint32_t);
#define bqd_diag_any(e, A, phase, nqubits) RUN(bqd_diag_any, e, A, phase, nqubits)

/**
 * The dense 2^k x 2^k matrix M applied to the qubits of e: row-major, bit
 * k-1-i of an index is qubits[i], so gates[] order for k = 1. It is
 * sum_r [x_Q = r] . sum_c M[r][c] . e|_{x_Q = c}, from the operations above,
 * memoised on nodes. For a single qubit it computes what bqd_gate does, by a
 * longer route: bqd_gate combines the two cofactors at the qubit's level and
 * leaves the levels below alone, where this builds each restriction and
 * projection over the whole width first.
 */
BQD bqd_local_matvec(BQD e, const EVBDD_WGT *M, const uint32_t *qubits, uint32_t k,
                     uint32_t nqubits);

/* --- gates, called as the LIMDD's ----------------------------------------- */

/**
 * `gateid` (an index into gates[]) applied to qubit `target`. Any 2x2 matrix
 * in the table: the Hadamard, X, Y, the square roots, and the diagonal gates,
 * which go through bqd_apply_diagonal; X, Y and any other matrix with a zero
 * diagonal go through bqd_x after the diagonal gate that makes them X. As
 * limdd_gate.
 */
TASK_DECL_4(BQD, bqd_gate, BQD, uint32_t, uint32_t, uint32_t);
#define bqd_gate(e, gateid, target, nqubits) RUN(bqd_gate, e, gateid, target, nqubits)

/**
 * `gateid` on `target`, conditioned on every qubit in `control_mask` (bit c
 * for qubit c, as limdd_cgate) being |1>. A diagonal gate whose u00 is one is
 * a monomial and may have its controls anywhere. Any other gate needs every
 * control above the target, control_mask < 2^target, as limdd_cgate does, and
 * the call exits with a message otherwise. A gate with a zero diagonal and one
 * control is bqd_perm between two monomials, the CX itself Perm alone.
 */
TASK_DECL_5(BQD, bqd_cgate, BQD, uint32_t, uint64_t, uint32_t, uint32_t);
#define bqd_cgate(e, gateid, control_mask, target, nqubits) \
    RUN(bqd_cgate, e, gateid, control_mask, target, nqubits)

/**
 * One control on either side of the target, as limdd_cgate_either: a diagonal
 * gate directly, and CX, CY and any gate with a zero diagonal by bqd_perm with
 * the control below the target, BQD_PERM_CX_UP, between two monomials. Any
 * other gate with the control below the target sets *ok to false and returns e.
 */
BQD bqd_cgate_either(BQD e, uint32_t gateid, uint32_t control, uint32_t target,
                     uint32_t nqubits, bool *ok);

/** Exchange two qubits, bqd_perm with BQD_PERM_SWAP, as limdd_swap. */
BQD bqd_swap(BQD e, uint32_t a, uint32_t b, uint32_t nqubits);

#ifdef __cplusplus
}
#endif

#endif
