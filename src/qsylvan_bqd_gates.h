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
 * two new cofactors, again with the copy. Every step is itself such an
 * operation one level down, so the recursion ends at the terminal. A gate on
 * qubit q is the same walk down to level q, where the two new cofactors are
 * the gate's linear combination of the old ones: for the Hadamard,
 * (f_0 + f_1)/sqrt2 and (f_0 - f_1)/sqrt2. Results are canonical, the same
 * edge bqd_from_vector builds for the same vector.
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
 * Everything is for the scalar family, and checked; from a Lace worker.
 */

#ifndef QSYLVAN_BQD_GATES_H
#define QSYLVAN_BQD_GATES_H

#include "qsylvan_bqd.h"
#include "qsylvan_bqd_ops.h"

#ifdef __cplusplus
extern "C" {
#endif

/* --- states --------------------------------------------------------------- */

/** The basis state |x>, x a vector index (qubit q at bit n-1-q): n nodes. */
BQD bqd_basis_state(uint64_t x, uint32_t nqubits);

/** Whether e is nonzero at every point. Memoised on nodes, so O(|e|) at most. */
bool bqd_has_full_support(BQD e);

/* --- pointwise operations on any support ---------------------------------- */

/** f . g, on any supports. bqd_product is the O(|f||g|) case of full support. */
TASK_DECL_2(BQD, bqd_multiply, BQD, BQD);
#define bqd_multiply(f, g) RUN(bqd_multiply, f, g)

/** f + g. */
TASK_DECL_2(BQD, bqd_add, BQD, BQD);
#define bqd_add(f, g) RUN(bqd_add, f, g)

/**
 * The canonical edge, at level `var`, of the function whose cofactors are
 * f_0 = lo and f_1 = hi (edges at level var + 1).
 */
TASK_DECL_3(BQD, bqd_compose, uint32_t, BQD, BQD);
#define bqd_compose(var, lo, hi) RUN(bqd_compose, var, lo, hi)

/** The cofactor of e at its top variable = b, as an edge one level down. */
TASK_DECL_2(BQD, bqd_cofactor, BQD, int);
#define bqd_cofactor(e, b) RUN(bqd_cofactor, e, b)

/** c . e, and -e. */
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
 * which go through bqd_apply_diagonal. As limdd_gate.
 */
TASK_DECL_4(BQD, bqd_gate, BQD, uint32_t, uint32_t, uint32_t);
#define bqd_gate(e, gateid, target, nqubits) RUN(bqd_gate, e, gateid, target, nqubits)

/**
 * `gateid` on `target`, conditioned on every qubit in `control_mask` (bit c
 * for qubit c, as limdd_cgate) being |1>. A diagonal gate whose u00 is one is
 * a monomial and may have its controls anywhere. Any other gate needs every
 * control above the target, control_mask < 2^target, as limdd_cgate does, and
 * the call exits with a message otherwise.
 */
TASK_DECL_5(BQD, bqd_cgate, BQD, uint32_t, uint64_t, uint32_t, uint32_t);
#define bqd_cgate(e, gateid, control_mask, target, nqubits) \
    RUN(bqd_cgate, e, gateid, control_mask, target, nqubits)

/**
 * One control on either side of the target, as limdd_cgate_either: a diagonal
 * gate directly, and CX with the control below the target by
 * (H (x) H) CNOT_{a->b} (H (x) H) = CNOT_{b->a}. Any other gate with the
 * control below the target sets *ok to false and returns e.
 */
BQD bqd_cgate_either(BQD e, uint32_t gateid, uint32_t control, uint32_t target,
                     uint32_t nqubits, bool *ok);

/** Exchange two qubits, as three CNOTs, as limdd_swap. */
BQD bqd_swap(BQD e, uint32_t a, uint32_t b, uint32_t nqubits);

#ifdef __cplusplus
}
#endif

#endif
