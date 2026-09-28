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
 * Low-level operations on EVDD (QMDD) state vectors, with the same names and
 * meaning as the LIMDD's (qsylvan_limdd_ops.h) and the BQD's
 * (qsylvan_bqd_gates.h), so the three diagrams can be driven the same way:
 *
 *     times      pointwise product f.g
 *     plus       pointwise sum f + g (evbdd_plus, in sylvan_evbdd.h)
 *     scale      c.f, and negate, -f
 *     restrict   f|_{x_q = b}, the cofactor, as a function of the same
 *                variables that does not depend on x_q
 *     project    f . [x_q = b], the other half set to zero: measurement
 *                before renormalising
 *     local_matvec  M applied to the qubits Q of f, for a dense 2^k x 2^k M
 *
 * A vector variable is a qubit, qubit q at variable q, and the diagram may
 * skip a variable it does not depend on.
 *
 * local_matvec is the relational product (exists X: M(X', X) /\ V(X, Y))
 * [X' := X] with the sum for the quantifier, done without a diagram for M:
 *
 *     M v  =  sum_r [x_Q = r] . sum_c M[r][c] . v|_{x_Q = c}
 *
 * so it is 2^k restrictions, 4^k scaled sums and 2^k projections, all
 * memoised on nodes by the operations above; nothing is keyed on M. M is
 * row-major, and bit k-1-i of a row or column index is qubits[i], so for
 * k = 1 it is the order of gates[] (u00, u01, u10, u11). evbdd_matvec_mult
 * is the other route, for a matrix that is a diagram of its own.
 */

#ifndef QSYLVAN_EVDD_OPS_H
#define QSYLVAN_EVDD_OPS_H

#include <sylvan_int.h>
#include "sylvan_evbdd.h"

#ifdef __cplusplus
extern "C" {
#endif

TASK_DECL_2(EVBDD, evbdd_times, EVBDD, EVBDD);
#define evbdd_times(a, b) RUN(evbdd_times, a, b)

EVBDD evbdd_scale(EVBDD e, EVBDD_WGT c);
EVBDD evbdd_negate(EVBDD e);

TASK_DECL_3(EVBDD, evbdd_restrict, EVBDD, BDDVAR, int);
#define evbdd_restrict(e, q, b) RUN(evbdd_restrict, e, q, b)

TASK_DECL_3(EVBDD, evbdd_project, EVBDD, BDDVAR, int);
#define evbdd_project(e, q, b) RUN(evbdd_project, e, q, b)

/** From a Lace worker. `nqubits` is the vector's width. */
EVBDD evbdd_local_matvec(EVBDD e, const EVBDD_WGT *M, const uint32_t *qubits,
                         uint32_t k, uint32_t nqubits);

#ifdef __cplusplus
}
#endif

#endif
