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
 * Binary quotient diagrams, on the LIMDD's tables.
 *
 * A Shannon node stores the two cofactors of the function it denotes. A
 * quotient node stores the low cofactor and the pointwise RATIO of the two,
 *
 *     f_0 = l_0 . g_0,        f_1 = l_1 . (g_0 (.) g_1),
 *
 * where g_0, g_1 are the children, l_0, l_1 the edge labels, and (.) is the
 * pointwise product extended by a copy where g_0 is zero: there is nothing to
 * divide by there, and copying the high value is the one total choice
 * (Laarman, "Binary Quotient Diagrams", def:bqd). On a stabiliser state the
 * ratio of two cofactors is a derivative of the phase polynomial, which lowers
 * its degree, so the recursion stops after k steps at level k of the Clifford
 * hierarchy. That is where the succinctness comes from.
 *
 * Three label families, each a diagram of its own:
 *
 *     BQD_FAMILY_SCALAR   labels c            the BQD
 *     BQD_FAMILY_X        labels c X^t        the X-BQD, a scalar and a
 *                                             translation of the domain
 *     BQD_FAMILY_PAULI    labels c Z^s X^t    the Pauli-BQD
 *
 * A label acts on a function of the child's variables by
 *
 *     (c Z^s X^t . g)(y) = c (-1)^{s.y} g(y xor t).
 *
 * That is exactly a LIM. So a BQD is a different READING of the LIMDD's
 * tables, not different storage: labels are LIMs in the LIM table, a node is
 * (var, low target, high label, high target) interned by limdd_makenode with
 * the LIMDD's packing, an edge is a LIMDD edge, limdd_protect keeps a root
 * alive across a collection, and limdd_width, limdd_level_counts and
 * limdd_fprintdot in qsylvan_limdd_inspect.h read a BQD unchanged. What this
 * module adds is the quotient semantics: the canonical form, the decoder, and
 * the amplitude query. The one thing to never do is hand a BQD edge to a
 * LIMDD operation, which would read the same triple under the Shannon rule.
 *
 * Canonical form. Every stored node is the REPRESENTATIVE of the function it
 * denotes: for the X and Pauli families translated so that the least point
 * of its support is 0; scaled so that the value at the least point is 1;
 * and, for the Pauli family, multiplied by the sign pattern that gives the
 * value at every pivot an argument in [0, pi) (def:rep, def:prep). The scalar
 * family has no translation, so there only the scale is divided out and the
 * least point stays where it is. The representative is a function of the
 * vector, so the reduced diagram in a fixed order is unique and equality is a
 * comparison of two edges. Building a node is a table lookup, with no search over the label
 * orbit -- the point of the design, and the difference from the LIMDD, whose
 * canonical form is an orbit computation. The price is that a function and
 * its translate can be two nodes, since the least point of a translated
 * support is not the translate of the least point: the diagram is reduced in
 * the sense that no two nodes denote one function, and not in the LIMDD's
 * sense that no two nodes are equal up to a label. Storing representatives
 * also means the low edge of every node carries the identity label, or is the
 * zero edge where the low cofactor is zero, which happens in the scalar
 * family only: in the other two, 0 is in the support of every stored node.
 * That is what lets a node pack exactly like a LIMDD node, whose low edge is
 * the identity or zero too. It is a consequence of def:rep, not a theorem.
 *
 * Conventions, fixed once. Qubit 0 is the top variable and is decided first.
 * A vector of length 2^n is indexed by x with qubit q at bit (n-1-q), so the
 * top qubit is the most significant bit and the first half of a vector is the
 * cofactor at x_top = 0. Translations and sign patterns are masks in that
 * same convention. The terminal is reached at level n only: levels are never
 * skipped and a redundant node is kept, which is the paper's convention and
 * NOT the fully reduced LIMDD's. One difference from the paper's counts
 * remains: the paper counts the function 0 as a node at every level where it
 * occurs (the +2 of thm:pcoset is the code-state node and 0), and here zero
 * is an edge, never a node. limdd_level_counts is therefore the paper's count
 * minus one at every level with a zero cofactor, and a comparison with the
 * paper's tables has to add it back.
 *
 * Scope. Here a diagram is built from an amplitude vector, exponential in n,
 * which is what the paper's evaluation does. qsylvan_bqd_ops.h has the two
 * operations the paper proves, and qsylvan_bqd_gates.h the states and gates
 * of a circuit simulator, the Hadamard among them, correct on every function
 * and without the bound the paper leaves open. *
 * Parallelism follows the rest of the LIMDD code. The builder and the decoder
 * are Lace tasks that spawn the low cofactor and compute the high one, and
 * every representative is computed in per-call storage; the only shared state
 * is the lock-free tables. Weight arithmetic goes through the operation
 * cache, so everything here runs on a Lace worker.
 */

#ifndef QSYLVAN_BQD_H
#define QSYLVAN_BQD_H

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#include <sylvan_int.h>
#include "sylvan_edge_weights.h"
#include "qsylvan_limdd_lim.h"
#include "qsylvan_limdd_node.h"

#ifdef __cplusplus
extern "C" {
#endif

/** A BQD edge is a LIMDD edge: the label's LIM index above the target. */
typedef LIMDD BQD;

typedef enum {
    BQD_FAMILY_SCALAR = 0,
    BQD_FAMILY_X      = 1,
    BQD_FAMILY_PAULI  = 2,
} bqd_family_t;

/**
 * Create the tables through limdd_nodes_init, and fix the label family. A
 * vector is indexed by a 64-bit word, so `nqubits` may not exceed 63. Must be
 * called from a Lace worker; the edge weight table must already exist and
 * must hold complex weights. Not to be combined with limdd_nodes_init in one
 * session: the two share the tables.
 */
void bqd_init(bqd_family_t family, size_t nqubits, size_t node_tablesize,
              size_t pauli_tablesize, size_t lim_tablesize, size_t stab_tablesize);
void bqd_quit(void);

bqd_family_t bqd_family(void);
const char  *bqd_family_name(bqd_family_t f);

/* --- labels --------------------------------------------------------------- */

/**
 * The label c Z^s X^t, with s and t as vector-index masks of an `nqubits`
 * vector. For the scalar family s and t must be 0; for the X family s must
 * be 0. Interned in the LIM table; the zero weight gives the zero label.
 */
LIMDD_LIM bqd_lim_make(EVBDD_WGT c, uint64_t s, uint64_t t, uint32_t nqubits);

/** Read a label back as (c, s, t). The zero label reads as (0, 0, 0). */
void bqd_lim_masks(LIMDD_LIM lim, uint32_t nqubits, EVBDD_WGT *c, uint64_t *s, uint64_t *t);

/* --- the canonical form --------------------------------------------------- */

/**
 * The representative of a vector of length `len` in the current family, and
 * how to get the vector back:  g = (c0 (-1)^{s.p}  Z^s X^p) . rep.
 *
 * `rep` may alias `g`. For the scalar family p and s are 0; for the X family
 * s is 0. A zero vector is its own representative with c0 = 1. Per-call
 * storage only.
 */
void bqd_representative(const EVBDD_WGT *g, uint64_t len, EVBDD_WGT *rep,
                        EVBDD_WGT *c0, uint64_t *s, uint64_t *p);

/**
 * The canonical diagram of a vector of length 2^nqubits, with EVBDD_ZERO for
 * a zero amplitude. Exponential in nqubits by nature. A Lace task: the low
 * cofactor is spawned, the high one computed, and the two joined.
 */
TASK_DECL_2(BQD, bqd_from_vector, const EVBDD_WGT *, uint32_t);
#define bqd_from_vector(f, n) RUN(bqd_from_vector, f, n)

/** Decode a diagram back into a vector of length 2^nqubits. A Lace task. */
VOID_TASK_DECL_3(bqd_to_vector, BQD, uint32_t, EVBDD_WGT *);
#define bqd_to_vector(e, n, out) RUN(bqd_to_vector, e, n, out)

/**
 * One amplitude. Under the quotient rule this is NOT one path: a high edge
 * needs both children at the same point, so the walk forks at every high
 * edge it takes and costs 2^{|x|} visits (sec:ops:queries). From a worker.
 */
EVBDD_WGT bqd_eval(BQD e, uint32_t nqubits, uint64_t x);

/* Node access, width, level counts and dot output are the LIMDD's:
 * limdd_node_var/low/high, limdd_width, limdd_level_counts, limdd_nodecount,
 * limdd_fprintdot. */

#ifdef __cplusplus
}
#endif

#endif
