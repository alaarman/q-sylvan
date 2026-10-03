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
 * That is the copy rule. A session runs rule SM by default, which has no
 * copy: a node whose cofactors have misaligned supports stores both, as a
 * Shannon node, and says so in a tag (bqd_zero_rule_t; qsylvan_bqd_sm.h in
 * the scalar family, qsylvan_bqd_xp_sm.h in the translation and Pauli
 * families), and the copy rule is there on request. The two rules are one
 * diagram on full support. The rest of this header describes the copy rule.
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
 * Levels are skipped. A node whose function does not depend on its variable
 * is redundant and is not stored: an edge may point at a node any number of
 * levels below, and denotes that node's function extended to the levels it
 * skips (skip:def:skip, in the note "Binary Quotient Diagrams with Level
 * Skipping", bqd-skip.tex, whose labels are cited here as skip:...). That is
 * the LIMDD's convention and NOT the paper's, which keeps every level, but
 * the redundant node is not the LIMDD's. Under the quotient rule two equal
 * edges do not make a node redundant: the redundant node has the identity to
 * a node M as its low edge and the identity to ind M, the indicator of M's
 * support, as its high edge (skip:thm:redundant), which on full support is
 * the constant one (skip:cor:full). The builder skips where the
 * representative has two equal cofactors, and an operation where the two
 * cofactors of its result are equal, which for two canonical edges is one
 * comparison. A label on an edge acts on the levels it skips too: it never
 * translates them, and in the Pauli family it may sign them (skip:lem:normal).
 * The diagram is the paper's with every redundant node contracted into its
 * low child (skip:thm:canon), as unique as the paper's, so equality stays
 * one comparison.
 *
 * Conventions, fixed once. Qubit 0 is the top variable and is decided first.
 * A vector of length 2^n is indexed by x with qubit q at bit (n-1-q), so the
 * top qubit is the most significant bit and the first half of a vector is the
 * cofactor at x_top = 0. Translations and sign patterns are masks in that
 * same convention, and the levels an edge skips are the leading bits of the
 * index, over which the vector of its node repeats. Two differences from the
 * paper's counts: the paper keeps the redundant nodes, which this diagram
 * contracts, so it has at least as many nodes at every level; and it counts
 * the function 0 as a node at every level where it occurs (the +2 of
 * thm:pcoset is the code-state node and 0), where here zero is an edge, never
 * a node. limdd_level_counts is therefore the paper's count less the
 * redundant nodes, and less one at every level with a zero cofactor, and a
 * comparison with the paper's tables has to add both back.
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
 * Create the tables through limdd_nodes_init, and fix the label family, with
 * the default zero rule, SM, in every family (bqd_init_rule). A vector is
 * indexed by a 64-bit word, so `nqubits` may not exceed 63. Must be
 * called from a Lace worker; the edge weight table must already exist and
 * must hold complex weights. Not to be combined with limdd_nodes_init in one
 * session: the two share the tables. In the translation and Pauli families
 * pauli_tablesize may not exceed 2^32, since their memo keys hold the index
 * of a Pauli word in 32 bits. A session starts with an empty operation cache
 * (limdd_nodes_init), so one session's memo cannot answer for the next.
 */
void bqd_init(bqd_family_t family, size_t nqubits, size_t node_tablesize,
              size_t pauli_tablesize, size_t lim_tablesize, size_t stab_tablesize);
void bqd_quit(void);

bqd_family_t bqd_family(void);
const char  *bqd_family_name(bqd_family_t f);

/**
 * What a node does where its low cofactor is zero, the zero rule of the
 * session. BQD_ZERO_COPY is def:bqd: a quotient node everywhere, its ratio a
 * copy of the high cofactor where the low one is zero. BQD_ZERO_SM is rule SM
 * of skip:def:sm (qsylvan_bqd_sm.h, and qsylvan_bqd_xp_sm.h for the labels of
 * skip:def:smxp): a quotient node where the support of the high cofactor,
 * moved by its least point in the translation and Pauli families, lies inside
 * that of the low one, and a Shannon node, which stores the two cofactors,
 * where it does not; there is no copy. The two give
 * one diagram on every function of full support and on every function whose
 * support is an affine subspace (skip:cor:smd0), and different ones
 * elsewhere, so a session has one rule, fixed with its tables.
 */
typedef enum {
    BQD_ZERO_SM   = 0,
    BQD_ZERO_COPY = 1,
} bqd_zero_rule_t;

/**
 * bqd_init with the zero rule named, in any family: SM of the scalar family
 * is qsylvan_bqd_sm.h, SM of the translation and Pauli families
 * qsylvan_bqd_xp_sm.h, and the copy rule is the rest of the BQD's code.
 * bqd_init is this call with SM; a session that wants the copy rule, as the
 * paper's def:bqd has it, asks for it here.
 */
void bqd_init_rule(bqd_family_t family, bqd_zero_rule_t rule, size_t nqubits,
                   size_t node_tablesize, size_t pauli_tablesize, size_t lim_tablesize,
                   size_t stab_tablesize);

bqd_zero_rule_t bqd_zero_rule(void);
const char     *bqd_zero_rule_name(bqd_zero_rule_t r);

/* Written once by bqd_init_rule and read afterwards, by bqd_sm and bqd_xp_sm. */
extern bool bqd_sm_session;
extern bool bqd_xp_sm_session;

/**
 * Whether the session is the scalar family's under rule SM. The public entry
 * points ask this first and go to qsylvan_bqd_sm.h when it holds.
 */
static inline bool
bqd_sm(void)
{
    return bqd_sm_session;
}

/**
 * Whether the session is the translation or the Pauli family's under rule SM.
 * The public entry points ask this second and go to qsylvan_bqd_xp_sm.h when
 * it holds; the copy rule's code is reached only where neither does.
 */
static inline bool
bqd_xp_sm(void)
{
    return bqd_xp_sm_session;
}

/* --- labels --------------------------------------------------------------- */

/**
 * The label c Z^s X^t, with s and t as vector-index masks of an `nqubits`
 * vector. For the scalar family s and t must be 0; for the X family s must
 * be 0. Interned in the LIM table; the zero weight gives the zero label.
 */
LIMDD_LIM bqd_lim_make(EVBDD_WGT c, uint64_t s, uint64_t t, uint32_t nqubits);

/** Read a label back as (c, s, t). The zero label reads as (0, 0, 0). */
void bqd_lim_masks(LIMDD_LIM lim, uint32_t nqubits, EVBDD_WGT *c, uint64_t *s, uint64_t *t);

/**
 * The label c Z^s X^t with s and t as LIM masks, bit q for qubit q, which is
 * what the recursions hold. Interned through a per-worker, direct-mapped memo
 * in front of the LIM table, so that a label met again costs no shared
 * access; a scalar label is interned in one table lookup, without its word.
 */
LIMDD_LIM bqd_lim_word(EVBDD_WGT c, uint64_t s, uint64_t t);

/* --- counts, for diagnostics and the size of the cache -------------------- */

/**
 * What the operations did, counted per worker on a padded line of its own and
 * summed when read, so that counting shares nothing: the high cofactors made
 * (Cof1, and Join in the translation and Pauli families), the misses of the
 * four pointwise operations, of Canon, of the recursions that avoid the high
 * cofactor (Perm, Pair, X, CanonT) and of the phase multiplications. A test
 * bounds a gate's work by these rather than by its time, and bqd_cache_fit
 * sizes the operation cache by them.
 */
typedef enum {
    BQD_COUNT_COF1 = 0,
    BQD_COUNT_APPLY,
    BQD_COUNT_CANON,
    BQD_COUNT_PERM,
    BQD_COUNT_PHASEMUL,
    BQD_COUNTS
} bqd_count_t;

typedef struct {
    uint64_t n[BQD_COUNTS];
    char pad[SYLVAN_SHARING_PAD - (BQD_COUNTS * sizeof(uint64_t)) % SYLVAN_SHARING_PAD];
} bqd_count_line_t;

extern bqd_count_line_t *bqd_count_lines;

/** Count one k on `worker`'s line, LACE_WORKER_ID from inside a task. */
static inline void
bqd_count(unsigned worker, bqd_count_t k)
{
    bqd_count_lines[worker].n[k]++;
}

/** Zero the counts, and read their sums over the workers; between operations only. */
void bqd_counts_reset(void);
void bqd_counts_read(uint64_t out[BQD_COUNTS]);

/**
 * Between gates: fit Sylvan's operation cache to what the gates since the last
 * call made, and start the counts above afresh. Their sum is the number of
 * results the memoised operations made. Where it reaches a sixteenth of the
 * cache, the cache grows to the least power of two of which it is less than a
 * sixteenth, up to the maximum sylvan_set_sizes was given, and keeps what it
 * holds (cache_grow). A cache of fixed size, minimum and maximum equal, is left
 * alone. The maximum is reserved from the start and costs nothing until the
 * size in use reaches into it: a gate, and a collection that empties the
 * cache, touch the size in use and not the maximum. The runner calls it after
 * every gate; it needs every worker idle, as a collection does.
 */
void bqd_cache_fit(void);

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
 * Whether arg(w) lies in [0, pi), for w nonzero: the test of def:prep at a
 * pivot, which the Pauli family's Compose repeats (qsylvan_bqd_xp.h). Exact
 * for the algebraic weights, and only as good as the rounding for floats.
 */
bool bqd_arg_in_upper(EVBDD_WGT w);

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
