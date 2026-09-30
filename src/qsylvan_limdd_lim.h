/*
 * Copyright 2026 System Verification Lab, LIACS, Leiden University
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
 * Interned LIMs: the labels on LIMDD edges.
 *
 * A Pauli-LIMDD labels each edge with a Linearly Invertible Map, which for the
 * Pauli group is a complex scalar times a tensor product of Pauli matrices.
 * Carrying that inline would make an edge far too wide -- the Pauli word alone
 * is two machine words -- so both halves are interned and an edge carries one
 * index, exactly as an EVBDD edge carries an index into the edge weight table.
 *
 * TWO TABLES
 *
 * Both are ordinary Sylvan llmsset instances, because both keys are exactly
 * the 16 bytes an llmsset bucket holds:
 *
 *   Pauli table:  (x mask, z mask)          -> LIMDD_PAULI_REF
 *   LIM table:    (Pauli ref, EVBDD_WGT)    -> LIMDD_LIM
 *
 * so lookup is Sylvan's existing lock-free CAS insert, with no custom leaf
 * callbacks, no side allocation and nothing to free. The scalar half reuses
 * Q-Sylvan's edge weight table, so LIMDD amplitudes and EVBDD/QMDD amplitudes
 * agree on when two weights are the same.
 *
 * This is what keeps a LIMDD edge the same shape as an EVBDD edge: a single
 * 64-bit [label index | node index] pair.
 *
 * CANONICITY
 *
 * Two rules are enforced when a LIM is made, so that equal maps always get
 * equal indices:
 *
 *   - A LIM with scalar zero is the zero map whatever its Pauli word is, so
 *     every such LIM collapses to LIMDD_LIM_ZERO. Without this the same zero
 *     map would get as many indices as there are Pauli words, and node lookup
 *     would stop finding existing nodes. (paulilim documents that its MakeEdge
 *     leaves the Pauli word unspecified when the scalar is zero, so this case
 *     genuinely arises.)
 *   - The phase from multiplying Pauli words is folded into the scalar rather
 *     than stored, so there is one representation of each map, not four.
 *
 * CONCURRENCY
 *
 * limdd_lim_t values are plain 64-bit indices and limdd_pauli_t values are
 * plain pairs of words: they are computed in registers, passed by value, and
 * never point at anything. The only shared state is the two tables and the
 * edge weight table, all of which are lock-free. Nothing here holds per-worker
 * mutable state, so there is no line for two workers to contend over.
 *
 * Interning requires a Lace worker, because llmsset claims its buckets from a
 * per-worker region. Call these from inside a task.
 */

#ifndef QSYLVAN_LIMDD_LIM_H
#define QSYLVAN_LIMDD_LIM_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "sylvan_edge_weights.h"
#include "qsylvan_limdd_pauli.h"

#ifdef __cplusplus
extern "C" {
#endif

/** An interned Pauli word. */
typedef uint64_t LIMDD_PAULI_REF;

/** An interned LIM: a scalar together with a Pauli word. */
typedef uint64_t LIMDD_LIM;

/**
 * The identity map, 1 * II..I. This is by far the most common label, since
 * every low edge of a normalised node carries it.
 */
extern LIMDD_LIM LIMDD_LIM_IDENTITY;

/**
 * The zero map. Every LIM whose scalar is zero collapses to this one.
 */
/**
 * Width of the LIM index where an edge carries it, and hence the largest
 * number of LIMs the table may hold. It lives here rather than with the edge
 * layout because it bounds this table: an index at or above LIMDD_LIM_MAX
 * cannot be put in an edge at all.
 */
#define LIMDD_LIM_BITS  32
#define LIMDD_LIM_MAX   ((uint64_t)1 << LIMDD_LIM_BITS)

extern LIMDD_LIM LIMDD_LIM_ZERO;

/**
 * Create the Pauli and LIM tables.
 *
 * `nqubits` is the width of the Pauli words that will be stored, used to check
 * the words handed in; it may not exceed LIMDD_MAX_QUBITS. The edge weight
 * table must already exist (sylvan_init_edge_weights), and the weights must be
 * complex, since a Pauli product can produce a factor of i.
 *
 * Must be called from a Lace worker.
 */
void limdd_lims_init(size_t nqubits, size_t pauli_tablesize, size_t lim_tablesize);

/**
 * Destroy both tables. The edge weight table is not touched.
 */
void limdd_lims_quit(void);

/** The qubit width passed to limdd_lims_init. */
size_t limdd_lims_nqubits(void);

/** Number of distinct Pauli words interned so far. Needs a running Lace. */
size_t limdd_pauli_table_count(void);

/** Capacity of the Pauli table, for deciding when to collect. */
size_t limdd_pauli_table_size(void);

/** Number of distinct LIMs interned so far. Needs a running Lace. */
size_t limdd_lim_table_count(void);

/** Capacity of the LIM table, for deciding when to collect. */
size_t limdd_lim_table_size(void);

/**
 * Intern `p`, returning its reference. Equal words always give equal
 * references. Aborts if the table is full.
 */
LIMDD_PAULI_REF limdd_pauli_intern(limdd_pauli_t p);

/** The Pauli word `ref` refers to. */
limdd_pauli_t limdd_pauli_deref(LIMDD_PAULI_REF ref);

/**
 * The LIM `w * p`. Equal maps always give equal LIMs; in particular any `p`
 * with `w == EVBDD_ZERO` gives LIMDD_LIM_ZERO.
 */
LIMDD_LIM limdd_lim_make(limdd_pauli_t p, EVBDD_WGT w);

/**
 * The LIM `w * I`, the one limdd_lim_make gives for the identity word, in one
 * table lookup: the identity word is not built and not looked up.
 */
LIMDD_LIM limdd_lim_scalar(EVBDD_WGT w);

/**
 * The generation of the LIM indices, which moves on at every collection and
 * every new pair of tables: a memo of LIMs kept outside the table holds for
 * one generation only, since a swept bucket is handed out again and a
 * collection may renumber the weights a LIM holds.
 */
uint64_t limdd_lim_generation(void);

/** The Pauli word of `lim`. For LIMDD_LIM_ZERO this is the identity. */
limdd_pauli_t limdd_lim_pauli(LIMDD_LIM lim);

/**
 * The reference of the Pauli word of `lim`, as limdd_pauli_intern returned it:
 * one read of the LIM's bucket, where limdd_lim_pauli also reads the word.
 */
LIMDD_PAULI_REF limdd_lim_pauli_ref(LIMDD_LIM lim);

/** The scalar of `lim`. */
EVBDD_WGT limdd_lim_weight(LIMDD_LIM lim);

/** True iff `lim` is the zero map. */
static inline bool
limdd_lim_is_zero(LIMDD_LIM lim)
{
    return lim == LIMDD_LIM_ZERO;
}

/** True iff `lim` is the identity map. */
static inline bool
limdd_lim_is_identity(LIMDD_LIM lim)
{
    return lim == LIMDD_LIM_IDENTITY;
}

/**
 * The composition `a * b`, as matrices, in that order.
 *
 * Pauli matrices do not commute, so this is not symmetric: limdd_lim_mul(a,b)
 * and limdd_lim_mul(b,a) can differ by a sign. The power of i that the Pauli
 * product produces is folded into the scalar here, so the result is already
 * canonical.
 */
LIMDD_LIM limdd_lim_mul(LIMDD_LIM a, LIMDD_LIM b);

/**
 * The inverse of `lim`, i.e. the LIM `m` with `lim * m` the identity map.
 *
 * A Pauli word is its own inverse, so only the scalar is inverted. Inverting
 * the zero map is not defined; that returns LIMDD_LIM_ZERO.
 */
LIMDD_LIM limdd_lim_inverse(LIMDD_LIM lim);

/**
 * The edge weight i^k, for k in 0..3.
 *
 * Applying a Pauli word to a basis state produces such a factor, so callers
 * that evaluate a diagram need them; they are interned once at init.
 */
EVBDD_WGT limdd_wgt_i_pow(unsigned k);

/* --- garbage collection --------------------------------------------------- */

/**
 * Mark `lim` live, and with it the Pauli word it names.
 *
 * The two tables are collected together because nothing else refers to a Pauli
 * word: it is reachable only through the LIMs that use it.
 */
void limdd_gc_mark_lim(LIMDD_LIM lim);

/** Forget which buckets of both tables are occupied. */
void limdd_gc_clear_lims(void);

/** Rebuild both hash arrays, and re-mark the reserved LIMs. */
void limdd_gc_rehash_lims(void);

/**
 * Move the edge weights of every live LIM into a fresh weight table, and
 * discard the old one. Returns the number kept.
 *
 * Nothing else reclaims them. The weight table is collected by copying, which
 * assigns new indices, and Q-Sylvan's own collector rewrites every EVBDD to
 * match -- a path a LIMDD never takes, so limdd_nodes_init switches it off and
 * the table then only grows. Measured, 95% to 99.8% of it is dead.
 *
 * Copying is safe here for a reason particular to llmsset: a bucket keeps its
 * index across a collection. A LIM is interned under (Pauli, weight), so a new
 * weight index changes the bucket's CONTENTS but not the index it sits at, and
 * every node referring to that LIM stays valid. Only the LIM table's hash has
 * to be rebuilt -- which limdd_gc_rehash_lims already does, so this must run
 * between marking and that rehash, and nowhere else.
 */
size_t limdd_gc_remap_weights(void);

/** Print `lim` as "(re,im) * IXYZ" to `out`, without a trailing newline. */
void limdd_lim_fprint(FILE *out, LIMDD_LIM lim);

#ifdef __cplusplus
}
#endif

#endif // QSYLVAN_LIMDD_LIM_H
