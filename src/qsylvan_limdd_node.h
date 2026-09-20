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
 * LIMDD nodes and edges.
 *
 * EDGES
 *
 * A LIMDD edge is one 64-bit word holding an interned LIM and a node index,
 * laid out exactly like an EVBDD edge with the edge weight generalised to a
 * whole LIM:
 *
 *      1 bit  unused
 *     23 bits index of the LIM labelling this edge
 *     40 bits index of the node this edge points at
 *
 * Interning the LIM (see qsylvan_limdd_lim.h) is what makes that fit: the
 * scalar and the Pauli word together are three machine words and could never
 * be carried inline.
 *
 * NODES
 *
 * A node is 16 bytes, one llmsset bucket, as two 64-bit words:
 *
 *   low word:   1 bit  unused
 *              16 bits variable (qubit) of this node
 *               1 bit  set iff the low edge is the zero map
 *               6 bits unused
 *              40 bits low edge target
 *
 *   high word:  1 bit  mark, for garbage collection
 *              23 bits LIM labelling the high edge
 *              40 bits high edge target
 *
 * Only the high edge carries a LIM. That is not a space trick, it is what
 * LIMDD normalisation gives you: the low edge's label is factored out to the
 * parent, leaving the low edge labelled with the identity, or with zero when
 * that branch vanishes. One bit distinguishes those two cases.
 *
 * FULLY REDUCED: LEVELS MAY BE SKIPPED
 *
 * An edge does not record what level it is at; that comes from context. The
 * root edge is at level 0, a node's child edges are at the node's level plus
 * one, and the terminal counts as level nqubits (see limdd_level). An edge
 * read at level k that points at a node of level k' > k SKIPS levels k..k'-1,
 * and each skipped level denotes the unnormalised |0>+|1>: the edge means
 *
 *     L ( (|0>+|1>)_k (x) ... (x) (|0>+|1>)_{k'-1} (x) |v> )
 *
 * with L the edge's LIM and |v> the target's state. L is a word over all n
 * qubits, so it has an entry at every skipped level too, and that entry acts
 * on the |0>+|1> there: I and X leave it, Z turns it into |0>-|1>, and Y into
 * -i(|0>-|1>). That is all a skipped level can carry -- the Pauli orbit of
 * |0>+|1> -- and it is why a level may be skipped only where the node would
 * have had two children that are one and the same state.
 *
 * Consequences worth keeping in mind:
 *
 *   - A node's two children may sit at different levels.
 *   - `var` is stored in the node so that an edge can tell how many levels it
 *     skips. It is the one thing context cannot supply.
 *   - An edge into the terminal denotes the scalar 1 only when read at level
 *     nqubits. Read at level 0, the same word (identity, terminal) is the
 *     unnormalised |+>^n. Nothing about the edge distinguishes the two.
 *   - Reading an edge at a level ABOVE its true one is a silent error: the
 *     reader synthesises |0>+|1> on levels the label already acts on. A level
 *     below the target's is caught by an assertion; a level below the edge's
 *     own is not detectable from the edge alone.
 *
 * Skipping is limdd_makeedge's decision, not limdd_makenode's. makenode stores
 * exactly what it is handed, and `if (low == high) return low` is still wrong
 * here: a bare node index does not say at what level it is read, and only the
 * parent, which knows the level and holds the label, can put the skipped
 * level's Pauli where it belongs. The EVBDD version of that rule is sound
 * because an EVBDD skips a level exactly when the amplitude ignores the
 * qubit; the LIMDD version has to account for the LIM as well.
 *
 * VARIABLE ORDER
 *
 * Variable k is qubit k, variable 0 is the root level, and bit k of a basis
 * state is qubit k. Note this is the opposite of QolDDer, which puts qubit
 * n-1 at the top; anything converting between the two has to flip.
 *
 * CONCURRENCY
 *
 * Edges and node indices are plain 64-bit values. The only shared state is the
 * node table, the two LIM tables and the edge weight table, all lock-free.
 * Nothing here is per-worker and mutable, so there is no line to false-share.
 * Making a node claims a bucket from a per-worker region, so these must be
 * called from a Lace worker.
 */

#ifndef QSYLVAN_LIMDD_NODE_H
#define QSYLVAN_LIMDD_NODE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "sylvan_edge_weights.h"
#include "qsylvan_limdd_lim.h"
#include "qsylvan_limdd_pauli.h"

#ifdef __cplusplus
extern "C" {
#endif

/** A LIMDD edge: an interned LIM together with a target node. */
typedef uint64_t LIMDD;

/** A LIMDD node index. */
typedef uint64_t LIMDD_TARG;

/**
 * The terminal node: the empty tensor product, at level nqubits. An edge into
 * it read at level k denotes the label applied to the unnormalised |+>^(n-k),
 * which is the scalar 1 only when k == nqubits.
 */
#define LIMDD_TERMINAL ((LIMDD_TARG)1)

/** Widths of the two fields of an edge. */
#define LIMDD_LIM_BITS  23
#define LIMDD_TARG_BITS 40

#define LIMDD_TARG_MASK (((uint64_t)1 << LIMDD_TARG_BITS) - 1)
#define LIMDD_LIM_MAX   ((uint64_t)1 << LIMDD_LIM_BITS)
#define LIMDD_TARG_MAX  LIMDD_TARG_MASK

/**
 * Create the LIMDD node table, and the LIM, Pauli and stabiliser tables under it.
 *
 * The edge weight table must already exist and must hold complex weights.
 * Must be called from a Lace worker.
 */
void limdd_nodes_init(size_t nqubits, size_t node_tablesize,
                      size_t pauli_tablesize, size_t lim_tablesize,
                      size_t stab_tablesize);

/** Destroy the node, LIM, Pauli and stabiliser tables. Weights are untouched. */
void limdd_nodes_quit(void);

/** Number of nodes currently in the table, terminal excluded. */
size_t limdd_node_table_count(void);

/** Pack a LIM and a target into an edge. */
static inline LIMDD
limdd_bundle(LIMDD_LIM lim, LIMDD_TARG target)
{
    return (lim << LIMDD_TARG_BITS) | target;
}

/** The LIM labelling `e`. */
static inline LIMDD_LIM
limdd_label(LIMDD e)
{
    return e >> LIMDD_TARG_BITS;
}

/** The node `e` points at. */
static inline LIMDD_TARG
limdd_target(LIMDD e)
{
    return e & LIMDD_TARG_MASK;
}

/**
 * The zero edge. Its label is LIMDD_LIM_ZERO and it points at the terminal;
 * pinning the target matters, because otherwise the same zero vector could be
 * written in as many ways as there are nodes to point at.
 */
LIMDD limdd_zero_edge(void);

/**
 * The edge labelled with the identity, pointing at the terminal. Read at level
 * nqubits it is the scalar 1; read at level k it is the unnormalised |+>^(n-k).
 */
LIMDD limdd_one_edge(void);

/** True iff `e` denotes the zero vector. */
static inline bool
limdd_edge_is_zero(LIMDD e)
{
    return limdd_lim_is_zero(limdd_label(e));
}

/**
 * Find or create the node (var, low, high), returning its index.
 *
 * `low` must be labelled with the identity or be the zero edge: the caller is
 * responsible for having factored any other low label out to the parent. It is
 * an error to pass anything else, and an assertion catches it. Both children
 * must point strictly below `var`, i.e. var < limdd_level(target); they need
 * not point at var+1, since an edge may skip levels.
 *
 * The node is interned structurally: equal (var, low, high) always give the
 * same index. That makes the STORAGE canonical. It does not by itself make the
 * diagram canonical -- that is limdd_makeedge's job, and so is deciding
 * whether this level should exist at all.
 *
 * Must be called from a Lace worker.
 */
LIMDD_TARG limdd_makenode(uint32_t var, LIMDD low, LIMDD high);

/**
 * As limdd_makenode, but `*created` reports whether the bucket was fresh.
 *
 * Anything caching data alongside a node needs this: after a collection a
 * bucket can be handed out again, and a cache slot left over from its previous
 * occupant would otherwise be read as that node's.
 */
LIMDD_TARG limdd_makenode_ex(uint32_t var, LIMDD low, LIMDD high, int *created);

/** The variable of node `p`. `p` must not be the terminal. */
uint32_t limdd_node_var(LIMDD_TARG p);

/**
 * The level an edge into `p` descends to: `p`'s variable, or nqubits for the
 * terminal. An edge read at level k skips levels k..limdd_level(target)-1.
 */
uint32_t limdd_level(LIMDD_TARG p);

/** The low edge of node `p`. */
LIMDD limdd_node_low(LIMDD_TARG p);

/** The high edge of node `p`. */
LIMDD limdd_node_high(LIMDD_TARG p);

/**
 * One machine word of storage per node, where the canonical form caches that
 * node's stabiliser group. Zero means "not computed yet"; the group handles
 * themselves are never zero.
 *
 * It lives here rather than in the node because the bucket is full, and it is
 * a plain array rather than a table because the key is already an index.
 *
 * Concurrency: entries are written with a relaxed atomic store, and two
 * workers racing on the same node both compute the same value, so the race is
 * benign and needs no lock. Neighbouring entries do share a cache line, but
 * llmsset hands each worker a 512-bucket region at a time, so concurrent
 * writers are normally far further apart than that.
 */
uint64_t limdd_node_stab_raw(LIMDD_TARG p);
void limdd_node_set_stab_raw(LIMDD_TARG p, uint64_t v);

/* --- garbage collection ---------------------------------------------------
 *
 * These are the node table's half of limdd_gc; see qsylvan_limdd_gc.h for what
 * drives them. llmsset keeps data at a fixed index across a collection, which
 * is what lets the stabiliser cache above survive one.
 */

/** Mark bucket `p` as live. Returns 1 if this call marked it, 0 if already. */
int limdd_gc_mark_node(LIMDD_TARG p);

/** Forget which buckets are occupied, so marking can rebuild the set. */
void limdd_gc_clear_nodes(void);

/** Rebuild the hash array over the marked buckets. */
void limdd_gc_rehash_nodes(void);

/**
 * Drop the cached stabiliser group of every bucket that marking did not reach.
 *
 * Their buckets are about to be handed out again, and a group left behind
 * would then be read as the new occupant's.
 */
void limdd_gc_purge_stab_cache(void);

/**
 * The amplitude that `e` assigns to the basis state `bits`, where bits[k] is
 * qubit k. `nqubits` entries are read.
 *
 * This is the semantic meaning of a diagram, and exists so that structural
 * changes -- normalisation rules, rewrites -- can be checked against what the
 * diagram actually denotes rather than against its shape.
 */
EVBDD_WGT limdd_eval(LIMDD e, const bool *bits, size_t nqubits);

#ifdef __cplusplus
}
#endif

#endif // QSYLVAN_LIMDD_NODE_H
