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
 * QUASI-REDUCED, NOT FULLY REDUCED
 *
 * Every path from the root to the terminal visits every level. In particular
 * limdd_makenode does NOT collapse a node whose children coincide.
 *
 * This is the one place where copying EVBDD would be wrong. An EVBDD skips a
 * level when the amplitude does not depend on that qubit and synthesises it
 * again on the way down, so `if (low == high) return low` is sound there. For
 * a LIMDD it is not: a node with low == high denotes |0>(x)v + |1>(x)v, which
 * is a state on one more qubit than v, not v. Dropping it would silently
 * change the state. The LIMDD normalisation rules are also stated over a
 * quasi-reduced structure, so keeping every level is what lets them apply.
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

/** The terminal node, representing the scalar 1. */
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

/** The edge labelled with the identity, pointing at the terminal. */
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
 * an error to pass anything else, and an assertion catches it.
 *
 * The node is interned structurally: equal (var, low, high) always give the
 * same index. That makes the STORAGE canonical. It does not by itself make the
 * diagram canonical -- that additionally needs the LIM normalisation rules,
 * which are not implemented yet.
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
