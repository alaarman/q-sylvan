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

/*
 * THE ISOMORPHISM-CLASS TABLE of the lazy canonical form (--canon=lazy).
 *
 * The lazy form keeps what canonicity needs -- one stored node per class of
 * states that are LIM images of each other -- and drops what it does not: the
 * stored node need not be the minimal one of its class, only the first one in.
 * A node with two different children is then stored as makeedge factors it,
 * without the coset search and without its stabiliser group, unless a node of
 * its class may already be stored. Whether one may is decided by a key that
 * every member of the class shares:
 *
 *     (var, {target(low), target(high)})        an unordered pair
 *
 * It is invariant because the children are class representatives themselves:
 * a LIM mu (P_top (x) P') maps the node's children to images of themselves
 * (P_top in {I, Z}) or swaps them (P_top in {X, Y}), and an image of a stored
 * class representative is that representative again. The high label's Pauli
 * word is left out on purpose: the children's stabiliser groups move it within
 * a double coset, so the word is not invariant.
 * Its scalar is invariant up to a unit. On exact weights every member is
 * stored with its children in index order (makeedge_lazy), so two members of one class are
 * related by a LIM with P_top in {I, Z}: the low branch gives mu R |v0> = |v0>,
 * so mu is an eigenvalue of the Pauli word R, +-1; the high branch then gives
 * lambda' P' |v1> = +-mu lambda R P |v1>, and R P and P'^-1 R P are Pauli words
 * up to a factor in {+-1, +-i}, the second with eigenvalue +-1 on v1. So
 * lambda' = lambda u for a unit u in {+-1, +-i}. On exact weights a hash of
 * the scalar up to units (limdd_class_scalar_key) picks the slot a class's
 * scan starts from, and filters the members exactly before the full
 * comparison. On floats rounding could break that invariant, so the scan
 * starts at slot 0, and the scalar's modulus (shared by the class, a unit
 * having modulus 1) filters with a tolerance far wider than any rounding
 * (see qsylvan_limdd_canon.c).
 *
 * Nodes with one child node twice (where the skip rule needs the child's
 * group, so the canonical form is computed anyway and stored, as the eager
 * form stores it) and nodes with a dead branch (one per child, so unique by
 * construction) are not in this table.
 *
 * The members of a class sit in slots s, s+1, s+2, ... (mod 2^20) of the key,
 * s the scalar key (0 on floats), each slot an llmsset bucket (var, lo, hi,
 * slot -> member) whose equality ignores the member, so the first worker to
 * ask for a slot owns it. Slots are shared by keys whose runs overlap, which
 * costs a comparison and nothing else. Every worker scans a class's slots in
 * the same order and compares against each member before it claims the next
 * free one; two isomorphic nodes therefore meet at the first slot either of
 * them takes, and only one of them is ever a member. No locks: the
 * llmsset's compare-and-swap on the hash bucket is the only synchronisation.
 *
 * Collection: the table holds node indices, and a collection recycles the
 * buckets of dead nodes, so limdd_class_gc_rebuild refills it from the marked
 * nodes after every LIMDD collection (stop the world), each from its own start
 * slot on, with no comparisons -- every marked node is a member, since makeedge never hands
 * out a node that lost its class to another.
 */

#ifndef QSYLVAN_LIMDD_CLASS_H
#define QSYLVAN_LIMDD_CLASS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "qsylvan_limdd_node.h"
#include "sylvan_edge_weights.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Creates the table with room for `tablesize` members (the node table's size). */
void limdd_class_init(size_t tablesize);

/**
 * Asks for the lazy form (limdd_set_canon_policy does). The table is created
 * with the node table's size as soon as both exist (limdd_nodes_attach_class),
 * and again for every later session until it is asked off.
 */
void limdd_class_request(bool on);
bool limdd_class_requested(void);

/** Whether the table exists, i.e. the lazy canonical form was turned on. */
bool limdd_class_on(void);

/**
 * The member in slot start + j (mod 2^20) of the class (var, {lo, hi}), lo < hi;
 * if the slot is free, `t` becomes its member and *claimed is set. Exits when
 * the table or the slot range (2^20 members per key) is full, like the other
 * LIMDD tables.
 */
LIMDD_TARG limdd_class_slot(uint32_t var, LIMDD_TARG lo, LIMDD_TARG hi, uint64_t start,
                            uint32_t j, LIMDD_TARG t, int *claimed);

/**
 * The high scalar of a member, stored with its children in index order, up to
 * the units +-1 and +-i, as a hash: equal for every member of a class. 0 unless
 * the weights are exact (qisq2).
 */
uint64_t limdd_class_scalar_key(LIMDD_LIM high);

/**
 * Whether the class of a node at `var` with the two different children v0 and
 * v1 goes through this table. Not when both children's groups are known
 * trivial (their inline bit, read at var+1 with no level skipped): then the
 * canonical label is a minimisation over two trivial groups, with no
 * elimination, and the node's own group is trivial too, so the canonical form
 * costs less than a lookup here and is applied as the eager form applies it.
 * A function of the two children alone, which every member of a class shares,
 * so each class is handled one way only. Only on exact weights: there the
 * canonical node is also the one makeedge_lazy would store (lazy_key), while
 * floats store the factored node unnormalised.
 */
static inline bool
limdd_class_wanted(uint32_t var, LIMDD_TARG v0, LIMDD_TARG v1)
{
    /* Floats keep the factored node as it comes (makeedge_lazy), and the
     * canonical path's swap would round it differently. */
    if (sylvan_get_edge_weight_type() != WGT_QISQ2) return true;
    return !(limdd_node_stab_is_trivial(v0) && limdd_level(v0) == var + 1 &&
             limdd_node_stab_is_trivial(v1) && limdd_level(v1) == var + 1);
}

/** Refills the table from the marked nodes; called by the LIMDD collector. */
void limdd_class_gc_rebuild(void);

/** Frees the table, which turns the lazy form off; called by limdd_nodes_quit. */
void limdd_class_quit(void);

#ifdef __cplusplus
}
#endif

#endif // QSYLVAN_LIMDD_CLASS_H
