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
 * The LIMDD canonical form.
 *
 * limdd_makeedge is what limdd_makenode is not: it does not store the node it
 * is handed, it stores the canonical representative of every node denoting a
 * state that one is related to by a LIM, and returns the label that makes up
 * the difference. Two calls whose arguments denote LIM-related states come
 * back pointing at the SAME node. That is the whole content of the data
 * structure -- it is why a LIMDD is smaller than a QMDD, and everything above
 * this line relies on it.
 *
 * WHAT IS FREE, AND WHAT EACH FREEDOM COSTS
 *
 * Given |0>(x)A|v0> + |1>(x)B|v1>, four things can change without changing the
 * state beyond a LIM, and the canonical form spends one rule on each:
 *
 *   the low label       A is divided out and handed to the parent, so the
 *                       stored low edge is always the identity. This is why a
 *                       node needs only one label.
 *
 *   the high label      What is left, A^-1 B, is replaced by the least member
 *                       of its class under Stab(v0) and Stab(v1) (see
 *                       limdd_stab_min_coset). Without this two nodes that are
 *                       the same state in disguise keep different labels and
 *                       never merge.
 *
 *   which child is which  X at this level swaps the branches, so the node with
 *                       the children exchanged is equally valid and must be
 *                       considered; the smaller of the two is kept and X goes
 *                       to the parent. Skipping this loses exactly the merges
 *                       that make stabiliser states compact.
 *
 *   a dead branch's side  |0>(x)v and |1>(x)v are related by X too, so a node
 *                       with one zero child is always stored with the zero on
 *                       the high side. |0..0> and |10..0> are the same node.
 *
 * WHAT IS NOT DONE HERE
 *
 * The LABEL on the returned edge is not itself canonical: it may be multiplied
 * by any element of the result node's stabiliser group and still denote the
 * same state. That does not affect node merging -- the next level divides the
 * label out again and minimises over the same coset -- so diagrams stay
 * canonical and stay small. It does mean two edges denoting one state can
 * differ, so comparing states by edge equality needs the label reduced against
 * the target's group first.
 *
 * Weight normalisation follows from the rules above rather than being a rule
 * of its own: the low edge carries the identity, so the node's weights are
 * (1, scalar of the high label) and all scale sits on the parent. The LIMDD
 * paper instead keeps |a0|^2 + |a1|^2 = 1, which needs a weight and a Pauli
 * stored separately per edge; here the two are fused into one interned LIM so
 * that an edge is 64 bits, and there is no room for a second weight. Both are
 * canonical; theirs is the better conditioned.
 *
 * CONCURRENCY
 *
 * No locks and no shared mutable state beyond the lock-free tables. Each
 * node's stabiliser group is computed on demand and cached; two workers racing
 * on one node both compute the same group, so the duplicate work is wasted but
 * never wrong.
 */

#ifndef QSYLVAN_LIMDD_CANON_H
#define QSYLVAN_LIMDD_CANON_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "qsylvan_limdd_node.h"
#include "qsylvan_limdd_stab.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * The canonical edge denoting |0>(x)low + |1>(x)high, where both are edges at
 * level `var`+1.
 *
 * Must be called from a Lace worker.
 */
LIMDD limdd_makeedge(uint32_t var, LIMDD low, LIMDD high);

/**
 * The stabiliser group of node `p`, computed on first use and cached.
 *
 * The terminal has the trivial group. Nodes must have been built by
 * limdd_makeedge, since limdd_stab_of_node needs canonical children.
 */
LIMDD_STAB limdd_node_stab(LIMDD_TARG p);

/** The group of the state an edge denotes, i.e. its target's group conjugated
 *  by its label. The zero edge has the trivial group. */
LIMDD_STAB limdd_edge_stab(LIMDD e);

#ifdef __cplusplus
}
#endif

#endif // QSYLVAN_LIMDD_CANON_H
