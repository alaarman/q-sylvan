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
 *   the level itself    When the two canonical children coincide the node
 *                       would denote (|0>+|1>)(x)w, and it is not stored at
 *                       all: the edge points straight at w, and the level is
 *                       read as |0>+|1> (see qsylvan_limdd_node.h). The sign
 *                       rule folds (|0>-|1>)(x)w into the same case, with a Z
 *                       on the parent edge, so |+> and |-> registers cost no
 *                       nodes. |0>+i|1> is a product state too but no Pauli
 *                       image of |0>+|1>, and keeps its level.
 *
 * THE LABEL TOO
 *
 * The label on the returned edge is reduced against the target's stabiliser
 * group, which is the one remaining freedom: L and L G denote the same state
 * for every G that fixes the target. Node merging does not need it -- the next
 * level divides the label out again -- but state equality does, and without it
 * two edges for one state can differ.
 *
 * LEVELS
 *
 * An edge does not know what level it is read at, and once levels can be
 * skipped the group it is reduced against depends on that: an edge read at
 * level k whose target sits at level k' is also fixed by X_j for every
 * skipped j in [k, k'). So limdd_edge_canonical and limdd_edge_stab take the
 * level. makeedge uses `var` for the edge it returns and `var + 1` for the
 * children it was handed; a root is at 0. Too LARGE a level is caught by an
 * assertion against the target; too SMALL a level is not detectable, and
 * silently reduces the label at levels the edge does not own.
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
 * Turn the high-determinism rule on or off. On by default.
 *
 * With it off, makeedge stops choosing a canonical high label: it divides the
 * low label out, applies the zero-edge rule, and interns what is left. The
 * coset minimisation, the child swap and the stabiliser groups are all
 * skipped, so making a node costs a LIM division instead of O(n^3) linear
 * algebra.
 *
 * What is lost is merging. Two nodes denoting the same state up to a LIM keep
 * whatever labels they were handed and stay distinct, so the diagram is bigger
 * -- for stabiliser states, potentially exponentially so, which is the whole
 * reason LIMDDs exist. Soundness is unaffected: the edges still denote the
 * states they are supposed to.
 *
 * Node stabiliser groups are not computed while this is off, so
 * limdd_node_stab will derive them on demand and its cache stays cold.
 * QolDDer calls the same switch high_determinism.
 */
void limdd_set_high_determinism(bool on);
bool limdd_get_high_determinism(void);

/**
 * When the canonical form is applied.
 *
 * ALWAYS is the original behaviour: every makeedge does the coset search, so
 * the diagram is canonical after every operation. The rest defer it, running
 * operations with high determinism off and rebuilding later. Deferring trades
 * a smaller number of expensive canonicalisations for a larger diagram in
 * between -- and the diagram in between is about the size a QMDD would be,
 * since without the orbit search nothing LIM-related merges.
 */
typedef enum {
    LIMDD_CANON_ALWAYS,    /* canonicalise inside every makeedge */
    LIMDD_CANON_MANUAL,    /* never automatically; the caller decides */
    LIMDD_CANON_OPS,       /* when `interval` nodes have been built */
    LIMDD_CANON_ADAPTIVE,  /* interval tuned by how much the last one helped */
} limdd_canon_policy_t;

/**
 * Choose the policy. `interval` counts NODE-BUILDING operations (makeedge
 * calls), not gates, so the trigger follows the work a circuit does rather
 * than how many gates it is written with. It is the threshold for OPS and the
 * starting threshold for ADAPTIVE, and is ignored otherwise.
 *
 * Anything but ALWAYS turns high determinism off for ordinary operations;
 * limdd_canonize turns it back on for the duration of a rebuild.
 */
void limdd_set_canon_policy(limdd_canon_policy_t policy, uint64_t interval);

/**
 * Count one unit of work towards the next rebuild.
 *
 * Safe from inside a parallel operation: each worker has its own counter on
 * its own cache line, so incrementing never invalidates another's. A single
 * shared counter here would be read-modify-written once per unit of work by
 * every worker, which costs far more than the rebuild it is scheduling.
 */
void limdd_canon_count(void);

/**
 * Sum the per-worker counters and say whether a rebuild is due.
 *
 * Called by whoever drives the operations, BETWEEN them -- never inside one.
 * limdd_canonize rebuilds nodes bottom-up, so it needs a root that is not in
 * the middle of being built. Summing here rather than at each increment is
 * what keeps the counting off the hot path.
 */
bool limdd_canon_due(void);

/**
 * The canonical diagram denoting the same state as `e`.
 *
 * A bottom-up rebuild: each node's children are canonicalised first, then
 * makeedge is applied to them, which is exactly the precondition makeedge
 * needs. The label that makeedge hands back is composed into the edge above,
 * so the state is unchanged.
 *
 * Building with the canonical form off and then calling this gives the SAME
 * edge as building with it on throughout.
 *
 * Must be called from a Lace worker.
 */
LIMDD limdd_canonize(LIMDD e);

/** Nodes reached by the last limdd_canonize, before and after. */
void limdd_canon_last(uint64_t *before, uint64_t *after);

/**
 * The canonical edge denoting |0>(x)low + |1>(x)high, where both are edges
 * read at level `var`+1 -- they may point below that, and so may the edge
 * returned, which is read at `var` and skips `var` itself when the two
 * children turn out to be one state.
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

/**
 * `e`, read at `level`, with its label replaced by the least member of
 * L * Stab(target at level).
 *
 * Every member of that coset denotes the same state, so an edge is only
 * determined up to it; reducing makes two edges for one state identical, which
 * is what comparing states by edge equality needs. makeedge applies this to
 * everything it returns, so edges from it are already reduced.
 *
 * The group includes X_j for each level j the edge skips, and because the X
 * columns lead the RREF that clears the label's X bit at every such level: a
 * skipped level ends up carrying I or Z, with a Y folded into -i times Z.
 */
LIMDD limdd_edge_canonical(uint32_t level, LIMDD e);

/**
 * The group of the state an edge read at `level` denotes: its target's group,
 * extended over the levels the edge skips, conjugated by its label. A skipped
 * level contributes +X_j, or -X_j where the label has Z or Y there, since
 * that level then holds |0>-|1>. The zero edge has the trivial group.
 */
LIMDD_STAB limdd_edge_stab(uint32_t level, LIMDD e);

#ifdef __cplusplus
}
#endif

#endif // QSYLVAN_LIMDD_CANON_H
