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
 * Garbage collection for LIMDDs.
 *
 * Four tables have to be collected together, because each is reachable only
 * through the one above it and through nothing else:
 *
 *     nodes  -> the LIM on each high edge, and the node's stabiliser group
 *     groups -> a list cell per generator, each naming a LIM
 *     LIMs   -> a Pauli word and an edge weight
 *     Paulis -> nothing
 *
 * Collecting any one alone would either sweep something still in use or keep
 * everything. So marking walks from the protected roots through all four, and
 * all four are swept and rehashed in one pass.
 *
 * WHY INDICES SURVIVE
 *
 * llmsset leaves bucket contents where they are and only rebuilds the hash
 * array, so a surviving entry keeps its index. Every LIMDD edge, label and
 * group handle is such an index, which is what makes a mark-and-sweep possible
 * at all here: nothing has to be rewritten. The per-node stabiliser cache
 * survives for the same reason -- it is indexed by bucket -- and only the
 * entries of swept buckets are dropped, so that a recycled bucket does not
 * inherit its predecessor's group.
 *
 * WHAT IS NOT COLLECTED: EDGE WEIGHTS
 *
 * The scalar of a LIM lives in Q-Sylvan's edge weight table, which is shared
 * with EVBDD and QMDD and is collected by copying: surviving weights are moved
 * to a fresh table and every EVBDD is rewritten with their new indices. That
 * cannot be done to a LIMDD in place. A LIM is interned under (Pauli, weight),
 * so new weight indices mean new LIM indices, which mean new node contents and
 * so new node indices -- the whole forest would have to be rebuilt.
 *
 * limdd_gc therefore leaves the weight table alone, and limdd_nodes_init turns
 * OFF the automatic weight-table collection that evbdd_makenode would
 * otherwise trigger. Without that, an EVBDD operation elsewhere in the program
 * could silently invalidate every live LIMDD. The cost is that dead weights
 * accumulate; for Pauli-LIMDDs the scalars are mostly fourth roots of unity,
 * so the table grows slowly, but a long run can still exhaust it.
 *
 * WHEN IT MAY RUN
 *
 * This is a stop-the-world collection with explicit roots: nothing but the
 * protected edges survives, and intermediate results inside a half-finished
 * makeedge are not protected. It must not run while another worker is building
 * a diagram. Sylvan's own collector suspends every worker at a Lace frame
 * boundary to get that guarantee; wiring LIMDD operations into that needs them
 * to be Lace tasks carrying a reference stack, which they are not yet.
 */

#ifndef QSYLVAN_LIMDD_GC_H
#define QSYLVAN_LIMDD_GC_H

#include <stddef.h>

#include <lace.h>

#include "qsylvan_limdd_canon.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Keep the edge stored at `*a` across collections, and keep tracking it: the
 * variable is read at collection time, so assigning to it afterwards protects
 * the new value instead.
 */
void limdd_protect(LIMDD *a);

/** Stop protecting `*a`. */
void limdd_unprotect(LIMDD *a);

/** How many variables are currently protected. */
size_t limdd_count_protected(void);

/* --- per-worker reference stacks ------------------------------------------
 *
 * limdd_protect is for the few long-lived roots a caller holds -- a
 * simulator's state edge. It is a shared set, so protecting an intermediate
 * there would have every worker writing to one structure on the hot path.
 * These stacks are the opposite: private to a worker, and the mechanism that
 * lets a collection run WHILE an operation is in flight.
 *
 * Use is the same discipline as Sylvan's own, which these mirror:
 *
 *   limdd_refs_push(e)    keep this value across a collection
 *   limdd_refs_pop(n)     drop the last n
 *   limdd_refs_pushptr(&v) keep whatever the variable holds, re-read at
 *                         collection time, so assigning to v is safe
 *   limdd_refs_popptr(n)
 *   limdd_refs_spawn(SPAWN(...))  keep a branch's result while it runs
 *   e = limdd_refs_sync(SYNC(...))
 *
 * Must be called from a Lace worker; the stack lives in thread-local storage
 * and is padded to its own cache line, since the cursors are written on every
 * push and pop.
 */
/**
 * True when a table is close enough to full that an operation should stop and
 * collect. Cheap: two counter reads.
 */
bool limdd_gc_wanted(void);

void  limdd_refs_init(void);
LIMDD limdd_refs_push(LIMDD e);
void  limdd_refs_pop(long amount);
void  limdd_refs_pushptr(const LIMDD *ptr);
void  limdd_refs_popptr(size_t amount);
void  limdd_refs_spawn(Task *t);
LIMDD limdd_refs_sync(LIMDD result);

/**
 * Collect the four LIMDD tables, keeping only what the protected edges reach.
 *
 * A Lace task, because collecting has to happen in a Lace frame of its own.
 * Clearing a table resets every worker's allocation cursor for it with
 * TOGETHER, and TOGETHER only reaches all the workers if they are in the
 * collector's frame; run from an ordinary task, some workers keep a cursor
 * into a region the collector has just declared free, two of them then
 * allocate from the same region, and the table quietly corrupts.
 *
 * Use RUN(limdd_gc) off a worker thread, CALL(limdd_gc) from inside a task.
 */
VOID_TASK_DECL_0(limdd_gc);

#ifdef __cplusplus
}
#endif

#endif // QSYLVAN_LIMDD_GC_H
