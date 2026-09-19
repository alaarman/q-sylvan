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

#include <assert.h>
#include <stdlib.h>

#include <sylvan_int.h>

#include "qsylvan_limdd_gc.h"
#include "sylvan_refs.h"

static refs_table_t limdd_protected;
static int limdd_protected_created = 0;

void
limdd_protect(LIMDD *a)
{
    if (!limdd_protected_created) {
        protect_create(&limdd_protected, 4096);
        limdd_protected_created = 1;
    }
    protect_up(&limdd_protected, (size_t)a);
}

void
limdd_unprotect(LIMDD *a)
{
    if (limdd_protected.refs_table != NULL) protect_down(&limdd_protected, (size_t)a);
}

size_t
limdd_count_protected(void)
{
    if (!limdd_protected_created) return 0;
    return protect_count(&limdd_protected);
}

/**
 * Mark everything reachable from `e`.
 *
 * The node's own bucket is what stops the recursion: marking it returns 0 the
 * second time, and a diagram is a DAG, so without that check a shared
 * subdiagram would be walked once per path into it -- which for the diagrams
 * LIMDDs are good at is exponential.
 */
VOID_TASK_1(limdd_gc_mark_edge, LIMDD, e)
{
    limdd_gc_mark_lim(limdd_label(e));

    const LIMDD_TARG t = limdd_target(e);
    if (t == LIMDD_TERMINAL) return;
    if (limdd_gc_mark_node(t) == 0) return;

    /*
     * The cached group is read here, before the cache is purged, and is the
     * only reference to its list cells. A node whose group was never computed
     * has nothing to mark.
     */
    const uint64_t s = limdd_node_stab_raw(t);
    if (s != 0) limdd_gc_mark_stab(s);

    SPAWN(limdd_gc_mark_edge, limdd_node_high(t));
    CALL(limdd_gc_mark_edge, limdd_node_low(t));
    SYNC(limdd_gc_mark_edge);
}

VOID_TASK_0(limdd_gc_mark_roots)
{
    if (!limdd_protected_created) return;

    uint64_t *it = protect_iter(&limdd_protected, 0, limdd_protected.refs_size);
    while (it != NULL) {
        LIMDD *p = (LIMDD *)protect_next(&limdd_protected, &it, limdd_protected.refs_size);
        if (p != NULL) CALL(limdd_gc_mark_edge, *p);
    }
}

VOID_TASK_0(limdd_gc_go)
{
    limdd_gc_clear_nodes();
    limdd_gc_clear_stabs();
    limdd_gc_clear_lims();

    /*
     * The identity and the zero map sit in ordinary buckets and are held in
     * globals rather than by any edge, so nothing would mark them. Marking the
     * identity also covers the identity Pauli word, which every other reserved
     * value shares.
     */
    limdd_gc_mark_lim(LIMDD_LIM_IDENTITY);
    limdd_gc_mark_lim(LIMDD_LIM_ZERO);

    CALL(limdd_gc_mark_roots);

    limdd_gc_purge_stab_cache();

    limdd_gc_rehash_lims();
    limdd_gc_rehash_stabs();
    limdd_gc_rehash_nodes();
}

VOID_TASK_IMPL_0(limdd_gc)
{
    NEWFRAME(limdd_gc_go);
}
