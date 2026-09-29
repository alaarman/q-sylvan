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

/*
 * What the last collection left behind. A collection is expensive -- it
 * begins by clearing the operation cache, so the recursion loses every
 * memoised result and redoes the work, which allocates again. If the previous
 * one freed almost nothing, running another immediately is a loop that gets
 * slower each time round rather than a garbage collector.
 */
static _Atomic(uint64_t) limdd_gc_floor_lim   = 0;
static _Atomic(uint64_t) limdd_gc_floor_pauli = 0;
static _Atomic(uint64_t) limdd_gc_floor_node  = 0;
static _Atomic(uint64_t) limdd_gc_floor_stab  = 0;
static _Atomic(uint64_t) limdd_gc_floor_wgt   = 0;

/*
 * A new session: limdd_nodes_init has just made new tables, which hand out
 * the indices of the last session's again. The operation cache still holds
 * that session's entries, keyed on its indices, and a lookup in the new one
 * would return them -- in the BQD's translation and Pauli families, whose
 * keys are a label's Pauli word and a node, a wrong state and not a miss
 * (two families share their cache ids, and one family run twice meets its
 * own). So it is cleared here, for every arm, as a collection clears it. The
 * floors are what the last collection left in tables that no longer exist,
 * and a floor above a smaller new table would stop collection for good.
 */
void
limdd_gc_new_session(void)
{
    cache_clear();
    limdd_gc_floor_lim   = 0;
    limdd_gc_floor_pauli = 0;
    limdd_gc_floor_node  = 0;
    limdd_gc_floor_stab  = 0;
    limdd_gc_floor_wgt   = 0;
}

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

/* --- per-worker reference stacks ------------------------------------------
 *
 * A direct mirror of Sylvan's mtbdd_refs. Three stacks per worker: values
 * pushed by copy (r), addresses of variables re-read at collection time (p),
 * and spawned tasks whose results must survive while they run (s).
 */

typedef struct limdd_refs_task
{
    Task *t;
    void *f;
} *limdd_refs_task_t;

typedef struct limdd_refs_internal
{
    const LIMDD **pbegin, **pend, **pcur;
    LIMDD *rbegin, *rend, *rcur;
    limdd_refs_task_t sbegin, send, scur;
} *limdd_refs_internal_t;

static SYLVAN_TLS limdd_refs_internal_t limdd_refs_key;

VOID_TASK_2(limdd_refs_mark_p_par, const LIMDD**, begin, size_t, count)
{
    if (count < 32) {
        while (count) { CALL(limdd_gc_mark_edge, **(begin++)); count--; }
    } else {
        SPAWN(limdd_refs_mark_p_par, begin, count / 2);
        CALL(limdd_refs_mark_p_par, begin + (count / 2), count - count / 2);
        SYNC(limdd_refs_mark_p_par);
    }
}

VOID_TASK_2(limdd_refs_mark_r_par, LIMDD*, begin, size_t, count)
{
    if (count < 32) {
        while (count) { CALL(limdd_gc_mark_edge, *begin++); count--; }
    } else {
        SPAWN(limdd_refs_mark_r_par, begin, count / 2);
        CALL(limdd_refs_mark_r_par, begin + (count / 2), count - count / 2);
        SYNC(limdd_refs_mark_r_par);
    }
}

VOID_TASK_2(limdd_refs_mark_s_par, limdd_refs_task_t, begin, size_t, count)
{
    if (count < 32) {
        while (count > 0) {
            Task *t = begin->t;
            if (!TASK_IS_STOLEN(t)) return;
            if (t->f == begin->f && TASK_IS_COMPLETED(t)) {
                CALL(limdd_gc_mark_edge, *(LIMDD*)TASK_RESULT(t));
            }
            begin += 1;
            count -= 1;
        }
    } else {
        if (!TASK_IS_STOLEN(begin->t)) return;
        SPAWN(limdd_refs_mark_s_par, begin, count / 2);
        CALL(limdd_refs_mark_s_par, begin + (count / 2), count - count / 2);
        SYNC(limdd_refs_mark_s_par);
    }
}

VOID_TASK_0(limdd_refs_mark_task)
{
    if (limdd_refs_key == NULL) return;
    SPAWN(limdd_refs_mark_p_par, limdd_refs_key->pbegin,
          (size_t)(limdd_refs_key->pcur - limdd_refs_key->pbegin));
    SPAWN(limdd_refs_mark_r_par, limdd_refs_key->rbegin,
          (size_t)(limdd_refs_key->rcur - limdd_refs_key->rbegin));
    CALL(limdd_refs_mark_s_par, limdd_refs_key->sbegin,
         (size_t)(limdd_refs_key->scur - limdd_refs_key->sbegin));
    SYNC(limdd_refs_mark_r_par);
    SYNC(limdd_refs_mark_p_par);
}

VOID_TASK_0(limdd_refs_mark_all)
{
    TOGETHER(limdd_refs_mark_task);
}

VOID_TASK_0(limdd_refs_init_task)
{
    /* Padded: the cursors below are written by every push and pop, so two
     * workers sharing a line would invalidate each other continuously. */
    limdd_refs_internal_t s =
        (limdd_refs_internal_t)sylvan_alloc_padded(sizeof(struct limdd_refs_internal));
    if (s == NULL) {
        fprintf(stderr, "sylvan: out of memory allocating a LIMDD reference stack\n");
        exit(1);
    }
    s->pcur = s->pbegin = (const LIMDD**)malloc(sizeof(LIMDD*) * 1024);
    s->pend = s->pbegin + 1024;
    s->rcur = s->rbegin = (LIMDD*)malloc(sizeof(LIMDD) * 1024);
    s->rend = s->rbegin + 1024;
    s->scur = s->sbegin = (limdd_refs_task_t)malloc(sizeof(struct limdd_refs_task) * 1024);
    s->send = s->sbegin + 1024;
    limdd_refs_key = s;
}

void
limdd_refs_init(void)
{
    TOGETHER(limdd_refs_init_task);
}

VOID_TASK_0(limdd_refs_quit_task)
{
    limdd_refs_internal_t s = limdd_refs_key;
    if (s == NULL) return;
    free(s->pbegin);
    free(s->rbegin);
    free(s->sbegin);
    sylvan_free_padded(s);
    limdd_refs_key = NULL;
}

/*
 * Every init allocates a set per worker, and the key is thread-local, so a
 * set that is not freed here is unreachable once lace_stop ends its thread.
 * A program that opens several sessions, as test_bqd does one per label
 * family, leaked one set per worker per session. evbdd_refs_cleanup is the
 * same thing for the EVBDD.
 */
void
limdd_refs_quit(void)
{
    TOGETHER(limdd_refs_quit_task);
}

static void SYLVAN_NOINLINE
limdd_refs_ptrs_up(limdd_refs_internal_t refs)
{
    size_t cur = (size_t)(refs->pcur - refs->pbegin);
    size_t size = (size_t)(refs->pend - refs->pbegin);
    refs->pbegin = (const LIMDD**)realloc(refs->pbegin, sizeof(LIMDD*) * size * 2);
    refs->pcur = refs->pbegin + cur;
    refs->pend = refs->pbegin + (size * 2);
}

static LIMDD SYLVAN_NOINLINE
limdd_refs_refs_up(limdd_refs_internal_t refs, LIMDD res)
{
    size_t size = (size_t)(refs->rend - refs->rbegin);
    refs->rbegin = (LIMDD*)realloc(refs->rbegin, sizeof(LIMDD) * size * 2);
    refs->rcur = refs->rbegin + size;
    refs->rend = refs->rbegin + (size * 2);
    return res;
}

static void SYLVAN_NOINLINE
limdd_refs_tasks_up(limdd_refs_internal_t refs)
{
    size_t size = (size_t)(refs->send - refs->sbegin);
    refs->sbegin = (limdd_refs_task_t)realloc(refs->sbegin,
                                              sizeof(struct limdd_refs_task) * size * 2);
    refs->scur = refs->sbegin + size;
    refs->send = refs->sbegin + (size * 2);
}

LIMDD
limdd_refs_push(LIMDD e)
{
    *(limdd_refs_key->rcur++) = e;
    if (limdd_refs_key->rcur == limdd_refs_key->rend)
        return limdd_refs_refs_up(limdd_refs_key, e);
    return e;
}

void
limdd_refs_pop(long amount)
{
    limdd_refs_key->rcur -= amount;
}

void
limdd_refs_pushptr(const LIMDD *ptr)
{
    *limdd_refs_key->pcur++ = ptr;
    if (limdd_refs_key->pcur == limdd_refs_key->pend) limdd_refs_ptrs_up(limdd_refs_key);
}

void
limdd_refs_popptr(size_t amount)
{
    limdd_refs_key->pcur -= amount;
}

void
limdd_refs_spawn(Task *t)
{
    limdd_refs_key->scur->t = t;
    limdd_refs_key->scur->f = t->f;
    limdd_refs_key->scur += 1;
    if (limdd_refs_key->scur == limdd_refs_key->send) limdd_refs_tasks_up(limdd_refs_key);
}

LIMDD
limdd_refs_sync(LIMDD result)
{
    limdd_refs_key->scur -= 1;
    return result;
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
    /*
     * The operation cache first. Its entries are keyed on edges and hold
     * edges, and a collection frees buckets for reuse -- so an entry that
     * survives a sweep can hand back a node that now belongs to something
     * else. Sylvan's own collector clears it for the same reason. Nothing
     * about this is visible until a run is long enough to collect, which is
     * what let it through: every circuit small enough to finish without
     * collecting gave the right answer.
     */
    cache_clear();

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

    /* Every worker's own stack, which is what keeps a half-finished operation
     * alive across a collection that happens underneath it. */
    CALL(limdd_refs_mark_all);

    limdd_gc_purge_stab_cache();

    /*
     * Edge weights, between marking and the rehash below -- the only window
     * where the live set is known and the LIM buckets may still be rewritten.
     *
     * Gated on the weight table filling rather than done on every collection,
     * because a collection here is triggered by the LIM, node or stabiliser
     * table and says nothing about the weights, while wgt_table_gc_init_new
     * doubles the table each time it runs. Collecting when there is no
     * pressure would grow it for no reason.
     */
    const uint64_t wgt_size = sylvan_get_edge_weight_table_size();
    if (wgt_table_entries_estimate() > wgt_size / 2 &&
        getenv("LIMDD_NO_WGT_GC") == NULL) {
        const uint64_t before = wgt_table_entries_estimate();
        const size_t kept = limdd_gc_remap_weights();
        if (getenv("LIMDD_GC_VERBOSE"))
            fprintf(stderr, "[gc] weights: %llu -> %zu kept (table %llu)\n",
                    (unsigned long long)before, kept, (unsigned long long)wgt_size);
    }

    limdd_gc_rehash_lims();
    limdd_gc_rehash_stabs();
    limdd_gc_rehash_nodes();

    /* Remember what is left, so a collection that freed nothing is not
     * immediately repeated. The weight table has also just doubled if it was
     * collected, so its floor is measured against the new size. */
    limdd_gc_floor_lim   = limdd_lim_table_count();
    limdd_gc_floor_pauli = limdd_pauli_table_count();
    limdd_gc_floor_node  = limdd_node_table_count();
    limdd_gc_floor_stab  = limdd_stab_table_count();
    limdd_gc_floor_wgt   = wgt_table_entries_estimate();
}

static _Atomic(int) limdd_gc_in_progress = 0;

VOID_TASK_IMPL_0(limdd_gc)
{
    /*
     * One worker opens the frame, the rest yield into it.
     *
     * Without the compare-exchange two workers that both notice a full table
     * would each open a frame, and the second would collect a heap the first
     * had already rebuilt. The loser yields into the frame, which is how it
     * ends up helping rather than spinning.
     *
     * It waits on the flag and not on the frame, unlike sylvan_gc in
     * sylvan_common.c. A loser whose compare-exchange fails after the
     * winner's frame has ended, but before the winner has reset the flag --
     * one that took part in that frame from a steal loop and then asks for a
     * collection of its own -- would otherwise wait for a frame that never
     * comes. It also yields into any other frame that opens meanwhile, since
     * a worker that spins without yielding holds up every NEWFRAME and
     * TOGETHER.
     */
    int zero = 0;
    if (atomic_compare_exchange_strong(&limdd_gc_in_progress, &zero, 1)) {
        NEWFRAME(limdd_gc_go);
        limdd_gc_in_progress = 0;
    } else {
        while (atomic_load_explicit(&limdd_gc_in_progress, memory_order_acquire) != 0) {
            YIELD_NEWFRAME();
        }
    }
}

/**
 * Is a table close enough to full that the caller should collect?
 *
 * This is the whole trigger policy, in one place. It used to be duplicated in
 * the simulator's gate loop, which is what let the damping below sit dead: the
 * loop tested its own undamped thresholds and never consulted this.
 *
 * It cannot be asked from inside an operation, which is where a full table is
 * actually hit -- one gate on a large diagram mints millions of intermediates.
 * Collecting mid-flight needs the operation cache to survive a sweep, which it
 * does not; see the note on 736bc97. So this is a between-gate check, and the
 * headroom below is what has to absorb a whole gate.
 */
/** Wanted when a table is nearly full AND the last collection left room. */
static bool
table_wants_gc(uint64_t count, uint64_t size, uint64_t floor)
{
    if (size == 0) return false;
    /*
     * A quarter free, not an eighth. The headroom has to absorb a whole gate,
     * because nothing can collect in the middle of one, and one gate on a
     * large diagram mints millions of weights. Measured: at an eighth,
     * rand_n20_d700 dies with a full weight table at 2^24, where a quarter
     * finishes it in 50s.
     */
    if (count <= size - (size >> 2)) return false;      /* under 3/4 full */
    /* Re-collect only once a quarter of what the last one freed is used up
     * again. When it freed nothing the condition never holds, and the caller
     * runs out of table and reports it rather than spinning. */
    return count > floor + ((size - floor) >> 2);
}

bool
limdd_gc_wanted(void)
{
    /*
     * The weight table keeps a plain counter, so it is asked every time. The
     * other four walk an occupancy bitmap to count, which cost 17% of a run
     * on a 20-qubit Clifford+T circuit when asked every gate -- more than the
     * coset search -- so they are asked every 16th time. The quarter of
     * headroom in table_wants_gc is what absorbs those 16 gates.
     */
    if (table_wants_gc(wgt_table_entries_estimate(),
                       sylvan_get_edge_weight_table_size(), limdd_gc_floor_wgt))
        return true;

    static _Atomic(uint64_t) calls = 0;
    if ((atomic_fetch_add(&calls, 1) & 15u) != 0) return false;

    const uint64_t lims = limdd_lim_table_count();
    if (table_wants_gc(lims, limdd_lim_table_size(), limdd_gc_floor_lim)) return true;

    /*
     * The Pauli words too. A LIMDD and the BQD's scalar family intern few of
     * them, but in the BQD's Pauli family every label product can make a new
     * word, at about a third of the rate of new LIMs and far above that of
     * new nodes, and a Pauli table that fills unasked stops the run with the
     * other tables nearly empty. Every word the table holds is the word of a
     * LIM, since only limdd_lim_make interns one, so the table can be three
     * quarters full only when there are as many LIMs, and until then its
     * count, another walk of a bitmap, is not taken.
     */
    const uint64_t pauli_size = limdd_pauli_table_size();
    if (lims > pauli_size - (pauli_size >> 2) &&
        table_wants_gc(limdd_pauli_table_count(), pauli_size, limdd_gc_floor_pauli))
        return true;

    return table_wants_gc(limdd_node_table_count(), limdd_node_table_size(),
                          limdd_gc_floor_node)
        || table_wants_gc(limdd_stab_table_count(), limdd_stab_table_size(),
                          limdd_gc_floor_stab);
}
