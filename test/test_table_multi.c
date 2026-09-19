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
 * Several llmsset tables live at the same time.
 *
 * This used to be forbidden: every table shared one thread-local region
 * cursor, so a worker inserting into one table would move the cursor another
 * table was using, and buckets would be handed out twice. Each table now
 * claims an id and the cursor is per (worker, table).
 *
 * The interesting case is therefore concurrent insertion into DIFFERENT tables
 * from the same workers, which is what the LIMDD work needs: separate tables
 * for nodes, Pauli words and LIM labels.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lace.h>

#include "sylvan_platform.h"
#include "sylvan_table.h"
#include "test_assert.h"

#define NTABLES  3
#define NKEYS    4000

static llmsset_t tables[NTABLES];
static uint64_t  indices[NTABLES][NKEYS];

/**
 * Insert keys [lo,hi) into ALL tables, interleaved.
 *
 * The interleaving is the whole point. When one worker inserts into table A
 * then table B then table C in a tight loop, a single shared region cursor
 * gets carried from one table into the next: the worker would go on using
 * region R in table B without ever having claimed R in B's ownership bitmap,
 * so a second worker can claim that same region and both then allocate from
 * the same bitmap2 word. That word is updated with a plain load followed by
 * an atomic_fetch_or rather than a CAS, so both workers can pick the same
 * free bit and walk away with the same bucket index.
 *
 * A test where each task touches only one table does NOT catch this.
 */
VOID_TASK_3(fill_range, size_t, lo, size_t, hi, int, expect_created)
{
    if (hi - lo > 64) {
        size_t mid = lo + (hi - lo) / 2;
        SPAWN(fill_range, lo, mid, expect_created);
        CALL(fill_range, mid, hi, expect_created);
        SYNC(fill_range);
        return;
    }
    for (size_t i = lo; i < hi; i++) {
        for (int t = 0; t < NTABLES; t++) {
            int created = -1;
            uint64_t idx = llmsset_lookup(tables[t], (uint64_t)i, (uint64_t)t, &created);
            if (idx == 0) {
                fprintf(stderr, "table %d: insertion of key %zu failed (table full?)\n", t, i);
                exit(1);
            }
            if (expect_created >= 0 && created != expect_created) {
                fprintf(stderr, "table %d key %zu: created=%d, expected %d\n",
                        t, i, created, expect_created);
                exit(1);
            }
            indices[t][i] = idx;
        }
    }
}

VOID_TASK_1(fill_all, int, expect_created)
{
    CALL(fill_range, 0, NKEYS, expect_created);
}

static int
cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

int
check_distinct_within_table(void)
{
    for (int t = 0; t < NTABLES; t++) {
        uint64_t *sorted = malloc(sizeof(uint64_t) * NKEYS);
        test_assert(sorted != NULL);
        memcpy(sorted, indices[t], sizeof(uint64_t) * NKEYS);
        qsort(sorted, NKEYS, sizeof(uint64_t), cmp_u64);
        for (size_t i = 1; i < NKEYS; i++) {
            /* Two distinct keys sharing a bucket is the exact corruption a
             * shared region cursor used to cause. */
            if (sorted[i] == sorted[i - 1]) {
                fprintf(stderr, "table %d: bucket %llu handed out twice\n",
                        t, (unsigned long long)sorted[i]);
                free(sorted);
                return 1;
            }
        }
        free(sorted);
    }
    return 0;
}

/**
 * Every occupied bucket must lie in a region this table actually owns.
 *
 * This is the invariant a shared region cursor breaks, and checking it
 * directly is deterministic. Trying instead to catch the resulting race --
 * two workers allocating from one bitmap2 word with a load followed by
 * atomic_fetch_or, rather than a CAS, and picking the same free bit -- only
 * fires when regions are exhausted and cursors converge, so it makes for a
 * flaky test. The ownership check fires on the very first insertion that a
 * carried-over cursor performs.
 */
int
check_region_ownership(void)
{
    for (int t = 0; t < NTABLES; t++) {
        llmsset_t dbs = tables[t];
        /* buckets 0 and 1 are marked occupied at creation without any region
         * being claimed, so they are not evidence of anything */
        for (uint64_t b = 2; b < dbs->table_size; b++) {
            const uint64_t occupied_word = dbs->bitmap2[b / 64];
            const uint64_t occupied_mask = UINT64_C(0x8000000000000000) >> (b & 63);
            if (!(occupied_word & occupied_mask)) continue;

            const uint64_t region = b / 512;
            const uint64_t owned_word = dbs->bitmap1[region / 64];
            const uint64_t owned_mask = UINT64_C(0x8000000000000000) >> (region & 63);
            if (!(owned_word & owned_mask)) {
                fprintf(stderr,
                        "table %d: bucket %llu is in use but its region %llu "
                        "was never claimed -- the region cursor leaked in from "
                        "another table\n",
                        t, (unsigned long long)b, (unsigned long long)region);
                return 1;
            }
        }
    }
    return 0;
}

int
check_contents(void)
{
    /* Every key must still be there, with the same index and its own data. */
    for (int t = 0; t < NTABLES; t++) {
        for (size_t i = 0; i < NKEYS; i++) {
            uint64_t *data = (uint64_t *)llmsset_index_to_ptr(tables[t], indices[t][i]);
            test_assert(data[0] == (uint64_t)i);
            test_assert(data[1] == (uint64_t)t);
        }
    }
    return 0;
}

TASK_0(int, runtests)
{
    /* First pass: everything is new. */
    RUN(fill_all, 1);
    if (check_distinct_within_table()) return 1;
    if (check_region_ownership()) return 1;
    if (check_contents()) return 1;
    printf("multi-table insert (3 tables, %d keys each):   ok\n", NKEYS);
    printf("every bucket lies in a region its table owns:  ok\n");

    /* Second pass: nothing is new, and every key must come back to the same
     * bucket it got the first time. */
    uint64_t (*first)[NKEYS] = malloc(sizeof(indices));
    test_assert(first != NULL);
    memcpy(first, indices, sizeof(indices));

    RUN(fill_all, 0);
    test_assert(memcmp(first, indices, sizeof(indices)) == 0);
    free(first);
    printf("multi-table lookup is idempotent:             ok\n");

    if (check_contents()) return 1;
    printf("tables did not corrupt each other:            ok\n");

    return 0;
}

/**
 * Table ids must be returned on free, otherwise a long-running program that
 * creates and destroys tables exhausts LLMSSET_MAX_TABLES. Cycle more times
 * than there are ids to prove they are reused.
 *
 * This runs as a Lace task because llmsset_lookup reaches lace_get_worker()
 * when it claims a region, which is NULL on the thread that called
 * lace_start().
 */
TASK_0(int, test_id_reuse)
{
    for (int round = 0; round < LLMSSET_MAX_TABLES * 3; round++) {
        /* >= 4096 buckets: bitmap1 is max_size/4096 bytes and must be non-empty */
        llmsset_t tmp = llmsset_create(1LL << 13, 1LL << 13);
        test_assert(tmp != NULL);
        int created = -1;
        uint64_t idx = llmsset_lookup(tmp, 42, (uint64_t)round, &created);
        test_assert(idx != 0);
        test_assert(created == 1);
        llmsset_free(tmp);
    }
    printf("table ids are released and reused:            ok\n");
    return 0;
}

int
main(void)
{
    lace_start(8, 0);

    for (int t = 0; t < NTABLES; t++) {
        tables[t] = llmsset_create(1LL << 14, 1LL << 14);
        if (tables[t] == NULL) { fprintf(stderr, "llmsset_create failed\n"); exit(1); }
    }

    int res = RUN(runtests);

    for (int t = 0; t < NTABLES; t++) llmsset_free(tables[t]);

    if (res == 0) res = RUN(test_id_reuse);

    lace_stop();
    return res;
}
