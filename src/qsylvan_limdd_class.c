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
#include <stdio.h>
#include <stdlib.h>

#include <sylvan_int.h>

#include "qsylvan_limdd_class.h"
#include "qsylvan_limdd_lim.h"
#include "sylvan_edge_weights.h"
#include "edge_weight_storage/qisq2_map.h"
#include "edge_weight_storage/wgt_storage_interface.h"

/*
 * A bucket is (a, b) with
 *     a = lo | hi << 32
 *     b = var << 52 | slot << 32 | member
 * and the key is everything but the member, so hash and equality mask it out.
 */
#define MEMBER_MASK  UINT64_C(0x00000000ffffffff)
#define SLOT_SHIFT   32
#define SLOT_MAX     UINT32_C(0xfffff)
#define VAR_SHIFT    52

static llmsset_t class_table = NULL;
static bool class_requested = false;     /* read-mostly: set between sessions */

void
limdd_class_request(bool on)
{
    class_requested = on;
}

bool
limdd_class_requested(void)
{
    return class_requested;
}

static uint64_t
class_hash_cb(uint64_t a, uint64_t b, uint64_t seed)
{
    return sylvan_tabhash16(a, b & ~MEMBER_MASK, seed);
}

static int
class_equals_cb(uint64_t a, uint64_t b, uint64_t aa, uint64_t bb)
{
    return a == aa && (b & ~MEMBER_MASK) == (bb & ~MEMBER_MASK);
}

/* The data is two plain words: nothing to build or to free. */
static void class_create_cb(uint64_t *a, uint64_t *b) { (void)a; (void)b; }
static void class_destroy_cb(uint64_t a, uint64_t b) { (void)a; (void)b; }

void
limdd_class_init(size_t tablesize)
{
    if (class_table != NULL) return;
    class_table = llmsset_create(tablesize, tablesize);
    if (class_table == NULL) {
        fprintf(stderr, "sylvan: could not create the LIMDD class table\n");
        exit(1);
    }
    llmsset_set_custom(class_table, class_hash_cb, class_equals_cb,
                       class_create_cb, class_destroy_cb);
}

bool
limdd_class_on(void)
{
    return class_table != NULL;
}

void
limdd_class_quit(void)
{
    if (class_table == NULL) return;
    llmsset_free(class_table);
    class_table = NULL;
}

LIMDD_TARG
limdd_class_slot(uint32_t var, LIMDD_TARG lo, LIMDD_TARG hi, uint64_t start,
                 uint32_t j, LIMDD_TARG t, int *claimed)
{
    assert(class_table != NULL);
    assert(lo < hi);
    assert(var < (UINT32_C(1) << (64 - VAR_SHIFT)));
    if (j > SLOT_MAX) {
        fprintf(stderr, "sylvan: more than %u LIMDD nodes with the same two children\n",
                (unsigned)SLOT_MAX + 1);
        exit(1);
    }
    const uint32_t slot = (uint32_t)(start + j) & SLOT_MAX;
    const uint64_t a = (uint64_t)lo | (uint64_t)hi << 32;
    const uint64_t b = (uint64_t)var << VAR_SHIFT | (uint64_t)slot << SLOT_SHIFT | (uint64_t)t;
    const uint64_t idx = llmsset_lookupc(class_table, a, b, claimed);
    if (idx == 0) {
        fprintf(stderr, "sylvan: the LIMDD class table is full; raise --node-tab-size\n");
        exit(1);
    }
    if (*claimed) return t;
    const uint64_t *d = (const uint64_t *)llmsset_index_to_ptr(class_table, idx);
    return (LIMDD_TARG)(d[1] & MEMBER_MASK);
}

/*
 * The scalar part of the key. On exact weights a member's high scalar is
 * shared by its class up to a unit: a LIM between two nodes with the same two
 * children in the same order multiplies it by +-1 or +-i (qsylvan_limdd_class.h).
 * Multiplying a + b sqrt2 + i (c + d sqrt2) by i gives -(c + d sqrt2) + i (a +
 * b sqrt2), so the four units permute the pairs (a, b) and (c, d) and negate
 * them; a hash of each pair that ignores a joint negation, combined so that
 * the order of the pairs does not matter, is the same for all four. Read from
 * the stored, reduced rationals, so it costs no arithmetic. Floats round, and
 * an invariant that rounding can break would part two members of a class: they
 * get 0, the key they had before.
 */
static inline uint64_t
mix64(uint64_t h)
{
    h ^= h >> 33; h *= UINT64_C(0xff51afd7ed558ccd);
    h ^= h >> 33; h *= UINT64_C(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;
    return h;
}

static uint64_t
mpz_abs_hash(mpz_srcptr z, uint64_t h)
{
    const size_t n = mpz_size(z);
    const mp_limb_t *l = mpz_limbs_read(z);
    h = mix64(h ^ (uint64_t)n);
    for (size_t i = 0; i < n; i++) h = mix64(h ^ (uint64_t)l[i]);
    return h;
}

/* The same for (x, y) and (-x, -y). */
static uint64_t
pair_hash(mpq_srcptr x, mpq_srcptr y)
{
    uint64_t h = (uint64_t)(mpq_sgn(x) * mpq_sgn(y) + 2);
    h = mpz_abs_hash(mpq_numref(x), h);
    h = mpz_abs_hash(mpq_denref(x), h);
    h = mpz_abs_hash(mpq_numref(y), h);
    return mpz_abs_hash(mpq_denref(y), h);
}

uint64_t
limdd_class_scalar_key(LIMDD_LIM high)
{
    if (sylvan_get_edge_weight_type() != WGT_QISQ2) return 0;
    const qisq2_t *q = (const qisq2_t *)wgt_store_get(wgt_storage, limdd_lim_weight(high));
    const uint64_t f = pair_hash(q->a, q->b), g = pair_hash(q->c, q->d);
    return mix64((f < g ? f : g) * UINT64_C(0x9e3779b97f4a7c15) ^ (f < g ? g : f));
}

/*
 * Every marked node with two different live children whose class goes through
 * the table (limdd_class_wanted) goes back in, in index
 * order, each in the first free slot of its key. No comparison is needed: the
 * marked nodes are all members, so no two of them share a class. Serial, like
 * the purge of the stabiliser cache beside it, and both are a single pass over
 * the node table during a stop-the-world collection.
 */
void
limdd_class_gc_rebuild(void)
{
    if (class_table == NULL) return;
    llmsset_clear(class_table);
    const size_t n = limdd_node_table_size();
    for (uint64_t p = 2; p < n; p++) {
        if (!limdd_gc_node_is_marked((LIMDD_TARG)p)) continue;
        const LIMDD low = limdd_node_low((LIMDD_TARG)p);
        const LIMDD high = limdd_node_high((LIMDD_TARG)p);
        if (limdd_edge_is_zero(low) || limdd_edge_is_zero(high)) continue;
        const LIMDD_TARG v0 = limdd_target(low), v1 = limdd_target(high);
        if (v0 == v1) continue;
        const uint32_t var = limdd_node_var((LIMDD_TARG)p);
        if (!limdd_class_wanted(var, v0, v1)) continue;   /* stored canonical */
        const LIMDD_TARG lo = v0 < v1 ? v0 : v1, hi = v0 < v1 ? v1 : v0;
        /* On exact weights members are stored with their children in index
         * order (makeedge_lazy), which the scalar key assumes; on floats the
         * key is 0 and the order does not matter. */
        const uint64_t start = limdd_class_scalar_key(limdd_label(high));
        assert(v0 < v1 || start == 0);
        for (uint32_t j = 0;; j++) {
            int claimed;
            limdd_class_slot(var, lo, hi, start, j, (LIMDD_TARG)p, &claimed);
            if (claimed) break;
        }
    }
}
