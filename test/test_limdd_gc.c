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
 * Tests for limdd_gc.
 *
 * A collector can be wrong in three ways, and each needs its own evidence:
 *
 *   it sweeps something live      -- the protected diagrams are read back
 *                                    amplitude by amplitude afterwards, and
 *                                    the whole run is repeated with fresh
 *                                    allocation in between, so a bucket handed
 *                                    out twice shows up as a changed state.
 *
 *   it keeps everything           -- the tables are counted before and after,
 *                                    and unreachable entries must actually go.
 *
 *   it leaves the tables unusable -- after collecting, rebuilding a protected
 *                                    diagram must return the SAME edge. That
 *                                    only happens if the rehash put every
 *                                    surviving entry back where lookup can
 *                                    find it; a broken rehash silently
 *                                    duplicates instead of failing.
 *
 * The fourth hazard is subtler and specific to this design: a swept bucket is
 * handed out again, and the stabiliser cache is indexed by bucket. If the
 * cache is not purged, the new occupant inherits the old one's group -- a
 * wrong group, on a node that looks fine. So new diagrams built after a
 * collection have every node's group checked against brute force.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lace.h>
#include <sylvan.h>

#include "qsylvan_limdd_gc.h"
#include "sylvan_edge_weights_complex.h"
#include "test_assert.h"

#define NQUBITS 4
#define NBASIS  (1u << NQUBITS)

static uint64_t rng_state = UINT64_C(0x9E3779B97F4A7C15);
static uint64_t
rnd(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

typedef struct { double re, im; } cx;

static cx
cx_mul(cx a, cx b)
{
    cx r = { a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re };
    return r;
}

static bool
cx_eq(cx a, cx b)
{
    return fabs(a.re - b.re) < 1e-9 && fabs(a.im - b.im) < 1e-9;
}

static cx
eval_at(LIMDD e, uint64_t b, uint32_t level)
{
    static const cx ipow[4] = { {1,0}, {0,1}, {-1,0}, {0,-1} };
    if (limdd_edge_is_zero(e)) { cx z = {0,0}; return z; }

    const LIMDD_LIM lim = limdd_label(e);
    const limdd_pauli_t p = limdd_lim_pauli(lim);
    const uint64_t c = b ^ p.x;
    const unsigned k = ((unsigned)__builtin_popcountll(p.x & p.z)
                        + 2u * (unsigned)__builtin_popcountll(p.z & c)) & 3u;
    complex_t w;
    weight_value(limdd_lim_weight(lim), &w);
    const cx acc = cx_mul((cx){ w.r, w.i }, ipow[k]);

    const LIMDD_TARG t = limdd_target(e);
    if (t == LIMDD_TERMINAL) return acc;
    /* The edge may skip levels down to the target's own; it may never climb.
     * A climbing edge is malformed, and NaN makes any comparison fail. */
    const uint32_t var = limdd_node_var(t);
    if (var < level) { cx bad = {NAN, NAN}; return bad; }
    const LIMDD child = ((c >> var) & 1) ? limdd_node_high(t) : limdd_node_low(t);
    return cx_mul(acc, eval_at(child, c, var + 1));
}

static void
vector_of(LIMDD e, uint32_t level, cx *out)
{
    const unsigned n = 1u << (NQUBITS - level);
    for (unsigned j = 0; j < n; j++) out[j] = eval_at(e, (uint64_t)j << level, level);
}

/* --- building --------------------------------------------------------- */

static limdd_pauli_t
random_word_above(uint32_t var)
{
    const uint64_t all = (UINT64_C(1) << NQUBITS) - 1;
    const uint64_t above = (var + 1 >= 64) ? 0 : ~((UINT64_C(1) << (var + 1)) - 1);
    limdd_pauli_t p;
    p.x = rnd() & all & above;
    p.z = rnd() & all & above;
    return p;
}

static EVBDD_WGT
random_scalar(void)
{
    static const double re[5] = { 1.0, 0.0, -1.0,  0.0, 0.5 };
    static const double im[5] = { 0.0, 1.0,  0.0, -1.0, 0.0 };
    const unsigned w = rnd() % 5;
    return complex_lookup(re[w], im[w]);
}

static LIMDD
random_edge(uint32_t var)
{
    if (var == NQUBITS) return limdd_one_edge();

    LIMDD lo = random_edge(var + 1);
    LIMDD hi;
    const unsigned dice = rnd() % 8;
    if (dice < 3) hi = lo;
    else if (dice < 6 && !limdd_edge_is_zero(lo))
        hi = limdd_bundle(limdd_lim_mul(limdd_lim_make(random_word_above(var),
                                                       random_scalar()),
                                        limdd_label(lo)), limdd_target(lo));
    else hi = random_edge(var + 1);

    if ((rnd() % 12) == 0) hi = limdd_zero_edge();
    return limdd_makeedge(var, lo, hi);
}

/* --- brute force groups, for the recycled-bucket check ------------------ */

static void
apply_pauli(const cx *in, cx *out, uint64_t x, uint64_t z, bool neg, unsigned nb)
{
    static const cx ipow[4] = { {1,0}, {0,1}, {-1,0}, {0,-1} };
    const unsigned ny = (unsigned)__builtin_popcountll(x & z);
    for (unsigned b = 0; b < nb; b++) {
        const unsigned c = b ^ (unsigned)x;
        const unsigned k = (ny + 2u * (unsigned)__builtin_popcountll(z & c)) & 3u;
        cx r = cx_mul(in[c], ipow[k]);
        if (neg) { r.re = -r.re; r.im = -r.im; }
        out[b] = r;
    }
}

static int
cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

static int
check_stab_below(LIMDD_TARG t)
{
    if (t == LIMDD_TERMINAL) return 0;

    const uint32_t var = limdd_node_var(t);
    const LIMDD low = limdd_node_low(t), high = limdd_node_high(t);
    if (!limdd_edge_is_zero(low) && check_stab_below(limdd_target(low))) return 1;
    if (!limdd_edge_is_zero(high) && check_stab_below(limdd_target(high))) return 1;

    const unsigned nb = 1u << (NQUBITS - var);
    cx vec[NBASIS], img[NBASIS];
    vector_of(limdd_bundle(LIMDD_LIM_IDENTITY, t), var, vec);

    uint64_t want[2 * NBASIS * NBASIS];
    size_t nw = 0;
    for (uint64_t x = 0; x < nb; x++)
        for (uint64_t z = 0; z < nb; z++)
            for (int neg = 0; neg < 2; neg++) {
                apply_pauli(vec, img, x, z, neg, nb);
                bool fixes = true;
                for (unsigned b = 0; b < nb && fixes; b++)
                    if (!cx_eq(img[b], vec[b])) fixes = false;
                if (fixes) want[nw++] = (x << 33) | (z << 1) | (uint64_t)neg;
            }
    qsort(want, nw, sizeof(uint64_t), cmp_u64);

    const LIMDD_STAB s = limdd_node_stab(t);
    const size_t ng = (size_t)1 << limdd_stab_ngens(s);
    uint64_t got[1u << NQUBITS];
    for (size_t i = 0; i < ng; i++) {
        const LIMDD_LIM e = limdd_stab_element(s, i);
        const limdd_pauli_t p = limdd_lim_pauli(e);
        got[i] = ((p.x >> var) << 33) | ((p.z >> var) << 1)
               | (uint64_t)(limdd_lim_weight(e) != EVBDD_ONE);
    }
    qsort(got, ng, sizeof(uint64_t), cmp_u64);

    if (nw != ng || memcmp(want, got, sizeof(uint64_t) * nw) != 0) {
        fprintf(stderr, "node %llu at qubit %u: cached group has %zu elements, "
                        "the state is fixed by %zu -- a recycled bucket kept "
                        "its predecessor's group?\n",
                (unsigned long long)t, var, ng, nw);
        return 1;
    }
    return 0;
}

/* --- tests -------------------------------------------------------------- */

#define NKEEP 12

static LIMDD kept[NKEEP];
static cx    kept_vec[NKEEP][NBASIS];

static int
check_kept(const char *when)
{
    for (int i = 0; i < NKEEP; i++) {
        cx now[NBASIS];
        vector_of(kept[i], 0, now);
        for (unsigned b = 0; b < NBASIS; b++) {
            if (!cx_eq(now[b], kept_vec[i][b])) {
                fprintf(stderr, "%s: protected diagram %d changed at index %u: "
                                "(%g,%g) was (%g,%g)\n", when, i, b,
                        now[b].re, now[b].im, kept_vec[i][b].re, kept_vec[i][b].im);
                return 1;
            }
        }
    }
    return 0;
}

TASK_0(int, runtests)
{
    limdd_nodes_init(NQUBITS, 1LL << 18, 1LL << 18, 1LL << 18, 1LL << 18);

    /* Diagrams to keep, plus a great deal of garbage. */
    for (int i = 0; i < NKEEP; i++) {
        do { kept[i] = random_edge(0); } while (limdd_edge_is_zero(kept[i]));
        limdd_protect(&kept[i]);
        vector_of(kept[i], 0, kept_vec[i]);
    }
    for (int i = 0; i < 600; i++) (void)random_edge(0);

    const size_t n_before = limdd_node_table_count();
    const size_t l_before = limdd_lim_table_count();
    const size_t s_before = limdd_stab_table_count();
    test_assert(limdd_count_protected() == NKEEP);

    /* Remember where the kept diagrams live, to prove indices are stable. */
    LIMDD kept_edges[NKEEP];
    memcpy(kept_edges, kept, sizeof(kept));

    CALL(limdd_gc);

    if (check_kept("straight after gc")) return 1;
    printf("protected diagrams survive collection:       ok\n");

    for (int i = 0; i < NKEEP; i++) test_assert(kept[i] == kept_edges[i]);
    printf("indices are unchanged by collection:         ok\n");

    const size_t n_after = limdd_node_table_count();
    printf("nodes %zu -> %zu, LIMs %zu -> %zu, cells %zu -> %zu\n",
           n_before, n_after, l_before, limdd_lim_table_count(),
           s_before, limdd_stab_table_count());
    if (n_after >= n_before) {
        fprintf(stderr, "collection freed no nodes -- everything was kept\n");
        return 1;
    }
    printf("unreachable entries are actually freed:      ok\n");

    /*
     * Rebuilding must find the survivors rather than make copies. A rehash
     * that lost entries does not fail, it silently duplicates, so the node
     * count is watched too.
     */
    const size_t n_rebuild = limdd_node_table_count();
    for (int i = 0; i < NKEEP; i++) {
        const LIMDD_TARG t = limdd_target(kept[i]);
        if (t == LIMDD_TERMINAL) continue;
        const LIMDD again = limdd_makeedge(limdd_node_var(t),
                                           limdd_node_low(t), limdd_node_high(t));
        if (limdd_target(again) != t) {
            fprintf(stderr, "kept %d: node %llu at qubit %u rebuilt as %llu\n",
                    i, (unsigned long long)t, limdd_node_var(t),
                    (unsigned long long)limdd_target(again));
            fprintf(stderr, "  low  -> %llu  high = ", (unsigned long long)limdd_target(limdd_node_low(t)));
            limdd_lim_fprint(stderr, limdd_label(limdd_node_high(t)));
            fprintf(stderr, " -> %llu\n", (unsigned long long)limdd_target(limdd_node_high(t)));
            LIMDD_TARG t2 = limdd_target(again);
            fprintf(stderr, "  new: low -> %llu  high = ", (unsigned long long)limdd_target(limdd_node_low(t2)));
            limdd_lim_fprint(stderr, limdd_label(limdd_node_high(t2)));
            fprintf(stderr, " -> %llu\n", (unsigned long long)limdd_target(limdd_node_high(t2)));
            fprintf(stderr, "  stab(low child):\n");
            limdd_stab_fprint(stderr, limdd_node_stab(limdd_target(limdd_node_low(t))), "    ");
            return 1;
        }
    }
    /*
     * And rebuilt from the swapped presentation too. This is where choosing
     * between the two candidate nodes by interned LIM index went wrong: an
     * index says when a label was first seen, and a collection that sweeps one
     * and re-interns it elsewhere flips the comparison, so the node comes back
     * swapped.
     */
    for (int i = 0; i < NKEEP; i++) {
        const LIMDD_TARG t = limdd_target(kept[i]);
        if (t == LIMDD_TERMINAL) continue;
        const LIMDD swapped = limdd_makeedge(limdd_node_var(t),
                                             limdd_node_high(t), limdd_node_low(t));
        if (limdd_target(swapped) != t) {
            fprintf(stderr, "kept %d: node %llu rebuilds as %llu from the "
                            "swapped presentation\n", i,
                    (unsigned long long)t, (unsigned long long)limdd_target(swapped));
            return 1;
        }
    }
    test_assert(limdd_node_table_count() == n_rebuild);
    printf("rebuilding finds survivors, adds nothing:    ok\n");

    /*
     * Now refill the tables. Whatever is allocated must come out of the freed
     * buckets, and if any of them were still live the kept diagrams break.
     */
    for (int i = 0; i < 600; i++) (void)random_edge(0);
    if (check_kept("after reallocating over the freed space")) return 1;
    printf("reused buckets do not disturb what was kept: ok\n");

    /* Nodes in recycled buckets must not inherit a stale stabiliser group. */
    for (int trial = 0; trial < 25; trial++) {
        const LIMDD e = random_edge(0);
        if (limdd_edge_is_zero(e)) continue;
        if (check_stab_below(limdd_target(e))) return 1;
    }
    printf("recycled buckets get their own group:        ok\n");

    /* And it all still holds after collecting repeatedly. */
    for (int round = 0; round < 5; round++) {
        for (int i = 0; i < 200; i++) (void)random_edge(0);
        CALL(limdd_gc);
        if (check_kept("after repeated collection")) return 1;
    }
    printf("five more collections change nothing:        ok\n");

    for (int i = 0; i < NKEEP; i++) limdd_unprotect(&kept[i]);
    test_assert(limdd_count_protected() == 0);

    CALL(limdd_gc);
    printf("after unprotecting: %zu nodes remain\n", limdd_node_table_count());
    test_assert(limdd_node_table_count() == 0);
    printf("dropping the last root empties the table:    ok\n");

    limdd_nodes_quit();
    return 0;
}

int
main(void)
{
    lace_start(8, 0);
    sylvan_set_sizes(1LL << 16, 1LL << 16, 1LL << 16, 1LL << 16);
    sylvan_init_package();
    sylvan_init_edge_weights(1LL << 18, 1LL << 18, 1e-14, WGT_COMPLEX_128, COMP_HASHMAP);

    int res = RUN(runtests);

    sylvan_edge_weights_free();
    sylvan_quit();
    lace_stop();
    return res;
}
