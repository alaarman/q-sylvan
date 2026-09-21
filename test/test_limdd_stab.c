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
 * Tests for stabiliser generator sets.
 *
 * The property everything else rests on is that a handle identifies a GROUP,
 * not a generating set: the canonical form compares stabiliser groups by
 * comparing handles, so two generating sets that span the same group must
 * intern to the same handle, and two that do not must not.
 *
 * Group equality is therefore checked here the only way that cannot itself be
 * wrong: by enumerating all 2^k elements of both groups and comparing the sets
 * exactly. That is exponential, so the tests stay at small k -- but it means a
 * bug in the reduction cannot hide behind the same bug in the comparison.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lace.h>
#include <sylvan.h>

#include "qsylvan_limdd_stab.h"
#include "sylvan_edge_weights_complex.h"
#include "test_assert.h"

#define NQUBITS 5

static uint64_t rng_state = UINT64_C(0xFEEDFACECAFEBEEF);
static uint64_t
rnd(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

static LIMDD_LIM
gen(const char *s, bool neg)
{
    limdd_pauli_t p;
    if (!limdd_pauli_from_string(s, NQUBITS, &p)) {
        fprintf(stderr, "bad Pauli string %s\n", s);
        exit(1);
    }
    return limdd_lim_make(p, neg ? EVBDD_MIN_ONE : EVBDD_ONE);
}

/* --- enumerating a group, as an independent check ------------------------ */

static int
cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

/** All 2^k elements of `s`, sorted. Caller frees. */
static uint64_t *
elements(LIMDD_STAB s, size_t *count)
{
    const size_t k = limdd_stab_ngens(s);
    const size_t n = (size_t)1 << k;
    uint64_t *out = malloc(sizeof(uint64_t) * n);
    if (out == NULL) { fprintf(stderr, "oom\n"); exit(1); }
    for (size_t i = 0; i < n; i++) out[i] = limdd_stab_element(s, i);
    qsort(out, n, sizeof(uint64_t), cmp_u64);
    *count = n;
    return out;
}

/** True iff the two handles denote the same set of group elements. */
static bool
same_group(LIMDD_STAB a, LIMDD_STAB b)
{
    size_t na, nb;
    uint64_t *ea = elements(a, &na);
    uint64_t *eb = elements(b, &nb);
    const bool eq = (na == nb) && memcmp(ea, eb, sizeof(uint64_t) * na) == 0;
    free(ea);
    free(eb);
    return eq;
}

/* --- tests --------------------------------------------------------------- */

int
test_trivial(void)
{
    test_assert(limdd_stab_is_trivial(limdd_stab_make(NULL, 0)));
    test_assert(limdd_stab_ngens(LIMDD_STAB_TRIVIAL) == 0);

    /* The identity generates nothing, however many times it is offered. */
    LIMDD_LIM ids[3] = { LIMDD_LIM_IDENTITY, LIMDD_LIM_IDENTITY, LIMDD_LIM_IDENTITY };
    test_assert(limdd_stab_is_trivial(limdd_stab_make(ids, 3)));

    test_assert(limdd_stab_contains(LIMDD_STAB_TRIVIAL, LIMDD_LIM_IDENTITY));
    test_assert(!limdd_stab_contains(LIMDD_STAB_TRIVIAL, gen("ZIIII", false)));
    return 0;
}

int
test_known_groups(void)
{
    /*
     * |00000> is stabilised by Z on every qubit, independently. Five
     * generators, 32 elements, and every one of them must be a product of Zs
     * with a plus sign.
     */
    LIMDD_LIM z[NQUBITS];
    z[0] = gen("ZIIII", false); z[1] = gen("IZIII", false);
    z[2] = gen("IIZII", false); z[3] = gen("IIIZI", false);
    z[4] = gen("IIIIZ", false);
    LIMDD_STAB all_zero = limdd_stab_make(z, NQUBITS);
    test_assert(limdd_stab_ngens(all_zero) == NQUBITS);

    for (uint64_t k = 0; k < 32; k++) {
        LIMDD_LIM e = limdd_stab_element(all_zero, k);
        test_assert(limdd_stab_contains(all_zero, e));
        test_assert(limdd_lim_weight(e) == EVBDD_ONE);
    }
    /* X flips |0>, so it stabilises nothing here; nor does -Z. */
    test_assert(!limdd_stab_contains(all_zero, gen("XIIII", false)));
    test_assert(!limdd_stab_contains(all_zero, gen("ZIIII", true)));

    /*
     * The Bell state on qubits 0,1 (times |000>) is stabilised by XX and ZZ but
     * by neither X nor Z alone -- the case a per-qubit representation cannot
     * express, and the reason LIMDDs need groups rather than single operators.
     */
    LIMDD_LIM bell[2] = { gen("XXIII", false), gen("ZZIII", false) };
    LIMDD_STAB b = limdd_stab_make(bell, 2);
    test_assert(limdd_stab_ngens(b) == 2);
    test_assert(limdd_stab_contains(b, gen("XXIII", false)));
    test_assert(limdd_stab_contains(b, gen("ZZIII", false)));
    test_assert(!limdd_stab_contains(b, gen("XIIII", false)));
    test_assert(!limdd_stab_contains(b, gen("ZIIII", false)));

    /*
     * XX * ZZ = -YY, not +YY. Getting this wrong is invisible to a membership
     * test that ignores signs, and would silently merge two different states.
     */
    test_assert(limdd_stab_contains(b, gen("YYIII", true)));
    test_assert(!limdd_stab_contains(b, gen("YYIII", false)));
    return 0;
}

int
test_generating_set_invariance(void)
{
    /*
     * The heart of the matter: a group must get ONE handle no matter which
     * generating set produced it. Handles are compared with ==, so a group that
     * can land on two handles would make two identical nodes look different and
     * break canonicity outright.
     */
    LIMDD_LIM a = gen("XXIII", false);
    LIMDD_LIM b = gen("ZZIII", false);
    LIMDD_LIM ab = limdd_lim_mul(a, b);          /* -YYIII */

    LIMDD_LIM s1[2] = { a, b };
    LIMDD_LIM s2[2] = { b, a };                   /* reordered */
    LIMDD_LIM s3[2] = { a, ab };                  /* different basis */
    LIMDD_LIM s4[3] = { a, b, ab };               /* redundant */
    LIMDD_LIM s5[4] = { b, ab, a, LIMDD_LIM_IDENTITY };

    LIMDD_STAB h = limdd_stab_make(s1, 2);
    test_assert(limdd_stab_make(s2, 2) == h);
    test_assert(limdd_stab_make(s3, 2) == h);
    test_assert(limdd_stab_make(s4, 3) == h);
    test_assert(limdd_stab_make(s5, 4) == h);

    /* And the handles really do denote the same elements, checked the slow way. */
    test_assert(same_group(limdd_stab_make(s3, 2), h));

    /* A genuinely different group must not collide. */
    LIMDD_LIM other[2] = { gen("XXIII", false), gen("ZZIII", true) };
    test_assert(limdd_stab_make(other, 2) != h);
    test_assert(!same_group(limdd_stab_make(other, 2), h));
    return 0;
}

/**
 * Random groups, built from random commuting generators, then re-derived from
 * a random alternative basis. Reordering alone would not exercise much; taking
 * random products does, because it forces the reduction to recover the same
 * echelon form from bases that share no rows.
 */
int
test_random_invariance(void)
{
    for (int trial = 0; trial < 2000; trial++) {
        /* Build up to 4 pairwise-commuting generators by rejection. */
        LIMDD_LIM g[4];
        size_t k = 0;
        for (int attempt = 0; attempt < 64 && k < 4; attempt++) {
            const uint64_t mask = (UINT64_C(1) << NQUBITS) - 1;
            limdd_pauli_t p = limdd_pauli_identity();
            {
                const uint64_t xw = rnd() & mask, zw = rnd() & mask;
                for (size_t qq = 0; qq < NQUBITS; qq++) {
                    const unsigned op = 2u * ((unsigned)(xw >> qq) & 1u)
                                      +      ((unsigned)(zw >> qq) & 1u);
                    if (op) limdd_pauli_set(&p, qq, (limdd_pauli_op_t)op);
                }
            }
            if (limdd_pauli_is_identity(p)) continue;

            bool commutes = true;
            for (size_t i = 0; i < k; i++) {
                if (limdd_pauli_commutation_phase(p, limdd_lim_pauli(g[i])) != 0) {
                    commutes = false;
                    break;
                }
            }
            if (!commutes) continue;

            LIMDD_LIM cand = limdd_lim_make(p, (rnd() & 1) ? EVBDD_MIN_ONE : EVBDD_ONE);
            /* Skip anything already in the span: -I would otherwise appear. */
            LIMDD_STAB sofar = limdd_stab_make(g, k);
            if (limdd_stab_contains(sofar, cand)) continue;
            limdd_pauli_t np = p;
            LIMDD_LIM negcand = limdd_lim_make(np, (limdd_lim_weight(cand) == EVBDD_ONE)
                                                   ? EVBDD_MIN_ONE : EVBDD_ONE);
            if (limdd_stab_contains(sofar, negcand)) continue;

            g[k++] = cand;
        }
        if (k == 0) continue;

        LIMDD_STAB h = limdd_stab_make(g, k);
        test_assert(limdd_stab_ngens(h) == k);

        /*
         * A random alternative basis. The coefficient matrix is made unit
         * lower-triangular -- generator i is a random product of generators
         * 0..i-1 times generator i -- because that is invertible over GF(2)
         * and so spans the SAME group. Merely forcing bit i is not enough:
         * two such rows can coincide, and the result is then a subgroup, for
         * which a different handle would be the correct answer.
         */
        LIMDD_LIM alt[4];
        for (size_t i = 0; i < k; i++) {
            const uint64_t bits = (rnd() & (((uint64_t)1 << i) - 1)) | ((uint64_t)1 << i);
            alt[i] = limdd_stab_element(h, bits);
        }
        LIMDD_STAB h2 = limdd_stab_make(alt, k);
        if (h2 != h) {
            fprintf(stderr, "trial %d: same group, two handles (%llu vs %llu)\n",
                    trial, (unsigned long long)h, (unsigned long long)h2);
            fprintf(stderr, "  original:\n");
            limdd_stab_fprint(stderr, h, "    ");
            fprintf(stderr, "  rebuilt:\n");
            limdd_stab_fprint(stderr, h2, "    ");
            return 1;
        }

        /* Every element must be recognised, and no element of the coset -G. */
        for (uint64_t e = 0; e < ((uint64_t)1 << k); e++) {
            LIMDD_LIM el = limdd_stab_element(h, e);
            test_assert(limdd_stab_contains(h, el));
            LIMDD_LIM neg = limdd_lim_make(limdd_lim_pauli(el),
                                           limdd_lim_weight(el) == EVBDD_ONE
                                               ? EVBDD_MIN_ONE : EVBDD_ONE);
            test_assert(!limdd_stab_contains(h, neg));
        }
    }
    return 0;
}

/* --- concurrent construction --------------------------------------------- */

#define NSHARED 256
static LIMDD_LIM  shared_gens[NSHARED][2];
static LIMDD_STAB seen[NSHARED];

VOID_TASK_2(build_shared, size_t, lo, size_t, hi)
{
    if (hi - lo > 32) {
        size_t mid = lo + (hi - lo) / 2;
        SPAWN(build_shared, lo, mid);
        CALL(build_shared, mid, hi);
        SYNC(build_shared);
        return;
    }
    for (size_t i = lo; i < hi; i++) {
        LIMDD_STAB s = limdd_stab_make(shared_gens[i], 2);
        if (seen[i] != 0 && seen[i] != s) {
            fprintf(stderr, "group %zu interned as %llu and also %llu\n",
                    i, (unsigned long long)seen[i], (unsigned long long)s);
            exit(1);
        }
        seen[i] = s;
    }
}

/**
 * Reduction runs in a per-worker scratch buffer. If that buffer were shared,
 * two workers reducing at once would corrupt each other's rows and produce
 * handles that disagree between runs -- which this catches, since every worker
 * reduces the same groups over several rounds.
 */
int
test_extend_skipped(void)
{
    /*
     * Extending a group over skipped levels is a prepend with no elimination,
     * and must give the very handle that building the same generators the
     * slow way gives -- that is what makes the shortcut legitimate.
     */
    for (int trial = 0; trial < 500; trial++) {
        const uint32_t to = rnd() % (NQUBITS + 1);
        const uint32_t from = rnd() % (to + 1);

        /* A random group living on qubits to..n-1 only. */
        const uint64_t above = (to >= 64) ? 0 : ~((UINT64_C(1) << to) - 1);
        const uint64_t mask = ((UINT64_C(1) << NQUBITS) - 1) & above;
        LIMDD_LIM g[NQUBITS + 1];
        size_t k = 0;
        for (int attempt = 0; attempt < 32 && k + to < NQUBITS; attempt++) {
            limdd_pauli_t p = limdd_pauli_identity();
            {
                const uint64_t xw = rnd() & mask, zw = rnd() & mask;
                for (size_t qq = 0; qq < NQUBITS; qq++) {
                    const unsigned op = 2u * ((unsigned)(xw >> qq) & 1u)
                                      +      ((unsigned)(zw >> qq) & 1u);
                    if (op) limdd_pauli_set(&p, qq, (limdd_pauli_op_t)op);
                }
            }
            if (limdd_pauli_is_identity(p)) continue;
            bool ok = true;
            for (size_t i = 0; i < k; i++) {
                if (limdd_pauli_commutation_phase(p, limdd_lim_pauli(g[i])) != 0) ok = false;
            }
            if (!ok) continue;
            const LIMDD_LIM cand = limdd_lim_make(p, (rnd() & 1) ? EVBDD_MIN_ONE : EVBDD_ONE);
            const LIMDD_LIM negc = limdd_lim_make(p, wgt_neg(limdd_lim_weight(cand)));
            const LIMDD_STAB sofar = limdd_stab_make(g, k);
            if (limdd_stab_contains(sofar, cand) || limdd_stab_contains(sofar, negc)) continue;
            g[k++] = cand;
        }
        const LIMDD_STAB s = limdd_stab_make(g, k);

        const LIMDD_STAB fast = limdd_stab_extend_skipped(s, from, to);

        for (uint32_t j = from; j < to; j++) {
            limdd_pauli_t xj = limdd_pauli_single(j, LIMDD_PAULI_X);
            g[k++] = limdd_lim_make(xj, EVBDD_ONE);
        }
        const LIMDD_STAB slow = limdd_stab_make(g, k);

        test_assert(fast == slow);
        test_assert(limdd_stab_ngens(fast) == k);
        for (uint32_t j = from; j < to; j++) {
            limdd_pauli_t xj = limdd_pauli_single(j, LIMDD_PAULI_X);
            test_assert(limdd_stab_contains(fast, limdd_lim_make(xj, EVBDD_ONE)));
        }
    }
    return 0;
}

int
test_concurrent(void)
{
    const uint64_t mask = (UINT64_C(1) << NQUBITS) - 1;
    for (size_t i = 0; i < NSHARED; i++) {
        /* Z-type generators always commute, so every pair is a valid group. */
        limdd_pauli_t p = limdd_pauli_identity();
        { const uint64_t zw = (rnd() & mask) | 1;
          for (size_t qq = 0; qq < NQUBITS; qq++)
              if ((zw >> qq) & 1) limdd_pauli_set(&p, qq, LIMDD_PAULI_Z); }
        limdd_pauli_t q = limdd_pauli_identity();
        { const uint64_t zw = (rnd() & mask) | 2;
          for (size_t qq = 0; qq < NQUBITS; qq++)
              if ((zw >> qq) & 1) limdd_pauli_set(&q, qq, LIMDD_PAULI_Z); }
        shared_gens[i][0] = limdd_lim_make(p, EVBDD_ONE);
        shared_gens[i][1] = limdd_lim_make(q, EVBDD_ONE);
        seen[i] = 0;
    }

    for (int round = 0; round < 8; round++) RUN(build_shared, 0, NSHARED);

    for (size_t i = 0; i < NSHARED; i++) {
        test_assert(seen[i] != 0);
        test_assert(seen[i] == limdd_stab_make(shared_gens[i], 2));
    }
    return 0;
}

TASK_0(int, runtests)
{
    limdd_lims_init(NQUBITS, 1LL << 16, 1LL << 16);
    limdd_stab_init(1LL << 16);

    if (test_trivial()) return 1;
    printf("stab trivial group:                      ok\n");
    if (test_known_groups()) return 1;
    printf("stab known groups and signs:             ok\n");
    if (test_generating_set_invariance()) return 1;
    printf("stab handle depends only on the group:   ok\n");
    if (test_random_invariance()) return 1;
    printf("stab invariant over 2000 random bases:   ok\n");
    if (test_extend_skipped()) return 1;
    printf("stab extension over skipped levels:      ok\n");
    if (test_concurrent()) return 1;
    printf("stab construction agrees across workers: ok\n");

    printf("(%zu generator cells interned)\n", limdd_stab_table_count());

    limdd_stab_quit();
    limdd_lims_quit();
    return 0;
}

int
main(void)
{
    lace_start(8, 0);
    sylvan_set_sizes(1LL << 16, 1LL << 16, 1LL << 16, 1LL << 16);
    sylvan_init_package();
    sylvan_init_edge_weights(1LL << 16, 1LL << 16, 1e-14, WGT_COMPLEX_128, COMP_HASHMAP);

    int res = RUN(runtests);

    sylvan_edge_weights_free();
    sylvan_quit();
    lace_stop();
    return res;
}
