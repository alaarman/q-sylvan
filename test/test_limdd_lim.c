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
 * Tests for interned LIMs.
 *
 * The property that matters is canonicity: equal maps must get equal indices,
 * because LIMDD node lookup compares labels by index. Most of what follows is
 * therefore about establishing that two routes to the same map agree.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lace.h>
#include <sylvan.h>

#include "qsylvan_limdd_lim.h"
#include "sylvan_edge_weights_complex.h"
#include "test_assert.h"

#define NQUBITS 8

static uint64_t rng_state = UINT64_C(0x123456789ABCDEF);
static uint64_t
rnd(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

static limdd_pauli_t
random_pauli(void)
{
    const uint64_t mask = (UINT64_C(1) << NQUBITS) - 1;
    const uint64_t xw = rnd() & mask, zw = rnd() & mask;
    limdd_pauli_t p = limdd_pauli_identity();
    for (size_t q = 0; q < NQUBITS; q++) {
        const unsigned op = 2u * ((unsigned)(xw >> q) & 1u)
                          +      ((unsigned)(zw >> q) & 1u);
        if (op) limdd_pauli_set(&p, q, (limdd_pauli_op_t)op);
    }
    return p;
}

/* A small set of exactly representable weights, so equality is exact. */
static EVBDD_WGT
random_weight(void)
{
    const double vals[4] = { 1.0, -1.0, 0.5, -0.25 };
    return complex_lookup(vals[rnd() & 3], vals[rnd() & 3]);
}

int
test_interning(void)
{
    /* The same word must always come back as the same reference. */
    for (int i = 0; i < 200; i++) {
        limdd_pauli_t p = random_pauli();
        LIMDD_PAULI_REF a = limdd_pauli_intern(p);
        LIMDD_PAULI_REF b = limdd_pauli_intern(p);
        test_assert(a == b);
        test_assert(a != 0);
        test_assert(limdd_pauli_equals(limdd_pauli_deref(a), p));
    }

    /* Different words must not collide. */
    limdd_pauli_t p1, p2;
    test_assert(limdd_pauli_from_string("XIIIIIII", NQUBITS, &p1));
    test_assert(limdd_pauli_from_string("ZIIIIIII", NQUBITS, &p2));
    test_assert(limdd_pauli_intern(p1) != limdd_pauli_intern(p2));

    /* Likewise for whole LIMs. */
    for (int i = 0; i < 200; i++) {
        limdd_pauli_t p = random_pauli();
        EVBDD_WGT w = random_weight();
        LIMDD_LIM a = limdd_lim_make(p, w);
        LIMDD_LIM b = limdd_lim_make(p, w);
        test_assert(a == b);
        test_assert(limdd_pauli_equals(limdd_lim_pauli(a), p));
        test_assert(limdd_lim_weight(a) == w);
    }
    return 0;
}

int
test_zero_and_identity(void)
{
    test_assert(limdd_lim_is_identity(LIMDD_LIM_IDENTITY));
    test_assert(limdd_lim_is_zero(LIMDD_LIM_ZERO));
    test_assert(LIMDD_LIM_IDENTITY != LIMDD_LIM_ZERO);

    /* The identity must be reachable the ordinary way. */
    test_assert(limdd_lim_make(limdd_pauli_identity(), EVBDD_ONE) == LIMDD_LIM_IDENTITY);

    /*
     * Every zero-scalar LIM is the same map whatever its Pauli word, so they
     * must all collapse onto one index. If they did not, node lookup would
     * treat two identical nodes as different.
     */
    for (int i = 0; i < 100; i++) {
        test_assert(limdd_lim_make(random_pauli(), EVBDD_ZERO) == LIMDD_LIM_ZERO);
    }
    test_assert(limdd_lim_weight(LIMDD_LIM_ZERO) == EVBDD_ZERO);

    /* Zero absorbs under composition. */
    for (int i = 0; i < 50; i++) {
        LIMDD_LIM x = limdd_lim_make(random_pauli(), random_weight());
        test_assert(limdd_lim_mul(x, LIMDD_LIM_ZERO) == LIMDD_LIM_ZERO);
        test_assert(limdd_lim_mul(LIMDD_LIM_ZERO, x) == LIMDD_LIM_ZERO);
        test_assert(limdd_lim_mul(x, LIMDD_LIM_IDENTITY) == x);
        test_assert(limdd_lim_mul(LIMDD_LIM_IDENTITY, x) == x);
    }
    return 0;
}

int
test_multiplication(void)
{
    /*
     * X * Z = -iY and Z * X = +iY. Checking both directions pins down that the
     * phase is folded in with the right sign, and that composition is not
     * accidentally commutative.
     */
    limdd_pauli_t px, pz;
    test_assert(limdd_pauli_from_string("XIIIIIII", NQUBITS, &px));
    test_assert(limdd_pauli_from_string("ZIIIIIII", NQUBITS, &pz));

    LIMDD_LIM x = limdd_lim_make(px, EVBDD_ONE);
    LIMDD_LIM z = limdd_lim_make(pz, EVBDD_ONE);

    LIMDD_LIM xz = limdd_lim_mul(x, z);
    LIMDD_LIM zx = limdd_lim_mul(z, x);

    test_assert(xz != zx);                                   /* X and Z anticommute */
    test_assert(limdd_pauli_equals(limdd_lim_pauli(xz),
                                   limdd_lim_pauli(zx)));    /* both are Y... */
    test_assert(limdd_lim_weight(xz) != limdd_lim_weight(zx)); /* ...differing by a sign */

    char buf[LIMDD_MAX_QUBITS + 1];
    limdd_pauli_to_string(limdd_lim_pauli(xz), NQUBITS, buf);
    test_assert(strcmp(buf, "YIIIIIII") == 0);

    test_assert(limdd_lim_weight(xz) == complex_lookup(0.0, -1.0)); /*  X*Z = -iY */
    test_assert(limdd_lim_weight(zx) == complex_lookup(0.0,  1.0)); /*  Z*X = +iY */

    /* A LIM composed with its own inverse is the identity. */
    for (int i = 0; i < 200; i++) {
        LIMDD_LIM a = limdd_lim_make(random_pauli(), random_weight());
        LIMDD_LIM inv = limdd_lim_inverse(a);
        test_assert(limdd_lim_mul(a, inv) == LIMDD_LIM_IDENTITY);
        test_assert(limdd_lim_mul(inv, a) == LIMDD_LIM_IDENTITY);
    }
    return 0;
}

int
test_associativity(void)
{
    /*
     * (a*b)*c and a*(b*c) are the same map, so canonicity demands the same
     * index. This is the strongest cheap check available: it exercises the
     * phase accumulation, the scalar multiplication and the interning together,
     * and any of the three getting it wrong shows up here.
     *
     * Weights are restricted to powers of i so that the products stay exactly
     * representable and equality is not a floating point question.
     */
    const EVBDD_WGT units[4] = {
        EVBDD_ONE,
        complex_lookup(0.0, 1.0),
        EVBDD_MIN_ONE,
        complex_lookup(0.0, -1.0),
    };

    for (int trial = 0; trial < 2000; trial++) {
        LIMDD_LIM a = limdd_lim_make(random_pauli(), units[rnd() & 3]);
        LIMDD_LIM b = limdd_lim_make(random_pauli(), units[rnd() & 3]);
        LIMDD_LIM c = limdd_lim_make(random_pauli(), units[rnd() & 3]);

        test_assert(limdd_lim_mul(limdd_lim_mul(a, b), c)
                    == limdd_lim_mul(a, limdd_lim_mul(b, c)));
    }
    return 0;
}

/* --- concurrent interning ------------------------------------------------ */

#define NSHARED 512
static limdd_pauli_t shared_pauli[NSHARED];
static EVBDD_WGT     shared_wgt[NSHARED];
static LIMDD_LIM     seen[NSHARED];

/**
 * Every worker interns the same LIMs. They must all get the same indices:
 * that is the whole point of a lock-free find-or-put, and it exercises the
 * Pauli table, the LIM table and the edge weight table at once, from several
 * workers, which is the configuration the multi-table region cursors exist for.
 */
VOID_TASK_2(intern_shared, size_t, lo, size_t, hi)
{
    if (hi - lo > 32) {
        size_t mid = lo + (hi - lo) / 2;
        SPAWN(intern_shared, lo, mid);
        CALL(intern_shared, mid, hi);
        SYNC(intern_shared);
        return;
    }
    for (size_t i = lo; i < hi; i++) {
        LIMDD_LIM lim = limdd_lim_make(shared_pauli[i], shared_wgt[i]);
        if (seen[i] != 0 && seen[i] != lim) {
            fprintf(stderr, "LIM %zu interned as %llu and also %llu\n",
                    i, (unsigned long long)seen[i], (unsigned long long)lim);
            exit(1);
        }
        seen[i] = lim;
    }
}

int
test_concurrent_interning(void)
{
    for (size_t i = 0; i < NSHARED; i++) {
        shared_pauli[i] = random_pauli();
        shared_wgt[i] = random_weight();
        seen[i] = 0;
    }

    /* Several rounds, so workers race on both fresh and existing entries. */
    for (int round = 0; round < 8; round++) RUN(intern_shared, 0, NSHARED);

    /* And the single-threaded answer must agree with what the workers found. */
    for (size_t i = 0; i < NSHARED; i++) {
        test_assert(seen[i] != 0);
        test_assert(seen[i] == limdd_lim_make(shared_pauli[i], shared_wgt[i]));
    }
    return 0;
}

TASK_0(int, runtests)
{
    limdd_lims_init(NQUBITS, 1LL << 16, 1LL << 16);

    if (test_interning()) return 1;
    printf("lim interning is idempotent:      ok\n");
    if (test_zero_and_identity()) return 1;
    printf("lim zero and identity collapse:   ok\n");
    if (test_multiplication()) return 1;
    printf("lim multiplication and inverse:   ok\n");
    if (test_associativity()) return 1;
    printf("lim multiplication associates:    ok\n");
    if (test_concurrent_interning()) return 1;
    printf("lim interning agrees across workers: ok\n");

    printf("(%zu Pauli words, %zu LIMs interned)\n",
           limdd_pauli_table_count(), limdd_lim_table_count());

    limdd_lims_quit();
    return 0;
}

int
main(void)
{
    lace_start(8, 0);

    /* The node table and, more to the point, the operation cache: weight
     * arithmetic goes through the cache, so wgt_mul needs the package up. */
    sylvan_set_sizes(1LL << 16, 1LL << 16, 1LL << 16, 1LL << 16);
    sylvan_init_package();

    /* LIMs need complex weights: a Pauli product can produce a factor of i. */
    sylvan_init_edge_weights(1LL << 16, 1LL << 16, 1e-14, WGT_COMPLEX_128, COMP_HASHMAP);

    int res = RUN(runtests);

    sylvan_edge_weights_free();
    sylvan_quit();
    lace_stop();
    return res;
}
