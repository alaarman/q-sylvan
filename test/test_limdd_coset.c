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
 * Tests for limdd_stab_min_coset.
 *
 * The function replaces enumeration of a class with up to 2^(k0+k1+1) members
 * by linear algebra over the words. The test does the enumeration anyway, for
 * small groups, and insists the two agree -- both on which element is least
 * and on the witness that produced it.
 *
 * The property that actually matters is invariance: feeding in any member of
 * the class must give back the same representative. That is what merges two
 * nodes denoting the same state, so it is tested directly, on every member of
 * the class rather than on a sample.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lace.h>
#include <sylvan.h>

#include "qsylvan_limdd_stab.h"
#include "sylvan_edge_weights_complex.h"
#include "test_assert.h"

#define NQUBITS 4

static uint64_t rng_state = UINT64_C(0x243F6A8885A308D3);
static uint64_t
rnd(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

/* --- an order on Pauli words, written out independently ------------------ */

/**
 * Column c of the symplectic vector: the X bits first, then the Z bits. This
 * repeats the library's convention on purpose -- if the two disagreed the test
 * would be comparing minima under different orders, which is worth catching.
 */
static bool
bit_at(limdd_pauli_t p, size_t c)
{
    return (c < NQUBITS) ? ((p.x >> c) & 1) : ((p.z >> (c - NQUBITS)) & 1);
}

static int
word_cmp(limdd_pauli_t a, limdd_pauli_t b)
{
    for (size_t c = 0; c < 2 * NQUBITS; c++) {
        const bool x = bit_at(a, c), y = bit_at(b, c);
        if (x != y) return x ? 1 : -1;
    }
    return 0;
}

static bool
sign_ok(EVBDD_WGT w)
{
    complex_t c;
    weight_value(w, &c);
    if (c.i > 1e-14) return true;
    if (c.i < -1e-14) return false;
    return c.r >= 0.0;
}

/* --- brute force over the class ------------------------------------------ */

/**
 * The least element of { +-G b H : G in s0, H in s1 }, by enumeration.
 * `*count` receives the class size, so the two views of it can be compared.
 */
static LIMDD_LIM
brute_min(LIMDD_LIM b, LIMDD_STAB s0, LIMDD_STAB s1, size_t *count)
{
    const size_t k0 = limdd_stab_ngens(s0);
    const size_t k1 = limdd_stab_ngens(s1);

    LIMDD_LIM best = 0;
    size_t seen = 0;

    for (uint64_t i = 0; i < ((uint64_t)1 << k0); i++) {
        for (uint64_t j = 0; j < ((uint64_t)1 << k1); j++) {
            const LIMDD_LIM g = limdd_stab_element(s0, i);
            const LIMDD_LIM h = limdd_stab_element(s1, j);
            LIMDD_LIM e = limdd_lim_mul(limdd_lim_mul(g, b), h);

            for (int neg = 0; neg < 2; neg++) {
                if (neg) e = limdd_lim_make(limdd_lim_pauli(e),
                                            wgt_neg(limdd_lim_weight(e)));
                seen++;
                if (best == 0) { best = e; continue; }

                const int c = word_cmp(limdd_lim_pauli(e), limdd_lim_pauli(best));
                if (c < 0) best = e;
                else if (c == 0 && sign_ok(limdd_lim_weight(e))
                                && !sign_ok(limdd_lim_weight(best))) best = e;
            }
        }
    }
    *count = seen;
    return best;
}

/** True iff `e` lies in { +-G b H }. */
static bool
in_class(LIMDD_LIM e, LIMDD_LIM b, LIMDD_STAB s0, LIMDD_STAB s1)
{
    const size_t k0 = limdd_stab_ngens(s0);
    const size_t k1 = limdd_stab_ngens(s1);
    for (uint64_t i = 0; i < ((uint64_t)1 << k0); i++) {
        for (uint64_t j = 0; j < ((uint64_t)1 << k1); j++) {
            LIMDD_LIM c = limdd_lim_mul(limdd_lim_mul(limdd_stab_element(s0, i), b),
                                        limdd_stab_element(s1, j));
            if (c == e) return true;
            if (limdd_lim_make(limdd_lim_pauli(c), wgt_neg(limdd_lim_weight(c))) == e)
                return true;
        }
    }
    return false;
}

/* --- random stabiliser groups -------------------------------------------- */

static LIMDD_STAB
random_group(size_t want)
{
    LIMDD_LIM g[NQUBITS];
    size_t k = 0;
    for (int attempt = 0; attempt < 64 && k < want; attempt++) {
        const uint64_t mask = (UINT64_C(1) << NQUBITS) - 1;
        limdd_pauli_t p;
        p.x = rnd() & mask;
        p.z = rnd() & mask;
        if (limdd_pauli_is_identity(p)) continue;

        bool ok = true;
        for (size_t i = 0; i < k; i++) {
            if (limdd_pauli_commutation_phase(p, limdd_lim_pauli(g[i])) != 0) ok = false;
        }
        if (!ok) continue;

        LIMDD_LIM cand = limdd_lim_make(p, (rnd() & 1) ? EVBDD_MIN_ONE : EVBDD_ONE);
        LIMDD_LIM negc = limdd_lim_make(p, wgt_neg(limdd_lim_weight(cand)));
        LIMDD_STAB sofar = limdd_stab_make(g, k);
        if (limdd_stab_contains(sofar, cand) || limdd_stab_contains(sofar, negc)) continue;

        g[k++] = cand;
    }
    return limdd_stab_make(g, k);
}

static LIMDD_LIM
random_lim(void)
{
    const uint64_t mask = (UINT64_C(1) << NQUBITS) - 1;
    limdd_pauli_t p;
    p.x = rnd() & mask;
    p.z = rnd() & mask;
    static const double re[6] = { 1.0, 0.0, -1.0,  0.0, 0.5, -0.25 };
    static const double im[6] = { 0.0, 1.0,  0.0, -1.0, 0.0,  0.5  };
    const unsigned w = rnd() % 6;
    return limdd_lim_make(p, complex_lookup(re[w], im[w]));
}

/* --- tests --------------------------------------------------------------- */

int
test_trivial_groups(void)
{
    /* With no freedom at all the answer is b, up to the sign rule. */
    for (int i = 0; i < 200; i++) {
        LIMDD_LIM b = random_lim();
        LIMDD_LIM w; bool neg;
        LIMDD_LIM m = limdd_stab_min_coset(b, LIMDD_STAB_TRIVIAL, LIMDD_STAB_TRIVIAL,
                                           &w, &neg);
        test_assert(w == LIMDD_LIM_IDENTITY);
        test_assert(limdd_pauli_equals(limdd_lim_pauli(m), limdd_lim_pauli(b)));
        test_assert(sign_ok(limdd_lim_weight(m)));
        test_assert(neg == !sign_ok(limdd_lim_weight(b)));
    }
    return 0;
}

int
test_matches_enumeration(void)
{
    for (int trial = 0; trial < 3000; trial++) {
        LIMDD_STAB s0 = random_group(1 + (rnd() % 3));
        LIMDD_STAB s1 = random_group(1 + (rnd() % 3));
        LIMDD_LIM b = random_lim();

        LIMDD_LIM witness; bool neg;
        LIMDD_LIM got = limdd_stab_min_coset(b, s0, s1, &witness, &neg);

        size_t seen;
        LIMDD_LIM want = brute_min(b, s0, s1, &seen);

        if (got != want) {
            fprintf(stderr, "trial %d: coset minimum disagrees\n", trial);
            fprintf(stderr, "  b        = "); limdd_lim_fprint(stderr, b);       fprintf(stderr, "\n");
            fprintf(stderr, "  computed = "); limdd_lim_fprint(stderr, got);     fprintf(stderr, "\n");
            fprintf(stderr, "  brute    = "); limdd_lim_fprint(stderr, want);    fprintf(stderr, "\n");
            fprintf(stderr, "  s0:\n"); limdd_stab_fprint(stderr, s0, "    ");
            fprintf(stderr, "  s1:\n"); limdd_stab_fprint(stderr, s1, "    ");
            return 1;
        }
        test_assert(in_class(got, b, s0, s1));
    }
    return 0;
}

int
test_invariance(void)
{
    /*
     * Every member of the class must reduce to the same representative. This
     * is the property the canonical form is built on: two nodes that differ
     * only by a LIM must end up with identical high labels, or they will never
     * be merged.
     */
    for (int trial = 0; trial < 1500; trial++) {
        LIMDD_STAB s0 = random_group(1 + (rnd() % 3));
        LIMDD_STAB s1 = random_group(1 + (rnd() % 3));
        LIMDD_LIM b = random_lim();

        LIMDD_LIM ref = limdd_stab_min_coset(b, s0, s1, NULL, NULL);

        const size_t k0 = limdd_stab_ngens(s0);
        const size_t k1 = limdd_stab_ngens(s1);
        for (uint64_t i = 0; i < ((uint64_t)1 << k0); i++) {
            for (uint64_t j = 0; j < ((uint64_t)1 << k1); j++) {
                LIMDD_LIM alt = limdd_lim_mul(
                    limdd_lim_mul(limdd_stab_element(s0, i), b),
                    limdd_stab_element(s1, j));
                for (int neg = 0; neg < 2; neg++) {
                    if (neg) alt = limdd_lim_make(limdd_lim_pauli(alt),
                                                  wgt_neg(limdd_lim_weight(alt)));
                    if (limdd_stab_min_coset(alt, s0, s1, NULL, NULL) != ref) {
                        fprintf(stderr, "trial %d: class member maps elsewhere\n", trial);
                        fprintf(stderr, "  b   = "); limdd_lim_fprint(stderr, b);   fprintf(stderr, "\n");
                        fprintf(stderr, "  alt = "); limdd_lim_fprint(stderr, alt); fprintf(stderr, "\n");
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

int
test_witness(void)
{
    /*
     * The witness is what the parent edge is built from, so a wrong one would
     * change the state while leaving the node looking right. It must satisfy
     * min = +-G b H for some H in s1, with the sign that `negated` reported.
     */
    for (int trial = 0; trial < 2000; trial++) {
        LIMDD_STAB s0 = random_group(1 + (rnd() % 3));
        LIMDD_STAB s1 = random_group(1 + (rnd() % 3));
        LIMDD_LIM b = random_lim();

        LIMDD_LIM g; bool neg;
        LIMDD_LIM m = limdd_stab_min_coset(b, s0, s1, &g, &neg);

        test_assert(limdd_stab_contains(s0, g) || g == LIMDD_LIM_IDENTITY);

        /* Undo the sign and the witness: what is left must be b times an
         * element of s1. */
        LIMDD_LIM t = neg ? limdd_lim_make(limdd_lim_pauli(m),
                                           wgt_neg(limdd_lim_weight(m))) : m;
        LIMDD_LIM h = limdd_lim_mul(limdd_lim_inverse(b), limdd_lim_mul(g, t));
        if (!limdd_stab_contains(s1, h)) {
            fprintf(stderr, "trial %d: witness does not reconstruct the minimum\n", trial);
            fprintf(stderr, "  b = "); limdd_lim_fprint(stderr, b); fprintf(stderr, "\n");
            fprintf(stderr, "  g = "); limdd_lim_fprint(stderr, g); fprintf(stderr, "\n");
            fprintf(stderr, "  m = "); limdd_lim_fprint(stderr, m); fprintf(stderr, "\n");
            fprintf(stderr, "  h = "); limdd_lim_fprint(stderr, h); fprintf(stderr, "\n");
            return 1;
        }
    }
    return 0;
}

TASK_0(int, runtests)
{
    limdd_lims_init(NQUBITS, 1LL << 18, 1LL << 18);
    limdd_stab_init(1LL << 18);

    if (test_trivial_groups()) return 1;
    printf("coset min with no freedom is the sign rule:  ok\n");
    if (test_matches_enumeration()) return 1;
    printf("coset min matches enumeration (3000 cases):  ok\n");
    if (test_invariance()) return 1;
    printf("every class member maps to one representative: ok\n");
    if (test_witness()) return 1;
    printf("witness reconstructs the minimum:            ok\n");

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
    sylvan_init_edge_weights(1LL << 18, 1LL << 18, 1e-14, WGT_COMPLEX_128, COMP_HASHMAP);

    int res = RUN(runtests);

    sylvan_edge_weights_free();
    sylvan_quit();
    lace_stop();
    return res;
}
