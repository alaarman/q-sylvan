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
 * Tests for the packed Pauli words used as LIMDD edge labels.
 *
 * The full 4x4 multiplication table below is written out by hand from the
 * definitions of the Pauli matrices rather than generated, so that this test
 * is an independent statement of what the code should do. It was also
 * cross-checked against an external reference implementation during
 * development.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "qsylvan_limdd_pauli.h"
#include "test_assert.h"

/* A tiny deterministic PRNG, so failures reproduce exactly. */
static uint64_t rng_state = UINT64_C(0x243F6A8885A308D3);
static uint64_t
rnd(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

static limdd_pauli_t
random_word(size_t n)
{
    const uint64_t mask = (n >= 64) ? UINT64_MAX : ((UINT64_C(1) << n) - 1);
    const uint64_t xw = rnd() & mask, zw = rnd() & mask;
    limdd_pauli_t p = limdd_pauli_identity();
    for (size_t q = 0; q < n; q++) {
        const unsigned op = 2u * ((unsigned)(xw >> q) & 1u)
                          +      ((unsigned)(zw >> q) & 1u);
        if (op) limdd_pauli_set(&p, q, (limdd_pauli_op_t)op);
    }
    return p;
}

int
test_get_set(void)
{
    limdd_pauli_t p = limdd_pauli_identity();
    test_assert(limdd_pauli_is_identity(p));
    test_assert(limdd_pauli_weight(p) == 0);

    /* Setting one qubit must not disturb any other. */
    limdd_pauli_set(&p, 5, LIMDD_PAULI_Y);
    test_assert(limdd_pauli_get(p, 5) == LIMDD_PAULI_Y);
    test_assert(limdd_pauli_get(p, 4) == LIMDD_PAULI_I);
    test_assert(limdd_pauli_get(p, 6) == LIMDD_PAULI_I);
    test_assert(!limdd_pauli_is_identity(p));
    test_assert(limdd_pauli_weight(p) == 1);

    limdd_pauli_set(&p, 5, LIMDD_PAULI_X);
    test_assert(limdd_pauli_get(p, 5) == LIMDD_PAULI_X);
    limdd_pauli_set(&p, 5, LIMDD_PAULI_Z);
    test_assert(limdd_pauli_get(p, 5) == LIMDD_PAULI_Z);

    /* Setting back to I must clear both components, not just one. */
    limdd_pauli_set(&p, 5, LIMDD_PAULI_I);
    test_assert(limdd_pauli_is_identity(p));

    /* The encoding is 2*x + z; check it explicitly. */
    limdd_pauli_t q = limdd_pauli_identity();
    limdd_pauli_set(&q, 0, LIMDD_PAULI_X);
    test_assert(limdd_pauli_get(q, 0) == LIMDD_PAULI_X);
    limdd_pauli_set(&q, 0, LIMDD_PAULI_Z);
    test_assert(limdd_pauli_get(q, 0) == LIMDD_PAULI_Z);
    limdd_pauli_set(&q, 0, LIMDD_PAULI_Y);
    test_assert(limdd_pauli_get(q, 0) == LIMDD_PAULI_Y);

    /* The top qubit must work, not just low ones. */
    limdd_pauli_t top = limdd_pauli_identity();
    limdd_pauli_set(&top, 63, LIMDD_PAULI_Y);
    test_assert(limdd_pauli_get(top, 63) == LIMDD_PAULI_Y);
    test_assert(limdd_pauli_weight(top) == 1);
    test_assert(limdd_pauli_is_canonical(top, 64));
    test_assert(!limdd_pauli_is_canonical(top, 63));

    return 0;
}

int
test_multiplication_table(void)
{
    /*
     * rows = left operand, cols = right operand, order I, Z, X, Y.
     * Each entry is the resulting Pauli and the power of i, from
     *   ZX = iY   XY = iZ   YZ = iX
     *   XZ = -iY  YX = -iZ  ZY = -iX
     * and P*P = I, I*P = P*I = P.
     */
    static const limdd_pauli_op_t want_op[4][4] = {
        /*        I              Z              X              Y      */
        /* I */ { LIMDD_PAULI_I, LIMDD_PAULI_Z, LIMDD_PAULI_X, LIMDD_PAULI_Y },
        /* Z */ { LIMDD_PAULI_Z, LIMDD_PAULI_I, LIMDD_PAULI_Y, LIMDD_PAULI_X },
        /* X */ { LIMDD_PAULI_X, LIMDD_PAULI_Y, LIMDD_PAULI_I, LIMDD_PAULI_Z },
        /* Y */ { LIMDD_PAULI_Y, LIMDD_PAULI_X, LIMDD_PAULI_Z, LIMDD_PAULI_I },
    };
    static const unsigned want_phase[4][4] = {
        /*        I  Z  X  Y */
        /* I */ { 0, 0, 0, 0 },
        /* Z */ { 0, 0, 1, 3 },
        /* X */ { 0, 3, 0, 1 },
        /* Y */ { 0, 1, 3, 0 },
    };

    for (int a = 0; a < 4; a++) {
        for (int b = 0; b < 4; b++) {
            limdd_pauli_t pa = limdd_pauli_identity();
            limdd_pauli_t pb = limdd_pauli_identity();
            limdd_pauli_set(&pa, 0, (limdd_pauli_op_t)a);
            limdd_pauli_set(&pb, 0, (limdd_pauli_op_t)b);

            limdd_pauli_t r = pa;
            unsigned phase = limdd_pauli_rightmul(&r, pb, 1);
            test_assert(limdd_pauli_get(r, 0) == want_op[a][b]);
            test_assert(phase == want_phase[a][b]);

            /* Left-multiplying a by b is right-multiplying b by a. */
            limdd_pauli_t l = pb;
            unsigned lphase = limdd_pauli_leftmul(&l, pa, 1);
            test_assert(limdd_pauli_get(l, 0) == want_op[a][b]);
            test_assert(lphase == want_phase[a][b]);
        }
    }
    return 0;
}

int
test_multiplication_is_elementwise(void)
{
    /*
     * A product over many qubits must be the product qubit by qubit, with
     * the phases summed mod 4. This is what catches a bit-parallel phase
     * formula that is right for one qubit but wrong in aggregate.
     */
    for (size_t n = 1; n <= 64; n++) {
        for (int trial = 0; trial < 200; trial++) {
            limdd_pauli_t a = random_word(n);
            limdd_pauli_t b = random_word(n);

            unsigned want_phase = 0;
            for (size_t i = 0; i < n; i++) {
                limdd_pauli_t sa = limdd_pauli_identity();
                limdd_pauli_t sb = limdd_pauli_identity();
                limdd_pauli_set(&sa, 0, limdd_pauli_get(a, i));
                limdd_pauli_set(&sb, 0, limdd_pauli_get(b, i));
                want_phase = (want_phase + limdd_pauli_rightmul(&sa, sb, 1)) & 3u;
            }

            limdd_pauli_t got = a;
            unsigned got_phase = limdd_pauli_rightmul(&got, b, n);

            test_assert(got_phase == want_phase);
            {   /* the word is the symplectic sum, phase aside */
                limdd_pauli_t want = a;
                limdd_pauli_xor(&want, b);
                test_assert(limdd_pauli_equals(got, want));
            }
            test_assert(limdd_pauli_is_canonical(got, n));
        }
    }
    return 0;
}

int
test_self_inverse(void)
{
    /* P * P = I with no leftover phase; this is why there is no inverse(). */
    for (size_t n = 1; n <= 64; n++) {
        for (int trial = 0; trial < 50; trial++) {
            limdd_pauli_t p = random_word(n);
            limdd_pauli_t q = p;
            unsigned phase = limdd_pauli_rightmul(&q, p, n);
            test_assert(limdd_pauli_is_identity(q));
            test_assert(phase == 0);
        }
    }
    return 0;
}

int
test_commutation(void)
{
    /* Single-qubit facts: distinct non-identity Paulis anticommute. */
    const limdd_pauli_op_t ops[3] = { LIMDD_PAULI_X, LIMDD_PAULI_Y, LIMDD_PAULI_Z };
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            limdd_pauli_t a = limdd_pauli_identity(), b = limdd_pauli_identity();
            limdd_pauli_set(&a, 0, ops[i]);
            limdd_pauli_set(&b, 0, ops[j]);
            test_assert(limdd_pauli_commutes(a, b) == (i == j));
        }
    }

    /* Paulis on disjoint qubits always commute. */
    limdd_pauli_t p = limdd_pauli_identity(), q = limdd_pauli_identity();
    limdd_pauli_set(&p, 1, LIMDD_PAULI_X);
    limdd_pauli_set(&q, 2, LIMDD_PAULI_Z);
    test_assert(limdd_pauli_commutes(p, q));

    /* In general the reported phase must match doing both products. */
    for (size_t n = 1; n <= 64; n++) {
        for (int trial = 0; trial < 100; trial++) {
            limdd_pauli_t a = random_word(n);
            limdd_pauli_t b = random_word(n);

            limdd_pauli_t ab = a, ba = b;
            unsigned p_ab = limdd_pauli_rightmul(&ab, b, n);
            unsigned p_ba = limdd_pauli_rightmul(&ba, a, n);

            test_assert(limdd_pauli_equals(ab, ba));

            unsigned want = (p_ba + 4u - p_ab) & 3u;
            unsigned got = limdd_pauli_commutation_phase(a, b);
            test_assert(got == want);
            test_assert(got == 0 || got == 2);
            test_assert(limdd_pauli_commutes(a, b) == (got == 0));
        }
    }
    return 0;
}

int
test_local_multiplication(void)
{
    /* The one-qubit helpers must agree with the general routine. */
    for (size_t n = 1; n <= 64; n++) {
        for (int trial = 0; trial < 20; trial++) {
            limdd_pauli_t base = random_word(n);
            size_t i = (size_t)(rnd() % n);

            const limdd_pauli_op_t ops[3] = { LIMDD_PAULI_X, LIMDD_PAULI_Y, LIMDD_PAULI_Z };
            unsigned (*rmul[3])(limdd_pauli_t *, size_t) = {
                limdd_pauli_rightmul_x, limdd_pauli_rightmul_y, limdd_pauli_rightmul_z };
            unsigned (*lmul[3])(limdd_pauli_t *, size_t) = {
                limdd_pauli_leftmul_x, limdd_pauli_leftmul_y, limdd_pauli_leftmul_z };

            for (int c = 0; c < 3; c++) {
                limdd_pauli_t one = limdd_pauli_identity();
                limdd_pauli_set(&one, i, ops[c]);

                limdd_pauli_t want = base, got = base;
                unsigned wp = limdd_pauli_rightmul(&want, one, n);
                unsigned gp = rmul[c](&got, i);
                test_assert(wp == gp);
                test_assert(limdd_pauli_equals(want, got));

                want = base; got = base;
                wp = limdd_pauli_leftmul(&want, one, n);
                gp = lmul[c](&got, i);
                test_assert(wp == gp);
                test_assert(limdd_pauli_equals(want, got));
            }
        }
    }
    return 0;
}

int
test_strings(void)
{
    limdd_pauli_t p;
    char buf[LIMDD_MAX_QUBITS + 1];

    test_assert(limdd_pauli_from_string("IZXY", 4, &p));
    test_assert(limdd_pauli_get(p, 0) == LIMDD_PAULI_I);
    test_assert(limdd_pauli_get(p, 1) == LIMDD_PAULI_Z);
    test_assert(limdd_pauli_get(p, 2) == LIMDD_PAULI_X);
    test_assert(limdd_pauli_get(p, 3) == LIMDD_PAULI_Y);

    limdd_pauli_to_string(p, 4, buf);
    test_assert(strcmp(buf, "IZXY") == 0);

    /* Malformed input must be rejected, not silently accepted. */
    limdd_pauli_t unused;
    test_assert(!limdd_pauli_from_string("IZXQ", 4, &unused));   /* bad character */
    test_assert(!limdd_pauli_from_string("IZX", 4, &unused));    /* too short */
    test_assert(!limdd_pauli_from_string("IZXYZ", 4, &unused));  /* too long */
    test_assert(!limdd_pauli_from_string(NULL, 4, &unused));

    /* Round-trip at full width. */
    for (int trial = 0; trial < 100; trial++) {
        limdd_pauli_t r = random_word(64), back;
        limdd_pauli_to_string(r, 64, buf);
        test_assert(limdd_pauli_from_string(buf, 64, &back));
        test_assert(limdd_pauli_equals(r, back));
    }
    return 0;
}

int
runtests(void)
{
    if (test_get_set()) return 1;
    printf("pauli get/set:                 ok\n");
    if (test_multiplication_table()) return 1;
    printf("pauli multiplication table:    ok\n");
    if (test_multiplication_is_elementwise()) return 1;
    printf("pauli multiplication (wide):   ok\n");
    if (test_self_inverse()) return 1;
    printf("pauli self-inverse:            ok\n");
    if (test_commutation()) return 1;
    printf("pauli commutation:             ok\n");
    if (test_local_multiplication()) return 1;
    printf("pauli single-qubit multiply:   ok\n");
    if (test_strings()) return 1;
    printf("pauli string conversion:       ok\n");
    return 0;
}

int
main(void)
{
    /*
     * No Lace and no Sylvan tables: packed Pauli words are plain values with
     * no shared state, so this test needs no runtime at all.
     */
    return runtests();
}
