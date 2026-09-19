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
#include <string.h>

#include <sylvan_platform.h>

#include "qsylvan_limdd_pauli.h"

/**
 * Mask with the low `nqubits` bits set. Written to avoid the undefined
 * behaviour of shifting a 64-bit value by 64.
 */
static inline uint64_t
qubit_mask(size_t nqubits)
{
    if (nqubits >= 64) return UINT64_MAX;
    return (UINT64_C(1) << nqubits) - 1;
}

bool
limdd_pauli_is_canonical(limdd_pauli_t p, size_t nqubits)
{
    if (nqubits > LIMDD_MAX_QUBITS) return false;
    const uint64_t mask = qubit_mask(nqubits);
    return (p.x & ~mask) == 0 && (p.z & ~mask) == 0;
}

size_t
limdd_pauli_weight(limdd_pauli_t p)
{
    return popcnt_uint64(p.x | p.z);
}

/**
 * The phase exponent of a * b, where `a` and `b` are the packed words.
 *
 * This is the bit-parallel form of the `g` function of the rowsum subroutine
 * in Aaronson and Gottesman, "Improved Simulation of Stabilizer Circuits"
 * (https://scottaaronson.com/papers/chp5.pdf). Per qubit, g contributes
 * -1, 0 or +1 to the exponent:
 *
 *      a_i = Z  ->  x_b(1 - 2 z_b)   +1 on  X    -1 on  Y
 *      a_i = X  ->  z_b(2 x_b - 1)   +1 on  Y    -1 on  Z
 *      a_i = Y  ->  z_b - x_b        +1 on  Z    -1 on  X
 *      a_i = I  ->  0
 *
 * Rather than loop, build one mask of the qubits contributing +1 and one of
 * those contributing -1, then take the difference of their popcounts mod 4.
 */
static inline unsigned
mul_phase(limdd_pauli_t a, limdd_pauli_t b)
{
    const uint64_t a_is_z = ~a.x & a.z;
    const uint64_t a_is_x = a.x & ~a.z;
    const uint64_t a_is_y = a.x & a.z;

    const uint64_t plus = (a_is_z & b.x & ~b.z)   // Z * X = +i Y
                        | (a_is_x & b.z & b.x)    // X * Y = +i Z
                        | (a_is_y & b.z & ~b.x);  // Y * Z = +i X

    const uint64_t minus = (a_is_z & b.x & b.z)   // Z * Y = -i X
                         | (a_is_x & b.z & ~b.x)  // X * Z = -i Y
                         | (a_is_y & b.x & ~b.z); // Y * X = -i Z

    const unsigned p = popcnt_uint64(plus);
    const unsigned m = popcnt_uint64(minus);

    // (p - m) mod 4. Use 3m rather than -m so the intermediate stays
    // non-negative in unsigned arithmetic; 3 == -1 (mod 4).
    return (p + 3u * m) & 3u;
}

/**
 * a * b, without the canonicity assertions. Callers that construct `b`
 * themselves know it is well formed.
 */
static inline unsigned
rightmul_raw(limdd_pauli_t *a, limdd_pauli_t b)
{
    const unsigned phase = mul_phase(*a, b);
    a->x ^= b.x;
    a->z ^= b.z;
    return phase;
}

/**
 * b * a, without the canonicity assertions.
 */
static inline unsigned
leftmul_raw(limdd_pauli_t *a, limdd_pauli_t b)
{
    // b is now the left operand, so it comes first in the phase function.
    const unsigned phase = mul_phase(b, *a);
    a->x ^= b.x;
    a->z ^= b.z;
    return phase;
}

unsigned
limdd_pauli_rightmul(limdd_pauli_t *a, limdd_pauli_t b, size_t nqubits)
{
    assert(limdd_pauli_is_canonical(*a, nqubits));
    assert(limdd_pauli_is_canonical(b, nqubits));
    (void)nqubits; // only used for the assertions above

    return rightmul_raw(a, b);
}

unsigned
limdd_pauli_leftmul(limdd_pauli_t *a, limdd_pauli_t b, size_t nqubits)
{
    assert(limdd_pauli_is_canonical(*a, nqubits));
    assert(limdd_pauli_is_canonical(b, nqubits));
    (void)nqubits;

    return leftmul_raw(a, b);
}

unsigned
limdd_pauli_commutation_phase(limdd_pauli_t a, limdd_pauli_t b)
{
    // Symplectic inner product: a and b anticommute iff it is odd.
    const unsigned ip = popcnt_uint64(a.x & b.z) + popcnt_uint64(a.z & b.x);
    return (ip & 1u) ? 2u : 0u;
}

/**
 * A word that is `op` on `index` and the identity everywhere else.
 *
 * The one-qubit operations below go through the raw multiply rather than the
 * checked public one: the word built here is canonical by construction, and
 * these helpers are not told how many qubits `*p` spans, so they are in no
 * position to check it.
 */
static inline limdd_pauli_t
single(size_t index, limdd_pauli_op_t op)
{
    assert(index < LIMDD_MAX_QUBITS);
    limdd_pauli_t p = limdd_pauli_identity();
    limdd_pauli_set(&p, index, op);
    return p;
}

unsigned
limdd_pauli_rightmul_x(limdd_pauli_t *p, size_t index)
{
    return rightmul_raw(p, single(index, LIMDD_PAULI_X));
}

unsigned
limdd_pauli_rightmul_y(limdd_pauli_t *p, size_t index)
{
    return rightmul_raw(p, single(index, LIMDD_PAULI_Y));
}

unsigned
limdd_pauli_rightmul_z(limdd_pauli_t *p, size_t index)
{
    return rightmul_raw(p, single(index, LIMDD_PAULI_Z));
}

unsigned
limdd_pauli_leftmul_x(limdd_pauli_t *p, size_t index)
{
    return leftmul_raw(p, single(index, LIMDD_PAULI_X));
}

unsigned
limdd_pauli_leftmul_y(limdd_pauli_t *p, size_t index)
{
    return leftmul_raw(p, single(index, LIMDD_PAULI_Y));
}

unsigned
limdd_pauli_leftmul_z(limdd_pauli_t *p, size_t index)
{
    return leftmul_raw(p, single(index, LIMDD_PAULI_Z));
}

bool
limdd_pauli_from_string(const char *str, size_t nqubits, limdd_pauli_t *out)
{
    if (str == NULL || out == NULL || nqubits > LIMDD_MAX_QUBITS) return false;
    if (strlen(str) != nqubits) return false;

    limdd_pauli_t p = limdd_pauli_identity();
    for (size_t i = 0; i < nqubits; i++) {
        switch (str[i]) {
            case 'I': break;
            case 'Z': limdd_pauli_set(&p, i, LIMDD_PAULI_Z); break;
            case 'X': limdd_pauli_set(&p, i, LIMDD_PAULI_X); break;
            case 'Y': limdd_pauli_set(&p, i, LIMDD_PAULI_Y); break;
            default: return false;
        }
    }
    *out = p;
    return true;
}

void
limdd_pauli_to_string(limdd_pauli_t p, size_t nqubits, char *buf)
{
    static const char chars[4] = { 'I', 'Z', 'X', 'Y' };
    assert(nqubits <= LIMDD_MAX_QUBITS);
    for (size_t i = 0; i < nqubits; i++) {
        buf[i] = chars[limdd_pauli_get(p, i)];
    }
    buf[nqubits] = '\0';
}

void
limdd_pauli_print(FILE *out, limdd_pauli_t p, size_t nqubits)
{
    char buf[LIMDD_MAX_QUBITS + 1];
    limdd_pauli_to_string(p, nqubits, buf);
    fputs(buf, out);
}
