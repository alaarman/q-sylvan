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
 * Packed Pauli words, for LIMDD edge labels.
 *
 * A LIMDD labels each edge with a Linearly Invertible Map; for Pauli-LIMDDs
 * that map is a scalar times a tensor product of Pauli matrices. This file
 * implements the Pauli factor. The scalar is an EVBDD_WGT and lives in the
 * edge weight table, so it is not represented here.
 *
 * REPRESENTATION
 *
 * A Pauli word is held in symplectic form as two bit masks, one word each:
 * bit i of `x` is the X-component of qubit i and bit i of `z` its
 * Z-component. The Pauli on qubit i is then
 *
 *     2 * x_i + z_i   in   { I = 0, Z = 1, X = 2, Y = 3 }
 *
 * This packs one qubit into two bits instead of the two bytes a
 * `bool[2n]` array would use, and turns multiplication into an XOR of the
 * masks plus a phase correction computed with popcount. Both matter: the
 * word is part of every LIMDD edge label, so it is copied and compared on
 * every node lookup.
 *
 * The phase of a Pauli word is NOT stored here. Multiplication returns the
 * power of i that the caller must fold into the accompanying scalar, which
 * is what keeps this type a plain value: two words, no allocation, no
 * ownership, safe to copy between Lace workers. Nothing in this file has
 * mutable shared state, so it introduces no false sharing.
 *
 * LIMITS
 *
 * One word per component caps a Pauli word at LIMDD_MAX_QUBITS (64) qubits.
 * Bits at or above `n` are required to be zero for every operation here;
 * `limdd_pauli_is_canonical()` checks that invariant.
 */

#ifndef QSYLVAN_LIMDD_PAULI_H
#define QSYLVAN_LIMDD_PAULI_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Maximum number of qubits a single packed Pauli word can describe.
 */
#define LIMDD_MAX_QUBITS 64

/**
 * The four Pauli matrices, encoded as 2*x + z.
 *
 * NOTE the ordering: Z is 1 and X is 2, not the other way round. This
 * matches the symplectic packing directly and is worth keeping in mind when
 * reading code that switches on these values.
 */
typedef enum {
    LIMDD_PAULI_I = 0,
    LIMDD_PAULI_Z = 1,
    LIMDD_PAULI_X = 2,
    LIMDD_PAULI_Y = 3,
} limdd_pauli_op_t;

/**
 * A Pauli word on at most LIMDD_MAX_QUBITS qubits, without phase.
 */
typedef struct {
    uint64_t x; // bit i set iff qubit i carries an X component (X or Y)
    uint64_t z; // bit i set iff qubit i carries a Z component (Z or Y)
} limdd_pauli_t;

/**
 * The all-identity word "II..I".
 */
static inline limdd_pauli_t
limdd_pauli_identity(void)
{
    limdd_pauli_t p = { 0, 0 };
    return p;
}

/**
 * True iff `p` is the identity on every qubit.
 */
static inline bool
limdd_pauli_is_identity(limdd_pauli_t p)
{
    return (p.x | p.z) == 0;
}

/**
 * True iff `a` and `b` are the same Pauli word (ignoring phase).
 */
static inline bool
limdd_pauli_equals(limdd_pauli_t a, limdd_pauli_t b)
{
    return a.x == b.x && a.z == b.z;
}

/**
 * True iff no bit at or above `nqubits` is set, i.e. `p` is a well-formed
 * word on `nqubits` qubits. Every other function in this file assumes this.
 */
bool limdd_pauli_is_canonical(limdd_pauli_t p, size_t nqubits);

/**
 * The Pauli on qubit `index`, as a limdd_pauli_op_t.
 */
static inline limdd_pauli_op_t
limdd_pauli_get(limdd_pauli_t p, size_t index)
{
    return (limdd_pauli_op_t)(2u * ((p.x >> index) & 1u) + ((p.z >> index) & 1u));
}

/**
 * Set the Pauli on qubit `index` to `op`, leaving every other qubit alone.
 */
static inline void
limdd_pauli_set(limdd_pauli_t *p, size_t index, limdd_pauli_op_t op)
{
    const uint64_t bit = UINT64_C(1) << index;
    p->x = (op & 2u) ? (p->x | bit) : (p->x & ~bit);
    p->z = (op & 1u) ? (p->z | bit) : (p->z & ~bit);
}

/**
 * Number of qubits on which `p` is not the identity.
 */
size_t limdd_pauli_weight(limdd_pauli_t p);

/**
 * Right-multiply: replaces `*a` with the word `c` where a * b = i^k * c,
 * and returns k in {0,1,2,3}.
 *
 * The caller is responsible for multiplying the accompanying scalar by i^k;
 * this function does not track phase itself.
 */
unsigned limdd_pauli_rightmul(limdd_pauli_t *a, limdd_pauli_t b, size_t nqubits);

/**
 * Left-multiply: replaces `*a` with the word `c` where b * a = i^k * c,
 * and returns k in {0,1,2,3}.
 */
unsigned limdd_pauli_leftmul(limdd_pauli_t *a, limdd_pauli_t b, size_t nqubits);

/**
 * The power of i by which `a` and `b` fail to commute: returns 0 when they
 * commute and 2 when they anticommute, so that b * a = i^k * (a * b).
 *
 * This is the symplectic inner product and is much cheaper than performing
 * both multiplications.
 */
unsigned limdd_pauli_commutation_phase(limdd_pauli_t a, limdd_pauli_t b);

/**
 * True iff `a` and `b` commute.
 */
static inline bool
limdd_pauli_commutes(limdd_pauli_t a, limdd_pauli_t b)
{
    return limdd_pauli_commutation_phase(a, b) == 0;
}

/*
 * There is deliberately no inverse function. Every Pauli matrix is its own
 * inverse and squares to I with no leftover phase, and that carries over to
 * tensor products, so the inverse of a packed word is the word itself. Only
 * the accompanying scalar needs inverting.
 */

/**
 * Right-multiply `*p` by a single X, Y or Z on qubit `index`, returning the
 * power of i as limdd_pauli_rightmul does. Gate application conjugates one
 * qubit at a time, so these are the common case.
 */
unsigned limdd_pauli_rightmul_x(limdd_pauli_t *p, size_t index);
unsigned limdd_pauli_rightmul_y(limdd_pauli_t *p, size_t index);
unsigned limdd_pauli_rightmul_z(limdd_pauli_t *p, size_t index);

/**
 * Left-multiply `*p` by a single X, Y or Z on qubit `index`.
 */
unsigned limdd_pauli_leftmul_x(limdd_pauli_t *p, size_t index);
unsigned limdd_pauli_leftmul_y(limdd_pauli_t *p, size_t index);
unsigned limdd_pauli_leftmul_z(limdd_pauli_t *p, size_t index);

/**
 * Parse `str`, which must be `nqubits` characters from {I,X,Y,Z}, into a
 * packed word. Character 0 of the string is qubit 0.
 *
 * Returns false and leaves `*out` untouched if the string is malformed.
 */
bool limdd_pauli_from_string(const char *str, size_t nqubits, limdd_pauli_t *out);

/**
 * Write `p` as `nqubits` characters from {I,X,Y,Z} plus a NUL into `buf`,
 * which must hold at least nqubits+1 bytes. Character 0 is qubit 0.
 */
void limdd_pauli_to_string(limdd_pauli_t p, size_t nqubits, char *buf);

/**
 * Print `p` to `out` in I/X/Y/Z notation, without a trailing newline.
 */
void limdd_pauli_print(FILE *out, limdd_pauli_t p, size_t nqubits);

#ifdef __cplusplus
}
#endif

#endif // QSYLVAN_LIMDD_PAULI_H
