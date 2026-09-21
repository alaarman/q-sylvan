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

#include <sylvan_platform.h>   /* popcnt_uint64 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Maximum number of qubits a single packed Pauli word can describe.
 */
#ifndef LIMDD_PAULI_WORDS
/* 64-bit words per component. Set with -DLIMDD_PAULI_WORDS=k, or through the
 * CMake cache variable of the same name. One is the historical width. */
#define LIMDD_PAULI_WORDS 1
#endif

#if LIMDD_PAULI_WORDS < 1
#error "LIMDD_PAULI_WORDS must be at least 1"
#endif

#define LIMDD_MAX_QUBITS (64 * LIMDD_PAULI_WORDS)

/** Which word of a component qubit `i` lives in, and which bit of it. */
#define LIMDD_PAULI_LANE(i) ((LIMDD_PAULI_WORDS == 1) ? 0u : (unsigned)((i) >> 6))
#define LIMDD_PAULI_BIT(i)  (UINT64_C(1) << ((i) & 63))

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
    // Lane i holds qubits 64i..64i+63. Bit j of x is set iff that qubit
    // carries an X component (X or Y); likewise z for a Z component.
    uint64_t x[LIMDD_PAULI_WORDS];
    uint64_t z[LIMDD_PAULI_WORDS];
} limdd_pauli_t;

/**
 * The all-identity word "II..I".
 */
static inline limdd_pauli_t
limdd_pauli_identity(void)
{
    limdd_pauli_t p;
    for (unsigned i = 0; i < LIMDD_PAULI_WORDS; i++) { p.x[i] = 0; p.z[i] = 0; }
    return p;
}

/**
 * True iff `p` is the identity on every qubit.
 */
static inline bool
limdd_pauli_is_identity(limdd_pauli_t p)
{
    uint64_t acc = 0;
    for (unsigned i = 0; i < LIMDD_PAULI_WORDS; i++) acc |= p.x[i] | p.z[i];
    return acc == 0;
}

/**
 * True iff `a` and `b` are the same Pauli word (ignoring phase).
 */
static inline bool
limdd_pauli_equals(limdd_pauli_t a, limdd_pauli_t b)
{
    uint64_t neq = 0;
    for (unsigned i = 0; i < LIMDD_PAULI_WORDS; i++)
        neq |= (a.x[i] ^ b.x[i]) | (a.z[i] ^ b.z[i]);
    return neq == 0;
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
    const unsigned w = LIMDD_PAULI_LANE(index);
    const uint64_t b = LIMDD_PAULI_BIT(index);
    return (limdd_pauli_op_t)(2u * ((p.x[w] & b) != 0) + ((p.z[w] & b) != 0));
}

/**
 * Set the Pauli on qubit `index` to `op`, leaving every other qubit alone.
 */
static inline void
limdd_pauli_set(limdd_pauli_t *p, size_t index, limdd_pauli_op_t op)
{
    const unsigned w = LIMDD_PAULI_LANE(index);
    const uint64_t bit = LIMDD_PAULI_BIT(index);
    p->x[w] = (op & 2u) ? (p->x[w] | bit) : (p->x[w] & ~bit);
    p->z[w] = (op & 1u) ? (p->z[w] | bit) : (p->z[w] & ~bit);
}

/**
 * Number of qubits on which `p` is not the identity.
 */
size_t limdd_pauli_weight(limdd_pauli_t p);

/* --- the idioms that would otherwise reach for .x and .z ------------------
 *
 * Everything below exists so that no other file has to know how a Pauli word
 * is laid out. Written over the current one-word-per-component form, they
 * compile to exactly the expressions they replace; when the components become
 * arrays these are the only bodies that change.
 */

/** A Pauli that is `op` on `index` and the identity elsewhere. */
static inline limdd_pauli_t
limdd_pauli_single(size_t index, limdd_pauli_op_t op)
{
    limdd_pauli_t p = limdd_pauli_identity();
    limdd_pauli_set(&p, index, op);
    return p;
}

/** a ^= b, word by word. The Pauli product's phase is NOT applied; this is
 *  the symplectic sum, which is what GF(2) elimination needs. */
static inline void
limdd_pauli_xor(limdd_pauli_t *a, limdd_pauli_t b)
{
    for (unsigned i = 0; i < LIMDD_PAULI_WORDS; i++) {
        a->x[i] ^= b.x[i];
        a->z[i] ^= b.z[i];
    }
}

/**
 * Bit `c` of the symplectic vector: the X half occupies columns 0..nqubits-1
 * and the Z half the rest. Elimination orders columns this way, so the X half
 * leading is what makes a reduced label carry I or Z on a skipped level
 * rather than X or Y.
 */
static inline bool
limdd_pauli_column(limdd_pauli_t p, size_t c, size_t nqubits)
{
    const size_t k = (c < nqubits) ? c : c - nqubits;
    const uint64_t w = (c < nqubits) ? p.x[LIMDD_PAULI_LANE(k)]
                                     : p.z[LIMDD_PAULI_LANE(k)];
    return (w & LIMDD_PAULI_BIT(k)) != 0;
}

/** True iff `p` acts on any qubit below `level`. */
static inline bool
limdd_pauli_acts_below(limdd_pauli_t p, size_t level)
{
    if (level == 0) return false;
    uint64_t acc = 0;
    for (unsigned i = 0; i < LIMDD_PAULI_WORDS; i++) {
        const size_t base = (size_t)i * 64;
        if (level <= base) break;
        const size_t k = level - base;
        const uint64_t m = (k >= 64) ? ~UINT64_C(0) : (UINT64_C(1) << k) - 1;
        acc |= (p.x[i] | p.z[i]) & m;
    }
    return acc != 0;
}

/**
 * Split `p` at `level`: the part acting below it is returned, and `p` keeps
 * the rest. The two have disjoint support, so neither the split nor putting
 * them back together picks up a phase.
 */
static inline limdd_pauli_t
limdd_pauli_split_below(limdd_pauli_t *p, size_t level)
{
    limdd_pauli_t out = limdd_pauli_identity();
    if (level == 0) return out;
    for (unsigned i = 0; i < LIMDD_PAULI_WORDS; i++) {
        const size_t base = (size_t)i * 64;
        if (level <= base) break;
        const size_t k = level - base;
        const uint64_t m = (k >= 64) ? ~UINT64_C(0) : (UINT64_C(1) << k) - 1;
        out.x[i] = p->x[i] & m;
        out.z[i] = p->z[i] & m;
        p->x[i] &= ~m;
        p->z[i] &= ~m;
    }
    return out;
}

/**
 * A total order on Pauli words, by value rather than by any interned index.
 * Comparing indices would make the order depend on what the program happened
 * to intern first, and on whether a collection has since moved anything.
 */
static inline int
limdd_pauli_cmp(limdd_pauli_t a, limdd_pauli_t b)
{
    /* Most significant lane first, so widening does not change the order. */
    for (unsigned i = LIMDD_PAULI_WORDS; i-- > 0; )
        if (a.x[i] != b.x[i]) return a.x[i] < b.x[i] ? -1 : 1;
    for (unsigned i = LIMDD_PAULI_WORDS; i-- > 0; )
        if (a.z[i] != b.z[i]) return a.z[i] < b.z[i] ? -1 : 1;
    return 0;
}

/**
 * Apply `p` to the basis state `*b`, returning the power of i it contributes.
 *
 * P sends |c> to i^(x&z) (-1)^(z&c) |c XOR x>, so `*b` is advanced to the
 * basis state that actually reaches the target and the phase comes back as a
 * power of i, mod 4.
 */
static inline unsigned
limdd_pauli_apply_basis(limdd_pauli_t p, uint64_t *b)
{
    /* The basis state is a single word, so this serves the evaluators only,
     * which are exponential and never run above a few tens of qubits. */
    const uint64_t c = *b ^ p.x[0];
    *b = c;
    return (unsigned)(popcnt_uint64(p.x[0] & p.z[0])
                      + 2u * popcnt_uint64(p.z[0] & c)) & 3u;
}

/** Serialise into, and load from, the two words of a table entry. */
static inline void
limdd_pauli_store(limdd_pauli_t p, uint64_t *w)
{
    for (unsigned i = 0; i < LIMDD_PAULI_WORDS; i++) {
        w[2 * i]     = p.x[i];
        w[2 * i + 1] = p.z[i];
    }
}

static inline limdd_pauli_t
limdd_pauli_load(const uint64_t *w)
{
    limdd_pauli_t p;
    for (unsigned i = 0; i < LIMDD_PAULI_WORDS; i++) {
        p.x[i] = w[2 * i];
        p.z[i] = w[2 * i + 1];
    }
    return p;
}

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
