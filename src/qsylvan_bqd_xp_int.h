/*
 * Copyright 2025 Q-Sylvan contributors
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
 * The labels and labelled edges of the translation and Pauli families, taken
 * out of the copy rule's operations (qsylvan_bqd_xp.c) so that another file of
 * operations on labelled edges can share them: a label c Z^s X^t read out of
 * its LIM, the label product and inverse of skip:sec:xp:edges, Norm, the top
 * of a labelled edge, a scale, the key word of a labelled edge in a memo, the
 * least point of a support, which reads the support shadow and no stored
 * ratio, and the conjugation of a label by a permutation of two qubits.
 * Private to the files of those families, and none of it reads the zero
 * rule. Every function is static inline but minpoint, which is static and not
 * inline: the copy rule's file called it as a static function before this
 * header took it, and inlining it there would change that rule's object code.
 */

#ifndef QSYLVAN_BQD_XP_INT_H
#define QSYLVAN_BQD_XP_INT_H

#include "qsylvan_bqd.h"
#include "qsylvan_bqd_gates.h"

/* --- labels ------------------------------------------------------------------ */

/**
 * A label c Z^s X^t, read out of its LIM. s is the LIM's z word and t its x
 * word, bit q for qubit q, in the first lane: a BQD has at most 63 qubits.
 */
typedef struct {
    EVBDD_WGT c;
    uint64_t  s, t;
} label_t;

static inline bool parity(uint64_t v) { return __builtin_parityll(v) != 0; }

/*
 * -c, as the product with -1, which the operation cache remembers. wgt_neg
 * works out the value and looks it up in the weight table every time, and the
 * signs of the Pauli family -- of a cofactor, of a label product, the test
 * b = -a and the sign repair -- took a fifth of the time of its gates that way.
 */
static inline EVBDD_WGT
neg(EVBDD_WGT c)
{
    return wgt_mul(EVBDD_MIN_ONE, c);
}

/** The qubits above variable v, which are the bits below bit v. */
static inline uint64_t
above(uint32_t v)
{
    return (UINT64_C(1) << v) - 1;
}

static inline label_t
label_read(LIMDD_LIM lim)
{
    label_t l = { EVBDD_ONE, 0, 0 };
    if (lim == LIMDD_LIM_IDENTITY) return l;
    if (limdd_lim_is_zero(lim)) { l.c = EVBDD_ZERO; return l; }
    const limdd_pauli_t p = limdd_lim_pauli(lim);
    l.c = limdd_lim_weight(lim);
    l.s = p.z[0];
    l.t = p.x[0];
    return l;
}

/** The LIM of a label, interned through the per-worker memo; the zero weight gives the zero label. */
static inline LIMDD_LIM
label_lim(label_t l)
{
    return bqd_lim_word(l.c, l.s, l.t);
}

/**
 * The label a . b, which acts as b and then a. X^t Z^s = (-1)^{s.t} Z^s X^t,
 * so moving b's Z past a's X costs the sign of s_b . t_a.
 */
static inline label_t
label_mul(label_t a, label_t b)
{
    label_t r = { wgt_mul(a.c, b.c), a.s ^ b.s, a.t ^ b.t };
    if (parity(b.s & a.t)) r.c = neg(r.c);
    return r;
}

/** The inverse: (c Z^s X^t)(c' Z^s X^t) = c c' (-1)^{s.t}, so c' = c^-1 (-1)^{s.t}. */
static inline label_t
label_inv(label_t l)
{
    label_t r = { wgt_div(EVBDD_ONE, l.c), l.s, l.t };
    if (parity(l.s & l.t)) r.c = neg(r.c);
    return r;
}

/* the label (1, 0, t), the translation X^t */
static inline label_t
xlabel(uint64_t t)
{
    const label_t l = { EVBDD_ONE, 0, t };
    return l;
}

/* --- labelled edges ---------------------------------------------------------- */

static inline BQD
unit(LIMDD_TARG t)
{
    return limdd_bundle(LIMDD_LIM_IDENTITY, t);
}

/**
 * The labelled edge l on N, normal: an X above the level of N acts on a
 * function that does not depend on that level, so Norm clears it
 * (skip:alg:xpedges). The zero scalar gives the zero edge.
 */
static inline BQD
xedge(label_t l, LIMDD_TARG N)
{
    if (l.c == EVBDD_ZERO) return limdd_zero_edge();
    l.t &= ~above(limdd_level(N));
    return limdd_bundle(label_lim(l), N);
}

/** The top of (l, N) as a variable: N's level, or the highest qubit with a Z above it. */
static inline uint32_t
top_var(LIMDD_TARG N, uint64_t s)
{
    const uint32_t v = limdd_level(N);
    if (s == 0) return v;
    const uint32_t z = (uint32_t)__builtin_ctzll(s);
    return z < v ? z : v;
}

/** The scalar of a labelled edge's label. */
static inline EVBDD_WGT
scalar_of(BQD e)
{
    return limdd_lim_weight(limdd_label(e));
}

/** The node of an edge, or 0 for the zero edge. */
static inline LIMDD_TARG
node_or_zero(BQD e)
{
    return limdd_edge_is_zero(e) ? 0 : limdd_target(e);
}

/** k . e, which is canonical when e is (skip:lem:xscale). */
static inline BQD
xscale(EVBDD_WGT k, BQD e)
{
    if (k == EVBDD_ZERO || limdd_edge_is_zero(e)) return limdd_zero_edge();
    if (k == EVBDD_ONE) return e;
    label_t l = label_read(limdd_label(e));
    l.c = wgt_mul(k, l.c);
    /* on floats a product of two small scalars can merge with zero, and a zero
     * label on a live node is not the zero edge: a division by it later gives
     * a NaN (as lim_times_edge in the LIMDD) */
    if (l.c == EVBDD_ZERO) return limdd_zero_edge();
    return limdd_bundle(label_lim(l), limdd_target(e));
}

/** c on the terminal, the constant c; zero for c = 0. */
static inline BQD
constant(EVBDD_WGT c)
{
    const label_t l = { c, 0, 0 };
    return xedge(l, LIMDD_TERMINAL);
}

/** label . e, normal and not canonical: the label product, which Canon makes canonical. */
static inline BQD
relabel(label_t l, BQD e)
{
    if (limdd_edge_is_zero(e)) return e;
    return xedge(label_mul(l, label_read(limdd_label(e))), limdd_target(e));
}

/**
 * A labelled edge in a memo key, without its scalar: the index of its label's
 * Pauli word and its node, in one word (bqd_init keeps the Pauli table within
 * 32 bits). Read from the LIM's bucket, so a key interns nothing. The word is
 * the edge's own, before Norm: an X above the node gives a second key for the
 * same function, which costs a miss and not a wrong result, and the edges the
 * recursion makes are normal already. An operation that misses pushes its
 * operands onto the refs stack, since it reads their labels again after a
 * point where a collection can run, and a caller may hand it an edge nothing
 * else holds (the reversed CX of bqd_cgate_either passes one gate's result
 * to the next); that keeps the word of the key as well.
 */
static inline uint64_t
key_edge(BQD e)
{
    return (limdd_lim_pauli_ref(limdd_label(e)) << 32) | limdd_target(e);
}

/**
 * The result r of an operation on edges with the scalar c taken out of the
 * key, as the memo keeps it: r / c, canonical when r is (skip:lem:xscale).
 * The recursion runs on the edges as they are, scalar and all, and so returns
 * the result to hand back without a scale; the memo pays the one scale
 * instead, and no key interns an edge without its scalar.
 */
static inline uint64_t
memo_of(BQD r, EVBDD_WGT c)
{
    return (uint64_t)xscale(wgt_div(EVBDD_ONE, c), r);
}

/* --- the least point --------------------------------------------------------- */

/**
 * The least point of supp(ext [N]) ^ tau, a mask of qubits (skip:lem:minpoint).
 * The support of a node is {0} x supp[N_0] together with {1} x (supp[N_1] ^ t_1),
 * the support shadow (lem:shadow), so the least point has a 0 at the node's
 * level exactly when the cofactor x = tau there is nonzero, and the walk goes
 * on in the cofactor it takes, moved by the translation that cofactor carries.
 * At a level an edge skips the support is closed under flipping the bit and
 * the least point has a 0 there, which the walk leaves in place. One path,
 * at most one step per level, and no memo.
 */
static uint64_t
minpoint(LIMDD_TARG N, uint64_t tau)
{
    uint64_t y = 0;
    while (N != LIMDD_TERMINAL) {
        const uint64_t bit = UINT64_C(1) << limdd_node_var(N);
        const LIMDD lo = limdd_node_low(N), hi = limdd_node_high(N);
        const bool tl = (tau & bit) != 0;
        const bool up = tl ? !limdd_edge_is_zero(hi) : limdd_edge_is_zero(lo);  /* x = 1 */
        if (up != tl) y |= bit;
        tau &= ~bit;
        if (up) {
            tau ^= label_read(limdd_label(hi)).t;
            N = limdd_target(hi);
        } else {
            N = limdd_target(lo);
        }
    }
    return y;
}

/* --- permutations of two qubits ---------------------------------------------- */

/** pi on a mask of qubits, (pi m) for the permutation `kind` of qa < qb. */
static inline uint64_t
perm_mask(uint32_t kind, uint32_t qa, uint32_t qb, uint64_t m)
{
    const uint64_t a = (m >> qa) & 1, b = (m >> qb) & 1;
    uint64_t na = a, nb = b;
    switch (kind) {
    case BQD_PERM_SWAP:    na = b; nb = a; break;
    case BQD_PERM_CX_DOWN: nb = b ^ a;     break;     /* control qa, target qb */
    default:               na = a ^ b;     break;     /* control qb, target qa */
    }
    m &= ~((UINT64_C(1) << qa) | (UINT64_C(1) << qb));
    return m | (na << qa) | (nb << qb);
}

/**
 * The conjugated label l^pi = c Z^{pi^T s} X^{pi t}, pi . (l . g) = l^pi . pi g:
 * s . pi x = pi^T s . x, and pi x xor t = pi (x xor pi t). The swap is its own
 * transpose, and the CX with the control above is the transpose of the one
 * with the control below.
 */
static inline label_t
perm_label(uint32_t kind, uint32_t qa, uint32_t qb, label_t l)
{
    const uint32_t tk = kind == BQD_PERM_SWAP ? kind
                      : kind == BQD_PERM_CX_DOWN ? (uint32_t)BQD_PERM_CX_UP : (uint32_t)BQD_PERM_CX_DOWN;
    l.s = perm_mask(tk, qa, qb, l.s);
    l.t = perm_mask(kind, qa, qb, l.t);
    return l;
}

#endif
