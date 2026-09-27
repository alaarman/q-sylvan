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

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "qsylvan_bqd.h"
#include "edge_weight_storage/qisq2_map.h"

/* Written once by bqd_init and read afterwards. */
static bqd_family_t family  = BQD_FAMILY_SCALAR;
static size_t       nqubits = 0;

void
bqd_init(bqd_family_t f, size_t n, size_t node_tablesize,
         size_t pauli_tablesize, size_t lim_tablesize, size_t stab_tablesize)
{
    if (n == 0 || n > 63) {
        fprintf(stderr, "sylvan: a BQD indexes its vectors by a 64-bit word, so "
                        "it holds 1 to 63 qubits, not %zu\n", n);
        exit(1);
    }
    family  = f;
    nqubits = n;
    limdd_nodes_init(n, node_tablesize, pauli_tablesize, lim_tablesize, stab_tablesize);
}

void
bqd_quit(void)
{
    limdd_nodes_quit();
    nqubits = 0;
}

bqd_family_t bqd_family(void) { return family; }

const char *
bqd_family_name(bqd_family_t f)
{
    switch (f) {
    case BQD_FAMILY_SCALAR: return "BQD";
    case BQD_FAMILY_X:      return "X-BQD";
    case BQD_FAMILY_PAULI:  return "Pauli-BQD";
    }
    return "?";
}

/* --- small helpers --------------------------------------------------------- */

static inline bool parity(uint64_t v) { return __builtin_parityll(v) != 0; }

/** Index of the least point of the support, or -1 for the zero vector. */
static int64_t
lexmin(const EVBDD_WGT *g, uint64_t len)
{
    for (uint64_t y = 0; y < len; y++) if (g[y] != EVBDD_ZERO) return (int64_t)y;
    return -1;
}

/**
 * Whether arg(w) lies in [0, pi): Im w > 0, or Im w = 0 and Re w > 0.
 *
 * Exact for the algebraic backend. A weight there is a + b sqrt2 + i(c + d
 * sqrt2) with rational a, b, c, d, and the sign of A + B sqrt2 is decided
 * without a root: when A and B agree in sign it is theirs, and otherwise it
 * is the sign of whichever of A^2 and 2B^2 is larger, which cannot tie for
 * rationals since sqrt2 is irrational.
 */
static int
sign_sqrt2(const mpq_t A, const mpq_t B)
{
    const int sa = mpq_sgn(A), sb = mpq_sgn(B);
    if (sa == 0 && sb == 0) return 0;
    if (sa >= 0 && sb >= 0) return 1;
    if (sa <= 0 && sb <= 0) return -1;
    mpq_t a2, b2;
    mpq_init(a2); mpq_init(b2);
    mpq_mul(a2, A, A);
    mpq_mul(b2, B, B);
    mpq_add(b2, b2, b2);                     /* 2 B^2 */
    const int c = mpq_cmp(a2, b2);
    mpq_clear(a2); mpq_clear(b2);
    return c > 0 ? sa : sb;
}

static bool
arg_in_upper(EVBDD_WGT w)
{
    if (sylvan_get_edge_weight_type() == WGT_QISQ2) {
        const qisq2_t *q = (const qisq2_t *)qisq2_map_get(wgt_storage, (uint64_t)w);
        const int im = sign_sqrt2(q->c, q->d);
        if (im != 0) return im > 0;
        return sign_sqrt2(q->a, q->b) > 0;
    }
    const complex_t z = weight_as_complex(w);
    return z.i > 0.0 || (z.i == 0.0 && z.r > 0.0);
}

/* --- labels --------------------------------------------------------------- */

LIMDD_LIM
bqd_lim_make(EVBDD_WGT c, uint64_t s, uint64_t t, uint32_t n)
{
    assert(n <= nqubits);
    assert(family != BQD_FAMILY_SCALAR || (s == 0 && t == 0));
    assert(family != BQD_FAMILY_X || s == 0);
    if (c == EVBDD_ZERO) return LIMDD_LIM_ZERO;
    limdd_pauli_t p;
    memset(&p, 0, sizeof(p));
    for (uint32_t q = 0; q < n; q++) {
        const uint32_t b = n - 1 - q;                 /* qubit q at bit n-1-q */
        if ((t >> b) & 1) p.x[LIMDD_PAULI_LANE(q)] |= LIMDD_PAULI_BIT(q);
        if ((s >> b) & 1) p.z[LIMDD_PAULI_LANE(q)] |= LIMDD_PAULI_BIT(q);
    }
    return limdd_lim_make(p, c);
}

void
bqd_lim_masks(LIMDD_LIM lim, uint32_t n, EVBDD_WGT *c, uint64_t *s, uint64_t *t)
{
    *s = *t = 0;
    if (limdd_lim_is_zero(lim)) { *c = EVBDD_ZERO; return; }
    *c = limdd_lim_weight(lim);
    const limdd_pauli_t p = limdd_lim_pauli(lim);
    for (uint32_t q = 0; q < n; q++) {
        const uint32_t b = n - 1 - q;
        if (p.x[LIMDD_PAULI_LANE(q)] & LIMDD_PAULI_BIT(q)) *t |= (uint64_t)1 << b;
        if (p.z[LIMDD_PAULI_LANE(q)] & LIMDD_PAULI_BIT(q)) *s |= (uint64_t)1 << b;
    }
}

/** The label that recovers  g = (c0 Z^s X^t) . rep  read as c Z^s X^t. */
static LIMDD_LIM
label_of(EVBDD_WGT c0, uint64_t s, uint64_t t, uint32_t n)
{
    const EVBDD_WGT c = parity(s & t) ? wgt_neg(c0) : c0;
    return bqd_lim_make(c, s, t, n);
}

/* --- the representative --------------------------------------------------- */

/**
 * The sign vector of def:prep, by back substitution.
 *
 * The pivots of g (which has 0 in its support) are one point per level down
 * the low edges: at the level deciding bit b, the least point of the high
 * cofactor with that bit set. The pivot of a level is the only pivot with a 1
 * at that bit, so reading the levels from the deepest up fixes one bit of s
 * at a time from the bits already fixed below it.
 */
static uint64_t
sign_vector(const EVBDD_WGT *g2, uint64_t len)
{
    uint64_t s = 0;
    int nb = 0;
    while ((UINT64_C(1) << nb) < len) nb++;           /* len = 2^nb */
    int64_t piv[64];
    for (int b = 0; b < nb; b++) {
        const uint64_t h = UINT64_C(1) << b;
        const int64_t m = lexmin(g2 + h, h);           /* high half of the prefix [0, 2h) */
        piv[b] = (m < 0) ? -1 : (int64_t)(h | (uint64_t)m);
    }
    for (int b = 0; b < nb; b++) {                     /* deepest level first */
        if (piv[b] < 0) continue;
        const uint64_t q = (uint64_t)piv[b];
        EVBDD_WGT v = g2[q];
        if (parity(s & q)) v = wgt_neg(v);
        if (!arg_in_upper(v)) s |= UINT64_C(1) << b;
    }
    return s;
}

void
bqd_representative(const EVBDD_WGT *g, uint64_t len, EVBDD_WGT *rep,
                   EVBDD_WGT *c0, uint64_t *s, uint64_t *p)
{
    const int64_t lm = lexmin(g, len);
    if (lm < 0) {
        for (uint64_t y = 0; y < len; y++) rep[y] = EVBDD_ZERO;
        *c0 = EVBDD_ONE; *s = 0; *p = 0;
        return;
    }
    const uint64_t pp = (family != BQD_FAMILY_SCALAR) ? (uint64_t)lm : 0;

    /* translate, so that the least point of the support sits at 0 */
    EVBDD_WGT *tmp = (rep == g) ? malloc(len * sizeof(EVBDD_WGT)) : rep;
    if (tmp == NULL) { fprintf(stderr, "sylvan: out of memory in a BQD build\n"); exit(1); }
    for (uint64_t y = 0; y < len; y++) tmp[y] = g[y ^ pp];

    /* scale, so that the first value is 1 */
    const EVBDD_WGT cc = tmp[lexmin(tmp, len)];
    for (uint64_t y = 0; y < len; y++)
        if (tmp[y] != EVBDD_ZERO) tmp[y] = wgt_div(tmp[y], cc);

    /* and for the Pauli family, the sign pattern that fixes every pivot */
    const uint64_t ss = (family == BQD_FAMILY_PAULI) ? sign_vector(tmp, len) : 0;
    if (ss != 0)
        for (uint64_t y = 0; y < len; y++)
            if (tmp[y] != EVBDD_ZERO && parity(ss & y)) tmp[y] = wgt_neg(tmp[y]);

    if (tmp != rep) { memcpy(rep, tmp, len * sizeof(EVBDD_WGT)); free(tmp); }
    *c0 = cc; *s = ss; *p = pp;
}

/* --- building ------------------------------------------------------------- */

/**
 * The node of a representative `g` at level `var`, of length 2^(n-var).
 *
 * The three steps of sec:canon. The low cofactor of a representative is its
 * own representative (its least support point is 0 with value 1, and its
 * pivots are among the parent's), so the low edge carries the identity. The
 * high edge gets the translation that aligns the least points of the two
 * cofactor supports, then the pointwise quotient with the copy where the low
 * cofactor is zero, then the representative of that, whose recovering label
 * is what the edge carries.
 */
TASK_3(LIMDD_TARG, bqd_build, const EVBDD_WGT *, g, uint32_t, var, uint32_t, n)
{
    if (var == n) return LIMDD_TERMINAL;             /* g = [1] */

    const uint64_t h = UINT64_C(1) << (n - var - 1);
    const EVBDD_WGT *f0 = g, *f1 = g + h;

    const int64_t s0 = lexmin(f0, h);
    const int64_t s1 = lexmin(f1, h);
    const bool low_zero = (s0 < 0);

    /* the low child and the high child are independent: overlap them */
    if (!low_zero) SPAWN(bqd_build, f0, var + 1, n);

    LIMDD high = limdd_zero_edge();
    if (s1 >= 0) {
        /* step 2: align the least support points, when the family can move one */
        const uint64_t t = (family != BQD_FAMILY_SCALAR && !low_zero)
                         ? (uint64_t)(s1 ^ s0) : 0;

        /* step 3: the pointwise quotient, copying where there is no divisor */
        EVBDD_WGT *raw = malloc(h * sizeof(EVBDD_WGT));
        if (raw == NULL) { fprintf(stderr, "sylvan: out of memory in a BQD build\n"); exit(1); }
        for (uint64_t y = 0; y < h; y++) {
            const EVBDD_WGT hy = f1[y ^ t];
            raw[y] = (f0[y] == EVBDD_ZERO || hy == EVBDD_ZERO) ? hy : wgt_div(hy, f0[y]);
        }

        EVBDD_WGT c0; uint64_t s, p;
        bqd_representative(raw, h, raw, &c0, &s, &p);
        /* the quotient is nonzero at 0 once the least points are aligned, and
         * the scalar family may not translate, so nothing is left to move */
        assert(p == 0);

        const LIMDD_TARG ht = CALL(bqd_build, raw, var + 1, n);
        high = limdd_bundle(label_of(c0, s, t, n), ht);
        free(raw);
    }

    const LIMDD low = low_zero ? limdd_zero_edge()
                               : limdd_bundle(LIMDD_LIM_IDENTITY, SYNC(bqd_build));
    return limdd_makenode(var, low, high);
}

TASK_IMPL_2(BQD, bqd_from_vector, const EVBDD_WGT *, f, uint32_t, n)
{
    assert(n >= 1 && n <= nqubits);
    const uint64_t len = UINT64_C(1) << n;
    if (lexmin(f, len) < 0) return limdd_zero_edge();

    EVBDD_WGT *rep = malloc(len * sizeof(EVBDD_WGT));
    if (rep == NULL) { fprintf(stderr, "sylvan: out of memory in a BQD build\n"); exit(1); }
    EVBDD_WGT c0; uint64_t s, p;
    bqd_representative(f, len, rep, &c0, &s, &p);
    const LIMDD_TARG root = CALL(bqd_build, rep, 0, n);
    free(rep);
    return limdd_bundle(label_of(c0, s, p, n), root);
}

/* --- decoding ------------------------------------------------------------- */

static void
fill_zero(EVBDD_WGT *out, uint64_t len)
{
    for (uint64_t y = 0; y < len; y++) out[y] = EVBDD_ZERO;
}

/** f_1 = (c Z^s X^t) . (g_0 (.) g_1), written into out[h..2h) from g_0 in out[0..h). */
static void
combine_high(EVBDD_WGT *out, const EVBDD_WGT *g1, uint64_t h, LIMDD_LIM lim, uint32_t n)
{
    EVBDD_WGT c; uint64_t s, t;
    bqd_lim_masks(lim, n, &c, &s, &t);
    for (uint64_t y = 0; y < h; y++) {
        const uint64_t yt = y ^ t;
        const EVBDD_WGT a = out[yt], b = g1[yt];
        EVBDD_WGT v = (a == EVBDD_ZERO) ? b : (b == EVBDD_ZERO ? EVBDD_ZERO : wgt_mul(a, b));
        if (v != EVBDD_ZERO) {
            if (parity(s & y)) v = wgt_neg(v);
            v = wgt_mul(c, v);
        }
        out[h + y] = v;
    }
}

VOID_TASK_4(bqd_decode, LIMDD_TARG, p, uint32_t, var, uint32_t, n, EVBDD_WGT *, out)
{
    if (var == n) { out[0] = EVBDD_ONE; return; }
    assert(p != LIMDD_TERMINAL && "a BQD never skips a level");
    assert(limdd_node_var(p) == var);

    const uint64_t h = UINT64_C(1) << (n - var - 1);
    const LIMDD low = limdd_node_low(p), high = limdd_node_high(p);
    const bool lz = limdd_edge_is_zero(low), hz = limdd_edge_is_zero(high);

    if (!lz) SPAWN(bqd_decode, limdd_target(low), var + 1, n, out);

    EVBDD_WGT *g1 = NULL;
    if (!hz) {
        g1 = malloc(h * sizeof(EVBDD_WGT));
        if (g1 == NULL) { fprintf(stderr, "sylvan: out of memory decoding a BQD\n"); exit(1); }
        CALL(bqd_decode, limdd_target(high), var + 1, n, g1);
    }

    if (!lz) SYNC(bqd_decode); else fill_zero(out, h);

    if (hz) fill_zero(out + h, h);
    else { combine_high(out, g1, h, limdd_label(high), n); free(g1); }
}

VOID_TASK_IMPL_3(bqd_to_vector, BQD, e, uint32_t, n, EVBDD_WGT *, out)
{
    assert(n >= 1 && n <= nqubits);
    const uint64_t len = UINT64_C(1) << n;
    if (limdd_edge_is_zero(e)) { fill_zero(out, len); return; }

    EVBDD_WGT *rep = malloc(len * sizeof(EVBDD_WGT));
    if (rep == NULL) { fprintf(stderr, "sylvan: out of memory decoding a BQD\n"); exit(1); }
    CALL(bqd_decode, limdd_target(e), 0, n, rep);

    EVBDD_WGT c; uint64_t s, t;
    bqd_lim_masks(limdd_label(e), n, &c, &s, &t);
    for (uint64_t y = 0; y < len; y++) {
        EVBDD_WGT v = rep[y ^ t];
        if (v != EVBDD_ZERO) {
            if (parity(s & y)) v = wgt_neg(v);
            v = wgt_mul(c, v);
        }
        out[y] = v;
    }
    free(rep);
}

/* --- one amplitude -------------------------------------------------------- */

static EVBDD_WGT
eval_node(LIMDD_TARG p, uint32_t var, uint32_t n, uint64_t y)
{
    if (var == n) return EVBDD_ONE;
    const uint64_t h = UINT64_C(1) << (n - var - 1);
    const LIMDD low = limdd_node_low(p);
    if ((y & h) == 0)
        return limdd_edge_is_zero(low) ? EVBDD_ZERO : eval_node(limdd_target(low), var + 1, n, y);

    const LIMDD high = limdd_node_high(p);
    if (limdd_edge_is_zero(high)) return EVBDD_ZERO;
    EVBDD_WGT c; uint64_t s, t;
    bqd_lim_masks(limdd_label(high), n, &c, &s, &t);

    const uint64_t yy = y & (h - 1), yt = yy ^ t;
    /* the fork: the high cofactor at yy is the product of both children at yt */
    const EVBDD_WGT a = limdd_edge_is_zero(low) ? EVBDD_ZERO
                                                : eval_node(limdd_target(low), var + 1, n, yt);
    const EVBDD_WGT b = eval_node(limdd_target(high), var + 1, n, yt);
    EVBDD_WGT v = (a == EVBDD_ZERO) ? b : (b == EVBDD_ZERO ? EVBDD_ZERO : wgt_mul(a, b));
    if (v == EVBDD_ZERO) return v;
    if (parity(s & yy)) v = wgt_neg(v);
    return wgt_mul(c, v);
}

EVBDD_WGT
bqd_eval(BQD e, uint32_t n, uint64_t x)
{
    assert(n >= 1 && n <= nqubits);
    if (limdd_edge_is_zero(e)) return EVBDD_ZERO;
    EVBDD_WGT c; uint64_t s, t;
    bqd_lim_masks(limdd_label(e), n, &c, &s, &t);
    EVBDD_WGT v = eval_node(limdd_target(e), 0, n, x ^ t);
    if (v == EVBDD_ZERO) return v;
    if (parity(s & x)) v = wgt_neg(v);
    return wgt_mul(c, v);
}
