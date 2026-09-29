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

#include "qsylvan_bqd_xp.h"
#include "qsylvan_bqd_gates.h"
#include "qsylvan_gates.h"
#include "qsylvan_limdd_gc.h"

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

/** The LIM of a label, interned; the zero weight gives the zero label. */
static inline LIMDD_LIM
label_lim(label_t l)
{
    if (l.c == EVBDD_ZERO) return LIMDD_LIM_ZERO;
    if (l.c == EVBDD_ONE && l.s == 0 && l.t == 0) return LIMDD_LIM_IDENTITY;
    limdd_pauli_t p = limdd_pauli_identity();
    p.x[0] = l.t;
    p.z[0] = l.s;
    return limdd_lim_make(p, l.c);
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

uint32_t
bqd_xp_top(BQD e)
{
    if (limdd_edge_is_zero(e)) return limdd_level(LIMDD_TERMINAL);
    return top_var(limdd_target(e), label_read(limdd_label(e)).s);
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

BQD
bqd_xp_scale(BQD e, EVBDD_WGT c)
{
    return xscale(c, e);
}

/** c on the terminal, the constant c; zero for c = 0. */
static inline BQD
constant(EVBDD_WGT c)
{
    const label_t l = { c, 0, 0 };
    return xedge(l, LIMDD_TERMINAL);
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

/* --- cofactors and the least point ------------------------------------------- */

/**
 * The canonical edge of [N_0] (.) [N_1], memoised on N. Neither family has a
 * zero low edge, since 0 is in the support of every stored node, so the copy
 * everywhere is here only for completeness. The high label is not applied to
 * the high child inside the product, as the scalar family's Cof1 does: a
 * scalar and a sign pattern commute with (.), and a translation does not.
 */
TASK_IMPL_1(BQD, bqd_xp_join, LIMDD_TARG, N)
{
    const BQD lo = limdd_node_low(N), hi = limdd_node_high(N);
    if (limdd_edge_is_zero(hi)) return limdd_zero_edge();
    if (limdd_edge_is_zero(lo)) return unit(limdd_target(hi));
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XP_JOIN, N, 0, 0, &hit)) return (BQD)hit;
    const BQD r = CALL(bqd_xp_apply_op, BQD_OP_XPROD, lo, unit(limdd_target(hi)));
    cache_put3(CACHE_BQD_XP_JOIN, N, 0, 0, (uint64_t)r);
    return r;
}

/**
 * skip:lem:xcof. With g = ext[N] read at var v and a point (b, y'),
 *
 *     (c Z^s X^t g)(b, y') = c (-1)^{b s_v} (-1)^{s'.y'} g(b ^ t_v, y' ^ t'),
 *
 * which is l' = c (-1)^{b s_v} Z^{s'} X^{t'} applied to the cofactor b ^ t_v
 * of g, s' and t' being s and t without bit v. Where N is below v that
 * cofactor is N itself, and the Z at v is all that made the edge depend on
 * x_v. Otherwise it is the low edge, or the high label times Join(N).
 */
TASK_IMPL_3(BQD, bqd_xp_cofactor, BQD, e, uint32_t, v, int, b)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG N = limdd_target(e);
    label_t l = label_read(limdd_label(e));
    const uint64_t t = l.t & ~above(limdd_level(N));
    const uint32_t top = top_var(N, l.s);
    assert(v <= top && "a labelled edge is read at or above its top");
    if (top > v) return t == l.t ? e : xedge(l, N);    /* e skips v */

    const uint64_t bit = UINT64_C(1) << v;
    const bool sv = (l.s & bit) != 0, tv = (t & bit) != 0;
    l.s &= ~bit;
    l.t = t & ~bit;
    if (sv && b) l.c = neg(l.c);
    if (limdd_level(N) > v) return xedge(l, N);         /* a Z above the node */
    if ((b != 0) == tv) {
        const LIMDD lo = limdd_node_low(N);
        return limdd_edge_is_zero(lo) ? lo : xedge(l, limdd_target(lo));
    }
    const LIMDD hi = limdd_node_high(N);
    if (limdd_edge_is_zero(hi)) return hi;
    const BQD J = CALL(bqd_xp_join, N);
    const label_t h = label_mul(label_mul(l, label_read(limdd_label(hi))),
                                label_read(limdd_label(J)));
    return xedge(h, limdd_target(J));
}

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

/* --- Compose and Canon ------------------------------------------------------- */

/**
 * skip:prop:xcompose, as skip:alg:xcompose. a and b are the canonical edges
 * of f_0 and f_1.
 *
 *   (i)   a = b: f does not depend on x_v, and its edge is a, which skips v.
 *   (ii)  Pauli, b = -a: the sign pattern of def:prep turns f into ext[a],
 *         so the level is skipped with Z_v on a's label (skip:lem:normal (b)).
 *   (iii) a zero: the least point of supp f is (1, t_b), and translating by
 *         it moves f_1 to the low side: the node (v; (id, N_b), 0) under
 *         c_b Z^{s_b} X^{t_b + e_v}, with no sign since b's masks have no bit v.
 *   (iv)  b zero: the node (v; (id, N_a), 0) under a's label.
 *   (v)   otherwise L = l_a^{-1} makes the low cofactor [N_a], 1 at 0, and
 *         the high one u = L . b. With t the least point of supp u the node
 *         is (v; (id, N_a), X^t . k), k = E([N_a] (/) X^t u) and (/) the
 *         xquot: 0 is in supp k, so k's translation is 0 and k(0) = c_k is
 *         the value at the pivot (1, t) of the level. In the Pauli family a
 *         pivot needs an argument in [0, pi), and where c_k has none the
 *         node is multiplied by (-1)^{x_v}: the high scalar is negated and
 *         Z_v goes on the root label.
 *
 * X^t . k is the label product, c_k (-1)^{s_k.t} Z^{s_k} X^t. Its sign is
 * always +1, and so is that of l_a^{-1}: a Z of a canonical label sits at a
 * pivot or at a level the function ignores, and a 1 of a least point at
 * neither, so s_k.t = s_a.t_a = 0 (test_bqd_skip checks s.t = 0 on every
 * edge). The products are taken as they are all the same. The nodes of (iii)
 * to (v) are not redundant.
 */
TASK_IMPL_3(BQD, bqd_xp_compose, uint32_t, v, BQD, a, BQD, b)
{
    if (a == b) return a;
    const uint64_t bit = UINT64_C(1) << v;
    if (limdd_edge_is_zero(a)) {
        label_t lb = label_read(limdd_label(b));
        lb.t |= bit;
        const LIMDD_TARG node = limdd_makenode(v, unit(limdd_target(b)), limdd_zero_edge());
        return limdd_bundle(label_lim(lb), node);
    }
    const LIMDD_TARG na = limdd_target(a);
    if (limdd_edge_is_zero(b))
        return limdd_bundle(limdd_label(a), limdd_makenode(v, unit(na), limdd_zero_edge()));

    const bool pauli = bqd_family() == BQD_FAMILY_PAULI;
    const LIMDD_TARG nb = limdd_target(b);
    label_t la = label_read(limdd_label(a));
    const label_t lb = label_read(limdd_label(b));
    if (pauli && na == nb && la.s == lb.s && la.t == lb.t && lb.c == neg(la.c)) {
        la.s |= bit;
        return limdd_bundle(label_lim(la), na);
    }

    label_t u = label_mul(label_inv(la), lb);
    u.t &= ~above(limdd_level(nb));
    const label_t xt = { EVBDD_ONE, 0, minpoint(nb, u.t) };
    const BQD h = limdd_refs_push(xedge(label_mul(xt, u), nb));
    const BQD k = CALL(bqd_xp_apply_op, BQD_OP_XQUOT, h, unit(na));
    limdd_refs_pop(1);
    label_t lk = label_read(limdd_label(k));
    assert(!limdd_edge_is_zero(k) && lk.t == 0 && "the ratio is nonzero at 0");
    if (pauli && !bqd_arg_in_upper(lk.c)) {
        lk.c = neg(lk.c);
        la.s |= bit;
    }
    const LIMDD_TARG node = limdd_makenode(v, unit(na), xedge(label_mul(xt, lk), limdd_target(k)));
    return limdd_bundle(label_lim(la), node);
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

/**
 * The canonical edge of a labelled edge (skip:prop:xcanon): Compose of the
 * canonical edges of its two cofactors at its top. A label with no Z and no X
 * is a scalar on a stored node, which is canonical as it is. The memo is on
 * the edge without its scalar, which goes back on afterwards, since the
 * representative does not see a scale (skip:lem:xscale).
 */
TASK_IMPL_1(BQD, bqd_xp_canon, BQD, e)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG N = limdd_target(e);
    const label_t l = label_read(limdd_label(e));
    const uint64_t t = l.t & ~above(limdd_level(N));
    if (l.s == 0 && t == 0) return t == l.t ? e : xedge(l, N);

    uint64_t hit;
    if (cache_get3(CACHE_BQD_XP_CANON, 0, key_edge(e), 0, &hit)) return xscale(l.c, (BQD)hit);

    const uint32_t v = top_var(N, l.s);
    limdd_refs_push(e);
    const BQD e0 = limdd_refs_push(CALL(bqd_xp_cofactor, e, v, 0));
    limdd_refs_spawn(SPAWN(bqd_xp_canon, e0));
    const BQD e1 = limdd_refs_push(CALL(bqd_xp_cofactor, e, v, 1));
    const BQD hi = limdd_refs_push(CALL(bqd_xp_canon, e1));
    const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xp_canon)));
    const BQD r = limdd_refs_push(CALL(bqd_xp_compose, v, lo, hi));

    cache_put3(CACHE_BQD_XP_CANON, 0, key_edge(e), 0, memo_of(r, l.c));
    limdd_refs_pop(6);
    return r;
}

/* --- Apply ------------------------------------------------------------------- */

static uint64_t
op_cache_id(int op)
{
    switch (op) {
    case BQD_OP_MUL:   return CACHE_BQD_XP_MUL;
    case BQD_OP_ADD:   return CACHE_BQD_XP_ADD;
    case BQD_OP_XPROD: return CACHE_BQD_XP_XPROD;
    default:           return CACHE_BQD_XP_XQUOT;
    }
}

static EVBDD_WGT
op_scalar(int op, EVBDD_WGT a, EVBDD_WGT b)
{
    switch (op) {
    case BQD_OP_MUL:   return wgt_mul(a, b);
    case BQD_OP_ADD:   return wgt_add(a, b);
    case BQD_OP_XPROD: return (a == EVBDD_ZERO) ? b : wgt_mul(a, b);
    default:           return (b == EVBDD_ZERO) ? a : wgt_div(a, b);
    }
}

/**
 * Apply of skip:alg:xapply. The scalars are taken out as far as the operation
 * allows, and what is left of them stays in the key:
 *
 *     mul    f.g            = c_f c_g . ((f/c_f)(g/c_g))      nothing left
 *     add    f + g          = c_f . (f/c_f + g/c_f)           c_g/c_f left
 *     xprod  xprod(f, g)    = c_g . xprod(f, g/c_g)           c_f left (the copy keeps it)
 *     xquot  xquot(f, g)    = c_f . xquot(f/c_f, g)           c_g left
 *
 * The key is that scalar and the two labelled edges without theirs, each the
 * index of its Pauli word and its node (key_edge), and no level: a labelled
 * edge denotes one function at every level from its top up (skip:lem:ext),
 * and its masks name levels, not positions below a node.
 *
 * The recursion is at the higher of the two tops, where an operand whose top
 * is lower skips the level and is both of its cofactors, and it ends when
 * both are constants: a scalar on the terminal, with no Z. An operand on the
 * terminal with a Z is a character, which is not a constant. Where the scalar
 * family returns an operand, a zero case or a constant one, the operand goes
 * through Canon, since a labelled edge need not be canonical. The cofactors
 * are those of f and g, scalars and all, since every pointwise operation
 * commutes with a cofactor whatever the scalars, so the result comes back
 * scaled, and the memo takes the scale out (memo_of).
 */
TASK_IMPL_3(BQD, bqd_xp_apply_op, int, op, BQD, f, BQD, g)
{
    const bool fz = limdd_edge_is_zero(f), gz = limdd_edge_is_zero(g);
    switch (op) {
    case BQD_OP_MUL:   if (fz || gz) return limdd_zero_edge(); break;
    case BQD_OP_ADD:   if (fz) return CALL(bqd_xp_canon, g);
                       if (gz) return CALL(bqd_xp_canon, f); break;
    case BQD_OP_XPROD: if (fz) return CALL(bqd_xp_canon, g);
                       if (gz) return limdd_zero_edge(); break;
    default:           if (gz) return CALL(bqd_xp_canon, f);
                       if (fz) return limdd_zero_edge(); break;
    }

    const LIMDD_TARG nf = limdd_target(f), ng = limdd_target(g);
    const label_t lf = label_read(limdd_label(f)), lg = label_read(limdd_label(g));
    const bool kf = nf == LIMDD_TERMINAL && lf.s == 0;
    const bool kg = ng == LIMDD_TERMINAL && lg.s == 0;
    if (kf && kg) return constant(op_scalar(op, lf.c, lg.c));
    /* a constant operand is a scale of the other, which is Canon of it */
    if (kf && (op == BQD_OP_MUL || op == BQD_OP_XPROD)) return xscale(lf.c, CALL(bqd_xp_canon, g));
    if (kg && op == BQD_OP_MUL) return xscale(lg.c, CALL(bqd_xp_canon, f));
    if (kg && op == BQD_OP_XQUOT) return xscale(wgt_div(EVBDD_ONE, lg.c), CALL(bqd_xp_canon, f));

    EVBDD_WGT outer, kept;
    switch (op) {
    case BQD_OP_MUL:   outer = wgt_mul(lf.c, lg.c); kept = 0;                   break;
    case BQD_OP_ADD:   outer = lf.c;                kept = wgt_div(lg.c, lf.c); break;
    case BQD_OP_XPROD: outer = lg.c;                kept = lf.c;                break;
    default:           outer = lf.c;                kept = lg.c;                break;
    }
    assert(kept < (UINT64_C(1) << 40) && "a weight index fits under the op id");
    const uint64_t kF = key_edge(f), kG = key_edge(g);

    uint64_t hit;
    if (cache_get3(op_cache_id(op), kept, kF, kG, &hit)) return xscale(outer, (BQD)hit);

    const uint32_t vf = top_var(nf, lf.s), vg = top_var(ng, lg.s);
    const uint32_t v = vf < vg ? vf : vg;
    limdd_refs_push(f);
    limdd_refs_push(g);
    const BQD f0 = limdd_refs_push(CALL(bqd_xp_cofactor, f, v, 0));
    const BQD g0 = limdd_refs_push(CALL(bqd_xp_cofactor, g, v, 0));

    /* the two cofactors of the result are independent, and the high
     * cofactors of the operands, which may need a Join, overlap the low side */
    limdd_refs_spawn(SPAWN(bqd_xp_apply_op, op, f0, g0));
    const BQD f1 = limdd_refs_push(CALL(bqd_xp_cofactor, f, v, 1));
    const BQD g1 = limdd_refs_push(CALL(bqd_xp_cofactor, g, v, 1));
    const BQD hi = limdd_refs_push(CALL(bqd_xp_apply_op, op, f1, g1));
    const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xp_apply_op)));
    const BQD r = limdd_refs_push(CALL(bqd_xp_compose, v, lo, hi));

    cache_put3(op_cache_id(op), kept, kF, kG, memo_of(r, outer));
    limdd_refs_pop(9);
    return r;
}

/* --- restriction, projection and gates --------------------------------------- */

/** The scalar of a labelled edge's label. */
static inline EVBDD_WGT
scalar_of(BQD e)
{
    return limdd_lim_weight(limdd_label(e));
}

/**
 * skip:alg:xrestrict, as bqd_restrict does it in the scalar family, at the
 * top of e instead of at its node. At or below q the restriction is the
 * cofactor at q, which is e itself where e skips q, and a labelled edge
 * either way, so it goes through Canon. Above q the two cofactors are
 * restricted and composed again. Linear, so the memo is on e without its
 * scalar, with q and b in the first word.
 */
TASK_IMPL_3(BQD, bqd_xp_restrict, BQD, e, uint32_t, q, int, b)
{
    if (limdd_edge_is_zero(e)) return e;
    const uint32_t v = bqd_xp_top(e);
    if (v > q) return CALL(bqd_xp_canon, e);
    const EVBDD_WGT c = scalar_of(e);
    const uint64_t key = ((uint64_t)q << 1) | (uint64_t)(b & 1);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XP_RESTRICT, key, key_edge(e), 0, &hit)) return xscale(c, (BQD)hit);

    limdd_refs_push(e);
    BQD r;
    if (v == q) {
        const BQD f = limdd_refs_push(CALL(bqd_xp_cofactor, e, q, b));
        r = CALL(bqd_xp_canon, f);
        limdd_refs_pop(1);
    } else {
        const BQD e0 = limdd_refs_push(CALL(bqd_xp_cofactor, e, v, 0));
        limdd_refs_spawn(SPAWN(bqd_xp_restrict, e0, q, b));
        const BQD e1 = limdd_refs_push(CALL(bqd_xp_cofactor, e, v, 1));
        const BQD hi = limdd_refs_push(CALL(bqd_xp_restrict, e1, q, b));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xp_restrict)));
        r = CALL(bqd_xp_compose, v, lo, hi);
        limdd_refs_pop(4);
    }
    limdd_refs_push(r);
    cache_put3(CACHE_BQD_XP_RESTRICT, key, key_edge(e), 0, memo_of(r, c));
    limdd_refs_pop(2);
    return r;
}

/**
 * skip:alg:xrestrict: e . [x_q = b], as bqd_project does it in the scalar
 * family. At q the chosen cofactor, through Canon, is composed with zero on
 * the other side; where e skips q that cofactor is e itself. Above q the two
 * cofactors are projected and composed again.
 */
TASK_IMPL_3(BQD, bqd_xp_project, BQD, e, uint32_t, q, int, b)
{
    if (limdd_edge_is_zero(e)) return e;
    const uint32_t v = bqd_xp_top(e);
    if (v > q) {
        const BQD g = limdd_refs_push(CALL(bqd_xp_canon, e));
        const BQD r = b ? CALL(bqd_xp_compose, q, limdd_zero_edge(), g)
                        : CALL(bqd_xp_compose, q, g, limdd_zero_edge());
        limdd_refs_pop(1);
        return r;
    }
    const EVBDD_WGT c = scalar_of(e);
    const uint64_t key = ((uint64_t)q << 1) | (uint64_t)(b & 1);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XP_PROJECT, key, key_edge(e), 0, &hit)) return xscale(c, (BQD)hit);

    limdd_refs_push(e);
    BQD r;
    if (v == q) {
        const BQD f = limdd_refs_push(CALL(bqd_xp_cofactor, e, q, b));
        const BQD g = limdd_refs_push(CALL(bqd_xp_canon, f));
        r = b ? CALL(bqd_xp_compose, q, limdd_zero_edge(), g)
              : CALL(bqd_xp_compose, q, g, limdd_zero_edge());
        limdd_refs_pop(2);
    } else {
        const BQD e0 = limdd_refs_push(CALL(bqd_xp_cofactor, e, v, 0));
        limdd_refs_spawn(SPAWN(bqd_xp_project, e0, q, b));
        const BQD e1 = limdd_refs_push(CALL(bqd_xp_cofactor, e, v, 1));
        const BQD hi = limdd_refs_push(CALL(bqd_xp_project, e1, q, b));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xp_project)));
        r = CALL(bqd_xp_compose, v, lo, hi);
        limdd_refs_pop(4);
    }
    limdd_refs_push(r);
    cache_put3(CACHE_BQD_XP_PROJECT, key, key_edge(e), 0, memo_of(r, c));
    limdd_refs_pop(2);
    return r;
}

/**
 * Gate of skip:alg:xapply, which is bqd_cgate_rec of the scalar family at
 * the top of e instead of at its node: U on qubit q where the controls in
 * `cmask` still pending, all above q, are 1. p is the first variable that
 * matters, the highest pending control or q.
 *
 *   e skips p      a control: (e, gate(e)) at p; the target: the row sums of
 *                  U times e at q; e through Canon either way
 *   top above p    the gate on both cofactors
 *   top at a control  the low cofactor through Canon, the gate without that
 *                  control on the high one
 *   top at q       the rows of U applied to the two cofactors, two sums
 *
 * Linear, so the memo is on e without its scalar, the gate's id and cmask
 * with bit q set, as in the scalar family: unambiguous since every control
 * is above q, and free of levels since a control leaves cmask once passed.
 */
TASK_IMPL_4(BQD, bqd_xp_cgate_rec, BQD, e, uint32_t, gid, uint64_t, cmask, uint32_t, q)
{
    if (limdd_edge_is_zero(e)) return e;
    const uint32_t v = bqd_xp_top(e);
    const EVBDD_WGT c = scalar_of(e);
    const uint64_t key = cmask | (UINT64_C(1) << q);
    const uint32_t p = (uint32_t)__builtin_ctzll(key);

    uint64_t hit;
    if (cache_get3(CACHE_BQD_XP_GATE, gid, key_edge(e), key, &hit)) return xscale(c, (BQD)hit);

    limdd_refs_push(e);
    const EVBDD_WGT u00 = gates[gid][0], u01 = gates[gid][1];
    const EVBDD_WGT u10 = gates[gid][2], u11 = gates[gid][3];
    BQD r;
    if (v > p && p < q) {
        limdd_refs_spawn(SPAWN(bqd_xp_canon, e));
        const BQD hi = limdd_refs_push(CALL(bqd_xp_cgate_rec, e, gid,
                                            cmask & ~(UINT64_C(1) << p), q));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xp_canon)));
        r = CALL(bqd_xp_compose, p, lo, hi);
        limdd_refs_pop(2);
    } else if (v > p) {
        const BQD C = limdd_refs_push(CALL(bqd_xp_canon, e));
        const BQD lo = limdd_refs_push(xscale(wgt_add(u00, u01), C));
        r = CALL(bqd_xp_compose, q, lo, xscale(wgt_add(u10, u11), C));
        limdd_refs_pop(2);
    } else if (v < p) {
        const BQD e0 = limdd_refs_push(CALL(bqd_xp_cofactor, e, v, 0));
        limdd_refs_spawn(SPAWN(bqd_xp_cgate_rec, e0, gid, cmask, q));
        const BQD e1 = limdd_refs_push(CALL(bqd_xp_cofactor, e, v, 1));
        const BQD hi = limdd_refs_push(CALL(bqd_xp_cgate_rec, e1, gid, cmask, q));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xp_cgate_rec)));
        r = CALL(bqd_xp_compose, v, lo, hi);
        limdd_refs_pop(4);
    } else if (v < q) {
        const BQD e0 = limdd_refs_push(CALL(bqd_xp_cofactor, e, v, 0));
        limdd_refs_spawn(SPAWN(bqd_xp_canon, e0));
        const BQD e1 = limdd_refs_push(CALL(bqd_xp_cofactor, e, v, 1));
        const BQD hi = limdd_refs_push(CALL(bqd_xp_cgate_rec, e1, gid,
                                            cmask & ~(UINT64_C(1) << v), q));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xp_canon)));
        r = CALL(bqd_xp_compose, v, lo, hi);
        limdd_refs_pop(4);
    } else {
        /* the two new cofactors are independent sums of the old ones */
        const BQD f0 = limdd_refs_push(CALL(bqd_xp_cofactor, e, q, 0));
        const BQD f1 = limdd_refs_push(CALL(bqd_xp_cofactor, e, q, 1));
        const BQD a0 = limdd_refs_push(xscale(u00, f0));
        const BQD a1 = limdd_refs_push(xscale(u01, f1));
        limdd_refs_spawn(SPAWN(bqd_xp_apply_op, BQD_OP_ADD, a0, a1));
        const BQD b0 = limdd_refs_push(xscale(u10, f0));
        const BQD b1 = limdd_refs_push(xscale(u11, f1));
        const BQD hi = limdd_refs_push(CALL(bqd_xp_apply_op, BQD_OP_ADD, b0, b1));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xp_apply_op)));
        r = CALL(bqd_xp_compose, q, lo, hi);
        limdd_refs_pop(8);
    }
    limdd_refs_push(r);
    cache_put3(CACHE_BQD_XP_GATE, gid, key_edge(e), key, memo_of(r, c));
    limdd_refs_pop(2);
    return r;
}

/* --- states and monomials ---------------------------------------------------- */

/**
 * |x> is X^x on |0...0> (skip:alg:xrestrict): the least point of its support
 * is x, and the translate is |0...0>, a node at every level with the chain
 * below it on the low side and zero on the high side, whose sign vector is 0
 * since it has no pivot. So every basis state is one chain under a different
 * label.
 */
BQD
bqd_xp_basis_state(uint64_t x, uint32_t n)
{
    LIMDD_TARG N = LIMDD_TERMINAL;
    label_t l = { EVBDD_ONE, 0, 0 };
    for (uint32_t v = n; v-- > 0; ) {
        N = limdd_makenode(v, unit(N), limdd_zero_edge());
        if ((x >> (n - 1 - v)) & 1) l.t |= UINT64_C(1) << v;   /* qubit v at bit n-1-v */
    }
    return limdd_bundle(label_lim(l), N);
}

/**
 * Mono of skip:alg:xrestrict: phase^{x_A} ignores every variable outside A,
 * and at the first variable a of A its cofactors are the constant one and
 * the monomial of the rest of A, which is below a. So it is Compose of those
 * two, from the last variable of A up. In the Pauli family a phase of -1
 * makes no node at the last variable, where the two cofactors are opposite:
 * (-1)^{x_a} is (Z_a, terminal).
 */
BQD
bqd_xp_monomial(uint64_t A, EVBDD_WGT phase, uint32_t n)
{
    BQD h = limdd_refs_push(constant(phase));
    if (phase != EVBDD_ONE) {
        /* qubit q is bit n-1-q, so the lowest variable is the lowest bit */
        for (uint64_t rest = A; rest != 0; rest &= rest - 1) {
            const uint32_t v = n - 1 - (uint32_t)__builtin_ctzll(rest);
            const BQD next = bqd_xp_compose(v, limdd_one_edge(), h);
            limdd_refs_pop(1);
            h = limdd_refs_push(next);
        }
    }
    limdd_refs_pop(1);
    return h;
}

/* --- the diagonal walk of the Pauli family ----------------------------------- */

/**
 * The canonical edge of the function whose cofactors at var a are [L] and
 * [L] . (h, H), for L and H nodes of a diagram of full support and (h, H) a
 * labelled edge with no X: Compose of skip:alg:xcompose, where the least
 * point is 0 and the ratio is known, so that nothing is divided. A constant
 * ratio of 1 skips the level, one of -1 skips it with Z_a, and a ratio whose
 * value at 0 has no argument in [0, pi) gets the sign repair.
 */
static BQD
pauli_node(uint32_t a, LIMDD_TARG L, label_t h, LIMDD_TARG H)
{
    const label_t za = { EVBDD_ONE, UINT64_C(1) << a, 0 };
    if (H == LIMDD_TERMINAL && h.s == 0) {
        if (h.c == EVBDD_ONE) return unit(L);
        if (h.c == EVBDD_MIN_ONE) return limdd_bundle(label_lim(za), L);
    }
    const bool repair = !bqd_arg_in_upper(h.c);
    if (repair) h.c = neg(h.c);
    const LIMDD_TARG node = limdd_makenode(a, unit(L), limdd_bundle(label_lim(h), H));
    return repair ? limdd_bundle(label_lim(za), node) : unit(node);
}

/**
 * Diag of skip:alg:xdiag: the canonical edge of [N] . phase^{x_A}, N a node
 * of a diagram of full support, A a vector-index mask. On full support the
 * least point of every function is 0, so no label has an X, the pivot of a
 * level is the point with only that level's bit set, and a Z toggles its sign
 * and no other: the canonical edge of Z^s g is Z^s times that of g, as a
 * label product. So the walk is the scalar family's (qsylvan_bqd_ops.c,
 * skip:prop:diag), one path, with a label product where that walk multiplies
 * scalars, and a level of A where the ratio changes goes through pauli_node:
 *
 *   N below a       the node at a over N with ratio phase^{x_rest}, the
 *                   monomial of the rest of A, or phase at the last of A
 *   N above a       the ratio of N is unchanged, so its high edge is kept and
 *                   its node needs no repair; the low child becomes the walk
 *                   on it, whose label, with every Z it carries, goes on top
 *   N at a          the low child is kept, and the ratio is the high label
 *                   times the walk on the high child, or phase times the high
 *                   edge at the last of A
 *
 * Outside A the ratio is N's own, which is neither 1 nor -1 on a canonical N,
 * so only a level of A can be skipped or repaired, and at a level of A before
 * the last the ratio's value at 0 is the high scalar of N, whose argument is
 * in [0, pi) already. So the repair is at the last level of A at most.
 */
static BQD
pauli_diag(LIMDD_TARG N, uint64_t A, EVBDD_WGT phase, uint32_t n, uint32_t *visits)
{
    /* a, the first variable of A, is its highest bit, since qubit q is bit n-1-q */
    const uint32_t hb = 63 - (uint32_t)__builtin_clzll(A);
    const uint64_t rest = A & ~(UINT64_C(1) << hb);
    const uint32_t a = n - 1 - hb;
    const uint32_t var = limdd_level(N);
    if (var > a) {
        if (rest == 0) {
            const label_t h = { phase, 0, 0 };
            return pauli_node(a, N, h, LIMDD_TERMINAL);
        }
        const BQD m = limdd_refs_push(bqd_xp_monomial(rest, phase, n));
        const BQD r = pauli_node(a, N, label_read(limdd_label(m)), limdd_target(m));
        limdd_refs_pop(1);
        return r;
    }

    (*visits)++;
    const LIMDD lo = limdd_node_low(N), hi = limdd_node_high(N);
    if (var < a) {
        const BQD d = limdd_refs_push(pauli_diag(limdd_target(lo), A, phase, n, visits));
        const LIMDD_TARG node = limdd_makenode(var, unit(limdd_target(d)), hi);
        limdd_refs_pop(1);
        return limdd_bundle(limdd_label(d), node);
    }
    label_t h = label_read(limdd_label(hi));
    if (rest == 0) {
        h.c = wgt_mul(phase, h.c);
        return pauli_node(a, limdd_target(lo), h, limdd_target(hi));
    }
    const BQD d = limdd_refs_push(pauli_diag(limdd_target(hi), rest, phase, n, visits));
    const BQD r = pauli_node(a, limdd_target(lo), label_mul(h, label_read(limdd_label(d))),
                             limdd_target(d));
    limdd_refs_pop(1);
    return r;
}

BQD
bqd_xp_diagonal(BQD e, uint64_t A, EVBDD_WGT phase, uint32_t n, uint32_t *visits)
{
    assert(bqd_family() == BQD_FAMILY_PAULI && bqd_has_full_support(e));
    uint32_t count = 0;
    BQD r = e;
    if (A == 0 || phase == EVBDD_ONE) {
        r = xscale(phase, e);
    } else if (!limdd_edge_is_zero(e)) {
        /* e is walked across the collections Mono can meet */
        limdd_refs_push(e);
        const BQD d = limdd_refs_push(pauli_diag(limdd_target(e), A, phase, n, &count));
        const label_t l = label_mul(label_read(limdd_label(e)), label_read(limdd_label(d)));
        r = (l.c == EVBDD_ZERO) ? limdd_zero_edge() : limdd_bundle(label_lim(l), limdd_target(d));
        limdd_refs_pop(2);
    }
    if (visits != NULL) *visits = count;
    return r;
}
