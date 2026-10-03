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
#include "qsylvan_bqd_xp_int.h"
#include "qsylvan_bqd_exp.h"
#include "qsylvan_bqd_gates.h"
#include "qsylvan_gates.h"
#include "qsylvan_limdd_gc.h"

/* The labels, the labelled edges, the least point and the conjugation of a
 * label by a permutation are in qsylvan_bqd_xp_int.h, which the translation
 * and Pauli families share with their rule SM. */

/* --- labelled edges ---------------------------------------------------------- */

uint32_t
bqd_xp_top(BQD e)
{
    if (limdd_edge_is_zero(e)) return limdd_level(LIMDD_TERMINAL);
    return top_var(limdd_target(e), label_read(limdd_label(e)).s);
}

BQD
bqd_xp_scale(BQD e, EVBDD_WGT c)
{
    return xscale(c, e);
}

/* --- cofactors --------------------------------------------------------------- */

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
    bqd_count(LACE_WORKER_ID, BQD_COUNT_COF1);
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

    /* the translation family: (c X^t, N) is canonical when t is the least
     * point of its support, and CanonT moves the node otherwise, on its
     * stored edges where the least point stays (skip:alg:xratio) */
    if (bqd_family() == BQD_FAMILY_X) {
        if (minpoint(N, t) == t) return t == l.t ? e : xedge(l, N);
        return xscale(l.c, CALL(bqd_xp_canon_t, N, t));
    }

    uint64_t hit;
    if (cache_get3(CACHE_BQD_XP_CANON, 0, key_edge(e), 0, &hit)) return xscale(l.c, (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_CANON);

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
    bqd_count(LACE_WORKER_ID, BQD_COUNT_APPLY);

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

/* --- gates without the high cofactor (skip:sec:ratio:xp) -------------------- */

/*
 * In these families the high edge (l_1, N_1) of a node is not its ratio: the
 * high cofactor is l_1 . ([N_0] (.) [N_1]), and l_1 may translate. What
 * carries over from the scalar family (qsylvan_bqd_gates.c) is that a
 * permutation of two qubits commutes with (.) and conjugates a label,
 *
 *     pi . (c Z^s X^t . g) = c Z^{pi^T s} X^{pi t} . pi g,
 *
 * so PermN maps a node through its two stored edges, and ProdHigh makes the
 * node of the result without a product wherever the least point of the high
 * cofactor's support stays where it was; the translation family's Canon does
 * the same for a translation (CanonT), and so X; and the phase
 * multiplications carry the labels along (skip:alg:xphase). Pair recurses on
 * labelled cofactors, since its two operands carry different high labels.
 */

/**
 * Ind of skip:alg:constructors in these families: the indicator of supp [t]
 * has the X parts of t's labels, since the high cofactor's support is the
 * high child's moved by t_1 (the support shadow), and none of their scalars
 * or signs, which the indicator does not see. It is 1 at 0, so its low edge is
 * the identity; a node whose two edges are then equal skips its level. The
 * high edge is normal: a bit of t_1 at a level the indicator of N_1 skips
 * would be a least point with a 1 where the support is closed under flipping
 * it, which it is not. Memoised on t.
 */
TASK_IMPL_1(LIMDD_TARG, bqd_xp_ind, LIMDD_TARG, t)
{
    if (t == LIMDD_TERMINAL) return t;
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XP_IND, t, 0, 0, &hit)) return (LIMDD_TARG)hit;
    const LIMDD lo = limdd_node_low(t), hi = limdd_node_high(t);
    const LIMDD_TARG n0 = node_or_zero(lo), n1 = node_or_zero(hi);
    if (n0 != 0) SPAWN(bqd_xp_ind, n0);
    const LIMDD_TARG i1 = (n1 != 0) ? CALL(bqd_xp_ind, n1) : 0;
    const LIMDD_TARG i0 = (n0 != 0) ? SYNC(bqd_xp_ind) : 0;
    const BQD e0 = i0 ? unit(i0) : limdd_zero_edge();
    const BQD e1 = i1 ? xedge(xlabel(label_read(limdd_label(hi)).t), i1) : limdd_zero_edge();
    const LIMDD_TARG r = (i0 != 0 && e0 == e1) ? i0 : limdd_makenode(limdd_node_var(t), e0, e1);
    cache_put3(CACHE_BQD_XP_IND, t, 0, 0, r);
    return r;
}

/** iota(t): the indicator of supp [t] as a 0/1 exponent, the high child's translated by t_1. */
TASK_IMPL_1(BQD_EXP, bqd_xp_iota, LIMDD_TARG, t)
{
    if (t == 0) return mtbdd_int64(0);
    if (t == LIMDD_TERMINAL) return mtbdd_int64(1);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XP_IOTA, t, 0, 0, &hit)) return (BQD_EXP)hit;
    const LIMDD hi = limdd_node_high(t);
    mtbdd_refs_spawn(SPAWN(bqd_xp_iota, node_or_zero(limdd_node_low(t))));
    BQD_EXP i1 = mtbdd_refs_push(CALL(bqd_xp_iota, node_or_zero(hi)));
    if (!limdd_edge_is_zero(hi)) {
        const uint64_t t1 = label_read(limdd_label(hi)).t;
        const BQD_EXP moved = CALL(bqd_exp_translate, i1, t1);
        mtbdd_refs_pop(1);
        i1 = mtbdd_refs_push(moved);
    }
    const BQD_EXP i0 = mtbdd_refs_push(mtbdd_refs_sync(SYNC(bqd_xp_iota)));
    const BQD_EXP r = bqd_exp_node(limdd_node_var(t), i0, i1);
    mtbdd_refs_pop(2);
    cache_put3(CACHE_BQD_XP_IOTA, t, 0, 0, r);
    return r;
}

/* --- phase multiplications ----------------------------------------------------- */

/**
 * PhaseMulX of skip:alg:xphase, the translation family: the canonical edge of
 * [N] . beta^eps, r the order of beta. Where N skips the level it is read
 * there as its virtual node (N, (id, Ind N)). The low call gives (alpha, U)
 * with alpha = beta^{eps_0(0)}, since 0 is the least point of every stored
 * node; the high cofactor keeps its support, so its least point is still t_1
 * and the quotient of case (v) of skip:prop:xcompose is c_1 [N_1] times
 * beta to X^{t_1} eps_1 - eps_0 on supp [N_0] and X^{t_1} eps_1 - v off it,
 * one call with that exponent. No Cof1, no Apply and no Canon. Memoised on N,
 * eps and beta; the two calls depend on the input alone, and overlap.
 */
TASK_IMPL_4(BQD, bqd_xp_phase_mul_x, LIMDD_TARG, N, BQD_EXP, eps, EVBDD_WGT, beta, uint32_t, r)
{
    if (bqd_exp_is_const(eps)) return xscale(bqd_exp_power(beta, bqd_exp_value(eps)), unit(N));
    if (N == LIMDD_TERMINAL) return CALL(bqd_exp_state, eps, beta, r);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XP_PHASEMUL, N, eps, beta, &hit)) return (BQD)hit;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PHASEMUL);

    const uint32_t lN = limdd_level(N), le = bqd_exp_var(eps);
    const uint32_t l = lN < le ? lN : le;
    const LIMDD_TARG N0 = lN > l ? N : limdd_target(limdd_node_low(N));
    const BQD hi = lN > l ? unit(CALL(bqd_xp_ind, N)) : limdd_node_high(N);
    const BQD_EXP e0 = bqd_exp_cof(eps, l, 0), e1 = bqd_exp_cof(eps, l, 1);
    BQD res;
    if (limdd_edge_is_zero(hi)) {
        const BQD lo = CALL(bqd_xp_phase_mul_x, N0, e0, beta, r);
        res = limdd_bundle(limdd_label(lo), limdd_makenode(l, unit(limdd_target(lo)), hi));
    } else {
        const label_t l1 = label_read(limdd_label(hi));
        const int64_t v = bqd_exp_eval(e0, 0);
        const BQD_EXP io = mtbdd_refs_push(CALL(bqd_xp_iota, N0));
        const BQD_EXP sel = mtbdd_refs_push(CALL(bqd_exp_sel, io, e0, v));
        const BQD_EXP moved = mtbdd_refs_push(CALL(bqd_exp_translate, e1, l1.t));
        const BQD_EXP e2 = mtbdd_refs_push(CALL(bqd_exp_sub, moved, sel, r));
        limdd_refs_spawn(SPAWN(bqd_xp_phase_mul_x, N0, e0, beta, r));
        const BQD k = limdd_refs_push(xscale(l1.c, CALL(bqd_xp_phase_mul_x, limdd_target(hi),
                                                        e2, beta, r)));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xp_phase_mul_x)));
        const LIMDD_TARG U = limdd_target(lo);
        if (l1.t == 0 && limdd_label(k) == LIMDD_LIM_IDENTITY
            && limdd_target(k) == CALL(bqd_xp_ind, U)) {
            res = lo;                                   /* f_1 = f_0: the level is skipped */
        } else {
            const BQD h = xedge(label_mul(xlabel(l1.t), label_read(limdd_label(k))), limdd_target(k));
            res = limdd_bundle(limdd_label(lo), limdd_makenode(l, unit(U), h));
        }
        limdd_refs_pop(2);
        mtbdd_refs_pop(4);
    }
    cache_put3(CACHE_BQD_XP_PHASEMUL, N, eps, beta, (uint64_t)res);
    return res;
}

/**
 * PhaseMulP of skip:alg:xphase, the Pauli family: the canonical edge of
 * [N] . w^eps, w = e^{i pi/4} and eps modulo 8. As PhaseMulX, with the sign
 * patterns carried along: the low result is (c_0 Z^{s_0}, U) with
 * [N_0] = [U] w^delta, delta = v - eps_0 + P_{s_0}, and moving X^{t_1} past
 * Z^w, w = s_0 xor s_1, costs (-1)^{w . t_1}, so the quotient is c_1 (-1)^{w.t_1}
 * [N_1] w^{eps'} with eps' = X^{t_1} eps_1 - v + P_w, plus delta on
 * supp [N_0]. Equal cofactors skip the level, opposite ones skip it with Z_l,
 * and a pivot without an argument in [0, pi) gets the sign repair. The high
 * call needs the low result's sign pattern, so the two run one after the
 * other. The terminal is its own virtual node, (T, (id, T)). Memoised on N and
 * eps.
 */
TASK_IMPL_2(BQD, bqd_xp_phase_mul_p, LIMDD_TARG, N, BQD_EXP, eps)
{
    if (bqd_exp_is_const(eps)) {
        const EVBDD_WGT w = gates[GATEID_T][3];
        return xscale(bqd_exp_power(w, bqd_exp_value(eps)), unit(N));
    }
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XP_PHASEMULP, N, eps, 0, &hit)) return (BQD)hit;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PHASEMUL);

    const uint32_t lN = limdd_level(N), le = bqd_exp_var(eps);
    const uint32_t l = lN < le ? lN : le;
    const uint64_t bit = UINT64_C(1) << l;
    const LIMDD_TARG N0 = lN > l ? N : limdd_target(limdd_node_low(N));
    const BQD hi = lN > l ? unit(CALL(bqd_xp_ind, N)) : limdd_node_high(N);
    const BQD_EXP e0 = bqd_exp_cof(eps, l, 0), e1 = bqd_exp_cof(eps, l, 1);
    const BQD lo = limdd_refs_push(CALL(bqd_xp_phase_mul_p, N0, e0));
    const label_t l0 = label_read(limdd_label(lo));
    const LIMDD_TARG U = limdd_target(lo);
    BQD res;
    if (limdd_edge_is_zero(hi)) {
        res = limdd_bundle(limdd_label(lo), limdd_makenode(l, unit(U), hi));
    } else {
        const label_t l1 = label_read(limdd_label(hi));
        const int64_t v = bqd_exp_eval(e0, 0);
        const uint64_t w = l0.s ^ l1.s;
        /* delta = v - eps_0 + P_{s_0}, and eps' = X^{t_1} eps_1 - v + Sel(iota N_0, delta, 0) + P_w */
        const BQD_EXP cv = mtbdd_refs_push(bqd_exp_const(v, 8));
        const BQD_EXP d1 = mtbdd_refs_push(CALL(bqd_exp_sub, cv, e0, 8));
        const BQD_EXP ps0 = mtbdd_refs_push(bqd_exp_parity4(l0.s));
        const BQD_EXP delta = mtbdd_refs_push(CALL(bqd_exp_add, d1, ps0, 8));
        const BQD_EXP io = mtbdd_refs_push(CALL(bqd_xp_iota, N0));
        const BQD_EXP sel = mtbdd_refs_push(CALL(bqd_exp_sel, io, delta, 0));
        const BQD_EXP moved = mtbdd_refs_push(CALL(bqd_exp_translate, e1, l1.t));
        const BQD_EXP a1 = mtbdd_refs_push(CALL(bqd_exp_sub, moved, cv, 8));
        const BQD_EXP a2 = mtbdd_refs_push(CALL(bqd_exp_add, a1, sel, 8));
        const BQD_EXP pw = mtbdd_refs_push(bqd_exp_parity4(w));
        const BQD_EXP e2 = mtbdd_refs_push(CALL(bqd_exp_add, a2, pw, 8));
        const EVBDD_WGT gam = parity(w & l1.t) ? neg(l1.c) : l1.c;
        const BQD k = limdd_refs_push(xscale(gam, CALL(bqd_xp_phase_mul_p, limdd_target(hi), e2)));
        mtbdd_refs_pop(11);
        label_t lk = label_read(limdd_label(k));
        const LIMDD_TARG ind = (l1.t == 0 && lk.s == 0 && lk.t == 0) ? CALL(bqd_xp_ind, U) : 0;
        if (ind != 0 && limdd_target(k) == ind && lk.c == EVBDD_ONE) {
            res = lo;                                   /* f_1 = f_0 */
        } else if (ind != 0 && limdd_target(k) == ind && lk.c == EVBDD_MIN_ONE) {
            label_t z = l0;                             /* f_1 = -f_0: skipped with Z_l */
            z.s ^= bit;
            res = limdd_bundle(label_lim(z), U);
        } else {
            label_t root = l0;
            if (!bqd_arg_in_upper(lk.c)) {              /* the sign repair */
                lk.c = neg(lk.c);
                root.s ^= bit;
            }
            const BQD h = xedge(label_mul(xlabel(l1.t), lk), limdd_target(k));
            res = limdd_bundle(label_lim(root), limdd_makenode(l, unit(U), h));
        }
        limdd_refs_pop(1);
    }
    limdd_refs_pop(1);
    cache_put3(CACHE_BQD_XP_PHASEMULP, N, eps, 0, (uint64_t)res);
    return res;
}

/**
 * MulOffX: the canonical edge of [k] . c^{[x not in supp [U]]} for a
 * canonical k = (c' X^t, V) of the translation family, which is
 * c' X^t . PhaseMulX(V, X^t (1 - iota U), c): the multiplication keeps the
 * support, so t stays its least point and the label is canonical.
 */
TASK_IMPL_3(BQD, bqd_xp_mul_off, BQD, k, EVBDD_WGT, c, LIMDD_TARG, U)
{
    if (limdd_edge_is_zero(k) || c == EVBDD_ONE || CALL(bqd_xp_ind, U) == LIMDD_TERMINAL) return k;
    const uint32_t r = bqd_exp_order(c);
    const label_t lk = label_read(limdd_label(k));
    const BQD_EXP one = mtbdd_refs_push(bqd_exp_const(1, r));
    const BQD_EXP io = mtbdd_refs_push(CALL(bqd_xp_iota, U));
    const BQD_EXP off = mtbdd_refs_push(CALL(bqd_exp_sub, one, io, r));
    const BQD_EXP eps = mtbdd_refs_push(CALL(bqd_exp_translate, off, lk.t));
    const BQD p = CALL(bqd_xp_phase_mul_x, limdd_target(k), eps, c, r);
    mtbdd_refs_pop(4);
    return xedge(label_mul(lk, label_read(limdd_label(p))), limdd_target(p));
}

/* --- Canon of a translation ---------------------------------------------------- */

/**
 * CanonT of skip:alg:xratio, the translation family: the canonical edge of
 * X^tau [N]. Where tau has a bit at N's level the two cofactors trade places,
 * and that is one step of the generic Canon. Otherwise the low cofactor is
 * X^tau [N_0], whose canonical edge (c_0 X^{t_0}, U) CanonT gives, and the
 * high one c_1 X^{t_1 xor t_0} (c_0 [U] (.) X^delta [N_1]), delta = tau xor
 * t_0; where the least point of its support is still t_1 the node's high edge
 * is X^{t_1} . c_1 MulOffX(CanonT(N_1, delta), 1/c_0, U), and otherwise it is
 * Compose of the two. No skip: a translation is a bijection, and this family
 * has no sign repair. Memoised on N and tau; the high call needs t_0.
 */
TASK_IMPL_2(BQD, bqd_xp_canon_t, LIMDD_TARG, N, uint64_t, tau)
{
    tau &= ~above(limdd_level(N));
    if (tau == 0) return unit(N);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XP_CANONT, N, tau, 0, &hit)) return (BQD)hit;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    const uint32_t l = limdd_node_var(N);
    BQD r;
    if ((tau >> l) & 1) {
        const BQD e = limdd_refs_push(xedge(xlabel(tau), N));
        const BQD c0 = limdd_refs_push(CALL(bqd_xp_cofactor, e, l, 0));
        limdd_refs_spawn(SPAWN(bqd_xp_canon, c0));
        const BQD c1 = limdd_refs_push(CALL(bqd_xp_cofactor, e, l, 1));
        const BQD hi = limdd_refs_push(CALL(bqd_xp_canon, c1));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xp_canon)));
        r = CALL(bqd_xp_compose, l, lo, hi);
        limdd_refs_pop(5);
    } else {
        const BQD lo = limdd_refs_push(CALL(bqd_xp_canon_t, limdd_target(limdd_node_low(N)), tau));
        const LIMDD_TARG U = limdd_target(lo);
        const label_t l0 = label_read(limdd_label(lo));
        const LIMDD high = limdd_node_high(N);
        if (limdd_edge_is_zero(high)) {
            r = limdd_bundle(limdd_label(lo), limdd_makenode(l, unit(U), high));
        } else {
            const label_t l1 = label_read(limdd_label(high));
            const LIMDD_TARG N1 = limdd_target(high);
            const uint64_t delta = tau ^ l0.t;
            if (minpoint(N1, delta ^ l1.t) == l1.t) {
                const BQD m = limdd_refs_push(CALL(bqd_xp_canon_t, N1, delta));
                const BQD k = xscale(l1.c, CALL(bqd_xp_mul_off, m, wgt_div(EVBDD_ONE, l0.c), U));
                const BQD h = xedge(label_mul(xlabel(l1.t), label_read(limdd_label(k))), limdd_target(k));
                r = limdd_bundle(limdd_label(lo), limdd_makenode(l, unit(U), h));
                limdd_refs_pop(1);
            } else {
                const BQD f1 = limdd_refs_push(CALL(bqd_xp_cofactor, unit(N), l, 1));
                const BQD moved = limdd_refs_push(xedge(label_mul(xlabel(tau), label_read(limdd_label(f1))),
                                                        limdd_target(f1)));
                const BQD hi = limdd_refs_push(CALL(bqd_xp_canon, moved));
                r = CALL(bqd_xp_compose, l, lo, hi);
                limdd_refs_pop(3);
            }
        }
        limdd_refs_pop(1);
    }
    cache_put3(CACHE_BQD_XP_CANONT, N, tau, 0, (uint64_t)r);
    return r;
}

/* --- permutations of two qubits ------------------------------------------------ */

/**
 * Pair of skip:alg:perm on labelled edges: [x_b = 0] A|_{x_b = s} +
 * [x_b = 1] B|_{x_b = t}, recursing on the labelled cofactors at the higher
 * of the two tops (the Shannon level), with Canon on the edges it composes.
 * Linear in the pair, so the scalar g of the first nonzero operand comes out,
 * and the memo is on the two edges without it, B's scalar over g, s, t and b:
 * the scalar and s and b in the first word, below 2^40 since a weight index
 * is below 2^33 and b below 64, and t in the id.
 */
TASK_IMPL_5(BQD, bqd_xp_pair, BQD, A, int, s, BQD, B, int, t, uint32_t, b)
{
    const bool az = limdd_edge_is_zero(A), bz = limdd_edge_is_zero(B);
    if (az && bz) return A;
    const EVBDD_WGT g = az ? scalar_of(B) : scalar_of(A);
    const EVBDD_WGT gi = wgt_div(EVBDD_ONE, g);
    const BQD A1 = az ? A : xscale(gi, A);
    const BQD B1 = bz ? B : xscale(gi, B);
    const uint64_t kappa = bz ? 0 : (uint64_t)scalar_of(B1);
    assert(kappa < (UINT64_C(1) << 33) && "a weight index is below 2^33");
    const uint64_t k0 = kappa | ((uint64_t)(s & 1) << 33) | ((uint64_t)b << 34);
    const uint64_t kA = az ? 0 : key_edge(A1), kB = bz ? 0 : key_edge(B1);
    const uint64_t id = t ? CACHE_BQD_XP_PAIR1 : CACHE_BQD_XP_PAIR0;
    uint64_t hit;
    if (cache_get3(id, k0, kA, kB, &hit)) return xscale(g, (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    limdd_refs_push(A1);
    limdd_refs_push(B1);
    const uint32_t ta = bqd_xp_top(A1), tb = bqd_xp_top(B1);
    const uint32_t l = ta < tb ? ta : tb;
    BQD r;
    if (l > b) {
        const BQD ca = limdd_refs_push(CALL(bqd_xp_canon, A1));
        const BQD cb = limdd_refs_push(CALL(bqd_xp_canon, B1));
        r = CALL(bqd_xp_compose, b, ca, cb);
        limdd_refs_pop(2);
    } else if (l == b) {
        const BQD fa = limdd_refs_push(CALL(bqd_xp_cofactor, A1, b, s));
        const BQD fb = limdd_refs_push(CALL(bqd_xp_cofactor, B1, b, t));
        limdd_refs_spawn(SPAWN(bqd_xp_canon, fa));
        const BQD cb = limdd_refs_push(CALL(bqd_xp_canon, fb));
        const BQD ca = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xp_canon)));
        r = CALL(bqd_xp_compose, b, ca, cb);
        limdd_refs_pop(4);
    } else {
        const BQD a0 = limdd_refs_push(CALL(bqd_xp_cofactor, A1, l, 0));
        const BQD b0 = limdd_refs_push(CALL(bqd_xp_cofactor, B1, l, 0));
        limdd_refs_spawn(SPAWN(bqd_xp_pair, a0, s, b0, t, b));
        const BQD a1 = limdd_refs_push(CALL(bqd_xp_cofactor, A1, l, 1));
        const BQD b1 = limdd_refs_push(CALL(bqd_xp_cofactor, B1, l, 1));
        const BQD hi = limdd_refs_push(CALL(bqd_xp_pair, a1, s, b1, t, b));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xp_pair)));
        r = CALL(bqd_xp_compose, l, lo, hi);
        limdd_refs_pop(6);
    }
    cache_put3(id, k0, kA, kB, (uint64_t)r);
    limdd_refs_pop(2);
    return xscale(g, r);
}

/**
 * ProdHigh of skip:alg:xratio: Compose(j, lo, Canon(lam . ([lo] (.) [W]))),
 * for lo and W canonical edges of functions with 0 in their support and the
 * value 1 there, whose cofactors are not equal or opposite. Where lo carries
 * the identity and W a sign pattern at most, lam . ([U] (.) Z^{s_W} [V]) is
 * c Z^s X^t . ([U] (.) [V]) with c Z^s X^t = lam Z^{s_W}, and where t is the
 * least point of the support of that, the ratio of case (v) of
 * skip:prop:xcompose is X^t c Z^s X^t [V]: no product and no quotient, only a
 * Canon of a label without X part, which the translation family does not
 * even need, and the sign repair of the Pauli family. Otherwise the product
 * is built.
 */
TASK_6(BQD, bqd_xp_prod_high, uint32_t, j, BQD, lo, BQD, W, EVBDD_WGT, lc, uint64_t, ls,
       uint64_t, lt)
{
    const label_t lam = { lc, ls, lt };
    const label_t lw = label_read(limdd_label(W));
    if (limdd_label(lo) == LIMDD_LIM_IDENTITY && lw.c == EVBDD_ONE && lw.t == 0) {
        const label_t zw = { EVBDD_ONE, lw.s, 0 };
        const label_t cst = label_mul(lam, zw);
        const LIMDD_TARG V = limdd_target(W);
        if (minpoint(V, cst.t) == cst.t) {
            const label_t xt = xlabel(cst.t);
            const BQD k = CALL(bqd_xp_canon, xedge(label_mul(xt, cst), V));
            label_t lk = label_read(limdd_label(k));
            label_t root = { EVBDD_ONE, 0, 0 };
            if (bqd_family() == BQD_FAMILY_PAULI && !bqd_arg_in_upper(lk.c)) {
                lk.c = neg(lk.c);
                root.s = UINT64_C(1) << j;
            }
            const BQD h = xedge(label_mul(xt, lk), limdd_target(k));
            return limdd_bundle(label_lim(root), limdd_makenode(j, unit(limdd_target(lo)), h));
        }
    }
    const BQD p = limdd_refs_push(CALL(bqd_xp_apply_op, BQD_OP_XPROD, lo, W));
    const BQD m = limdd_refs_push(relabel(lam, p));
    const BQD hi = limdd_refs_push(CALL(bqd_xp_canon, m));
    const BQD r = CALL(bqd_xp_compose, j, lo, hi);
    limdd_refs_pop(3);
    return r;
}

/**
 * PermN of skip:alg:xratio: the canonical edge of pi [N], memoised on N, the
 * kind and the two qubits. Below qa, [N] does not depend on x_qa; at qa the
 * cofactors are the labelled ones, and the new ones Pairs at qb, or the low
 * one and X_qb of the high one; above qa the node's two children are mapped,
 * since pi [N] has the cofactors pi [N_0] and l_1^pi . (pi [N_0] (.) pi [N_1]),
 * and ProdHigh makes the node.
 */
TASK_4(BQD, bqd_xp_perm_node, LIMDD_TARG, N, uint32_t, kind, uint32_t, qa, uint32_t, qb)
{
    const uint32_t var = limdd_level(N);
    if (var > qa && kind == BQD_PERM_CX_UP) return unit(N);
    const uint64_t key = (uint64_t)kind | ((uint64_t)qa << 8) | ((uint64_t)qb << 16);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XP_PERM, N, key, 0, &hit)) return (BQD)hit;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    const BQD E = unit(N);
    BQD r;
    if (var > qa) {
        BQD g0, g1;
        if (kind == BQD_PERM_CX_DOWN) {
            g0 = limdd_refs_push(E);
            g1 = limdd_refs_push(CALL(bqd_xp_canon, xedge(xlabel(UINT64_C(1) << qb), N)));
        } else {
            limdd_refs_spawn(SPAWN(bqd_xp_restrict, E, qb, 1));
            g0 = limdd_refs_push(CALL(bqd_xp_restrict, E, qb, 0));
            g1 = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xp_restrict)));
        }
        r = CALL(bqd_xp_compose, qa, g0, g1);
        limdd_refs_pop(2);
    } else if (var == qa) {
        const BQD f0 = limdd_refs_push(CALL(bqd_xp_cofactor, E, qa, 0));
        const BQD f1 = limdd_refs_push(CALL(bqd_xp_cofactor, E, qa, 1));
        BQD g0, g1;
        if (kind == BQD_PERM_CX_DOWN) {
            limdd_refs_spawn(SPAWN(bqd_xp_canon, f0));
            g1 = limdd_refs_push(CALL(bqd_xp_canon, relabel(xlabel(UINT64_C(1) << qb), f1)));
            g0 = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xp_canon)));
        } else {
            const bool sw = kind == BQD_PERM_SWAP;
            limdd_refs_spawn(SPAWN(bqd_xp_pair, sw ? f0 : f1, sw ? 1 : 0, sw ? f1 : f0, 1, qb));
            g0 = limdd_refs_push(CALL(bqd_xp_pair, f0, 0, f1, sw ? 0 : 1, qb));
            g1 = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xp_pair)));
        }
        r = CALL(bqd_xp_compose, qa, g0, g1);
        limdd_refs_pop(4);
    } else {
        const LIMDD high = limdd_node_high(N);
        const LIMDD_TARG N0 = limdd_target(limdd_node_low(N));
        if (limdd_edge_is_zero(high)) {
            const BQD lo = limdd_refs_push(CALL(bqd_xp_perm_node, N0, kind, qa, qb));
            r = CALL(bqd_xp_compose, var, lo, high);
            limdd_refs_pop(1);
        } else {
            limdd_refs_spawn(SPAWN(bqd_xp_perm_node, N0, kind, qa, qb));
            const BQD W = limdd_refs_push(CALL(bqd_xp_perm_node, limdd_target(high), kind, qa, qb));
            const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xp_perm_node)));
            const label_t lam = perm_label(kind, qa, qb, label_read(limdd_label(high)));
            r = CALL(bqd_xp_prod_high, var, lo, W, lam.c, lam.s, lam.t);
            limdd_refs_pop(2);
        }
    }
    cache_put3(CACHE_BQD_XP_PERM, N, key, 0, (uint64_t)r);
    return r;
}

TASK_IMPL_4(BQD, bqd_xp_perm, BQD, e, uint32_t, kind, uint32_t, qa, uint32_t, qb)
{
    if (limdd_edge_is_zero(e)) return e;
    const label_t le = perm_label(kind, qa, qb, label_read(limdd_label(e)));
    const BQD p = limdd_refs_push(CALL(bqd_xp_perm_node, limdd_target(e), kind, qa, qb));
    const BQD m = limdd_refs_push(relabel(le, p));
    const BQD r = CALL(bqd_xp_canon, m);
    limdd_refs_pop(2);
    return r;
}

/** X on qubit q: the label product X_q . e, made canonical by Canon, CanonT in the translation family. */
TASK_IMPL_2(BQD, bqd_xp_x, BQD, e, uint32_t, q)
{
    if (limdd_edge_is_zero(e)) return e;
    const BQD m = limdd_refs_push(relabel(xlabel(UINT64_C(1) << q), e));
    const BQD r = CALL(bqd_xp_canon, m);
    limdd_refs_pop(1);
    return r;
}

int
bqd_xp_w8_log(EVBDD_WGT beta)
{
    const EVBDD_WGT w = gates[GATEID_T][3];
    EVBDD_WGT p = EVBDD_ONE;
    for (int k = 0; k < 8; k++, p = wgt_mul(p, w)) if (p == beta) return k;
    return -1;
}

/**
 * e . beta^eps on any support in these families: c Z^s X^t . PhaseMulX(N,
 * X^t eps, beta) in the translation family, and c Z^s X^t . PhaseMulP(N,
 * X^t (m eps)) in the Pauli family for beta = w^m; the product label is
 * canonical, since the multiplication keeps the support, and with it the
 * least point and the pivots. eps is protected by the caller.
 */
TASK_IMPL_3(BQD, bqd_xp_phase_mul, BQD, e, BQD_EXP, eps, EVBDD_WGT, beta)
{
    if (limdd_edge_is_zero(e) || beta == EVBDD_ONE) return e;
    const label_t l = label_read(limdd_label(e));
    const LIMDD_TARG N = limdd_target(e);
    BQD p;
    if (bqd_family() == BQD_FAMILY_X) {
        const uint32_t r = bqd_exp_order(beta);
        const BQD_EXP moved = mtbdd_refs_push(CALL(bqd_exp_translate, eps, l.t));
        p = CALL(bqd_xp_phase_mul_x, N, moved, beta, r);
        mtbdd_refs_pop(1);
    } else {
        const int m = bqd_xp_w8_log(beta);
        if (m < 0) {
            fprintf(stderr, "sylvan: the Pauli-BQD multiplies by a power of w_8 only\n");
            exit(1);
        }
        /* m eps, the monomial's exponent in base w */
        const BQD_EXP zero = mtbdd_refs_push(bqd_exp_const(0, 8));
        BQD_EXP me = mtbdd_refs_push(zero);
        for (int k = 0; k < m; k++) {
            const BQD_EXP next = CALL(bqd_exp_add, me, eps, 8);
            mtbdd_refs_pop(1);
            me = mtbdd_refs_push(next);
        }
        const BQD_EXP moved = mtbdd_refs_push(CALL(bqd_exp_translate, me, l.t));
        p = CALL(bqd_xp_phase_mul_p, N, moved);
        mtbdd_refs_pop(3);
    }
    return xedge(label_mul(l, label_read(limdd_label(p))), limdd_target(p));
}
