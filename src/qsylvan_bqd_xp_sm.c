/*
 * Copyright 2026 Q-Sylvan contributors
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

#include "qsylvan_bqd_xp_sm.h"
#include "qsylvan_bqd_xp_int.h"
#include "qsylvan_bqd_xp.h"
#include "qsylvan_bqd_exp.h"
#include "qsylvan_bqd_sm.h"
#include "qsylvan_gates.h"
#include "qsylvan_limdd_gc.h"

/* --- masks and the sign repair ------------------------------------------------ */

static inline uint64_t
bit(uint32_t v)
{
    return UINT64_C(1) << v;
}

/** The qubits below variable v, which are the bits above bit v. */
static inline uint64_t
lower(uint32_t v)
{
    return ~((UINT64_C(2) << v) - 1);
}

/**
 * Whether a pivot value needs the Pauli family's sign repair (def:prep): an
 * argument outside [0, pi). Never in the translation family, which has none.
 */
static inline bool
repair(EVBDD_WGT c)
{
    return bqd_family() == BQD_FAMILY_PAULI && !bqd_arg_in_upper(c);
}

/**
 * The labelled edge e with the scalar c in place of its own, by its label and
 * no arithmetic: the operations recurse on their operands with the scalar
 * that leaves the memo key taken out this way, as the scalar family's do, and
 * put it back on the result once. Carried down the recursion instead, it
 * would multiply into every label below and be divided out of every memo
 * entry: on exact weights that division, a gcd of two large numbers, took
 * three quarters of the time of a Hadamard on clifford_T_circuit_20_700.
 */
static inline BQD
with_scalar(BQD e, EVBDD_WGT c)
{
    if (limdd_edge_is_zero(e)) return e;
    label_t l = label_read(limdd_label(e));
    l.c = c;
    return limdd_bundle(label_lim(l), limdd_target(e));
}

/* --- the support shadow --------------------------------------------------------- */

/*
 * The support of a node is {0} x supp [A] and {1} x (supp [R] ^ t_1) for either
 * tag (skip:lem:smxpshadow), so the walks below read the targets of the two
 * stored edges and the X part of the high label, and nothing else: no scalar,
 * no sign, no product. A node below the level a walk is at is both of its
 * cofactors there.
 */

/**
 * Ind: the node of the indicator of supp [N]. Its cofactors are the indicator
 * of A and that of R moved by t_1, and its tag N's, since the supports are
 * N's: a Q node's moved high support is inside its low one, and an S node's
 * is not. The indicator is 1 at 0, its own representative, so the high edge
 * is X^{t_1} on the indicator of R, and a Q indicator skips its level where
 * the two cofactors are one node with no translation between them. The
 * terminal on full support, from the flag. Memoised on N; the two calls
 * overlap.
 */
TASK_IMPL_1(LIMDD_TARG, bqd_xpsm_ind, LIMDD_TARG, N)
{
    if (N == 0) return 0;
    if (bqd_sm_full(N)) return LIMDD_TERMINAL;
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XPSM_IND, N, 0, 0, &hit)) return (LIMDD_TARG)hit;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_CANON);
    const BQD lo = limdd_node_low(N), hi = limdd_node_high(N);
    assert(!limdd_edge_is_zero(lo) && "0 is in the support of every stored node");
    SPAWN(bqd_xpsm_ind, limdd_target(lo));
    const LIMDD_TARG B = CALL(bqd_xpsm_ind, node_or_zero(hi));
    const LIMDD_TARG A = SYNC(bqd_xpsm_ind);
    const uint32_t var = limdd_node_var(N);
    const bool s = bqd_sm_is_s(N);
    LIMDD_TARG r;
    if (B == 0) {
        r = bqd_sm_mk(var, false, unit(A), limdd_zero_edge());
    } else {
        const uint64_t t1 = label_read(limdd_label(hi)).t;
        r = (!s && t1 == 0 && A == B) ? A : bqd_sm_mk(var, s, unit(A), xedge(xlabel(t1), B));
    }
    cache_put3(CACHE_BQD_XPSM_IND, N, 0, 0, r);
    return r;
}

/**
 * The support of the cofactor x_v = c of X^tau [N], N a node or the terminal
 * read at v: a node, returned, and the translation that moves it, in *tr. A
 * node below v is both of its cofactors; at v the low one is A, and the high
 * one R moved by t_1, or nothing where the high edge is zero.
 */
static inline LIMDD_TARG
shadow(LIMDD_TARG N, uint64_t tau, uint32_t v, int c, uint64_t *tr)
{
    *tr = tau & lower(v);
    if (limdd_level(N) > v) return N;
    if (c == 0) return limdd_target(limdd_node_low(N));
    const BQD hi = limdd_node_high(N);
    if (limdd_edge_is_zero(hi)) { *tr = 0; return 0; }
    *tr ^= label_read(limdd_label(hi)).t;
    return limdd_target(hi);
}

/**
 * Sub: whether supp [X] ^ tau lies inside supp [Y], at the higher of the two
 * nodes, on both halves: the half x_v = b of the left side is the cofactor
 * b ^ tau_v of X moved by the rest of tau. Only the bits of tau at and below
 * the lower node matter, since that node's support is closed under flipping
 * any bit above it, so the others are dropped before the memo, which is on
 * (X, Y, tau). Free where X is zero or Y full, false where Y is zero or X is
 * full and Y is not. The halves are independent, and one is spawned.
 */
TASK_IMPL_3(int, bqd_xpsm_sub, LIMDD_TARG, X, LIMDD_TARG, Y, uint64_t, tau)
{
    if (X == 0 || bqd_sm_full(Y)) return 1;
    if (Y == 0 || bqd_sm_full(X)) return 0;
    const uint32_t lx = limdd_level(X), ly = limdd_level(Y);
    tau &= ~above(lx > ly ? lx : ly);
    if (X == Y && tau == 0) return 1;
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XPSM_SUBSET, X, Y, tau, &hit)) return (int)hit - 1;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_CANON);
    const uint32_t v = lx < ly ? lx : ly;
    const int tv = (int)((tau >> v) & 1);
    uint64_t tx0, ty0, tx1, ty1;
    const LIMDD_TARG x0 = shadow(X, tau, v, tv, &tx0), y0 = shadow(Y, 0, v, 0, &ty0);
    const LIMDD_TARG x1 = shadow(X, tau, v, 1 - tv, &tx1), y1 = shadow(Y, 0, v, 1, &ty1);
    SPAWN(bqd_xpsm_sub, x0, y0, tx0 ^ ty0);
    const int s1 = CALL(bqd_xpsm_sub, x1, y1, tx1 ^ ty1);
    const int s0 = SYNC(bqd_xpsm_sub);
    const int r = s0 && s1;
    cache_put3(CACHE_BQD_XPSM_SUBSET, X, Y, tau, (uint64_t)r + 1);   /* 0 would read as a miss */
    return r;
}

/**
 * Dep: whether [N] as the diagram stores it depends on x_q, that is, whether a
 * node of variable q or a Z_q on a high label is reachable from N; an X_q is
 * a translation, which keeps a function free of x_q. The high label first,
 * then the low child, and the high one only where neither has it. Memoised on
 * (N, q).
 */
TASK_IMPL_2(int, bqd_xpsm_dep, LIMDD_TARG, N, uint32_t, q)
{
    if (N == 0 || N == LIMDD_TERMINAL) return 0;
    const uint32_t var = limdd_node_var(N);
    if (var > q) return 0;
    if (var == q) return 1;
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XPSM_DEP, N, q, 0, &hit)) return (int)hit - 1;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_CANON);
    const BQD hi = limdd_node_high(N);
    const int r = (!limdd_edge_is_zero(hi) && ((label_read(limdd_label(hi)).s >> q) & 1))
               || CALL(bqd_xpsm_dep, limdd_target(limdd_node_low(N)), q)
               || CALL(bqd_xpsm_dep, node_or_zero(hi), q);
    cache_put3(CACHE_BQD_XPSM_DEP, N, q, 0, (uint64_t)r + 1);
    return r;
}

/* --- cofactors ------------------------------------------------------------------ */

/**
 * Cof1: the canonical edge of [A] . [R] for a Q node, the product of its two
 * children with the identity, memoised on N. The high label is not applied to
 * R inside the product: a translation does not commute with it, and the
 * cofactor puts the label on outside, as the copy rule's Join.
 */
TASK_IMPL_1(BQD, bqd_xpsm_cof1, LIMDD_TARG, N)
{
    const BQD lo = limdd_node_low(N), hi = limdd_node_high(N);
    assert(!bqd_sm_is_s(N) && !limdd_edge_is_zero(hi));
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XPSM_COF1, N, 0, 0, &hit)) return (BQD)hit;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_COF1);
    const BQD r = CALL(bqd_xpsm_apply, BQD_SM_MUL, lo, unit(limdd_target(hi)));
    cache_put3(CACHE_BQD_XPSM_COF1, N, 0, 0, (uint64_t)r);
    return r;
}

/**
 * skip:lem:xcof, as bqd_xp_cofactor: the cofactor b of (c Z^s X^t, N) read
 * at v is l' = c (-1)^{b s_v} Z^{s'} X^{t'} applied to the cofactor b ^ t_v of
 * N, s' and t' being s and t without bit v, and N itself where N is below v.
 * The high side of a node at v is the stored edge for an S node, with no
 * product, and the high label times Cof1 for a Q node.
 */
TASK_IMPL_3(BQD, bqd_xpsm_cofactor, BQD, e, uint32_t, v, int, b)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG N = limdd_target(e);
    label_t l = label_read(limdd_label(e));
    const uint64_t t = l.t & ~above(limdd_level(N));
    const uint32_t top = top_var(N, l.s);
    assert(v <= top && "a labelled edge is read at or above its top");
    if (top > v) return t == l.t ? e : xedge(l, N);    /* e skips v */

    const bool sv = (l.s & bit(v)) != 0, tv = (t & bit(v)) != 0;
    l.s &= ~bit(v);
    l.t = t & ~bit(v);
    if (sv && b) l.c = neg(l.c);
    if (limdd_level(N) > v) return xedge(l, N);         /* a Z above the node */
    if ((b != 0) == tv) {
        const BQD lo = limdd_node_low(N);
        return limdd_edge_is_zero(lo) ? lo : xedge(l, limdd_target(lo));
    }
    const BQD hi = limdd_node_high(N);
    if (limdd_edge_is_zero(hi)) return hi;
    if (bqd_sm_is_s(N)) return relabel(l, hi);
    const BQD P = CALL(bqd_xpsm_cof1, N);
    const label_t h = label_mul(label_mul(l, label_read(limdd_label(hi))),
                                label_read(limdd_label(P)));
    return xedge(h, limdd_target(P));
}

/* --- constructors ----------------------------------------------------------------- */

/**
 * Case (v) of Compose (skip:alg:smxpcons), the two cofactors nonzero and
 * neither equal nor opposite, given a, the canonical edge of f_0, and
 * u = l_a^{-1} . f_1, a labelled edge: the node of l_a^{-1} f has the low
 * cofactor [N_a] and the high one u. With t the least point of supp u, the
 * shadow decides whether h = X^t u has its support inside supp N_a, with no
 * arithmetic: then the node is Q and its ratio k = h / [N_a], the RATIO of
 * Apply, and otherwise it is S and k is the canonical edge of h, one Canon,
 * which is free where h's label has no Z and no X. Either way 0 is the least
 * point of supp k, so k's label has no X, and its scalar is the value at the
 * pivot (1, t) of the level, which in the Pauli family gets the sign repair:
 * negated, with Z_v on the root label. The high edge is X^t . k.
 */
TASK_3(BQD, bqd_xpsm_compose_u, uint32_t, v, BQD, a, BQD, u)
{
    const label_t xt = xlabel(minpoint(limdd_target(u), label_read(limdd_label(u)).t));
    const BQD h = limdd_refs_push(relabel(xt, u));
    const LIMDD_TARG Na = limdd_target(a);
    const bool s = !CALL(bqd_xpsm_sub, limdd_target(h), Na, label_read(limdd_label(h)).t);
    const BQD k = s ? CALL(bqd_xpsm_canon, h) : CALL(bqd_xpsm_apply, BQD_SM_RATIO, unit(Na), h);
    limdd_refs_pop(1);
    label_t lk = label_read(limdd_label(k));
    assert(!limdd_edge_is_zero(k) && lk.t == 0 && "the ratio is nonzero at 0");
    label_t root = label_read(limdd_label(a));
    if (repair(lk.c)) {
        lk.c = neg(lk.c);
        root.s |= bit(v);
    }
    const LIMDD_TARG node = bqd_sm_mk(v, s, unit(Na), xedge(label_mul(xt, lk), limdd_target(k)));
    return limdd_bundle(label_lim(root), node);
}

/**
 * Compose of skip:alg:smxpcons, skip:prop:xcompose read under rule SM. a and
 * b are the canonical edges of f_0 and f_1.
 *
 *   (i)   a = b: f does not depend on x_v, and its edge is a, which skips v.
 *   (ii)  Pauli, b = -a: skipped with Z_v on a's label.
 *   (iii) a zero: the least point of supp f is (1, t_b), and translating by
 *         it moves f_1 to the low side: the Q node (v; (id, N_b), 0) under
 *         c_b Z^{s_b} X^{t_b + e_v}.
 *   (iv)  b zero: the Q node (v; (id, N_a), 0) under a's label.
 *   (v)   otherwise the node of l_a^{-1} f, compose_u.
 *
 * None of the nodes of (iii) to (v) is redundant.
 */
TASK_IMPL_3(BQD, bqd_xpsm_compose, uint32_t, v, BQD, a, BQD, b)
{
    if (a == b) return a;
    if (limdd_edge_is_zero(a)) {
        label_t lb = label_read(limdd_label(b));
        lb.t |= bit(v);
        const LIMDD_TARG node = bqd_sm_mk(v, false, unit(limdd_target(b)), limdd_zero_edge());
        return limdd_bundle(label_lim(lb), node);
    }
    const LIMDD_TARG na = limdd_target(a);
    if (limdd_edge_is_zero(b))
        return limdd_bundle(limdd_label(a), bqd_sm_mk(v, false, unit(na), limdd_zero_edge()));

    const LIMDD_TARG nb = limdd_target(b);
    label_t la = label_read(limdd_label(a));
    const label_t lb = label_read(limdd_label(b));
    if (bqd_family() == BQD_FAMILY_PAULI && na == nb && la.s == lb.s && la.t == lb.t
        && lb.c == neg(la.c)) {
        la.s |= bit(v);
        return limdd_bundle(label_lim(la), na);
    }
    const BQD u = limdd_refs_push(xedge(label_mul(label_inv(la), lb), nb));
    const BQD r = CALL(bqd_xpsm_compose_u, v, a, u);
    limdd_refs_pop(1);
    return r;
}

/**
 * VIEW of skip:lem:smxpview. With F = (c Z^s X^tau, N) read at v, lam the
 * label without bit v and N a Q node at v with tau_v = 0, the low cofactor is
 * lam . [A], and the high one is (-1)^{s_v} lam . (l_1 . ([A] . [R])), which
 * is X^{t_1} (F_0 . rho) for rho = c_rho Z^{s_1} X^{tau'} [R], tau' the X part
 * of lam: moving Z^{s_1} past X^{tau'}, and Z^{s'} Z^{s_1} past X^{t_1}, are
 * the signs c_rho takes on beside c_1 and (-1)^{s_v}. supp rho = supp R ^ tau'
 * lies inside supp A ^ tau' = supp F_0. A node below v is the virtual Q node
 * (N, Ind N): F_1 = (-1)^{s_v} F_0, whose ratio is the sign on the indicator
 * of supp F_0. An S node, or an X at v, has no view.
 */
TASK_IMPL_5(int, bqd_xpsm_view, BQD, F, uint32_t, v, BQD *, lo, uint64_t *, t, BQD *, rho)
{
    assert(!limdd_edge_is_zero(F));
    const label_t l = label_read(limdd_label(F));
    const LIMDD_TARG N = limdd_target(F);
    assert(top_var(N, l.s) >= v && "a labelled edge is read at or above its top");
    const bool sv = ((l.s >> v) & 1) != 0;
    const label_t lam = { l.c, l.s & lower(v), l.t & lower(v) };
    if (limdd_level(N) > v) {
        const label_t ri = { sv ? EVBDD_MIN_ONE : EVBDD_ONE, 0, lam.t };
        *lo = xedge(lam, N);
        *t = 0;
        *rho = xedge(ri, CALL(bqd_xpsm_ind, N));
        return 1;
    }
    if (((l.t >> v) & 1) || bqd_sm_is_s(N)) return 0;
    *lo = xedge(lam, limdd_target(limdd_node_low(N)));
    const BQD hi = limdd_node_high(N);
    if (limdd_edge_is_zero(hi)) {
        *t = 0;
        *rho = hi;
        return 1;
    }
    const label_t l1 = label_read(limdd_label(hi));
    const bool sign = sv ^ parity(l1.s & lam.t) ^ parity((lam.s ^ l1.s) & l1.t);
    const label_t r = { sign ? neg(l1.c) : l1.c, l1.s, lam.t };
    *t = l1.t;
    *rho = xedge(r, limdd_target(hi));
    return 1;
}

/**
 * NEST of skip:alg:smxpcons: the canonical edge of the function with the
 * cofactors [lo] and X^t ([lo] . [rho]). With lo = L . [N], L = c Z^s X^p,
 * L^{-1} f has the cofactors [N] and sigma X^t ([N] . X^p rho), with
 * sigma = (-1)^{s.t}. Where the least point of the support of that is t, its
 * move by t is sigma [N] . X^p rho, whose support is inside supp N, so the
 * node is Q and its ratio sigma X^p rho, canonical once X^p rho is, with no
 * product and no quotient; 0 is then its least point, the Pauli sign repair
 * goes on, and the level is skipped where the ratio is the identity on the
 * indicator of N, which in the Pauli family after the repair is f_1 = -f_0
 * with Z_v. Where the least point moved, Compose of the product decides the
 * node. A zero ratio is a node with no high cofactor.
 */
TASK_IMPL_5(BQD, bqd_xpsm_nest, uint32_t, v, BQD, lo, uint64_t, t, BQD, rho, int, rho_canon)
{
    if (limdd_edge_is_zero(lo)) return lo;
    if (limdd_edge_is_zero(rho)) return CALL(bqd_xpsm_compose, v, lo, rho);
    const label_t L = label_read(limdd_label(lo));
    const LIMDD_TARG N = limdd_target(lo);
    const BQD rp = limdd_refs_push(L.t != 0 ? CALL(bqd_xpsm_canon, relabel(xlabel(L.t), rho))
                                 : rho_canon ? rho : CALL(bqd_xpsm_canon, rho));
    label_t k = label_read(limdd_label(rp));
    const LIMDD_TARG R = limdd_target(rp);
    BQD r;
    if (minpoint(R, k.t ^ t) == t) {
        assert(k.t == 0 && "the ratio is nonzero at 0");
        if (parity(L.s & t)) k.c = neg(k.c);
        label_t root = L;
        if (repair(k.c)) {
            k.c = neg(k.c);
            root.s |= bit(v);
        }
        /* an indicator's node is never above the node it is of */
        if (t == 0 && k.c == EVBDD_ONE && k.s == 0 && limdd_level(R) >= limdd_level(N)
            && R == CALL(bqd_xpsm_ind, N)) {
            r = limdd_bundle(label_lim(root), N);
        } else {
            const BQD h = xedge(label_mul(xlabel(t), k), R);
            r = limdd_bundle(label_lim(root), bqd_sm_mk(v, false, unit(N), h));
        }
    } else {
        const label_t xt = xlabel(t);
        const BQD a = limdd_refs_push(relabel(xt, lo));
        const BQD b = limdd_refs_push(relabel(xt, rho));
        const BQD m = limdd_refs_push(CALL(bqd_xpsm_apply, BQD_SM_MUL, a, b));
        r = CALL(bqd_xpsm_compose, v, lo, m);
        limdd_refs_pop(3);
    }
    limdd_refs_pop(1);
    return r;
}

/**
 * RebuildS of skip:alg:smxpcons. The image of an S node under an affine
 * bijection y -> M y + w of the levels below v moves both supports by it, so
 * the misalignment stays, but the tag is decided after the move by the least
 * point, and the least point of the image is the image of the old one only
 * where it is M t_1 (skip:lem:smxptag). There the node is S again, with the
 * canonical edge of L^{-1} f_1, L the root label, and the Pauli repair of its
 * pivot; elsewhere Compose's case (v) decides, which may make a Q node.
 */
TASK_IMPL_4(BQD, bqd_xpsm_rebuild_s, uint32_t, v, BQD, lo, BQD, f1, uint64_t, tt)
{
    assert(!limdd_edge_is_zero(lo) && !limdd_edge_is_zero(f1));
    label_t root = label_read(limdd_label(lo));
    const BQD k = limdd_refs_push(CALL(bqd_xpsm_canon, relabel(label_inv(root), f1)));
    label_t lk = label_read(limdd_label(k));
    BQD r;
    if (lk.t != tt) {
        r = CALL(bqd_xpsm_compose_u, v, lo, k);
    } else {
        if (repair(lk.c)) {
            lk.c = neg(lk.c);
            root.s |= bit(v);
        }
        const LIMDD_TARG node = bqd_sm_mk(v, true, unit(limdd_target(lo)), xedge(lk, limdd_target(k)));
        r = limdd_bundle(label_lim(root), node);
    }
    limdd_refs_pop(1);
    return r;
}

/**
 * Canon of skip:alg:smxpcons: the canonical edge of a labelled edge. A label
 * with no Z and no X is a scalar on a stored node, canonical as it is, and in
 * the translation family so is c X^t where t is the least point of the
 * support (CanonT's test). Otherwise the memo, on the edge without its
 * scalar, which goes back on (skip:lem:xscale), and at the top v:
 *
 *   N at v, no X at v, Q    NEST of the view, with Canon of its low edge:
 *                           the node moved by the label, no high cofactor
 *   N at v, no X at v, S    RebuildS of the two cofactors, whose high one is
 *                           the stored edge, under the label's translation
 *   otherwise               Compose of the two canonical cofactors (an X at v
 *                           trades them, a Z above N signs one)
 *
 * The recursion runs on e without its scalar (with_scalar), and the scalar
 * goes on the result.
 */
TASK_IMPL_1(BQD, bqd_xpsm_canon, BQD, e)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG N = limdd_target(e);
    const label_t l = label_read(limdd_label(e));
    const uint64_t t = l.t & ~above(limdd_level(N));
    if (l.s == 0 && t == 0) return t == l.t ? e : xedge(l, N);
    if (bqd_family() == BQD_FAMILY_X && minpoint(N, t) == t) return t == l.t ? e : xedge(l, N);

    const uint64_t ke = key_edge(e);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XPSM_CANON, 0, ke, 0, &hit)) return xscale(l.c, (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_CANON);

    const uint32_t v = top_var(N, l.s);
    e = limdd_refs_push(with_scalar(e, EVBDD_ONE));
    BQD r;
    if (limdd_level(N) == v && !((t >> v) & 1) && !bqd_sm_is_s(N)) {
        BQD F0, rho;
        uint64_t t1;
        CALL(bqd_xpsm_view, e, v, &F0, &t1, &rho);
        limdd_refs_push(F0);
        limdd_refs_push(rho);
        const BQD lo = limdd_refs_push(CALL(bqd_xpsm_canon, F0));
        const int rc = limdd_edge_is_zero(rho) || label_read(limdd_label(rho)).t == 0;
        r = CALL(bqd_xpsm_nest, v, lo, t1, rho, rc);
        limdd_refs_pop(3);
    } else if (limdd_level(N) == v && !((t >> v) & 1)) {
        const BQD f0 = limdd_refs_push(CALL(bqd_xpsm_cofactor, e, v, 0));
        const BQD f1 = limdd_refs_push(CALL(bqd_xpsm_cofactor, e, v, 1));
        const BQD lo = limdd_refs_push(CALL(bqd_xpsm_canon, f0));
        r = CALL(bqd_xpsm_rebuild_s, v, lo, f1, label_read(limdd_label(limdd_node_high(N))).t);
        limdd_refs_pop(3);
    } else {
        const BQD e0 = limdd_refs_push(CALL(bqd_xpsm_cofactor, e, v, 0));
        limdd_refs_spawn(SPAWN(bqd_xpsm_canon, e0));
        const BQD e1 = limdd_refs_push(CALL(bqd_xpsm_cofactor, e, v, 1));
        const BQD hi = limdd_refs_push(CALL(bqd_xpsm_canon, e1));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_canon)));
        r = CALL(bqd_xpsm_compose, v, lo, hi);
        limdd_refs_pop(4);
    }
    cache_put3(CACHE_BQD_XPSM_CANON, 0, ke, 0, (uint64_t)r);
    limdd_refs_pop(1);
    return xscale(l.c, r);
}

/* --- pointwise operations ----------------------------------------------------------- */

static uint64_t
op_cache_id(int op)
{
    switch (op) {
    case BQD_SM_MUL:   return CACHE_BQD_XPSM_MUL;
    case BQD_SM_RATIO: return CACHE_BQD_XPSM_RATIO;
    default:           return CACHE_BQD_XPSM_ADD;
    }
}

/*
 * The ratio b / a of two scalars as b . (1 / a), as the scalar family's Apply
 * takes it at the terminal: a is the scalar of the ratio's denominator, a
 * node's ratio scalar that recurs, so its inverse is found in the weight
 * operation cache, where the quotient of the pair, on exact weights a gcd of
 * two large numbers, mostly is not. On clifford_T_circuit_20_700 it took the
 * Hadamards of the translation family from 1.2 times the scalar family's time
 * to 0.8 times, for the ratio of two constants and the scalar of the ratio
 * outside its memo together.
 */
static EVBDD_WGT
op_scalar(int op, EVBDD_WGT a, EVBDD_WGT b)
{
    switch (op) {
    case BQD_SM_MUL:   return wgt_mul(a, b);
    case BQD_SM_RATIO: return wgt_mul(wgt_div(EVBDD_ONE, a), b);
    default:           return wgt_add(a, b);
    }
}

/**
 * Apply of skip:alg:smxpapply. The scalars come out as far as the operation
 * allows, and what is left of them stays in the key, beside the two labelled
 * edges without theirs (key_edge):
 *
 *     mul     f.g     = c_f c_g . ((f/c_f)(g/c_g))      nothing left, the pair
 *                                                       ordered, since it commutes
 *     ratio   g/f     = (c_g/c_f) . ((g/c_g)/(f/c_f))   nothing left
 *     add     f + g   = c_f . (f/c_f + g/c_f)           c_g/c_f left
 *
 * The recursion is at v, the higher of the two tops, where an operand that
 * skips v is its own virtual node. It ends at two constants, and a constant
 * operand of a product, or the denominator of a ratio, is a scale of Canon of
 * the other. Before the cofactors, the views of skip:lem:smxpview. Where both
 * operands have one at v, with one translation t (or a zero ratio, which
 * fixes none), the high cofactors of a product are X^t (F_0 G_0 rho sigma), and
 * of a ratio X^t ((G_0/F_0)(sigma/rho)), so the result is NEST of the product,
 * or ratio, of the low edges and of the ratios: the paper's prop:prodscalar on
 * nested nodes, with no high cofactor and no quotient where the least point
 * stays. A sum of two views with one translation and one ratio rho is
 * X^t ((F_0 + G_0) rho), and its ratio is rho restricted to the support of the
 * new low cofactor. Every other pair, an S node on either side among them,
 * takes the labelled cofactors and Compose. The recursion runs on f / c_f
 * and on g with what is left of its scalar in the key (with_scalar), the
 * scalar family's F and G, and the outer scalar goes on the result; the two
 * calls of each case overlap.
 */
TASK_IMPL_3(BQD, bqd_xpsm_apply, int, op, BQD, f, BQD, g)
{
    const bool fz = limdd_edge_is_zero(f), gz = limdd_edge_is_zero(g);
    if (op == BQD_SM_ADD) {
        if (fz) return CALL(bqd_xpsm_canon, g);
        if (gz) return CALL(bqd_xpsm_canon, f);
    } else if (fz || gz) {
        return limdd_zero_edge();
    }

    const LIMDD_TARG nf = limdd_target(f), ng = limdd_target(g);
    const label_t lf = label_read(limdd_label(f)), lg = label_read(limdd_label(g));
    const bool kf = nf == LIMDD_TERMINAL && lf.s == 0;
    const bool kg = ng == LIMDD_TERMINAL && lg.s == 0;
    if (kf && kg) return constant(op_scalar(op, lf.c, lg.c));
    if (kf && op == BQD_SM_MUL) return xscale(lf.c, CALL(bqd_xpsm_canon, g));
    if (kg && op == BQD_SM_MUL) return xscale(lg.c, CALL(bqd_xpsm_canon, f));
    if (kf && op == BQD_SM_RATIO) return xscale(wgt_div(EVBDD_ONE, lf.c), CALL(bqd_xpsm_canon, g));

    EVBDD_WGT outer, kept = 0;
    switch (op) {
    case BQD_SM_MUL:   outer = wgt_mul(lf.c, lg.c);                             break;
    case BQD_SM_RATIO: outer = op_scalar(op, lf.c, lg.c);                      break;
    default:           outer = lf.c;               kept = wgt_div(lg.c, lf.c);  break;
    }
    assert(kept < (UINT64_C(1) << 40) && "a weight index fits under the op id");
    uint64_t kF = key_edge(f), kG = key_edge(g);
    if (op == BQD_SM_MUL && kF > kG) { const uint64_t k = kF; kF = kG; kG = k; }

    uint64_t hit;
    if (cache_get3(op_cache_id(op), kept, kF, kG, &hit)) return xscale(outer, (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_APPLY);

    const uint32_t tf = top_var(nf, lf.s), tg = top_var(ng, lg.s);
    const uint32_t v = tf < tg ? tf : tg;
    f = limdd_refs_push(with_scalar(f, EVBDD_ONE));
    g = limdd_refs_push(with_scalar(g, op == BQD_SM_ADD ? kept : EVBDD_ONE));
    BQD r, F0 = 0, rF = 0, G0 = 0, rG = 0;
    uint64_t tF = 0, tG = 0;
    const uint64_t xf = lf.t & ~above(limdd_level(nf)), xg = lg.t & ~above(limdd_level(ng));
    if (op == BQD_SM_RATIO && nf == ng && lf.s == lg.s && xf == xg) {
        /* g = (c_g/c_f) f: the ratio is that scalar on the indicator of supp f */
        const BQD I = limdd_refs_push(xedge(xlabel(xf), CALL(bqd_xpsm_ind, nf)));
        r = CALL(bqd_xpsm_canon, I);
        limdd_refs_pop(1);
    } else if (CALL(bqd_xpsm_view, f, v, &F0, &tF, &rF) && CALL(bqd_xpsm_view, g, v, &G0, &tG, &rG)
               && (op != BQD_SM_ADD ? limdd_edge_is_zero(rF) || limdd_edge_is_zero(rG) || tF == tG
                                    : tF == tG && rF == rG)) {
        limdd_refs_push(F0);
        limdd_refs_push(rF);
        limdd_refs_push(G0);
        limdd_refs_push(rG);
        if (op != BQD_SM_ADD) {
            const uint64_t t = limdd_edge_is_zero(rF) ? tG : tF;
            limdd_refs_spawn(SPAWN(bqd_xpsm_apply, op, F0, G0));
            const BQD hi = limdd_refs_push(CALL(bqd_xpsm_apply, op, rF, rG));
            const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_apply)));
            r = CALL(bqd_xpsm_nest, v, lo, t, hi, 1);
            limdd_refs_pop(2);
        } else {
            const BQD L = limdd_refs_push(CALL(bqd_xpsm_apply, BQD_SM_ADD, F0, G0));
            if (limdd_edge_is_zero(L)) {
                r = L;
            } else {
                /* the restriction is canonical where it is a product, and rF is
                 * where its translation is 0 */
                const BQD rr = limdd_refs_push(CALL(bqd_xpsm_restr, rF, L));
                const int rc = rr != rF || limdd_edge_is_zero(rr)
                            || label_read(limdd_label(rF)).t == 0;
                r = CALL(bqd_xpsm_nest, v, L, tF, rr, rc);
                limdd_refs_pop(1);
            }
            limdd_refs_pop(1);
        }
        limdd_refs_pop(4);
    } else {
        const BQD f0 = limdd_refs_push(CALL(bqd_xpsm_cofactor, f, v, 0));
        const BQD g0 = limdd_refs_push(CALL(bqd_xpsm_cofactor, g, v, 0));
        limdd_refs_spawn(SPAWN(bqd_xpsm_apply, op, f0, g0));
        const BQD f1 = limdd_refs_push(CALL(bqd_xpsm_cofactor, f, v, 1));
        const BQD g1 = limdd_refs_push(CALL(bqd_xpsm_cofactor, g, v, 1));
        const BQD hi = limdd_refs_push(CALL(bqd_xpsm_apply, op, f1, g1));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_apply)));
        r = CALL(bqd_xpsm_compose, v, lo, hi);
        limdd_refs_pop(6);
    }
    cache_put3(op_cache_id(op), kept, kF, kG, (uint64_t)r);
    limdd_refs_pop(2);
    return xscale(outer, r);
}

/**
 * Restr: rho . [Y != 0], for Y nonzero. Nothing to do where supp rho lies
 * inside supp Y already, which the shadow decides; otherwise the product
 * with the indicator of supp Y, X^{t_Y} on the indicator of its node.
 */
TASK_IMPL_2(BQD, bqd_xpsm_restr, BQD, rho, BQD, Y)
{
    assert(!limdd_edge_is_zero(Y));
    if (limdd_edge_is_zero(rho)) return rho;
    const uint64_t ty = label_read(limdd_label(Y)).t;
    if (CALL(bqd_xpsm_sub, limdd_target(rho), limdd_target(Y), label_read(limdd_label(rho)).t ^ ty))
        return rho;
    const BQD I = limdd_refs_push(xedge(xlabel(ty), CALL(bqd_xpsm_ind, limdd_target(Y))));
    const BQD r = CALL(bqd_xpsm_apply, BQD_SM_MUL, rho, I);
    limdd_refs_pop(1);
    return r;
}

/* --- selections ------------------------------------------------------------------- */

/*
 * Restrict, Project, Pair, X and Perm of skip:alg:smxpsel take each value of
 * their result from one value of their operands, so they commute with the
 * view of a Q node (skip:lem:smxpview) wherever the view's translation leaves
 * the qubits they act on alone: the images of F_0 and of the ratio are the
 * view of the image, and NEST makes the node, with no high cofactor and no
 * quotient where the least point stays. Elsewhere, an S node, an X at the
 * level, or a translation on one of the qubits, they recurse on the labelled
 * cofactors and Compose decides the node, as the copy rule's do. Each is
 * linear, so the memo is on the labelled edge without its scalar, which the
 * recursion runs without (with_scalar) and which goes on the result.
 */

/**
 * Restrict of skip:alg:smxpsel: e|_{x_q = b}, canonical, a function that
 * ignores x_q. At or below q the cofactor at q, which is e itself where e
 * skips q, through Canon. Above q, at a stored Q node whose view has no X at
 * q, NEST of the restrictions of the low edge and of the ratio, which still
 * nest; otherwise Compose of the restrictions of the two cofactors. Memoised
 * on e without its scalar, q and b in the first word; the two calls overlap.
 */
TASK_IMPL_3(BQD, bqd_xpsm_restrict, BQD, e, uint32_t, q, int, b)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG N = limdd_target(e);
    const label_t l = label_read(limdd_label(e));
    const uint32_t v = top_var(N, l.s);
    if (v > q) return CALL(bqd_xpsm_canon, e);
    const uint64_t key = ((uint64_t)q << 1) | (uint64_t)(b & 1), ke = key_edge(e);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XPSM_RESTRICT, key, ke, 0, &hit)) return xscale(l.c, (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    e = limdd_refs_push(with_scalar(e, EVBDD_ONE));
    BQD r, F0, rho;
    uint64_t t;
    if (v == q) {
        const BQD f = limdd_refs_push(CALL(bqd_xpsm_cofactor, e, q, b));
        r = CALL(bqd_xpsm_canon, f);
        limdd_refs_pop(1);
    } else if (limdd_level(N) == v && CALL(bqd_xpsm_view, e, v, &F0, &t, &rho) && !((t >> q) & 1)) {
        limdd_refs_push(F0);
        limdd_refs_push(rho);
        limdd_refs_spawn(SPAWN(bqd_xpsm_restrict, rho, q, b));
        const BQD lo = limdd_refs_push(CALL(bqd_xpsm_restrict, F0, q, b));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_restrict)));
        r = CALL(bqd_xpsm_nest, v, lo, t, hi, 1);
        limdd_refs_pop(4);
    } else {
        const BQD e0 = limdd_refs_push(CALL(bqd_xpsm_cofactor, e, v, 0));
        limdd_refs_spawn(SPAWN(bqd_xpsm_restrict, e0, q, b));
        const BQD e1 = limdd_refs_push(CALL(bqd_xpsm_cofactor, e, v, 1));
        const BQD hi = limdd_refs_push(CALL(bqd_xpsm_restrict, e1, q, b));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_restrict)));
        r = CALL(bqd_xpsm_compose, v, lo, hi);
        limdd_refs_pop(4);
    }
    cache_put3(CACHE_BQD_XPSM_RESTRICT, key, ke, 0, (uint64_t)r);
    limdd_refs_pop(1);
    return xscale(l.c, r);
}

/**
 * Project of skip:alg:smxpsel: e . [x_q = b], canonical. Where e skips q,
 * Canon of e composed with zero on the other side, before the memo; at q the
 * chosen cofactor through Canon, composed so; above q as Restrict, NEST of
 * the projections of a Q node's low edge and ratio, or Compose of those of
 * the two cofactors.
 */
TASK_IMPL_3(BQD, bqd_xpsm_project, BQD, e, uint32_t, q, int, b)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG N = limdd_target(e);
    const label_t l = label_read(limdd_label(e));
    const uint32_t v = top_var(N, l.s);
    if (v > q) {
        const BQD g = limdd_refs_push(CALL(bqd_xpsm_canon, e));
        const BQD r = b ? CALL(bqd_xpsm_compose, q, limdd_zero_edge(), g)
                        : CALL(bqd_xpsm_compose, q, g, limdd_zero_edge());
        limdd_refs_pop(1);
        return r;
    }
    const uint64_t key = ((uint64_t)q << 1) | (uint64_t)(b & 1), ke = key_edge(e);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XPSM_PROJECT, key, ke, 0, &hit)) return xscale(l.c, (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    e = limdd_refs_push(with_scalar(e, EVBDD_ONE));
    BQD r, F0, rho;
    uint64_t t;
    if (v == q) {
        const BQD f = limdd_refs_push(CALL(bqd_xpsm_cofactor, e, q, b));
        const BQD g = limdd_refs_push(CALL(bqd_xpsm_canon, f));
        r = b ? CALL(bqd_xpsm_compose, q, limdd_zero_edge(), g)
              : CALL(bqd_xpsm_compose, q, g, limdd_zero_edge());
        limdd_refs_pop(2);
    } else if (limdd_level(N) == v && CALL(bqd_xpsm_view, e, v, &F0, &t, &rho) && !((t >> q) & 1)) {
        limdd_refs_push(F0);
        limdd_refs_push(rho);
        limdd_refs_spawn(SPAWN(bqd_xpsm_project, rho, q, b));
        const BQD lo = limdd_refs_push(CALL(bqd_xpsm_project, F0, q, b));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_project)));
        r = CALL(bqd_xpsm_nest, v, lo, t, hi, 1);
        limdd_refs_pop(4);
    } else {
        const BQD e0 = limdd_refs_push(CALL(bqd_xpsm_cofactor, e, v, 0));
        limdd_refs_spawn(SPAWN(bqd_xpsm_project, e0, q, b));
        const BQD e1 = limdd_refs_push(CALL(bqd_xpsm_cofactor, e, v, 1));
        const BQD hi = limdd_refs_push(CALL(bqd_xpsm_project, e1, q, b));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_project)));
        r = CALL(bqd_xpsm_compose, v, lo, hi);
        limdd_refs_pop(4);
    }
    cache_put3(CACHE_BQD_XPSM_PROJECT, key, ke, 0, (uint64_t)r);
    limdd_refs_pop(1);
    return xscale(l.c, r);
}

/**
 * X on qubit q: Canon of the label product X_q . e. Canon moves a Q node by
 * NEST and an S node by RebuildS, which keeps the tag only where the least
 * point stays (skip:lem:smxptag); in the translation family a label whose
 * translation is the least point already is free.
 */
TASK_IMPL_2(BQD, bqd_xpsm_x, BQD, e, uint32_t, q)
{
    if (limdd_edge_is_zero(e)) return e;
    const BQD m = limdd_refs_push(relabel(xlabel(bit(q)), e));
    const BQD r = CALL(bqd_xpsm_canon, m);
    limdd_refs_pop(1);
    return r;
}

/**
 * Pair of skip:alg:smxpsel on labelled edges: [x_b = 0] A|_{x_b = s} +
 * [x_b = 1] B|_{x_b = t}, at the higher of the two tops, l. Below b both are
 * their own restrictions and the result is their Compose at b; at b Compose
 * of the two chosen cofactors. Above b, where both operands have a view at l
 * (a zero one has the view (0, 0, 0)) with one translation among the nonzero
 * ratios and no X at b in it, the pair of the low edges and the pair of the
 * ratios are the view of the result (skip:lem:smpair (d) with labels), and
 * NEST makes it; otherwise the pairs of the labelled cofactors and Compose.
 * Linear in the pair, so the scalar g of the first nonzero operand comes out
 * and the memo is the copy rule's: the two edges without it, B's scalar over
 * g, s and b in the first word, below 2^40, and t in the id.
 */
TASK_IMPL_5(BQD, bqd_xpsm_pair, BQD, A, int, s, BQD, B, int, t, uint32_t, b)
{
    const bool az = limdd_edge_is_zero(A), bz = limdd_edge_is_zero(B);
    if (az && bz) return A;
    const EVBDD_WGT g = az ? scalar_of(B) : scalar_of(A);
    const BQD A1 = with_scalar(A, EVBDD_ONE);
    const BQD B1 = bz ? B : az ? with_scalar(B, EVBDD_ONE) : with_scalar(B, wgt_div(scalar_of(B), g));
    const uint64_t kappa = bz ? 0 : (uint64_t)scalar_of(B1);
    assert(kappa < (UINT64_C(1) << 33) && "a weight index is below 2^33");
    const uint64_t k0 = kappa | ((uint64_t)(s & 1) << 33) | ((uint64_t)b << 34);
    const uint64_t kA = az ? 0 : key_edge(A1), kB = bz ? 0 : key_edge(B1);
    const uint64_t id = t ? CACHE_BQD_XPSM_PAIR1 : CACHE_BQD_XPSM_PAIR0;
    uint64_t hit;
    if (cache_get3(id, k0, kA, kB, &hit)) return xscale(g, (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    limdd_refs_push(A1);
    limdd_refs_push(B1);
    const uint32_t ta = bqd_xp_top(A1), tb = bqd_xp_top(B1);
    const uint32_t l = ta < tb ? ta : tb;
    BQD r;
    if (l > b) {
        const BQD ca = limdd_refs_push(CALL(bqd_xpsm_canon, A1));
        const BQD cb = limdd_refs_push(CALL(bqd_xpsm_canon, B1));
        r = CALL(bqd_xpsm_compose, b, ca, cb);
        limdd_refs_pop(2);
    } else if (l == b) {
        const BQD fa = limdd_refs_push(CALL(bqd_xpsm_cofactor, A1, b, s));
        const BQD fb = limdd_refs_push(CALL(bqd_xpsm_cofactor, B1, b, t));
        limdd_refs_spawn(SPAWN(bqd_xpsm_canon, fa));
        const BQD cb = limdd_refs_push(CALL(bqd_xpsm_canon, fb));
        const BQD ca = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_canon)));
        r = CALL(bqd_xpsm_compose, b, ca, cb);
        limdd_refs_pop(4);
    } else {
        BQD a0 = limdd_zero_edge(), ra = a0, b0 = a0, rb = a0;
        uint64_t sa = 0, sb = 0;
        const bool viewed = (az || CALL(bqd_xpsm_view, A1, l, &a0, &sa, &ra))
                         && (bz || CALL(bqd_xpsm_view, B1, l, &b0, &sb, &rb));
        const bool za = limdd_edge_is_zero(ra), zb = limdd_edge_is_zero(rb);
        const uint64_t tt = za ? sb : sa;                /* 0 where both ratios are zero */
        if (viewed && (za || zb || sa == sb) && !((tt >> b) & 1)) {
            limdd_refs_push(a0);
            limdd_refs_push(ra);
            limdd_refs_push(b0);
            limdd_refs_push(rb);
            limdd_refs_spawn(SPAWN(bqd_xpsm_pair, ra, s, rb, t, b));
            const BQD lo = limdd_refs_push(CALL(bqd_xpsm_pair, a0, s, b0, t, b));
            const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_pair)));
            r = CALL(bqd_xpsm_nest, l, lo, tt, hi, 1);
            limdd_refs_pop(6);
        } else {
            const BQD a0c = limdd_refs_push(CALL(bqd_xpsm_cofactor, A1, l, 0));
            const BQD b0c = limdd_refs_push(CALL(bqd_xpsm_cofactor, B1, l, 0));
            limdd_refs_spawn(SPAWN(bqd_xpsm_pair, a0c, s, b0c, t, b));
            const BQD a1c = limdd_refs_push(CALL(bqd_xpsm_cofactor, A1, l, 1));
            const BQD b1c = limdd_refs_push(CALL(bqd_xpsm_cofactor, B1, l, 1));
            const BQD hi = limdd_refs_push(CALL(bqd_xpsm_pair, a1c, s, b1c, t, b));
            const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_pair)));
            r = CALL(bqd_xpsm_compose, l, lo, hi);
            limdd_refs_pop(6);
        }
    }
    cache_put3(id, k0, kA, kB, (uint64_t)r);
    limdd_refs_pop(2);
    return xscale(g, r);
}

/**
 * PermN of skip:alg:smxpsel: the canonical edge of pi [N], memoised on N, the
 * kind and the two qubits. Below qa, [N] does not depend on x_qa; at qa the
 * new cofactors are Pairs at qb, or the low one and X_qb of the high one, as
 * the copy rule's. Above qa pi is a bijection of the levels below, and the
 * children are mapped: a Q node's cofactors become pi [A] and
 * l_1^pi . (pi [A] . pi [R]), which is X^tau ([lo] . (X^tau l_1^pi) [W]) for
 * tau = pi t_1, so NEST of the two images, with no product where the least
 * point stays; an S node's high cofactor becomes l_1^pi . [W], and RebuildS
 * keeps its tag where the least point is still tau.
 */
TASK_4(BQD, bqd_xpsm_perm_node, LIMDD_TARG, N, uint32_t, kind, uint32_t, qa, uint32_t, qb)
{
    const uint32_t var = limdd_level(N);
    if (var > qa && kind == BQD_PERM_CX_UP) return unit(N);
    const uint64_t key = (uint64_t)kind | ((uint64_t)qa << 8) | ((uint64_t)qb << 16);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XPSM_PERM, N, key, 0, &hit)) return (BQD)hit;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    const BQD E = unit(N);
    BQD r;
    if (var > qa) {
        BQD g0, g1;
        if (kind == BQD_PERM_CX_DOWN) {
            g0 = limdd_refs_push(E);
            g1 = limdd_refs_push(CALL(bqd_xpsm_x, E, qb));
        } else {
            limdd_refs_spawn(SPAWN(bqd_xpsm_restrict, E, qb, 1));
            g0 = limdd_refs_push(CALL(bqd_xpsm_restrict, E, qb, 0));
            g1 = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_restrict)));
        }
        r = CALL(bqd_xpsm_compose, qa, g0, g1);
        limdd_refs_pop(2);
    } else if (var == qa) {
        const BQD f0 = limdd_refs_push(CALL(bqd_xpsm_cofactor, E, qa, 0));
        const BQD f1 = limdd_refs_push(CALL(bqd_xpsm_cofactor, E, qa, 1));
        BQD g0, g1;
        if (kind == BQD_PERM_CX_DOWN) {
            limdd_refs_spawn(SPAWN(bqd_xpsm_canon, f0));
            g1 = limdd_refs_push(CALL(bqd_xpsm_x, f1, qb));
            g0 = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_canon)));
        } else {
            /* SW: (Pair(f0,0; f1,0), Pair(f0,1; f1,1)); CX up: (Pair(f0,0; f1,1), Pair(f1,0; f0,1)) */
            const bool sw = kind == BQD_PERM_SWAP;
            limdd_refs_spawn(SPAWN(bqd_xpsm_pair, sw ? f0 : f1, sw ? 1 : 0, sw ? f1 : f0, 1, qb));
            g0 = limdd_refs_push(CALL(bqd_xpsm_pair, f0, 0, f1, sw ? 0 : 1, qb));
            g1 = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_pair)));
        }
        r = CALL(bqd_xpsm_compose, qa, g0, g1);
        limdd_refs_pop(4);
    } else {
        const BQD high = limdd_node_high(N);
        const LIMDD_TARG A = limdd_target(limdd_node_low(N));
        if (limdd_edge_is_zero(high)) {
            const BQD lo = limdd_refs_push(CALL(bqd_xpsm_perm_node, A, kind, qa, qb));
            r = CALL(bqd_xpsm_compose, var, lo, high);
            limdd_refs_pop(1);
        } else {
            limdd_refs_spawn(SPAWN(bqd_xpsm_perm_node, A, kind, qa, qb));
            const BQD W = limdd_refs_push(CALL(bqd_xpsm_perm_node, limdd_target(high), kind, qa, qb));
            const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_perm_node)));
            const label_t L = perm_label(kind, qa, qb, label_read(limdd_label(high)));
            if (bqd_sm_is_s(N)) {
                const BQD f1 = limdd_refs_push(relabel(L, W));
                r = CALL(bqd_xpsm_rebuild_s, var, lo, f1, L.t);
            } else {
                /* X^tau L has no X part: the ratio's label, which NEST puts through Canon */
                const BQD rho = limdd_refs_push(relabel(label_mul(xlabel(L.t), L), W));
                r = CALL(bqd_xpsm_nest, var, lo, L.t, rho, 0);
            }
            limdd_refs_pop(3);
        }
    }
    cache_put3(CACHE_BQD_XPSM_PERM, N, key, 0, (uint64_t)r);
    return r;
}

/** Perm: Canon(l^pi . PermN(N)) for e = (l, N), the conjugated label on the image of the node. */
TASK_IMPL_4(BQD, bqd_xpsm_perm, BQD, e, uint32_t, kind, uint32_t, qa, uint32_t, qb)
{
    if (limdd_edge_is_zero(e)) return e;
    limdd_refs_push(e);
    const label_t le = perm_label(kind, qa, qb, label_read(limdd_label(e)));
    const BQD p = limdd_refs_push(CALL(bqd_xpsm_perm_node, limdd_target(e), kind, qa, qb));
    const BQD m = limdd_refs_push(relabel(le, p));
    const BQD r = CALL(bqd_xpsm_canon, m);
    limdd_refs_pop(3);
    return r;
}

/* --- gates ------------------------------------------------------------------------- */

/**
 * Whether a ratio rho does not depend on x_q, so that a gate on q passes
 * through it: no Z_q on its label, and no node of q and no Z_q on a high label
 * reachable from its node (Dep). A zero ratio depends on nothing.
 */
TASK_2(int, bqd_xpsm_free_of, BQD, rho, uint32_t, q)
{
    if (limdd_edge_is_zero(rho)) return 1;
    if ((label_read(limdd_label(rho)).s >> q) & 1) return 0;
    return !CALL(bqd_xpsm_dep, limdd_target(rho), q);
}

/**
 * Side: whether supp [N] lies in {x_q = 0}, the scalar family's walk on the
 * shadow of skip:lem:smxpshadow. 0 is in the support of every stored node, so
 * that is the only side it can lie on, and an edge (l, N) lies on the side t_q
 * of its label where N does. The support is {0} x supp [A] and
 * {1} x (supp [R] ^ t_1), so it lies there where A's does and, for a nonzero
 * high edge, R's does and t_1 has no X_q. A node of variable q lies there
 * exactly where its high edge is zero; a node below q, the terminal and every
 * full node do not. No scalar, no sign and no product: the targets and the X
 * part of the high label. The low child first, and the high one only where
 * the low one lies there; memoised on (N, q).
 */
TASK_IMPL_2(int, bqd_xpsm_side, LIMDD_TARG, N, uint32_t, q)
{
    if (N == LIMDD_TERMINAL || limdd_node_var(N) > q || bqd_sm_full(N)) return 0;
    const BQD hi = limdd_node_high(N);
    if (limdd_node_var(N) == q) return limdd_edge_is_zero(hi);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XPSM_SIDE, N, q, 0, &hit)) return (int)hit - 1;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_CANON);
    int r = CALL(bqd_xpsm_side, limdd_target(limdd_node_low(N)), q);
    if (r && !limdd_edge_is_zero(hi))
        r = !((label_read(limdd_label(hi)).t >> q) & 1) && CALL(bqd_xpsm_side, limdd_target(hi), q);
    cache_put3(CACHE_BQD_XPSM_SIDE, N, q, 0, (uint64_t)r + 1);    /* 0 would read as a miss */
    return r;
}

/**
 * Gate of skip:alg:smxpapply, the copy rule's Gate at the top v of a labelled
 * edge, with the scalar family's one-sided gate and pass-through
 * (skip:lem:smpass): U = gates[gid] on qubit q where the pending controls of
 * cmask, all above q, are 1, and p the first of them, or q.
 *
 *   e skips p      a control: (Canon e, gate(e)) at p; the target: the row
 *                  sums of U times Canon e at q
 *   v above p      Side, where no control is pending and supp e lies in
 *                  {x_q = s}: e's node in {x_q = 0}, moved to s = t_q by
 *                  the label. Ue is u_0s g at x_q = 0 and u_1s g at x_q = 1
 *                  for the restriction g = e|_{x_q = s}, so u_0s g .
 *                  beta^{x_q} with beta = u_1s / u_0s, PhaseMul of the
 *                  exponent x_q, which keeps every ratio of g; where
 *                  u_0s u_1s is not zero, and in the Pauli family where
 *                  beta is a power of w_8, the phases that family
 *                  multiplies by. g is Canon of the label without q on the
 *                  restriction of the bare node, which carries no
 *                  translation down and shares its memo across labels.
 *                  Otherwise,
 *                  where e has a view at v whose translation has no X at q
 *                  or at a pending control, and whose ratio does not depend
 *                  on x_q, the gate commutes with the translation and with
 *                  the product by the ratio: Y = gate(F_0), and NEST of Y and
 *                  the ratio restricted to supp Y (Restr), with no high
 *                  cofactor and no sum; otherwise the gate on both
 *                  labelled cofactors and Compose
 *   v at a control the low cofactor through Canon, the gate without that
 *                  control on the high one
 *   v at q         the rows of U applied to the two cofactors, two sums
 *
 * Linear, so the memo is on e without its scalar, the gate's id first and
 * cmask with bit q set last, as the copy rule's, and the recursion runs
 * without the scalar.
 */
TASK_IMPL_4(BQD, bqd_xpsm_gate, BQD, e, uint32_t, gid, uint64_t, cmask, uint32_t, q)
{
    if (limdd_edge_is_zero(e)) return e;
    const label_t l = label_read(limdd_label(e));
    const uint32_t v = top_var(limdd_target(e), l.s);
    const uint64_t key = cmask | bit(q);
    const uint32_t p = (uint32_t)__builtin_ctzll(key);
    const uint64_t ke = key_edge(e);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XPSM_GATE, gid, ke, key, &hit)) return xscale(l.c, (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    e = limdd_refs_push(with_scalar(e, EVBDD_ONE));
    const EVBDD_WGT u00 = gates[gid][0], u01 = gates[gid][1];
    const EVBDD_WGT u10 = gates[gid][2], u11 = gates[gid][3];
    /* Side: the side s of x_q that supp e lies on, 2 where it does not apply */
    int s = 2;
    EVBDD_WGT beta = EVBDD_ONE;
    if (v < p && p == q && CALL(bqd_xpsm_side, limdd_target(e), q)) {
        const int side = (int)((l.t >> q) & 1);
        const EVBDD_WGT u0s = gates[gid][side], u1s = gates[gid][2 + side];
        if (u0s != EVBDD_ZERO && u1s != EVBDD_ZERO) {
            beta = wgt_div(u1s, u0s);
            if (bqd_family() == BQD_FAMILY_X || bqd_xp_w8_log(beta) >= 0) s = side;
        }
    }
    BQD r, F0, rho;
    uint64_t t;
    if (v > p && p < q) {
        limdd_refs_spawn(SPAWN(bqd_xpsm_canon, e));
        const BQD hi = limdd_refs_push(CALL(bqd_xpsm_gate, e, gid, cmask & ~bit(p), q));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_canon)));
        r = CALL(bqd_xpsm_compose, p, lo, hi);
        limdd_refs_pop(2);
    } else if (v > p) {
        const BQD C = limdd_refs_push(CALL(bqd_xpsm_canon, e));
        const BQD lo = limdd_refs_push(xscale(wgt_add(u00, u01), C));
        r = CALL(bqd_xpsm_compose, q, lo, xscale(wgt_add(u10, u11), C));
        limdd_refs_pop(2);
    } else if (s != 2) {
        /* e|_{x_q = s} = l' . [N]|_{x_q = 0}, l' the label Z^{s_e} X^{t_e} of e
         * without bit q and with the sign of Z_q at x_q = s */
        label_t lr = l;
        lr.c = (((l.s >> q) & 1) && s) ? EVBDD_MIN_ONE : EVBDD_ONE;
        lr.s &= ~bit(q);
        lr.t &= ~bit(q);
        const BQD g0 = limdd_refs_push(CALL(bqd_xpsm_restrict, unit(limdd_target(e)), q, 0));
        const BQD g1 = limdd_refs_push(relabel(lr, g0));
        const BQD g = limdd_refs_push(CALL(bqd_xpsm_canon, g1));
        const uint32_t ord = bqd_family() == BQD_FAMILY_PAULI ? 8 : bqd_exp_order(beta);
        const BQD_EXP eps = mtbdd_refs_push(bqd_exp_monomial(bit(q), 1, ord));
        r = xscale(gates[gid][s], CALL(bqd_xpsm_phase_mul, g, eps, beta));
        mtbdd_refs_pop(1);
        limdd_refs_pop(3);
    } else if (v < p && CALL(bqd_xpsm_view, e, v, &F0, &t, &rho) && (t & key) == 0
               && CALL(bqd_xpsm_free_of, rho, q)) {
        limdd_refs_push(F0);
        limdd_refs_push(rho);
        const BQD Y = limdd_refs_push(CALL(bqd_xpsm_gate, F0, gid, cmask, q));
        if (limdd_edge_is_zero(Y)) {
            r = Y;
        } else {
            /* the restriction is canonical where it is a product, and rho is
             * where its translation is 0 */
            const BQD rr = limdd_refs_push(CALL(bqd_xpsm_restr, rho, Y));
            const int rc = rr != rho || limdd_edge_is_zero(rr) || label_read(limdd_label(rho)).t == 0;
            r = CALL(bqd_xpsm_nest, v, Y, t, rr, rc);
            limdd_refs_pop(1);
        }
        limdd_refs_pop(3);
    } else if (v < p) {
        const BQD e0 = limdd_refs_push(CALL(bqd_xpsm_cofactor, e, v, 0));
        limdd_refs_spawn(SPAWN(bqd_xpsm_gate, e0, gid, cmask, q));
        const BQD e1 = limdd_refs_push(CALL(bqd_xpsm_cofactor, e, v, 1));
        const BQD hi = limdd_refs_push(CALL(bqd_xpsm_gate, e1, gid, cmask, q));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_gate)));
        r = CALL(bqd_xpsm_compose, v, lo, hi);
        limdd_refs_pop(4);
    } else if (v < q) {
        const BQD e0 = limdd_refs_push(CALL(bqd_xpsm_cofactor, e, v, 0));
        limdd_refs_spawn(SPAWN(bqd_xpsm_canon, e0));
        const BQD e1 = limdd_refs_push(CALL(bqd_xpsm_cofactor, e, v, 1));
        const BQD hi = limdd_refs_push(CALL(bqd_xpsm_gate, e1, gid, cmask & ~bit(v), q));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_canon)));
        r = CALL(bqd_xpsm_compose, v, lo, hi);
        limdd_refs_pop(4);
    } else {
        /* the two new cofactors are independent sums of the old ones */
        const BQD f0 = limdd_refs_push(CALL(bqd_xpsm_cofactor, e, q, 0));
        const BQD f1 = limdd_refs_push(CALL(bqd_xpsm_cofactor, e, q, 1));
        const BQD a0 = limdd_refs_push(xscale(u00, f0));
        const BQD a1 = limdd_refs_push(xscale(u01, f1));
        limdd_refs_spawn(SPAWN(bqd_xpsm_apply, BQD_SM_ADD, a0, a1));
        const BQD b0 = limdd_refs_push(xscale(u10, f0));
        const BQD b1 = limdd_refs_push(xscale(u11, f1));
        const BQD hi = limdd_refs_push(CALL(bqd_xpsm_apply, BQD_SM_ADD, b0, b1));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_apply)));
        r = CALL(bqd_xpsm_compose, q, lo, hi);
        limdd_refs_pop(8);
    }
    cache_put3(CACHE_BQD_XPSM_GATE, gid, ke, key, (uint64_t)r);
    limdd_refs_pop(1);
    return xscale(l.c, r);
}

/* --- phases and diagonal gates ------------------------------------------------------ */

/*
 * A phase multiplies by a function with no zero, which keeps every support,
 * and so every least point and every tag (skip:alg:smxpphase): the root of a
 * stored node keeps its translation 0 and its high edge its t_1, a Q node's
 * ratio is multiplied by the phase's quotient and an S node's high cofactor by
 * the phase itself, and the low result carries the scalar (and the sign
 * pattern, in the Pauli family) that the high one is divided by. No Sel, no
 * indicator exponent, no least point and no subset test, which the copy rule's
 * PhaseMul needs to multiply its copies off.
 */

/**
 * N read at level l, (tag, A, l_1, R): its own tag, low target, high label and
 * high target (0 for a zero high edge) where N is at l; where it is below, the
 * terminal included, the virtual Q node (N, Ind N) with the identity, whose
 * R the caller makes.
 */
static inline bool
split(LIMDD_TARG N, uint32_t l, bool *s, LIMDD_TARG *A, label_t *l1, LIMDD_TARG *R)
{
    const label_t id = { EVBDD_ONE, 0, 0 };
    *l1 = id;
    *s = false;
    *A = N;
    *R = 0;
    if (limdd_level(N) > l) return false;
    const BQD hi = limdd_node_high(N);
    *s = bqd_sm_is_s(N);
    *A = limdd_target(limdd_node_low(N));
    *R = node_or_zero(hi);
    if (*R != 0) *l1 = label_read(limdd_label(hi));
    return true;
}

/**
 * PhaseMulX of skip:alg:smxpphase, the translation family: the canonical edge
 * of [N] . beta^eps, r the order of beta, at l, the higher of N's level and
 * eps's. The low call gives (alpha, U), 0 being the least point of [A] and so
 * of U. With E = X^{t_1} eps_1, a Q node's ratio becomes c_1 [R] beta^{E -
 * eps_0} and an S node's high cofactor X^{t_1} (c_1 / alpha) [R] beta^E,
 * moved by the same t_1; a Q node whose new ratio is the indicator of U with
 * no translation is skipped. The terminal is beta^eps, Exp of the scalar
 * family, whose diagram this family shares on full support. Memoised on N,
 * eps and beta; the two calls depend on the input alone, and overlap.
 */
TASK_IMPL_4(BQD, bqd_xpsm_phase_mul_x, LIMDD_TARG, N, BQD_EXP, eps, EVBDD_WGT, beta, uint32_t, r)
{
    if (bqd_exp_is_const(eps)) return xscale(bqd_exp_power(beta, bqd_exp_value(eps)), unit(N));
    if (N == LIMDD_TERMINAL) return CALL(bqd_sm_exp, eps, beta, r);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XPSM_PHASEMUL, N, eps, beta, &hit)) return (BQD)hit;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PHASEMUL);

    const uint32_t lN = limdd_level(N), le = bqd_exp_var(eps);
    const uint32_t l = lN < le ? lN : le;
    const BQD_EXP e0 = bqd_exp_cof(eps, l, 0), e1 = bqd_exp_cof(eps, l, 1);
    bool s;
    LIMDD_TARG A, R;
    label_t l1;
    if (!split(N, l, &s, &A, &l1, &R)) R = CALL(bqd_xpsm_ind, N);
    BQD res;
    if (R == 0) {
        const BQD lo = CALL(bqd_xpsm_phase_mul_x, A, e0, beta, r);
        res = limdd_bundle(limdd_label(lo),
                           bqd_sm_mk(l, false, unit(limdd_target(lo)), limdd_zero_edge()));
    } else {
        const BQD_EXP E = mtbdd_refs_push(CALL(bqd_exp_translate, e1, l1.t));
        const BQD_EXP d = s ? E : mtbdd_refs_push(CALL(bqd_exp_sub, E, e0, r));
        limdd_refs_spawn(SPAWN(bqd_xpsm_phase_mul_x, A, e0, beta, r));
        const BQD K = limdd_refs_push(CALL(bqd_xpsm_phase_mul_x, R, d, beta, r));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_phase_mul_x)));
        const LIMDD_TARG U = limdd_target(lo);
        const BQD k = xscale(s ? wgt_div(l1.c, scalar_of(lo)) : l1.c, K);
        if (!s && l1.t == 0 && limdd_label(k) == LIMDD_LIM_IDENTITY
            && limdd_target(k) == CALL(bqd_xpsm_ind, U)) {
            res = lo;                                   /* f_1 = f_0: the level is skipped */
        } else {
            const BQD h = xedge(label_mul(xlabel(l1.t), label_read(limdd_label(k))), limdd_target(k));
            res = limdd_bundle(limdd_label(lo), bqd_sm_mk(l, s, unit(U), h));
        }
        limdd_refs_pop(2);
        mtbdd_refs_pop(s ? 1 : 2);
    }
    cache_put3(CACHE_BQD_XPSM_PHASEMUL, N, eps, beta, (uint64_t)res);
    return res;
}

/**
 * PhaseMulP of skip:alg:smxpphase, the Pauli family: the canonical edge of
 * [N] . w^eps, w = e^{i pi/4} and eps modulo 8. The low result is
 * (c_0 Z^{s_0}, U), and moving X^{t_1} past Z^w, w = s_0 xor s_1, costs
 * (-1)^{w . t_1}: a Q node's ratio becomes c_1 (-1)^{w.t_1} [R] w^{E - eps_0 +
 * P_{s_1}}, independent of the low result, and an S node's high cofactor
 * (c_1 / c_0) (-1)^{w.t_1} [R] w^{E + P_w}, which needs its sign pattern, so
 * the low call comes first there. A Q node whose new ratio is plus or minus
 * the indicator of U, with no translation, is skipped, with Z_l for minus,
 * and a pivot without an argument in [0, pi) gets the sign repair. The
 * terminal is its own virtual node, (T, (id, T)). Memoised on N and eps.
 */
TASK_IMPL_2(BQD, bqd_xpsm_phase_mul_p, LIMDD_TARG, N, BQD_EXP, eps)
{
    if (bqd_exp_is_const(eps)) {
        const EVBDD_WGT w = gates[GATEID_T][3];
        return xscale(bqd_exp_power(w, bqd_exp_value(eps)), unit(N));
    }
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XPSM_PHASEMULP, N, eps, 0, &hit)) return (BQD)hit;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PHASEMUL);

    const uint32_t lN = limdd_level(N), le = bqd_exp_var(eps);
    const uint32_t l = lN < le ? lN : le;
    const BQD_EXP e0 = bqd_exp_cof(eps, l, 0), e1 = bqd_exp_cof(eps, l, 1);
    bool s;
    LIMDD_TARG A, R;
    label_t l1;
    if (!split(N, l, &s, &A, &l1, &R)) R = CALL(bqd_xpsm_ind, N);
    BQD res;
    if (R == 0) {
        const BQD lo = CALL(bqd_xpsm_phase_mul_p, A, e0);
        res = limdd_bundle(limdd_label(lo),
                           bqd_sm_mk(l, false, unit(limdd_target(lo)), limdd_zero_edge()));
    } else {
        const BQD_EXP E = mtbdd_refs_push(CALL(bqd_exp_translate, e1, l1.t));
        BQD lo, K;
        label_t l0;
        EVBDD_WGT gam;
        if (!s) {
            const BQD_EXP d1 = mtbdd_refs_push(CALL(bqd_exp_sub, E, e0, 8));
            const BQD_EXP ps = mtbdd_refs_push(bqd_exp_parity4(l1.s));
            const BQD_EXP d = mtbdd_refs_push(CALL(bqd_exp_add, d1, ps, 8));
            limdd_refs_spawn(SPAWN(bqd_xpsm_phase_mul_p, A, e0));
            K = limdd_refs_push(CALL(bqd_xpsm_phase_mul_p, R, d));
            lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xpsm_phase_mul_p)));
            l0 = label_read(limdd_label(lo));
            gam = l1.c;
        } else {
            lo = limdd_refs_push(CALL(bqd_xpsm_phase_mul_p, A, e0));
            l0 = label_read(limdd_label(lo));
            const BQD_EXP pw = mtbdd_refs_push(bqd_exp_parity4(l0.s ^ l1.s));
            const BQD_EXP d = mtbdd_refs_push(CALL(bqd_exp_add, E, pw, 8));
            K = limdd_refs_push(CALL(bqd_xpsm_phase_mul_p, R, d));
            gam = wgt_div(l1.c, l0.c);
        }
        if (parity((l0.s ^ l1.s) & l1.t)) gam = neg(gam);
        const BQD k = xscale(gam, K);
        label_t lk = label_read(limdd_label(k));
        assert(lk.t == 0 && "0 is the least point of a node's support");
        const LIMDD_TARG U = limdd_target(lo);
        label_t root = { l0.c, l0.s, 0 };
        if (!s && l1.t == 0 && lk.s == 0 && (lk.c == EVBDD_ONE || lk.c == EVBDD_MIN_ONE)
            && limdd_target(k) == CALL(bqd_xpsm_ind, U)) {
            /* f_1 = f_0, skipped; f_1 = -f_0, skipped with Z_l */
            if (lk.c == EVBDD_ONE) {
                res = lo;
            } else {
                root.s |= bit(l);
                res = limdd_bundle(label_lim(root), U);
            }
        } else {
            if (repair(lk.c)) {
                lk.c = neg(lk.c);
                root.s |= bit(l);
            }
            const BQD h = xedge(label_mul(xlabel(l1.t), lk), limdd_target(k));
            res = limdd_bundle(label_lim(root), bqd_sm_mk(l, s, unit(U), h));
        }
        limdd_refs_pop(2);
        mtbdd_refs_pop(s ? 3 : 4);
    }
    cache_put3(CACHE_BQD_XPSM_PHASEMULP, N, eps, 0, (uint64_t)res);
    return res;
}

/** m eps modulo 8, by m - 1 sums: the exponent of w_8 that beta^eps is for beta = w_8^m. */
TASK_2(BQD_EXP, bqd_xpsm_exp_times, BQD_EXP, eps, int, m)
{
    BQD_EXP me = mtbdd_refs_push(bqd_exp_const(0, 8));
    for (int k = 0; k < m; k++) {
        const BQD_EXP next = CALL(bqd_exp_add, me, eps, 8);
        mtbdd_refs_pop(1);
        me = mtbdd_refs_push(next);
    }
    mtbdd_refs_pop(1);
    return me;
}

/** The exponent of w_8 that beta is, or an exit: the Pauli family multiplies by those only. */
static int
w8_exponent(EVBDD_WGT beta)
{
    const int m = bqd_xp_w8_log(beta);
    if (m < 0) {
        fprintf(stderr, "sylvan: the Pauli-BQD multiplies by a power of w_8 only\n");
        exit(1);
    }
    return m;
}

/**
 * e . beta^eps for a labelled edge e = (l, N): l . PhaseMulX(N, X^t eps,
 * beta) in the translation family, and l . PhaseMulP(N, X^t (m eps)) in the
 * Pauli family for beta = w^m. The product label is canonical where e is,
 * since the multiplication keeps the support, and with it the least point
 * and the pivots. eps is protected by the caller.
 */
TASK_IMPL_3(BQD, bqd_xpsm_phase_mul, BQD, e, BQD_EXP, eps, EVBDD_WGT, beta)
{
    if (limdd_edge_is_zero(e) || beta == EVBDD_ONE) return e;
    const label_t l = label_read(limdd_label(e));
    const LIMDD_TARG N = limdd_target(e);
    BQD p;
    if (bqd_family() == BQD_FAMILY_X) {
        const uint32_t r = bqd_exp_order(beta);
        const BQD_EXP moved = mtbdd_refs_push(CALL(bqd_exp_translate, eps, l.t));
        p = CALL(bqd_xpsm_phase_mul_x, N, moved, beta, r);
        mtbdd_refs_pop(1);
    } else {
        const BQD_EXP me = mtbdd_refs_push(CALL(bqd_xpsm_exp_times, eps, w8_exponent(beta)));
        const BQD_EXP moved = mtbdd_refs_push(CALL(bqd_exp_translate, me, l.t));
        p = CALL(bqd_xpsm_phase_mul_p, N, moved);
        mtbdd_refs_pop(2);
    }
    return xedge(label_mul(l, label_read(limdd_label(p))), limdd_target(p));
}

/**
 * beta^eps, of full support: Exp of the scalar family in the translation
 * family, whose diagram on full support is that family's, and PhaseMulP on
 * the terminal of m eps in the Pauli family, for beta = w^m.
 */
TASK_IMPL_3(BQD, bqd_xpsm_exp_state, BQD_EXP, eps, EVBDD_WGT, beta, uint32_t, r)
{
    if (bqd_family() == BQD_FAMILY_X) return CALL(bqd_sm_exp, eps, beta, r);
    const BQD_EXP me = mtbdd_refs_push(CALL(bqd_xpsm_exp_times, eps, w8_exponent(beta)));
    const BQD res = CALL(bqd_xpsm_phase_mul_p, LIMDD_TERMINAL, me);
    mtbdd_refs_pop(1);
    return res;
}

/**
 * MulOffX, the translation family: [k] on supp [U], c [k] off it, for a
 * canonical k = (c' X^t, V): k where c is 1, U has full support, or the shadow
 * puts supp [k] inside supp [U]; c k where U is zero; and otherwise
 * c' X^t . PhaseMulX(V, X^t (1 - iota U), c), whose label is canonical since
 * the multiplication keeps the support. Only the public entry point calls it:
 * no operation of this rule multiplies off a support.
 */
TASK_IMPL_3(BQD, bqd_xpsm_mul_off, BQD, k, EVBDD_WGT, c, LIMDD_TARG, U)
{
    if (limdd_edge_is_zero(k) || c == EVBDD_ONE) return k;
    if (U == 0) return xscale(c, k);
    const label_t lk = label_read(limdd_label(k));
    if (bqd_sm_full(U) || CALL(bqd_xpsm_sub, limdd_target(k), U, lk.t)) return k;
    const uint32_t r = bqd_exp_order(c);
    const BQD_EXP one = mtbdd_refs_push(bqd_exp_const(1, r));
    const BQD_EXP io = mtbdd_refs_push(CALL(bqd_xp_iota, U));
    const BQD_EXP off = mtbdd_refs_push(CALL(bqd_exp_sub, one, io, r));
    const BQD_EXP eps = mtbdd_refs_push(CALL(bqd_exp_translate, off, lk.t));
    const BQD p = CALL(bqd_xpsm_phase_mul_x, limdd_target(k), eps, c, r);
    mtbdd_refs_pop(4);
    return xedge(label_mul(lk, label_read(limdd_label(p))), limdd_target(p));
}

/**
 * The canonical edge of the function whose cofactors at var a are [L] and
 * [L] . (h, H), for L and H nodes of full support and (h, H) a labelled edge
 * with no X: Compose where the least point is 0 and the ratio is known, so
 * that nothing is divided. A constant ratio of 1 skips the level, one of -1
 * skips it with Z_a, and a ratio whose value at 0 has no argument in [0, pi)
 * gets the sign repair. The node is a Q node of full support, the copy
 * rule's, made by bqd_sm_mk.
 */
static BQD
pauli_node(uint32_t a, LIMDD_TARG L, label_t h, LIMDD_TARG H)
{
    const label_t za = { EVBDD_ONE, bit(a), 0 };
    if (H == LIMDD_TERMINAL && h.s == 0) {
        if (h.c == EVBDD_ONE) return unit(L);
        if (h.c == EVBDD_MIN_ONE) return limdd_bundle(label_lim(za), L);
    }
    const bool rep = !bqd_arg_in_upper(h.c);
    if (rep) h.c = neg(h.c);
    const LIMDD_TARG node = bqd_sm_mk(a, false, unit(L), limdd_bundle(label_lim(h), H));
    return rep ? limdd_bundle(label_lim(za), node) : unit(node);
}

/**
 * The diagonal walk of the Pauli family on full support (skip:alg:xdiag),
 * bqd_xp_diagonal's: on full support the two rules are one diagram, every
 * node a Q node, so the walk is the copy rule's path with bqd_sm_mk and this
 * rule's Mono. One path, n steps at most, no memo; *visits counts the nodes
 * of e on it.
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
        const BQD m = limdd_refs_push(bqd_xpsm_monomial(rest, phase, n));
        const BQD r = pauli_node(a, N, label_read(limdd_label(m)), limdd_target(m));
        limdd_refs_pop(1);
        return r;
    }

    assert(bqd_sm_full(N));
    (*visits)++;
    const BQD lo = limdd_node_low(N), hi = limdd_node_high(N);
    if (var < a) {
        const BQD d = limdd_refs_push(pauli_diag(limdd_target(lo), A, phase, n, visits));
        const LIMDD_TARG node = bqd_sm_mk(var, false, unit(limdd_target(d)), hi);
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
bqd_xpsm_diagonal(BQD e, uint64_t A, EVBDD_WGT phase, uint32_t n, uint32_t *visits)
{
    assert(bqd_family() == BQD_FAMILY_PAULI && (limdd_edge_is_zero(e) || bqd_has_full_support(e)));
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

/* --- building ----------------------------------------------------------------------- */

static inline bool
all_zero(const EVBDD_WGT *g, uint64_t len)
{
    for (uint64_t y = 0; y < len; y++) if (g[y] != EVBDD_ZERO) return false;
    return true;
}

/** The index of the least point of the support, or -1 for the zero vector. */
static int64_t
lexmin(const EVBDD_WGT *g, uint64_t len)
{
    for (uint64_t y = 0; y < len; y++) if (g[y] != EVBDD_ZERO) return (int64_t)y;
    return -1;
}

/**
 * The label that recovers g = (c0 Z^s X^t) . rep read as c Z^s X^t, from
 * vector-index masks, as bqd_build's.
 */
static inline LIMDD_LIM
label_of(EVBDD_WGT c0, uint64_t s, uint64_t t, uint32_t n)
{
    return bqd_lim_make(parity(s & t) ? neg(c0) : c0, s, t, n);
}

/**
 * The node of a representative g at level var, of length 2^(n-var): Build
 * of skip:alg:smxpcons, skip:def:smxp read as code, with vector-index masks.
 * Equal cofactors skip the level; in the Pauli family the sign pattern of the
 * representative has already turned f_1 = -f_0 into two equal ones. The low
 * cofactor of a representative is one itself, so the low edge is the
 * identity on its node. The high cofactor is moved by its least point t, and
 * the tag is S where it then has a point outside supp f_0: the high edge is
 * then the representative's edge of the moved cofactor, and otherwise that of
 * its ratio to f_0, zero off supp f_0; X^t goes on the label either way. The
 * low child is spawned and the high one computed, as in the copy rule's.
 */
TASK_3(LIMDD_TARG, bqd_xpsm_build, const EVBDD_WGT *, g, uint32_t, var, uint32_t, n)
{
    if (var == n) return LIMDD_TERMINAL;                 /* g = [1] */

    const uint64_t h = UINT64_C(1) << (n - var - 1);
    const EVBDD_WGT *f0 = g, *f1 = g + h;
    if (memcmp(f0, f1, h * sizeof(EVBDD_WGT)) == 0) return CALL(bqd_xpsm_build, f0, var + 1, n);

    SPAWN(bqd_xpsm_build, f0, var + 1, n);

    bool s = false;
    LIMDD high = limdd_zero_edge();
    const int64_t m = lexmin(f1, h);
    if (m >= 0) {
        const uint64_t t = (uint64_t)m;
        EVBDD_WGT *raw = malloc(h * sizeof(EVBDD_WGT));
        if (raw == NULL) { fprintf(stderr, "sylvan: out of memory in a BQD build\n"); exit(1); }
        for (uint64_t y = 0; y < h; y++) {
            raw[y] = f1[y ^ t];
            if (raw[y] != EVBDD_ZERO && f0[y] == EVBDD_ZERO) s = true;
        }
        if (!s)
            for (uint64_t y = 0; y < h; y++)
                if (raw[y] != EVBDD_ZERO) raw[y] = wgt_div(raw[y], f0[y]);
        EVBDD_WGT c0; uint64_t zs, zp;
        bqd_representative(raw, h, raw, &c0, &zs, &zp);
        assert(zp == 0 && "0 is the least point once the high cofactor is moved");
        const LIMDD_TARG ht = CALL(bqd_xpsm_build, raw, var + 1, n);
        high = limdd_bundle(label_of(c0, zs, t, n), ht);
        free(raw);
    }

    const LIMDD low = unit(SYNC(bqd_xpsm_build));
    return bqd_sm_mk(var, s, low, high);
}

TASK_IMPL_2(BQD, bqd_xpsm_from_vector, const EVBDD_WGT *, f, uint32_t, n)
{
    const uint64_t len = UINT64_C(1) << n;
    if (all_zero(f, len)) return limdd_zero_edge();
    EVBDD_WGT *rep = malloc(len * sizeof(EVBDD_WGT));
    if (rep == NULL) { fprintf(stderr, "sylvan: out of memory in a BQD build\n"); exit(1); }
    EVBDD_WGT c0; uint64_t s, p;
    bqd_representative(f, len, rep, &c0, &s, &p);
    const LIMDD_TARG root = CALL(bqd_xpsm_build, rep, 0, n);
    free(rep);
    return limdd_bundle(label_of(c0, s, p, n), root);
}

/* --- decoding ----------------------------------------------------------------------- */

static void
fill_zero(EVBDD_WGT *out, uint64_t len)
{
    for (uint64_t y = 0; y < len; y++) out[y] = EVBDD_ZERO;
}

/**
 * The function of p read at level var, into out[0..2^(n-var)), as the copy
 * rule's decoder, with the high half by the tag: f_1(y) = c (-1)^{s.y} times
 * [A](y ^ t) [R](y ^ t) for a Q node and [R](y ^ t) for an S node, the plain
 * product, for the high label c Z^s X^t. The levels an edge skips are the
 * leading bits of the index, over which p's block repeats.
 */
VOID_TASK_4(bqd_xpsm_decode, LIMDD_TARG, p, uint32_t, var, uint32_t, n, EVBDD_WGT *, out)
{
    const uint32_t v = (p == LIMDD_TERMINAL) ? n : limdd_node_var(p);
    assert(var <= v && v <= n);
    const uint64_t len = UINT64_C(1) << (n - v);

    if (v == n) {
        out[0] = EVBDD_ONE;
    } else {
        const uint64_t h = len >> 1;
        const LIMDD low = limdd_node_low(p), high = limdd_node_high(p);
        const bool lz = limdd_edge_is_zero(low), hz = limdd_edge_is_zero(high);

        if (!lz) SPAWN(bqd_xpsm_decode, limdd_target(low), v + 1, n, out);

        EVBDD_WGT *g1 = NULL;
        if (!hz) {
            g1 = malloc(h * sizeof(EVBDD_WGT));
            if (g1 == NULL) { fprintf(stderr, "sylvan: out of memory decoding a BQD\n"); exit(1); }
            CALL(bqd_xpsm_decode, limdd_target(high), v + 1, n, g1);
        }

        if (!lz) SYNC(bqd_xpsm_decode); else fill_zero(out, h);

        if (hz) {
            fill_zero(out + h, h);
        } else {
            EVBDD_WGT c; uint64_t zs, t;
            bqd_lim_masks(limdd_label(high), n, &c, &zs, &t);
            const bool s = bqd_sm_is_s(p);
            for (uint64_t y = 0; y < h; y++) {
                const EVBDD_WGT a = out[y ^ t], b = g1[y ^ t];
                EVBDD_WGT w = EVBDD_ZERO;
                if (b != EVBDD_ZERO && (s || a != EVBDD_ZERO)) {
                    w = s ? b : wgt_mul(a, b);
                    if (parity(zs & y)) w = neg(w);
                    w = wgt_mul(c, w);
                }
                out[h + y] = w;
            }
            free(g1);
        }
    }

    for (uint64_t y = len; y < (UINT64_C(1) << (n - var)); y += len)
        memcpy(out + y, out, len * sizeof(EVBDD_WGT));
}

VOID_TASK_IMPL_3(bqd_xpsm_to_vector, BQD, e, uint32_t, n, EVBDD_WGT *, out)
{
    const uint64_t len = UINT64_C(1) << n;
    if (limdd_edge_is_zero(e)) { fill_zero(out, len); return; }

    EVBDD_WGT *rep = malloc(len * sizeof(EVBDD_WGT));
    if (rep == NULL) { fprintf(stderr, "sylvan: out of memory decoding a BQD\n"); exit(1); }
    CALL(bqd_xpsm_decode, limdd_target(e), 0, n, rep);

    EVBDD_WGT c; uint64_t s, t;
    bqd_lim_masks(limdd_label(e), n, &c, &s, &t);
    for (uint64_t y = 0; y < len; y++) {
        EVBDD_WGT w = rep[y ^ t];
        if (w != EVBDD_ZERO) {
            if (parity(s & y)) w = neg(w);
            w = wgt_mul(c, w);
        }
        out[y] = w;
    }
    free(rep);
}

/* --- one amplitude ------------------------------------------------------------------ */

/**
 * The value of p's function at y. The low half is one path, and so is the
 * high half of an S node; that of a Q node needs both children at the same
 * point, unless the low one is zero there.
 */
static EVBDD_WGT
eval_node(LIMDD_TARG p, uint32_t n, uint64_t y)
{
    if (p == LIMDD_TERMINAL) return EVBDD_ONE;
    const uint64_t h = UINT64_C(1) << (n - limdd_node_var(p) - 1);
    const uint64_t yy = y & (h - 1);
    const LIMDD low = limdd_node_low(p);
    if ((y & h) == 0)
        return limdd_edge_is_zero(low) ? EVBDD_ZERO : eval_node(limdd_target(low), n, yy);

    const LIMDD high = limdd_node_high(p);
    if (limdd_edge_is_zero(high)) return EVBDD_ZERO;
    EVBDD_WGT c; uint64_t s, t;
    bqd_lim_masks(limdd_label(high), n, &c, &s, &t);
    const uint64_t yt = yy ^ t;
    EVBDD_WGT a = EVBDD_ONE;
    if (!bqd_sm_is_s(p)) {
        a = limdd_edge_is_zero(low) ? EVBDD_ZERO : eval_node(limdd_target(low), n, yt);
        if (a == EVBDD_ZERO) return a;
    }
    const EVBDD_WGT b = eval_node(limdd_target(high), n, yt);
    if (b == EVBDD_ZERO) return b;
    EVBDD_WGT w = (a == EVBDD_ONE) ? b : wgt_mul(a, b);
    if (parity(s & yy)) w = neg(w);
    return wgt_mul(c, w);
}

EVBDD_WGT
bqd_xpsm_eval(BQD e, uint32_t n, uint64_t x)
{
    if (limdd_edge_is_zero(e)) return EVBDD_ZERO;
    EVBDD_WGT c; uint64_t s, t;
    bqd_lim_masks(limdd_label(e), n, &c, &s, &t);
    EVBDD_WGT w = eval_node(limdd_target(e), n, x ^ t);
    if (w == EVBDD_ZERO) return w;
    if (parity(s & x)) w = neg(w);
    return wgt_mul(c, w);
}

/* --- states ------------------------------------------------------------------------- */

/**
 * |x> is X^x on |0...0>: the least point of its support is x, and the
 * translate is |0...0>, a Q node at every level with the chain below it on
 * the low side and no high cofactor, with no pivot and so no sign. Every
 * basis state is one chain under a different label, as under the copy rule.
 */
BQD
bqd_xpsm_basis_state(uint64_t x, uint32_t n)
{
    LIMDD_TARG N = LIMDD_TERMINAL;
    label_t l = { EVBDD_ONE, 0, 0 };
    for (uint32_t v = n; v-- > 0; ) {
        N = bqd_sm_mk(v, false, unit(N), limdd_zero_edge());
        if ((x >> (n - 1 - v)) & 1) l.t |= bit(v);           /* qubit v at bit n-1-v */
    }
    return limdd_bundle(label_lim(l), N);
}

/**
 * Mono in the Pauli family, bqd_xp_monomial with this rule's Compose: the
 * constant one at the low side of each variable of A, from the last up, and
 * the monomial of the rest of A on the high side; a phase of -1 at the last
 * variable is a Z on the label and no node. Every node is a Q node of full
 * support.
 */
BQD
bqd_xpsm_monomial(uint64_t A, EVBDD_WGT phase, uint32_t n)
{
    BQD h = limdd_refs_push(constant(phase));
    if (phase != EVBDD_ONE) {
        /* qubit q is bit n-1-q, so the lowest variable is the lowest bit */
        for (uint64_t rest = A; rest != 0; rest &= rest - 1) {
            const uint32_t v = n - 1 - (uint32_t)__builtin_ctzll(rest);
            const BQD next = bqd_xpsm_compose(v, limdd_one_edge(), h);
            limdd_refs_pop(1);
            h = limdd_refs_push(next);
        }
    }
    limdd_refs_pop(1);
    return h;
}
