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

#include "qsylvan_bqd_sm.h"
#include "qsylvan_bqd_gates.h"
#include "qsylvan_gates.h"
#include "qsylvan_limdd_gc.h"

/* --- edges ---------------------------------------------------------------- */

static inline EVBDD_WGT
scalar_of(BQD e)
{
    return limdd_lim_weight(limdd_label(e));
}

/** c . [t], the zero edge for c = 0, with the label bqd_lim_word interns. */
static inline BQD
edge(EVBDD_WGT c, LIMDD_TARG t)
{
    if (c == EVBDD_ZERO) return limdd_zero_edge();
    return limdd_bundle(bqd_lim_word(c, 0, 0), t);
}

static inline BQD
unit(LIMDD_TARG t)
{
    return limdd_bundle(LIMDD_LIM_IDENTITY, t);
}

static inline BQD
scale(EVBDD_WGT c, BQD e)
{
    if (c == EVBDD_ZERO || limdd_edge_is_zero(e)) return limdd_zero_edge();
    if (c == EVBDD_ONE) return e;
    return edge(wgt_mul(c, scalar_of(e)), limdd_target(e));
}

/** The node of an edge, or 0 for the zero edge: 0 is no bucket llmsset hands out. */
static inline LIMDD_TARG
node_or_zero(BQD e)
{
    return limdd_edge_is_zero(e) ? 0 : limdd_target(e);
}

/** Whether [N] read at level var is nested there: a Q node, or below var, or the terminal. */
static inline bool
nested_at(LIMDD_TARG N, uint32_t var)
{
    return limdd_level(N) > var || !bqd_sm_is_s(N);
}

/* --- nodes ---------------------------------------------------------------- */

LIMDD_TARG
bqd_sm_mk(uint32_t var, bool s, BQD lo, BQD hi)
{
    const bool lz = limdd_edge_is_zero(lo), hz = limdd_edge_is_zero(hi);
    assert(s || !lz);                                            /* a Q node has a low cofactor */
    assert(!s || !hz);                                           /* an S node a high one */
    assert(!lz || limdd_label(hi) == LIMDD_LIM_IDENTITY);        /* and is its representative */
    const bool full = !s && !lz && !hz
                   && bqd_sm_full(limdd_target(lo)) && bqd_sm_full(limdd_target(hi));
    return limdd_makenode_flags(var, lo, hi, (s ? BQD_SM_S : 0) | (full ? BQD_SM_FULL : 0));
}

/* --- the support shadow ---------------------------------------------------- */

/*
 * The support of a node is {0} x supp [A] and {1} x supp [R] for either tag
 * (skip:lem:smshadow (a)), so the three walks below read the targets of the
 * two stored edges and nothing else: no scalar, no product. A node below the
 * level a walk is at is both of its cofactors there.
 */

/**
 * Ind: the node of the indicator of supp [N]. Its cofactors are the
 * indicators of the two targets, and its tag N's: a Q node's high support is
 * inside its low one, and an S node's is not, which the indicator keeps. The
 * ratio of a Q indicator is the indicator of its high cofactor, so the level
 * is skipped where the two are one node. The terminal on full support, from
 * the flag. Memoised on N.
 */
TASK_IMPL_1(LIMDD_TARG, bqd_sm_ind, LIMDD_TARG, N)
{
    if (N == 0) return 0;
    if (bqd_sm_full(N)) return LIMDD_TERMINAL;
    uint64_t hit;
    if (cache_get3(CACHE_BQD_SM_IND, N, 0, 0, &hit)) return (LIMDD_TARG)hit;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_CANON);
    const LIMDD_TARG n0 = node_or_zero(limdd_node_low(N)), n1 = node_or_zero(limdd_node_high(N));
    if (n0 != 0) SPAWN(bqd_sm_ind, n0);
    const LIMDD_TARG i1 = (n1 != 0) ? CALL(bqd_sm_ind, n1) : 0;
    const LIMDD_TARG i0 = (n0 != 0) ? SYNC(bqd_sm_ind) : 0;
    const bool s = bqd_sm_is_s(N);
    const LIMDD_TARG r = (!s && i0 == i1)
        ? i0
        : bqd_sm_mk(limdd_node_var(N), s, i0 ? unit(i0) : limdd_zero_edge(),
                    i1 ? unit(i1) : limdd_zero_edge());
    cache_put3(CACHE_BQD_SM_IND, N, 0, 0, r);
    return r;
}

/**
 * Sub: whether supp [X] lies inside supp [Y], at the higher of the two levels
 * on both halves. Free where X is zero, the two are one node, or Y is full,
 * and false where Y is zero or X is full and Y is not; memoised on the pair
 * otherwise. The halves are independent, and one is spawned.
 */
TASK_IMPL_2(int, bqd_sm_sub, LIMDD_TARG, X, LIMDD_TARG, Y)
{
    if (X == 0 || X == Y || bqd_sm_full(Y)) return 1;
    if (Y == 0 || bqd_sm_full(X)) return 0;
    uint64_t hit;
    if (cache_get3(CACHE_BQD_SM_SUBSET, X, Y, 0, &hit)) return (int)hit - 1;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_CANON);
    const uint32_t lx = limdd_node_var(X), ly = limdd_node_var(Y);
    const uint32_t v = lx < ly ? lx : ly;
    const LIMDD_TARG x0 = lx > v ? X : node_or_zero(limdd_node_low(X));
    const LIMDD_TARG x1 = lx > v ? X : node_or_zero(limdd_node_high(X));
    const LIMDD_TARG y0 = ly > v ? Y : node_or_zero(limdd_node_low(Y));
    const LIMDD_TARG y1 = ly > v ? Y : node_or_zero(limdd_node_high(Y));
    SPAWN(bqd_sm_sub, x0, y0);
    const int s1 = CALL(bqd_sm_sub, x1, y1);
    const int s0 = SYNC(bqd_sm_sub);
    const int r = s0 && s1;
    cache_put3(CACHE_BQD_SM_SUBSET, X, Y, 0, (uint64_t)r + 1);   /* 0 would read as a miss */
    return r;
}

/**
 * Dep: whether a node of variable q is reachable from N, that is, whether
 * some part of [N] that the diagram stores depends on x_q. The low child
 * first, and the high one only where it has none. Memoised on (N, q).
 */
TASK_IMPL_2(int, bqd_sm_dep, LIMDD_TARG, N, uint32_t, q)
{
    if (N == 0 || N == LIMDD_TERMINAL) return 0;
    const uint32_t var = limdd_node_var(N);
    if (var > q) return 0;
    if (var == q) return 1;
    uint64_t hit;
    if (cache_get3(CACHE_BQD_SM_DEP, N, q, 0, &hit)) return (int)hit - 1;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_CANON);
    const int r = CALL(bqd_sm_dep, node_or_zero(limdd_node_low(N)), q)
               || CALL(bqd_sm_dep, node_or_zero(limdd_node_high(N)), q);
    cache_put3(CACHE_BQD_SM_DEP, N, q, 0, (uint64_t)r + 1);
    return r;
}

/* --- constructors ----------------------------------------------------------- */

/**
 * mk_nested of skip:alg:smcons. With lo = a [U], the high cofactor is
 * a [U] . rho, and rho is zero off supp U, so the node of f / a is the Q node
 * (U, rho) as it stands: no quotient and no subset test. The level is skipped
 * where the two cofactors agree, which is where rho is the identity on the
 * node of supp U's indicator. The indicator is asked for only when the label
 * is the identity and rho's node is not above U, since an indicator's node is
 * never above the node it is of. A zero rho is the Q node with no high
 * cofactor.
 */
TASK_IMPL_3(BQD, bqd_sm_mk_nested, uint32_t, var, BQD, lo, BQD, rho)
{
    if (limdd_edge_is_zero(lo)) return lo;
    const LIMDD_TARG U = limdd_target(lo);
    /* the zero edge's label is not the identity */
    if (limdd_label(rho) == LIMDD_LIM_IDENTITY && limdd_level(limdd_target(rho)) >= limdd_level(U)
        && limdd_target(rho) == CALL(bqd_sm_ind, U)) return lo;
    return edge(scalar_of(lo), bqd_sm_mk(var, false, unit(U), rho));
}

/**
 * Compose of skip:alg:smcons. Equal cofactors are a function that does not
 * depend on x_var, whose canonical edge is theirs read at var: the level is
 * skipped, in one comparison of canonical edges. A zero low cofactor makes an
 * S node whose high edge is hi's node, with hi's scale on top. Otherwise the
 * node is f / a for lo = a [U], with h = hi / a: a Q node with the ratio of h
 * to [U] where supp h lies inside supp U, which the support shadow decides
 * with no arithmetic, and an S node with h itself where it does not, which
 * costs no quotient at all. The Q node's ratio is never the indicator of
 * supp U, since h would then be [U] = lo / a and lo == hi.
 */
TASK_IMPL_3(BQD, bqd_sm_compose, uint32_t, var, BQD, lo, BQD, hi)
{
    if (lo == hi) return lo;
    if (limdd_edge_is_zero(lo))
        return edge(scalar_of(hi), bqd_sm_mk(var, true, lo, unit(limdd_target(hi))));
    const EVBDD_WGT a = scalar_of(lo);
    const LIMDD_TARG U = limdd_target(lo);
    const BQD h = scale(wgt_div(EVBDD_ONE, a), hi);
    if (limdd_edge_is_zero(h) || CALL(bqd_sm_sub, limdd_target(h), U)) {
        const BQD r = CALL(bqd_sm_apply, BQD_SM_RATIO, unit(U), h);
        return edge(a, bqd_sm_mk(var, false, unit(U), r));
    }
    return edge(a, bqd_sm_mk(var, true, unit(U), h));
}

BQD
bqd_sm_rebuild_s(uint32_t var, BQD lo, BQD hi)
{
    assert(!limdd_edge_is_zero(hi));
    if (limdd_edge_is_zero(lo))
        return edge(scalar_of(hi), bqd_sm_mk(var, true, lo, unit(limdd_target(hi))));
    const EVBDD_WGT a = scalar_of(lo);
    return edge(a, bqd_sm_mk(var, true, unit(limdd_target(lo)), scale(wgt_div(EVBDD_ONE, a), hi)));
}

/* --- cofactors -------------------------------------------------------------- */

VOID_TASK_IMPL_4(bqd_sm_nest, LIMDD_TARG, N, uint32_t, var, BQD *, lo, BQD *, ratio)
{
    if (limdd_level(N) > var) {
        *lo = unit(N);
        *ratio = unit(CALL(bqd_sm_ind, N));
        return;
    }
    assert(!bqd_sm_is_s(N));
    *lo = limdd_node_low(N);
    *ratio = limdd_node_high(N);
}

/**
 * Cof1: the high cofactor of node N. An S node stores it, and so does a Q node
 * with no high cofactor, where it is zero; otherwise it is the product of the
 * low edge and the ratio, memoised on N.
 */
TASK_IMPL_1(BQD, bqd_sm_cof1, LIMDD_TARG, N)
{
    const BQD hi = limdd_node_high(N);
    if (bqd_sm_is_s(N) || limdd_edge_is_zero(hi)) return hi;
    uint64_t hit;
    if (cache_get3(CACHE_BQD_SM_COF1, N, 0, 0, &hit)) return (BQD)hit;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_COF1);
    const BQD r = CALL(bqd_sm_apply, BQD_SM_MUL, limdd_node_low(N), hi);
    cache_put3(CACHE_BQD_SM_COF1, N, 0, 0, (uint64_t)r);
    return r;
}

/*
 * An edge that skips var denotes a function that does not depend on x_var,
 * so it is both of its cofactors there (skip:lem:virtual).
 */
TASK_IMPL_3(BQD, bqd_sm_cof, BQD, e, uint32_t, var, int, b)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG N = limdd_target(e);
    assert(var <= limdd_level(N) && "an edge is read at or above its node");
    if (limdd_level(N) > var) return e;
    return scale(scalar_of(e), b ? CALL(bqd_sm_cof1, N) : limdd_node_low(N));
}

/* --- pointwise operations ---------------------------------------------------- */

/**
 * Apply of skip:alg:smapply, op(f, g) for two edges read at one level. The
 * scalars come out first, so that the memo is on the two nodes and what is
 * left of the scalars:
 *
 *     mul     c.F  d.G  =  cd . (F G)              key (min, max of the nodes)
 *     ratio   c.F  d.G  =  (d/c) . (G / F)         key (F, G)
 *     add     c.F  d.G  =  c . (F + (d/c) G)       key (F, G, d/c)
 *
 * The recursion is at the higher of the two nodes' levels, var. Where both are
 * nested there, the cofactors of a product are the products of the cofactors,
 * and with f_1 = f_0 rho and g_1 = g_0 sigma so is the ratio: (fg)_1 =
 * (f_0 g_0)(rho sigma), and rho sigma is zero off supp f_0 g_0. So the node of
 * a product or a ratio of two nested operands is the Q node of the two results
 * on the stored edges, with no high cofactor and no quotient: the paper's
 * prop:prodscalar, extended from full support to nested nodes. A sum of two
 * nested operands with one ratio rho is (f_0 + g_0) rho, whose ratio is rho
 * restricted to the support of the new low cofactor. Every other pair takes
 * the cofactors and Compose: an S node on either side, or a sum with two
 * ratios. An operand below var takes part as its virtual node, which is
 * nested, and the terminal, the constant one, is the unit of the product and
 * of the ratio's denominator. Both recursions overlap their two calls.
 */
TASK_IMPL_3(BQD, bqd_sm_apply, int, op, BQD, f, BQD, g)
{
    const bool fz = limdd_edge_is_zero(f), gz = limdd_edge_is_zero(g);
    const LIMDD_TARG Nf = limdd_target(f), Ng = limdd_target(g);
    EVBDD_WGT outer, k = EVBDD_ONE;
    uint64_t id, k0 = Nf, k1 = Ng, k2 = 0;
    switch (op) {
    case BQD_SM_MUL:
        if (fz || gz) return limdd_zero_edge();
        if (Nf == LIMDD_TERMINAL) return scale(scalar_of(f), g);
        if (Ng == LIMDD_TERMINAL) return scale(scalar_of(g), f);
        outer = wgt_mul(scalar_of(f), scalar_of(g));
        if (Nf > Ng) { k0 = Ng; k1 = Nf; }               /* it commutes */
        id = CACHE_BQD_SM_MUL;
        break;
    case BQD_SM_RATIO:
        if (fz || gz) return limdd_zero_edge();
        if (Nf == LIMDD_TERMINAL) return scale(wgt_div(EVBDD_ONE, scalar_of(f)), g);
        if (Nf == Ng) return edge(wgt_div(scalar_of(g), scalar_of(f)), CALL(bqd_sm_ind, Nf));
        outer = wgt_div(scalar_of(g), scalar_of(f));
        id = CACHE_BQD_SM_RATIO;
        break;
    default:
        if (fz) return g;
        if (gz) return f;
        if (Nf == Ng) return edge(wgt_add(scalar_of(f), scalar_of(g)), Nf);
        k = wgt_div(scalar_of(g), scalar_of(f));
        outer = scalar_of(f);
        k2 = k;
        id = CACHE_BQD_SM_ADD;
        break;
    }

    uint64_t hit;
    if (cache_get3(id, k0, k1, k2, &hit)) return scale(outer, (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_APPLY);

    const BQD F = unit(Nf), G = edge(k, Ng);
    const uint32_t lf = limdd_level(Nf), lg = limdd_level(Ng);
    const uint32_t var = lf < lg ? lf : lg;
    BQD r;
    bool done = false;
    if (nested_at(Nf, var) && nested_at(Ng, var)) {
        BQD a, rho, b, sigma;
        CALL(bqd_sm_nest, Nf, var, &a, &rho);
        CALL(bqd_sm_nest, Ng, var, &b, &sigma);
        if (op != BQD_SM_ADD) {
            limdd_refs_spawn(SPAWN(bqd_sm_apply, op, a, b));
            const BQD hi = limdd_refs_push(CALL(bqd_sm_apply, op, rho, sigma));
            const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_sm_apply)));
            r = CALL(bqd_sm_mk_nested, var, lo, hi);
            limdd_refs_pop(2);
            done = true;
        } else if (rho == sigma) {
            /* the ratio's scale is not G's, so only the low edge takes k */
            const BQD L = limdd_refs_push(CALL(bqd_sm_apply, BQD_SM_ADD, a, scale(k, b)));
            if (limdd_edge_is_zero(L)) {
                r = L;
            } else {
                const BQD rr = limdd_refs_push(CALL(bqd_sm_restr, rho, L));
                r = CALL(bqd_sm_mk_nested, var, L, rr);
                limdd_refs_pop(1);
            }
            limdd_refs_pop(1);
            done = true;
        }
    }
    if (!done) {
        const BQD f0 = CALL(bqd_sm_cof, F, var, 0), g0 = CALL(bqd_sm_cof, G, var, 0);
        const BQD f1 = CALL(bqd_sm_cof, F, var, 1), g1 = CALL(bqd_sm_cof, G, var, 1);
        limdd_refs_spawn(SPAWN(bqd_sm_apply, op, f0, g0));
        const BQD hi = limdd_refs_push(CALL(bqd_sm_apply, op, f1, g1));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_sm_apply)));
        r = CALL(bqd_sm_compose, var, lo, hi);
        limdd_refs_pop(2);
    }
    cache_put3(id, k0, k1, k2, (uint64_t)r);
    return scale(outer, r);
}

/**
 * Restr: rho on the support of Y, zero off it. Nothing to do where supp rho
 * lies inside supp Y already, which the shadow decides; otherwise the product
 * with Y's indicator.
 */
TASK_IMPL_2(BQD, bqd_sm_restr, BQD, rho, BQD, Y)
{
    assert(!limdd_edge_is_zero(Y));
    if (limdd_edge_is_zero(rho) || CALL(bqd_sm_sub, limdd_target(rho), limdd_target(Y))) return rho;
    const LIMDD_TARG I = CALL(bqd_sm_ind, limdd_target(Y));
    return CALL(bqd_sm_apply, BQD_SM_MUL, rho, unit(I));
}

/* --- selections ----------------------------------------------------------------- */

/**
 * Rebuild of skip:alg:smsel: the node of tag S (s) or Q from the images lo and
 * hi of a node's two stored edges under a map that keeps what the tag says,
 * the Shannon rebuild for S and mk_nested for Q.
 */
TASK_4(BQD, bqd_sm_rebuild, int, s, uint32_t, var, BQD, lo, BQD, hi)
{
    if (s) return bqd_sm_rebuild_s(var, lo, hi);
    return CALL(bqd_sm_mk_nested, var, lo, hi);
}

/*
 * Restrict and Project of skip:alg:smsel, both linear and so memoised on the
 * node, q and b, with the scale outside. An edge that skips q is its own
 * restriction, and its projection is the node at q with the edge on one side
 * and zero on the other; both before the lookup. At q the restriction is the
 * chosen cofactor, which skips q as an edge read at q, and the projection
 * composes it with zero. Above q each value of the result comes from one
 * value of the operand at the same x_var, so a Q node maps to mk_nested of the
 * images of its low edge and its ratio, whose image is still zero off the
 * image of the low support (skip:lem:smpair (d)). An S node maps to Compose of
 * the images of its two cofactors, which may come out nested, or equal.
 */
TASK_IMPL_3(BQD, bqd_sm_restrict, BQD, e, uint32_t, q, int, b)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG N = limdd_target(e);
    const uint32_t var = limdd_level(N);
    if (var > q) return e;
    const uint64_t key = ((uint64_t)q << 1) | (uint64_t)(b & 1);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_SM_RESTRICT, N, key, 0, &hit)) return scale(scalar_of(e), (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    BQD r;
    if (var == q) {
        r = b ? CALL(bqd_sm_cof1, N) : limdd_node_low(N);
    } else {
        limdd_refs_spawn(SPAWN(bqd_sm_restrict, limdd_node_high(N), q, b));
        const BQD lo = limdd_refs_push(CALL(bqd_sm_restrict, limdd_node_low(N), q, b));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_sm_restrict)));
        r = bqd_sm_is_s(N) ? CALL(bqd_sm_compose, var, lo, hi) : CALL(bqd_sm_mk_nested, var, lo, hi);
        limdd_refs_pop(2);
    }
    cache_put3(CACHE_BQD_SM_RESTRICT, N, key, 0, (uint64_t)r);
    return scale(scalar_of(e), r);
}

TASK_IMPL_3(BQD, bqd_sm_project, BQD, e, uint32_t, q, int, b)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG N = limdd_target(e);
    const uint32_t var = limdd_level(N);
    if (var > q) return b ? CALL(bqd_sm_compose, q, limdd_zero_edge(), e)
                          : CALL(bqd_sm_compose, q, e, limdd_zero_edge());
    const uint64_t key = ((uint64_t)q << 1) | (uint64_t)(b & 1);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_SM_PROJECT, N, key, 0, &hit)) return scale(scalar_of(e), (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    BQD r;
    if (var == q && b) {
        const BQD f1 = limdd_refs_push(CALL(bqd_sm_cof1, N));
        r = CALL(bqd_sm_compose, q, limdd_zero_edge(), f1);
        limdd_refs_pop(1);
    } else if (var == q) {
        r = CALL(bqd_sm_compose, q, limdd_node_low(N), limdd_zero_edge());
    } else {
        limdd_refs_spawn(SPAWN(bqd_sm_project, limdd_node_high(N), q, b));
        const BQD lo = limdd_refs_push(CALL(bqd_sm_project, limdd_node_low(N), q, b));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_sm_project)));
        r = bqd_sm_is_s(N) ? CALL(bqd_sm_compose, var, lo, hi) : CALL(bqd_sm_mk_nested, var, lo, hi);
        limdd_refs_pop(2);
    }
    cache_put3(CACHE_BQD_SM_PROJECT, N, key, 0, (uint64_t)r);
    return scale(scalar_of(e), r);
}

/**
 * XQ of skip:alg:smsel, X on qubit q, memoised on the node and q. Above q a
 * bijection of the points below, which keeps a nested pair nested and a
 * misaligned one misaligned, so the node keeps its tag and is rebuilt from
 * the images of its stored edges; at q the two cofactors trade places, which
 * takes the high one and Compose. An edge that skips q is its own image.
 */
TASK_IMPL_2(BQD, bqd_sm_xq, BQD, e, uint32_t, q)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG N = limdd_target(e);
    const uint32_t var = limdd_level(N);
    if (var > q) return e;
    uint64_t hit;
    if (cache_get3(CACHE_BQD_SM_XQ, N, q, 0, &hit)) return scale(scalar_of(e), (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    BQD r;
    if (var == q) {
        const BQD f1 = limdd_refs_push(CALL(bqd_sm_cof1, N));
        r = CALL(bqd_sm_compose, q, f1, limdd_node_low(N));
        limdd_refs_pop(1);
    } else {
        limdd_refs_spawn(SPAWN(bqd_sm_xq, limdd_node_high(N), q));
        const BQD lo = limdd_refs_push(CALL(bqd_sm_xq, limdd_node_low(N), q));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_sm_xq)));
        r = CALL(bqd_sm_rebuild, bqd_sm_is_s(N), var, lo, hi);
        limdd_refs_pop(2);
    }
    cache_put3(CACHE_BQD_SM_XQ, N, q, 0, (uint64_t)r);
    return scale(scalar_of(e), r);
}

/**
 * Pair of skip:alg:smsel, [x_b = 0] A|_{x_b = s} + [x_b = 1] B|_{x_b = t},
 * linear in the pair, so the scale g of the first nonzero operand comes out
 * and the memo is on the two nodes, B's scale over g, s, t and b, as the copy
 * rule's. Below b both operands are their own restrictions and the result is
 * their Compose at b; at b it composes the two chosen cofactors. Above b it
 * is a selection at the higher of the two levels, l: where both operands are
 * nested there, or zero, whose pair is (0, 0), the pair of the low cofactors
 * and the pair of the ratios are the result's low cofactor and ratio, B's
 * scale on its low cofactor alone (skip:lem:smpair (d)); otherwise the pairs
 * of the cofactors and Compose.
 */
TASK_IMPL_5(BQD, bqd_sm_pair, BQD, A, int, s, BQD, B, int, t, uint32_t, b)
{
    const bool az = limdd_edge_is_zero(A), bz = limdd_edge_is_zero(B);
    if (az && bz) return A;
    const EVBDD_WGT g = az ? scalar_of(B) : scalar_of(A);
    const BQD A1 = az ? A : unit(limdd_target(A));
    const BQD B1 = bz ? B : az ? unit(limdd_target(B)) : scale(wgt_div(EVBDD_ONE, g), B);
    const uint64_t kappa = bz ? 0 : (uint64_t)scalar_of(B1);
    const uint64_t key = ((uint64_t)node_or_zero(B1) << 16) | ((uint64_t)b << 2)
                       | ((uint64_t)(s & 1) << 1) | (uint64_t)(t & 1);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_SM_PAIR, node_or_zero(A1), key, kappa, &hit)) return scale(g, (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    const LIMDD_TARG NA = limdd_target(A1), NB = limdd_target(B1);   /* the terminal for zero */
    const uint32_t la = limdd_level(NA), lb = limdd_level(NB);
    const uint32_t l = la < lb ? la : lb;
    BQD r;
    if (l > b) {
        r = CALL(bqd_sm_compose, b, A1, B1);
    } else if (l == b) {
        const BQD a0 = limdd_refs_push(CALL(bqd_sm_cof, A1, b, s));
        const BQD b0 = limdd_refs_push(CALL(bqd_sm_cof, B1, b, t));
        r = CALL(bqd_sm_compose, b, a0, b0);
        limdd_refs_pop(2);
    } else if ((az || nested_at(NA, l)) && (bz || nested_at(NB, l))) {
        BQD u = limdd_zero_edge(), rho = u, w = u, sigma = u;
        if (!az) CALL(bqd_sm_nest, NA, l, &u, &rho);
        limdd_refs_push(rho);
        if (!bz) {
            CALL(bqd_sm_nest, NB, l, &w, &sigma);
            w = scale(scalar_of(B1), w);                 /* the ratio carries no scale */
        }
        limdd_refs_push(w);
        limdd_refs_push(sigma);
        limdd_refs_spawn(SPAWN(bqd_sm_pair, rho, s, sigma, t, b));
        const BQD lo = limdd_refs_push(CALL(bqd_sm_pair, u, s, w, t, b));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_sm_pair)));
        r = CALL(bqd_sm_mk_nested, l, lo, hi);
        limdd_refs_pop(5);
    } else {
        const BQD a0 = limdd_refs_push(CALL(bqd_sm_cof, A1, l, 0));
        const BQD a1 = limdd_refs_push(CALL(bqd_sm_cof, A1, l, 1));
        const BQD b0 = limdd_refs_push(CALL(bqd_sm_cof, B1, l, 0));
        const BQD b1 = limdd_refs_push(CALL(bqd_sm_cof, B1, l, 1));
        limdd_refs_spawn(SPAWN(bqd_sm_pair, a1, s, b1, t, b));
        const BQD lo = limdd_refs_push(CALL(bqd_sm_pair, a0, s, b0, t, b));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_sm_pair)));
        r = CALL(bqd_sm_compose, l, lo, hi);
        limdd_refs_pop(6);
    }
    cache_put3(CACHE_BQD_SM_PAIR, node_or_zero(A1), key, kappa, (uint64_t)r);
    return scale(g, r);
}

/**
 * Perm of skip:alg:smsel, the permutation `kind` of the qubits qa < qb,
 * (pi f)(x) = f(pi x), memoised on the node, the kind and the two qubits.
 * Above qa a bijection of the points below, so the node keeps its tag and is
 * rebuilt from the images of its stored edges; at qa the new cofactors are
 * Pairs at qb, or for the CX with the control above f_0 and X_qb f_1; an edge
 * that skips qa is its own image under the CX with the control below, and has
 * the cofactors (f, X_qb f) or the two restrictions at qb otherwise.
 */
TASK_IMPL_4(BQD, bqd_sm_perm, BQD, e, uint32_t, kind, uint32_t, qa, uint32_t, qb)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG N = limdd_target(e);
    const uint32_t var = limdd_level(N);
    if (var > qa && kind == BQD_PERM_CX_UP) return e;
    const uint64_t key = (uint64_t)kind | ((uint64_t)qa << 8) | ((uint64_t)qb << 16);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_SM_PERM, N, key, 0, &hit)) return scale(scalar_of(e), (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    const BQD E = unit(N);
    BQD r;
    if (var > qa) {
        BQD g0, g1;
        if (kind == BQD_PERM_CX_DOWN) {
            g0 = limdd_refs_push(E);
            g1 = limdd_refs_push(CALL(bqd_sm_xq, E, qb));
        } else {
            limdd_refs_spawn(SPAWN(bqd_sm_restrict, E, qb, 1));
            g0 = limdd_refs_push(CALL(bqd_sm_restrict, E, qb, 0));
            g1 = limdd_refs_push(limdd_refs_sync(SYNC(bqd_sm_restrict)));
        }
        r = CALL(bqd_sm_compose, qa, g0, g1);
        limdd_refs_pop(2);
    } else if (var < qa) {
        limdd_refs_spawn(SPAWN(bqd_sm_perm, limdd_node_high(N), kind, qa, qb));
        const BQD lo = limdd_refs_push(CALL(bqd_sm_perm, limdd_node_low(N), kind, qa, qb));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_sm_perm)));
        r = CALL(bqd_sm_rebuild, bqd_sm_is_s(N), var, lo, hi);
        limdd_refs_pop(2);
    } else {
        const BQD f0 = limdd_node_low(N);
        const BQD f1 = limdd_refs_push(CALL(bqd_sm_cof1, N));
        BQD g0, g1;
        if (kind == BQD_PERM_CX_DOWN) {
            g0 = limdd_refs_push(f0);
            g1 = limdd_refs_push(CALL(bqd_sm_xq, f1, qb));
        } else {
            /* SW: (Pair(f0,0; f1,0), Pair(f0,1; f1,1)); CX up: (Pair(f0,0; f1,1), Pair(f1,0; f0,1)) */
            const bool sw = kind == BQD_PERM_SWAP;
            limdd_refs_spawn(SPAWN(bqd_sm_pair, sw ? f0 : f1, sw ? 1 : 0, sw ? f1 : f0, 1, qb));
            g0 = limdd_refs_push(CALL(bqd_sm_pair, f0, 0, f1, sw ? 0 : 1, qb));
            g1 = limdd_refs_push(limdd_refs_sync(SYNC(bqd_sm_pair)));
        }
        r = CALL(bqd_sm_compose, qa, g0, g1);
        limdd_refs_pop(3);
    }
    cache_put3(CACHE_BQD_SM_PERM, N, key, 0, (uint64_t)r);
    return scale(scalar_of(e), r);
}

/* --- gates -------------------------------------------------------------------- */

#define SIDE_NONE 3      /* a zero edge, which constrains nothing */

/**
 * Side: the support of a node is {0} x supp [A] and {1} x supp [R]
 * (skip:lem:smshadow), so it lies on one side of x_q exactly where both
 * targets that are not zero lie on that side, and at q itself where one of
 * the two edges is zero: a zero low edge is an S node with f_0 = 0, and a
 * zero high edge a Q node with f_1 = 0. An edge that skips q, and every full
 * node, has points on both sides. The low child first, and the high one only
 * where the low one has not decided it; memoised on (N, q).
 */
TASK_IMPL_2(int, bqd_sm_side, LIMDD_TARG, N, uint32_t, q)
{
    if (N == LIMDD_TERMINAL || limdd_node_var(N) > q || bqd_sm_full(N)) return 2;
    const BQD low = limdd_node_low(N), high = limdd_node_high(N);
    const bool lz = limdd_edge_is_zero(low), hz = limdd_edge_is_zero(high);
    if (limdd_node_var(N) == q) return lz ? 1 : (hz ? 0 : 2);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_SM_SIDE, N, q, 0, &hit)) return (int)hit - 1;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_CANON);
    const int s0 = lz ? SIDE_NONE : CALL(bqd_sm_side, limdd_target(low), q);
    int r = 2;
    if (s0 != 2) {
        const int s1 = hz ? SIDE_NONE : CALL(bqd_sm_side, limdd_target(high), q);
        r = (s0 == SIDE_NONE) ? s1 : ((s1 == SIDE_NONE || s1 == s0) ? s0 : 2);
    }
    cache_put3(CACHE_BQD_SM_SIDE, N, q, 0, (uint64_t)r + 1);    /* 0 would read as a miss */
    return r;
}

/**
 * Gate of skip:alg:smgate: U = gates[gid] on qubit q, conditioned on the
 * controls in cmask that are still pending, all above q. Linear, so memoised
 * on the node, gid and cmask | 1 << q, the scale outside, which is the copy
 * rule's key: a control is dropped once passed, so the key names the
 * controls to come and holds no level. The work is at p, the highest
 * pending control, or q when none is left.
 *
 * An edge that skips p is both of its cofactors there: a skipped control
 * gives (e, gate(e)), and a skipped target the row sums of U times e. Above p
 * there are two shortcuts before the recursion on both cofactors, which
 * share the memo entry with it (skip:lem:smpass). Where no control is pending
 * and supp f lies in {x_q = s}, Uf is u_0s g at x_q = 0 and u_1s g at x_q = 1
 * for the restriction g, a phase on x_q, whenever u_0s u_1s is not zero. And
 * where a Q node's ratio rho does not depend on x_q, which holds when no node
 * of q is reachable from it, G(f_0 rho) = (G f_0) rho, so the gate acts on
 * the low cofactor alone and the ratio stays, restricted to the new low
 * support, with no high cofactor and no sum. At a control the low cofactor is
 * kept and the gate passes, without that control, to the high one, and at q
 * the new cofactors are the rows of U applied to the old ones.
 */
TASK_IMPL_4(BQD, bqd_sm_gate, BQD, e, uint32_t, gid, uint64_t, cmask, uint32_t, q)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG N = limdd_target(e);
    const uint32_t var = limdd_level(N);
    const uint64_t key = cmask | (UINT64_C(1) << q);
    const uint32_t p = (uint32_t)__builtin_ctzll(key);

    uint64_t hit;
    if (cache_get3(CACHE_BQD_SM_GATE, N, gid, key, &hit)) return scale(scalar_of(e), (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    const EVBDD_WGT u00 = gates[gid][0], u01 = gates[gid][1];
    const EVBDD_WGT u10 = gates[gid][2], u11 = gates[gid][3];
    BQD r;
    if (var > p) {
        const BQD E = unit(N);
        if (p < q) {
            const BQD hi = limdd_refs_push(CALL(bqd_sm_gate, E, gid, cmask & ~(UINT64_C(1) << p), q));
            r = CALL(bqd_sm_compose, p, E, hi);
            limdd_refs_pop(1);
        } else {
            r = CALL(bqd_sm_compose, q, scale(wgt_add(u00, u01), E), scale(wgt_add(u10, u11), E));
        }
    } else if (var < p) {
        const int s = (p == q) ? CALL(bqd_sm_side, N, q) : 2;
        if (s != 2 && gates[gid][s] != EVBDD_ZERO && gates[gid][2 + s] != EVBDD_ZERO) {
            const EVBDD_WGT u0s = gates[gid][s], u1s = gates[gid][2 + s];
            const BQD g = limdd_refs_push(CALL(bqd_sm_restrict, unit(N), q, s));
            r = scale(u0s, CALL(bqd_sm_diag, g, UINT64_C(1) << q, wgt_div(u1s, u0s)));
            limdd_refs_pop(1);
        } else if (!bqd_sm_is_s(N) && !CALL(bqd_sm_dep, node_or_zero(limdd_node_high(N)), q)) {
            const BQD Y = limdd_refs_push(CALL(bqd_sm_gate, limdd_node_low(N), gid, cmask, q));
            if (limdd_edge_is_zero(Y)) {
                r = Y;
            } else {
                const BQD rho = limdd_refs_push(CALL(bqd_sm_restr, limdd_node_high(N), Y));
                r = CALL(bqd_sm_mk_nested, var, Y, rho);
                limdd_refs_pop(1);
            }
            limdd_refs_pop(1);
        } else {
            const BQD f1 = limdd_refs_push(CALL(bqd_sm_cof1, N));
            limdd_refs_spawn(SPAWN(bqd_sm_gate, limdd_node_low(N), gid, cmask, q));
            const BQD hi = limdd_refs_push(CALL(bqd_sm_gate, f1, gid, cmask, q));
            const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_sm_gate)));
            r = CALL(bqd_sm_compose, var, lo, hi);
            limdd_refs_pop(3);
        }
    } else if (var < q) {
        const BQD f1 = limdd_refs_push(CALL(bqd_sm_cof1, N));
        const BQD hi = limdd_refs_push(CALL(bqd_sm_gate, f1, gid, cmask & ~(UINT64_C(1) << var), q));
        r = CALL(bqd_sm_compose, var, limdd_node_low(N), hi);
        limdd_refs_pop(2);
    } else {
        const BQD f0 = limdd_node_low(N);
        const BQD f1 = limdd_refs_push(CALL(bqd_sm_cof1, N));
        limdd_refs_spawn(SPAWN(bqd_sm_apply, BQD_SM_ADD, scale(u10, f0), scale(u11, f1)));
        const BQD lo = limdd_refs_push(CALL(bqd_sm_apply, BQD_SM_ADD, scale(u00, f0), scale(u01, f1)));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_sm_apply)));
        r = CALL(bqd_sm_compose, q, lo, hi);
        limdd_refs_pop(3);
    }
    cache_put3(CACHE_BQD_SM_GATE, N, gid, key, (uint64_t)r);
    return scale(scalar_of(e), r);
}

/* --- diagonal gates and phases ------------------------------------------------- */

/*
 * A diagonal gate multiplies by a function with no zero, which keeps every
 * support, so a node keeps its tag (skip:lem:smpair (e)): above the gate's
 * variables the ratio of a Q node is kept and the low cofactor multiplied,
 * and an S node has both cofactors multiplied. No support is multiplied off,
 * which is what the copy rule's PhaseMul needed Sel and the least point for.
 */

/**
 * Mono of skip:alg:smdiag, phase^{prod of x_v, v in A}: from the bottom up a
 * Q node at each variable of A, its low cofactor the constant one, since the
 * monomial is 0 there, and its ratio the monomial of the variables of A below
 * it, until the last, where it is the phase. Every node has full support.
 */
BQD
bqd_sm_monomial(uint64_t A, EVBDD_WGT phase)
{
    BQD h = edge(phase, LIMDD_TERMINAL);
    if (phase == EVBDD_ONE) return h;
    for (uint64_t rest = A; rest != 0; ) {
        const uint32_t v = 63 - (uint32_t)__builtin_clzll(rest);
        rest &= ~(UINT64_C(1) << v);
        h = unit(bqd_sm_mk(v, false, limdd_one_edge(), h));
    }
    return h;
}

/**
 * Diag of skip:alg:smdiag on any support, memoised on the node, A and the
 * phase, with a = the top variable of A. Where e skips a its cofactors there
 * are f and f phase^{x_rest}, so the node at a is the Q node of f and the
 * indicator of its support times the rest of the monomial. At a, the low
 * cofactor is kept and the high one, or the ratio, multiplied by the rest.
 * Above a, a Q node keeps its ratio, since the phase does not depend on x_var,
 * and an S node has both cofactors multiplied.
 */
TASK_IMPL_3(BQD, bqd_sm_diag, BQD, e, uint64_t, A, EVBDD_WGT, phase)
{
    if (limdd_edge_is_zero(e) || phase == EVBDD_ONE) return e;
    if (A == 0) return scale(phase, e);
    const LIMDD_TARG N = limdd_target(e);
    if (N == LIMDD_TERMINAL) return scale(scalar_of(e), bqd_sm_monomial(A, phase));
    uint64_t hit;
    if (cache_get3(CACHE_BQD_SM_DIAG, N, A, phase, &hit)) return scale(scalar_of(e), (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PHASEMUL);

    const uint32_t a = (uint32_t)__builtin_ctzll(A), var = limdd_level(N);
    const uint64_t rest = A & (A - 1);
    BQD r;
    if (var > a) {
        const BQD I = limdd_refs_push(unit(CALL(bqd_sm_ind, N)));
        const BQD hi = limdd_refs_push(rest ? CALL(bqd_sm_diag, I, rest, phase) : scale(phase, I));
        r = CALL(bqd_sm_mk_nested, a, unit(N), hi);
        limdd_refs_pop(2);
    } else if (var == a) {
        const BQD high = limdd_node_high(N);
        const BQD hi = limdd_refs_push(rest ? CALL(bqd_sm_diag, high, rest, phase) : scale(phase, high));
        r = CALL(bqd_sm_rebuild, bqd_sm_is_s(N), var, limdd_node_low(N), hi);
        limdd_refs_pop(1);
    } else if (!bqd_sm_is_s(N)) {
        const BQD lo = limdd_refs_push(CALL(bqd_sm_diag, limdd_node_low(N), A, phase));
        r = CALL(bqd_sm_mk_nested, var, lo, limdd_node_high(N));
        limdd_refs_pop(1);
    } else {
        limdd_refs_spawn(SPAWN(bqd_sm_diag, limdd_node_high(N), A, phase));
        const BQD lo = limdd_refs_push(CALL(bqd_sm_diag, limdd_node_low(N), A, phase));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_sm_diag)));
        r = bqd_sm_rebuild_s(var, lo, hi);
        limdd_refs_pop(2);
    }
    cache_put3(CACHE_BQD_SM_DIAG, N, A, phase, (uint64_t)r);
    return scale(scalar_of(e), r);
}

/** A Q node on full support, or u itself where the ratio is the constant one. */
static LIMDD_TARG
node_on_full(uint32_t var, LIMDD_TARG u, BQD r)
{
    if (r == limdd_one_edge()) return u;
    return bqd_sm_mk(var, false, unit(u), r);
}

/*
 * The walk of skip:prop:diag, the copy rule's diag_rec (qsylvan_bqd_ops.c)
 * with SM's constructors: on full support every node is a Q node with two
 * full children, the copy rule's node, so the walk is the same path, down the
 * low edge above a and down the high edge at a variable of A, and stops at the
 * last of A or where e skips the next one. Every node of e and of the
 * monomial is 1 at the point 0, and so is the product, so nothing but e's
 * scale comes back up.
 */
static BQD
diag_walk(BQD e, uint64_t A, EVBDD_WGT phase, uint32_t *visits)
{
    const EVBDD_WGT c = scalar_of(e);
    const LIMDD_TARG N = limdd_target(e);
    const uint32_t a = (uint32_t)__builtin_ctzll(A), var = limdd_level(N);
    const uint64_t rest = A & (A - 1);
    if (var > a) return edge(c, node_on_full(a, N, bqd_sm_monomial(rest, phase)));

    assert(bqd_sm_full(N));
    if (visits != NULL) (*visits)++;
    const BQD low = limdd_node_low(N), high = limdd_node_high(N);
    if (var < a) {
        const BQD lo = diag_walk(unit(limdd_target(low)), A, phase, visits);
        assert(scalar_of(lo) == EVBDD_ONE);
        return edge(c, node_on_full(var, limdd_target(lo), high));
    }
    BQD up;
    if (rest == 0) {
        up = edge(wgt_mul(scalar_of(high), phase), limdd_target(high));
    } else {
        const BQD hi = diag_walk(unit(limdd_target(high)), rest, phase, visits);
        assert(scalar_of(hi) == EVBDD_ONE);
        up = limdd_bundle(limdd_label(high), limdd_target(hi));
    }
    return edge(c, node_on_full(var, limdd_target(low), up));
}

BQD
bqd_sm_diag_walk(BQD e, uint64_t A, EVBDD_WGT phase, uint32_t *visits)
{
    if (limdd_edge_is_zero(e)) return e;
    if (A == 0) return scale(phase, e);
    return diag_walk(e, A, phase, visits);
}

/**
 * iota of skip:alg:smdiag, the indicator of supp [N] as a 0/1 exponent, from
 * the shadow: the exponent node of the indicators of the two targets, for
 * either tag. 1 on full support, from the flag. Memoised on N.
 */
TASK_IMPL_1(BQD_EXP, bqd_sm_iota, LIMDD_TARG, N)
{
    if (N == 0) return mtbdd_int64(0);
    if (bqd_sm_full(N)) return mtbdd_int64(1);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_SM_IOTA, N, 0, 0, &hit)) return (BQD_EXP)hit;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_CANON);
    mtbdd_refs_spawn(SPAWN(bqd_sm_iota, node_or_zero(limdd_node_low(N))));
    const BQD_EXP i1 = mtbdd_refs_push(CALL(bqd_sm_iota, node_or_zero(limdd_node_high(N))));
    const BQD_EXP i0 = mtbdd_refs_sync(SYNC(bqd_sm_iota));
    mtbdd_refs_pop(1);
    const BQD_EXP r = bqd_exp_node(limdd_node_var(N), i0, i1);
    cache_put3(CACHE_BQD_SM_IOTA, N, 0, 0, r);
    return r;
}

/**
 * PhaseMul of skip:alg:smdiag: [k] . beta^eps on any support, r the order of
 * beta, memoised on the node, eps and beta, the scale outside. The recursion
 * is at the higher of the node's and the exponent's tops, l. Where the node is
 * nested at l, f_1 beta^{eps_1} = (f_0 beta^{eps_0}) (rho beta^{eps_1 - eps_0}),
 * so the low cofactor and the ratio are multiplied, by eps_0 and by the
 * difference; an S node has its two cofactors multiplied, by eps_0 and eps_1.
 * Neither changes a support, so neither needs the least point or Sel, which
 * the copy rule's needs to multiply off one. The two calls overlap.
 */
TASK_IMPL_4(BQD, bqd_sm_phase_mul, BQD, k, BQD_EXP, eps, EVBDD_WGT, beta, uint32_t, r)
{
    if (limdd_edge_is_zero(k)) return k;
    if (bqd_exp_is_const(eps)) return scale(bqd_exp_power(beta, bqd_exp_value(eps)), k);
    const LIMDD_TARG N = limdd_target(k);
    if (N == LIMDD_TERMINAL) return scale(scalar_of(k), CALL(bqd_sm_exp, eps, beta, r));
    uint64_t hit;
    if (cache_get3(CACHE_BQD_SM_PHASEMUL, N, eps, beta, &hit)) return scale(scalar_of(k), (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PHASEMUL);

    const uint32_t lN = limdd_level(N), le = bqd_exp_var(eps);
    const uint32_t l = lN < le ? lN : le;
    const BQD_EXP e0 = bqd_exp_cof(eps, l, 0), e1 = bqd_exp_cof(eps, l, 1);
    BQD res;
    if (nested_at(N, l)) {
        BQD u, rho;
        CALL(bqd_sm_nest, N, l, &u, &rho);
        limdd_refs_push(rho);
        const BQD_EXP d = mtbdd_refs_push(CALL(bqd_exp_sub, e1, e0, r));
        limdd_refs_spawn(SPAWN(bqd_sm_phase_mul, rho, d, beta, r));
        const BQD lo = limdd_refs_push(CALL(bqd_sm_phase_mul, u, e0, beta, r));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_sm_phase_mul)));
        res = CALL(bqd_sm_mk_nested, l, lo, hi);
        limdd_refs_pop(3);
        mtbdd_refs_pop(1);
    } else {
        limdd_refs_spawn(SPAWN(bqd_sm_phase_mul, limdd_node_high(N), e1, beta, r));
        const BQD lo = limdd_refs_push(CALL(bqd_sm_phase_mul, limdd_node_low(N), e0, beta, r));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_sm_phase_mul)));
        res = bqd_sm_rebuild_s(l, lo, hi);
        limdd_refs_pop(2);
    }
    cache_put3(CACHE_BQD_SM_PHASEMUL, N, eps, beta, (uint64_t)res);
    return scale(scalar_of(k), res);
}

/**
 * Exp of skip:alg:smdiag, the canonical edge of beta^eps, of full support:
 * at the exponent's top v the low cofactor is beta^{eps_0} and the ratio
 * beta^{eps_1 - eps_0}, and the level is skipped where that is the constant
 * one, as mk_nested would on full support. Memoised on beta, eps and r, as
 * the copy rule's; the two calls overlap.
 */
TASK_IMPL_3(BQD, bqd_sm_exp, BQD_EXP, eps, EVBDD_WGT, beta, uint32_t, r)
{
    if (bqd_exp_is_const(eps)) return edge(bqd_exp_power(beta, bqd_exp_value(eps)), LIMDD_TERMINAL);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_SM_EXP, beta, eps, r, &hit)) return (BQD)hit;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PHASEMUL);
    const uint32_t v = bqd_exp_var(eps);
    const BQD_EXP e0 = bqd_exp_cof(eps, v, 0), e1 = bqd_exp_cof(eps, v, 1);
    const BQD_EXP d = mtbdd_refs_push(CALL(bqd_exp_sub, e1, e0, r));
    limdd_refs_spawn(SPAWN(bqd_sm_exp, e0, beta, r));
    const BQD h = limdd_refs_push(CALL(bqd_sm_exp, d, beta, r));
    const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_sm_exp)));
    const BQD res = (h == limdd_one_edge())
                  ? lo : edge(scalar_of(lo), bqd_sm_mk(v, false, unit(limdd_target(lo)), h));
    limdd_refs_pop(2);
    mtbdd_refs_pop(1);
    cache_put3(CACHE_BQD_SM_EXP, beta, eps, r, (uint64_t)res);
    return res;
}

/**
 * MulOff of skip:alg:smdiag, [k] . c^{[x not in supp [u]]}: k where c is 1, u
 * has full support or supp [k] lies inside supp [u], c k where u is zero, and
 * otherwise the phase multiplication by c^{1 - iota(u)}. Only the public entry
 * point calls it: no operation of SM multiplies off a support.
 */
TASK_IMPL_3(BQD, bqd_sm_mul_off, BQD, k, EVBDD_WGT, c, LIMDD_TARG, u)
{
    if (limdd_edge_is_zero(k) || c == EVBDD_ONE) return k;
    if (u == 0) return scale(c, k);
    if (bqd_sm_full(u) || CALL(bqd_sm_sub, limdd_target(k), u)) return k;
    const uint32_t r = bqd_exp_order(c);
    const BQD_EXP one = mtbdd_refs_push(bqd_exp_const(1, r));
    const BQD_EXP io = mtbdd_refs_push(CALL(bqd_sm_iota, u));
    const BQD_EXP eps = mtbdd_refs_push(CALL(bqd_exp_sub, one, io, r));
    const BQD res = CALL(bqd_sm_phase_mul, k, eps, c, r);
    mtbdd_refs_pop(3);
    return res;
}

/* --- building ----------------------------------------------------------------- */

static bool
all_zero(const EVBDD_WGT *g, uint64_t len)
{
    for (uint64_t y = 0; y < len; y++) if (g[y] != EVBDD_ZERO) return false;
    return true;
}

/**
 * The node of a representative g at level var, of length 2^(n-var): Build of
 * skip:alg:smbuild, skip:def:sm read as code. Equal cofactors skip the level.
 * A zero low cofactor makes an S node over the high cofactor, which is then a
 * representative itself, since g's least point is in it. Otherwise the low
 * cofactor is a representative, and the low edge the identity on its node;
 * the high edge is the representative's edge of the ratio, zero off the low
 * support, where the high support lies inside the low one, and of the high
 * cofactor itself where it does not. The low child is spawned and the high
 * one computed, as in the copy rule's builder.
 */
TASK_3(LIMDD_TARG, bqd_sm_build, const EVBDD_WGT *, g, uint32_t, var, uint32_t, n)
{
    if (var == n) return LIMDD_TERMINAL;                 /* g = [1] */

    const uint64_t h = UINT64_C(1) << (n - var - 1);
    const EVBDD_WGT *f0 = g, *f1 = g + h;
    if (memcmp(f0, f1, h * sizeof(EVBDD_WGT)) == 0) return CALL(bqd_sm_build, f0, var + 1, n);
    if (all_zero(f0, h)) {
        const LIMDD_TARG R = CALL(bqd_sm_build, f1, var + 1, n);
        return bqd_sm_mk(var, true, limdd_zero_edge(), unit(R));
    }

    SPAWN(bqd_sm_build, f0, var + 1, n);

    bool s = false;
    LIMDD high = limdd_zero_edge();
    if (!all_zero(f1, h)) {
        bool nested = true;
        for (uint64_t y = 0; y < h && nested; y++)
            if (f0[y] == EVBDD_ZERO && f1[y] != EVBDD_ZERO) nested = false;
        EVBDD_WGT *raw = malloc(h * sizeof(EVBDD_WGT));
        if (raw == NULL) { fprintf(stderr, "sylvan: out of memory in a BQD build\n"); exit(1); }
        for (uint64_t y = 0; y < h; y++) {
            if (!nested) raw[y] = f1[y];
            else raw[y] = (f1[y] == EVBDD_ZERO) ? EVBDD_ZERO : wgt_div(f1[y], f0[y]);
        }
        EVBDD_WGT c0; uint64_t zs, zp;
        bqd_representative(raw, h, raw, &c0, &zs, &zp);
        assert(zs == 0 && zp == 0);                      /* the scalar family's representative */
        const LIMDD_TARG ht = CALL(bqd_sm_build, raw, var + 1, n);
        high = edge(c0, ht);
        free(raw);
        s = !nested;
    }

    const LIMDD low = unit(SYNC(bqd_sm_build));
    return bqd_sm_mk(var, s, low, high);
}

TASK_IMPL_2(BQD, bqd_sm_from_vector, const EVBDD_WGT *, f, uint32_t, n)
{
    const uint64_t len = UINT64_C(1) << n;
    if (all_zero(f, len)) return limdd_zero_edge();
    EVBDD_WGT *rep = malloc(len * sizeof(EVBDD_WGT));
    if (rep == NULL) { fprintf(stderr, "sylvan: out of memory in a BQD build\n"); exit(1); }
    EVBDD_WGT c0; uint64_t s, p;
    bqd_representative(f, len, rep, &c0, &s, &p);
    const LIMDD_TARG root = CALL(bqd_sm_build, rep, 0, n);
    free(rep);
    return edge(c0, root);
}

/* --- decoding --------------------------------------------------------------- */

static void
fill_zero(EVBDD_WGT *out, uint64_t len)
{
    for (uint64_t y = 0; y < len; y++) out[y] = EVBDD_ZERO;
}

/**
 * The function of p read at level var, into out[0..2^(n-var)), as the copy
 * rule's decoder, with the high half by the tag: c_1 [A] [R] for a Q node and
 * c_1 [R] for an S node, the plain product in both. The levels an edge skips
 * are the leading bits of the index, over which p's block repeats.
 */
VOID_TASK_4(bqd_sm_decode, LIMDD_TARG, p, uint32_t, var, uint32_t, n, EVBDD_WGT *, out)
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

        if (!lz) SPAWN(bqd_sm_decode, limdd_target(low), v + 1, n, out);

        EVBDD_WGT *g1 = NULL;
        if (!hz) {
            g1 = malloc(h * sizeof(EVBDD_WGT));
            if (g1 == NULL) { fprintf(stderr, "sylvan: out of memory decoding a BQD\n"); exit(1); }
            CALL(bqd_sm_decode, limdd_target(high), v + 1, n, g1);
        }

        if (!lz) SYNC(bqd_sm_decode); else fill_zero(out, h);

        if (hz) {
            fill_zero(out + h, h);
        } else {
            const EVBDD_WGT c1 = scalar_of(high);
            const bool s = bqd_sm_is_s(p);
            for (uint64_t y = 0; y < h; y++) {
                const EVBDD_WGT a = out[y], b = g1[y];
                EVBDD_WGT w = EVBDD_ZERO;
                if (b != EVBDD_ZERO && (s || a != EVBDD_ZERO)) w = wgt_mul(c1, s ? b : wgt_mul(a, b));
                out[h + y] = w;
            }
            free(g1);
        }
    }

    for (uint64_t y = len; y < (UINT64_C(1) << (n - var)); y += len)
        memcpy(out + y, out, len * sizeof(EVBDD_WGT));
}

VOID_TASK_IMPL_3(bqd_sm_to_vector, BQD, e, uint32_t, n, EVBDD_WGT *, out)
{
    const uint64_t len = UINT64_C(1) << n;
    if (limdd_edge_is_zero(e)) { fill_zero(out, len); return; }
    CALL(bqd_sm_decode, limdd_target(e), 0, n, out);
    const EVBDD_WGT c = scalar_of(e);
    if (c == EVBDD_ONE) return;
    for (uint64_t y = 0; y < len; y++)
        if (out[y] != EVBDD_ZERO) out[y] = wgt_mul(c, out[y]);
}

/* --- one amplitude ------------------------------------------------------------ */

/**
 * The value of p's function at y. The low half is one path; the high half of
 * an S node is one path too, and that of a Q node needs both children at the
 * same point, unless the low one is zero there.
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
    EVBDD_WGT v = EVBDD_ONE;
    if (!bqd_sm_is_s(p)) {
        v = eval_node(limdd_target(low), n, yy);
        if (v == EVBDD_ZERO) return v;
    }
    const EVBDD_WGT b = eval_node(limdd_target(high), n, yy);
    if (b == EVBDD_ZERO) return b;
    return wgt_mul(scalar_of(high), wgt_mul(v, b));
}

EVBDD_WGT
bqd_sm_eval(BQD e, uint32_t n, uint64_t x)
{
    if (limdd_edge_is_zero(e)) return EVBDD_ZERO;
    const EVBDD_WGT v = eval_node(limdd_target(e), n, x);
    return v == EVBDD_ZERO ? v : wgt_mul(scalar_of(e), v);
}

/* --- states ------------------------------------------------------------------- */

/*
 * At a 1 the low cofactor is zero, an S node over the rest; at a 0 the high
 * one is, a Q node with no ratio. Both depend on their variable, so nothing
 * is skipped.
 */
BQD
bqd_sm_basis_state(uint64_t x, uint32_t n)
{
    LIMDD_TARG t = LIMDD_TERMINAL;
    for (uint32_t v = n; v-- > 0; ) {
        const bool one = (x >> (n - 1 - v)) & 1;
        t = one ? bqd_sm_mk(v, true, limdd_zero_edge(), unit(t))
                : bqd_sm_mk(v, false, unit(t), limdd_zero_edge());
    }
    return unit(t);
}
