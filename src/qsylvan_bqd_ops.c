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

#include "qsylvan_bqd_ops.h"
#include "qsylvan_limdd_gc.h"

static inline EVBDD_WGT
scalar_of(BQD e)
{
    return limdd_lim_weight(limdd_label(e));
}

static inline BQD
with_scalar(EVBDD_WGT c, LIMDD_TARG t, uint32_t n)
{
    return limdd_bundle(bqd_lim_make(c, 0, 0, n), t);
}

/**
 * prop:prodscalar. The node functions of the two operands are
 * representatives, so their product is memoised on the pair of targets and
 * the two edge scalars are multiplied on afterwards. What the memo holds is
 * the canonical edge of that product: its scalar is the value of the product
 * at 0, which the parent needs, and its node is a representative.
 *
 * With full support the copy never fires, so the high child of the product
 * is the product of the two high children: the ratio of the products is the
 * product of the ratios, with no cross term, which is the identity the
 * quotient rule has and the moment rule does not.
 */
TASK_IMPL_3(BQD, bqd_product, BQD, f, BQD, g, uint32_t, var)
{
    assert(bqd_family() == BQD_FAMILY_SCALAR && "prop:prodscalar is for scalar labels");
    if (limdd_edge_is_zero(f) || limdd_edge_is_zero(g)) return limdd_zero_edge();

    const uint32_t  n  = (uint32_t)limdd_lims_nqubits();
    const EVBDD_WGT c  = wgt_mul(scalar_of(f), scalar_of(g));
    const LIMDD_TARG tf = limdd_target(f), tg = limdd_target(g);

    /*
     * The terminal is the base case, not a level count: the operands may be
     * narrower than the table, and the table's width says nothing about where
     * a given vector's bottom is. Both reach it together, since both are at
     * the same level.
     */
    if (tf == LIMDD_TERMINAL || tg == LIMDD_TERMINAL) {
        assert(tf == LIMDD_TERMINAL && tg == LIMDD_TERMINAL);
        return with_scalar(c, LIMDD_TERMINAL, n);
    }
    assert(limdd_node_var(tf) == var && limdd_node_var(tg) == var);

    uint64_t hit;
    if (cache_get3(CACHE_BQD_PROD, tf, tg, var, &hit)) {
        const BQD r = (BQD)hit;
        return with_scalar(wgt_mul(c, scalar_of(r)), limdd_target(r), n);
    }

    const LIMDD fl = limdd_node_low(tf), fh = limdd_node_high(tf);
    const LIMDD gl = limdd_node_low(tg), gh = limdd_node_high(tg);
    assert(!limdd_edge_is_zero(fl) && !limdd_edge_is_zero(fh) &&
           !limdd_edge_is_zero(gl) && !limdd_edge_is_zero(gh) &&
           "the product recursion needs full support");

    /* the high children without their scalars: those go on the result edge */
    const LIMDD v = limdd_bundle(LIMDD_LIM_IDENTITY, limdd_target(fh));
    const LIMDD q = limdd_bundle(LIMDD_LIM_IDENTITY, limdd_target(gh));

    /* the two child products are independent: overlap them, as limdd_plus does */
    limdd_refs_spawn(SPAWN(bqd_product, fl, gl, var + 1));
    const BQD hi = limdd_refs_push(CALL(bqd_product, v, q, var + 1));
    const BQD lo = limdd_refs_sync(SYNC(bqd_product));
    limdd_refs_pop(1);

    /* u.p = a U  and  v.q = b Q, both canonical; the node of the product is
     * U below with ratio c1 d1 b Q, and the product itself is a times that */
    const EVBDD_WGT a = scalar_of(lo), b = scalar_of(hi);
    const EVBDD_WGT ratio = wgt_mul(wgt_mul(scalar_of(fh), scalar_of(gh)), b);
    const LIMDD_TARG node = limdd_makenode(var,
        limdd_bundle(LIMDD_LIM_IDENTITY, limdd_target(lo)),
        with_scalar(ratio, limdd_target(hi), n));

    cache_put3(CACHE_BQD_PROD, tf, tg, var, (uint64_t)with_scalar(a, node, n));
    return with_scalar(wgt_mul(c, a), node, n);
}

LIMDD_TARG
bqd_ones(uint32_t var, uint32_t n)
{
    if (var == n) return LIMDD_TERMINAL;
    const LIMDD e = limdd_bundle(LIMDD_LIM_IDENTITY, bqd_ones(var + 1, n));
    return limdd_makenode(var, e, e);        /* both cofactors one, ratio one */
}

/**
 * At a variable outside the monomial the two cofactors agree, so the ratio is
 * the constant one. At a variable inside it the low cofactor is one, since
 * the monomial is then zero, and the ratio is the monomial on what is left,
 * whose value at 0 is one, so it is its own representative with scalar one,
 * until the last variable of the monomial, where it is the phase itself.
 */
static LIMDD_TARG
monomial_node(uint64_t A, EVBDD_WGT phase, uint32_t var, uint32_t n)
{
    if (var == n) return LIMDD_TERMINAL;
    const uint64_t bit = UINT64_C(1) << (n - 1 - var);
    if ((A & bit) == 0) {
        const LIMDD low = limdd_bundle(LIMDD_LIM_IDENTITY, monomial_node(A, phase, var + 1, n));
        return limdd_makenode(var, low, limdd_bundle(LIMDD_LIM_IDENTITY, bqd_ones(var + 1, n)));
    }
    const uint64_t rest = A & ~bit;
    const LIMDD low = limdd_bundle(LIMDD_LIM_IDENTITY, bqd_ones(var + 1, n));
    const LIMDD high = (rest == 0)
        ? with_scalar(phase, bqd_ones(var + 1, n), n)
        : limdd_bundle(LIMDD_LIM_IDENTITY, monomial_node(rest, phase, var + 1, n));
    return limdd_makenode(var, low, high);
}

BQD
bqd_monomial(uint64_t A, EVBDD_WGT phase, uint32_t n)
{
    assert(bqd_family() == BQD_FAMILY_SCALAR);
    if (A == 0) return with_scalar(phase, bqd_ones(0, n), n);
    return limdd_bundle(LIMDD_LIM_IDENTITY, monomial_node(A, phase, 0, n));
}

BQD
bqd_apply_diagonal(BQD e, uint64_t A, EVBDD_WGT phase, uint32_t n)
{
    return bqd_product(e, bqd_monomial(A, phase, n), 0);
}
