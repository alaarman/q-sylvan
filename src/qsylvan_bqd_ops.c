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

#include "qsylvan_bqd_ops.h"
#include "qsylvan_bqd_gates.h"
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

/*
 * The two hypotheses of prop:prodscalar, checked in every build and not by
 * assert. Outside them the recursion does not fail, it returns a diagram of
 * the wrong function: with a zero cofactor the copy fires and the ratio of
 * the products is no longer the product of the ratios, and a translation or
 * sign label is read as a scalar. The entry points test full support of the
 * whole diagram and take the general product of qsylvan_bqd_gates.h when it
 * fails; the per-node check below is what is left, a guard that fires only
 * if that routing is wrong.
 */
static void
require_scalar_family(const char *op)
{
    if (bqd_family() != BQD_FAMILY_SCALAR) {
        fprintf(stderr, "sylvan: %s is proved for scalar labels only (prop:prodscalar), "
                        "and this session holds a %s\n", op, bqd_family_name(bqd_family()));
        exit(1);
    }
}

static void
require_full_support(LIMDD low, LIMDD high, const char *op)
{
    if (limdd_edge_is_zero(low) || limdd_edge_is_zero(high)) {
        fprintf(stderr, "sylvan: %s needs operands of full support (prop:prodscalar), "
                        "and one has a zero cofactor\n", op);
        exit(1);
    }
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
    require_scalar_family("bqd_product");
    if (limdd_edge_is_zero(f) || limdd_edge_is_zero(g)) return limdd_zero_edge();
    /* at the root, the hypothesis for the whole recursion: a zero anywhere
     * below is where the ratio of the products stops being the product */
    if (var == 0 && !(bqd_has_full_support(f) && bqd_has_full_support(g)))
        return bqd_multiply(f, g);

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
    require_full_support(fl, fh, "bqd_product");
    require_full_support(gl, gh, "bqd_product");

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
    require_scalar_family("bqd_monomial");
    if (A == 0) return with_scalar(phase, bqd_ones(0, n), n);
    return limdd_bundle(LIMDD_LIM_IDENTITY, monomial_node(A, phase, 0, n));
}

/**
 * prop:diag, as the paper's proof runs it: the product of e with the
 * monomial phase^{x_A}, where the gate's operand is never built.
 *
 * In the recursion of prop:prodscalar the gate's node has the constant one as
 * its high child at a variable outside A and as its low child at a variable
 * in A, and a product with the constant one returns the other operand. So of
 * the two child products at each node one is free, and what is left is a
 * single path: down the low edge outside A, down the high edge in A, and it
 * stops at the last variable of A, where the ratio is the phase times one.
 * bqd_product cannot see this, since it has no way to tell the constant one
 * from any other node, and would pair every node of e with it: |e| memo
 * entries for what is n steps here.
 *
 * With the node of e being L below and h.H above, the node of the product is
 *
 *     outside A:           U below, h.H above       (U = L x_A, the ratio kept)
 *     in A, more to come:  L below, h.V above       (V = H x_rest)
 *     the last of A:       L below, (h phase).H above
 *
 * No scalar comes back from below. Every node here is a representative, 1
 * at the point 0, and so is x_A while A is not empty, so the product of the
 * two is 1 at 0 as well and its canonical edge carries the identity. Each
 * node made is a representative for the same reason.
 */
static BQD
diag_rec(BQD e, uint64_t A, EVBDD_WGT phase, uint32_t n, uint32_t *visits)
{
    if (limdd_edge_is_zero(e)) return e;
    const EVBDD_WGT c = scalar_of(e);
    const LIMDD_TARG node = limdd_target(e);
    if (A == 0) return with_scalar(wgt_mul(c, phase), node, n);   /* x_A = 1 */

    /* A is a set of qubits of this vector, so the walk reaches the last of
     * them before the terminal */
    assert(node != LIMDD_TERMINAL);
    const uint32_t var = limdd_node_var(node);
    const uint64_t bit = UINT64_C(1) << (n - 1 - var);
    const LIMDD low = limdd_node_low(node), high = limdd_node_high(node);
    require_full_support(low, high, "bqd_apply_diagonal");
    if (visits != NULL) (*visits)++;

    if ((A & bit) == 0) {
        const BQD lo = diag_rec(limdd_bundle(LIMDD_LIM_IDENTITY, limdd_target(low)),
                                A, phase, n, visits);
        assert(scalar_of(lo) == EVBDD_ONE);
        const LIMDD_TARG r = limdd_makenode(var,
            limdd_bundle(LIMDD_LIM_IDENTITY, limdd_target(lo)), high);
        return with_scalar(c, r, n);
    }
    const uint64_t rest = A & ~bit;
    LIMDD up;
    if (rest == 0) {
        up = with_scalar(wgt_mul(scalar_of(high), phase), limdd_target(high), n);
    } else {
        const BQD hi = diag_rec(limdd_bundle(LIMDD_LIM_IDENTITY, limdd_target(high)),
                                rest, phase, n, visits);
        assert(scalar_of(hi) == EVBDD_ONE);
        up = limdd_bundle(limdd_label(high), limdd_target(hi));
    }
    return with_scalar(c, limdd_makenode(var, low, up), n);
}

/*
 * The walk needs full support of all of e, not only of the path it takes:
 * outside A it keeps e's high edge, and that edge is the new ratio only where
 * the low child has no zero, since at a zero the copy fires and returns the
 * old high value without the phase. A zero off the path is therefore a wrong
 * result and not a failed check. So the whole diagram is tested, which is a
 * walk memoised on nodes, and a state with a zero takes the general product.
 */
BQD
bqd_apply_diagonal_counted(BQD e, uint64_t A, EVBDD_WGT phase, uint32_t n, uint32_t *visits)
{
    require_scalar_family("bqd_apply_diagonal");
    if (visits != NULL) *visits = 0;
    if (!bqd_has_full_support(e)) {
        if (limdd_edge_is_zero(e)) return e;
        return bqd_multiply(e, bqd_monomial(A, phase, n));
    }
    return diag_rec(e, A, phase, n, visits);
}

BQD
bqd_apply_diagonal(BQD e, uint64_t A, EVBDD_WGT phase, uint32_t n)
{
    return bqd_apply_diagonal_counted(e, A, phase, n, NULL);
}
