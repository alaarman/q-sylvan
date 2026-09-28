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
 * The node at level var with low child u and high edge r, the ratio, on full
 * support: FromRatio of skip:alg:constructors with its flag set. The node is
 * redundant exactly when r is the constant one (skip:cor:full), since the
 * indicator of u's support is then the terminal, and its function is then
 * u's, extended to var: the level is skipped and u is returned. One
 * comparison, and no walk of u.
 */
static LIMDD_TARG
node_on_full(uint32_t var, LIMDD_TARG u, LIMDD r)
{
    if (r == limdd_one_edge()) return u;
    return limdd_makenode(var, limdd_bundle(LIMDD_LIM_IDENTITY, u), r);
}

/**
 * The low and high edge of node t at level var: its own where it sits there,
 * and where it is below, the virtual node of skip:lem:virtual, whose low edge
 * is t and whose ratio is the constant one. That is the virtual node on full
 * support, where the indicator of t's support is the terminal.
 */
static void
edges_at(LIMDD_TARG t, uint32_t var, LIMDD *low, LIMDD *high, const char *op)
{
    if (limdd_node_var(t) == var) {
        *low = limdd_node_low(t);
        *high = limdd_node_high(t);
        require_full_support(*low, *high, op);
    } else {
        *low = limdd_bundle(LIMDD_LIM_IDENTITY, t);
        *high = limdd_one_edge();
    }
}

/**
 * prop:prodscalar, as Prod of skip:alg:diag. The node functions of the two
 * operands are representatives, so their product is memoised on the pair of
 * targets and the two edge scalars are multiplied on afterwards. What the
 * memo holds is the canonical edge of that product: its scalar is the value
 * of the product at 0, which the parent needs, and its node is a
 * representative. The pair is unordered, since the product commutes, and the
 * key holds no level (skip:lem:ext).
 *
 * With full support the copy never fires, so the high child of the product
 * is the product of the two high children: the ratio of the products is the
 * product of the ratios, with no cross term, which is the identity the
 * quotient rule has and the moment rule does not.
 *
 * The recursion is at the higher of the two nodes' levels. An operand whose
 * node is lower skips that level and takes part as its virtual node, and the
 * terminal, the constant one, is the unit of the product. That is the base
 * case, not a level count: the operands may be narrower than the table, and
 * the table's width says nothing about where a given vector's bottom is.
 */
TASK_2(BQD, bqd_product_rec, BQD, f, BQD, g)
{
    if (limdd_edge_is_zero(f) || limdd_edge_is_zero(g)) return limdd_zero_edge();

    const uint32_t  n  = (uint32_t)limdd_lims_nqubits();
    const EVBDD_WGT c  = wgt_mul(scalar_of(f), scalar_of(g));
    LIMDD_TARG tf = limdd_target(f), tg = limdd_target(g);
    if (tf == LIMDD_TERMINAL) return with_scalar(c, tg, n);
    if (tg == LIMDD_TERMINAL) return with_scalar(c, tf, n);
    if (tf > tg) { const LIMDD_TARG t = tf; tf = tg; tg = t; }

    uint64_t hit;
    if (cache_get3(CACHE_BQD_PROD, tf, tg, 0, &hit)) {
        const BQD r = (BQD)hit;
        return with_scalar(wgt_mul(c, scalar_of(r)), limdd_target(r), n);
    }

    const uint32_t lf = limdd_node_var(tf), lg = limdd_node_var(tg);
    const uint32_t var = lf < lg ? lf : lg;
    LIMDD fl, fh, gl, gh;
    edges_at(tf, var, &fl, &fh, "bqd_product");
    edges_at(tg, var, &gl, &gh, "bqd_product");

    /* the high children without their scalars: those go on the result edge */
    const LIMDD v = limdd_bundle(LIMDD_LIM_IDENTITY, limdd_target(fh));
    const LIMDD q = limdd_bundle(LIMDD_LIM_IDENTITY, limdd_target(gh));

    /* the two child products are independent: overlap them, as limdd_plus does */
    limdd_refs_spawn(SPAWN(bqd_product_rec, fl, gl));
    const BQD hi = limdd_refs_push(CALL(bqd_product_rec, v, q));
    const BQD lo = limdd_refs_sync(SYNC(bqd_product_rec));
    limdd_refs_pop(1);

    /* u.p = a U  and  v.q = b Q, both canonical; the node of the product is
     * U below with ratio c1 d1 b Q, and the product itself is a times that.
     * Where that ratio is the constant one the product does not depend on
     * the variable, and its node is U */
    const EVBDD_WGT a = scalar_of(lo), b = scalar_of(hi);
    const EVBDD_WGT ratio = wgt_mul(wgt_mul(scalar_of(fh), scalar_of(gh)), b);
    const LIMDD_TARG node = node_on_full(var, limdd_target(lo),
                                         with_scalar(ratio, limdd_target(hi), n));

    cache_put3(CACHE_BQD_PROD, tf, tg, 0, (uint64_t)with_scalar(a, node, n));
    return with_scalar(wgt_mul(c, a), node, n);
}

TASK_IMPL_2(BQD, bqd_product, BQD, f, BQD, g)
{
    require_scalar_family("bqd_product");
    if (limdd_edge_is_zero(f) || limdd_edge_is_zero(g)) return limdd_zero_edge();
    /* the hypothesis for the whole recursion, tested once: a zero anywhere
     * below is where the ratio of the products stops being the product */
    if (!(bqd_has_full_support(f) && bqd_has_full_support(g))) return bqd_multiply(f, g);
    return CALL(bqd_product_rec, f, g);
}

/*
 * The constant one depends on no variable, so every level is skipped and its
 * node, from any level down, is the terminal (skip:def:fr).
 */
LIMDD_TARG
bqd_ones(uint32_t var, uint32_t n)
{
    assert(var <= n);
    (void)var; (void)n;
    return LIMDD_TERMINAL;
}

/**
 * Mono of skip:alg:diag, the monomial phase^{x_A} as an edge: one node at
 * each variable of A, and none elsewhere, since outside A the two cofactors
 * agree and the level is skipped. At a variable of A the low cofactor is the
 * constant one, since the monomial is then zero, and the ratio is the
 * monomial on the variables of A below it, whose value at 0 is one, so it is
 * its own representative with scalar one, until the last variable of A,
 * where it is the phase itself. So the nodes are made from the bottom up,
 * each the high child of the next. A phase of one makes the constant one,
 * which is the terminal, whatever A is.
 */
static BQD
monomial_edge(uint64_t A, EVBDD_WGT phase, uint32_t n)
{
    BQD h = with_scalar(phase, LIMDD_TERMINAL, n);
    if (phase == EVBDD_ONE) return h;
    /* qubit q is bit n-1-q, so the lowest variable is the lowest bit */
    for (uint64_t rest = A; rest != 0; rest &= rest - 1) {
        const uint32_t var = n - 1 - (uint32_t)__builtin_ctzll(rest);
        h = limdd_bundle(LIMDD_LIM_IDENTITY, limdd_makenode(var, limdd_one_edge(), h));
    }
    return h;
}

BQD
bqd_monomial(uint64_t A, EVBDD_WGT phase, uint32_t n)
{
    require_scalar_family("bqd_monomial");
    return monomial_edge(A, phase, n);
}

/**
 * prop:diag, as the paper's proof runs it: the product of e with the
 * monomial phase^{x_A}, where the gate's operand is never built. It is Diag
 * of skip:alg:diag.
 *
 * In the recursion of prop:prodscalar the gate's node has the constant one as
 * its high child at a variable outside A, where the monomial skips and takes
 * part as its virtual node, and as its low child at a variable in A, and a
 * product with the constant one returns the other operand. So of the two
 * child products at each node one is free, and what is left is a single
 * path: down the low edge outside A, down the high edge in A, and it stops at
 * the last variable of A, where the ratio is the phase times one. bqd_product
 * with the gate's diagram takes the same path, since the terminal is its
 * unit, and one node pair more where the walk stops at a node; but it builds
 * the diagram first and goes through the memo at every step, where this is
 * n steps at most and nothing else.
 *
 * With the node of e being L below and h.H above, the node of the product is
 *
 *     outside A:           U below, h.H above       (U = L x_A, the ratio kept)
 *     in A, more to come:  L below, h.V above       (V = H x_rest)
 *     the last of A:       L below, (h phase).H above
 *
 * Where e skips the first variable a of A that is left, its node N, or the
 * terminal, is below a, and the product has a node at a with N below and the
 * monomial of the rest of A above: the low cofactor of x_A is zero there, so
 * the product's is e itself, and the ratio of the two cofactors is x_rest on
 * full support. The walk stops there, and the node of e it stops at is not
 * counted as visited.
 *
 * Every node is made by node_on_full, since a ratio can come out as the
 * constant one, and the product then does not depend on the variable: in A
 * where h.V or (h phase).H is one, and at the stop where the phase is. Outside
 * A the ratio is e's own, which is not the constant one on a fully reduced e,
 * and the test there costs a comparison.
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

    /* a, the first variable of A, is its highest bit, since qubit q is bit n-1-q */
    const uint32_t hb = 63 - (uint32_t)__builtin_clzll(A);
    const uint64_t bit = UINT64_C(1) << hb;
    const uint32_t a = n - 1 - hb;
    const uint32_t var = limdd_level(node);
    if (var > a) {
        const LIMDD_TARG r = node_on_full(a, node, monomial_edge(A & ~bit, phase, n));
        return with_scalar(c, r, n);
    }

    const LIMDD low = limdd_node_low(node), high = limdd_node_high(node);
    require_full_support(low, high, "bqd_apply_diagonal");
    if (visits != NULL) (*visits)++;

    if (var < a) {
        const BQD lo = diag_rec(limdd_bundle(LIMDD_LIM_IDENTITY, limdd_target(low)),
                                A, phase, n, visits);
        assert(scalar_of(lo) == EVBDD_ONE);
        return with_scalar(c, node_on_full(var, limdd_target(lo), high), n);
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
    return with_scalar(c, node_on_full(var, limdd_target(low), up), n);
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
