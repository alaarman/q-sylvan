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

#include "qsylvan_bqd_exp.h"

/* --- scalars --------------------------------------------------------------- */

uint32_t
bqd_exp_order(EVBDD_WGT beta)
{
    EVBDD_WGT p = beta;
    for (uint32_t k = 1; k <= 8; k++) {
        if (p == EVBDD_ONE) return k;
        p = wgt_mul(p, beta);
    }
    return 0;
}

/* By squaring, through wgt_mul, whose results the operation cache keeps. */
EVBDD_WGT
bqd_exp_power(EVBDD_WGT beta, int64_t v)
{
    uint64_t k = (v < 0) ? (uint64_t)(-(v + 1)) + 1 : (uint64_t)v;
    EVBDD_WGT r = EVBDD_ONE, b = beta;
    while (k != 0) {
        if (k & 1) r = wgt_mul(r, b);
        k >>= 1;
        if (k != 0) b = wgt_mul(b, b);
    }
    return (v < 0) ? wgt_div(EVBDD_ONE, r) : r;
}

/* --- constants and nodes ----------------------------------------------------- */

static inline int64_t
reduce(int64_t v, uint32_t r)
{
    if (r == 0) return v;
    const int64_t m = v % (int64_t)r;
    return m < 0 ? m + (int64_t)r : m;
}

BQD_EXP
bqd_exp_const(int64_t v, uint32_t r)
{
    return mtbdd_int64(reduce(v, r));
}

int64_t
bqd_exp_value(BQD_EXP e)
{
    return mtbdd_getint64(e);
}

uint32_t
bqd_exp_var(BQD_EXP e)
{
    return mtbdd_isleaf(e) ? UINT32_MAX : mtbdd_getvar(e);
}

BQD_EXP
bqd_exp_cof(BQD_EXP e, uint32_t var, int b)
{
    if (mtbdd_isleaf(e) || mtbdd_getvar(e) > var) return e;
    assert(mtbdd_getvar(e) == var && "an exponent is read at or above its root");
    return b ? mtbdd_gethigh(e) : mtbdd_getlow(e);
}

int64_t
bqd_exp_eval(BQD_EXP e, uint64_t y)
{
    while (!mtbdd_isleaf(e)) e = ((y >> mtbdd_getvar(e)) & 1) ? mtbdd_gethigh(e) : mtbdd_getlow(e);
    return mtbdd_getint64(e);
}

/*
 * For the tests: a Sylvan collection at every k-th exponent node that the
 * worker which set the knob makes, so that collections happen inside the
 * operations, where they would otherwise wait for Sylvan's table to fill. One
 * worker only, since two that ask at once can wait on each other in
 * sylvan_gc. Off, 0, by default; the count is the worker's own.
 */
static uint64_t collect_every = 0;
static unsigned collect_worker = 0;
static SYLVAN_TLS uint64_t nodes_made = 0;

void
bqd_exp_collect_every(uint64_t k)
{
    WorkerP *w = lace_get_worker();
    collect_worker = (w != NULL) ? (unsigned)w->worker : 0;
    collect_every = k;
}

BQD_EXP
bqd_exp_node(uint32_t var, BQD_EXP lo, BQD_EXP hi)
{
    if (collect_every != 0 && lo != hi && ++nodes_made % collect_every == 0) {
        WorkerP *w = lace_get_worker();
        if (w != NULL && (unsigned)w->worker == collect_worker) {
            mtbdd_refs_push(lo);
            mtbdd_refs_push(hi);
            sylvan_gc();
            mtbdd_refs_pop(2);
        }
    }
    /* mtbdd_makenode reduces lo == hi, and protects both across a collection */
    return mtbdd_makenode(var, lo, hi);
}

/* --- leafwise operations ----------------------------------------------------- */

/*
 * a - b and a + b, the recursion of mtbdd_apply with the reduction modulo r at
 * the leaves, on the higher of the two roots. Memoised on (a, b, r); the two
 * cofactor results are independent and overlap. An exponent is never a
 * complemented edge, so the first key word is an index below 2^40.
 */
TASK_IMPL_3(BQD_EXP, bqd_exp_sub, BQD_EXP, a, BQD_EXP, b, uint32_t, r)
{
    const bool la = mtbdd_isleaf(a), lb = mtbdd_isleaf(b);
    if (lb && mtbdd_getint64(b) == 0) return a;
    if (a == b) return bqd_exp_const(0, r);
    if (la && lb) return bqd_exp_const(mtbdd_getint64(a) - mtbdd_getint64(b), r);

    uint64_t hit;
    if (cache_get3(CACHE_BQD_EXP_SUB, a, b, r, &hit)) return (BQD_EXP)hit;
    const uint32_t va = bqd_exp_var(a), vb = bqd_exp_var(b);
    const uint32_t v = va < vb ? va : vb;
    mtbdd_refs_spawn(SPAWN(bqd_exp_sub, bqd_exp_cof(a, v, 0), bqd_exp_cof(b, v, 0), r));
    const BQD_EXP hi = mtbdd_refs_push(CALL(bqd_exp_sub, bqd_exp_cof(a, v, 1),
                                            bqd_exp_cof(b, v, 1), r));
    const BQD_EXP lo = mtbdd_refs_sync(SYNC(bqd_exp_sub));
    mtbdd_refs_pop(1);
    const BQD_EXP res = bqd_exp_node(v, lo, hi);
    cache_put3(CACHE_BQD_EXP_SUB, a, b, r, res);
    return res;
}

TASK_IMPL_3(BQD_EXP, bqd_exp_add, BQD_EXP, a, BQD_EXP, b, uint32_t, r)
{
    const bool la = mtbdd_isleaf(a), lb = mtbdd_isleaf(b);
    if (lb && mtbdd_getint64(b) == 0) return a;
    if (la && mtbdd_getint64(a) == 0) return b;
    if (la && lb) return bqd_exp_const(mtbdd_getint64(a) + mtbdd_getint64(b), r);
    if (a > b) { const BQD_EXP t = a; a = b; b = t; }      /* it commutes */

    uint64_t hit;
    if (cache_get3(CACHE_BQD_EXP_ADD, a, b, r, &hit)) return (BQD_EXP)hit;
    const uint32_t va = bqd_exp_var(a), vb = bqd_exp_var(b);
    const uint32_t v = va < vb ? va : vb;
    mtbdd_refs_spawn(SPAWN(bqd_exp_add, bqd_exp_cof(a, v, 0), bqd_exp_cof(b, v, 0), r));
    const BQD_EXP hi = mtbdd_refs_push(CALL(bqd_exp_add, bqd_exp_cof(a, v, 1),
                                            bqd_exp_cof(b, v, 1), r));
    const BQD_EXP lo = mtbdd_refs_sync(SYNC(bqd_exp_add));
    mtbdd_refs_pop(1);
    const BQD_EXP res = bqd_exp_node(v, lo, hi);
    cache_put3(CACHE_BQD_EXP_ADD, a, b, r, res);
    return res;
}

/*
 * The constant v is made where iota is 0 and x is not already v, so it is
 * built inside the recursion and not handed in: an argument the caller made
 * would have to be protected, and one made here is used at once.
 */
TASK_IMPL_3(BQD_EXP, bqd_exp_sel, BQD_EXP, iota, BQD_EXP, x, int64_t, v)
{
    if (mtbdd_isleaf(iota)) return mtbdd_getint64(iota) ? x : mtbdd_int64(v);
    if (mtbdd_isleaf(x) && mtbdd_getint64(x) == v) return x;

    uint64_t hit;
    if (cache_get3(CACHE_BQD_EXP_SEL, iota, x, (uint64_t)v, &hit)) return (BQD_EXP)hit;
    const uint32_t vi = bqd_exp_var(iota), vx = bqd_exp_var(x);
    const uint32_t var = vi < vx ? vi : vx;
    mtbdd_refs_spawn(SPAWN(bqd_exp_sel, bqd_exp_cof(iota, var, 0), bqd_exp_cof(x, var, 0), v));
    const BQD_EXP hi = mtbdd_refs_push(CALL(bqd_exp_sel, bqd_exp_cof(iota, var, 1),
                                            bqd_exp_cof(x, var, 1), v));
    const BQD_EXP lo = mtbdd_refs_sync(SYNC(bqd_exp_sel));
    mtbdd_refs_pop(1);
    const BQD_EXP res = bqd_exp_node(var, lo, hi);
    cache_put3(CACHE_BQD_EXP_SEL, iota, x, (uint64_t)v, res);
    return res;
}

/*
 * The bits of t at the root's variable and below are the ones that act; the
 * key holds t with the others cleared, so that two translations that differ
 * only above the root share an entry.
 */
TASK_IMPL_2(BQD_EXP, bqd_exp_translate, BQD_EXP, e, uint64_t, t)
{
    if (mtbdd_isleaf(e)) return e;
    const uint32_t v = mtbdd_getvar(e);
    t &= ~((UINT64_C(1) << v) - 1);
    if (t == 0) return e;

    uint64_t hit;
    if (cache_get3(CACHE_BQD_EXP_TRANSLATE, e, t, 0, &hit)) return (BQD_EXP)hit;
    mtbdd_refs_spawn(SPAWN(bqd_exp_translate, mtbdd_getlow(e), t));
    const BQD_EXP hi = mtbdd_refs_push(CALL(bqd_exp_translate, mtbdd_gethigh(e), t));
    const BQD_EXP lo = mtbdd_refs_sync(SYNC(bqd_exp_translate));
    mtbdd_refs_pop(1);
    const BQD_EXP res = ((t >> v) & 1) ? bqd_exp_node(v, hi, lo) : bqd_exp_node(v, lo, hi);
    cache_put3(CACHE_BQD_EXP_TRANSLATE, e, t, 0, res);
    return res;
}

/* --- the two shapes the gates need ------------------------------------------- */

/*
 * From the last qubit of vars up: at each the low child is 0, since the
 * monomial is 0 where x_v is, and the high child the monomial of the qubits
 * below it. The partial result is protected across the next node.
 */
BQD_EXP
bqd_exp_monomial(uint64_t vars, int64_t gamma, uint32_t r)
{
    BQD_EXP m = mtbdd_refs_push(bqd_exp_const(gamma, r));
    if (reduce(gamma, r) != 0) {
        const BQD_EXP zero = mtbdd_refs_push(bqd_exp_const(0, r));
        for (uint64_t rest = vars; rest != 0; rest &= ~(UINT64_C(1) << (63 - __builtin_clzll(rest)))) {
            const uint32_t v = 63 - (uint32_t)__builtin_clzll(rest);     /* the lowest qubit left */
            const BQD_EXP next = bqd_exp_node(v, zero, m);
            mtbdd_refs_pop(2);
            m = mtbdd_refs_push(next);
            mtbdd_refs_push(zero);
        }
        mtbdd_refs_pop(1);
    }
    mtbdd_refs_pop(1);
    return m;
}

/*
 * The parity of the qubits of s below v and its complement, 4 - it, built
 * together from the last qubit up: at a qubit of s the low child is the
 * parity of the rest and the high child its complement.
 */
BQD_EXP
bqd_exp_parity4(uint64_t s)
{
    BQD_EXP p = mtbdd_refs_push(bqd_exp_const(0, 8));
    BQD_EXP q = mtbdd_refs_push(bqd_exp_const(4, 8));
    for (uint64_t rest = s; rest != 0; rest &= ~(UINT64_C(1) << (63 - __builtin_clzll(rest)))) {
        const uint32_t v = 63 - (uint32_t)__builtin_clzll(rest);
        const BQD_EXP np = mtbdd_refs_push(bqd_exp_node(v, p, q));
        const BQD_EXP nq = bqd_exp_node(v, q, p);
        mtbdd_refs_pop(3);
        p = mtbdd_refs_push(np);
        q = mtbdd_refs_push(nq);
    }
    mtbdd_refs_pop(2);
    return p;
}
