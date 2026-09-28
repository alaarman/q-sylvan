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

#include <stdlib.h>

#include "qsylvan_evdd_ops.h"
#include "sylvan_evbdd_int.h"

static inline EVBDD
zero_edge(void)
{
    return evbdd_bundle(EVBDD_TERMINAL, EVBDD_ZERO);
}

static inline bool
is_zero(EVBDD e)
{
    return EVBDD_WEIGHT(e) == EVBDD_ZERO;
}

/** The variable of e's node, or UINT32_MAX at the terminal. */
static inline BDDVAR
var_of(EVBDD e)
{
    if (EVBDD_TARGET(e) == EVBDD_TERMINAL) return UINT32_MAX;
    return evbddnode_getvar(EVBDD_GETNODE(EVBDD_TARGET(e)));
}

EVBDD
evbdd_scale(EVBDD e, EVBDD_WGT c)
{
    if (is_zero(e) || c == EVBDD_ZERO) return zero_edge();
    if (c == EVBDD_ONE) return e;
    return evbdd_bundle(EVBDD_TARGET(e), wgt_mul(c, EVBDD_WEIGHT(e)));
}

EVBDD
evbdd_negate(EVBDD e)
{
    return evbdd_scale(e, EVBDD_MIN_ONE);
}

/**
 * f.g. The root weights factor out, so the memo is on the two nodes, in a
 * fixed order since the product commutes. Each cofactor is taken at the
 * higher of the two top variables, with evbdd_get_topvar reading a skipped
 * variable as a node whose two children are the same.
 */
TASK_IMPL_2(EVBDD, evbdd_times, EVBDD, a, EVBDD, b)
{
    if (is_zero(a) || is_zero(b)) return zero_edge();
    sylvan_gc_test();

    const EVBDD_WGT w = wgt_mul(EVBDD_WEIGHT(a), EVBDD_WEIGHT(b));
    EVBDD na = evbdd_bundle(EVBDD_TARGET(a), EVBDD_ONE);
    EVBDD nb = evbdd_bundle(EVBDD_TARGET(b), EVBDD_ONE);
    if (EVBDD_TARGET(a) == EVBDD_TERMINAL && EVBDD_TARGET(b) == EVBDD_TERMINAL)
        return evbdd_bundle(EVBDD_TERMINAL, w);
    if (na > nb) { const EVBDD t = na; na = nb; nb = t; }

    EVBDD res;
    if (cache_get3(CACHE_EVBDD_TIMES, sylvan_false, na, nb, &res)) return evbdd_scale(res, w);

    const BDDVAR va = var_of(na), vb = var_of(nb);
    BDDVAR top;
    EVBDD la, ha, lb, hb;
    evbdd_get_topvar(na, vb, &top, &la, &ha);
    evbdd_get_topvar(nb, va, &top, &lb, &hb);

    evbdd_refs_push(la); evbdd_refs_push(ha);
    evbdd_refs_push(lb); evbdd_refs_push(hb);
    evbdd_refs_spawn(SPAWN(evbdd_times, ha, hb));
    const EVBDD lo = evbdd_refs_push(CALL(evbdd_times, la, lb));
    const EVBDD hi = evbdd_refs_sync(SYNC(evbdd_times));
    evbdd_refs_pop(5);

    res = evbdd_makenode(top, lo, hi);
    cache_put3(CACHE_EVBDD_TIMES, sylvan_false, na, nb, res);
    return evbdd_scale(res, w);
}

/**
 * f|_{x_q = b}. Linear, so the memo is on the node. A node below q, or the
 * terminal, does not depend on x_q and is its own restriction; at q the
 * restriction is the child, which then skips q; above q both children are
 * restricted and joined again.
 */
TASK_IMPL_3(EVBDD, evbdd_restrict, EVBDD, e, BDDVAR, q, int, b)
{
    if (is_zero(e)) return e;
    const BDDVAR v = var_of(e);
    if (v == UINT32_MAX || v > q) return e;
    sylvan_gc_test();

    const EVBDD n = evbdd_bundle(EVBDD_TARGET(e), EVBDD_ONE);
    const uint64_t key = ((uint64_t)q << 1) | (uint64_t)(b & 1);
    EVBDD res;
    if (cache_get3(CACHE_EVBDD_RESTRICT, sylvan_false, n, key, &res))
        return evbdd_scale(res, EVBDD_WEIGHT(e));

    BDDVAR top;
    EVBDD lo, hi;
    evbdd_get_topvar(n, v, &top, &lo, &hi);
    if (v == q) {
        res = b ? hi : lo;
    } else {
        evbdd_refs_push(lo); evbdd_refs_push(hi);
        evbdd_refs_spawn(SPAWN(evbdd_restrict, hi, q, b));
        const EVBDD l = evbdd_refs_push(CALL(evbdd_restrict, lo, q, b));
        const EVBDD h = evbdd_refs_sync(SYNC(evbdd_restrict));
        evbdd_refs_pop(3);
        res = evbdd_makenode(v, l, h);
    }
    cache_put3(CACHE_EVBDD_RESTRICT, sylvan_false, n, key, res);
    return evbdd_scale(res, EVBDD_WEIGHT(e));
}

/**
 * f . [x_q = b]. A node below q, or the terminal, is the same function on
 * both sides of x_q, so it becomes a node at q with it on one side and zero
 * on the other; at q the other child is dropped; above q both are projected.
 */
TASK_IMPL_3(EVBDD, evbdd_project, EVBDD, e, BDDVAR, q, int, b)
{
    if (is_zero(e)) return e;
    const BDDVAR v = var_of(e);
    const EVBDD n = evbdd_bundle(EVBDD_TARGET(e), EVBDD_ONE);
    if (v == UINT32_MAX || v > q) {
        const EVBDD r = b ? evbdd_makenode(q, zero_edge(), n) : evbdd_makenode(q, n, zero_edge());
        return evbdd_scale(r, EVBDD_WEIGHT(e));
    }
    sylvan_gc_test();

    const uint64_t key = ((uint64_t)q << 1) | (uint64_t)(b & 1);
    EVBDD res;
    if (cache_get3(CACHE_EVBDD_PROJECT, sylvan_false, n, key, &res))
        return evbdd_scale(res, EVBDD_WEIGHT(e));

    BDDVAR top;
    EVBDD lo, hi;
    evbdd_get_topvar(n, v, &top, &lo, &hi);
    if (v == q) {
        res = b ? evbdd_makenode(q, zero_edge(), hi) : evbdd_makenode(q, lo, zero_edge());
    } else {
        evbdd_refs_push(lo); evbdd_refs_push(hi);
        evbdd_refs_spawn(SPAWN(evbdd_project, hi, q, b));
        const EVBDD l = evbdd_refs_push(CALL(evbdd_project, lo, q, b));
        const EVBDD h = evbdd_refs_sync(SYNC(evbdd_project));
        evbdd_refs_pop(3);
        res = evbdd_makenode(v, l, h);
    }
    cache_put3(CACHE_EVBDD_PROJECT, sylvan_false, n, key, res);
    return evbdd_scale(res, EVBDD_WEIGHT(e));
}

EVBDD
evbdd_local_matvec(EVBDD e, const EVBDD_WGT *M, const uint32_t *qubits, uint32_t k,
                   uint32_t n)
{
    (void)n;
    const uint32_t dim = 1u << k;
    EVBDD *vc = malloc(dim * sizeof(EVBDD));
    if (vc == NULL) { fprintf(stderr, "sylvan: out of memory in evbdd_local_matvec\n"); exit(1); }

    /* v|_{x_Q = c} for every c, each held on the reference stack, since
     * every operation below may join a collection */
    for (uint32_t c = 0; c < dim; c++) {
        EVBDD r = e;
        for (uint32_t i = 0; i < k; i++) r = evbdd_restrict(r, qubits[i], (int)((c >> (k - 1 - i)) & 1));
        vc[c] = evbdd_refs_push(r);
    }
    EVBDD out = evbdd_refs_push(zero_edge());
    for (uint32_t r = 0; r < dim; r++) {
        EVBDD w = zero_edge();
        for (uint32_t c = 0; c < dim; c++) {
            const EVBDD_WGT m = M[(size_t)r * dim + c];
            if (m == EVBDD_ZERO) continue;
            evbdd_refs_push(w);
            w = evbdd_plus(w, evbdd_scale(vc[c], m));
            evbdd_refs_pop(1);
        }
        for (uint32_t i = 0; i < k; i++) {
            evbdd_refs_push(w);
            w = evbdd_project(w, qubits[i], (int)((r >> (k - 1 - i)) & 1));
            evbdd_refs_pop(1);
        }
        evbdd_refs_push(w);
        const EVBDD next = evbdd_plus(out, w);
        evbdd_refs_pop(2);                 /* w and the old out */
        out = evbdd_refs_push(next);
    }
    evbdd_refs_pop(1 + dim);
    free(vc);
    return out;
}
