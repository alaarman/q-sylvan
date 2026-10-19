/*
 * Copyright 2026 System Verification Lab, LIACS, Leiden University
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

/*
 * The exact weights' arithmetic by index (qisq2_wgt_*: operands read in place,
 * per-thread temporaries, the inverse by its closed form) against the
 * weight_qisq2_* functions it replaced, which work on copies. Both intern
 * their result, and the table holds one entry per value, so the two must
 * return the same index. The pool is grown by the new operations themselves,
 * so it reaches long rationals, and the checks run on several workers at once.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <sylvan.h>
#include <sylvan_int.h>

#include "qsylvan_simulator.h"
#include "sylvan_edge_weights_qisq2.h"
#include "test_assert.h"

#define NPOOL 256
static EVBDD_WGT pool[NPOOL];
static uint64_t rng = UINT64_C(0x9E3779B97F4A7C15);

static uint64_t
rnd(void)
{
    rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
    return rng;
}

/* The old path: copies, weight_qisq2_<op>, intern. */
static EVBDD_WGT
ref_binary(EVBDD_WGT a, EVBDD_WGT b, void (*op)(qisq2_t *, qisq2_t *))
{
    qisq2_t x, y;
    weight_value(a, &x);
    weight_value(b, &y);
    op(&x, &y);
    const EVBDD_WGT r = weight_qisq2_lookup(&x);
    qisq2_clear(&x);
    qisq2_clear(&y);
    return r;
}

static EVBDD_WGT
ref_unary(EVBDD_WGT a, void (*op)(qisq2_t *))
{
    qisq2_t x;
    weight_value(a, &x);
    op(&x);
    const EVBDD_WGT r = weight_qisq2_lookup(&x);
    qisq2_clear(&x);
    return r;
}

static void
check(const char *what, EVBDD_WGT a, EVBDD_WGT b, EVBDD_WGT got, EVBDD_WGT want)
{
    if (got == want) return;
    fprintf(stderr, "%s of %llu and %llu: new %llu, old %llu\n", what,
            (unsigned long long)a, (unsigned long long)b,
            (unsigned long long)got, (unsigned long long)want);
    wgt_fprint(stderr, got);
    wgt_fprint(stderr, want);
    exit(1);
}

/* Every operation on the pairs [lo, hi) of a fixed list. */
VOID_TASK_3(check_range, uint64_t, seed, size_t, lo, size_t, hi)
{
    if (hi - lo > 64) {
        const size_t mid = lo + (hi - lo) / 2;
        SPAWN(check_range, seed, lo, mid);
        CALL(check_range, seed, mid, hi);
        SYNC(check_range);
        return;
    }
    for (size_t k = lo; k < hi; k++) {
        uint64_t s = seed ^ (k * UINT64_C(0xD1B54A32D192ED03));
        s ^= s >> 31;
        const EVBDD_WGT a = pool[s % NPOOL], b = pool[(s >> 20) % NPOOL];
        if (a == EVBDD_ZERO || b == EVBDD_ZERO) continue;
        check("mul", a, b, qisq2_wgt_mul(a, b), ref_binary(a, b, weight_qisq2_mul));
        check("div", a, b, qisq2_wgt_div(a, b), ref_binary(a, b, weight_qisq2_div));
        check("inv", EVBDD_ONE, b, qisq2_wgt_inv(b), ref_binary(EVBDD_ONE, b, weight_qisq2_div));
        check("add", a, b, qisq2_wgt_add(a, b, false), ref_binary(a, b, weight_qisq2_add));
        check("sub", a, b, qisq2_wgt_add(a, b, true), ref_binary(a, b, weight_qisq2_sub));
        check("neg", a, a, qisq2_wgt_neg(a), ref_unary(a, weight_qisq2_neg));
        check("conj", a, a, qisq2_wgt_conj(a), ref_unary(a, weight_qisq2_conj));
        /* and through the dispatch, which adds the cache and the special cases */
        check("wgt_mul", a, b, wgt_mul(a, b), ref_binary(a, b, weight_qisq2_mul));
        check("wgt_div", a, b, wgt_div(a, b), ref_binary(a, b, weight_qisq2_div));
        double re, im;
        qisq2_wgt_as_doubles(a, &re, &im);
        qisq2_t x;
        weight_value(a, &x);
        const double r2 = mpq_get_d(x.a) + mpq_get_d(x.b) * sqrt(2.0);
        const double i2 = mpq_get_d(x.c) + mpq_get_d(x.d) * sqrt(2.0);
        qisq2_clear(&x);
        if (re != r2 || im != i2) {
            fprintf(stderr, "as doubles of %llu: new (%.17g, %.17g), old (%.17g, %.17g)\n",
                    (unsigned long long)a, re, im, r2, i2);
            exit(1);
        }
    }
}

TASK_0(int, run)
{
    /* the units, and a few values with nontrivial rationals in every slot */
    size_t n = 0;
    pool[n++] = EVBDD_ONE;
    pool[n++] = EVBDD_MIN_ONE;
    pool[n++] = qisq2_lookup(0, 1, 0, 1, 1, 1, 0, 1);      /*  i */
    pool[n++] = qisq2_lookup(0, 1, 0, 1, -1, 1, 0, 1);     /* -i */
    pool[n++] = qisq2_lookup(0, 1, 1, 2, 0, 1, 0, 1);      /* sqrt2/2 */
    pool[n++] = qisq2_lookup(1, 2, 0, 1, 1, 2, 0, 1);      /* (1+i)/2 */
    pool[n++] = qisq2_lookup(3, 7, -2, 5, 0, 1, 4, 9);
    pool[n++] = qisq2_lookup(-1, 3, 0, 1, 5, 11, -1, 2);
    pool[n++] = qisq2_lookup(0, 1, 7, 3, -2, 13, 0, 1);
    pool[n++] = qisq2_lookup(1, 1, 1, 1, 0, 1, 0, 1);      /* 1 + sqrt2 */
    pool[n++] = EVBDD_ZERO;
    /* grown by the new operations, so long rationals appear */
    while (n < NPOOL) {
        const EVBDD_WGT a = pool[rnd() % n], b = pool[rnd() % n];
        EVBDD_WGT c;
        switch (rnd() % 4) {
        case 0:  c = (a == EVBDD_ZERO || b == EVBDD_ZERO) ? EVBDD_ZERO : qisq2_wgt_mul(a, b); break;
        case 1:  c = (a == EVBDD_ZERO || b == EVBDD_ZERO) ? EVBDD_ZERO : qisq2_wgt_div(a, b); break;
        case 2:  c = qisq2_wgt_add(a, b, false); break;
        default: c = qisq2_wgt_add(a, b, true); break;
        }
        pool[n++] = c;
    }
    size_t longest = 0;
    for (size_t i = 0; i < NPOOL; i++) {
        const size_t bits = qisq2_size(wgt_storage, pool[i]);
        if (bits > longest) longest = bits;
    }
    for (int round = 0; round < 4; round++) CALL(check_range, rnd(), 0, 20000);
    printf("qisq2 arithmetic by index equals the copying path (4 x 20000 pairs, up to %zu-bit integers): ok\n",
           longest);
    return 0;
}

int
main(void)
{
    lace_start(getenv("LL_WORKERS") ? (unsigned)atoi(getenv("LL_WORKERS")) : 4, 0);
    sylvan_set_sizes(1LL << 20, 1LL << 20, 1LL << 20, 1LL << 20);
    sylvan_init_package();
    qsylvan_init_simulator(1LL << 22, 1LL << 22, -1, QISQ2_MAP, NORM_LOW);
    const int res = RUN(run);
    sylvan_quit();
    lace_stop();
    return res;
}
