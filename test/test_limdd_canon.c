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

/**
 * Tests for limdd_makeedge.
 *
 * Two things have to hold, and they pull against each other.
 *
 * SOUNDNESS. The returned edge must denote |0>(x)low + |1>(x)high. Every rule
 * in the canonical form moves something from the node to the label, and any
 * of them dropping a factor of i or a sign would leave a diagram that is
 * beautifully canonical and means the wrong thing. So every result is read
 * back out amplitude by amplitude and compared with what went in.
 *
 * CANONICITY. Inputs denoting states related by a LIM must come back pointing
 * at the SAME node. This is tested by transforming an input in each of the
 * ways the state is invariant under -- a common LIM on both children, the
 * children exchanged, a stabiliser element on a label -- and insisting the
 * target does not move.
 *
 * A canonical form can fail either way round: too coarse merges states that
 * differ, too fine leaves merges on the table. Soundness catches the first and
 * canonicity the second, so both are needed.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lace.h>
#include <sylvan.h>

#include "qsylvan_limdd_canon.h"
#include "sylvan_edge_weights_complex.h"
#include "test_assert.h"

#define NQUBITS 4
#define NBASIS  (1u << NQUBITS)

static uint64_t rng_state = UINT64_C(0xB5026F5AA96619E9);
static uint64_t
rnd(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

/* --- reading a diagram back out ------------------------------------------ */

typedef struct { double re, im; } cx;

static cx
cx_mul(cx a, cx b)
{
    cx r = { a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re };
    return r;
}

static bool
cx_eq(cx a, cx b)
{
    return fabs(a.re - b.re) < 1e-9 && fabs(a.im - b.im) < 1e-9;
}

/** Amplitude of `b` under `e`, which sits at `level`. Bit k of `b` is qubit k. */
static cx
eval_at(LIMDD e, uint64_t b, uint32_t level)
{
    static const cx ipow[4] = { {1,0}, {0,1}, {-1,0}, {0,-1} };
    if (limdd_edge_is_zero(e)) { cx z = {0,0}; return z; }

    const LIMDD_LIM lim = limdd_label(e);
    const limdd_pauli_t p = limdd_lim_pauli(lim);
    const uint64_t c = b ^ p.x;
    const unsigned k = ((unsigned)__builtin_popcountll(p.x & p.z)
                        + 2u * (unsigned)__builtin_popcountll(p.z & c)) & 3u;

    complex_t w;
    weight_value(limdd_lim_weight(lim), &w);
    const cx acc = cx_mul((cx){ w.r, w.i }, ipow[k]);

    const LIMDD_TARG t = limdd_target(e);
    if (t == LIMDD_TERMINAL) return acc;

    const LIMDD child = ((c >> level) & 1) ? limdd_node_high(t) : limdd_node_low(t);
    return cx_mul(acc, eval_at(child, c, level + 1));
}

/** The vector `e` denotes, over the 2^(n-level) basis states at and below it. */
static void
vector_of(LIMDD e, uint32_t level, cx *out)
{
    const unsigned n = 1u << (NQUBITS - level);
    for (unsigned j = 0; j < n; j++) out[j] = eval_at(e, (uint64_t)j << level, level);
}

/* --- building random inputs ---------------------------------------------- */

static limdd_pauli_t
random_word_above(uint32_t var)
{
    const uint64_t all = (UINT64_C(1) << NQUBITS) - 1;
    const uint64_t above = (var + 1 >= 64) ? 0 : ~((UINT64_C(1) << (var + 1)) - 1);
    const uint64_t mask = all & above;
    limdd_pauli_t p;
    p.x = rnd() & mask;
    p.z = rnd() & mask;
    return p;
}

static EVBDD_WGT
random_scalar(void)
{
    static const double re[6] = { 1.0, 0.0, -1.0,  0.0, 0.5, -0.25 };
    static const double im[6] = { 0.0, 1.0,  0.0, -1.0, 0.0,  0.5  };
    const unsigned w = rnd() % 6;
    return complex_lookup(re[w], im[w]);
}

/** `e` with its label multiplied on the left by `l`. */
static LIMDD
relabel(LIMDD e, LIMDD_LIM l)
{
    if (limdd_edge_is_zero(e)) return e;
    return limdd_bundle(limdd_lim_mul(l, limdd_label(e)), limdd_target(e));
}

/**
 * A random canonical edge at `var`, built with makeedge all the way down.
 *
 * Children are reused most of the time, and when they are not the second is
 * usually a relabelling of the first: both produce the isomorphic children
 * that the swap rule and the coset minimisation exist for. Drawing two
 * unrelated children would leave those paths mostly cold.
 */
static LIMDD
random_edge(uint32_t var)
{
    if (var == NQUBITS) return limdd_one_edge();

    LIMDD lo = random_edge(var + 1);
    LIMDD hi;
    const unsigned dice = rnd() % 8;
    if (dice < 3)      hi = lo;
    else if (dice < 6) hi = relabel(lo, limdd_lim_make(random_word_above(var), random_scalar()));
    else               hi = random_edge(var + 1);

    if ((rnd() % 12) == 0) hi = limdd_zero_edge();
    else if ((rnd() % 12) == 0) lo = limdd_zero_edge();

    return limdd_makeedge(var, lo, hi);
}

/* --- soundness ------------------------------------------------------------ */

static int
check_sound(uint32_t var, LIMDD lo, LIMDD hi, LIMDD got, const char *what)
{
    cx vlo[NBASIS], vhi[NBASIS], vgot[NBASIS];
    vector_of(lo, var + 1, vlo);
    vector_of(hi, var + 1, vhi);
    vector_of(got, var, vgot);

    const unsigned n = 1u << (NQUBITS - var);
    for (unsigned j = 0; j < n; j++) {
        const cx want = (j & 1) ? vhi[j >> 1] : vlo[j >> 1];
        if (!cx_eq(vgot[j], want)) {
            fprintf(stderr,
                    "%s: makeedge at qubit %u changed the state at index %u: "
                    "got (%g,%g), expected (%g,%g)\n",
                    what, var, j, vgot[j].re, vgot[j].im, want.re, want.im);
            return 1;
        }
    }
    return 0;
}

int
test_soundness(void)
{
    for (int trial = 0; trial < 3000; trial++) {
        const uint32_t var = rnd() % NQUBITS;
        LIMDD lo = random_edge(var + 1);
        LIMDD hi = ((rnd() & 1)) ? lo : random_edge(var + 1);
        if ((rnd() % 10) == 0) hi = limdd_zero_edge();
        if (limdd_edge_is_zero(lo) && limdd_edge_is_zero(hi)) continue;

        char label[48];
        snprintf(label, sizeof(label), "trial %d", trial);
        if (check_sound(var, lo, hi, limdd_makeedge(var, lo, hi), label)) return 1;
    }
    return 0;
}

/* --- canonicity ----------------------------------------------------------- */

/** Report a canonicity failure with enough context to act on. */
static int
differ(const char *what, uint32_t var, LIMDD a, LIMDD b)
{
    fprintf(stderr, "%s at qubit %u: node %llu vs %llu -- LIM-related states "
                    "must share a node\n", what, var,
            (unsigned long long)limdd_target(a), (unsigned long long)limdd_target(b));
    return 1;
}

int
test_canonical_under_lim(void)
{
    /*
     * Multiplying both children by one LIM multiplies the whole state by it,
     * so the node must not move. This is what the low-factoring and the coset
     * minimisation are for, together.
     */
    for (int trial = 0; trial < 2000; trial++) {
        const uint32_t var = rnd() % NQUBITS;
        LIMDD lo = random_edge(var + 1);
        LIMDD hi = (rnd() & 1) ? lo : random_edge(var + 1);
        if (limdd_edge_is_zero(lo) && limdd_edge_is_zero(hi)) continue;

        const LIMDD_LIM g = limdd_lim_make(random_word_above(var), random_scalar());
        const LIMDD ref = limdd_makeedge(var, lo, hi);
        const LIMDD alt = limdd_makeedge(var, relabel(lo, g), relabel(hi, g));

        if (limdd_target(ref) != limdd_target(alt))
            return differ("common LIM on both children", var, ref, alt);

        /* And it really is the same state up to a LIM: soundness of the alt. */
        if (check_sound(var, relabel(lo, g), relabel(hi, g), alt, "relabelled")) return 1;
    }
    return 0;
}

int
test_canonical_under_swap(void)
{
    /*
     * X at this level exchanges the branches, so (low, high) and (high, low)
     * are the same node. Without the swap rule these stay distinct and every
     * merge it would have bought is lost.
     */
    for (int trial = 0; trial < 2000; trial++) {
        const uint32_t var = rnd() % NQUBITS;
        LIMDD lo = random_edge(var + 1);
        LIMDD hi = (rnd() & 1) ? lo : random_edge(var + 1);
        if (limdd_edge_is_zero(lo) && limdd_edge_is_zero(hi)) continue;

        const LIMDD ref = limdd_makeedge(var, lo, hi);
        const LIMDD alt = limdd_makeedge(var, hi, lo);
        if (limdd_target(ref) != limdd_target(alt))
            return differ("children exchanged", var, ref, alt);
        if (check_sound(var, hi, lo, alt, "swapped")) return 1;
    }
    return 0;
}

int
test_canonical_under_stabiliser(void)
{
    /*
     * A stabiliser element on a child's label does not change that child's
     * state at all, so it certainly must not change the node. This is the
     * freedom the coset minimisation quotients out.
     */
    for (int trial = 0; trial < 2000; trial++) {
        const uint32_t var = rnd() % NQUBITS;
        LIMDD lo = random_edge(var + 1);
        LIMDD hi = random_edge(var + 1);
        if (limdd_edge_is_zero(lo) || limdd_edge_is_zero(hi)) continue;

        const LIMDD_STAB sh = limdd_edge_stab(hi);
        const size_t k = limdd_stab_ngens(sh);
        if (k == 0) continue;
        const LIMDD_LIM h = limdd_stab_element(sh, rnd() & (((uint64_t)1 << k) - 1));

        const LIMDD ref = limdd_makeedge(var, lo, hi);
        const LIMDD alt = limdd_makeedge(var, lo, relabel(hi, h));
        if (limdd_target(ref) != limdd_target(alt))
            return differ("stabiliser on the high label", var, ref, alt);
    }
    return 0;
}

int
test_zero_branch_merges(void)
{
    /*
     * |0>(x)v and |1>(x)v differ by X, so they are one node. The visible
     * consequence is that |0000> and |1000> share their whole diagram.
     */
    for (int trial = 0; trial < 500; trial++) {
        const uint32_t var = rnd() % NQUBITS;
        const LIMDD v = random_edge(var + 1);
        if (limdd_edge_is_zero(v)) continue;

        const LIMDD a = limdd_makeedge(var, v, limdd_zero_edge());
        const LIMDD b = limdd_makeedge(var, limdd_zero_edge(), v);
        if (limdd_target(a) != limdd_target(b))
            return differ("dead branch on either side", var, a, b);
        if (check_sound(var, limdd_zero_edge(), v, b, "dead low branch")) return 1;
    }

    /* |0000> and |1000>, end to end. */
    LIMDD zero = limdd_one_edge(), one = limdd_one_edge();
    for (int q = NQUBITS - 1; q >= 0; q--) {
        zero = limdd_makeedge((uint32_t)q, zero, limdd_zero_edge());
        one  = limdd_makeedge((uint32_t)q, limdd_zero_edge(), one);
    }
    test_assert(limdd_target(zero) == limdd_target(one));
    return 0;
}

int
test_bell_shares_a_node(void)
{
    /*
     * The Bell state needs one node per level, not two: the |1> branch is the
     * |0> branch with X applied, and the canonical form is what notices.
     */
    const size_t before = limdd_node_table_count();

    LIMDD k0 = limdd_one_edge();
    for (int q = NQUBITS - 1; q >= 1; q--) k0 = limdd_makeedge((uint32_t)q, k0, limdd_zero_edge());

    LIMDD k1 = limdd_one_edge();
    for (int q = NQUBITS - 1; q >= 2; q--) k1 = limdd_makeedge((uint32_t)q, k1, limdd_zero_edge());
    k1 = limdd_makeedge(1, limdd_zero_edge(), k1);

    /* |00..> and |01..> on qubits 1.. must already be the same node. */
    test_assert(limdd_target(k0) == limdd_target(k1));

    const LIMDD bell = limdd_makeedge(0, k0, k1);
    test_assert(!limdd_edge_is_zero(bell));

    /*
     * The property, stated directly: the top node's two branches land on ONE
     * node, differing only by the label. Counting nodes would not say this --
     * earlier tests have already filled the table, so the delta can be zero
     * whatever happens here.
     */
    const LIMDD_TARG top = limdd_target(bell);
    test_assert(limdd_target(limdd_node_low(top)) == limdd_target(limdd_node_high(top)));
    test_assert(limdd_node_low(top) != limdd_node_high(top));
    test_assert(!limdd_lim_is_identity(limdd_label(limdd_node_high(top))));

    /* And it is the Bell state: amplitudes equal on |00..> and |11..>, zero
     * on the other two. */
    cx v[NBASIS];
    vector_of(bell, 0, v);
    test_assert(!cx_eq(v[0], (cx){0,0}));
    test_assert(cx_eq(v[0], v[3]) || cx_eq(v[0], (cx){-v[3].re, -v[3].im}));
    test_assert(cx_eq(v[1], (cx){0,0}));
    test_assert(cx_eq(v[2], (cx){0,0}));

    (void)before;
    return 0;
}

/* --- the cached stabiliser groups ---------------------------------------- */

static void
apply_pauli(const cx *in, cx *out, uint64_t x, uint64_t z, bool neg, unsigned nb)
{
    static const cx ipow[4] = { {1,0}, {0,1}, {-1,0}, {0,-1} };
    const unsigned ny = (unsigned)__builtin_popcountll(x & z);
    for (unsigned b = 0; b < nb; b++) {
        const unsigned c = b ^ (unsigned)x;
        const unsigned k = (ny + 2u * (unsigned)__builtin_popcountll(z & c)) & 3u;
        cx r = cx_mul(in[c], ipow[k]);
        if (neg) { r.re = -r.re; r.im = -r.im; }
        out[b] = r;
    }
}

static int
cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

/**
 * Check limdd_node_stab against the definition, for every node below `t`.
 *
 * limdd_stab_of_node only computes the right group when its children have
 * already been merged by isomorphism -- a precondition the raw makenode tests
 * had to arrange by hand, by rejecting children that were isomorphic but
 * stored apart. Diagrams built through makeedge are supposed to satisfy it for
 * free. This is where that claim is actually tested rather than assumed.
 */
static int
check_stab_below(LIMDD_TARG t, uint32_t var)
{
    if (t == LIMDD_TERMINAL) return 0;

    const LIMDD low = limdd_node_low(t);
    const LIMDD high = limdd_node_high(t);
    if (!limdd_edge_is_zero(low) && check_stab_below(limdd_target(low), var + 1)) return 1;
    if (!limdd_edge_is_zero(high) && check_stab_below(limdd_target(high), var + 1)) return 1;

    const unsigned nq = NQUBITS - var;
    const unsigned nb = 1u << nq;

    cx vec[NBASIS], img[NBASIS];
    vector_of(limdd_bundle(LIMDD_LIM_IDENTITY, t), var, vec);

    uint64_t want[2 * NBASIS * NBASIS];
    size_t nw = 0;
    for (uint64_t x = 0; x < nb; x++) {
        for (uint64_t z = 0; z < nb; z++) {
            for (int neg = 0; neg < 2; neg++) {
                apply_pauli(vec, img, x, z, neg, nb);
                bool fixes = true;
                for (unsigned b = 0; b < nb && fixes; b++) {
                    if (!cx_eq(img[b], vec[b])) fixes = false;
                }
                if (fixes) want[nw++] = (x << 33) | (z << 1) | (uint64_t)neg;
            }
        }
    }
    qsort(want, nw, sizeof(uint64_t), cmp_u64);

    const LIMDD_STAB s = limdd_node_stab(t);
    const size_t k = limdd_stab_ngens(s);
    const size_t ng = (size_t)1 << k;
    uint64_t got[1u << NQUBITS];
    for (size_t i = 0; i < ng; i++) {
        const LIMDD_LIM e = limdd_stab_element(s, i);
        const limdd_pauli_t p = limdd_lim_pauli(e);
        got[i] = ((p.x >> var) << 33) | ((p.z >> var) << 1)
               | (uint64_t)(limdd_lim_weight(e) != EVBDD_ONE);
    }
    qsort(got, ng, sizeof(uint64_t), cmp_u64);

    if (nw != ng || memcmp(want, got, sizeof(uint64_t) * nw) != 0) {
        fprintf(stderr, "node %llu at qubit %u: cached group has %zu elements, "
                        "the state is fixed by %zu\n",
                (unsigned long long)t, var, ng, nw);
        limdd_stab_fprint(stderr, s, "    ");
        return 1;
    }
    return 0;
}

int
test_stab_of_canonical_nodes(void)
{
    for (int trial = 0; trial < 300; trial++) {
        const LIMDD e = random_edge(0);
        if (limdd_edge_is_zero(e)) continue;
        if (check_stab_below(limdd_target(e), 0)) {
            fprintf(stderr, "  (trial %d)\n", trial);
            return 1;
        }
    }
    return 0;
}

/* --- concurrency ---------------------------------------------------------- */

#define NSHARED 128
static LIMDD shared_lo[NSHARED], shared_hi[NSHARED];
static LIMDD_TARG seen[NSHARED];

VOID_TASK_2(build_shared, size_t, lo, size_t, hi)
{
    if (hi - lo > 16) {
        size_t mid = lo + (hi - lo) / 2;
        SPAWN(build_shared, lo, mid);
        CALL(build_shared, mid, hi);
        SYNC(build_shared);
        return;
    }
    for (size_t i = lo; i < hi; i++) {
        const LIMDD_TARG t = limdd_target(limdd_makeedge(0, shared_lo[i], shared_hi[i]));
        if (seen[i] != 0 && seen[i] != t) {
            fprintf(stderr, "edge %zu canonicalised to %llu and also %llu\n",
                    i, (unsigned long long)seen[i], (unsigned long long)t);
            exit(1);
        }
        seen[i] = t;
    }
}

/**
 * Workers race on the same edges. Each node's stabiliser group is cached
 * lazily, so a worker can read a slot another is still filling and recompute
 * it; that is meant to be harmless, and this is what says so.
 */
int
test_concurrent(void)
{
    for (size_t i = 0; i < NSHARED; i++) {
        shared_lo[i] = random_edge(1);
        shared_hi[i] = (rnd() & 1) ? shared_lo[i] : random_edge(1);
        seen[i] = 0;
    }
    for (int round = 0; round < 8; round++) RUN(build_shared, 0, NSHARED);
    for (size_t i = 0; i < NSHARED; i++) {
        test_assert(seen[i] != 0);
        test_assert(seen[i] == limdd_target(limdd_makeedge(0, shared_lo[i], shared_hi[i])));
    }
    return 0;
}

TASK_0(int, runtests)
{
    limdd_nodes_init(NQUBITS, 1LL << 18, 1LL << 18, 1LL << 18, 1LL << 18);

    if (test_soundness()) return 1;
    printf("makeedge preserves the state (3000 cases):   ok\n");
    if (test_canonical_under_lim()) return 1;
    printf("a common LIM on both children keeps the node: ok\n");
    if (test_canonical_under_swap()) return 1;
    printf("exchanging the children keeps the node:       ok\n");
    if (test_canonical_under_stabiliser()) return 1;
    printf("a stabiliser on a label keeps the node:       ok\n");
    if (test_zero_branch_merges()) return 1;
    printf("|0000> and |1000> share every node:           ok\n");
    if (test_bell_shares_a_node()) return 1;
    printf("the Bell state needs one node per level:      ok\n");
    if (test_stab_of_canonical_nodes()) return 1;
    printf("cached groups match the states they fix:      ok\n");
    if (test_concurrent()) return 1;
    printf("makeedge agrees across workers:               ok\n");

    printf("(%zu nodes, %zu LIMs, %zu generator cells)\n",
           limdd_node_table_count(), limdd_lim_table_count(), limdd_stab_table_count());

    limdd_nodes_quit();
    return 0;
}

int
main(void)
{
    lace_start(8, 0);
    sylvan_set_sizes(1LL << 16, 1LL << 16, 1LL << 16, 1LL << 16);
    sylvan_init_package();
    sylvan_init_edge_weights(1LL << 18, 1LL << 18, 1e-14, WGT_COMPLEX_128, COMP_HASHMAP);

    int res = RUN(runtests);

    sylvan_edge_weights_free();
    sylvan_quit();
    lace_stop();
    return res;
}
