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
 * Tests for limdd_stab_of_node.
 *
 * The group a node has is a fact about the state the node denotes, so that is
 * what it is checked against. For each diagram the amplitude vector is read out
 * with limdd_eval, every signed Pauli word is applied to it, and the ones that
 * fix it are collected: that IS the stabiliser group, by definition, and it
 * shares nothing with the code under test -- no reduction, no intersection, no
 * case analysis. A slip in the derivation behind limdd_stab_of_node cannot hide
 * from it.
 *
 * Brute force costs 2 * 4^n vector comparisons, so the qubit count stays small.
 * The structures generated are not: random diagrams deliberately reuse children
 * so that the anti-diagonal case, the one that needs both branches swapped, is
 * exercised rather than assumed unreachable.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lace.h>
#include <sylvan.h>

#include "qsylvan_limdd_stab.h"
#include "sylvan_edge_weights_complex.h"
#include "test_assert.h"

#define NQUBITS 4
#define NBASIS  (1u << NQUBITS)

static uint64_t rng_state = UINT64_C(0x5DEECE66D1234567);
static uint64_t
rnd(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

/* --- the amplitude vector of a diagram ----------------------------------- */

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

/**
 * The amplitude of basis state `b` under `e`, where `e` sits at `level` and bit
 * k of `b` is qubit k.
 *
 * limdd_eval only starts at the root, but every node in a diagram has a state
 * of its own and each one is worth checking, so this walks from an arbitrary
 * level. It repeats the amplitude formula rather than calling into the library,
 * which is the point: the oracle should not share code with what it judges.
 */
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
    cx acc = cx_mul((cx){ w.r, w.i }, ipow[k]);

    const LIMDD_TARG t = limdd_target(e);
    if (t == LIMDD_TERMINAL) return acc;

    /* The edge may skip levels down to the target's own; it may never climb.
     * A climbing edge is malformed, and NaN makes any comparison fail. */
    const uint32_t var = limdd_node_var(t);
    if (var < level) { cx bad = {NAN, NAN}; return bad; }
    const LIMDD child = ((c >> var) & 1) ? limdd_node_high(t) : limdd_node_low(t);
    return cx_mul(acc, eval_at(child, c, var + 1));
}

/**
 * The state of node `t`, which lives at `var`, as a vector over the 2^(n-var)
 * basis states of qubits var..n-1. Index j holds qubit var+i in bit i, so the
 * vector and the shifted Pauli words below line up.
 */
static void
state_of_node(LIMDD_TARG t, uint32_t var, cx *vec)
{
    const unsigned nb = 1u << (NQUBITS - var);
    for (unsigned j = 0; j < nb; j++) {
        vec[j] = eval_at(limdd_bundle(LIMDD_LIM_IDENTITY, t),
                         (uint64_t)j << var, var);
    }
}

/* --- brute force stabiliser group ---------------------------------------- */

/**
 * Apply the signed Pauli word (x, z, neg) to `in`, writing `out`.
 *
 * P maps |c> to a phase times |c XOR x>, so the amplitude arriving at b comes
 * from b XOR x; the phase is i per Y and -1 per Z meeting a set bit.
 */
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

/** Every signed Pauli word that fixes `vec`, as (x, z, neg) packed triples. */
static uint64_t *
brute_stab(const cx *vec, unsigned nq, size_t *count)
{
    const unsigned nb = 1u << nq;
    uint64_t *found = malloc(sizeof(uint64_t) * 2 * (size_t)nb * nb);
    if (found == NULL) { fprintf(stderr, "oom\n"); exit(1); }

    size_t n = 0;
    cx img[NBASIS];
    for (uint64_t x = 0; x < nb; x++) {
        for (uint64_t z = 0; z < nb; z++) {
            for (int neg = 0; neg < 2; neg++) {
                apply_pauli(vec, img, x, z, neg, nb);
                bool fixes = true;
                for (unsigned b = 0; b < nb && fixes; b++) {
                    if (!cx_eq(img[b], vec[b])) fixes = false;
                }
                if (fixes) found[n++] = (x << 33) | (z << 1) | (uint64_t)neg;
            }
        }
    }
    *count = n;
    return found;
}

static int
cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

/** The elements of `s`, in the same packed form brute_stab produces. */
static uint64_t *
stab_elements(LIMDD_STAB s, uint32_t var, size_t *count)
{
    const size_t k = limdd_stab_ngens(s);
    const size_t n = (size_t)1 << k;
    uint64_t *out = malloc(sizeof(uint64_t) * n);
    if (out == NULL) { fprintf(stderr, "oom\n"); exit(1); }

    for (size_t i = 0; i < n; i++) {
        const LIMDD_LIM e = limdd_stab_element(s, i);
        const limdd_pauli_t p = limdd_lim_pauli(e);
        const bool neg = (limdd_lim_weight(e) != EVBDD_ONE);
        out[i] = ((p.x >> var) << 33) | ((p.z >> var) << 1) | (uint64_t)neg;
    }
    qsort(out, n, sizeof(uint64_t), cmp_u64);
    *count = n;
    return out;
}

/* --- recursive construction of a diagram's groups ------------------------ */

/**
 * Stab of node `t`, computed bottom-up, checking EVERY node against brute
 * force on the way back up.
 *
 * Checking only the root would say a group is wrong without saying where; the
 * deepest node to fail is the one with the bug, since its parents were handed
 * a wrong group to start from. `*bad` names it.
 */
static LIMDD_STAB
stab_checked(LIMDD_TARG t, int *bad, const char *what)
{
    if (t == LIMDD_TERMINAL) return LIMDD_STAB_TRIVIAL;

    const uint32_t var = limdd_node_var(t);
    const LIMDD low = limdd_node_low(t);
    const LIMDD high = limdd_node_high(t);
    /* limdd_stab_of_node wants each child's group as read at var+1, i.e.
     * extended over any levels the child edge skips. */
    const LIMDD_STAB s0 = limdd_edge_is_zero(low)  ? LIMDD_STAB_TRIVIAL
                          : limdd_stab_extend_skipped(
                                stab_checked(limdd_target(low), bad, what),
                                var + 1, limdd_level(limdd_target(low)));
    const LIMDD_STAB s1 = limdd_edge_is_zero(high) ? LIMDD_STAB_TRIVIAL
                          : limdd_stab_extend_skipped(
                                stab_checked(limdd_target(high), bad, what),
                                var + 1, limdd_level(limdd_target(high)));
    if (*bad) return LIMDD_STAB_TRIVIAL;

    const LIMDD_STAB s = limdd_stab_of_node(var, low, high, s0, s1);

    cx vec[NBASIS];
    state_of_node(t, var, vec);

    size_t nb, nc;
    uint64_t *want = brute_stab(vec, NQUBITS - var, &nb);
    qsort(want, nb, sizeof(uint64_t), cmp_u64);
    uint64_t *got = stab_elements(s, var, &nc);

    if (nb != nc || memcmp(want, got, sizeof(uint64_t) * nb) != 0) {
        fprintf(stderr,
                "%s: node %llu at qubit %u -- brute force found %zu elements, "
                "limdd_stab_of_node %zu\n",
                what, (unsigned long long)t, var, nb, nc);
        fprintf(stderr, "  low  = "); limdd_lim_fprint(stderr, limdd_label(low));
        fprintf(stderr, " -> %llu\n", (unsigned long long)limdd_target(low));
        fprintf(stderr, "  high = "); limdd_lim_fprint(stderr, limdd_label(high));
        fprintf(stderr, " -> %llu\n", (unsigned long long)limdd_target(high));
        fprintf(stderr, "  Stab(low child):\n");  limdd_stab_fprint(stderr, s0, "    ");
        fprintf(stderr, "  Stab(high child):\n"); limdd_stab_fprint(stderr, s1, "    ");
        fprintf(stderr, "  computed:\n");         limdd_stab_fprint(stderr, s, "    ");
        *bad = 1;
    }
    free(want);
    free(got);
    return s;
}

static int
check_node(LIMDD_TARG t, const char *what)
{
    int bad = 0;
    stab_checked(t, &bad, what);
    return bad;
}

/* --- hand-built states --------------------------------------------------- */

static limdd_pauli_t
pauli(const char *s)
{
    limdd_pauli_t p;
    if (!limdd_pauli_from_string(s, NQUBITS, &p)) {
        fprintf(stderr, "bad Pauli string %s\n", s);
        exit(1);
    }
    return p;
}

/** |0> on every qubit from `var` down. */
static LIMDD_TARG
ket0_from(uint32_t var)
{
    LIMDD_TARG t = LIMDD_TERMINAL;
    for (int q = NQUBITS - 1; q >= (int)var; q--) {
        t = limdd_makenode((uint32_t)q, limdd_bundle(LIMDD_LIM_IDENTITY, t),
                           limdd_zero_edge());
    }
    return t;
}

int
test_known_states(void)
{
    /* |0000>: stabilised by all 16 products of Zs. */
    if (check_node(ket0_from(0), "|0000>")) return 1;

    /*
     * The Bell state on qubits 0,1 over |00>, built the LIMDD way: both edges
     * of the top node point at the SAME node, with X on the high edge. This is
     * the case the anti-diagonal branch exists for -- without it the computed
     * group would be half the size.
     */
    LIMDD_TARG rest = ket0_from(1);
    LIMDD_TARG bell = limdd_makenode(0, limdd_bundle(LIMDD_LIM_IDENTITY, rest),
                                     limdd_bundle(limdd_lim_make(pauli("IXII"), EVBDD_ONE),
                                                  rest));
    if (check_node(bell, "bell")) return 1;

    /* |1> at the top: the low-is-zero branch, where -Z belongs at this level. */
    LIMDD_TARG one = limdd_makenode(0, limdd_zero_edge(),
                                    limdd_bundle(LIMDD_LIM_IDENTITY, rest));
    if (check_node(one, "|1000>")) return 1;

    /* A high edge with scalar i, so the anti-diagonal case picks Y not X. */
    LIMDD_TARG iphase = limdd_makenode(0, limdd_bundle(LIMDD_LIM_IDENTITY, rest),
                                       limdd_bundle(limdd_lim_make(pauli("IXII"),
                                                    complex_lookup(0.0, 1.0)), rest));
    if (check_node(iphase, "(|0>+i|1>)/bell-like")) return 1;

    /* A high edge with a scalar that squares to neither +1 nor -1: no
     * anti-diagonal element exists, and claiming one would be wrong. */
    LIMDD_TARG half = limdd_makenode(0, limdd_bundle(LIMDD_LIM_IDENTITY, rest),
                                     limdd_bundle(limdd_lim_make(pauli("IXII"),
                                                  complex_lookup(0.5, 0.0)), rest));
    if (check_node(half, "unequal branch weights")) return 1;
    return 0;
}

/* --- random diagrams ----------------------------------------------------- */

/**
 * True iff the nodes `a` and `b` denote states related by a LIM.
 *
 * limdd_stab_of_node assumes canonical children, and canonical children means
 * isomorphic ones have already been merged into a single node -- that is what
 * makes comparing node indices a valid test for isomorphism. makenode does not
 * merge yet (canonicalisation is built ON this function), so the generator
 * below has to avoid producing two distinct-but-isomorphic children itself, or
 * it would be testing against a precondition it has broken.
 *
 * Brute force again: try every Pauli word and ask whether it carries one state
 * onto a multiple of the other.
 */
static bool
iso(LIMDD_TARG a, LIMDD_TARG b, uint32_t var)
{
    if (a == b) return true;

    const unsigned nb = 1u << (NQUBITS - var);
    cx va[NBASIS], vb[NBASIS], img[NBASIS];
    state_of_node(a, var, va);
    state_of_node(b, var, vb);

    for (uint64_t x = 0; x < nb; x++) {
        for (uint64_t z = 0; z < nb; z++) {
            apply_pauli(vb, img, x, z, false, nb);

            /* Proportional? Fix the ratio on the first entry that is nonzero
             * in either vector, then require it to hold everywhere. */
            bool ok = true;
            bool have = false;
            cx ratio = { 0, 0 };
            for (unsigned j = 0; j < nb && ok; j++) {
                const bool za = (fabs(va[j].re) < 1e-9 && fabs(va[j].im) < 1e-9);
                const bool zi = (fabs(img[j].re) < 1e-9 && fabs(img[j].im) < 1e-9);
                if (za != zi) { ok = false; break; }
                if (za) continue;
                const double d = img[j].re * img[j].re + img[j].im * img[j].im;
                const cx r = { (va[j].re * img[j].re + va[j].im * img[j].im) / d,
                               (va[j].im * img[j].re - va[j].re * img[j].im) / d };
                if (!have) { ratio = r; have = true; }
                else if (!cx_eq(ratio, r)) ok = false;
            }
            if (ok && have) return true;
        }
    }
    return false;
}


/**
 * A random quasi-reduced diagram.
 *
 * Children are reused with high probability, because that is what puts the two
 * edges of a node on the same target and so reaches the anti-diagonal case;
 * drawing both children independently would leave it almost untested.
 */
static LIMDD_TARG
random_node(uint32_t var)
{
    if (var == NQUBITS) return LIMDD_TERMINAL;

    LIMDD_TARG a = random_node(var + 1);

    /*
     * Either the same child -- which is what reaches the anti-diagonal case --
     * or one that is genuinely not isomorphic to it. Anything in between would
     * violate the canonical-children precondition; see iso().
     */
    LIMDD_TARG b = a;
    if ((rnd() & 3) == 0) {
        for (int attempt = 0; attempt < 8; attempt++) {
            LIMDD_TARG cand = random_node(var + 1);
            if (!iso(a, cand, var + 1)) { b = cand; break; }
        }
    }

    const uint64_t above = (var + 1 >= 64) ? 0
                         : ~((UINT64_C(1) << (var + 1)) - 1);
    const uint64_t mask = ((UINT64_C(1) << NQUBITS) - 1) & above;

    limdd_pauli_t p;
    p.x = rnd() & mask;
    p.z = rnd() & mask;

    /* Scalars from the fourth roots of unity, plus 1/2 so that the case with
     * no anti-diagonal element occurs too. */
    static const double re[5] = { 1.0, 0.0, -1.0,  0.0, 0.5 };
    static const double im[5] = { 0.0, 1.0,  0.0, -1.0, 0.0 };
    const unsigned w = rnd() % 5;
    const EVBDD_WGT scalar = complex_lookup(re[w], im[w]);

    LIMDD low = limdd_bundle(LIMDD_LIM_IDENTITY, a);
    LIMDD high = limdd_bundle(limdd_lim_make(p, scalar), b);

    /* Occasionally kill a branch, to reach the two zero-edge cases. */
    const unsigned dice = rnd() % 16;
    if (dice == 0) high = limdd_zero_edge();
    else if (dice == 1) low = limdd_zero_edge();

    return limdd_makenode(var, low, high);
}

int
test_random_diagrams(void)
{
    for (int trial = 0; trial < 400; trial++) {
        char label[64];
        snprintf(label, sizeof(label), "random diagram %d", trial);
        if (check_node(random_node(0), label)) return 1;
    }
    return 0;
}

TASK_0(int, runtests)
{
    limdd_nodes_init(NQUBITS, 1LL << 18, 1LL << 16, 1LL << 18, 1LL << 18);

    if (test_known_states()) return 1;
    printf("stab of hand-built nodes matches brute force: ok\n");
    if (test_random_diagrams()) return 1;
    printf("stab of 400 random diagrams matches brute force: ok\n");

    printf("(%zu nodes, %zu generator cells)\n",
           limdd_node_table_count(), limdd_stab_table_count());

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
