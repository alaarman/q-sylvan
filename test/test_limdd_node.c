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
 * Tests for LIMDD nodes and edges.
 *
 * Every test here is anchored on limdd_eval: a diagram is built by hand, and
 * the amplitude vector it denotes is compared against the state it is supposed
 * to be. Checking shape instead would only confirm that the code does what it
 * does; checking amplitudes is what makes eval usable as an oracle for the
 * normalisation rules later on.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lace.h>
#include <sylvan.h>

#include "qsylvan_limdd_node.h"
#include "sylvan_edge_weights_complex.h"
#include "test_assert.h"

#define NQUBITS 3

static const double SQRT2_INV = 0.70710678118654752440;

/* Amplitudes are compared as complex numbers, not as weight indices: two
 * different indices can hold the same value, and it is the value that the
 * state is. */
static bool
amp_is(EVBDD_WGT w, double re, double im)
{
    complex_t c;
    weight_value(w, &c);
    return fabs(c.r - re) < 1e-12 && fabs(c.i - im) < 1e-12;
}

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

/* --- edge encoding ------------------------------------------------------- */

int
test_edge_encoding(void)
{
    /* Round-tripping matters because the two fields share one word: a wrong
     * shift shows up as a target that silently carries label bits. */
    for (uint64_t lim = 0; lim < 64; lim++) {
        for (uint64_t t = 0; t < 64; t++) {
            LIMDD e = limdd_bundle(lim, t);
            test_assert(limdd_label(e) == lim);
            test_assert(limdd_target(e) == t);
        }
    }

    /* The extremes of both fields, where an off-by-one bit would overlap. */
    LIMDD e = limdd_bundle(LIMDD_LIM_MAX - 1, LIMDD_TARG_MAX);
    test_assert(limdd_label(e) == LIMDD_LIM_MAX - 1);
    test_assert(limdd_target(e) == LIMDD_TARG_MAX);

    test_assert(limdd_edge_is_zero(limdd_zero_edge()));
    test_assert(!limdd_edge_is_zero(limdd_one_edge()));
    test_assert(limdd_target(limdd_one_edge()) == LIMDD_TERMINAL);
    test_assert(limdd_label(limdd_one_edge()) == LIMDD_LIM_IDENTITY);
    return 0;
}

/* --- node interning ------------------------------------------------------ */

int
test_node_interning(void)
{
    LIMDD one = limdd_one_edge();
    LIMDD zero = limdd_zero_edge();

    /* Equal triples must give equal indices; that is what lets makenode share
     * subdiagrams at all. */
    LIMDD_TARG a = limdd_makenode(2, one, zero);
    LIMDD_TARG b = limdd_makenode(2, one, zero);
    test_assert(a == b);

    /* And differing in any one field must not. */
    test_assert(limdd_makenode(1, one, zero) != a);
    test_assert(limdd_makenode(2, zero, one) != a);
    test_assert(limdd_makenode(2, one, one) != a);

    test_assert(limdd_node_var(a) == 2);
    test_assert(limdd_node_low(a) == one);
    test_assert(limdd_edge_is_zero(limdd_node_high(a)));

    /*
     * A zero edge is the zero vector whatever it points at, so two nodes whose
     * dead branches dangle off different targets are the same node. If this
     * failed, the table would fill up with duplicates that never match.
     */
    LIMDD_TARG other = limdd_makenode(2, one, one);
    LIMDD zero_elsewhere = limdd_bundle(LIMDD_LIM_ZERO, other);
    test_assert(limdd_makenode(2, one, zero_elsewhere) == a);

    /*
     * Quasi-reduced: a node with identical children is a real node. Collapsing
     * it would drop a qubit from the state, which is exactly what the |+++>
     * test below would then get wrong.
     */
    LIMDD_TARG twin = limdd_makenode(2, one, one);
    test_assert(twin != LIMDD_TERMINAL);
    test_assert(limdd_node_low(twin) == limdd_node_high(twin));
    return 0;
}

/* --- evaluation ---------------------------------------------------------- */

/** The node for |0> on the last qubit: low = 1, high = 0. */
static LIMDD_TARG
make_ket0(uint32_t var)
{
    return limdd_makenode(var, limdd_one_edge(), limdd_zero_edge());
}

/** The node for |1> on the last qubit: low = 0, high = 1. */
static LIMDD_TARG
make_ket1(uint32_t var)
{
    return limdd_makenode(var, limdd_zero_edge(), limdd_one_edge());
}

static void
set_bits(bool *bits, unsigned v)
{
    for (int k = 0; k < NQUBITS; k++) bits[k] = (v >> k) & 1;
}

int
test_eval_basis_states(void)
{
    bool bits[NQUBITS];

    /* |000>: every level is |0>, so exactly one basis state has amplitude 1. */
    LIMDD e = limdd_one_edge();
    for (int q = NQUBITS - 1; q >= 0; q--) {
        e = limdd_bundle(LIMDD_LIM_IDENTITY, limdd_makenode(q, e, limdd_zero_edge()));
    }
    for (unsigned v = 0; v < (1u << NQUBITS); v++) {
        set_bits(bits, v);
        test_assert(amp_is(limdd_eval(e, bits, NQUBITS), v == 0 ? 1.0 : 0.0, 0.0));
    }

    /* |101>, built structurally: qubit 0 and 2 are |1>, qubit 1 is |0>. */
    LIMDD f = limdd_one_edge();
    f = limdd_bundle(LIMDD_LIM_IDENTITY, make_ket1(2));
    f = limdd_bundle(LIMDD_LIM_IDENTITY, limdd_makenode(1, f, limdd_zero_edge()));
    f = limdd_bundle(LIMDD_LIM_IDENTITY, limdd_makenode(0, limdd_zero_edge(), f));
    for (unsigned v = 0; v < (1u << NQUBITS); v++) {
        set_bits(bits, v);
        test_assert(amp_is(limdd_eval(f, bits, NQUBITS), v == 0x5 ? 1.0 : 0.0, 0.0));
    }
    return 0;
}

int
test_eval_lim_acts(void)
{
    bool bits[NQUBITS];

    /* |000> again, to act on. */
    LIMDD ket000 = limdd_one_edge();
    for (int q = NQUBITS - 1; q >= 0; q--) {
        ket000 = limdd_bundle(LIMDD_LIM_IDENTITY,
                              limdd_makenode(q, ket000, limdd_zero_edge()));
    }

    /*
     * X on qubit 0 must turn |000> into |100>, WITHOUT touching a single node.
     * This is the property the whole representation exists for: a label change
     * is a state change, so evaluation has to honour it.
     */
    LIMDD x0 = limdd_bundle(limdd_lim_make(pauli("XII"), EVBDD_ONE),
                            limdd_target(ket000));
    for (unsigned v = 0; v < (1u << NQUBITS); v++) {
        set_bits(bits, v);
        test_assert(amp_is(limdd_eval(x0, bits, NQUBITS), v == 0x1 ? 1.0 : 0.0, 0.0));
    }

    /* Z fixes |000>; on |0> it is the identity. */
    LIMDD z0 = limdd_bundle(limdd_lim_make(pauli("ZII"), EVBDD_ONE),
                            limdd_target(ket000));
    for (unsigned v = 0; v < (1u << NQUBITS); v++) {
        set_bits(bits, v);
        test_assert(amp_is(limdd_eval(z0, bits, NQUBITS), v == 0 ? 1.0 : 0.0, 0.0));
    }

    /* Y|0> = i|1>, so the i must survive into the amplitude. */
    LIMDD y0 = limdd_bundle(limdd_lim_make(pauli("YII"), EVBDD_ONE),
                            limdd_target(ket000));
    set_bits(bits, 0x1);
    test_assert(amp_is(limdd_eval(y0, bits, NQUBITS), 0.0, 1.0));
    set_bits(bits, 0x0);
    test_assert(amp_is(limdd_eval(y0, bits, NQUBITS), 0.0, 0.0));

    /*
     * ZX|000> = Z|100> = -|100>. Built as ONE Pauli word, with the leftover
     * power of i folded into the scalar, which is the same bookkeeping
     * limdd_lim_mul does -- so this checks that eval and composition agree on
     * where the phase lives.
     */
    limdd_pauli_t zx = pauli("ZII");
    const unsigned k = limdd_pauli_rightmul(&zx, pauli("XII"), NQUBITS);
    LIMDD zx000 = limdd_bundle(limdd_lim_make(zx, limdd_wgt_i_pow(k)),
                               limdd_target(ket000));
    set_bits(bits, 0x1);
    test_assert(amp_is(limdd_eval(zx000, bits, NQUBITS), -1.0, 0.0));
    set_bits(bits, 0x0);
    test_assert(amp_is(limdd_eval(zx000, bits, NQUBITS), 0.0, 0.0));
    return 0;
}

int
test_eval_superposition(void)
{
    bool bits[NQUBITS];

    /*
     * |+++>: every level is a node whose two children are the SAME edge --
     * exactly the node a fully-reduced diagram would collapse away. All 8
     * amplitudes must come out equal, which they cannot if the level is lost.
     *
     * The 1/sqrt(2) per qubit is not on the child edges: a low edge carries no
     * label, so each level's scalar is factored up, and the root edge ends up
     * with the whole (1/sqrt(2))^n.
     */
    LIMDD_TARG t = LIMDD_TERMINAL;
    for (int q = NQUBITS - 1; q >= 0; q--) {
        LIMDD child = limdd_bundle(LIMDD_LIM_IDENTITY, t);
        t = limdd_makenode(q, child, child);
        test_assert(limdd_node_low(t) == limdd_node_high(t));
    }

    double amp = 1.0;
    for (int q = 0; q < NQUBITS; q++) amp *= SQRT2_INV;
    LIMDD e = limdd_bundle(limdd_lim_make(limdd_pauli_identity(),
                                          complex_lookup(amp, 0.0)), t);

    for (unsigned v = 0; v < (1u << NQUBITS); v++) {
        set_bits(bits, v);
        test_assert(amp_is(limdd_eval(e, bits, NQUBITS), amp, 0.0));
    }
    return 0;
}

int
test_bell_state_two_ways(void)
{
    /*
     * The Bell state on qubits 0 and 1, with qubit 2 left in |0>, built twice:
     * once with two distinct subdiagrams, and once with ONE subdiagram whose
     * second use is relabelled with X. The second is the LIMDD; both must
     * denote the same vector, and the second must use fewer nodes.
     *
     * This is the smallest test that shows a LIM buying something, and it is
     * the case that a plain EVBDD cannot share.
     */
    bool bits[NQUBITS];
    const EVBDD_WGT inv_sqrt2 = complex_lookup(SQRT2_INV, 0.0);

    LIMDD_TARG q2_zero = make_ket0(2);

    /* --- explicit: |00> and |11> get their own qubit-1 nodes --- */
    LIMDD_TARG v0 = limdd_makenode(1, limdd_bundle(LIMDD_LIM_IDENTITY, q2_zero),
                                      limdd_zero_edge());
    LIMDD_TARG v1 = limdd_makenode(1, limdd_zero_edge(),
                                      limdd_bundle(LIMDD_LIM_IDENTITY, q2_zero));
    LIMDD_TARG root_a = limdd_makenode(0, limdd_bundle(LIMDD_LIM_IDENTITY, v0),
                                          limdd_bundle(LIMDD_LIM_IDENTITY, v1));
    LIMDD bell_a = limdd_bundle(limdd_lim_make(limdd_pauli_identity(), inv_sqrt2), root_a);

    /* --- shared: v1 == X_1 * v0, so the high edge reuses v0 --- */
    LIMDD_TARG root_b = limdd_makenode(0, limdd_bundle(LIMDD_LIM_IDENTITY, v0),
                                          limdd_bundle(limdd_lim_make(pauli("IXI"), EVBDD_ONE), v0));
    LIMDD bell_b = limdd_bundle(limdd_lim_make(limdd_pauli_identity(), inv_sqrt2), root_b);

    for (unsigned v = 0; v < (1u << NQUBITS); v++) {
        set_bits(bits, v);
        EVBDD_WGT a = limdd_eval(bell_a, bits, NQUBITS);
        EVBDD_WGT b = limdd_eval(bell_b, bits, NQUBITS);
        const bool nonzero = (v == 0x0) || (v == 0x3);
        test_assert(amp_is(a, nonzero ? SQRT2_INV : 0.0, 0.0));
        test_assert(amp_is(b, nonzero ? SQRT2_INV : 0.0, 0.0));
    }

    /* The shared form really is smaller: it never needed v1. */
    test_assert(root_a != root_b);
    test_assert(limdd_target(limdd_node_high(root_b)) == v0);
    test_assert(limdd_target(limdd_node_high(root_a)) == v1);
    return 0;
}

/* --- concurrent node creation -------------------------------------------- */

#define NSHARED 512
static LIMDD_TARG seen[NSHARED];

/**
 * Every worker builds the same nodes. They must agree on the indices: node
 * lookup compares children by index, so two workers disagreeing would produce
 * diagrams that are structurally equal and never compare equal.
 */
VOID_TASK_2(build_shared, size_t, lo, size_t, hi)
{
    if (hi - lo > 32) {
        size_t mid = lo + (hi - lo) / 2;
        SPAWN(build_shared, lo, mid);
        CALL(build_shared, mid, hi);
        SYNC(build_shared);
        return;
    }
    for (size_t i = lo; i < hi; i++) {
        limdd_pauli_t p = { i & 7, (i >> 3) & 7 };
        LIMDD high = limdd_bundle(limdd_lim_make(p, EVBDD_ONE), LIMDD_TERMINAL);
        LIMDD_TARG n = limdd_makenode(NQUBITS - 1, limdd_one_edge(), high);
        if (seen[i] != 0 && seen[i] != n) {
            fprintf(stderr, "node %zu built as %llu and also %llu\n",
                    i, (unsigned long long)seen[i], (unsigned long long)n);
            exit(1);
        }
        seen[i] = n;
    }
}

int
test_concurrent_makenode(void)
{
    memset(seen, 0, sizeof(seen));
    for (int round = 0; round < 8; round++) RUN(build_shared, 0, NSHARED);
    for (size_t i = 0; i < NSHARED; i++) test_assert(seen[i] != 0);
    return 0;
}

TASK_0(int, runtests)
{
    limdd_nodes_init(NQUBITS, 1LL << 16, 1LL << 16, 1LL << 16);

    if (test_edge_encoding()) return 1;
    printf("limdd edge fields round-trip:            ok\n");
    if (test_node_interning()) return 1;
    printf("limdd nodes intern, quasi-reduced:       ok\n");
    if (test_eval_basis_states()) return 1;
    printf("limdd eval on basis states:              ok\n");
    if (test_eval_lim_acts()) return 1;
    printf("limdd eval honours edge labels:          ok\n");
    if (test_eval_superposition()) return 1;
    printf("limdd eval on |+++>:                     ok\n");
    if (test_bell_state_two_ways()) return 1;
    printf("limdd shares a Bell state through a LIM: ok\n");
    if (test_concurrent_makenode()) return 1;
    printf("limdd makenode agrees across workers:    ok\n");

    printf("(%zu nodes, %zu LIMs)\n", limdd_node_table_count(), limdd_lim_table_count());

    limdd_nodes_quit();
    return 0;
}

int
main(void)
{
    lace_start(8, 0);

    /* Weight arithmetic goes through the operation cache, so the package has
     * to be up before any wgt_mul happens. */
    sylvan_set_sizes(1LL << 16, 1LL << 16, 1LL << 16, 1LL << 16);
    sylvan_init_package();
    sylvan_init_edge_weights(1LL << 16, 1LL << 16, 1e-14, WGT_COMPLEX_128, COMP_HASHMAP);

    int res = RUN(runtests);

    sylvan_edge_weights_free();
    sylvan_quit();
    lace_stop();
    return res;
}
