/*
 * Corner cases, checked by shape rather than by state.
 *
 * test_limdd_structure asserts the size claims on random families. This file
 * does the opposite: a handful of states whose diagram can be worked out on
 * paper, including the degenerate ones that the general families never
 * produce -- an empty diagram, a diagram that is only skips, a node whose low
 * edge is zero, a state that is not a stabiliser state at all.
 *
 * Everything here runs on the EXACT backend. The shapes below are statements
 * about the state, and under float a lost merge turns any of them into a
 * different but still correct diagram, which would make the assertions
 * measure the rounding rather than the reduction.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sylvan.h>
#include <sylvan_int.h>

#include "qsylvan_limdd_node.h"
#include "qsylvan_limdd_lim.h"
#include "qsylvan_limdd_ops.h"
#include "qsylvan_limdd_canon.h"
#include "qsylvan_limdd_inspect.h"
#include "qsylvan_gates.h"
#include "qsylvan_simulator.h"

#include "test_assert.h"

#define NQ 8

static int failures = 0;

static void
expect(bool ok, const char *what, const char *detail)
{
    if (!ok) { fprintf(stderr, "FAIL %s: %s\n", what, detail); failures++; }
}

static void
expect_shape(LIMDD e, size_t want_nodes, size_t want_width, const char *what)
{
    const size_t nodes = limdd_nodecount(e, NQ);
    const size_t width = limdd_width(e, NQ);
    if (nodes != want_nodes || width != want_width) {
        fprintf(stderr, "FAIL %s: %zu nodes / width %zu, expected %zu / %zu\n",
                what, nodes, width, want_nodes, want_width);
        size_t counts[NQ];
        limdd_level_counts(e, counts, NQ);
        for (unsigned i = 0; i < NQ; i++)
            fprintf(stderr, "    level %u: %zu\n", i, counts[i]);
        failures++;
    }
}

/* --- builders ------------------------------------------------------------ */

static LIMDD
all_zero(void) { return limdd_all_zero_state(NQ); }

/** H on every qubit: |+>^n. */
static LIMDD
plus_all(void)
{
    LIMDD e = all_zero();
    for (uint32_t q = 0; q < NQ; q++) e = limdd_gate(e, GATEID_H, q, NQ);
    return e;
}

/** |->^n, which is |+>^n with a Z on every qubit. */
static LIMDD
minus_all(void)
{
    LIMDD e = plus_all();
    for (uint32_t q = 0; q < NQ; q++) e = limdd_gate(e, GATEID_Z, q, NQ);
    return e;
}

/** GHZ: H on qubit 0, then a CNOT ladder. */
static LIMDD
ghz(void)
{
    LIMDD e = all_zero();
    e = limdd_gate(e, GATEID_H, 0, NQ);
    for (uint32_t q = 1; q < NQ; q++) {
        bool ok = true;
        e = limdd_cgate_either(e, GATEID_X, 0, q, NQ, &ok);
    }
    return e;
}

/** |1>^n: X on every qubit. Every node then has a ZERO low edge. */
static LIMDD
ones(void)
{
    LIMDD e = all_zero();
    for (uint32_t q = 0; q < NQ; q++) e = limdd_gate(e, GATEID_X, q, NQ);
    return e;
}

/** Alternating |0>|+>|0>|+>..., so skipped and unskipped levels interleave. */
static LIMDD
alternating(void)
{
    LIMDD e = all_zero();
    for (uint32_t q = 1; q < NQ; q += 2) e = limdd_gate(e, GATEID_H, q, NQ);
    return e;
}

/**
 * The W state on NQ qubits, which is NOT a stabiliser state.
 *
 * Built by hand from amplitudes rather than by a circuit: the rotations a W
 * state needs are outside Clifford+T, and the point here is the shape, not
 * how it was reached.
 */
static LIMDD
w_state(void)
{
    /* |W> = sum_k |0..1..0>, unnormalised: a sum of NQ basis states, each
     * reached from |0..0> by a single X. */
    LIMDD acc = limdd_zero_edge();
    for (uint32_t k = 0; k < NQ; k++) {
        const LIMDD term = limdd_gate(all_zero(), GATEID_X, k, NQ);
        acc = limdd_plus(acc, term, 0);
    }
    return acc;
}

/* --- tests --------------------------------------------------------------- */

static void
test_degenerate(void)
{
    /* The zero vector: no nodes at all, and the edge itself is zero. */
    const LIMDD z = limdd_zero_edge();
    expect(limdd_edge_is_zero(z), "zero edge", "limdd_edge_is_zero was false");
    expect_shape(z, 0, 0, "the zero vector");

    /*
     * |+>^n is every level skipped, so the diagram is the terminal and
     * nothing else. This is the case a node COUNT cannot describe -- an
     * n-qubit state in zero nodes -- and the reason width is the invariant
     * the other tests assert.
     */
    expect_shape(plus_all(), 0, 0, "|+>^n");

    /* |->^n is the same diagram; the Z's ride on the root label. */
    const LIMDD m = minus_all();
    expect_shape(m, 0, 0, "|->^n");
    expect(limdd_target(m) == limdd_target(plus_all()),
           "|->^n and |+>^n", "different targets; they differ only by a LIM");
    expect(limdd_label(m) != limdd_label(plus_all()),
           "|->^n and |+>^n", "identical labels; |-> is not |+>");
}

static void
test_basis_states(void)
{
    /* |0>^n: no level can skip, because one child is zero and the other is
     * not, so there is a node per level and each is used once. */
    expect_shape(all_zero(), NQ, 1, "|0>^n");

    /* |1>^n: the mirror image, with a zero LOW edge at every level. */
    expect_shape(ones(), NQ, 1, "|1>^n");

    /*
     * And they are the SAME node, reached by a different label: |1..1> is
     * X^(x)n |0..0>, and identifying states that differ by a local Pauli is
     * precisely what a LIMDD does. An EVDD needs two diagrams here.
     */
    expect(limdd_target(all_zero()) == limdd_target(ones()),
           "|0>^n vs |1>^n", "different targets, but they differ by X^(x)n");
    expect(limdd_label(all_zero()) != limdd_label(ones()),
           "|0>^n vs |1>^n", "identical labels; |0..0> is not |1..1>");
}

static void
test_skips_interleave(void)
{
    /*
     * |0>|+>|0>|+>... Only the |0> levels survive; each |+> level has equal
     * children and is skipped. With NQ even and qubits 1,3,5,... put into |+>,
     * that leaves NQ/2 nodes.
     */
    expect_shape(alternating(), NQ / 2, 1, "|0>|+>|0>|+>...");
}

static void
test_ghz_is_a_tower(void)
{
    /* GHZ is a stabiliser state: one node per level, all NQ of them, because
     * no level can skip -- the two children are Pauli-equivalent but the node
     * is still needed to hold the branch. */
    expect_shape(ghz(), NQ, 1, "GHZ");
}

static void
test_pauli_images_share_every_node(void)
{
    /*
     * The whole point of the data structure: a state and its image under a
     * local Pauli are the SAME node, reached by a different label. If this
     * ever stops holding, the diagram silently doubles and every state-based
     * test still passes.
     */
    const LIMDD g = ghz();
    for (uint32_t q = 0; q < NQ; q++) {
        for (unsigned gate = 0; gate < 3; gate++) {
            const uint32_t id = (gate == 0) ? GATEID_X
                              : (gate == 1) ? GATEID_Y : GATEID_Z;
            const LIMDD img = limdd_gate(g, id, q, NQ);
            if (limdd_edge_is_zero(img)) continue;
            char buf[96];
            snprintf(buf, sizeof(buf), "%c_%u on GHZ",
                     gate == 0 ? 'X' : gate == 1 ? 'Y' : 'Z', q);
            expect(limdd_target(img) == limdd_target(g), buf,
                   "a Pauli image landed on a different node");
            expect_shape(img, NQ, 1, buf);
        }
    }
}

static void
test_global_phase_shares_the_node(void)
{
    /* i|psi> is |psi> with a different root weight, never a different node. */
    const LIMDD g = ghz();
    LIMDD p = g;
    for (int i = 0; i < 4; i++) p = limdd_gate(p, GATEID_S, 0, NQ);  /* S^4 = I */
    expect(limdd_target(p) == limdd_target(g), "S^4 on GHZ",
           "a global phase changed the node");
}

static void
test_w_state_is_not_a_tower(void)
{
    /*
     * A negative control. The W state is not a stabiliser state, so it must
     * NOT come out as a tower -- if it did, the tower tests in
     * test_limdd_structure would be passing vacuously.
     */
    const LIMDD w = w_state();
    expect(!limdd_edge_is_zero(w), "W state", "came out as the zero vector");
    const size_t width = limdd_width(w, NQ);
    expect(width > 1, "W state",
           "width 1, i.e. a tower -- a non-stabiliser state cannot be one");
    if (width > 1) printf("  (W state width %zu, as expected > 1)\n", width);
}

TASK_0(int, runtests)
{
    limdd_nodes_init(NQ, 1LL << 18, 1LL << 18, 1LL << 20, 1LL << 20);

    test_degenerate();          printf("degenerate diagrams (zero, |+>^n, |->^n):  ok\n");
    test_basis_states();        printf("basis states and zero low/high edges:      ok\n");
    test_skips_interleave();    printf("skipped and unskipped levels interleaved:  ok\n");
    test_ghz_is_a_tower();      printf("GHZ is one node per level:                 ok\n");
    test_pauli_images_share_every_node();
                                printf("every Pauli image shares every node:       ok\n");
    test_global_phase_shares_the_node();
                                printf("a global phase does not make a new node:   ok\n");
    test_w_state_is_not_a_tower();
                                printf("the W state is NOT a tower:                ok\n");

    limdd_nodes_quit();
    return failures != 0;
}

int
main(void)
{
    lace_start(1, 0);
    sylvan_set_sizes(1LL << 18, 1LL << 18, 1LL << 18, 1LL << 18);
    sylvan_init_package();
    /* Exact weights: see the note at the top. */
    qsylvan_init_simulator(1LL << 20, 1LL << 20, -1, QISQ2_MAP, NORM_LOW);

    const int res = RUN(runtests);

    sylvan_quit();
    lace_stop();
    return res;
}
