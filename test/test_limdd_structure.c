/*
 * Structural tests: the shape a LIMDD is supposed to have, not just the state
 * it denotes.
 *
 * The rest of the suite checks that a diagram means the right vector. That
 * cannot fail when the canonical form quietly stops reducing -- a correct but
 * exponentially large diagram denotes the same state as a small one. These
 * tests assert the size claims instead, which is where the data structure
 * earns its place:
 *
 *   Vinkhuijzen et al., LIMDD: a stabiliser state is a Tower-LIMDD, so at most
 *   one node per level and hence width 1.
 *
 *   Quist, Coopmans and Laarman, Thm (m independent Pauli stabilisers): an
 *   n-qubit state with m <= n linearly independent Pauli-string stabilisers
 *   has a LIMDD of width at most 2^(n-m).  Corollary: a state reached from
 *   |0..0> by Cliffords and t T gates has width at most 2^t.
 *
 * Note on node COUNTS: "one node per level" is an upper bound, not an
 * equality. This implementation is fully reduced, so a level whose two
 * children are Pauli-equivalent carries no node at all -- |+>^n is 0 nodes,
 * not n. Width is the invariant that survives that, which is why it is what
 * is asserted here.
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

#define MAXQ 24

static uint64_t rng_state = UINT64_C(0x243F6A8885A308D3);

static uint64_t
rnd(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

static uint32_t rnd_below(uint32_t n) { return (uint32_t)(rnd() % n); }

/* --- state constructions ------------------------------------------------- */

/**
 * A random stabiliser state: |0..0> under a random Clifford circuit.
 *
 * H, S and CX generate the Clifford group, so a long enough random word in
 * them reaches a state that is stabiliser but otherwise unconstrained -- which
 * is what makes this a stronger test than any fixed family.
 */
static LIMDD
random_stabiliser_state(uint32_t n, uint32_t depth)
{
    LIMDD e = limdd_all_zero_state(MAXQ);
    for (uint32_t i = 0; i < depth; i++) {
        switch (rnd_below(3)) {
        case 0: e = limdd_gate(e, GATEID_H, rnd_below(n), MAXQ); break;
        case 1: e = limdd_gate(e, GATEID_S, rnd_below(n), MAXQ); break;
        default: {
            if (n < 2) break;
            uint32_t c = rnd_below(n), t = rnd_below(n);
            while (t == c) t = rnd_below(n);
            bool ok = true; e = limdd_cgate_either(e, GATEID_X, c, t, MAXQ, &ok);
            break;
        }
        }
    }
    return e;
}

/**
 * The graph state of a random graph: H on every qubit, then CZ on each edge.
 *
 * Graph states are the family the LIMDD literature leans on, and they are the
 * ones a QMDD represents worst, so they are the sharpest tower test available.
 */
static LIMDD
graph_state(uint32_t n, unsigned edge_pct)
{
    LIMDD e = limdd_all_zero_state(MAXQ);
    for (uint32_t q = 0; q < n; q++) e = limdd_gate(e, GATEID_H, q, MAXQ);
    for (uint32_t a = 0; a < n; a++)
        for (uint32_t b = a + 1; b < n; b++)
            if (rnd() % 100 < edge_pct)
                { bool ok = true; e = limdd_cgate_either(e, GATEID_Z, a, b, MAXQ, &ok); }
    return e;
}

/**
 * A coset state: the uniform superposition over an affine subspace of F_2^n.
 *
 * Built as the image of |0..0> under H on a set of `k` free qubits, CNOTs that
 * make the remaining qubits linear functions of those, and X gates for the
 * offset. The result is sum over x in (C + v) of |x>, with C a k-dimensional
 * subspace -- a stabiliser state, so a tower, and a family whose LIMDD is
 * genuinely non-trivial because the linear dependencies differ per instance.
 */
static LIMDD
coset_state(uint32_t n, uint32_t k)
{
    LIMDD e = limdd_all_zero_state(MAXQ);
    if (k > n) k = n;
    for (uint32_t q = 0; q < k; q++) e = limdd_gate(e, GATEID_H, q, MAXQ);
    /* every dependent qubit is an XOR of a random non-empty set of free ones */
    for (uint32_t q = k; q < n; q++) {
        bool any = false;
        for (uint32_t f = 0; f < k; f++) {
            if (rnd() % 2) { bool ok = true; e = limdd_cgate_either(e, GATEID_X, f, q, MAXQ, &ok); any = true; }
        }
        if (!any && k > 0) { bool ok = true; e = limdd_cgate_either(e, GATEID_X, 0, q, MAXQ, &ok); }
    }
    for (uint32_t q = 0; q < n; q++)          /* the offset v */
        if (rnd() % 2) e = limdd_gate(e, GATEID_X, q, MAXQ);
    return e;
}

/* --- checks -------------------------------------------------------------- */

static int
check_tower(LIMDD e, uint32_t n, const char *what)
{
    size_t counts[MAXQ];
    const size_t total = limdd_level_counts(e, counts, MAXQ);
    for (uint32_t q = 0; q < MAXQ; q++) {
        if (counts[q] > 1) {
            fprintf(stderr,
                    "%s on %u qubits: level %u has %zu nodes, a tower has at "
                    "most 1 (total %zu)\n", what, n, q, counts[q], total);
            for (uint32_t i = 0; i < MAXQ; i++)
                fprintf(stderr, "    level %2u: %zu\n", i, counts[i]);
            return 1;
        }
    }
    /* The state is built on MAXQ qubits with only the first n touched, so the
     * |0> tail contributes one node per level as well: a tower over MAXQ, not
     * over n. */
    if (total > MAXQ) {
        fprintf(stderr, "%s (%u active qubits): %zu nodes over %u levels, "
                "more than one per level\n", what, n, total, (unsigned)MAXQ);
        return 1;
    }
    return 0;
}

static int
check_width_at_most(LIMDD e, uint32_t n, size_t bound, const char *what)
{
    const size_t w = limdd_width(e, MAXQ);
    if (w > bound) {
        fprintf(stderr, "%s on %u qubits: width %zu exceeds the bound %zu\n",
                what, n, w, bound);
        return 1;
    }
    return 0;
}

/* --- tests --------------------------------------------------------------- */

static int
test_stabiliser_towers(void)
{
    for (uint32_t n = 2; n <= MAXQ; n += 2) {
        for (int rep = 0; rep < 8; rep++) {
            LIMDD e = random_stabiliser_state(n, 6 * n);
            if (limdd_edge_is_zero(e)) continue;
            if (check_tower(e, n, "random stabiliser state")) return 1;
        }
    }
    return 0;
}

static int
test_graph_states(void)
{
    for (uint32_t n = 2; n <= MAXQ; n += 2) {
        for (unsigned pct = 20; pct <= 80; pct += 30) {
            LIMDD e = graph_state(n, pct);
            if (check_tower(e, n, "graph state")) return 1;
        }
    }
    return 0;
}

static int
test_coset_states(void)
{
    for (uint32_t n = 2; n <= MAXQ; n += 2) {
        for (uint32_t k = 1; k <= n; k += (n > 8 ? n / 4 : 1)) {
            LIMDD e = coset_state(n, k);
            if (check_tower(e, n, "coset state")) return 1;
        }
    }
    return 0;
}

/**
 * Cliffords plus t T gates: width at most 2^t.
 *
 * The T gates go on qubits chosen at random and are separated by Clifford
 * layers, so the t of them are not trivially cancelling or acting on the same
 * qubit in sequence, which would leave the nullity below t and make the bound
 * pass for the wrong reason.
 */
static int
test_nullity_width_bound(void)
{
    for (uint32_t t = 0; t <= 4; t++) {
        const size_t bound = (size_t)1 << t;
        for (uint32_t n = 4; n <= 12; n += 4) {
            for (int rep = 0; rep < 6; rep++) {
                LIMDD e = limdd_all_zero_state(MAXQ);
                for (uint32_t i = 0; i <= t; i++) {
                    for (uint32_t g = 0; g < 3 * n; g++) {
                        switch (rnd_below(3)) {
                        case 0: e = limdd_gate(e, GATEID_H, rnd_below(n), MAXQ); break;
                        case 1: e = limdd_gate(e, GATEID_S, rnd_below(n), MAXQ); break;
                        default: {
                            uint32_t c = rnd_below(n), tg = rnd_below(n);
                            while (tg == c) tg = rnd_below(n);
                            bool ok = true; e = limdd_cgate_either(e, GATEID_X, c, tg, MAXQ, &ok);
                            break;
                        }
                        }
                    }
                    if (i < t) e = limdd_gate(e, GATEID_T, rnd_below(n), MAXQ);
                }
                if (limdd_edge_is_zero(e)) continue;
                char what[64];
                snprintf(what, sizeof(what), "Clifford + %u T", t);
                if (check_width_at_most(e, n, bound, what)) return 1;
            }
        }
    }
    return 0;
}

/**
 * A tower is not an accident of the node count: the state must also still be
 * right. A diagram that reduced too far would pass every width check above.
 */
static int
test_towers_still_denote_the_state(void)
{
    const uint32_t n = 6;
    for (int rep = 0; rep < 20; rep++) {
        LIMDD e = random_stabiliser_state(n, 6 * n);
        if (limdd_edge_is_zero(e)) continue;
        /* a stabiliser state is normalised: the squared amplitudes sum to 1 */
        double sum = 0.0;
        for (uint64_t b = 0; b < (1ULL << n); b++) {
            bool bits[MAXQ];
            memset(bits, 0, sizeof(bits));
            for (uint32_t q = 0; q < n; q++) bits[q] = (b >> (n - 1 - q)) & 1;
            const complex_t c = weight_as_complex(limdd_eval(e, bits, MAXQ));
            sum += c.r * c.r + c.i * c.i;
        }
        if (sum < 1.0 - 1e-9 || sum > 1.0 + 1e-9) {
            fprintf(stderr, "tower denotes an unnormalised state: |psi|^2 = %.12f\n", sum);
            return 1;
        }
    }
    return 0;
}

/** The dot writer has to produce something, and not crash on a tower. */
static int
test_dot_output(void)
{
    LIMDD e = graph_state(5, 50);
    char *buf = NULL;
    size_t len = 0;
    FILE *f = open_memstream(&buf, &len);
    if (f == NULL) { fprintf(stderr, "open_memstream failed\n"); return 1; }
    limdd_fprintdot(f, e, MAXQ);
    fclose(f);
    if (len == 0 || strstr(buf, "digraph limdd") == NULL ||
        strstr(buf, "term") == NULL) {
        fprintf(stderr, "dot output does not look like a graph (%zu bytes)\n", len);
        free(buf);
        return 1;
    }
    free(buf);
    return 0;
}

TASK_0(int, runtests)
{
    limdd_nodes_init(MAXQ, 1LL << 21, 1LL << 21, 1LL << 23, 1LL << 23);

    if (test_stabiliser_towers()) return 1;
    printf("random stabiliser states are towers (2..24 qubits):  ok\n");
    if (test_graph_states()) return 1;
    printf("graph states are towers (2..24 qubits, 3 densities): ok\n");
    if (test_coset_states()) return 1;
    printf("coset states are towers (affine subspaces):          ok\n");
    if (test_towers_still_denote_the_state()) return 1;
    printf("towers still denote a normalised state:              ok\n");
    if (test_nullity_width_bound()) return 1;
    printf("Clifford + t T gates: width <= 2^t (t = 0..4):       ok\n");
    if (test_dot_output()) return 1;
    printf("dot output is a graph:                               ok\n");

    printf("(%zu nodes, %zu LIMs)\n",
           limdd_node_table_count(), limdd_lim_table_count());
    limdd_nodes_quit();
    return 0;
}

static int
run_with(int backend, const char *name)
{
    printf("== LIMDD structure with %s edge weights ==\n", name);
    /*
     * One worker, deliberately.
     *
     * With four, this test is not reproducible: the same fixed seed gave
     * widths of 91, 19, 5, 62 and 23 on five consecutive runs against a bound
     * of 16. The diagram is not racy -- makeedge agreeing across workers is
     * checked in test_limdd_canon -- it is the FLOAT weight table. Which
     * weights merge under an absolute tolerance depends on the order they are
     * interned, that order depends on which worker gets there first, and a
     * lost merge is a node that fails to be shared. The size claims are then
     * schedule-dependent, which is a precision result rather than a structural
     * one and does not belong in a pass/fail test.
     *
     * The exact backend has no such dependence: equality is equality however
     * the table fills.
     */
    lace_start(1, 0);
    sylvan_set_sizes(1LL << 18, 1LL << 18, 1LL << 18, 1LL << 18);
    sylvan_init_package();
    /*
     * Through the simulator's initialiser, not sylvan_init_edge_weights: it is
     * what populates the gate table AND sets the normalisation strategy. Going
     * around it leaves the strategy at its default, and a LIMDD built under the
     * wrong one is not the diagram the theory is about -- which is exactly how
     * an earlier draft of this file produced a width of 19 where the bound is
     * 16, and blamed the implementation for it.
     */
    qsylvan_init_simulator(1LL << 22, 1LL << 22, -1, backend, NORM_LOW);

    const int res = RUN(runtests);

    sylvan_quit();
    lace_stop();
    return res;
}

int
main(void)
{
    if (run_with(COMP_HASHMAP, "complex")) return 1;
    /* And with exact coefficients, where no merge can be lost to rounding:
     * the size claims are statements about the state, so they must hold in
     * both, and a difference between the two is a precision failure. */
    if (run_with(QISQ2_MAP, "exact (Q[i,sqrt2])")) return 1;
    return 0;
}
