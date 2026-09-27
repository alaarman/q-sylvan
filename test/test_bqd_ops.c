/*
 * The BQD operations, checked the way the paper checked them: node for node
 * against the canonical diagram built from the amplitude vector.
 *
 *   prop:prodscalar   the product of two full-support phase states equals the
 *                     canonical build of the pointwise product vector, as an
 *                     edge, so the recursion is canonical and not merely right
 *   prop:diag         a monomial gate is the same, and its operand has one node
 *                     per level
 *   thm:size          an IQP circuit -- Hadamards then Z, S, T, CZ, CS, CCZ --
 *                     keeps every intermediate state within
 *                     sum_{i<3} C(v, i) + 1 nodes at the level with v decided
 *                     variables, and the final amplitudes are the phase
 *                     polynomial
 *
 * Scalar family, exact weights in Q(w_8).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sylvan.h>
#include <sylvan_int.h>

#include "qsylvan_bqd.h"
#include "qsylvan_bqd_ops.h"
#include "qsylvan_limdd_inspect.h"
#include "qsylvan_simulator.h"
#include "sylvan_edge_weights_qisq2.h"

#define NQ 12
#define MAXV (1u << NQ)

static int failures = 0;
static EVBDD_WGT pw[8];

static uint64_t rng_state = UINT64_C(0xC0FFEE1234567);
static uint64_t rnd(void)
{
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17;
    return rng_state;
}
static uint64_t rnd_below(uint64_t n) { return rnd() % n; }

static void
expect(bool ok, const char *what, const char *detail)
{
    if (!ok) { fprintf(stderr, "FAIL %s: %s\n", what, detail); failures++; }
}

static uint64_t
binom(unsigned n, unsigned k)
{
    if (k > n) return 0;
    uint64_t r = 1;
    for (unsigned i = 1; i <= k; i++) r = r * (n - k + i) / i;
    return r;
}

typedef struct { uint64_t W; unsigned c; } term_t;

/** A level-k phase polynomial over Z_8, as in test_bqd. */
static unsigned
rand_level(unsigned n, unsigned k, term_t *terms)
{
    unsigned nt = 0;
    for (unsigned d = 1; d <= k; d++) {
        const unsigned step = 1u << (3 - k + d - 1);
        if (step > 8) continue;
        for (uint64_t W = 1; W < (UINT64_C(1) << n); W++) {
            if (__builtin_popcountll(W) != (int)d) continue;
            const unsigned c = (unsigned)rnd_below(8 / step) * step;
            if (c) { terms[nt].W = W; terms[nt].c = c; nt++; }
        }
    }
    return nt;
}

static unsigned
poly_val(const term_t *terms, unsigned nt, uint64_t x)
{
    unsigned v = 0;
    for (unsigned i = 0; i < nt; i++) if ((x & terms[i].W) == terms[i].W) v += terms[i].c;
    return v & 7;
}

/** A full-support phase state with a random global scalar. */
static void
phase_state(unsigned n, unsigned k, EVBDD_WGT *f)
{
    term_t *terms = malloc(4096 * sizeof(term_t));
    const unsigned nt = rand_level(n, k, terms);
    const EVBDD_WGT scale = pw[rnd_below(8)];
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
        f[x] = wgt_mul(scale, pw[poly_val(terms, nt, x)]);
    free(terms);
}

static void
check_product(void)
{
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *g = malloc(MAXV * sizeof(EVBDD_WGT)),
              *h = malloc(MAXV * sizeof(EVBDD_WGT));
    unsigned bad = 0, tot = 0;
    for (unsigned n = 1; n <= 8; n++) for (int rep = 0; rep < 30; rep++) {
        phase_state(n, 2 + (unsigned)rnd_below(2), f);
        phase_state(n, 2 + (unsigned)rnd_below(2), g);
        for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) h[x] = wgt_mul(f[x], g[x]);
        const BQD F = bqd_from_vector(f, n), G = bqd_from_vector(g, n);
        const BQD want = bqd_from_vector(h, n);
        const BQD got = bqd_product(F, G, 0);
        tot++;
        if (got != want) bad++;
    }
    char buf[128];
    snprintf(buf, sizeof(buf), "%u of %u products differ from the canonical build of f.g", bad, tot);
    expect(bad == 0, "prop:prodscalar", buf);
    free(f); free(g); free(h);
}

static void
check_diagonal(void)
{
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *h = malloc(MAXV * sizeof(EVBDD_WGT));
    size_t counts[NQ];
    unsigned bad = 0, badshape = 0, tot = 0;
    for (unsigned n = 1; n <= 8; n++) for (int rep = 0; rep < 30; rep++) {
        phase_state(n, 2 + (unsigned)rnd_below(2), f);
        const unsigned arity = 1 + (unsigned)rnd_below(n < 3 ? n : 3);
        uint64_t A = 0;
        while (__builtin_popcountll(A) < (int)arity) A |= UINT64_C(1) << rnd_below(n);
        const EVBDD_WGT phase = pw[1 + rnd_below(7)];

        /* the gate's own diagram is one chain of monomial nodes and the
         * constant-one chain they point at: at most two nodes per level */
        const BQD gate = bqd_monomial(A, phase, n);
        limdd_level_counts(gate, counts, n);
        for (unsigned v = 0; v < n; v++) if (counts[v] < 1 || counts[v] > 2) { badshape++; break; }

        for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
            h[x] = ((x & A) == A) ? wgt_mul(f[x], phase) : f[x];
        const BQD want = bqd_from_vector(h, n);
        const BQD got = bqd_apply_diagonal(bqd_from_vector(f, n), A, phase, n);
        tot++;
        if (got != want) bad++;
    }
    char buf[128];
    snprintf(buf, sizeof(buf), "%u of %u gate applications differ from the canonical build", bad, tot);
    expect(bad == 0, "prop:diag", buf);
    snprintf(buf, sizeof(buf), "%u monomial diagrams have more than two nodes at a level", badshape);
    expect(badshape == 0, "monomial shape", buf);
    free(f); free(h);
}

/** An IQP circuit: |+>^n, then random diagonal gates of level at most three. */
static void
check_iqp(void)
{
    size_t counts[NQ];
    unsigned states = 0, bound_bad = 0, amp_bad = 0;
    size_t max_width = 0;
    term_t *poly = malloc(4096 * sizeof(term_t));
    for (unsigned n = 2; n <= 10; n++) {
        BQD e = limdd_bundle(LIMDD_LIM_IDENTITY, bqd_ones(0, n));     /* |+>^n up to scale */
        unsigned npoly = 0;
        for (unsigned gate = 0; gate < 6 * n; gate++) {
            /* Z, S, T, CZ, CS, CCZ: a monomial of arity d with exponent e in Z_8 */
            static const unsigned arity[6] = { 1, 1, 1, 2, 2, 3 };
            static const unsigned expo[6]  = { 4, 2, 1, 4, 2, 4 };
            const unsigned gi = (unsigned)rnd_below(n >= 3 ? 6 : 5);
            uint64_t A = 0;
            while (__builtin_popcountll(A) < (int)arity[gi]) A |= UINT64_C(1) << rnd_below(n);
            e = bqd_apply_diagonal(e, A, pw[expo[gi]], n);
            poly[npoly].W = A; poly[npoly].c = expo[gi]; npoly++;

            states++;
            limdd_level_counts(e, counts, n);
            for (unsigned v = 0; v < n; v++) {
                if (counts[v] > max_width) max_width = counts[v];
                const uint64_t bound = binom(v, 0) + binom(v, 1) + binom(v, 2) + 1;
                if (counts[v] > bound) { bound_bad++; break; }
            }
        }
        /* the amplitudes are w_8^{P(x)}, since the start was all ones */
        for (int q = 0; q < 32; q++) {
            const uint64_t x = rnd_below(UINT64_C(1) << n);
            if (bqd_eval(e, n, x) != pw[poly_val(poly, npoly, x)]) { amp_bad++; break; }
        }
    }
    char buf[160];
    printf("  IQP: %u intermediate states on 2..10 qubits, max width %zu\n", states, max_width);
    snprintf(buf, sizeof(buf), "%u of %u intermediate states exceed sum_{i<3} C(v,i) + 1", bound_bad, states);
    expect(bound_bad == 0, "thm:size on an IQP circuit", buf);
    snprintf(buf, sizeof(buf), "%u circuits have a final amplitude that is not w_8^P(x)", amp_bad);
    expect(amp_bad == 0, "IQP amplitudes", buf);
    free(poly);
}

TASK_0(int, runtests)
{
    bqd_init(BQD_FAMILY_SCALAR, NQ, 1LL << 20, 1LL << 20, 1LL << 22, 1LL << 20);
    pw[0] = EVBDD_ONE;
    const EVBDD_WGT w = qisq2_lookup(0, 1, 1, 2, 0, 1, 1, 2);
    for (int e = 1; e < 8; e++) pw[e] = wgt_mul(pw[e - 1], w);

    check_product();          printf("product equals the canonical build (prop:prodscalar): %s\n", failures ? "FAILED" : "ok");
    const int f1 = failures;
    check_diagonal();         printf("diagonal gates equal the canonical build (prop:diag):  %s\n", failures > f1 ? "FAILED" : "ok");
    const int f2 = failures;
    check_iqp();              printf("IQP circuits stay within thm:size:                     %s\n", failures > f2 ? "FAILED" : "ok");

    bqd_quit();
    return failures != 0;
}

int
main(void)
{
    lace_start(4, 0);
    sylvan_set_sizes(1LL << 20, 1LL << 20, 1LL << 20, 1LL << 20);
    sylvan_init_package();
    qsylvan_init_simulator(1LL << 22, 1LL << 22, -1, QISQ2_MAP, NORM_LOW);
    const int res = RUN(runtests);
    sylvan_quit();
    lace_stop();
    return res;
}
