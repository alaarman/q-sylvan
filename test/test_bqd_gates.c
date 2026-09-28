/*
 * Gates on a BQD, and the pointwise operations they are made of, checked
 * against the vector: every result must be the canonical diagram of what a
 * dense simulation gives, as an edge, so a result that is right in value and
 * wrong in shape fails too. The dense simulation is the paper's own Hadamard
 * (tab:ops: rebuild the vector, O(n^2 2^n)), which is what the recursion has
 * to agree with.
 *
 *   basis states      bqd_basis_state(x) is the diagram of e_x
 *   compose           a diagram is the composition of its two cofactors, and
 *                     (f_0, 0) and (0, f_1) compose to what they denote
 *   pointwise         f.g and f + g on random vectors with zeros, with
 *                     algebraic values, and with cancellation to zero
 *   gates             every 2x2 gate in the table on every qubit, including
 *                     the projectors; H H = I
 *   controlled        controls above the target for any gate, anywhere for
 *                     a phase; cgate_either and swap
 *   circuits          random Clifford+T circuits from |0...0>, Hadamards
 *                     anywhere, the diagram compared after every gate
 *
 * Scalar family, exact weights only. On floats these operations are not
 * reliable, and no tolerance makes them so; this is measured, not assumed.
 * A cancellation that should give 0 leaves a residue, and under the copy
 * clause a residue is not a small error: it turns a copy point into a ratio
 * point, so the value stored there is read as a ratio and multiplied by the
 * residue, and an amplitude of size 1 comes out 0. At the default tolerance,
 * 1e-14, every run of these checks failed, in a different place each time,
 * since the schedule decides which values merge. At 1e-10 or 1e-8 the
 * residues snap to 0 and single gates pass, but one circuit of 67 gates on 7
 * qubits then went wrong at its last gate by 0.044, which it did not at
 * 1e-14: a tolerance coarse enough to catch the residues merges values that
 * differ. The same circuit is exact, and canonical, on the algebraic backend.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sylvan.h>
#include <sylvan_int.h>

#include "qsylvan_bqd.h"
#include "qsylvan_bqd_ops.h"
#include "qsylvan_bqd_gates.h"
#include "qsylvan_limdd_inspect.h"
#include "qsylvan_simulator.h"
#include "sylvan_edge_weights_qisq2.h"

#define NQ 10
#define MAXV (1u << NQ)

static int failures = 0;
static EVBDD_WGT pw[8];

static uint64_t rng_state = UINT64_C(0xB0D1E5EED5);
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

static void
report(const char *what, unsigned bad, unsigned tot)
{
    char buf[128];
    snprintf(buf, sizeof(buf), "%u of %u differ", bad, tot);
    expect(bad == 0, what, buf);
    printf("  %-44s %s (%u)\n", what, bad ? "FAILED" : "ok", tot);
}

/* --- vectors --------------------------------------------------------------- */

/** Roots of unity and zeros, or, with `alg`, small algebraic values too. */
static void
rand_vector(unsigned n, double pzero, bool alg, EVBDD_WGT *f)
{
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) {
        if ((double)rnd_below(1000) / 1000.0 < pzero) { f[x] = EVBDD_ZERO; continue; }
        if (alg && rnd_below(2)) {
            long a, b, c, d;
            do {
                a = (long)rnd_below(5) - 2; b = (long)rnd_below(5) - 2;
                c = (long)rnd_below(5) - 2; d = (long)rnd_below(5) - 2;
            } while (a == 0 && b == 0 && c == 0 && d == 0);
            f[x] = qisq2_lookup(a, 1, b, 1, c, 1, d, 1);
        } else {
            f[x] = pw[rnd_below(8)];
        }
    }
}

/** U on qubit q of a vector (qubit q at bit n-1-q), where the qubits of cmask are 1. */
static void
dense_gate(EVBDD_WGT *v, unsigned n, uint32_t gid, uint64_t cmask, unsigned q)
{
    const uint64_t qb = UINT64_C(1) << (n - 1 - q);
    uint64_t cm = 0;
    for (unsigned c = 0; c < n; c++) if ((cmask >> c) & 1) cm |= UINT64_C(1) << (n - 1 - c);
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) {
        if ((x & qb) || (x & cm) != cm) continue;
        const EVBDD_WGT a = v[x], b = v[x | qb];
        v[x]      = wgt_add(wgt_mul(gates[gid][0], a), wgt_mul(gates[gid][1], b));
        v[x | qb] = wgt_add(wgt_mul(gates[gid][2], a), wgt_mul(gates[gid][3], b));
    }
}

static void
dense_swap(EVBDD_WGT *v, unsigned n, unsigned a, unsigned b)
{
    const uint64_t ab = UINT64_C(1) << (n - 1 - a), bb = UINT64_C(1) << (n - 1 - b);
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) {
        if ((x & ab) && !(x & bb)) {
            const uint64_t y = (x & ~ab) | bb;
            const EVBDD_WGT t = v[x]; v[x] = v[y]; v[y] = t;
        }
    }
}

/** Whether e is the canonical diagram of v: one edge comparison. */
static bool
denotes(BQD e, const EVBDD_WGT *v, unsigned n)
{
    return e == bqd_from_vector(v, n);
}

/* --- checks ---------------------------------------------------------------- */

static void
check_basis_and_compose(void)
{
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *g = malloc(MAXV * sizeof(EVBDD_WGT));
    unsigned bad_basis = 0, tb = 0, bad_comp = 0, tc = 0;
    for (unsigned n = 1; n <= 8; n++) {
        for (int rep = 0; rep < 8; rep++) {
            const uint64_t x = rnd_below(UINT64_C(1) << n);
            memset(f, 0, sizeof(EVBDD_WGT) << n);
            for (uint64_t y = 0; y < (UINT64_C(1) << n); y++) f[y] = (y == x) ? EVBDD_ONE : EVBDD_ZERO;
            tb++;
            if (bqd_basis_state(x, n) != bqd_from_vector(f, n)) bad_basis++;
        }
        const double dens[3] = { 0.0, 0.3, 0.6 };
        for (unsigned di = 0; di < 3; di++) for (int rep = 0; rep < 12; rep++) {
            rand_vector(n, dens[di], true, f);
            const uint64_t h = UINT64_C(1) << (n - 1);
            const BQD e = bqd_from_vector(f, n);
            tc++;
            if (limdd_edge_is_zero(e)) continue;
            const BQD c0 = bqd_cofactor(e, 0, 0), c1 = bqd_cofactor(e, 0, 1);
            bool ok = bqd_compose(0, c0, c1) == e;
            /* (f_0, 0) and (0, f_1) */
            memcpy(g, f, sizeof(EVBDD_WGT) << n);
            for (uint64_t y = h; y < 2 * h; y++) g[y] = EVBDD_ZERO;
            ok = ok && bqd_compose(0, c0, limdd_zero_edge()) == bqd_from_vector(g, n);
            memcpy(g, f, sizeof(EVBDD_WGT) << n);
            for (uint64_t y = 0; y < h; y++) g[y] = EVBDD_ZERO;
            ok = ok && bqd_compose(0, limdd_zero_edge(), c1) == bqd_from_vector(g, n);
            if (!ok) bad_comp++;
        }
    }
    report("basis states", bad_basis, tb);
    report("compose of the two cofactors", bad_comp, tc);
    free(f); free(g);
}

static void
check_pointwise(void)
{
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *g = malloc(MAXV * sizeof(EVBDD_WGT)),
              *h = malloc(MAXV * sizeof(EVBDD_WGT));
    unsigned bad_mul = 0, bad_add = 0, tot = 0, cancelled = 0;
    for (unsigned n = 1; n <= 7; n++) for (int rep = 0; rep < 30; rep++) {
        const double pz = (double)(rep % 3) * 0.3;
        rand_vector(n, pz, rep & 1, f);
        rand_vector(n, pz, rep & 1, g);
        /* cancellation: g = -f on a random set, so f + g has zeros f and g lack */
        for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
            if (rnd_below(3) == 0) { g[x] = wgt_neg(f[x]); if (f[x] != EVBDD_ZERO) cancelled++; }
        const BQD F = bqd_from_vector(f, n), G = bqd_from_vector(g, n);
        tot++;
        for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) h[x] = wgt_mul(f[x], g[x]);
        if (!denotes(bqd_multiply(F, G), h, n)) bad_mul++;
        for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) h[x] = wgt_add(f[x], g[x]);
        if (!denotes(bqd_add(F, G), h, n)) bad_add++;
    }
    printf("  %u vector pairs, %u points where f + g cancels to zero\n", tot, cancelled);
    report("f . g, any support", bad_mul, tot);
    report("f + g, with cancellation", bad_add, tot);
    free(f); free(g); free(h);
}

static const uint32_t all_gates[] = {
    GATEID_I, GATEID_X, GATEID_Y, GATEID_Z, GATEID_H, GATEID_S, GATEID_Sdag,
    GATEID_T, GATEID_Tdag, GATEID_sqrtX, GATEID_sqrtXdag, GATEID_sqrtY,
    GATEID_sqrtYdag, GATEID_proj0, GATEID_proj1,
};
#define N_ALL_GATES (sizeof(all_gates) / sizeof(all_gates[0]))

static void
check_gates(void)
{
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *v = malloc(MAXV * sizeof(EVBDD_WGT));
    unsigned bad = 0, tot = 0, bad_hh = 0, thh = 0;
    for (unsigned n = 1; n <= 7; n++) for (int rep = 0; rep < 6; rep++) {
        rand_vector(n, (double)(rep % 3) * 0.3, rep & 1, f);
        const BQD F = bqd_from_vector(f, n);
        for (unsigned gi = 0; gi < N_ALL_GATES; gi++) for (unsigned q = 0; q < n; q++) {
            memcpy(v, f, sizeof(EVBDD_WGT) << n);
            dense_gate(v, n, all_gates[gi], 0, q);
            tot++;
            if (!denotes(bqd_gate(F, all_gates[gi], q, n), v, n)) {
                if (bad < 3) fprintf(stderr, "  gate %u on qubit %u of %u\n", all_gates[gi], q, n);
                bad++;
            }
        }
        for (unsigned q = 0; q < n; q++) {
            thh++;
            const BQD HH = bqd_gate(bqd_gate(F, GATEID_H, q, n), GATEID_H, q, n);
            if (HH != F) bad_hh++;
        }
    }
    report("every table gate on every qubit", bad, tot);
    report("H H = I", bad_hh, thh);
    free(f); free(v);
}

static void
check_controlled(void)
{
    static const uint32_t cg[] = { GATEID_X, GATEID_Y, GATEID_Z, GATEID_H, GATEID_S,
                                   GATEID_T, GATEID_sqrtX, GATEID_Tdag };
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *v = malloc(MAXV * sizeof(EVBDD_WGT));
    unsigned bad = 0, tot = 0, bad_phase = 0, tp = 0, bad_either = 0, te = 0;
    for (unsigned n = 2; n <= 7; n++) for (int rep = 0; rep < 40; rep++) {
        rand_vector(n, (double)(rep % 3) * 0.3, rep & 1, f);
        const BQD F = bqd_from_vector(f, n);
        const uint32_t gid = cg[rnd_below(sizeof(cg) / sizeof(cg[0]))];

        /* controls above the target, any gate */
        const unsigned q = 1 + (unsigned)rnd_below(n - 1);
        uint64_t cm = 0;
        while (cm == 0) cm = rnd_below(UINT64_C(1) << q);
        memcpy(v, f, sizeof(EVBDD_WGT) << n);
        dense_gate(v, n, gid, cm, q);
        tot++;
        if (!denotes(bqd_cgate(F, gid, cm, q, n), v, n)) bad++;

        /* a phase with its controls anywhere */
        static const uint32_t ph[] = { GATEID_Z, GATEID_S, GATEID_T, GATEID_Sdag };
        const uint32_t pg = ph[rnd_below(4)];
        const unsigned t = (unsigned)rnd_below(n);
        uint64_t pm = 0;
        while (pm == 0) pm = rnd_below(UINT64_C(1) << n) & ~(UINT64_C(1) << t);
        memcpy(v, f, sizeof(EVBDD_WGT) << n);
        dense_gate(v, n, pg, pm, t);
        tp++;
        if (!denotes(bqd_cgate(F, pg, pm, t, n), v, n)) bad_phase++;

        /* one control on either side, and a swap */
        const unsigned a = (unsigned)rnd_below(n);
        unsigned b = (unsigned)rnd_below(n - 1); if (b >= a) b++;
        const uint32_t eg = rnd_below(2) ? GATEID_X : GATEID_Z;
        memcpy(v, f, sizeof(EVBDD_WGT) << n);
        dense_gate(v, n, eg, UINT64_C(1) << a, b);
        bool ok;
        te++;
        if (!denotes(bqd_cgate_either(F, eg, a, b, n, &ok), v, n) || !ok) bad_either++;
        memcpy(v, f, sizeof(EVBDD_WGT) << n);
        dense_swap(v, n, a, b);
        te++;
        if (!denotes(bqd_swap(F, a, b, n), v, n)) bad_either++;
    }
    report("controlled, controls above the target", bad, tot);
    report("controlled phase, controls anywhere", bad_phase, tp);
    report("cgate_either and swap", bad_either, te);
    free(f); free(v);
}

/**
 * Random Clifford+T circuits from |0...0>, the gate set of the qasm runner,
 * Hadamards anywhere. After every gate the diagram must be the vector's.
 */
static void
check_circuits(unsigned max_n, unsigned per_n)
{
    EVBDD_WGT *v = malloc(MAXV * sizeof(EVBDD_WGT));
    unsigned bad = 0, states = 0, not_full = 0;
    size_t max_nodes = 0;
    for (unsigned n = 1; n <= max_n; n++) for (unsigned c = 0; c < per_n; c++) {
        BQD e = bqd_basis_state(0, n);
        for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) v[x] = x ? EVBDD_ZERO : EVBDD_ONE;
        bool diverged = false;
        unsigned log_k[128], log_g[128], log_a[128], log_b[128];
        for (unsigned gi = 0; gi < 12 * n && !diverged; gi++) {
            static const uint32_t one[] = { GATEID_H, GATEID_H, GATEID_X, GATEID_Y, GATEID_Z,
                                            GATEID_S, GATEID_Sdag, GATEID_T, GATEID_Tdag };
            const unsigned kind = (n >= 2) ? (unsigned)rnd_below(4) : 0;
            const unsigned a = (unsigned)rnd_below(n);
            unsigned b = (n >= 2) ? (unsigned)rnd_below(n - 1) : 0; if (n >= 2 && b >= a) b++;
            log_k[gi] = kind; log_a[gi] = a; log_b[gi] = b; log_g[gi] = 0;
            if (kind <= 1) {
                const uint32_t g = one[rnd_below(sizeof(one) / sizeof(one[0]))];
                log_g[gi] = g;
                e = bqd_gate(e, g, a, n);
                dense_gate(v, n, g, 0, a);
            } else if (kind == 2) {
                const uint32_t g = rnd_below(2) ? GATEID_X : GATEID_Z;
                log_g[gi] = g;
                bool ok;
                e = bqd_cgate_either(e, g, a, b, n, &ok);
                dense_gate(v, n, g, UINT64_C(1) << a, b);
            } else {
                e = bqd_swap(e, a, b, n);
                dense_swap(v, n, a, b);
            }
            states++;
            if (!bqd_has_full_support(e)) not_full++;
            const size_t nodes = limdd_nodecount(e, n);
            if (nodes > max_nodes) max_nodes = nodes;
            if (!denotes(e, v, n)) {
                bad++; diverged = true;
                if (getenv("BQD_DEBUG")) {
                    fprintf(stderr, "  circuit on %u qubits:", n);
                    for (unsigned j = 0; j <= gi; j++)
                        fprintf(stderr, " %u:%u:%u:%u", log_k[j], log_g[j], log_a[j], log_b[j]);
                    fprintf(stderr, "\n");
                }
            }
        }
    }
    printf("  %u states on 1..%u qubits, %u without full support, at most %zu nodes\n",
           states, max_n, not_full, max_nodes);
    report("circuits, every intermediate state", bad, states);
    free(v);
}

/* --- harness -------------------------------------------------------------- */

TASK_0(int, run_exact)
{
    bqd_init(BQD_FAMILY_SCALAR, NQ, 1LL << 20, 1LL << 20, 1LL << 22, 1LL << 20);
    pw[0] = EVBDD_ONE;
    const EVBDD_WGT w = qisq2_lookup(0, 1, 1, 2, 0, 1, 1, 2);
    for (int k = 1; k < 8; k++) pw[k] = wgt_mul(pw[k - 1], w);
    const int before = failures;
    check_basis_and_compose();
    check_pointwise();
    check_gates();
    check_controlled();
    check_circuits(8, 12);
    bqd_quit();
    return failures != before;
}

int
main(void)
{
    printf("== scalar BQD, exact weights ==\n");
    lace_start(4, 0);
    sylvan_set_sizes(1LL << 20, 1LL << 20, 1LL << 20, 1LL << 20);
    sylvan_init_package();
    qsylvan_init_simulator(1LL << 22, 1LL << 22, -1, QISQ2_MAP, NORM_LOW);
    const int res = RUN(run_exact);
    sylvan_quit();
    lace_stop();
    return res;
}
