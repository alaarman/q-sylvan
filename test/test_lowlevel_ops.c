/*
 * The low-level operations of the three diagrams, held to one dense oracle:
 * the EVDD (qsylvan_evdd_ops.h), the LIMDD (qsylvan_limdd_ops.h) and the BQD
 * (qsylvan_bqd_gates.h) must give the same vector for
 *
 *   times, plus, scale, negate
 *   restrict   f|_{x_q = b}, not depending on x_q
 *   project    f . [x_q = b]
 *   local_matvec  a dense 2^k x 2^k matrix on k qubits, k = 1, 2, 3, in any
 *              order of the qubits, with zeros and algebraic entries
 *
 * and for a single-qubit gate the local matvec must be what the diagram's own
 * gate gives, as an edge. Each result is also compared as an edge with the
 * diagram built from the oracle's vector, so a right value in the wrong shape
 * fails. Exact weights, so equality is exact.
 *
 * A LIMDD's recursions end at the width of its tables, so each width runs in
 * a session of its own; the BQD and the LIMDD share tables and run apart.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sylvan.h>
#include <sylvan_int.h>

#include "qsylvan_simulator.h"
#include "qsylvan_evdd_ops.h"
#include "qsylvan_limdd_ops.h"
#include "qsylvan_limdd_canon.h"
#include "qsylvan_bqd.h"
#include "qsylvan_bqd_gates.h"
#include "sylvan_edge_weights_qisq2.h"
#include "sylvan_evbdd_int.h"

#define MAXN 7
#define MAXV (1u << MAXN)

static int failures = 0;
static uint64_t rng_state = UINT64_C(0x10E1E7E1);
static uint64_t rnd(void)
{
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17;
    return rng_state;
}
static uint64_t rnd_below(uint64_t n) { return rnd() % n; }

typedef enum { EVDD = 0, LIMDD_K = 1, BQD_K = 2 } kind_t;
static const char *kind_name[3] = { "EVDD", "LIMDD", "BQD" };

/* --- vectors: qubit q is bit n-1-q of an index --------------------------- */

static EVBDD_WGT
rand_weight(void)
{
    long a, b, c, d;
    do {
        a = (long)rnd_below(5) - 2; b = (long)rnd_below(3) - 1;
        c = (long)rnd_below(5) - 2; d = (long)rnd_below(3) - 1;
    } while (a == 0 && b == 0 && c == 0 && d == 0);
    return qisq2_lookup(a, 1, b, 2, c, 1, d, 2);
}

static void
rand_vector(unsigned n, double pzero, EVBDD_WGT *f)
{
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
        f[x] = ((double)rnd_below(1000) / 1000.0 < pzero) ? EVBDD_ZERO : rand_weight();
}

/* --- the three diagrams, built from and read back into vectors ----------- */

static EVBDD
ev_build(const EVBDD_WGT *f, uint32_t var, uint32_t n)
{
    if (var == n) return evbdd_bundle(EVBDD_TERMINAL, f[0]);
    const uint64_t h = UINT64_C(1) << (n - var - 1);
    return evbdd_makenode(var, ev_build(f, var + 1, n), ev_build(f + h, var + 1, n));
}

static LIMDD
li_build(const EVBDD_WGT *f, uint32_t var, uint32_t n)
{
    if (var == n) {
        if (f[0] == EVBDD_ZERO) return limdd_zero_edge();
        return limdd_bundle(limdd_lim_make(limdd_pauli_identity(), f[0]), LIMDD_TERMINAL);
    }
    const uint64_t h = UINT64_C(1) << (n - var - 1);
    return limdd_makeedge(var, li_build(f, var + 1, n), li_build(f + h, var + 1, n));
}

static uint64_t
build(kind_t k, const EVBDD_WGT *f, unsigned n)
{
    switch (k) {
    case EVDD:    return ev_build(f, 0, n);
    case LIMDD_K: return li_build(f, 0, n);
    default:      return bqd_from_vector(f, n);
    }
}

static void
read_back(kind_t k, uint64_t e, unsigned n, EVBDD_WGT *out)
{
    if (k == BQD_K) { bqd_to_vector(e, n, out); return; }
    bool bits[MAXN];
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) {
        for (unsigned q = 0; q < n; q++) bits[q] = (x >> (n - 1 - q)) & 1;
        out[x] = (k == EVDD) ? evbdd_getvalue(e, bits) : limdd_eval(e, bits, n);
    }
}

/* --- the operations, one name per diagram -------------------------------- */

static uint64_t op_times(kind_t k, uint64_t a, uint64_t b)
{ return k == EVDD ? evbdd_times(a, b) : k == LIMDD_K ? limdd_times(a, b, 0) : bqd_multiply(a, b); }
static uint64_t op_plus(kind_t k, uint64_t a, uint64_t b)
{ return k == EVDD ? evbdd_plus(a, b) : k == LIMDD_K ? limdd_plus(a, b, 0) : bqd_add(a, b); }
static uint64_t op_scale(kind_t k, uint64_t a, EVBDD_WGT c)
{ return k == EVDD ? evbdd_scale(a, c) : k == LIMDD_K ? limdd_scale(a, c) : bqd_scale(a, c); }
static uint64_t op_negate(kind_t k, uint64_t a)
{ return k == EVDD ? evbdd_negate(a) : k == LIMDD_K ? limdd_negate(a) : bqd_negate(a); }
static uint64_t op_restrict(kind_t k, uint64_t a, unsigned q, int b)
{ return k == EVDD ? evbdd_restrict(a, q, b) : k == LIMDD_K ? limdd_restrict(a, q, b, 0) : bqd_restrict(a, q, b); }
static uint64_t op_project(kind_t k, uint64_t a, unsigned q, int b)
{ return k == EVDD ? evbdd_project(a, q, b) : k == LIMDD_K ? limdd_project(a, q, b, 0) : bqd_project(a, q, b); }
static uint64_t op_matvec(kind_t k, uint64_t a, const EVBDD_WGT *M, const uint32_t *qs, uint32_t kk, unsigned n)
{
    return k == EVDD ? evbdd_local_matvec(a, M, qs, kk, n)
         : k == LIMDD_K ? limdd_local_matvec(a, M, qs, kk, n)
         : bqd_local_matvec(a, M, qs, kk, n);
}
static uint64_t op_gate(kind_t k, uint64_t a, uint32_t gid, unsigned q, unsigned n)
{
    return k == EVDD ? qmdd_gate(a, gid, q) : k == LIMDD_K ? limdd_gate(a, gid, q, n)
         : bqd_gate(a, gid, q, n);
}

/* --- checking ------------------------------------------------------------- */

typedef struct { unsigned tot, bad_value, bad_edge; } tally_t;
enum { T_TIMES, T_PLUS, T_SCALE, T_NEGATE, T_RESTRICT, T_PROJECT, T_MATVEC, T_GATE, N_T };
static const char *tname[N_T] = { "times", "plus", "scale", "negate", "restrict",
                                  "project", "local_matvec", "matvec = gate" };
static tally_t tally[3][N_T];

static void
check(kind_t k, int t, uint64_t e, const EVBDD_WGT *want, unsigned n)
{
    EVBDD_WGT got[MAXV];
    read_back(k, e, n, got);
    tally[k][t].tot++;
    if (memcmp(got, want, sizeof(EVBDD_WGT) << n) != 0) tally[k][t].bad_value++;
    else if (build(k, want, n) != e) {
        tally[k][t].bad_edge++;
        if (getenv("LL_DEBUG")) {
            const uint64_t b2 = build(k, want, n);
            fprintf(stderr, "  %s %s n=%u: got target %llu label %llu, built target %llu label %llu; nonzeros",
                    kind_name[k], tname[t], n,
                    (unsigned long long)limdd_target(e), (unsigned long long)limdd_label(e),
                    (unsigned long long)limdd_target(b2), (unsigned long long)limdd_label(b2));
            unsigned nz = 0; for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) nz += want[x] != EVBDD_ZERO;
            fprintf(stderr, " %u; built again equal: %d\n", nz, build(k, want, n) == b2);
        }
    }
}

/** out[x] = sum_c M[r][c] v[x with Q := c], r = x on Q */
static void
dense_matvec(const EVBDD_WGT *v, const EVBDD_WGT *M, const uint32_t *qs, uint32_t kk,
             unsigned n, EVBDD_WGT *out)
{
    const uint32_t dim = 1u << kk;
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) {
        uint32_t r = 0;
        uint64_t base = x;
        for (uint32_t i = 0; i < kk; i++) {
            const uint64_t bit = UINT64_C(1) << (n - 1 - qs[i]);
            if (x & bit) r |= 1u << (kk - 1 - i);
            base &= ~bit;
        }
        EVBDD_WGT s = EVBDD_ZERO;
        for (uint32_t c = 0; c < dim; c++) {
            uint64_t y = base;
            for (uint32_t i = 0; i < kk; i++)
                if ((c >> (kk - 1 - i)) & 1) y |= UINT64_C(1) << (n - 1 - qs[i]);
            s = wgt_add(s, wgt_mul(M[(size_t)r * dim + c], v[y]));
        }
        out[x] = s;
    }
}

static uint64_t op_gate(kind_t k, uint64_t a, uint32_t gid, unsigned q, unsigned n);

/** |0..0> and 4n random gates from H, T, S, X and CX, either way round. */
static uint64_t
circuit_state(kind_t k, unsigned n)
{
    EVBDD_WGT z[MAXV];
    memset(z, 0, sizeof(z));
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) z[x] = x ? EVBDD_ZERO : EVBDD_ONE;
    uint64_t e = build(k, z, n);
    static const uint32_t one[] = { GATEID_H, GATEID_H, GATEID_T, GATEID_S, GATEID_X };
    for (unsigned i = 0; i < 4 * n; i++) {
        const uint32_t a = (uint32_t)rnd_below(n);
        if (n > 1 && rnd_below(3) == 0) {
            uint32_t b = (uint32_t)rnd_below(n - 1); if (b >= a) b++;
            bool ok;
            e = k == EVDD ? qmdd_cgate(e, GATEID_X, a, b, n)
              : k == LIMDD_K ? limdd_cgate_either(e, GATEID_X, a, b, n, &ok)
              : bqd_cgate_either(e, GATEID_X, a, b, n, &ok);
        } else {
            e = op_gate(k, e, one[rnd_below(5)], a, n);
        }
    }
    return e;
}

static void
run_checks(kind_t k, unsigned n, int reps)
{
    EVBDD_WGT f[MAXV], g[MAXV], want[MAXV];
    const uint64_t len = UINT64_C(1) << n;
    for (int rep = 0; rep < reps; rep++) {
        const double pz = (double)(rep % 3) * 0.3;
        uint64_t F, G;
        if (rep & 1) {
            /* states of random Clifford+T circuits, whose stabilisers are
             * rich where random vectors have none: that is where a LIMDD
             * label can be left unreduced */
            F = circuit_state(k, n);
            G = circuit_state(k, n);
            read_back(k, F, n, f);
            read_back(k, G, n, g);
        } else {
            rand_vector(n, pz, f);
            rand_vector(n, pz, g);
            F = build(k, f, n);
            G = build(k, g, n);
        }

        for (uint64_t x = 0; x < len; x++) want[x] = wgt_mul(f[x], g[x]);
        check(k, T_TIMES, op_times(k, F, G), want, n);
        for (uint64_t x = 0; x < len; x++) want[x] = wgt_add(f[x], g[x]);
        check(k, T_PLUS, op_plus(k, F, G), want, n);
        const EVBDD_WGT c = rand_weight();
        for (uint64_t x = 0; x < len; x++) want[x] = wgt_mul(c, f[x]);
        check(k, T_SCALE, op_scale(k, F, c), want, n);
        for (uint64_t x = 0; x < len; x++) want[x] = wgt_neg(f[x]);
        check(k, T_NEGATE, op_negate(k, F), want, n);

        for (unsigned q = 0; q < n; q++) for (int b = 0; b < 2; b++) {
            const uint64_t qb = UINT64_C(1) << (n - 1 - q);
            for (uint64_t x = 0; x < len; x++) want[x] = f[b ? (x | qb) : (x & ~qb)];
            check(k, T_RESTRICT, op_restrict(k, F, q, b), want, n);
            for (uint64_t x = 0; x < len; x++) want[x] = (((x & qb) != 0) == b) ? f[x] : EVBDD_ZERO;
            check(k, T_PROJECT, op_project(k, F, q, b), want, n);
        }

        /* a dense matrix on k random qubits in a random order */
        const uint32_t kk = 1 + (uint32_t)rnd_below(n < 3 ? n : 3);
        uint32_t qs[3];
        for (uint32_t i = 0; i < kk; i++) {
            bool again;
            do {
                qs[i] = (uint32_t)rnd_below(n);
                again = false;
                for (uint32_t j = 0; j < i; j++) if (qs[j] == qs[i]) again = true;
            } while (again);
        }
        EVBDD_WGT M[64];
        for (uint32_t i = 0; i < (1u << (2 * kk)); i++) M[i] = rnd_below(3) ? rand_weight() : EVBDD_ZERO;
        dense_matvec(f, M, qs, kk, n, want);
        check(k, T_MATVEC, op_matvec(k, F, M, qs, kk, n), want, n);

        /* for one qubit the matvec is the gate */
        static const uint32_t gs[] = { GATEID_H, GATEID_X, GATEID_Y, GATEID_S, GATEID_T, GATEID_sqrtX };
        const uint32_t gid = gs[rnd_below(6)];
        const uint32_t q = (uint32_t)rnd_below(n);
        const EVBDD_WGT U[4] = { gates[gid][0], gates[gid][1], gates[gid][2], gates[gid][3] };
        tally[k][T_GATE].tot++;
        if (op_matvec(k, F, U, &q, 1, n) != op_gate(k, F, gid, q, n)) tally[k][T_GATE].bad_edge++;
    }
}

static void
report(kind_t k)
{
    for (int t = 0; t < N_T; t++) {
        const tally_t *s = &tally[k][t];
        if (s->tot == 0) continue;
        const bool ok = s->bad_value == 0 && s->bad_edge == 0;
        printf("  %-6s %-14s %s (%u", kind_name[k], tname[t], ok ? "ok" : "FAILED", s->tot);
        if (!ok) printf("; %u wrong value, %u right value in another shape", s->bad_value, s->bad_edge);
        printf(")\n");
        if (!ok) failures++;
    }
}

/* --- sessions -------------------------------------------------------------- */

static unsigned sess_n;

TASK_0(int, run_ev_li)
{
    limdd_nodes_init(sess_n, 1LL << 20, 1LL << 20, 1LL << 22, 1LL << 20);
    limdd_set_canon_policy(LIMDD_CANON_ALWAYS, 0);
    run_checks(EVDD, sess_n, 30);
    run_checks(LIMDD_K, sess_n, 30);
    limdd_nodes_quit();
    return 0;
}

TASK_0(int, run_bqd)
{
    bqd_init(BQD_FAMILY_SCALAR, MAXN, 1LL << 20, 1LL << 20, 1LL << 22, 1LL << 20);
    for (unsigned n = 1; n <= MAXN; n++) run_checks(BQD_K, n, 30);
    bqd_quit();
    return 0;
}

static void
session(int (*go)(void))
{
    lace_start(getenv("LL_WORKERS") ? (unsigned)atoi(getenv("LL_WORKERS")) : 4, 0);
    sylvan_set_sizes(1LL << 20, 1LL << 20, 1LL << 20, 1LL << 20);
    sylvan_init_package();
    qsylvan_init_simulator(1LL << 22, 1LL << 22, -1, QISQ2_MAP, NORM_LOW);
    go();
    sylvan_quit();
    lace_stop();
}

static int go_ev_li(void) { return RUN(run_ev_li); }
static int go_bqd(void) { return RUN(run_bqd); }

int
main(void)
{
    static const unsigned widths[] = { 1, 2, 3, 5, 7 };
    for (unsigned i = 0; i < sizeof(widths) / sizeof(widths[0]); i++) {
        sess_n = widths[i];
        session(go_ev_li);
    }
    session(go_bqd);
    printf("exact weights, EVDD and LIMDD at widths 1, 2, 3, 5, 7, BQD at 1..7:\n");
    report(EVDD);
    report(LIMDD_K);
    report(BQD_K);
    return failures != 0;
}
