/*
 * The binary quotient diagram, in its three label families, checked against
 * what the paper proves and what its prototype measured.
 *
 *   round trip        decode(build(f)) = f, and the amplitude query agrees,
 *                     on random vectors at several zero densities
 *   canonicity        building a vector twice gives the same edge, and the
 *                     representative is idempotent
 *   thm:coset         on a level-k coset state the X-BQD has at most
 *                     sum_{i<k} C(v, i) + 2 nodes at the level with v decided
 *                     variables, in random orders, and the copy clause never
 *                     fires there
 *   thm:pcoset        the Pauli-BQD has at most sum_{i<k-1} C(v, i) + 2, so at
 *                     most 3 per level on every stabiliser state
 *   def:prep          the Pauli representative is invariant under scalars
 *                     and sign patterns wherever the pivots span the support,
 *                     which is every node of a coset state
 *   full support      the X-BQD carries no translation, so it is the BQD
 *   scalar labels     on a coset state with a proper support the copy fires,
 *                     which is the degradation the translation label removes
 *
 * Everything is exact, in Q(w_8): coset states have amplitudes that are
 * eighth roots of unity or zero, which is what the paper's prototype stores.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sylvan.h>
#include <sylvan_int.h>

#include "qsylvan_bqd.h"
#include "qsylvan_limdd_inspect.h"
#include "qsylvan_simulator.h"
#include "sylvan_edge_weights_qisq2.h"

#define NQ 12
#define MAXV (1u << NQ)

static int failures = 0;
static EVBDD_WGT pw[8];                    /* w_8^e */

static uint64_t rng_state = UINT64_C(0x9E3779B97F4A7C15);
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

/* --- states --------------------------------------------------------------- */

/** A random linear code of dimension k in F_2^n, as a membership table. */
static void
rand_code(unsigned n, unsigned k, uint8_t *in_code)
{
    for (;;) {
        memset(in_code, 0, (size_t)1 << n);
        in_code[0] = 1;
        uint64_t size = 1;
        for (unsigned i = 0; i < k; i++) {
            const uint64_t b = 1 + rnd_below((UINT64_C(1) << n) - 1);
            if (in_code[b]) continue;
            for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
                if (in_code[x] && !in_code[x ^ b]) { in_code[x ^ b] = 2; }
            for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
                if (in_code[x] == 2) { in_code[x] = 1; size++; }
        }
        if (size == (UINT64_C(1) << k)) return;
    }
}

typedef struct { uint64_t W; unsigned c; } term_t;

/**
 * A level-k phase polynomial over Z_8, as the prototype draws it: the
 * degree-d coefficient is a multiple of 2^{3-k+d-1}, so a level-2 polynomial
 * has even linear and quadratic-in-4 coefficients, and a level-3 one is any
 * cubic with quadratic coefficients even and cubic ones in 4.
 */
static unsigned
rand_level(unsigned n, unsigned k, term_t *terms)
{
    unsigned nt = 0;
    for (unsigned d = 1; d <= k; d++) {
        const unsigned step = 1u << (3 - k + d - 1);
        if (step > 8) continue;
        /* every subset of size d */
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

/** The level-k coset state  [x in a + C] w_8^{P(x)}. */
static void
coset_state(unsigned n, const uint8_t *in_code, uint64_t a, const term_t *terms, unsigned nt,
            EVBDD_WGT *f)
{
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
        f[x] = in_code[x ^ a] ? pw[poly_val(terms, nt, x)] : EVBDD_ZERO;
}

/** Read f in another variable order: bit k of the new index is bit perm[k] of the old. */
static void
reorder(const EVBDD_WGT *f, unsigned n, const unsigned *perm, EVBDD_WGT *out)
{
    for (uint64_t idx = 0; idx < (UINT64_C(1) << n); idx++) {
        uint64_t o = 0;
        for (unsigned k = 0; k < n; k++) if ((idx >> k) & 1) o |= UINT64_C(1) << perm[k];
        out[idx] = f[o];
    }
}

static void
rand_perm(unsigned n, unsigned *perm)
{
    for (unsigned i = 0; i < n; i++) perm[i] = i;
    for (unsigned i = n; i > 1; i--) {
        const unsigned j = (unsigned)rnd_below(i);
        const unsigned t = perm[i - 1]; perm[i - 1] = perm[j]; perm[j] = t;
    }
}

static void
rand_vector(unsigned n, double pzero, EVBDD_WGT *f)
{
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
        f[x] = ((double)rnd_below(1000) / 1000.0 < pzero) ? EVBDD_ZERO : pw[rnd_below(8)];
}

/* --- diagram walks used by the checks ------------------------------------- */

/**
 * The support of the function a node denotes, from the diagram alone:
 * supp(g_0 (.) g_1) = supp(g_1) (lem:shadow), so the high half is the
 * translated support of the high child.
 */
static void
support(LIMDD_TARG p, unsigned var, unsigned n, uint8_t *bits)
{
    const uint64_t len = UINT64_C(1) << (n - var);
    if (var == n) { bits[0] = 1; return; }
    const uint64_t h = len >> 1;
    const LIMDD low = limdd_node_low(p), high = limdd_node_high(p);
    if (limdd_edge_is_zero(low)) memset(bits, 0, h);
    else support(limdd_target(low), var + 1, n, bits);
    if (limdd_edge_is_zero(high)) { memset(bits + h, 0, h); return; }
    uint8_t *g1 = malloc(h);
    support(limdd_target(high), var + 1, n, g1);
    EVBDD_WGT c; uint64_t s, t;
    bqd_lim_masks(limdd_label(high), n, &c, &s, &t);
    for (uint64_t y = 0; y < h; y++) bits[h + y] = g1[y ^ t];
    free(g1);
}

/** Points where the copy clause fires: g_0 zero and g_1 nonzero. Over every path. */
static uint64_t
copies_below(LIMDD_TARG p, unsigned var, unsigned n)
{
    if (var == n) return 0;
    const uint64_t h = UINT64_C(1) << (n - var - 1);
    const LIMDD low = limdd_node_low(p), high = limdd_node_high(p);
    uint64_t fires = 0;
    if (!limdd_edge_is_zero(high)) {
        uint8_t *s0 = calloc(h, 1), *s1 = calloc(h, 1);
        if (!limdd_edge_is_zero(low)) support(limdd_target(low), var + 1, n, s0);
        support(limdd_target(high), var + 1, n, s1);
        for (uint64_t y = 0; y < h; y++) if (s1[y] && !s0[y]) fires++;
        free(s0); free(s1);
        fires += copies_below(limdd_target(high), var + 1, n);
    }
    if (!limdd_edge_is_zero(low)) fires += copies_below(limdd_target(low), var + 1, n);
    return fires;
}

/** Whether any high label below carries a translation. */
static bool
any_translation(LIMDD_TARG p, unsigned var, unsigned n)
{
    if (var == n) return false;
    const LIMDD low = limdd_node_low(p), high = limdd_node_high(p);
    if (!limdd_edge_is_zero(high)) {
        EVBDD_WGT c; uint64_t s, t;
        bqd_lim_masks(limdd_label(high), n, &c, &s, &t);
        if (t != 0) return true;
        if (any_translation(limdd_target(high), var + 1, n)) return true;
    }
    return !limdd_edge_is_zero(low) && any_translation(limdd_target(low), var + 1, n);
}

/* --- checks --------------------------------------------------------------- */

static void
check_round_trip(void)
{
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *g = malloc(MAXV * sizeof(EVBDD_WGT));
    unsigned bad = 0, badeval = 0, badcanon = 0, tot = 0;
    const double dens[4] = { 0.0, 0.2, 0.5, 0.8 };
    for (unsigned n = 1; n <= 8; n++) {
        for (unsigned di = 0; di < 4; di++) {
            for (int rep = 0; rep < 40; rep++) {
                rand_vector(n, dens[di], f);
                const BQD e = bqd_from_vector(f, n);
                bqd_to_vector(e, n, g);
                tot++;
                if (memcmp(f, g, ((size_t)1 << n) * sizeof(EVBDD_WGT)) != 0) bad++;
                for (int q = 0; q < 8; q++) {
                    const uint64_t x = rnd_below(UINT64_C(1) << n);
                    if (bqd_eval(e, n, x) != f[x]) { badeval++; break; }
                }
                if (bqd_from_vector(f, n) != e) badcanon++;
            }
        }
    }
    char buf[128];
    snprintf(buf, sizeof(buf), "%u of %u vectors do not decode to themselves", bad, tot);
    expect(bad == 0, "round trip", buf);
    snprintf(buf, sizeof(buf), "%u of %u vectors have an amplitude the query gets wrong", badeval, tot);
    expect(badeval == 0, "amplitude query", buf);
    snprintf(buf, sizeof(buf), "%u of %u vectors built twice gave two different edges", badcanon, tot);
    expect(badcanon == 0, "canonicity", buf);
    free(f); free(g);
}

static void
check_representative_idempotent(void)
{
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *r = malloc(MAXV * sizeof(EVBDD_WGT)),
              *r2 = malloc(MAXV * sizeof(EVBDD_WGT));
    unsigned bad = 0, tot = 0;
    for (unsigned n = 1; n <= 8; n++) for (int rep = 0; rep < 60; rep++) {
        rand_vector(n, 0.3, f);
        const uint64_t len = UINT64_C(1) << n;
        EVBDD_WGT c; uint64_t s, p;
        bqd_representative(f, len, r, &c, &s, &p);
        bqd_representative(r, len, r2, &c, &s, &p);
        tot++;
        if (memcmp(r, r2, len * sizeof(EVBDD_WGT)) != 0 || c != EVBDD_ONE || s != 0 || p != 0) bad++;
    }
    char buf[96];
    snprintf(buf, sizeof(buf), "%u of %u representatives are not fixed by taking the representative again", bad, tot);
    expect(bad == 0, "representative idempotent", buf);
    free(f); free(r); free(r2);
}

static void
check_coset_states(void)
{
    const bqd_family_t fam = bqd_family();
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *g = malloc(MAXV * sizeof(EVBDD_WGT));
    uint8_t *code = malloc(MAXV);
    term_t *terms = malloc(4096 * sizeof(term_t));
    unsigned perm[NQ];
    size_t counts[NQ];

    unsigned states = 0, bound_bad = 0, copy_bad = 0, copies_seen = 0, trans_bad = 0, prep_bad = 0;
    size_t max_width = 0;

    for (unsigned n = 3; n <= 8; n++) {
        for (unsigned k = 2; k <= 3; k++) {
            for (int rep = 0; rep < 25; rep++) {
                const unsigned dim = (unsigned)rnd_below(n + 1);        /* 0..n */
                rand_code(n, dim, code);
                const uint64_t a = rnd_below(UINT64_C(1) << n);
                const unsigned nt = rand_level(n, k, terms);
                coset_state(n, code, a, terms, nt, g);
                rand_perm(n, perm);
                reorder(g, n, perm, f);

                const BQD e = bqd_from_vector(f, n);
                states++;
                limdd_level_counts(e, counts, n);
                for (unsigned v = 0; v < n; v++) if (counts[v] > max_width) max_width = counts[v];

                if (fam == BQD_FAMILY_X || fam == BQD_FAMILY_PAULI) {
                    /* thm:coset / thm:pcoset, per level: v decided variables */
                    const unsigned kk = (fam == BQD_FAMILY_X) ? k : k - 1;
                    for (unsigned v = 0; v < n; v++) {
                        uint64_t bound = 2;
                        for (unsigned i = 0; i < kk; i++) bound += binom(v, i);
                        if (counts[v] > bound) {
                            bound_bad++;
                            fprintf(stderr, "  %s: n=%u k=%u level v=%u has %zu nodes, bound %llu\n",
                                    bqd_family_name(fam), n, k, v, counts[v],
                                    (unsigned long long)bound);
                            break;
                        }
                    }
                    /* the copy never fires on a coset state */
                    if (!limdd_edge_is_zero(e) && copies_below(limdd_target(e), 0, n) != 0) copy_bad++;
                } else {
                    if (!limdd_edge_is_zero(e) && copies_below(limdd_target(e), 0, n) != 0) copies_seen++;
                }

                /* full support: the X-BQD is the BQD, no translation anywhere */
                if (fam == BQD_FAMILY_X && dim == n && !limdd_edge_is_zero(e)
                    && any_translation(limdd_target(e), 0, n)) trans_bad++;

                /* def:prep is an invariant under scalars and sign patterns
                 * wherever the pivots span the support, i.e. on a coset state */
                if (fam == BQD_FAMILY_PAULI) {
                    const uint64_t len = UINT64_C(1) << n;
                    const EVBDD_WGT c = pw[rnd_below(8)];
                    const uint64_t s = rnd_below(len);
                    EVBDD_WGT *h = malloc(len * sizeof(EVBDD_WGT));
                    for (uint64_t y = 0; y < len; y++) {
                        EVBDD_WGT v = f[y];
                        if (v != EVBDD_ZERO) {
                            if (__builtin_parityll(s & y)) v = wgt_neg(v);
                            v = wgt_mul(c, v);
                        }
                        h[y] = v;
                    }
                    if (bqd_from_vector(h, n) != e && limdd_target(bqd_from_vector(h, n)) != limdd_target(e))
                        prep_bad++;
                    free(h);
                }
            }
        }
    }

    char buf[160];
    printf("  %u coset states (n = 3..8, levels 2 and 3, random orders), max width %zu\n",
           states, max_width);
    if (fam == BQD_FAMILY_X || fam == BQD_FAMILY_PAULI) {
        snprintf(buf, sizeof(buf), "%u of %u states exceed the per-level bound", bound_bad, states);
        expect(bound_bad == 0, fam == BQD_FAMILY_X ? "thm:coset" : "thm:pcoset", buf);
        snprintf(buf, sizeof(buf), "the copy clause fired on %u of %u coset states", copy_bad, states);
        expect(copy_bad == 0, "no copy on a coset state", buf);
    } else {
        snprintf(buf, sizeof(buf), "the copy never fired on any of %u coset states; the scalar "
                 "family should degrade on a proper support", states);
        expect(copies_seen > 0, "scalar labels degrade", buf);
        printf("  copy fired on %u of %u states (expected: it does, without a translation)\n",
               copies_seen, states);
    }
    if (fam == BQD_FAMILY_X) {
        snprintf(buf, sizeof(buf), "%u full-support states carried a translation", trans_bad);
        expect(trans_bad == 0, "full support is translation-free", buf);
    }
    if (fam == BQD_FAMILY_PAULI) {
        snprintf(buf, sizeof(buf), "%u of %u: (c Z^s) f and f got different nodes", prep_bad, states);
        expect(prep_bad == 0, "def:prep invariant under scalars and Z", buf);
    }
    free(f); free(g); free(code); free(terms);
}

/* --- harness -------------------------------------------------------------- */

TASK_1(int, runtests, int, fam)
{
    bqd_init((bqd_family_t)fam, NQ, 1LL << 20, 1LL << 20, 1LL << 22, 1LL << 20);

    /* w_8 = (1 + i)/sqrt2 = sqrt2/2 + i sqrt2/2, and its powers by multiplication */
    pw[0] = EVBDD_ONE;
    const EVBDD_WGT w = qisq2_lookup(0, 1, 1, 2, 0, 1, 1, 2);
    for (int e = 1; e < 8; e++) pw[e] = wgt_mul(pw[e - 1], w);
    if (pw[4] != EVBDD_MIN_ONE) { fprintf(stderr, "w_8^4 is not -1\n"); return 1; }

    const int before = failures;
    check_round_trip();
    check_representative_idempotent();
    check_coset_states();

    bqd_quit();
    return failures != before;
}

static int
run_family(bqd_family_t fam)
{
    printf("== %s ==\n", bqd_family_name(fam));
    lace_start(4, 0);
    sylvan_set_sizes(1LL << 20, 1LL << 20, 1LL << 20, 1LL << 20);
    sylvan_init_package();
    qsylvan_init_simulator(1LL << 22, 1LL << 22, -1, QISQ2_MAP, NORM_LOW);

    const int res = RUN(runtests, (int)fam);

    sylvan_quit();
    lace_stop();
    printf("  %s\n", res ? "FAILED" : "ok");
    return res;
}

int
main(void)
{
    int bad = 0;
    bad |= run_family(BQD_FAMILY_SCALAR);
    bad |= run_family(BQD_FAMILY_X);
    bad |= run_family(BQD_FAMILY_PAULI);
    return bad;
}
