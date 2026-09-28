/*
 * The binary quotient diagram, in its three label families, checked against
 * what the paper proves and what its prototype measured.
 *
 *   round trip        decode(build(f)) = f, and the amplitude query agrees,
 *                     on random vectors at several zero densities
 *   determinism       building a vector twice, on four workers, gives the
 *                     same edge
 *   def:rep           c.f and f have the same root node, in every family:
 *                     the representative divides the scale out, so the
 *                     orbit under scalars is one node; and taking the
 *                     representative again changes nothing
 *   thm:coset         on a level-k coset state the X-BQD has at most
 *                     sum_{i<k} C(v, i) + 1 nodes at the level with v decided
 *                     variables, the node 0 counted where it occurs as the
 *                     paper counts it, in random orders, and the copy clause
 *                     never fires there
 *   thm:pcoset        the Pauli-BQD has at most sum_{i<k-1} C(v, i) + 1
 *                     nodes that are not 0 at such a level. The theorem says
 *                     + 2, and one of the two is the node 0 (bqd2.tex, after
 *                     thm:pcoset), which this port keeps as an edge and never
 *                     counts. At k = 2 that is two per level.
 *   def:prep          the Pauli representative is invariant under scalars
 *                     and sign patterns wherever the pivots span the support,
 *                     which is every node of a coset state; and, on vectors
 *                     over Z[sqrt2, i] with both signs in a component, every
 *                     pivot value of the representative has its argument in
 *                     [0, pi), decided exactly (A + B sqrt2 with A, B of
 *                     opposite sign is the case sign_sqrt2 has to compare)
 *   full support      the X-BQD carries no translation, so it is the BQD
 *   scalar labels     on a coset state with a proper support the copy fires,
 *                     which is the degradation the translation label removes
 *   node counts       limdd_countnodes, the count the qasm runner takes after
 *                     every gate with -c, is limdd_nodecount on vectors of 9
 *                     to 12 qubits, whose diagrams outgrow the set of met
 *                     nodes it starts with
 *   float weights     the Pauli family on the float backend, whose sign
 *                     decision is a comparison of doubles: round trip,
 *                     the pivot arguments, and def:prep invariance
 *   skip:def:fr       levels are skipped: every stored node depends on its
 *                     variable (skip:thm:canon (ii)), on all of the vectors
 *                     above and on vectors that do not depend on random
 *                     qubits, which round-trip; |+>^n is the terminal and a
 *                     basis state n nodes; and the note's examples of a
 *                     skipping edge with a label, a Z in the Pauli family and
 *                     a scalar in the X family, come out as it says
 *                     (skip:lem:normal)
 *
 * Everything but the last is exact, in Q(w_8, sqrt2): coset states have
 * amplitudes that are eighth roots of unity or zero, which is what the
 * paper's prototype stores.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sylvan.h>
#include <sylvan_int.h>

#include "qsylvan_bqd.h"
#include "qsylvan_limdd_inspect.h"
#include "qsylvan_limdd_ops.h"
#include "qsylvan_simulator.h"
#include "sylvan_edge_weights_complex.h"
#include "sylvan_edge_weights_qisq2.h"
#include "edge_weight_storage/qisq2_map.h"

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

static int64_t
lexmin_nonzero(const EVBDD_WGT *f, uint64_t len)
{
    for (uint64_t x = 0; x < len; x++) if (f[x] != EVBDD_ZERO) return (int64_t)x;
    return -1;
}

/* --- diagram walks used by the checks ------------------------------------- */

/** The level a node sits at: its variable, or n for the terminal of this vector. */
static unsigned
level_of(LIMDD_TARG p, unsigned n)
{
    return p == LIMDD_TERMINAL ? n : limdd_node_var(p);
}

/**
 * The support of the function a node denotes read at level var, from the
 * diagram alone: supp(g_0 (.) g_1) = supp(g_1) (lem:shadow), so the high half
 * is the translated support of the high child. An edge that skips levels
 * denotes its node's function extended to them, so the support repeats.
 */
static void
support(LIMDD_TARG p, unsigned var, unsigned n, uint8_t *bits)
{
    const unsigned v = level_of(p, n);
    const uint64_t len = UINT64_C(1) << (n - v);
    if (v == n) {
        bits[0] = 1;
    } else {
        const uint64_t h = len >> 1;
        const LIMDD low = limdd_node_low(p), high = limdd_node_high(p);
        if (limdd_edge_is_zero(low)) memset(bits, 0, h);
        else support(limdd_target(low), v + 1, n, bits);
        if (limdd_edge_is_zero(high)) {
            memset(bits + h, 0, h);
        } else {
            uint8_t *g1 = malloc(h);
            support(limdd_target(high), v + 1, n, g1);
            EVBDD_WGT c; uint64_t s, t;
            bqd_lim_masks(limdd_label(high), n, &c, &s, &t);
            for (uint64_t y = 0; y < h; y++) bits[h + y] = g1[y ^ t];
            free(g1);
        }
    }
    for (uint64_t y = len; y < (UINT64_C(1) << (n - var)); y++) bits[y] = bits[y - len];
}

/** Points where the copy clause fires: g_0 zero and g_1 nonzero. Over every path. */
static uint64_t
copies_below(LIMDD_TARG p, unsigned n)
{
    if (p == LIMDD_TERMINAL) return 0;
    const unsigned var = limdd_node_var(p);
    const uint64_t h = UINT64_C(1) << (n - var - 1);
    const LIMDD low = limdd_node_low(p), high = limdd_node_high(p);
    uint64_t fires = 0;
    if (!limdd_edge_is_zero(high)) {
        uint8_t *s0 = calloc(h, 1), *s1 = calloc(h, 1);
        if (!limdd_edge_is_zero(low)) support(limdd_target(low), var + 1, n, s0);
        support(limdd_target(high), var + 1, n, s1);
        for (uint64_t y = 0; y < h; y++) if (s1[y] && !s0[y]) fires++;
        free(s0); free(s1);
        fires += copies_below(limdd_target(high), n);
    }
    if (!limdd_edge_is_zero(low)) fires += copies_below(limdd_target(low), n);
    return fires;
}

/**
 * The least variable of a node with a zero child, or n if none. The paper
 * counts the function 0 as a node at every level it occurs (thm:coset), and
 * it occurs at every level below such a node, where this port keeps an edge.
 */
static unsigned
min_zero_var(LIMDD_TARG p, unsigned n)
{
    if (p == LIMDD_TERMINAL) return n;
    const LIMDD low = limdd_node_low(p), high = limdd_node_high(p);
    if (limdd_edge_is_zero(low) || limdd_edge_is_zero(high)) return limdd_node_var(p);
    const unsigned a = min_zero_var(limdd_target(low), n);
    const unsigned b = min_zero_var(limdd_target(high), n);
    return a < b ? a : b;
}

/** Whether any high label below carries a translation. */
static bool
any_translation(LIMDD_TARG p, unsigned n)
{
    if (p == LIMDD_TERMINAL) return false;
    const LIMDD low = limdd_node_low(p), high = limdd_node_high(p);
    if (!limdd_edge_is_zero(high)) {
        EVBDD_WGT c; uint64_t s, t;
        bqd_lim_masks(limdd_label(high), n, &c, &s, &t);
        if (t != 0) return true;
        if (any_translation(limdd_target(high), n)) return true;
    }
    return !limdd_edge_is_zero(low) && any_translation(limdd_target(low), n);
}

/**
 * The nodes reachable from p whose function does not depend on their own
 * variable, from the decoded vector of the edge to each: its two halves at
 * the node's bit are equal. A fully reduced diagram has none
 * (skip:thm:canon (ii)). Over every path, which is at most 2^n here.
 */
static unsigned
redundant_below(LIMDD_TARG p, unsigned n, EVBDD_WGT *buf)
{
    if (p == LIMDD_TERMINAL) return 0;
    bqd_to_vector(limdd_bundle(LIMDD_LIM_IDENTITY, p), n, buf);
    const uint64_t bit = UINT64_C(1) << (n - 1 - limdd_node_var(p));
    unsigned r = 1;
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
        if (!(x & bit) && buf[x] != buf[x | bit]) { r = 0; break; }
    const LIMDD low = limdd_node_low(p), high = limdd_node_high(p);
    if (!limdd_edge_is_zero(low)) r += redundant_below(limdd_target(low), n, buf);
    if (!limdd_edge_is_zero(high)) r += redundant_below(limdd_target(high), n, buf);
    return r;
}

/** Whether the edge into p, read at `var`, or an edge below p skips a level. */
static bool
skips_below(LIMDD_TARG p, unsigned var, unsigned n)
{
    const unsigned v = level_of(p, n);
    if (v > var) return true;
    if (p == LIMDD_TERMINAL) return false;
    const LIMDD low = limdd_node_low(p), high = limdd_node_high(p);
    if (!limdd_edge_is_zero(low) && skips_below(limdd_target(low), v + 1, n)) return true;
    return !limdd_edge_is_zero(high) && skips_below(limdd_target(high), v + 1, n);
}

/* --- checks --------------------------------------------------------------- */

static void
check_round_trip(void)
{
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *g = malloc(MAXV * sizeof(EVBDD_WGT));
    unsigned bad = 0, badeval = 0, badcanon = 0, badscale = 0, redundant = 0, tot = 0;
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
                if (!limdd_edge_is_zero(e)) redundant += redundant_below(limdd_target(e), n, g);

                /* def:rep: c.f is f up to the root's label, so one node */
                const EVBDD_WGT c = pw[1 + rnd_below(7)];
                for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
                    g[x] = (f[x] == EVBDD_ZERO) ? EVBDD_ZERO : wgt_mul(c, f[x]);
                const BQD ec = bqd_from_vector(g, n);
                if (limdd_target(ec) != limdd_target(e)) badscale++;
                else if (!limdd_edge_is_zero(e) && limdd_label(ec) == limdd_label(e)) badscale++;
            }
        }
    }
    char buf[128];
    snprintf(buf, sizeof(buf), "%u of %u vectors do not decode to themselves", bad, tot);
    expect(bad == 0, "round trip", buf);
    snprintf(buf, sizeof(buf), "%u of %u vectors have an amplitude the query gets wrong", badeval, tot);
    expect(badeval == 0, "amplitude query", buf);
    snprintf(buf, sizeof(buf), "%u of %u vectors built twice gave two different edges", badcanon, tot);
    expect(badcanon == 0, "determinism", buf);
    snprintf(buf, sizeof(buf), "%u of %u: c.f and f do not share the root node, or share its label", badscale, tot);
    expect(badscale == 0, "def:rep, scalar orbit", buf);
    snprintf(buf, sizeof(buf), "%u nodes of %u diagrams do not depend on their variable", redundant, tot);
    expect(redundant == 0, "skip:thm:canon (ii), random vectors", buf);
    free(f); free(g);
}

/**
 * limdd_countnodes against limdd_nodecount. The first keeps the nodes it has
 * met in a set that starts at 1024 slots and doubles at half full, so a
 * diagram of more than 512 nodes takes it through a doubling at least once.
 */
static void
check_countnodes(void)
{
    const uint64_t saved_rng = rng_state;
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT));
    unsigned bad = 0, big = 0, tot = 0;
    for (unsigned n = 9; n <= NQ; n++) for (int rep = 0; rep < 3; rep++) {
        rand_vector(n, (rep == 2) ? 0.3 : 0.0, f);
        const BQD e = bqd_from_vector(f, n);
        const size_t want = limdd_nodecount(e, n);
        tot++;
        if (want > 512) big++;
        if (limdd_countnodes(e) != want) bad++;
    }
    char buf[128];
    snprintf(buf, sizeof(buf), "%u of %u diagrams (%u of more than 512 nodes)", bad, tot, big);
    expect(bad == 0, "limdd_countnodes", buf);
    expect(big > 0, "limdd_countnodes", "no diagram of more than 512 nodes");
    free(f);
    rng_state = saved_rng;
}

/**
 * A vector that depends only on the qubits of a random set, then, where the
 * family has them, times a random sign pattern and translated, so that the
 * Pauli family meets a Z at a level it skips (skip:lem:normal).
 */
static void
rand_skipping(unsigned n, double pzero, EVBDD_WGT *f)
{
    const bqd_family_t fam = bqd_family();
    const uint64_t len = UINT64_C(1) << n;
    uint64_t D = 0;
    for (unsigned b = 0; b < n; b++) if (rnd_below(3) != 0) D |= UINT64_C(1) << b;
    const uint64_t s = (fam == BQD_FAMILY_PAULI && rnd_below(2)) ? rnd_below(len) : 0;
    const uint64_t t = (fam != BQD_FAMILY_SCALAR && rnd_below(2)) ? rnd_below(len) : 0;
    EVBDD_WGT *g = malloc(len * sizeof(EVBDD_WGT));
    rand_vector(n, pzero, g);
    for (uint64_t x = 0; x < len; x++) {
        EVBDD_WGT v = g[(x ^ t) & D];
        if (v != EVBDD_ZERO && __builtin_parityll(s & x)) v = wgt_neg(v);
        f[x] = v;
    }
    free(g);
}

static void
check_skipping(void)
{
    /* draws of its own, so that the checks of the next family see the
     * vectors they saw before this one was added */
    const uint64_t saved_rng = rng_state;
    const bqd_family_t fam = bqd_family();
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *g = malloc(MAXV * sizeof(EVBDD_WGT));
    unsigned tot = 0, skipping = 0, bad = 0, badeval = 0, badcanon = 0, redundant = 0;
    for (unsigned n = 1; n <= 8; n++) for (int rep = 0; rep < 40; rep++) {
        rand_skipping(n, (rep & 1) ? 0.3 : 0.0, f);
        const BQD e = bqd_from_vector(f, n);
        if (limdd_edge_is_zero(e)) continue;
        tot++;
        if (skips_below(limdd_target(e), 0, n)) skipping++;
        bqd_to_vector(e, n, g);
        if (memcmp(f, g, ((size_t)1 << n) * sizeof(EVBDD_WGT)) != 0) bad++;
        for (int q = 0; q < 8; q++) {
            const uint64_t x = rnd_below(UINT64_C(1) << n);
            if (bqd_eval(e, n, x) != f[x]) { badeval++; break; }
        }
        if (bqd_from_vector(f, n) != e) badcanon++;
        redundant += redundant_below(limdd_target(e), n, g);
    }

    /* the constant one depends on nothing, and a basis state on everything */
    unsigned badplus = 0, badbasis = 0;
    for (unsigned n = 1; n <= 8; n++) {
        const uint64_t len = UINT64_C(1) << n;
        for (uint64_t x = 0; x < len; x++) f[x] = EVBDD_ONE;
        if (limdd_target(bqd_from_vector(f, n)) != LIMDD_TERMINAL) badplus++;
        const uint64_t y = rnd_below(len);
        for (uint64_t x = 0; x < len; x++) f[x] = (x == y) ? EVBDD_ONE : EVBDD_ZERO;
        if (limdd_nodecount(bqd_from_vector(f, n), n) != n) badbasis++;
    }

    /* the note's examples, indexed by (x_top x_bottom), qubit 0 on top */
    unsigned badex = 0;
    EVBDD_WGT c; uint64_t sm, tm;
    if (fam == BQD_FAMILY_PAULI) {
        /* (1, -1): the representative is (1, 1), and Z goes on the root edge */
        f[0] = EVBDD_ONE; f[1] = EVBDD_MIN_ONE;
        const BQD e = bqd_from_vector(f, 1);
        bqd_lim_masks(limdd_label(e), 1, &c, &sm, &tm);
        if (limdd_target(e) != LIMDD_TERMINAL || c != EVBDD_ONE || sm != 1 || tm != 0) badex++;
        /* (1, 0, 1, -1): the copy puts the -1 in the ratio (1, -1), and the
         * high edge is Z_1 to the terminal, skipping qubit 1 */
        f[0] = EVBDD_ONE; f[1] = EVBDD_ZERO; f[2] = EVBDD_ONE; f[3] = EVBDD_MIN_ONE;
        const LIMDD_TARG p = limdd_target(bqd_from_vector(f, 2));
        const LIMDD high = limdd_node_high(p);
        bqd_lim_masks(limdd_label(high), 2, &c, &sm, &tm);
        if (limdd_node_var(p) != 0 || limdd_target(high) != LIMDD_TERMINAL
            || c != EVBDD_ONE || sm != 1 || tm != 0) badex++;
    }
    if (fam == BQD_FAMILY_X) {
        /* (-1)^{x_1 + x_2}: the root's high edge is -1 to the terminal, with
         * no translation at the level it skips */
        f[0] = EVBDD_ONE; f[1] = EVBDD_MIN_ONE; f[2] = EVBDD_MIN_ONE; f[3] = EVBDD_ONE;
        const LIMDD_TARG p = limdd_target(bqd_from_vector(f, 2));
        const LIMDD high = limdd_node_high(p);
        bqd_lim_masks(limdd_label(high), 2, &c, &sm, &tm);
        if (limdd_node_var(p) != 0 || limdd_target(high) != LIMDD_TERMINAL
            || c != EVBDD_MIN_ONE || sm != 0 || tm != 0) badex++;
    }

    char buf[128];
    printf("  %u of %u vectors that skip qubits skip a level in the diagram\n", skipping, tot);
    snprintf(buf, sizeof(buf), "%u of %u vectors do not decode to themselves", bad, tot);
    expect(bad == 0, "round trip, skipping", buf);
    snprintf(buf, sizeof(buf), "%u of %u vectors have an amplitude the query gets wrong", badeval, tot);
    expect(badeval == 0, "amplitude query, skipping", buf);
    snprintf(buf, sizeof(buf), "%u of %u vectors built twice gave two different edges", badcanon, tot);
    expect(badcanon == 0, "determinism, skipping", buf);
    snprintf(buf, sizeof(buf), "%u nodes of %u diagrams do not depend on their variable", redundant, tot);
    expect(redundant == 0, "skip:thm:canon (ii)", buf);
    expect(skipping > 0, "skip:def:fr", "no diagram skipped a level");
    snprintf(buf, sizeof(buf), "%u of 8 constant vectors have a node", badplus);
    expect(badplus == 0, "|+>^n is the terminal", buf);
    snprintf(buf, sizeof(buf), "%u of 8 basis states are not n nodes", badbasis);
    expect(badbasis == 0, "a basis state skips nothing", buf);
    snprintf(buf, sizeof(buf), "%u of the note's examples differ", badex);
    expect(badex == 0, "skip:lem:normal examples", buf);
    free(f); free(g);
    rng_state = saved_rng;
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
    unsigned redundant = 0;
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
                if (!limdd_edge_is_zero(e)) redundant += redundant_below(limdd_target(e), n, g);
                limdd_level_counts(e, counts, n);
                for (unsigned v = 0; v < n; v++) if (counts[v] > max_width) max_width = counts[v];

                if (fam == BQD_FAMILY_X || fam == BQD_FAMILY_PAULI) {
                    /* thm:coset / thm:pcoset, per level: v decided variables,
                     * nonzero nodes only, since that is what this port counts */
                    const unsigned kk = (fam == BQD_FAMILY_X) ? k : k - 1;
                    /* thm:coset counts the node 0 within its + 1; thm:pcoset's
                     * + 2 is the code-state node and 0, so its nonzero nodes
                     * are within + 1 whether or not 0 occurs */
                    const unsigned zmin = limdd_edge_is_zero(e) ? n
                                        : min_zero_var(limdd_target(e), n);
                    for (unsigned v = 0; v < n; v++) {
                        uint64_t bound = 1;
                        for (unsigned i = 0; i < kk; i++) bound += binom(v, i);
                        const size_t zero_node = (fam == BQD_FAMILY_X && v > zmin) ? 1 : 0;
                        if (counts[v] + zero_node > bound) {
                            bound_bad++;
                            fprintf(stderr, "  %s: n=%u k=%u level v=%u has %zu nodes, bound %llu\n",
                                    bqd_family_name(fam), n, k, v, counts[v],
                                    (unsigned long long)bound);
                            break;
                        }
                    }
                    /* the copy never fires on a coset state */
                    if (!limdd_edge_is_zero(e) && copies_below(limdd_target(e), n) != 0) copy_bad++;
                } else {
                    if (!limdd_edge_is_zero(e) && copies_below(limdd_target(e), n) != 0) copies_seen++;
                }

                /* full support: the X-BQD is the BQD, no translation anywhere */
                if (fam == BQD_FAMILY_X && dim == n && !limdd_edge_is_zero(e)
                    && any_translation(limdd_target(e), n)) trans_bad++;

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
                    if (limdd_target(bqd_from_vector(h, n)) != limdd_target(e)) prep_bad++;
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
    snprintf(buf, sizeof(buf), "%u nodes of %u coset states do not depend on their variable", redundant, states);
    expect(redundant == 0, "skip:thm:canon (ii), coset states", buf);
    free(f); free(g); free(code); free(terms);
}

/** Whether arg(w) is in [0, pi), from the value as doubles. */
static bool
upper_by_double(EVBDD_WGT w)
{
    const complex_t z = weight_as_complex(w);
    return z.i > 0.0 || (z.i == 0.0 && z.r > 0.0);
}

/** Whether deciding the sign of A + B sqrt2 needs the comparison of A^2 with 2B^2. */
static bool
mixed(const mpq_t A, const mpq_t B)
{
    return (mpq_sgn(A) > 0 && mpq_sgn(B) < 0) || (mpq_sgn(A) < 0 && mpq_sgn(B) > 0);
}

/**
 * def:prep on values the eighth roots of unity never produce. Every entry is
 * a + b sqrt2 + i (c + d sqrt2) with small integers of either sign, so the
 * ratios at the pivots have components like 1 - sqrt2, whose sign is the
 * comparison sign_sqrt2 makes, and every one of these vectors has full
 * support, so the pivots are the unit vectors and span it. Checked: the
 * round trip, the argument at every pivot of the representative (against the
 * doubles, which cannot be wrong here: an exact zero is 0.0, and A + B sqrt2
 * with integers A, B not both 0 is at least 1/(|A| + |B| sqrt2) from it, far
 * above the rounding at these sizes), and invariance under c Z^s.
 */
static void
check_prep_algebraic(void)
{
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *g = malloc(MAXV * sizeof(EVBDD_WGT)),
              *r = malloc(MAXV * sizeof(EVBDD_WGT));
    unsigned tot = 0, badrt = 0, badarg = 0, badinv = 0;
    uint64_t pivots = 0, needed_compare = 0;
    for (unsigned n = 1; n <= 7; n++) for (int rep = 0; rep < 40; rep++) {
        const uint64_t len = UINT64_C(1) << n;
        for (uint64_t x = 0; x < len; x++) {
            long a, b, c, d;
            do {
                a = (long)rnd_below(7) - 3; b = (long)rnd_below(7) - 3;
                c = (long)rnd_below(7) - 3; d = (long)rnd_below(7) - 3;
            } while (a == 0 && b == 0 && c == 0 && d == 0);
            f[x] = qisq2_lookup(a, 1, b, 1, c, 1, d, 1);
        }
        tot++;
        const BQD e = bqd_from_vector(f, n);
        bqd_to_vector(e, n, g);
        if (memcmp(f, g, len * sizeof(EVBDD_WGT)) != 0) badrt++;

        EVBDD_WGT c0; uint64_t sv, p;
        bqd_representative(f, len, r, &c0, &sv, &p);
        for (unsigned b = 0; b < n; b++) {
            const EVBDD_WGT w = r[UINT64_C(1) << b];
            pivots++;
            if (!upper_by_double(w)) badarg++;
            const qisq2_t *q = (const qisq2_t *)qisq2_map_get(wgt_storage, (uint64_t)w);
            const bool im_zero = mpq_sgn(q->c) == 0 && mpq_sgn(q->d) == 0;
            if (mixed(q->c, q->d) || (im_zero && mixed(q->a, q->b))) needed_compare++;
        }

        const EVBDD_WGT cc = f[rnd_below(len)];
        const uint64_t s = rnd_below(len);
        for (uint64_t y = 0; y < len; y++) {
            EVBDD_WGT v = wgt_mul(cc, f[y]);
            g[y] = __builtin_parityll(s & y) ? wgt_neg(v) : v;
        }
        if (limdd_target(bqd_from_vector(g, n)) != limdd_target(e)) badinv++;
    }
    char buf[160];
    printf("  %u vectors over Z[sqrt2, i]: %llu pivot values, %llu of them decided by "
           "comparing A^2 with 2B^2\n", tot, (unsigned long long)pivots,
           (unsigned long long)needed_compare);
    snprintf(buf, sizeof(buf), "%u of %u vectors do not decode to themselves", badrt, tot);
    expect(badrt == 0, "round trip over Z[sqrt2, i]", buf);
    snprintf(buf, sizeof(buf), "%u of %llu pivot values have an argument outside [0, pi)",
             badarg, (unsigned long long)pivots);
    expect(badarg == 0, "def:prep pivot arguments", buf);
    snprintf(buf, sizeof(buf), "%u of %u: (c Z^s) f and f got different nodes", badinv, tot);
    expect(badinv == 0, "def:prep invariance over Z[sqrt2, i]", buf);
    expect(needed_compare > 0, "sign_sqrt2 exercised",
           "no pivot value needed the A^2 against 2B^2 comparison");
    free(f); free(g); free(r);
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
    if (fam == BQD_FAMILY_PAULI) check_prep_algebraic();
    check_skipping();
    check_countnodes();

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

/**
 * The Pauli family on float weights, where arg_in_upper compares doubles.
 * Amplitudes are eighth roots of unity and zero, which the weight table
 * holds as single entries, so -1 is stored with an imaginary part of exactly
 * 0.0 and the comparison sees what the exact one does. Values are compared
 * to 1e-9, not by index.
 */
TASK_0(int, run_float_pauli)
{
    bqd_init(BQD_FAMILY_PAULI, NQ, 1LL << 20, 1LL << 20, 1LL << 22, 1LL << 20);
    for (int e = 0; e < 8; e++) pw[e] = complex_lookup_angle((fl_t)e * 0.25 * M_PI, 1.0);
    pw[0] = EVBDD_ONE;
    pw[4] = EVBDD_MIN_ONE;

    const int before = failures;
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *g = malloc(MAXV * sizeof(EVBDD_WGT)),
              *r = malloc(MAXV * sizeof(EVBDD_WGT));
    unsigned tot = 0, badrt = 0, badarg = 0, badinv = 0, full = 0;
    for (unsigned n = 1; n <= 8; n++) for (int rep = 0; rep < 40; rep++) {
        const uint64_t len = UINT64_C(1) << n;
        rand_vector(n, (rep & 1) ? 0.0 : 0.3, f);
        tot++;
        const BQD e = bqd_from_vector(f, n);
        bqd_to_vector(e, n, g);
        for (uint64_t x = 0; x < len; x++) {
            const complex_t a = weight_as_complex(f[x]), b = weight_as_complex(g[x]);
            if (fabs(a.r - b.r) > 1e-9 || fabs(a.i - b.i) > 1e-9) { badrt++; break; }
        }
        if (lexmin_nonzero(f, len) != 0) continue;
        bool support_full = true;
        for (uint64_t x = 0; x < len; x++) if (f[x] == EVBDD_ZERO) support_full = false;
        if (!support_full) continue;
        full++;
        EVBDD_WGT c0; uint64_t sv, p;
        bqd_representative(f, len, r, &c0, &sv, &p);
        for (unsigned b = 0; b < n; b++) if (!upper_by_double(r[UINT64_C(1) << b])) { badarg++; break; }
        const EVBDD_WGT cc = pw[rnd_below(8)];
        const uint64_t s = rnd_below(len);
        for (uint64_t y = 0; y < len; y++) {
            const EVBDD_WGT v = wgt_mul(cc, f[y]);
            g[y] = __builtin_parityll(s & y) ? wgt_neg(v) : v;
        }
        if (limdd_target(bqd_from_vector(g, n)) != limdd_target(e)) badinv++;
    }
    char buf[128];
    printf("  %u vectors, %u of full support\n", tot, full);
    snprintf(buf, sizeof(buf), "%u of %u vectors do not decode to themselves", badrt, tot);
    expect(badrt == 0, "float round trip", buf);
    snprintf(buf, sizeof(buf), "%u of %u representatives have a pivot argument outside [0, pi)", badarg, full);
    expect(badarg == 0, "float def:prep pivot arguments", buf);
    snprintf(buf, sizeof(buf), "%u of %u: (c Z^s) f and f got different nodes", badinv, full);
    expect(badinv == 0, "float def:prep invariance", buf);
    free(f); free(g); free(r);
    bqd_quit();
    return failures != before;
}

int
main(void)
{
    int bad = 0;
    bad |= run_family(BQD_FAMILY_SCALAR);
    bad |= run_family(BQD_FAMILY_X);
    bad |= run_family(BQD_FAMILY_PAULI);

    printf("== Pauli-BQD, float weights ==\n");
    lace_start(4, 0);
    sylvan_set_sizes(1LL << 20, 1LL << 20, 1LL << 20, 1LL << 20);
    sylvan_init_package();
    qsylvan_init_simulator(1LL << 22, 1LL << 22, -1, COMP_HASHMAP, NORM_LOW);
    const int res = RUN(run_float_pauli);
    sylvan_quit();
    lace_stop();
    printf("  %s\n", res ? "FAILED" : "ok");
    bad |= res;
    return bad;
}
