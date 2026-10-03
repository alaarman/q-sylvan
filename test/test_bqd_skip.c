/*
 * Level skipping in the BQD, checked so that a wrong skip rule cannot pass:
 * a diagram with a level too many or too few fails here, whatever made it.
 * The rule is that of the note "Binary Quotient Diagrams with Level
 * Skipping" (bqd-skip.tex), whose labels are cited as skip:...
 *
 *   structure     every diagram below, node by node: no node whose function
 *                 does not depend on its variable (skip:thm:canon (ii)); a
 *                 low edge that is the identity or zero (skip:lem:normal
 *                 (c)); no X at a level an edge skips (skip:lem:normal (a))
 *                 and no label bit above the level an edge is read at; on a
 *                 Pauli root edge a Z at a skipped level exactly where the
 *                 vector changes sign across it (skip:lem:normal (b)); and
 *                 every node the one bqd_from_vector builds for its own
 *                 function; no label with a Z and an X at one qubit. In
 *                 every family also Cof and Compose (skip:alg:apply,
 *                 skip:alg:constructors, skip:alg:xcompose): at every node
 *                 the two cofactors are canonical and compose back to the
 *                 node, and at every level an edge skips the edge is both of
 *                 its cofactors and composes with itself to itself, or in
 *                 the Pauli family, where its label has a Z there, the two
 *                 are opposite and compose back to the edge
 *   build         decode(build(f)) = f, by bqd_to_vector and by a decoder
 *                 written here from skip:def:skip, and bqd_eval agrees at
 *                 every index: on every vector over {0, 1, -1, i} on 1 to 3
 *                 qubits, and at random on 4 to 8 qubits on Clifford+T
 *                 circuit states, phase states, product states and states
 *                 with zeros that ignore random qubits, in all three
 *                 families; |+>^n is the terminal under one label, and in
 *                 the Pauli family so is a character; a character has one
 *                 node per qubit it reads in the other two, and a basis
 *                 state n nodes in all three, bqd_basis_state's, which in the
 *                 other two is X^x on the chain of |0...0>
 *   operations    every family: every operation of qsylvan_bqd_gates.h and
 *                 qsylvan_bqd_ops.h returns the edge bqd_from_vector builds
 *                 for the dense result, and a diagram that passes the
 *                 structure checks. OP_XQUOT is not public and is checked
 *                 inside bqd_compose, and OP_XPROD inside the high side of
 *                 bqd_cofactor. The states ignore random qubits, so targets,
 *                 controls and restricted qubits fall on skipped levels as
 *                 often as not, and each such case is counted and must occur;
 *                 in the translation and Pauli families they are translated
 *                 and signed at random, and Gate, Restrict and Project also
 *                 get random labelled edges on their nodes, as their
 *                 recursions hand themselves
 *   operations    translation and Pauli families: the pointwise product,
 *   (X, Pauli)    sum, xprod and xquot, scale, negate, the cofactor at every
 *                 level down to an edge's top, Compose, and Canon of labelled
 *                 edges (skip:sec:xp), each result the edge bqd_from_vector
 *                 builds for the dense result and a diagram that passes the
 *                 structure checks: exhaustively on every vector over
 *                 {0, 1, -1, i} on 1 to 3 qubits for the one-operand ones and
 *                 Compose, on every pair on 1 and 2 qubits and a random
 *                 partner on 3 for the pointwise ones, and at random on 4 to
 *                 8 qubits; the operands canonical, and random labelled edges
 *                 on their nodes, as the recursion hands Apply and Canon;
 *                 every case of skip:prop:xcompose counted and met
 *   gates         translation and Pauli families: every gate of the table
 *   (X, Pauli)    and the dynamic one on every qubit, every set of controls
 *                 above every target, every phase with controls anywhere,
 *                 cgate_either with X and Z and swap and the local matvec on
 *                 every pair, restriction and projection on every qubit, the
 *                 monomial on every set with every phase, the product with
 *                 every partner and the full-support test, on every vector
 *                 over {0, 1, -1, i} on 1 and 2 qubits, and a random draw of
 *                 them on one vector in eight on 3; skipped targets and
 *                 controls, a root that skips, a zero low cofactor and a
 *                 cancellation each counted and met
 *   skip:prop:diag  the diagonal walk visits at most n - ctz(A) nodes, the
 *                 path the proposition names, counted here from the diagram,
 *                 and stops where the state skips a qubit of A; in every
 *                 family, the Pauli family's walk (skip:alg:xdiag) with its
 *                 sign repair and its skip at the last qubit of A met, and a
 *                 monomial there may have fewer nodes than qubits
 *   collections   every family, limdd_gc between operations in small
 *                 tables, so that swept buckets are soon built on again:
 *                 gate sequences on 2 to 8 qubits in a table two wider, the
 *                 state the edge bqd_from_vector builds after every gate; and
 *                 protected results on 8 and 6 qubits under the operations
 *                 above, each result checked as there, a 6-qubit one also
 *                 against the same operation on its edge read on 8 qubits,
 *                 which no memo key may tell apart, and after every
 *                 collection every kept edge decoded and built again, and
 *                 every node met after it checked from scratch; and in the
 *                 translation and Pauli families, from any basis state and on
 *                 translated and signed states, and the same kept results
 *                 under the operations of skip:sec:xp, on labelled operands
 *                 too
 *   sessions      two BQD sessions in one Sylvan package, of two families
 *                 or of one, the same gates in each and every state checked:
 *                 no memo entry outlives bqd_quit; every worker asking for a
 *                 collection at once, a different number of times each,
 *                 under an alarm; and a Pauli table smaller than the LIM
 *                 table, which only the trigger keeps from filling
 *
 * The scalar family runs the copy rule of def:bqd here, which bqd_init_rule
 * gives on request; its default, rule SM, is test_bqd_sm's. Exact weights, in
 * Q(w_8, sqrt2), throughout. Every check runs on the number
 * of workers in BQD_SKIP_WORKERS (default 4), BQD_SKIP_REPEAT times (default
 * 1), each run in a fresh session on the same draws; BQD_SKIP_FAMILY (0, 1
 * or 2) runs one family alone, with its collections, and 3 the sessions
 * alone.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <sylvan.h>
#include <sylvan_int.h>

#include "qsylvan_bqd.h"
#include "qsylvan_bqd_ops.h"
#include "qsylvan_bqd_gates.h"
#include "qsylvan_bqd_xp.h"
#include "qsylvan_limdd_gc.h"
#include "qsylvan_limdd_inspect.h"
#include "qsylvan_simulator.h"
#include "sylvan_edge_weights_qisq2.h"

#define NQ 8
#define MAXV (1u << NQ)

/*
 * Under a sanitizer every check runs some twenty times slower, and this test
 * came near ctest's timeout on a loaded machine. There the draws on three
 * qubits, the gate sequences and the collection pools are thinned by THIN;
 * every kind of check still runs, and every case it counts is still met.
 */
#if defined(__SANITIZE_ADDRESS__) || defined(__SANITIZE_THREAD__)
#define THIN 4
#elif defined(__has_feature)
#if __has_feature(address_sanitizer) || __has_feature(thread_sanitizer) || \
    __has_feature(memory_sanitizer)
#define THIN 4
#endif
#endif
#ifndef THIN
#define THIN 1
#endif

/*
 * The diagram here is the copy rule's, def:bqd, in every family: the
 * structure checks below are written from it. The scalar family's default is
 * rule SM, which test_bqd_sm checks, so every session here asks for the copy
 * rule by name.
 */
static void
init_copy(bqd_family_t fam, size_t n, size_t nodes, size_t paulis, size_t lims, size_t stabs)
{
    bqd_init_rule(fam, BQD_ZERO_COPY, n, nodes, paulis, lims, stabs);
}

static int failures = 0;
static EVBDD_WGT pw[8];                     /* w_8^e */
static EVBDD_WGT isq2;                      /* 1/sqrt2 */

/** w_8^e and 1/sqrt2, looked up again in every session's weight table. */
static void
exact_constants(void)
{
    pw[0] = EVBDD_ONE;
    const EVBDD_WGT w = qisq2_lookup(0, 1, 1, 2, 0, 1, 1, 2);
    for (int e = 1; e < 8; e++) pw[e] = wgt_mul(pw[e - 1], w);
    isq2 = qisq2_lookup(0, 1, 1, 2, 0, 1, 0, 1);
}

static uint64_t rng_state;
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
    printf("  %-56s %s (%u)\n", what, bad ? "FAILED" : "ok", tot);
}

/** A case the checks must have met, or they say nothing about it. */
static void
covered(const char *what, unsigned seen)
{
    expect(seen > 0, what, "never met");
    printf("  %-56s %u\n", what, seen);
}

static inline bool parity(uint64_t v) { return __builtin_parityll(v) != 0; }

/** Qubit q as a vector-index bit: qubit 0 is the most significant of n. */
static inline uint64_t
ibit(unsigned q, unsigned n)
{
    return UINT64_C(1) << (n - 1 - q);
}

/** A set of qubits, bit q for qubit q as in a control mask, as an index mask. */
static uint64_t
imask(uint64_t Q, unsigned n)
{
    uint64_t m = 0;
    for (unsigned q = 0; q < n; q++) if ((Q >> q) & 1) m |= ibit(q, n);
    return m;
}

static inline BQD
unit(LIMDD_TARG t)
{
    return limdd_bundle(LIMDD_LIM_IDENTITY, t);
}

/** The level of a node: its variable, or n for the terminal of an n-qubit vector. */
static inline unsigned
level_of(LIMDD_TARG p, unsigned n)
{
    return p == LIMDD_TERMINAL ? n : limdd_node_var(p);
}

/**
 * The top of an edge read on n qubits (skip:sec:xp:edges): the level of its
 * node, or in the Pauli family a Z above it, where the function starts to
 * depend on its variables; n for zero and for a constant.
 */
static inline unsigned
top_of(BQD e, unsigned n)
{
    if (limdd_edge_is_zero(e)) return n;
    const unsigned t = bqd_xp_top(e);
    return t < n ? t : n;
}

static void
zero_fill(EVBDD_WGT *v, uint64_t len)
{
    for (uint64_t y = 0; y < len; y++) v[y] = EVBDD_ZERO;
}

static bool
same_vector(const EVBDD_WGT *a, const EVBDD_WGT *b, unsigned n)
{
    return memcmp(a, b, ((size_t)1 << n) * sizeof(EVBDD_WGT)) == 0;
}

/*
 * A set of nodes, for walks that meet a node more than once. A slot is in
 * the set while its stamp is the current one, so clearing is one increment.
 * It is only a memo: when it fills up it is cleared, and a node is then
 * walked again.
 */
#define SET_BITS 18
typedef struct {
    LIMDD_TARG key[1u << SET_BITS];
    uint32_t   stamp[1u << SET_BITS];
    uint32_t   now, count;
} nodeset_t;

static nodeset_t *walked, *checked;

static void
set_clear(nodeset_t *s)
{
    s->now++;
    s->count = 0;
}

/** Add p, and say whether it was new. */
static bool
set_add(nodeset_t *s, LIMDD_TARG p)
{
    if (s->count >= (1u << SET_BITS) / 2) set_clear(s);
    uint64_t i = ((uint64_t)p * UINT64_C(0x9E3779B97F4A7C15)) >> (64 - SET_BITS);
    for (;;) {
        if (s->stamp[i] != s->now) { s->stamp[i] = s->now; s->key[i] = p; s->count++; return true; }
        if (s->key[i] == p) return false;
        i = (i + 1) & ((1u << SET_BITS) - 1);
    }
}

/* --- the vector of a diagram, from skip:def:skip -------------------------- */

/** The label c Z^s X^t applied to u, the 2^(n-var) values of an edge read at var. */
static void
apply_label(LIMDD_LIM lim, unsigned n, const EVBDD_WGT *u, uint64_t len, EVBDD_WGT *out)
{
    EVBDD_WGT c; uint64_t s, t;
    bqd_lim_masks(lim, n, &c, &s, &t);
    for (uint64_t y = 0; y < len; y++) {
        EVBDD_WGT v = u[(y ^ t) & (len - 1)];
        if (v != EVBDD_ZERO) {
            if (parity(s & y)) v = wgt_neg(v);
            v = wgt_mul(c, v);
        }
        out[y] = v;
    }
}

static void node_dense(LIMDD_TARG p, unsigned n, EVBDD_WGT *out);

/** ext of p's function to the levels var..n-1: its block, repeated. */
static void
ext_dense(LIMDD_TARG p, unsigned var, unsigned n, EVBDD_WGT *out)
{
    const uint64_t blk = UINT64_C(1) << (n - level_of(p, n));
    node_dense(p, n, out);
    for (uint64_t y = blk; y < (UINT64_C(1) << (n - var)); y++) out[y] = out[y - blk];
}

/**
 * The function of node p read at its own level v, 2^(n-v) values, by the
 * quotient rule with skipping edges: f_0 = l_0 . ext[N_0] and
 * f_1 = l_1 . (ext[N_0] (.) ext[N_1]), the product copying the high value
 * where the low one is zero. Written from the definition, not from
 * bqd_decode, so that the two check each other.
 */
static void
node_dense(LIMDD_TARG p, unsigned n, EVBDD_WGT *out)
{
    if (p == LIMDD_TERMINAL) { out[0] = EVBDD_ONE; return; }
    const unsigned v = limdd_node_var(p);
    const uint64_t h = UINT64_C(1) << (n - v - 1);
    const LIMDD low = limdd_node_low(p), high = limdd_node_high(p);
    EVBDD_WGT *u0 = malloc(h * sizeof(EVBDD_WGT)), *u1 = malloc(h * sizeof(EVBDD_WGT));
    if (limdd_edge_is_zero(low)) zero_fill(u0, h); else ext_dense(limdd_target(low), v + 1, n, u0);
    if (limdd_edge_is_zero(high)) zero_fill(u1, h); else ext_dense(limdd_target(high), v + 1, n, u1);
    apply_label(limdd_label(low), n, u0, h, out);
    for (uint64_t y = 0; y < h; y++)
        if (u0[y] != EVBDD_ZERO) u1[y] = (u1[y] == EVBDD_ZERO) ? EVBDD_ZERO : wgt_mul(u0[y], u1[y]);
    apply_label(limdd_label(high), n, u1, h, out + h);
    free(u0); free(u1);
}

/** The function the edge e denotes read at level var, 2^(n-var) values. */
static void
edge_dense(BQD e, unsigned var, unsigned n, EVBDD_WGT *out)
{
    const uint64_t len = UINT64_C(1) << (n - var);
    if (limdd_edge_is_zero(e)) { zero_fill(out, len); return; }
    EVBDD_WGT *u = malloc(len * sizeof(EVBDD_WGT));
    ext_dense(limdd_target(e), var, n, u);
    apply_label(limdd_label(e), n, u, len, out);
    free(u);
}

/* --- the structure of a diagram ------------------------------------------- */

typedef struct {
    unsigned nodes;         /* distinct nodes checked */
    unsigned redundant;     /* a node whose two cofactors are equal */
    unsigned low_label;     /* a low edge neither the identity nor zero */
    unsigned x_skipped;     /* an X at a level the edge skips */
    unsigned scope;         /* a label bit above the level the edge is read at */
    unsigned sign;          /* a Z at a skipped level of a root the vector does not show */
    unsigned not_canon;     /* a node that is not bqd_from_vector's for its function */
    unsigned cof;           /* Cof or Compose disagrees at a node or a skipped level */
    unsigned zx;            /* a Z and an X at one qubit of a label */
    unsigned z_skipped;     /* met: a Z at a level the edge skips */
    unsigned translated;    /* met: an edge with a translation */
    unsigned skipping;      /* met: an edge that skips a level */
} structure_t;

static structure_t sf;

/** The checks on an edge read at level r, which are about its label and the levels it skips. */
static void
check_edge(BQD e, unsigned r, unsigned n)
{
    if (limdd_edge_is_zero(e)) return;
    const LIMDD_TARG p = limdd_target(e);
    const unsigned v = level_of(p, n);
    EVBDD_WGT c; uint64_t s, t;
    bqd_lim_masks(limdd_label(e), n, &c, &s, &t);
    /* the variables r..n-1 are the index bits below n-r, and those the edge
     * skips, r..v-1, are the bits from n-v up */
    const uint64_t scope = (UINT64_C(1) << (n - r)) - 1;
    const uint64_t skipped = scope & ~((UINT64_C(1) << (n - v)) - 1);
    if ((s | t) & ~scope) sf.scope++;
    if (t & skipped) sf.x_skipped++;
    /* the translation is a least point of a support, and a 1 of it is at
     * neither a pivot nor a level the support ignores, so no Z shares it */
    if (s & t) sf.zx++;
    if (s & skipped) sf.z_skipped++;
    if (t != 0) sf.translated++;
    if (v <= r) return;
    sf.skipping++;
    /* an edge that skips m denotes a function that does not depend on x_m:
     * it is both of its cofactors there, and Compose of it with itself is
     * itself (skip:lem:virtual, Compose's first line). In the Pauli family a
     * Z at m makes the two cofactors opposite, and Compose of those puts the
     * Z back (skip:lem:normal (b), skip:prop:xcompose (ii)). The levels are
     * taken from the top down, each read on the low cofactor of the one
     * above, which has shed the Z bits above it. */
    BQD cur = e;
    for (unsigned m = r; m < v; m++) {
        bqd_lim_masks(limdd_label(cur), n, &c, &s, &t);
        const BQD c0 = bqd_cofactor(cur, m, 0), c1 = bqd_cofactor(cur, m, 1);
        const bool ok = (s & ibit(m, n)) ? c0 != cur && c1 == bqd_negate(c0)
                                         : c0 == cur && c1 == cur;
        if (!ok || bqd_compose(m, c0, c1) != cur) {
            sf.cof++;
            break;
        }
        cur = c0;
    }
}

/** The checks on node p and every node below it, each node once per section. */
static void
check_node(LIMDD_TARG p, unsigned n)
{
    if (p == LIMDD_TERMINAL || !set_add(checked, p)) return;
    sf.nodes++;
    const unsigned v = limdd_node_var(p);
    const LIMDD low = limdd_node_low(p), high = limdd_node_high(p);
    if (!limdd_edge_is_zero(low) && limdd_label(low) != LIMDD_LIM_IDENTITY) sf.low_label++;
    check_edge(low, v + 1, n);
    check_edge(high, v + 1, n);

    /* its function read at v, the first 2^(n-v) values, extended to all n qubits */
    const uint64_t len = UINT64_C(1) << n, h = UINT64_C(1) << (n - v - 1);
    EVBDD_WGT *g = malloc(len * sizeof(EVBDD_WGT));
    ext_dense(p, 0, n, g);
    if (memcmp(g, g + h, h * sizeof(EVBDD_WGT)) == 0) sf.redundant++;
    if (bqd_from_vector(g, n) != unit(p)) sf.not_canon++;

    /* its cofactors are the canonical edges of the two halves, extended, and
     * composing them gives the node back: a node stores a representative, so
     * its own edge carries the identity in every family */
    EVBDD_WGT *g0 = malloc(len * sizeof(EVBDD_WGT)), *g1 = malloc(len * sizeof(EVBDD_WGT));
    for (uint64_t x = 0; x < len; x++) { g0[x] = g[x & ~h]; g1[x] = g[x | h]; }
    const BQD c0 = bqd_cofactor(unit(p), v, 0), c1 = bqd_cofactor(unit(p), v, 1);
    if (c0 != bqd_from_vector(g0, n) || c1 != bqd_from_vector(g1, n)
        || bqd_compose(v, c0, c1) != unit(p)) sf.cof++;
    free(g0); free(g1);
    free(g);

    if (!limdd_edge_is_zero(low)) check_node(limdd_target(low), n);
    if (!limdd_edge_is_zero(high)) check_node(limdd_target(high), n);
}

/** The structure of e, the diagram of f on n qubits. */
static void
check_diagram(BQD e, const EVBDD_WGT *f, unsigned n)
{
    if (limdd_edge_is_zero(e)) return;
    check_edge(e, 0, n);
    const unsigned v = level_of(limdd_target(e), n);
    if (bqd_family() == BQD_FAMILY_PAULI && v > 0) {
        /* skip:lem:normal (b): across a level the root skips the vector is
         * the same, or, where the label has a Z, the same negated */
        EVBDD_WGT c; uint64_t s, t;
        bqd_lim_masks(limdd_label(e), n, &c, &s, &t);
        for (unsigned m = 0; m < v; m++) {
            const uint64_t b = ibit(m, n);
            bool ok = true;
            for (uint64_t x = 0; x < (UINT64_C(1) << n) && ok; x++)
                if (!(x & b)) ok = f[x | b] == ((s & b) ? wgt_neg(f[x]) : f[x]);
            if (!ok) { sf.sign++; break; }
        }
    }
    check_node(limdd_target(e), n);
}

static void
structure_begin(void)
{
    memset(&sf, 0, sizeof(sf));
    set_clear(checked);
}

static void
structure_report(const char *section)
{
    const struct { const char *what; unsigned bad; } k[] = {
        { "a node whose two cofactors are equal (skip:thm:canon (ii))", sf.redundant },
        { "a low edge neither the identity nor zero (skip:lem:normal (c))", sf.low_label },
        { "an X at a level the edge skips (skip:lem:normal (a))", sf.x_skipped },
        { "a label bit above the level its edge is read at", sf.scope },
        { "a Z on a root at a skipped level the vector does not show (skip:lem:normal (b))", sf.sign },
        { "a node that is not the one bqd_from_vector builds for its function", sf.not_canon },
        { "Cof or Compose disagrees at a node or a skipped level", sf.cof },
        { "a label with a Z and an X at one qubit", sf.zx },
    };
    unsigned bad = 0;
    char buf[160];
    for (size_t i = 0; i < sizeof(k) / sizeof(k[0]); i++) {
        if (k[i].bad == 0) continue;
        snprintf(buf, sizeof(buf), "%u times %s", k[i].bad, k[i].what);
        expect(false, section, buf);
        bad += k[i].bad;
    }
    snprintf(buf, sizeof(buf), "%s: structure", section);
    printf("  %-56s %s (%u nodes, %u skipping edges)\n", buf, bad ? "FAILED" : "ok",
           sf.nodes, sf.skipping);
}

/**
 * The qubits some edge of e skips, bit q for qubit q, and whether an edge
 * reaches the terminal from above the last level.
 */
static uint64_t
skipped_below(LIMDD_TARG p, unsigned r, unsigned n, bool *to_terminal)
{
    const unsigned v = level_of(p, n);
    uint64_t m = 0;
    for (unsigned q = r; q < v; q++) m |= UINT64_C(1) << q;
    if (p == LIMDD_TERMINAL) { if (v > r) *to_terminal = true; return m; }
    if (!set_add(walked, p)) return m;
    const LIMDD e[2] = { limdd_node_low(p), limdd_node_high(p) };
    for (int b = 0; b < 2; b++)
        if (!limdd_edge_is_zero(e[b])) m |= skipped_below(limdd_target(e[b]), v + 1, n, to_terminal);
    return m;
}

static uint64_t
skipped_qubits(BQD e, unsigned n, bool *to_terminal)
{
    *to_terminal = false;
    if (limdd_edge_is_zero(e)) return 0;
    set_clear(walked);
    return skipped_below(limdd_target(e), 0, n, to_terminal);
}

/** The qubits the vector depends on, bit q for qubit q. */
static uint64_t
depends_on(const EVBDD_WGT *f, unsigned n)
{
    uint64_t D = 0;
    for (unsigned q = 0; q < n; q++) {
        const uint64_t b = ibit(q, n);
        for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
            if (!(x & b) && f[x] != f[x | b]) { D |= UINT64_C(1) << q; break; }
    }
    return D;
}

/* --- vectors: qubit q is bit n-1-q of an index ----------------------------- */

static EVBDD_WGT
rand_alg(void)
{
    long a, b, c, d;
    do {
        a = (long)rnd_below(5) - 2; b = (long)rnd_below(5) - 2;
        c = (long)rnd_below(5) - 2; d = (long)rnd_below(5) - 2;
    } while (a == 0 && b == 0 && c == 0 && d == 0);
    return qisq2_lookup(a, 1, b, 2, c, 1, d, 2);
}

static EVBDD_WGT
rand_value(double pzero, bool alg)
{
    if ((double)rnd_below(1000) / 1000.0 < pzero) return EVBDD_ZERO;
    if (alg && rnd_below(2)) return rand_alg();
    return pw[rnd_below(8)];
}

/** An entry of a random exact matrix: zero a quarter of the time. */
static EVBDD_WGT
rand_entry(void)
{
    switch (rnd_below(8)) {
    case 0: case 1: return EVBDD_ZERO;
    case 2: case 3: return wgt_mul(isq2, pw[rnd_below(8)]);
    case 4: case 5: return pw[rnd_below(8)];
    default:        return rand_alg();
    }
}

/** A random set of qubits, each in it with probability 1/2. */
static uint64_t
rand_qubits(unsigned n)
{
    return rnd_below(UINT64_C(1) << n);
}

/** A random function of the qubits in D alone, repeated over the others. */
static void
rand_on(unsigned n, uint64_t D, double pzero, bool alg, EVBDD_WGT *f)
{
    const uint64_t dm = imask(D, n);
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
        f[x] = (x & ~dm) ? f[x & dm] : rand_value(pzero, alg);
}

/** c . w_8^{P(x)}, P a random polynomial over Z_8 of degree at most 3 in the qubits of D. */
static void
phase_on(unsigned n, uint64_t D, EVBDD_WGT *f)
{
    const uint64_t dm = imask(D, n);
    uint64_t W[96]; unsigned c[96], nt = 0;
    for (uint64_t w = 1; w < (UINT64_C(1) << n); w++)
        if ((w & ~dm) == 0 && __builtin_popcountll(w) <= 3 && rnd_below(2)) {
            W[nt] = w; c[nt] = (unsigned)rnd_below(8); nt++;
        }
    const EVBDD_WGT scale = pw[rnd_below(8)];
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) {
        unsigned e = 0;
        for (unsigned i = 0; i < nt; i++) if ((x & W[i]) == W[i]) e += c[i];
        f[x] = wgt_mul(scale, pw[e & 7]);
    }
}

/** U, row-major, on qubit q where the qubits of cmask (bit c for qubit c) are 1. */
static void
dense_gate(EVBDD_WGT *v, unsigned n, const EVBDD_WGT *U, uint64_t cmask, unsigned q)
{
    const uint64_t qb = ibit(q, n), cm = imask(cmask, n);
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) {
        if ((x & qb) || (x & cm) != cm) continue;
        const EVBDD_WGT a = v[x], b = v[x | qb];
        v[x]      = wgt_add(wgt_mul(U[0], a), wgt_mul(U[1], b));
        v[x | qb] = wgt_add(wgt_mul(U[2], a), wgt_mul(U[3], b));
    }
}

static void
dense_swap(EVBDD_WGT *v, unsigned n, unsigned a, unsigned b)
{
    const uint64_t ab = ibit(a, n), bb = ibit(b, n);
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) {
        if ((x & ab) && !(x & bb)) {
            const uint64_t y = (x & ~ab) | bb;
            const EVBDD_WGT t = v[x]; v[x] = v[y]; v[y] = t;
        }
    }
}

/** A random Clifford+T circuit from |0...0>, the gate set of the qasm runner. */
static void
circuit_state(unsigned n, EVBDD_WGT *f)
{
    static const uint32_t one[] = { GATEID_H, GATEID_H, GATEID_S, GATEID_T, GATEID_X,
                                    GATEID_Z, GATEID_Tdag };
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) f[x] = x ? EVBDD_ZERO : EVBDD_ONE;
    const unsigned len = n + (unsigned)rnd_below(6 * n);
    for (unsigned i = 0; i < len; i++) {
        const unsigned a = (unsigned)rnd_below(n);
        if (n >= 2 && rnd_below(3) == 0) {
            unsigned b = (unsigned)rnd_below(n - 1); if (b >= a) b++;
            dense_gate(f, n, gates[rnd_below(2) ? GATEID_X : GATEID_Z], UINT64_C(1) << a, b);
        } else {
            dense_gate(f, n, gates[one[rnd_below(sizeof(one) / sizeof(one[0]))]], 0, a);
        }
    }
}

/** A tensor product of random one-qubit states, |+> and |-> among them. */
static void
product_state(unsigned n, EVBDD_WGT *f)
{
    EVBDD_WGT a[NQ][2];
    for (unsigned q = 0; q < n; q++) {
        switch (rnd_below(6)) {
        case 0:  a[q][0] = EVBDD_ONE;  a[q][1] = EVBDD_ZERO;    break;
        case 1:  a[q][0] = EVBDD_ZERO; a[q][1] = EVBDD_ONE;     break;
        case 2:  a[q][0] = EVBDD_ONE;  a[q][1] = EVBDD_ONE;     break;
        case 3:  a[q][0] = EVBDD_ONE;  a[q][1] = EVBDD_MIN_ONE; break;
        case 4:  a[q][0] = EVBDD_ONE;  a[q][1] = pw[2];         break;
        default: a[q][0] = rand_value(0.0, true); a[q][1] = rand_value(0.0, true); break;
        }
    }
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) {
        EVBDD_WGT v = EVBDD_ONE;
        for (unsigned q = 0; q < n && v != EVBDD_ZERO; q++) {
            const EVBDD_WGT w = a[q][(x >> (n - 1 - q)) & 1];
            v = (w == EVBDD_ZERO) ? EVBDD_ZERO : wgt_mul(v, w);
        }
        f[x] = v;
    }
}

enum { S_ZEROS, S_ALGEBRAIC, S_PHASE, S_CIRCUIT, S_PRODUCT, S_BASIS_PLUS, S_KINDS };

/**
 * A state of one of six kinds: a function of random qubits with zeros, the
 * same over Z[sqrt2, i]/2 without zeros, a phase state of random qubits, a
 * Clifford+T circuit state, a product state, and a basis state on random
 * qubits with |+> on the rest.
 */
static void
rand_state(unsigned n, unsigned kind, EVBDD_WGT *f)
{
    switch (kind) {
    case S_ZEROS:     rand_on(n, rand_qubits(n), 0.2 * (double)(1 + rnd_below(3)), false, f); break;
    case S_ALGEBRAIC: rand_on(n, rand_qubits(n), 0.0, true, f); break;
    case S_PHASE:     phase_on(n, rand_qubits(n), f); break;
    case S_CIRCUIT:   circuit_state(n, f); break;
    case S_PRODUCT:   product_state(n, f); break;
    default: {
        const uint64_t D = imask(rand_qubits(n), n), x0 = rnd_below(UINT64_C(1) << n) & D;
        for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
            f[x] = ((x & D) == x0) ? EVBDD_ONE : EVBDD_ZERO;
    }
    }
}

/** In the families with them, a random sign pattern and translation: (-1)^{s.x} f(x ^ t). */
static void
twist(unsigned n, EVBDD_WGT *f)
{
    const bqd_family_t fam = bqd_family();
    const uint64_t len = UINT64_C(1) << n;
    const uint64_t s = (fam == BQD_FAMILY_PAULI) ? rnd_below(len) : 0;
    const uint64_t t = (fam != BQD_FAMILY_SCALAR) ? rnd_below(len) : 0;
    EVBDD_WGT *g = malloc(len * sizeof(EVBDD_WGT));
    memcpy(g, f, len * sizeof(EVBDD_WGT));
    for (uint64_t x = 0; x < len; x++) {
        EVBDD_WGT v = g[x ^ t];
        if (v != EVBDD_ZERO && parity(s & x)) v = wgt_neg(v);
        f[x] = v;
    }
    free(g);
}

/* --- building ------------------------------------------------------------- */

typedef struct {
    unsigned tot, decode, dense, eval, again;
    unsigned skip, root_skip, to_terminal;
} builds_t;

/** Build f, and check the build against f and the structure of the diagram. */
static BQD
build_checked(const EVBDD_WGT *f, unsigned n, builds_t *b)
{
    const uint64_t len = UINT64_C(1) << n;
    EVBDD_WGT *g = malloc(len * sizeof(EVBDD_WGT));
    const BQD e = bqd_from_vector(f, n);
    b->tot++;
    bqd_to_vector(e, n, g);
    if (!same_vector(f, g, n)) b->decode++;
    edge_dense(e, 0, n, g);
    if (!same_vector(f, g, n)) b->dense++;
    for (uint64_t x = 0; x < len; x++) if (bqd_eval(e, n, x) != f[x]) { b->eval++; break; }
    if (bqd_from_vector(f, n) != e) b->again++;
    check_diagram(e, f, n);
    bool to_terminal;
    if (skipped_qubits(e, n, &to_terminal) != 0) b->skip++;
    if (to_terminal) b->to_terminal++;
    if (!limdd_edge_is_zero(e) && level_of(limdd_target(e), n) > 0) b->root_skip++;
    free(g);
    return e;
}

static void
builds_report(const char *section, const builds_t *b)
{
    char what[96];
    snprintf(what, sizeof(what), "%s: decode(build(f)) = f", section);
    report(what, b->decode, b->tot);
    snprintf(what, sizeof(what), "%s: the decoder of skip:def:skip", section);
    report(what, b->dense, b->tot);
    snprintf(what, sizeof(what), "%s: bqd_eval at every index", section);
    report(what, b->eval, b->tot);
    snprintf(what, sizeof(what), "%s: built twice, one edge", section);
    report(what, b->again, b->tot);
    snprintf(what, sizeof(what), "%s: diagrams that skip a level", section);
    covered(what, b->skip);
    snprintf(what, sizeof(what), "%s: root edges that skip", section);
    covered(what, b->root_skip);
    snprintf(what, sizeof(what), "%s: edges that skip to the terminal", section);
    covered(what, b->to_terminal);
}

/**
 * The shapes the note fixes: |+>^n up to scale is the terminal under one
 * label, and so, in the Pauli family, is a character c (-1)^{s.x}, whose Z
 * sits on the root at levels it skips; in the other two a character has one
 * node at each qubit it reads; a basis state has a node at every level.
 */
static void
check_shapes(unsigned max_n)
{
    const bqd_family_t fam = bqd_family();
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT));
    size_t counts[NQ];
    unsigned bad_const = 0, bad_char = 0, bad_basis = 0, tc = 0, th = 0, tb = 0;
    for (unsigned n = 1; n <= max_n; n++) {
        const uint64_t len = UINT64_C(1) << n;
        const EVBDD_WGT cs[4] = { EVBDD_ONE, EVBDD_MIN_ONE, pw[2], wgt_mul(isq2, pw[3]) };
        for (int i = 0; i < 4; i++) {
            for (uint64_t x = 0; x < len; x++) f[x] = cs[i];
            tc++;
            if (bqd_from_vector(f, n) != limdd_bundle(bqd_lim_make(cs[i], 0, 0, n), LIMDD_TERMINAL))
                bad_const++;
        }
        for (int rep = 0; rep < 8; rep++) {
            const uint64_t s = rnd_below(len);
            const EVBDD_WGT c = pw[rnd_below(8)];
            for (uint64_t x = 0; x < len; x++) f[x] = parity(s & x) ? wgt_neg(c) : c;
            const BQD e = bqd_from_vector(f, n);
            th++;
            if (fam == BQD_FAMILY_PAULI) {
                if (e != limdd_bundle(bqd_lim_make(c, s, 0, n), LIMDD_TERMINAL)) bad_char++;
            } else {
                limdd_level_counts(e, counts, n);
                for (unsigned q = 0; q < n; q++)
                    if (counts[q] != ((s & ibit(q, n)) ? 1u : 0u)) { bad_char++; break; }
            }
            const uint64_t y = rnd_below(len);
            for (uint64_t x = 0; x < len; x++) f[x] = (x == y) ? EVBDD_ONE : EVBDD_ZERO;
            const BQD b = bqd_from_vector(f, n);
            tb++;
            limdd_level_counts(b, counts, n);
            bool ok = true;
            for (unsigned q = 0; q < n; q++) ok = ok && counts[q] == 1;
            ok = ok && b == bqd_basis_state(y, n);
            /* in the translation and Pauli families every basis state is X^x
             * on the one chain of |0...0> */
            if (fam != BQD_FAMILY_SCALAR)
                ok = ok && limdd_target(b) == limdd_target(bqd_basis_state(0, n));
            if (!ok) bad_basis++;
        }
    }
    report("|+>^n is the terminal under one label", bad_const, tc);
    report(fam == BQD_FAMILY_PAULI ? "a character is the terminal under c Z^s"
                                   : "a character has a node at each qubit it reads",
           bad_char, th);
    report(fam == BQD_FAMILY_SCALAR ? "a basis state: a node per level, bqd_basis_state"
                                    : "a basis state: X^x on the chain of |0...0>",
           bad_basis, tb);
    free(f);
}

/** Every vector over {0, 1, -1, i} on 1 to 3 qubits: 16, 256 and 65536 of them. */
static void
check_exhaustive(void)
{
    const EVBDD_WGT vals[4] = { EVBDD_ZERO, EVBDD_ONE, EVBDD_MIN_ONE, pw[2] };
    EVBDD_WGT f[8];
    builds_t b;
    memset(&b, 0, sizeof(b));
    structure_begin();
    for (unsigned n = 1; n <= 3; n++) {
        const uint64_t len = UINT64_C(1) << n;
        for (uint64_t idx = 0; idx < (UINT64_C(1) << (2 * len)); idx++) {
            for (uint64_t x = 0; x < len; x++) f[x] = vals[(idx >> (2 * x)) & 3];
            build_checked(f, n, &b);
        }
    }
    builds_report("n <= 3, all", &b);
    structure_report("n <= 3, all");
    if (bqd_family() != BQD_FAMILY_SCALAR)
        covered("n <= 3, all: edges with a translation", sf.translated);
    if (bqd_family() == BQD_FAMILY_PAULI)
        covered("n <= 3, all: a Z at a level its edge skips", sf.z_skipped);
}

static void
check_random_builds(void)
{
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT));
    builds_t b;
    memset(&b, 0, sizeof(b));
    structure_begin();
    for (unsigned n = 4; n <= NQ; n++) for (int rep = 0; rep < 8; rep++) {
        for (unsigned kind = 0; kind < S_KINDS; kind++) {
            rand_state(n, kind, f);
            if (rnd_below(2)) twist(n, f);
            build_checked(f, n, &b);
        }
    }
    builds_report("n = 4..8, random", &b);
    structure_report("n = 4..8, random");
    if (bqd_family() == BQD_FAMILY_PAULI)
        covered("n = 4..8, random: a Z at a level its edge skips", sf.z_skipped);
    check_shapes(NQ);
    free(f);
}

/* --- operations, scalar family --------------------------------------------- */

typedef struct { unsigned bad, tot; } tally_t;

/** Count got wrong unless it is the edge bqd_from_vector builds for h, and check its structure. */
static void
tally(tally_t *t, BQD got, const EVBDD_WGT *h, unsigned n)
{
    t->tot++;
    check_diagram(got, h, n);
    if (got != bqd_from_vector(h, n)) t->bad++;
}

/** Put U in the table as the dynamic gate. The gate memo is keyed on the id, so the cache goes. */
static void
set_dynamic(const EVBDD_WGT *U)
{
    for (int i = 0; i < 4; i++) gates[GATEID_dynamic][i] = U[i];
    sylvan_clear_cache();
}

static void
dense_restrict(const EVBDD_WGT *f, unsigned n, unsigned q, int b, EVBDD_WGT *out)
{
    const uint64_t qb = ibit(q, n);
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) out[x] = f[b ? (x | qb) : (x & ~qb)];
}

static void
dense_project(const EVBDD_WGT *f, unsigned n, unsigned q, int b, EVBDD_WGT *out)
{
    const uint64_t qb = ibit(q, n);
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
        out[x] = (((x & qb) != 0) == (b != 0)) ? f[x] : EVBDD_ZERO;
}

/** sum_c M[r][c] f|_{x_Q = c} at the points where x_Q = r, bit k-1-i of r and c being qubits[i]. */
static void
dense_matvec(const EVBDD_WGT *f, unsigned n, const EVBDD_WGT *M, const uint32_t *qs, unsigned k,
             EVBDD_WGT *out)
{
    const uint32_t dim = 1u << k;
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) {
        uint32_t r = 0;
        uint64_t base = x;
        for (unsigned i = 0; i < k; i++) {
            if (x & ibit(qs[i], n)) r |= 1u << (k - 1 - i);
            base &= ~ibit(qs[i], n);
        }
        EVBDD_WGT acc = EVBDD_ZERO;
        for (uint32_t c = 0; c < dim; c++) {
            const EVBDD_WGT m = M[(size_t)r * dim + c];
            if (m == EVBDD_ZERO) continue;
            uint64_t y = base;
            for (unsigned i = 0; i < k; i++) if ((c >> (k - 1 - i)) & 1) y |= ibit(qs[i], n);
            acc = wgt_add(acc, wgt_mul(m, f[y]));
        }
        out[x] = acc;
    }
}

enum {
    T_MUL, T_ADD, T_CANCEL, T_SCALE, T_NEGATE, T_COMPOSE, T_COF, T_RESTRICT, T_PROJECT,
    T_GATE, T_DYN, T_CGATE, T_CPHASE, T_EITHER, T_SWAP, T_MATVEC1, T_MATVEC2, T_PRODUCT,
    T_PRODFULL, T_FULL, T_COUNT
};
static const char *tally_name[T_COUNT] = {
    "bqd_multiply", "bqd_add", "bqd_add with cancellation", "bqd_scale", "bqd_negate",
    "bqd_compose (OP_XQUOT inside)", "bqd_cofactor (OP_XPROD on the high side)",
    "bqd_restrict, every qubit", "bqd_project, every qubit", "bqd_gate, table gates, every qubit",
    "bqd_gate, random exact 2x2, every qubit", "bqd_cgate, controls above the target",
    "bqd_cgate, phase with controls anywhere", "bqd_cgate_either", "bqd_swap",
    "bqd_local_matvec on one qubit", "bqd_local_matvec on two qubits", "bqd_product, any support",
    "bqd_product, full support (prop:prodscalar)", "bqd_has_full_support",
};

enum {
    C_ROOT, C_TERMINAL, C_LEVELS, C_COMPOSE, C_COF, C_RESTRICT, C_GATE, C_GATE_FREE, C_CONTROL,
    C_TARGET, C_BELOW, C_EITHER, C_MATVEC, C_PROD_LEVELS, C_PROD_SKIP, C_PROD_CANCEL,
    C_ZERO_LOW, C_CANCEL, C_COUNT
};
static const char *cover_name[C_COUNT] = {
    "states whose root edge skips", "states with an edge that skips to the terminal",
    "operand pairs at different root levels", "compose of two equal edges",
    "cofactor at a level the edge skips", "restrict or project on a skipped qubit",
    "gates on a skipped target", "gates on a qubit the state ignores",
    "controlled gates with a skipped control", "controlled gates with a skipped target",
    "phases with a control below the target", "cgate_either with the control below",
    "local matvecs on a skipped qubit", "full-support products at different root levels",
    "full-support products that skip a level",
    "full-support products that drop a qubit both read",
    "gates on a qubit whose low cofactor is zero", "gates that cancel an amplitude",
};

static LIMDD_LIM rand_label(unsigned n, bool alg);
static BQD relabel(BQD e, LIMDD_LIM lim, unsigned n, EVBDD_WGT *out);

/**
 * What a gate on qubit q took f to h through: a zero low cofactor, the
 * support x_q = 1 alone, and a cancellation, a zero of h where neither
 * cofactor of f is zero.
 */
static void
gate_cases(const EVBDD_WGT *f, const EVBDD_WGT *h, unsigned n, unsigned q, unsigned *zero_low,
           unsigned *cancels)
{
    const uint64_t qb = ibit(q, n);
    bool low_zero = true, nonzero = false, cancel = false;
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) {
        if (x & qb) continue;
        if (f[x] != EVBDD_ZERO) low_zero = false;
        if (f[x | qb] != EVBDD_ZERO) nonzero = true;
        if (f[x] != EVBDD_ZERO && f[x | qb] != EVBDD_ZERO
            && (h[x] == EVBDD_ZERO || h[x | qb] == EVBDD_ZERO)) cancel = true;
    }
    if (low_zero && nonzero) (*zero_low)++;
    if (cancel) (*cancels)++;
}

/**
 * Gate, Restrict and Project of qsylvan_bqd_xp.h on `draws` random labelled
 * edges on F's node, which is what their recursions hand themselves: a
 * cofactor is a label product, canonical or not, and each of the three
 * passes such an edge through Canon wherever the scalar family returns an
 * operand or a cofactor as it is. The gate has random controls above a
 * random target, or none, and is any of the table's or the dynamic one.
 * Tallied in t[0], t[1] and t[2].
 */
static void
xp_labelled_gates(BQD F, unsigned n, int draws, tally_t *t)
{
    static const uint32_t gs[] = { GATEID_H, GATEID_X, GATEID_Y, GATEID_Z, GATEID_S,
                                   GATEID_T, GATEID_sqrtX, GATEID_dynamic };
    const size_t vb = (UINT64_C(1) << n) * sizeof(EVBDD_WGT);
    EVBDD_WGT *lf = malloc(vb), *h = malloc(vb);
    for (int k = 0; k < draws; k++) {
        const BQD L = relabel(F, rand_label(n, false), n, lf);
        const uint32_t q = (uint32_t)rnd_below(n);
        const uint64_t cm = rnd_below(UINT64_C(1) << q);
        const uint32_t gid = gs[rnd_below(sizeof(gs) / sizeof(gs[0]))];
        memcpy(h, lf, vb);
        dense_gate(h, n, gates[gid], cm, q);
        tally(&t[0], bqd_xp_cgate_rec(L, gid, cm, q), h, n);
        const uint32_t r = (uint32_t)rnd_below(n);
        const int b = (int)rnd_below(2);
        dense_restrict(lf, n, r, b, h);
        tally(&t[1], bqd_xp_restrict(L, r, b), h, n);
        dense_project(lf, n, r, b, h);
        tally(&t[2], bqd_xp_project(L, r, b), h, n);
    }
    free(lf); free(h);
}

static void
check_operations(void)
{
    const bqd_family_t fam = bqd_family();
    static const uint32_t table[] = { GATEID_H, GATEID_X, GATEID_Y, GATEID_Z, GATEID_S,
                                      GATEID_T, GATEID_Tdag, GATEID_sqrtX, GATEID_sqrtY,
                                      GATEID_proj0, GATEID_proj1 };
    static const uint32_t ctl[] = { GATEID_X, GATEID_Y, GATEID_Z, GATEID_H, GATEID_S,
                                    GATEID_T, GATEID_sqrtX, GATEID_dynamic };
    static const uint32_t ph[] = { GATEID_Z, GATEID_S, GATEID_T, GATEID_Sdag };
    const size_t vb = MAXV * sizeof(EVBDD_WGT);
    EVBDD_WGT *f = malloc(vb), *g = malloc(vb), *h = malloc(vb), *g0 = malloc(vb), *g1 = malloc(vb);
    tally_t t[T_COUNT], lab[3];
    unsigned cov[C_COUNT];
    memset(t, 0, sizeof(t));
    memset(lab, 0, sizeof(lab));
    memset(cov, 0, sizeof(cov));
    structure_begin();

    for (unsigned n = 1; n <= 7; n++) for (int rep = 0; rep < 12; rep++) {
        const uint64_t len = UINT64_C(1) << n;
        rand_state(n, (unsigned)rnd_below(S_KINDS), f);
        if (fam != BQD_FAMILY_SCALAR && rnd_below(2)) twist(n, f);
        rand_state(n, (unsigned)rnd_below(S_KINDS), g);
        if (fam != BQD_FAMILY_SCALAR && rnd_below(2)) twist(n, g);
        const BQD F = bqd_from_vector(f, n), G = bqd_from_vector(g, n);
        check_diagram(F, f, n);
        check_diagram(G, g, n);
        /* full support is kept by extension, so a skipping edge changes nothing */
        bool full = true;
        for (uint64_t x = 0; x < len; x++) if (f[x] == EVBDD_ZERO) full = false;
        t[T_FULL].tot++;
        if (bqd_has_full_support(F) != full) t[T_FULL].bad++;
        bool to_terminal;
        const uint64_t skF = skipped_qubits(F, n, &to_terminal);
        const uint64_t dF = depends_on(f, n);
        if (to_terminal) cov[C_TERMINAL]++;
        if (!limdd_edge_is_zero(F) && level_of(limdd_target(F), n) > 0) cov[C_ROOT]++;

        /* pointwise, and the two that only touch the scalar */
        if (level_of(limdd_target(F), n) != level_of(limdd_target(G), n)) cov[C_LEVELS]++;
        for (uint64_t x = 0; x < len; x++) h[x] = wgt_mul(f[x], g[x]);
        tally(&t[T_MUL], bqd_multiply(F, G), h, n);
        for (uint64_t x = 0; x < len; x++) h[x] = wgt_add(f[x], g[x]);
        tally(&t[T_ADD], bqd_add(F, G), h, n);
        for (uint64_t x = 0; x < len; x++) g1[x] = rnd_below(3) ? wgt_neg(f[x]) : g[x];
        for (uint64_t x = 0; x < len; x++) h[x] = wgt_add(f[x], g1[x]);
        tally(&t[T_CANCEL], bqd_add(F, bqd_from_vector(g1, n)), h, n);
        const EVBDD_WGT c = rand_entry();
        for (uint64_t x = 0; x < len; x++) h[x] = wgt_mul(c, f[x]);
        tally(&t[T_SCALE], bqd_scale(F, c), h, n);
        for (uint64_t x = 0; x < len; x++) h[x] = wgt_neg(f[x]);
        tally(&t[T_NEGATE], bqd_negate(F), h, n);

        /* compose at var, of two functions of the qubits below it: equal, a
         * multiple, one zero, the other zero, or unrelated */
        for (int k = 0; k < 3; k++) {
            const unsigned var = (unsigned)rnd_below(n);
            const uint64_t below = ((UINT64_C(1) << n) - 1) & ~((UINT64_C(2) << var) - 1);
            rand_on(n, rand_qubits(n) & below, 0.3 * (double)rnd_below(2), rnd_below(2), g0);
            rand_on(n, rand_qubits(n) & below, 0.3 * (double)rnd_below(2), rnd_below(2), g1);
            const EVBDD_WGT m = rand_value(0.0, true);
            const unsigned rel = (unsigned)rnd_below(5);
            for (uint64_t x = 0; x < len; x++) {
                switch (rel) {
                case 0: g1[x] = g0[x]; break;
                case 1: g1[x] = wgt_mul(m, g0[x]); break;
                case 2: g1[x] = EVBDD_ZERO; break;
                case 3: g0[x] = EVBDD_ZERO; break;
                default: break;
                }
            }
            for (uint64_t x = 0; x < len; x++) h[x] = (x & ibit(var, n)) ? g1[x] : g0[x];
            const BQD e0 = bqd_from_vector(g0, n), e1 = bqd_from_vector(g1, n);
            if (e0 == e1) cov[C_COMPOSE]++;
            /* both read at var + 1, so a top at var or above is a build that
             * kept a level it should have skipped, and compose may not have it */
            if (top_of(e0, n) <= var || top_of(e1, n) <= var) {
                t[T_COMPOSE].tot++;
                t[T_COMPOSE].bad++;
                continue;
            }
            tally(&t[T_COMPOSE], bqd_compose(var, e0, e1), h, n);
        }

        /* the cofactors of F at every level down to its top */
        const unsigned lF = top_of(F, n);
        for (unsigned var = 0; var <= lF && var < n; var++) for (int b = 0; b < 2; b++) {
            if (var < lF && !limdd_edge_is_zero(F)) cov[C_COF]++;
            dense_restrict(f, n, var, b, h);
            tally(&t[T_COF], bqd_cofactor(F, var, b), h, n);
        }

        for (unsigned q = 0; q < n; q++) for (int b = 0; b < 2; b++) {
            if ((skF >> q) & 1) cov[C_RESTRICT]++;
            dense_restrict(f, n, q, b, h);
            tally(&t[T_RESTRICT], bqd_restrict(F, q, b), h, n);
            dense_project(f, n, q, b, h);
            tally(&t[T_PROJECT], bqd_project(F, q, b), h, n);
        }

        for (size_t gi = 0; gi < sizeof(table) / sizeof(table[0]); gi++) for (unsigned q = 0; q < n; q++) {
            if ((skF >> q) & 1) cov[C_GATE]++;
            if (!((dF >> q) & 1)) cov[C_GATE_FREE]++;
            memcpy(h, f, len * sizeof(EVBDD_WGT));
            dense_gate(h, n, gates[table[gi]], 0, q);
            gate_cases(f, h, n, q, &cov[C_ZERO_LOW], &cov[C_CANCEL]);
            tally(&t[T_GATE], bqd_gate(F, table[gi], q, n), h, n);
        }

        /* a random exact matrix, as the dynamic gate and as a local matvec */
        EVBDD_WGT U[4];
        for (int i = 0; i < 4; i++) U[i] = rand_entry();
        set_dynamic(U);
        for (uint32_t q = 0; q < n; q++) {
            memcpy(h, f, len * sizeof(EVBDD_WGT));
            dense_gate(h, n, U, 0, q);
            tally(&t[T_DYN], bqd_gate(F, GATEID_dynamic, q, n), h, n);
            if ((skF >> q) & 1) cov[C_MATVEC]++;
            tally(&t[T_MATVEC1], bqd_local_matvec(F, U, &q, 1, n), h, n);
        }
        if (fam != BQD_FAMILY_SCALAR) xp_labelled_gates(F, n, 8, lab);

        if (n >= 2) {
            for (int k = 0; k < 3; k++) {
                /* controls above the target, any gate */
                const unsigned q = 1 + (unsigned)rnd_below(n - 1);
                uint64_t cm = 0;
                while (cm == 0) cm = rnd_below(UINT64_C(1) << q);
                const uint32_t gid = ctl[rnd_below(sizeof(ctl) / sizeof(ctl[0]))];
                if (cm & skF) cov[C_CONTROL]++;
                if ((skF >> q) & 1) cov[C_TARGET]++;
                memcpy(h, f, len * sizeof(EVBDD_WGT));
                dense_gate(h, n, gates[gid], cm, q);
                tally(&t[T_CGATE], bqd_cgate(F, gid, cm, q, n), h, n);

                /* a phase, its controls anywhere */
                const uint32_t pg = ph[rnd_below(4)];
                const unsigned tq = (unsigned)rnd_below(n);
                uint64_t pm = 0;
                while (pm == 0) pm = rnd_below(UINT64_C(1) << n) & ~(UINT64_C(1) << tq);
                if (pm >> (tq + 1)) cov[C_BELOW]++;
                memcpy(h, f, len * sizeof(EVBDD_WGT));
                dense_gate(h, n, gates[pg], pm, tq);
                tally(&t[T_CPHASE], bqd_cgate(F, pg, pm, tq, n), h, n);
            }

            /* one control on either side, and a swap */
            const unsigned a = (unsigned)rnd_below(n);
            unsigned b = (unsigned)rnd_below(n - 1); if (b >= a) b++;
            const uint32_t eg = rnd_below(2) ? GATEID_X : GATEID_Z;
            if (a > b) cov[C_EITHER]++;
            memcpy(h, f, len * sizeof(EVBDD_WGT));
            dense_gate(h, n, gates[eg], UINT64_C(1) << a, b);
            bool ok;
            const BQD r = bqd_cgate_either(F, eg, a, b, n, &ok);
            tally(&t[T_EITHER], r, h, n);
            if (!ok) t[T_EITHER].bad++;
            memcpy(h, f, len * sizeof(EVBDD_WGT));
            dense_swap(h, n, a, b);
            tally(&t[T_SWAP], bqd_swap(F, a, b, n), h, n);

            /* two qubits, in either order */
            const uint32_t qs[2] = { a, b };
            EVBDD_WGT M[16];
            for (int i = 0; i < 16; i++) M[i] = rand_entry();
            if (skF & ((UINT64_C(1) << a) | (UINT64_C(1) << b))) cov[C_MATVEC]++;
            dense_matvec(f, n, M, qs, 2, h);
            tally(&t[T_MATVEC2], bqd_local_matvec(F, M, qs, 2, n), h, n);
        }

        /* products: the operands as drawn, and two phase states of full support */
        for (uint64_t x = 0; x < len; x++) h[x] = wgt_mul(f[x], g[x]);
        tally(&t[T_PRODUCT], bqd_product(F, G), h, n);
        if (bqd_product(G, F) != bqd_product(F, G)) t[T_PRODUCT].bad++;
        phase_on(n, rand_qubits(n), g0);
        phase_on(n, rand_qubits(n), g1);
        for (uint64_t x = 0; x < len; x++) h[x] = wgt_mul(g0[x], g1[x]);
        const BQD P0 = bqd_from_vector(g0, n), P1 = bqd_from_vector(g1, n);
        if (level_of(limdd_target(P0), n) != level_of(limdd_target(P1), n)) cov[C_PROD_LEVELS]++;
        const BQD P = bqd_product(P0, P1);
        if (skipped_qubits(P, n, &to_terminal) != 0) cov[C_PROD_SKIP]++;
        tally(&t[T_PRODFULL], P, h, n);
        if (bqd_product(P1, P0) != P) t[T_PRODFULL].bad++;
        /* and one whose second operand undoes the first but for a phase of
         * other qubits: the product ignores qubits both operands read, so its
         * ratio there is the constant one and the level goes (skip:cor:full) */
        phase_on(n, rand_qubits(n), g);
        for (uint64_t x = 0; x < len; x++) g1[x] = wgt_mul(wgt_div(EVBDD_ONE, g0[x]), g[x]);
        const uint64_t both = depends_on(g0, n) & depends_on(g1, n) & ~depends_on(g, n);
        if (both != 0) cov[C_PROD_CANCEL]++;
        const BQD Q = bqd_product(P0, bqd_from_vector(g1, n));
        tally(&t[T_PRODFULL], Q, g, n);
    }

    for (int i = 0; i < T_COUNT; i++) report(tally_name[i], t[i].bad, t[i].tot);
    if (fam != BQD_FAMILY_SCALAR) {
        report("Gate, labelled operands (bqd_xp_cgate_rec)", lab[0].bad, lab[0].tot);
        report("Restrict, labelled operands (bqd_xp_restrict)", lab[1].bad, lab[1].tot);
        report("Project, labelled operands (bqd_xp_project)", lab[2].bad, lab[2].tot);
    }
    for (int i = 0; i < C_COUNT; i++) covered(cover_name[i], cov[i]);
    structure_report("operations");
    free(f); free(g); free(h); free(g0); free(g1);
}

/* --- skip:prop:diag ------------------------------------------------------ */

/**
 * The nodes of e the diagonal walk visits, from the diagram alone
 * (skip:prop:diag): the path from the root that goes low at a node outside A
 * and high at a node in A, and ends after the node of the last qubit of A, or
 * at the first edge that skips the first qubit of A not yet passed, whose
 * node is not counted; *stopped says whether it ended there. A is a set of
 * qubits, bit q for qubit q.
 */
static uint32_t
diag_path(BQD e, uint64_t A, unsigned n, bool *stopped)
{
    uint32_t visits = 0;
    LIMDD_TARG t = limdd_target(e);
    *stopped = false;
    while (A != 0) {
        const unsigned a = (unsigned)__builtin_ctzll(A);    /* the top qubit of A left */
        const unsigned v = level_of(t, n);
        if (v > a) { *stopped = true; break; }
        visits++;
        if (v < a) { t = limdd_target(limdd_node_low(t)); continue; }
        A &= A - 1;
        if (A != 0) t = limdd_target(limdd_node_high(t));
    }
    return visits;
}

/** The Z mask of e's label, over the qubits of the vector-index mask A. */
static uint64_t
z_on(BQD e, uint64_t A, unsigned n)
{
    EVBDD_WGT c;
    uint64_t z, x;
    if (limdd_edge_is_zero(e)) return 0;
    bqd_lim_masks(limdd_label(e), n, &c, &z, &x);
    return z & A;
}

static void
check_diagonal(void)
{
    /* the walk serves every family: the scalar one's, which the translation
     * family's diagram of full support is node for node, and the Pauli
     * family's with label products, where a monomial may have fewer nodes */
    const bqd_family_t fam = bqd_family();
    const bool pauli = fam == BQD_FAMILY_PAULI;
    const size_t vb = MAXV * sizeof(EVBDD_WGT);
    EVBDD_WGT *f = malloc(vb), *h = malloc(vb);
    size_t counts[NQ];
    tally_t walk = { 0, 0 }, undo = { 0, 0 }, prod = { 0, 0 }, mono = { 0, 0 }, fallback = { 0, 0 };
    unsigned badvisits = 0, badshape = 0, badones = 0;
    unsigned at_bound = 0, stops = 0, outside = 0, root_skip = 0, to_term = 0, undone = 0;
    unsigned fewer = 0, signed_a = 0;
    structure_begin();

    for (unsigned n = 1; n <= NQ; n++) for (int rep = 0; rep < 40; rep++) {
        const uint64_t len = UINT64_C(1) << n;
        const uint64_t D = rand_qubits(n);
        phase_on(n, D, f);
        if (fam != BQD_FAMILY_SCALAR && rnd_below(2)) twist(n, f);
        const unsigned arity = 1 + (unsigned)rnd_below(n < 3 ? n : 3);
        uint64_t Aq = 0;
        while (__builtin_popcountll(Aq) < (int)arity) Aq |= UINT64_C(1) << rnd_below(n);
        const uint64_t A = imask(Aq, n);
        const EVBDD_WGT phase = pw[1 + rnd_below(7)];
        if (Aq & ~D) outside++;

        /* the monomial is one node at each qubit of A, and none elsewhere; in
         * the Pauli family at most one, and none at the last qubit of A for a
         * phase of -1, which is a Z there (skip:alg:xrestrict) */
        for (uint64_t x = 0; x < len; x++) h[x] = ((x & A) == A) ? phase : EVBDD_ONE;
        const BQD gate = bqd_monomial(A, phase, n);
        tally(&mono, gate, h, n);
        limdd_level_counts(gate, counts, n);
        unsigned nodes = 0;
        for (unsigned q = 0; q < n; q++) {
            const size_t in_A = (Aq >> q) & 1;
            nodes += (unsigned)counts[q];
            if (pauli ? counts[q] > in_A : counts[q] != in_A) { badshape++; break; }
        }
        if (nodes < (unsigned)__builtin_popcountll(Aq)) fewer++;
        if (bqd_monomial(A, EVBDD_ONE, n) != limdd_one_edge()) badones++;
        for (unsigned v = 0; v <= n; v++) if (bqd_ones(v, n) != LIMDD_TERMINAL) { badones++; break; }

        const BQD psi = bqd_from_vector(f, n);
        check_diagram(psi, f, n);
        bool term;
        (void)skipped_qubits(psi, n, &term);
        if (term) to_term++;
        if (level_of(limdd_target(psi), n) > 0) root_skip++;
        for (uint64_t x = 0; x < len; x++) h[x] = ((x & A) == A) ? wgt_mul(f[x], phase) : f[x];
        uint32_t visits;
        const BQD got = bqd_apply_diagonal_counted(psi, A, phase, n, &visits);
        tally(&walk, got, h, n);
        prod.tot++;
        if (bqd_product(psi, gate) != got) prod.bad++;
        /* a Z at a qubit of A that the state had not: the sign repair, or a
         * level of A skipped with a sign */
        if (z_on(got, A, n) & ~z_on(psi, A, n)) signed_a++;

        /* at most one node per level down to the last qubit of A, and
         * exactly the nodes of the path, which stops where psi skips a qubit
         * of A */
        bool stopped;
        const uint32_t path = diag_path(psi, Aq, n, &stopped);
        const uint32_t bound = n - (uint32_t)__builtin_ctzll(A);
        if (visits > bound || visits != path) badvisits++;
        if (stopped) stops++;
        if (visits == bound) at_bound++;

        /* the gate undoing a phase of the state on A: the result ignores the
         * qubits of A that only that phase read, and the walk has to drop
         * their levels (skip:cor:full) */
        for (uint64_t x = 0; x < len; x++) h[x] = ((x & A) == A) ? wgt_div(f[x], phase) : f[x];
        if (depends_on(h, n) & ~depends_on(f, n)) undone++;
        const BQD chi = bqd_from_vector(h, n);
        tally(&undo, bqd_apply_diagonal_counted(chi, A, phase, n, &visits), f, n);
        if (visits > bound || visits != diag_path(chi, Aq, n, &stopped)) badvisits++;

        /* a state with a zero takes the general product, and walks nothing */
        rand_on(n, rand_qubits(n), 0.3, true, f);
        if (fam != BQD_FAMILY_SCALAR && rnd_below(2)) twist(n, f);
        for (uint64_t x = 0; x < len; x++) h[x] = ((x & A) == A) ? wgt_mul(f[x], phase) : f[x];
        const BQD z = bqd_from_vector(f, n);
        tally(&fallback, bqd_apply_diagonal_counted(z, A, phase, n, &visits), h, n);
        if (!bqd_has_full_support(z) && visits != 0) fallback.bad++;
    }

    /* every vector over {1, -1, i, w_8} on 1 and 2 qubits, and one in 64 on
     * 3, under every monomial with every phase but 1: the walk against the
     * edge bqd_from_vector builds, and its visits against the path */
    tally_t every = { 0, 0 };
    const EVBDD_WGT fv[4] = { EVBDD_ONE, EVBDD_MIN_ONE, pw[2], pw[1] };
    for (unsigned n = 1; n <= 3; n++) {
        const uint64_t len = UINT64_C(1) << n, count = UINT64_C(1) << (2 * len);
        for (uint64_t idx = 0; idx < count; idx++) {
            if (n == 3 && rnd_below(64 * THIN) != 0) continue;
            for (uint64_t x = 0; x < len; x++) f[x] = fv[(idx >> (2 * x)) & 3];
            const BQD psi = bqd_from_vector(f, n);
            for (uint64_t Aq = 1; Aq < len; Aq++) for (int k = 1; k < 8; k++) {
                const uint64_t A = imask(Aq, n);
                for (uint64_t x = 0; x < len; x++) h[x] = ((x & A) == A) ? wgt_mul(f[x], pw[k]) : f[x];
                uint32_t visits;
                bool stopped;
                const BQD got = bqd_apply_diagonal_counted(psi, A, pw[k], n, &visits);
                tally(&every, got, h, n);
                if (visits != diag_path(psi, Aq, n, &stopped)) every.bad++;
                if (z_on(got, A, n) & ~z_on(psi, A, n)) signed_a++;
            }
        }
    }

    report("bqd_apply_diagonal_counted, full support", walk.bad, walk.tot);
    report("bqd_apply_diagonal_counted, every monomial on n <= 3", every.bad, every.tot);
    report("bqd_apply_diagonal_counted, the gate undoing a phase", undo.bad, undo.tot);
    report("the same edge as bqd_product with the monomial", prod.bad, prod.tot);
    report("bqd_apply_diagonal_counted, with zeros", fallback.bad, fallback.tot);
    report("bqd_monomial", mono.bad, mono.tot);
    report(pauli ? "the monomial is at most one node at each qubit of A"
                 : "the monomial is one node at each qubit of A", badshape, mono.tot);
    report("bqd_monomial of phase 1 and bqd_ones are the terminal", badones, mono.tot);
    report("visits: the path of skip:prop:diag, <= n - ctz(A)", badvisits, walk.tot + undo.tot);
    covered("walks that stop at a qubit of A the state skips", stops);
    covered("walks at the bound n - ctz(A)", at_bound);
    if (pauli) {
        covered("monomials with fewer nodes than qubits", fewer);
        covered("walks that put a Z on a qubit of A (repair or signed skip)", signed_a);
    }
    covered("monomials on a qubit the state ignores", outside);
    covered("gates that make the state ignore a qubit it read", undone);
    covered("walks on a state whose root edge skips", root_skip);
    covered("walks on a state with an edge that skips to the terminal", to_term);
    structure_report("diagonal");
    free(f); free(h);
}

/* --- operations, translation and Pauli families ---------------------------- */

/*
 * The operations of qsylvan_bqd_xp.h, through the entry points where there is
 * one: every result the edge bqd_from_vector builds for the dense result, and
 * a diagram that passes the structure checks. Apply and Canon are handed
 * labelled edges by the recursion, which need not be canonical or even
 * normal, so they get those here too: a random label of the family on the
 * node of a canonical edge.
 */

/**
 * A random label of the session's family on n qubits, vector-index masks,
 * whose scalar is a power of w_8, or a third of the time a random algebraic
 * number when `alg`.
 */
static LIMDD_LIM
rand_label(unsigned n, bool alg)
{
    const uint64_t len = UINT64_C(1) << n;
    const EVBDD_WGT c = (alg && rnd_below(3) == 0) ? rand_alg() : pw[rnd_below(8)];
    const uint64_t s = (bqd_family() == BQD_FAMILY_PAULI) ? rnd_below(len) : 0;
    return bqd_lim_make(c, s, rnd_below(len), n);
}

/** The labelled edge lim on the node of e, and in out the function it denotes. */
static BQD
relabel(BQD e, LIMDD_LIM lim, unsigned n, EVBDD_WGT *out)
{
    const uint64_t len = UINT64_C(1) << n;
    if (limdd_edge_is_zero(e)) { zero_fill(out, len); return e; }
    EVBDD_WGT *u = malloc(len * sizeof(EVBDD_WGT));
    ext_dense(limdd_target(e), 0, n, u);
    apply_label(lim, n, u, len, out);
    free(u);
    return limdd_bundle(lim, limdd_target(e));
}

static void
dense_pointwise(int op, const EVBDD_WGT *f, const EVBDD_WGT *g, unsigned n, EVBDD_WGT *h)
{
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) {
        const EVBDD_WGT a = f[x], b = g[x];
        switch (op) {
        case BQD_OP_MUL:   h[x] = wgt_mul(a, b); break;
        case BQD_OP_ADD:   h[x] = wgt_add(a, b); break;
        case BQD_OP_XPROD: h[x] = (a == EVBDD_ZERO) ? b : wgt_mul(a, b); break;
        default:           h[x] = (b == EVBDD_ZERO) ? a : wgt_div(a, b); break;
        }
    }
}

/** op by its entry point: bqd_multiply, bqd_add, and Apply for the two of the quotient rule. */
static BQD
pointwise(int op, BQD f, BQD g)
{
    switch (op) {
    case BQD_OP_MUL: return bqd_multiply(f, g);
    case BQD_OP_ADD: return bqd_add(f, g);
    default:         return bqd_xp_apply_op(op, f, g);
    }
}

enum {
    X_MUL, X_ADD, X_XPROD, X_XQUOT, X_LMUL, X_LADD, X_LXPROD, X_LXQUOT, X_CANCEL, X_SCALE,
    X_NEGATE, X_COF, X_COMPOSE, X_CANON, X_COUNT
};
static const char *xtally_name[X_COUNT] = {
    "bqd_multiply", "bqd_add", "xprod (bqd_xp_apply_op)", "xquot (bqd_xp_apply_op)",
    "bqd_multiply, labelled operands", "bqd_add, labelled operands",
    "xprod, labelled operands", "xquot, labelled operands", "bqd_add with cancellation",
    "bqd_scale", "bqd_negate", "bqd_cofactor, every level down to the top", "bqd_compose",
    "bqd_xp_canon of a labelled edge",
};

enum {
    XC_EQUAL, XC_OPPOSITE, XC_ZERO_LOW, XC_ZERO_HIGH, XC_GENERAL, XC_REPAIR, XC_SHIFT,
    XC_NONCANON, XC_ZTOP, XC_COF_SKIP, XC_COUNT
};
static const char *xcover_name[XC_COUNT] = {
    "compose (i): equal cofactors", "compose (ii): opposite cofactors",
    "compose (iii): a zero low cofactor", "compose (iv): a zero high cofactor",
    "compose (v): the general case", "compose (v) with the sign repair",
    "compose (v) with a translation on the high edge", "labelled edges that are not canonical",
    "labelled operands with a Z above their node", "cofactors at a level the edge skips",
};

/**
 * Which case of skip:prop:xcompose made r from e0 and e1 at var, read off
 * the result: (ii) keeps a's node, and (v) makes a node at var whose root
 * label has a Z at var where the sign was repaired.
 */
static void
compose_case(BQD e0, BQD e1, BQD r, unsigned var, unsigned n, unsigned *cov)
{
    if (e0 == e1) { cov[XC_EQUAL]++; return; }
    if (limdd_edge_is_zero(e0)) { cov[XC_ZERO_LOW]++; return; }
    if (limdd_edge_is_zero(e1)) { cov[XC_ZERO_HIGH]++; return; }
    if (limdd_target(r) == limdd_target(e0)) { cov[XC_OPPOSITE]++; return; }
    if (limdd_target(r) == LIMDD_TERMINAL || limdd_node_var(limdd_target(r)) != var) return;
    cov[XC_GENERAL]++;
    EVBDD_WGT c; uint64_t s, t;
    bqd_lim_masks(limdd_label(r), n, &c, &s, &t);
    if (s & ibit(var, n)) cov[XC_REPAIR]++;
    bqd_lim_masks(limdd_label(limdd_node_high(limdd_target(r))), n, &c, &s, &t);
    if (t != 0) cov[XC_SHIFT]++;
}

/**
 * The operations on one operand F, the diagram of f: scale, negate, the
 * cofactors at every level down to F's top, and Canon of labelled edges on
 * F's node, every label of the family with the scalars 1, -1 and w when
 * `all`, and `labels` random ones otherwise. Random algebraic scalars when
 * `alg`, and powers of w_8 otherwise, which keep the values of the
 * exhaustive runs few enough that their nodes repeat.
 */
static void
xp_unary(BQD F, const EVBDD_WGT *f, unsigned n, bool all, int labels, bool alg, tally_t *t,
         unsigned *cov)
{
    const uint64_t len = UINT64_C(1) << n;
    EVBDD_WGT *h = malloc(len * sizeof(EVBDD_WGT));
    const EVBDD_WGT c = alg ? rand_entry() : rnd_below(8) ? pw[rnd_below(8)] : EVBDD_ZERO;
    for (uint64_t x = 0; x < len; x++) h[x] = wgt_mul(c, f[x]);
    tally(&t[X_SCALE], bqd_scale(F, c), h, n);
    for (uint64_t x = 0; x < len; x++) h[x] = wgt_neg(f[x]);
    tally(&t[X_NEGATE], bqd_negate(F), h, n);

    /* read at var + 1 and so at 0 as well, a cofactor is the restriction */
    const unsigned top = bqd_xp_top(F);
    for (unsigned var = 0; var <= top && var < n; var++) for (int b = 0; b < 2; b++) {
        if (var < top && !limdd_edge_is_zero(F)) cov[XC_COF_SKIP]++;
        dense_restrict(f, n, var, b, h);
        tally(&t[X_COF], bqd_cofactor(F, var, b), h, n);
    }

    const bool pauli = bqd_family() == BQD_FAMILY_PAULI;
    const EVBDD_WGT cs[3] = { EVBDD_ONE, EVBDD_MIN_ONE, pw[1] };
    const uint64_t count = all ? 3 * len * (pauli ? len : 1) : (uint64_t)labels;
    for (uint64_t k = 0; k < count; k++) {
        const LIMDD_LIM lim = !all ? rand_label(n, alg)
                            : bqd_lim_make(cs[k % 3], pauli ? (k / 3) / len : 0, (k / 3) % len, n);
        const BQD L = relabel(F, lim, n, h);
        const BQD r = bqd_xp_canon(L);
        if (r != L) cov[XC_NONCANON]++;
        tally(&t[X_CANON], r, h, n);
    }
    free(h);
}

/**
 * The four pointwise operations on F and G, and when `labelled` on random
 * labelled edges on their nodes, with algebraic scalars when `alg`.
 */
static void
xp_binary(BQD F, const EVBDD_WGT *f, BQD G, const EVBDD_WGT *g, unsigned n, bool labelled,
          bool alg, tally_t *t, unsigned *cov)
{
    const size_t vb = (UINT64_C(1) << n) * sizeof(EVBDD_WGT);
    EVBDD_WGT *h = malloc(vb), *lf = malloc(vb), *lg = malloc(vb);
    for (int op = BQD_OP_MUL; op <= BQD_OP_XQUOT; op++) {
        dense_pointwise(op, f, g, n, h);
        tally(&t[X_MUL + op], pointwise(op, F, G), h, n);
    }
    if (!labelled) { free(h); free(lf); free(lg); return; }
    const BQD LF = relabel(F, rand_label(n, alg), n, lf), LG = relabel(G, rand_label(n, alg), n, lg);
    if (!limdd_edge_is_zero(LF) && bqd_xp_top(LF) < limdd_level(limdd_target(LF))) cov[XC_ZTOP]++;
    for (int op = BQD_OP_MUL; op <= BQD_OP_XQUOT; op++) {
        dense_pointwise(op, lf, lg, n, h);
        tally(&t[X_LMUL + op], pointwise(op, LF, LG), h, n);
    }
    free(h); free(lf); free(lg);
}

/**
 * Compose at a random var of two functions of the qubits below it: equal,
 * a multiple, opposite, one zero, the other zero, a translate with a sign
 * pattern, or unrelated.
 */
static void
xp_compose(unsigned n, tally_t *t, unsigned *cov)
{
    const uint64_t len = UINT64_C(1) << n;
    const size_t vb = len * sizeof(EVBDD_WGT);
    EVBDD_WGT *g0 = malloc(vb), *g1 = malloc(vb), *h = malloc(vb);
    const unsigned var = (unsigned)rnd_below(n);
    const uint64_t below = ((UINT64_C(1) << n) - 1) & ~((UINT64_C(2) << var) - 1);
    rand_on(n, rand_qubits(n) & below, 0.3 * (double)rnd_below(2), rnd_below(2), g0);
    rand_on(n, rand_qubits(n) & below, 0.3 * (double)rnd_below(2), rnd_below(2), g1);
    const EVBDD_WGT m = rand_value(0.0, true);
    const uint64_t im = imask(below, n);
    const uint64_t s = (bqd_family() == BQD_FAMILY_PAULI) ? rnd_below(len) & im : 0;
    const uint64_t tr = rnd_below(len) & im;
    const unsigned rel = (unsigned)rnd_below(7);
    for (uint64_t x = 0; x < len; x++) {
        switch (rel) {
        case 0: g1[x] = g0[x]; break;
        case 1: g1[x] = wgt_mul(m, g0[x]); break;
        case 2: g1[x] = wgt_neg(g0[x]); break;
        case 3: g1[x] = EVBDD_ZERO; break;
        case 4: g0[x] = EVBDD_ZERO; break;
        case 5: g1[x] = parity(s & x) ? wgt_neg(g0[x ^ tr]) : g0[x ^ tr]; break;
        default: break;
        }
    }
    for (uint64_t x = 0; x < len; x++) h[x] = (x & ibit(var, n)) ? g1[x] : g0[x];
    const BQD e0 = bqd_from_vector(g0, n), e1 = bqd_from_vector(g1, n);
    if ((!limdd_edge_is_zero(e0) && bqd_xp_top(e0) <= var)
        || (!limdd_edge_is_zero(e1) && bqd_xp_top(e1) <= var)) {
        /* a build that depends on a level it should not */
        t[X_COMPOSE].tot++;
        t[X_COMPOSE].bad++;
    } else {
        const BQD r = bqd_compose(var, e0, e1);
        compose_case(e0, e1, r, var, n, cov);
        tally(&t[X_COMPOSE], r, h, n);
    }
    free(g0); free(g1); free(h);
}

static void
check_xp_operations(void)
{
    const bool pauli = bqd_family() == BQD_FAMILY_PAULI;
    const EVBDD_WGT vals[4] = { EVBDD_ZERO, EVBDD_ONE, EVBDD_MIN_ONE, pw[2] };
    const size_t vb = MAXV * sizeof(EVBDD_WGT);
    EVBDD_WGT *f = malloc(vb), *g = malloc(vb), *h = malloc(vb), *g0 = malloc(vb), *g1 = malloc(vb);
    tally_t t[X_COUNT];
    unsigned cov[XC_COUNT];
    memset(t, 0, sizeof(t));
    memset(cov, 0, sizeof(cov));
    structure_begin();

    /* every vector over {0, 1, -1, i} on 1 to 3 qubits: the one-operand
     * operations on each; the pointwise ones on every pair on 1 and 2 qubits
     * and on a random partner on 3, and on labelled edges on the nodes of
     * every pair on 1 qubit, one pair in 16 on 2 and one in 4 on 3; and
     * Compose at var 0 of its two halves, which are then every pair of
     * functions on one qubit fewer */
    for (unsigned n = 1; n <= 3; n++) {
        const uint64_t len = UINT64_C(1) << n, half = len >> 1, count = UINT64_C(1) << (2 * len);
        for (uint64_t idx = 0; idx < count; idx++) {
            if (n == 3 && THIN > 1 && rnd_below(THIN) != 0) continue;
            for (uint64_t x = 0; x < len; x++) f[x] = vals[(idx >> (2 * x)) & 3];
            const BQD F = bqd_from_vector(f, n);
            xp_unary(F, f, n, n <= 2, 2, false, t, cov);
            for (uint64_t jdx = 0; jdx < count; jdx++) {
                const uint64_t j = (n <= 2) ? jdx : rnd_below(count);
                for (uint64_t x = 0; x < len; x++) g[x] = vals[(j >> (2 * x)) & 3];
                const bool lab = n == 1 || rnd_below(n == 2 ? 16 : 4) == 0;
                xp_binary(F, f, bqd_from_vector(g, n), g, n, lab, false, t, cov);
                if (n > 2) break;
            }
            for (uint64_t x = 0; x < len; x++) {
                g0[x] = f[x & (half - 1)];
                g1[x] = f[half | (x & (half - 1))];
            }
            const BQD e0 = bqd_from_vector(g0, n), e1 = bqd_from_vector(g1, n);
            const BQD r = bqd_compose(0, e0, e1);
            compose_case(e0, e1, r, 0, n, cov);
            tally(&t[X_COMPOSE], r, f, n);
        }
    }

    /* random states on 4 to 8 qubits, translated and signed at random */
    for (unsigned n = 4; n <= NQ; n++) for (int rep = 0; rep < 16; rep++) {
        const uint64_t len = UINT64_C(1) << n;
        rand_state(n, (unsigned)rnd_below(S_KINDS), f);
        if (rnd_below(2)) twist(n, f);
        rand_state(n, (unsigned)rnd_below(S_KINDS), g);
        if (rnd_below(2)) twist(n, g);
        const BQD F = bqd_from_vector(f, n), G = bqd_from_vector(g, n);
        check_diagram(F, f, n);
        check_diagram(G, g, n);
        xp_unary(F, f, n, false, 4, true, t, cov);
        xp_binary(F, f, G, g, n, true, true, t, cov);
        for (uint64_t x = 0; x < len; x++) g1[x] = rnd_below(3) ? wgt_neg(f[x]) : g[x];
        for (uint64_t x = 0; x < len; x++) h[x] = wgt_add(f[x], g1[x]);
        tally(&t[X_CANCEL], bqd_add(F, bqd_from_vector(g1, n)), h, n);
        for (int k = 0; k < 4; k++) xp_compose(n, t, cov);
    }

    for (int i = 0; i < X_COUNT; i++) report(xtally_name[i], t[i].bad, t[i].tot);
    for (int i = 0; i < XC_COUNT; i++) {
        if (!pauli && (i == XC_OPPOSITE || i == XC_REPAIR || i == XC_ZTOP)) continue;
        covered(xcover_name[i], cov[i]);
    }
    structure_report("operations");
    free(f); free(g); free(h); free(g0); free(g1);
}

/*
 * The operations of qsylvan_bqd_gates.h and qsylvan_bqd_ops.h in the
 * translation and Pauli families, where every one of them is the recursion
 * of qsylvan_bqd_xp.h on labelled edges: on every vector over {0, 1, -1, i}
 * on 1 and 2 qubits every gate, restriction, projection, local matvec,
 * controlled gate, swap, monomial and product there is, and on one vector in
 * eight on 3 qubits, at random, a random draw of them. check_operations and
 * check_diagonal run on larger random states in these families as well.
 */
enum {
    XE_GATE, XE_DYN, XE_RESTRICT, XE_PROJECT, XE_MATVEC1, XE_CGATE, XE_CPHASE, XE_EITHER,
    XE_SWAP, XE_MATVEC2, XE_DIAG, XE_PRODUCT, XE_FULL, XE_LGATE, XE_LRESTRICT, XE_LPROJECT,
    XE_COUNT
};
static const char *xe_name[XE_COUNT] = {
    "bqd_gate, table gates", "bqd_gate, random exact 2x2", "bqd_restrict", "bqd_project",
    "bqd_local_matvec on one qubit", "bqd_cgate, controls above the target",
    "bqd_cgate, phase with controls anywhere", "bqd_cgate_either, X and Z", "bqd_swap",
    "bqd_local_matvec on two qubits", "bqd_apply_diagonal", "bqd_product", "bqd_has_full_support",
    "Gate, labelled operands (bqd_xp_cgate_rec)", "Restrict, labelled operands (bqd_xp_restrict)",
    "Project, labelled operands (bqd_xp_project)",
};

enum {
    XEC_ROOT, XEC_TARGET, XEC_CONTROL, XEC_ZERO_LOW, XEC_CANCEL, XEC_BELOW, XEC_SKIP, XEC_COUNT
};
static const char *xec_name[XEC_COUNT] = {
    "states whose root edge skips", "gates on a skipped target",
    "controlled gates with a skipped control", "gates on a qubit whose low cofactor is zero",
    "gates that cancel an amplitude", "controls below the target (either, phase)",
    "results that skip a level",
};

/** Draw one of `count` things, or with `all` take them all: the range [*lo, *hi). */
static void
one_or_all(bool all, uint64_t count, uint64_t *lo, uint64_t *hi)
{
    *lo = all ? 0 : rnd_below(count);
    *hi = all ? count : *lo + 1;
}

/** tally, and count a result that skips a level. */
static void
xe_tally(tally_t *t, BQD got, const EVBDD_WGT *h, unsigned n, unsigned *cov)
{
    bool to_terminal;
    if (skipped_qubits(got, n, &to_terminal) != 0) cov[XEC_SKIP]++;
    tally(t, got, h, n);
}

/**
 * The operations on F, the diagram of f on n qubits, and when `all` every
 * gate on every qubit and every set of controls, and one random draw of each
 * kind otherwise. The dynamic gate is whatever set_dynamic put there.
 */
static void
xp_gates_on(BQD F, const EVBDD_WGT *f, unsigned n, bool all, tally_t *t, unsigned *cov)
{
    static const uint32_t table[] = { GATEID_H, GATEID_X, GATEID_Y, GATEID_Z, GATEID_S,
                                      GATEID_T, GATEID_Tdag, GATEID_sqrtX, GATEID_sqrtY,
                                      GATEID_proj0, GATEID_proj1 };
    static const uint32_t ctl[] = { GATEID_X, GATEID_Y, GATEID_Z, GATEID_H, GATEID_S,
                                    GATEID_T, GATEID_sqrtX, GATEID_dynamic };
    static const uint32_t ph[] = { GATEID_Z, GATEID_S, GATEID_T, GATEID_Sdag };
    const uint64_t len = UINT64_C(1) << n;
    EVBDD_WGT *h = malloc(len * sizeof(EVBDD_WGT));
    bool to_terminal;
    const uint64_t skF = skipped_qubits(F, n, &to_terminal);
    if (!limdd_edge_is_zero(F) && level_of(limdd_target(F), n) > 0) cov[XEC_ROOT]++;
    uint64_t lo, hi, glo, ghi;

    one_or_all(all, n, &lo, &hi);
    for (uint32_t q = (uint32_t)lo; q < hi; q++) {
        one_or_all(all, sizeof(table) / sizeof(table[0]), &glo, &ghi);
        for (uint64_t gi = glo; gi < ghi; gi++) {
            memcpy(h, f, len * sizeof(EVBDD_WGT));
            dense_gate(h, n, gates[table[gi]], 0, q);
            if ((skF >> q) & 1) cov[XEC_TARGET]++;
            gate_cases(f, h, n, q, &cov[XEC_ZERO_LOW], &cov[XEC_CANCEL]);
            xe_tally(&t[XE_GATE], bqd_gate(F, table[gi], q, n), h, n, cov);
        }
        memcpy(h, f, len * sizeof(EVBDD_WGT));
        dense_gate(h, n, gates[GATEID_dynamic], 0, q);
        gate_cases(f, h, n, q, &cov[XEC_ZERO_LOW], &cov[XEC_CANCEL]);
        xe_tally(&t[XE_DYN], bqd_gate(F, GATEID_dynamic, q, n), h, n, cov);
        xe_tally(&t[XE_MATVEC1], bqd_local_matvec(F, gates[GATEID_dynamic], &q, 1, n), h, n, cov);
        for (int b = 0; b < 2; b++) {
            dense_restrict(f, n, q, b, h);
            xe_tally(&t[XE_RESTRICT], bqd_restrict(F, q, b), h, n, cov);
            dense_project(f, n, q, b, h);
            xe_tally(&t[XE_PROJECT], bqd_project(F, q, b), h, n, cov);
        }
    }

    if (n >= 2) {
        /* controls above the target: every nonempty set of them, or one */
        one_or_all(all, n - 1, &lo, &hi);
        for (uint32_t q = (uint32_t)lo + 1; q < hi + 1; q++) {
            uint64_t clo, chi;
            one_or_all(all, (UINT64_C(1) << q) - 1, &clo, &chi);
            for (uint64_t cm = clo + 1; cm < chi + 1; cm++) {
                one_or_all(all, sizeof(ctl) / sizeof(ctl[0]), &glo, &ghi);
                for (uint64_t gi = glo; gi < ghi; gi++) {
                    if (cm & skF) cov[XEC_CONTROL]++;
                    if ((skF >> q) & 1) cov[XEC_TARGET]++;
                    memcpy(h, f, len * sizeof(EVBDD_WGT));
                    dense_gate(h, n, gates[ctl[gi]], cm, q);
                    xe_tally(&t[XE_CGATE], bqd_cgate(F, ctl[gi], cm, q, n), h, n, cov);
                }
            }
        }
        /* a phase, its controls anywhere but the target: every set, or one */
        one_or_all(all, n, &lo, &hi);
        for (uint32_t q = (uint32_t)lo; q < hi; q++) {
            const uint64_t others = ((UINT64_C(1) << n) - 1) & ~(UINT64_C(1) << q);
            uint64_t mlo, mhi;
            one_or_all(all, UINT64_C(1) << n, &mlo, &mhi);
            for (uint64_t m = mlo; m < mhi; m++) {
                const uint64_t pm = all ? m & others : others & (m | (UINT64_C(1) << rnd_below(n)));
                if (pm == 0 || (all && (m & ~others))) continue;
                if (pm >> (q + 1)) cov[XEC_BELOW]++;
                one_or_all(all, 4, &glo, &ghi);
                for (uint64_t gi = glo; gi < ghi; gi++) {
                    memcpy(h, f, len * sizeof(EVBDD_WGT));
                    dense_gate(h, n, gates[ph[gi]], pm, q);
                    xe_tally(&t[XE_CPHASE], bqd_cgate(F, ph[gi], pm, q, n), h, n, cov);
                }
            }
        }
        /* one control on either side, a swap, and a matvec, on every pair or one */
        one_or_all(all, (uint64_t)n * n, &lo, &hi);
        for (uint64_t ab = lo; ab < hi; ab++) {
            const unsigned a = (unsigned)(ab / n), b = (unsigned)(ab % n);
            if (a == b) continue;
            if (a > b) cov[XEC_BELOW]++;
            for (int k = 0; k < 2; k++) {
                const uint32_t eg = k ? GATEID_Z : GATEID_X;
                memcpy(h, f, len * sizeof(EVBDD_WGT));
                dense_gate(h, n, gates[eg], UINT64_C(1) << a, b);
                bool ok;
                const BQD r = bqd_cgate_either(F, eg, a, b, n, &ok);
                xe_tally(&t[XE_EITHER], r, h, n, cov);
                if (!ok) t[XE_EITHER].bad++;
            }
            memcpy(h, f, len * sizeof(EVBDD_WGT));
            dense_swap(h, n, a, b);
            xe_tally(&t[XE_SWAP], bqd_swap(F, a, b, n), h, n, cov);
            const uint32_t qs[2] = { a, b };
            EVBDD_WGT M[16];
            for (int i = 0; i < 16; i++) M[i] = rand_entry();
            dense_matvec(f, n, M, qs, 2, h);
            xe_tally(&t[XE_MATVEC2], bqd_local_matvec(F, M, qs, 2, n), h, n, cov);
        }
    }

    /* the monomials on every set of qubits with every phase but 1, or one */
    one_or_all(all, ((UINT64_C(1) << n) - 1) * 7, &lo, &hi);
    for (uint64_t k = lo; k < hi; k++) {
        const uint64_t A = imask(1 + k / 7, n);
        const EVBDD_WGT phase = pw[1 + k % 7];
        for (uint64_t x = 0; x < len; x++) h[x] = ((x & A) == A) ? wgt_mul(f[x], phase) : f[x];
        xe_tally(&t[XE_DIAG], bqd_apply_diagonal(F, A, phase, n), h, n, cov);
    }

    bool full = true;
    for (uint64_t x = 0; x < len; x++) if (f[x] == EVBDD_ZERO) full = false;
    t[XE_FULL].tot++;
    if (bqd_has_full_support(F) != full) t[XE_FULL].bad++;
    free(h);
}

static void
check_xp_exhaustive(void)
{
    const EVBDD_WGT vals[4] = { EVBDD_ZERO, EVBDD_ONE, EVBDD_MIN_ONE, pw[2] };
    EVBDD_WGT f[8], g[8], h[8];
    tally_t t[XE_COUNT];
    unsigned cov[XEC_COUNT];
    memset(t, 0, sizeof(t));
    memset(cov, 0, sizeof(cov));
    structure_begin();

    for (unsigned n = 1; n <= 3; n++) {
        const uint64_t len = UINT64_C(1) << n, count = UINT64_C(1) << (2 * len);
        EVBDD_WGT U[4];
        for (int i = 0; i < 4; i++) U[i] = rand_entry();
        set_dynamic(U);
        for (uint64_t idx = 0; idx < count; idx++) {
            if (n == 3 && rnd_below(8 * THIN) != 0) continue;     /* one vector in 8 on 3 qubits */
            for (uint64_t x = 0; x < len; x++) f[x] = vals[(idx >> (2 * x)) & 3];
            const BQD F = bqd_from_vector(f, n);
            xp_gates_on(F, f, n, n <= 2, t, cov);
            xp_labelled_gates(F, n, 2, &t[XE_LGATE]);
            /* the product with every partner on 1 and 2 qubits, a random one on 3 */
            for (uint64_t jdx = 0; jdx < count; jdx++) {
                const uint64_t j = (n <= 2) ? jdx : rnd_below(count);
                for (uint64_t x = 0; x < len; x++) g[x] = vals[(j >> (2 * x)) & 3];
                for (uint64_t x = 0; x < len; x++) h[x] = wgt_mul(f[x], g[x]);
                xe_tally(&t[XE_PRODUCT], bqd_product(F, bqd_from_vector(g, n)), h, n, cov);
                if (n > 2) break;
            }
        }
    }

    for (int i = 0; i < XE_COUNT; i++) report(xe_name[i], t[i].bad, t[i].tot);
    for (int i = 0; i < XEC_COUNT; i++) covered(xec_name[i], cov[i]);
    structure_report("n <= 3, gates and the rest");
}

/* --- collections, scalar family ------------------------------------------ */

/*
 * A collection sweeps what the protected edges do not reach, clears the
 * operation cache and hands the swept buckets out again. An unmarked node
 * that an edge reaches across skipped levels, a memo that outlives a sweep,
 * or a node built again in a recycled bucket then shows as an edge that is no
 * longer the one bqd_from_vector builds. The tables are small, so that a
 * collection is cheap and what it frees is soon built on again, and GC_WIDTH
 * is wider than any state here, so no diagram ends at the table's terminal
 * level. The weight table is not collected here: it is far from half full,
 * and a collection of it would move the weights the dense vectors hold.
 */
#define GC_WIDTH (NQ + 2)
#define POOL     5

static unsigned collections;

/**
 * Gate sequences as the qasm runner applies them, collecting between gates:
 * from |0...0> on 2 to NQ qubits, the runner's gates with a control on either
 * side, swaps, the diagonal walk, and a random exact 2x2, singular as often
 * as not, under a control above its target. After every gate the state must
 * be the edge bqd_from_vector builds for the dense state.
 */
VOID_TASK_0(check_gc_circuits)
{
    static const uint32_t one[] = { GATEID_H, GATEID_H, GATEID_X, GATEID_Y, GATEID_Z,
                                    GATEID_S, GATEID_T, GATEID_Tdag, GATEID_sqrtX };
    static const uint32_t two[] = { GATEID_X, GATEID_Z, GATEID_S, GATEID_T };
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT));
    tally_t t = { 0, 0 };
    BQD state = limdd_zero_edge();
    limdd_protect(&state);

    for (int run = 0; run < 40 / THIN; run++) {
        const unsigned n = 2 + (unsigned)rnd_below(NQ - 1);
        const uint64_t len = UINT64_C(1) << n;
        EVBDD_WGT U[4];
        for (int i = 0; i < 4; i++) U[i] = rand_entry();
        set_dynamic(U);
        /* from |0...0>, and in the translation and Pauli families from any
         * basis state, which is a translate of it */
        const uint64_t x0 = bqd_family() != BQD_FAMILY_SCALAR ? rnd_below(len) : 0;
        for (uint64_t x = 0; x < len; x++) f[x] = (x == x0) ? EVBDD_ONE : EVBDD_ZERO;
        state = bqd_basis_state(x0, n);
        for (int k = 0; k < 40; k++) {
            const unsigned q = (unsigned)rnd_below(n);
            unsigned c = (unsigned)rnd_below(n - 1); if (c >= q) c++;
            switch (rnd_below(6)) {
            case 0: case 1: {
                const uint32_t gid = one[rnd_below(sizeof(one) / sizeof(one[0]))];
                dense_gate(f, n, gates[gid], 0, q);
                state = bqd_gate(state, gid, q, n);
                break;
            }
            case 2: {
                const uint32_t gid = two[rnd_below(4)];
                bool ok;
                dense_gate(f, n, gates[gid], UINT64_C(1) << c, q);
                state = bqd_cgate_either(state, gid, c, q, n, &ok);
                if (!ok) t.bad++;
                break;
            }
            case 3:
                dense_swap(f, n, c, q);
                state = bqd_swap(state, c, q, n);
                break;
            case 4: {
                const unsigned ctl = c < q ? c : q, tgt = c < q ? q : c;
                dense_gate(f, n, U, UINT64_C(1) << ctl, tgt);
                state = bqd_cgate(state, GATEID_dynamic, UINT64_C(1) << ctl, tgt, n);
                break;
            }
            default: {
                const uint64_t A = imask(1 + rnd_below(len - 1), n);
                const EVBDD_WGT phase = pw[1 + rnd_below(7)];
                for (uint64_t x = 0; x < len; x++) if ((x & A) == A) f[x] = wgt_mul(f[x], phase);
                state = bqd_apply_diagonal(state, A, phase, n);
                break;
            }
            }
            t.tot++;
            if (state != bqd_from_vector(f, n)) t.bad++;
            if (limdd_edge_is_zero(state)) break;      /* a singular gate ended the run */
            if (rnd_below(15) == 0) { CALL(limdd_gc); collections++; }
        }
    }

    limdd_unprotect(&state);
    report("gate sequences, a collection every ~15 gates", t.bad, t.tot);
    free(f);
}

/** Kept results on n qubits: protected edges beside their dense vectors. */
typedef struct {
    unsigned  n;
    BQD       e[POOL];
    EVBDD_WGT v[POOL][MAXV];
} pool_t;

enum {
    G_MUL, G_ADD, G_GATE, G_CGATE, G_CPHASE, G_EITHER, G_RESTRICT, G_PROJECT, G_PRODUCT,
    G_PRODFULL, G_DIAG, G_SWAP, G_MATVEC, G_KINDS
};

/** An operation, drawn once so that it can be applied at two widths. */
typedef struct {
    unsigned  kind, q, c;
    int       b;
    uint64_t  cm;           /* the controls, or the qubits of the monomial */
    uint32_t  gid;
    EVBDD_WGT phase, M[16];
} gc_op_t;

static void
draw_op(gc_op_t *o, unsigned n)
{
    static const uint32_t table[] = { GATEID_H, GATEID_X, GATEID_Y, GATEID_Z, GATEID_S,
                                      GATEID_T, GATEID_sqrtX, GATEID_proj1 };
    static const uint32_t ctl[] = { GATEID_X, GATEID_Y, GATEID_H, GATEID_S, GATEID_sqrtX };
    static const uint32_t ph[] = { GATEID_Z, GATEID_S, GATEID_T, GATEID_Sdag };
    o->kind = (unsigned)rnd_below(G_KINDS);
    o->q = (unsigned)rnd_below(n);
    o->c = (unsigned)rnd_below(n - 1); if (o->c >= o->q) o->c++;
    o->b = (int)rnd_below(2);
    o->gid = table[rnd_below(sizeof(table) / sizeof(table[0]))];
    o->cm = 1 + rnd_below((UINT64_C(1) << n) - 1);
    o->phase = pw[1 + rnd_below(7)];
    for (int i = 0; i < 16; i++) o->M[i] = rand_entry();
    switch (o->kind) {
    case G_GATE:
        o->cm = 0;
        break;
    case G_CGATE:           /* controls above the target */
        o->q = 1 + (unsigned)rnd_below(n - 1);
        o->cm = 1 + rnd_below((UINT64_C(1) << o->q) - 1);
        o->gid = ctl[rnd_below(sizeof(ctl) / sizeof(ctl[0]))];
        break;
    case G_CPHASE:          /* controls anywhere but the target */
        o->gid = ph[rnd_below(4)];
        while ((o->cm &= ~(UINT64_C(1) << o->q)) == 0) o->cm = rnd_below(UINT64_C(1) << n);
        break;
    case G_EITHER:          /* one control, on either side */
        o->gid = rnd_below(2) ? GATEID_X : GATEID_Z;
        o->cm = UINT64_C(1) << o->c;
        break;
    default:
        break;
    }
}

static void
dense_op(const gc_op_t *o, const EVBDD_WGT *f, const EVBDD_WGT *g, unsigned n, EVBDD_WGT *h)
{
    const uint64_t len = UINT64_C(1) << n;
    const uint32_t qs[2] = { o->q, o->c };
    switch (o->kind) {
    case G_MUL: case G_PRODUCT: case G_PRODFULL:
        for (uint64_t x = 0; x < len; x++) h[x] = wgt_mul(f[x], g[x]);
        break;
    case G_ADD:
        for (uint64_t x = 0; x < len; x++) h[x] = wgt_add(f[x], g[x]);
        break;
    case G_GATE: case G_CGATE: case G_CPHASE: case G_EITHER:
        memcpy(h, f, len * sizeof(EVBDD_WGT));
        dense_gate(h, n, gates[o->gid], o->cm, o->q);
        break;
    case G_RESTRICT:
        dense_restrict(f, n, o->q, o->b, h);
        break;
    case G_PROJECT:
        dense_project(f, n, o->q, o->b, h);
        break;
    case G_DIAG: {
        const uint64_t A = imask(o->cm, n);
        for (uint64_t x = 0; x < len; x++) h[x] = ((x & A) == A) ? wgt_mul(f[x], o->phase) : f[x];
        break;
    }
    case G_SWAP:
        memcpy(h, f, len * sizeof(EVBDD_WGT));
        dense_swap(h, n, o->q, o->c);
        break;
    default:
        dense_matvec(f, n, o->M, qs, 2, h);
        break;
    }
}

/** The operation on F and G read on n qubits. */
static BQD
apply_op(const gc_op_t *o, BQD F, BQD G, unsigned n)
{
    const uint32_t qs[2] = { o->q, o->c };
    switch (o->kind) {
    case G_MUL:      return bqd_multiply(F, G);
    case G_ADD:      return bqd_add(F, G);
    case G_GATE:     return bqd_gate(F, o->gid, o->q, n);
    case G_CGATE: case G_CPHASE:
                     return bqd_cgate(F, o->gid, o->cm, o->q, n);
    case G_EITHER: {
        /* X and Z are taken on either side, so a refusal must fail */
        bool ok;
        const BQD r = bqd_cgate_either(F, o->gid, o->c, o->q, n, &ok);
        return ok ? r : limdd_zero_edge();
    }
    case G_RESTRICT: return bqd_restrict(F, o->q, o->b);
    case G_PROJECT:  return bqd_project(F, o->q, o->b);
    case G_PRODUCT: case G_PRODFULL:
                     return bqd_product(F, G);
    case G_DIAG:     return bqd_apply_diagonal(F, imask(o->cm, n), o->phase, n);
    case G_SWAP:     return bqd_swap(F, o->q, o->c, n);
    default:         return bqd_local_matvec(F, o->M, qs, 2, n);
    }
}

static void
pool_fresh(pool_t *p, unsigned i)
{
    rand_state(p->n, (unsigned)rnd_below(S_KINDS), p->v[i]);
    if (bqd_family() != BQD_FAMILY_SCALAR && rnd_below(2)) twist(p->n, p->v[i]);
    p->e[i] = bqd_from_vector(p->v[i], p->n);
}

/**
 * Results kept across collections, in two pools that share the tables, on NQ
 * and on NQ - 2 qubits. A random operation takes its operands from one pool
 * and puts its result back, and every one to four operations a collection
 * runs. After it every kept edge must still decode to its vector and still
 * be the edge bqd_from_vector builds, and the set of checked nodes starts
 * empty, since the indices it holds may have been handed out again: every
 * node a result reaches after a collection passes the structure checks anew.
 */
VOID_TASK_0(check_gc_pool)
{
    pool_t *pools = malloc(2 * sizeof(pool_t));
    EVBDD_WGT *h = malloc(MAXV * sizeof(EVBDD_WGT)), *g = malloc(MAXV * sizeof(EVBDD_WGT));
    tally_t t = { 0, 0 }, wide = { 0, 0 }, kept = { 0, 0 };
    unsigned skipping = 0, freed = 0;
    pools[0].n = NQ;
    pools[1].n = NQ - 2;
    for (int k = 0; k < 2; k++) for (unsigned i = 0; i < POOL; i++) {
        pool_fresh(&pools[k], i);
        limdd_protect(&pools[k].e[i]);
    }

    unsigned next = 1 + (unsigned)rnd_below(4);
    for (int it = 0; it < 250 / THIN; it++) {
        pool_t *p = &pools[rnd_below(3) == 0];
        const unsigned n = p->n;
        const unsigned a = (unsigned)rnd_below(POOL), b = (unsigned)rnd_below(POOL);
        gc_op_t o;
        draw_op(&o, n);
        if (o.kind == G_PRODFULL) {
            /* two phase states, so that the product contracts on skip:cor:full */
            phase_on(n, rand_qubits(n), p->v[a]);
            p->e[a] = bqd_from_vector(p->v[a], n);
            phase_on(n, rand_qubits(n), p->v[b]);
            p->e[b] = bqd_from_vector(p->v[b], n);
        }
        dense_op(&o, p->v[a], p->v[b], n, h);
        const BQD r = apply_op(&o, p->e[a], p->e[b], n);
        tally(&t, r, h, n);
        bool to_terminal;
        if (skipped_qubits(r, n, &to_terminal) != 0) skipping++;
        if (n < NQ) {
            /* the same edges read on NQ qubits, the last two of which they ignore */
            wide.tot++;
            if (apply_op(&o, p->e[a], p->e[b], NQ) != r) wide.bad++;
        }
        const unsigned d = (unsigned)rnd_below(POOL);
        p->e[d] = r;
        memcpy(p->v[d], h, (UINT64_C(1) << n) * sizeof(EVBDD_WGT));
        /* products and projections thin a pool out, so now and then a fresh state */
        if (rnd_below(8) == 0) pool_fresh(p, (unsigned)rnd_below(POOL));

        if (--next > 0) continue;
        next = 1 + (unsigned)rnd_below(4);
        const size_t before = limdd_node_table_count();
        CALL(limdd_gc);
        collections++;
        if (limdd_node_table_count() < before) freed++;
        set_clear(checked);
        for (int k = 0; k < 2; k++) for (unsigned i = 0; i < POOL; i++) {
            const pool_t *kp = &pools[k];
            kept.tot++;
            bqd_to_vector(kp->e[i], kp->n, g);
            if (!same_vector(g, kp->v[i], kp->n) || kp->e[i] != bqd_from_vector(kp->v[i], kp->n))
                kept.bad++;
        }
    }

    for (int k = 0; k < 2; k++) for (unsigned i = 0; i < POOL; i++) limdd_unprotect(&pools[k].e[i]);
    report("operations on kept results, collected every 1 to 4", t.bad, t.tot);
    report("the same operation on the edges read on NQ qubits", wide.bad, wide.tot);
    report("kept results after a collection, decoded and rebuilt", kept.bad, kept.tot);
    covered("results that skip a level", skipping);
    covered("collections that freed nodes", freed);
    free(pools); free(h); free(g);
}

TASK_0(int, run_collections)
{
    init_copy(BQD_FAMILY_SCALAR, GC_WIDTH, 1LL << 16, 1LL << 16, 1LL << 17, 1LL << 16);
    exact_constants();
    rng_state = UINT64_C(0x5EED6C000000);

    const int before = failures;
    collections = 0;
    structure_begin();
    CALL(check_gc_circuits);
    CALL(check_gc_pool);
    covered("collections", collections);
    structure_report("collections");
    bqd_quit();
    return failures != before;
}

/* --- collections, translation and Pauli families ---------------------------- */

enum {
    XG_MUL, XG_ADD, XG_XPROD, XG_XQUOT, XG_SCALE, XG_NEGATE, XG_COF, XG_FLIP, XG_SIGN,
    XG_CANON, XG_KINDS
};

/**
 * A fresh state of a kind whose values are in Z[w_8, 1/2], a phase, circuit
 * or basis state or one with values w_8^k and zeros, translated and signed
 * at random in the families that have those.
 */
static void
pool_fresh_xp(pool_t *p, unsigned i)
{
    static const unsigned kinds[] = { S_PHASE, S_CIRCUIT, S_BASIS_PLUS };
    if (rnd_below(4) == 0) rand_on(p->n, rand_qubits(p->n), 0.3, false, p->v[i]);
    else rand_state(p->n, kinds[rnd_below(3)], p->v[i]);
    twist(p->n, p->v[i]);
    p->e[i] = bqd_from_vector(p->v[i], p->n);
}

/**
 * As check_gc_pool, with the operations of the two families: the four
 * pointwise ones on kept results or on random labelled edges on their
 * nodes, scale and negate, a cofactor, an X and a Z at a level down to the
 * top made from the two cofactors by Compose, and Canon of a labelled edge.
 * Each result must be the edge bqd_from_vector builds, and every kept one
 * must be it still after each collection. The states are those of
 * pool_fresh_xp, scalars are powers of w_8, and a quotient is not kept:
 * kept products and quotients of algebraic numbers grow without bound in
 * exact arithmetic, and the test would time the rationals.
 */
VOID_TASK_0(check_gc_xp_pool)
{
    pool_t *pools = malloc(2 * sizeof(pool_t));
    EVBDD_WGT *h = malloc(MAXV * sizeof(EVBDD_WGT)), *g = malloc(MAXV * sizeof(EVBDD_WGT));
    EVBDD_WGT *la = malloc(MAXV * sizeof(EVBDD_WGT)), *lb = malloc(MAXV * sizeof(EVBDD_WGT));
    tally_t t = { 0, 0 }, kept = { 0, 0 };
    unsigned labelled = 0, freed = 0;
    pools[0].n = NQ;
    pools[1].n = NQ - 2;
    for (int k = 0; k < 2; k++) for (unsigned i = 0; i < POOL; i++) {
        pool_fresh_xp(&pools[k], i);
        limdd_protect(&pools[k].e[i]);
    }

    unsigned next = 1 + (unsigned)rnd_below(4);
    for (int it = 0; it < 400 / THIN; it++) {
        pool_t *p = &pools[rnd_below(3) == 0];
        const unsigned n = p->n;
        const uint64_t len = UINT64_C(1) << n;
        const unsigned a = (unsigned)rnd_below(POOL), b = (unsigned)rnd_below(POOL);
        BQD A = p->e[a], B = p->e[b];
        const unsigned top = bqd_xp_top(A);
        const unsigned v = (unsigned)rnd_below((top < n ? top : n - 1) + 1);
        bool keep = true;
        BQD r;
        switch (rnd_below(XG_KINDS)) {
        case XG_MUL: case XG_ADD: case XG_XPROD: case XG_XQUOT: {
            const int op = (int)rnd_below(4);
            keep = op != BQD_OP_XQUOT;
            memcpy(la, p->v[a], len * sizeof(EVBDD_WGT));
            memcpy(lb, p->v[b], len * sizeof(EVBDD_WGT));
            if (rnd_below(2)) {
                A = relabel(A, rand_label(n, false), n, la);
                B = relabel(B, rand_label(n, false), n, lb);
                labelled++;
            }
            dense_pointwise(op, la, lb, n, h);
            r = pointwise(op, A, B);
            break;
        }
        case XG_SCALE: {
            const EVBDD_WGT c = pw[rnd_below(8)];
            for (uint64_t x = 0; x < len; x++) h[x] = wgt_mul(c, p->v[a][x]);
            r = bqd_scale(A, c);
            break;
        }
        case XG_NEGATE:
            for (uint64_t x = 0; x < len; x++) h[x] = wgt_neg(p->v[a][x]);
            r = bqd_negate(A);
            break;
        case XG_COF: {
            const int bb = (int)rnd_below(2);
            dense_restrict(p->v[a], n, v, bb, h);
            r = bqd_cofactor(A, v, bb);
            break;
        }
        case XG_FLIP: case XG_SIGN: {
            const bool flip = rnd_below(2);
            memcpy(h, p->v[a], len * sizeof(EVBDD_WGT));
            dense_gate(h, n, gates[flip ? GATEID_X : GATEID_Z], 0, v);
            const BQD c0 = bqd_cofactor(A, v, 0), c1 = bqd_cofactor(A, v, 1);
            r = flip ? bqd_compose(v, c1, c0) : bqd_compose(v, c0, bqd_negate(c1));
            break;
        }
        default:
            r = bqd_xp_canon(relabel(A, rand_label(n, false), n, h));
            break;
        }
        tally(&t, r, h, n);
        if (keep) {
            const unsigned d = (unsigned)rnd_below(POOL);
            p->e[d] = r;
            memcpy(p->v[d], h, len * sizeof(EVBDD_WGT));
        }
        /* products and cofactors thin a pool out, so now and then a fresh state */
        if (rnd_below(8) == 0) pool_fresh_xp(p, (unsigned)rnd_below(POOL));

        if (--next > 0) continue;
        next = 1 + (unsigned)rnd_below(4);
        const size_t before = limdd_node_table_count();
        CALL(limdd_gc);
        collections++;
        if (limdd_node_table_count() < before) freed++;
        set_clear(checked);
        for (int k = 0; k < 2; k++) for (unsigned i = 0; i < POOL; i++) {
            const pool_t *kp = &pools[k];
            kept.tot++;
            bqd_to_vector(kp->e[i], kp->n, g);
            if (!same_vector(g, kp->v[i], kp->n) || kp->e[i] != bqd_from_vector(kp->v[i], kp->n))
                kept.bad++;
        }
    }

    for (int k = 0; k < 2; k++) for (unsigned i = 0; i < POOL; i++) limdd_unprotect(&pools[k].e[i]);
    report("operations on kept results, collected every 1 to 4", t.bad, t.tot);
    report("kept results after a collection, decoded and rebuilt", kept.bad, kept.tot);
    covered("operations on random labelled edges", labelled);
    covered("collections that freed nodes", freed);
    free(pools); free(h); free(g); free(la); free(lb);
}

TASK_1(int, run_collections_xp, int, fam)
{
    init_copy((bqd_family_t)fam, GC_WIDTH, 1LL << 16, 1LL << 16, 1LL << 17, 1LL << 16);
    exact_constants();
    rng_state = UINT64_C(0x5EED6C000000) + (uint64_t)fam;

    const int before = failures;
    collections = 0;
    structure_begin();
    CALL(check_gc_xp_pool);
    CALL(check_gc_circuits);
    CALL(check_gc_pool);
    covered("collections", collections);
    structure_report("collections");
    bqd_quit();
    return failures != before;
}

/* --- sessions and contention --------------------------------------------- */

/*
 * Two BQD sessions, bqd_init to bqd_quit, in one Sylvan package, each running
 * the same gates from |0...0> on 3 qubits and checked after every gate. The
 * second session's tables hand out the first one's indices again, so a memo
 * entry that outlived bqd_quit would answer for a different state; and the
 * translation and Pauli families share their cache ids, so an entry of one
 * would answer the other, whose Compose differs.
 */
static unsigned
session_gates(bqd_family_t fam)
{
    enum { SN = 3, SLEN = 1 << SN };
    static const uint32_t one[] = { GATEID_H, GATEID_X, GATEID_S, GATEID_T, GATEID_Y, GATEID_Z };
    init_copy(fam, SN, 1LL << 16, 1LL << 16, 1LL << 18, 1LL << 16);
    exact_constants();
    rng_state = UINT64_C(0x5E55105);
    EVBDD_WGT v[SLEN], w[SLEN];
    for (unsigned x = 0; x < SLEN; x++) v[x] = x ? EVBDD_ZERO : EVBDD_ONE;
    BQD s = bqd_basis_state(0, SN);
    unsigned bad = 0;
    for (int i = 0; i < 60; i++) {
        const unsigned a = (unsigned)rnd_below(SN);
        if (rnd_below(3) == 0 && a + 1 < SN) {
            s = bqd_cgate(s, GATEID_X, UINT64_C(1) << a, a + 1, SN);
            dense_gate(v, SN, gates[GATEID_X], UINT64_C(1) << a, a + 1);
        } else {
            const uint32_t g = one[rnd_below(6)];
            s = bqd_gate(s, g, a, SN);
            dense_gate(v, SN, gates[g], 0, a);
        }
        bqd_to_vector(s, SN, w);
        if (!same_vector(v, w, SN) || s != bqd_from_vector(v, SN)) bad++;
    }
    bqd_quit();
    return bad;
}

TASK_0(int, run_sessions)
{
    static const int pairs[][2] = {
        { BQD_FAMILY_X, BQD_FAMILY_PAULI }, { BQD_FAMILY_PAULI, BQD_FAMILY_X },
        { BQD_FAMILY_X, BQD_FAMILY_X }, { BQD_FAMILY_PAULI, BQD_FAMILY_PAULI },
        { BQD_FAMILY_SCALAR, BQD_FAMILY_SCALAR }, { BQD_FAMILY_SCALAR, BQD_FAMILY_PAULI },
    };
    const int before = failures;
    for (size_t k = 0; k < sizeof(pairs) / sizeof(pairs[0]); k++) {
        const unsigned first = session_gates((bqd_family_t)pairs[k][0]);
        const unsigned second = session_gates((bqd_family_t)pairs[k][1]);
        char what[96];
        snprintf(what, sizeof(what), "%s, then %s: gates wrong",
                 bqd_family_name((bqd_family_t)pairs[k][0]),
                 bqd_family_name((bqd_family_t)pairs[k][1]));
        report(what, first + second, 120);
    }
    return failures != before;
}

/*
 * Every worker asks for collections at once, each a different number of
 * times, so that some ask again right after a frame another one opened, and
 * go on asking after the others are done. A loser once waited for the
 * winner's frame rather than for the winner: one that asked after that frame
 * had ended, but before the winner had lowered its flag, waited for a frame
 * that no one would open. The alarm turns such a hang into a failure.
 */
VOID_TASK_0(ask_for_collections)
{
    const unsigned rounds = 4 + 4 * (unsigned)LACE_WORKER_ID;
    for (unsigned i = 0; i < rounds; i++) CALL(limdd_gc);
}

TASK_0(int, run_contention)
{
    enum { CN = 6, CLEN = 1 << CN };
    init_copy(BQD_FAMILY_PAULI, GC_WIDTH, 1LL << 16, 1LL << 16, 1LL << 17, 1LL << 16);
    exact_constants();
    rng_state = UINT64_C(0xC0117E57);
    EVBDD_WGT f[CLEN], g[CLEN];
    rand_state(CN, S_CIRCUIT, f);
    BQD psi = bqd_from_vector(f, CN);
    limdd_protect(&psi);

    const int before = failures;
    const unsigned rounds = 40;
    alarm(300);
    for (unsigned r = 0; r < rounds; r++) TOGETHER(ask_for_collections);
    alarm(0);
    bqd_to_vector(psi, CN, g);
    report("a state kept across collections every worker asked for",
           !same_vector(f, g, CN) || psi != bqd_from_vector(f, CN), 1);
    covered("rounds of collections asked for by every worker at once", rounds);

    limdd_unprotect(&psi);
    bqd_quit();
    return failures != before;
}

/*
 * The Pauli table among the tables whose filling asks for a collection. The
 * Pauli family makes a new Pauli word with many of its label products, far
 * more of them than new nodes, so in a Pauli table no larger than the node
 * table and a sixty-fourth of the LIM table only the Pauli table fills.
 * Clifford gates on 16 qubits as the runner applies them, asking
 * limdd_gc_wanted after each: without the Pauli table in the trigger the run
 * stops with that table full. The state after the last gate must be the one
 * the same gates make in large tables, which never collect.
 */
TASK_2(BQD, trigger_gates, bool, collect, unsigned*, collected)
{
    enum { TN = 16 };
    static const uint32_t one[] = { GATEID_H, GATEID_H, GATEID_X, GATEID_Y, GATEID_Z,
                                    GATEID_S, GATEID_Sdag };
    rng_state = UINT64_C(0x7E57);
    BQD s = bqd_basis_state(0, TN);
    limdd_protect(&s);
    for (int i = 0; i < 1200; i++) {
        if (rnd() % 2) {
            s = bqd_gate(s, one[rnd() % 7], (uint32_t)(rnd() % TN), TN);
        } else {
            const unsigned a = (unsigned)(rnd() % TN), b0 = (unsigned)(rnd() % TN);
            const unsigned b = a == b0 ? (a + 1) % TN : b0;
            if (rnd() % 2) {
                s = bqd_cgate(s, GATEID_Z, UINT64_C(1) << a, b, TN);
            } else {
                const unsigned c = a < b ? a : b, t = a < b ? b : a;
                s = bqd_cgate(s, GATEID_X, UINT64_C(1) << c, t, TN);
            }
        }
        if (collect && limdd_gc_wanted()) { CALL(limdd_gc); (*collected)++; }
    }
    limdd_unprotect(&s);
    return s;
}

TASK_0(int, run_pauli_trigger)
{
    enum { TN = 16, TLEN = 1 << TN };
    EVBDD_WGT *a = malloc(TLEN * sizeof(EVBDD_WGT)), *b = malloc(TLEN * sizeof(EVBDD_WGT));
    unsigned collected = 0;
    const int before = failures;

    init_copy(BQD_FAMILY_PAULI, TN, 1LL << 14, 1LL << 14, 1LL << 20, 1LL << 14);
    bqd_to_vector(CALL(trigger_gates, true, &collected), TN, a);
    bqd_quit();
    init_copy(BQD_FAMILY_PAULI, TN, 1LL << 20, 1LL << 22, 1LL << 22, 1LL << 20);
    bqd_to_vector(CALL(trigger_gates, false, &collected), TN, b);
    bqd_quit();

    report("Pauli table in the trigger: the state in small tables", !same_vector(a, b, TN), 1);
    covered("collections asked for between the gates", collected);
    free(a); free(b);
    return failures != before;
}

/* --- harness -------------------------------------------------------------- */

TASK_1(int, run_family, int, fam)
{
    init_copy((bqd_family_t)fam, NQ, 1LL << 21, 1LL << 20, 1LL << 22, 1LL << 20);
    exact_constants();
    rng_state = UINT64_C(0x5EED5C1B0000) + (uint64_t)fam;

    const int before = failures;
    check_exhaustive();
    check_random_builds();
    if (fam != BQD_FAMILY_SCALAR) {
        check_xp_operations();
        check_xp_exhaustive();
    }
    check_operations();
    check_diagonal();
    printf("  %zu nodes in the table\n", limdd_node_table_count());
    bqd_quit();
    return failures != before;
}

/** A fresh session, which every family and the collections run in, with 2^cache memo entries. */
static void
session_begin(unsigned workers, int cache)
{
    lace_start(workers, 0);
    sylvan_set_sizes(1LL << 20, 1LL << 20, 1LL << cache, 1LL << cache);
    sylvan_init_package();
    qsylvan_init_simulator(1LL << 22, 1LL << 22, -1, QISQ2_MAP, NORM_LOW);
}

static void
session_end(void)
{
    sylvan_quit();
    lace_stop();
}

int
main(void)
{
    /* a line at a time, so that a crash leaves the check it happened in */
    setvbuf(stdout, NULL, _IOLBF, 0);
    const unsigned workers = getenv("BQD_SKIP_WORKERS") ? (unsigned)atoi(getenv("BQD_SKIP_WORKERS")) : 4;
    const unsigned repeat = getenv("BQD_SKIP_REPEAT") ? (unsigned)atoi(getenv("BQD_SKIP_REPEAT")) : 1;
    /* one family alone, by its number, or 3 for the sessions, when chasing a failure */
    const int only = getenv("BQD_SKIP_FAMILY") ? atoi(getenv("BQD_SKIP_FAMILY")) : -1;
    walked = calloc(1, sizeof(nodeset_t));
    checked = calloc(1, sizeof(nodeset_t));
    if (walked == NULL || checked == NULL) { fprintf(stderr, "out of memory\n"); return 1; }
    int bad = 0;
    /* The collection tests keep weight indices across limdd_gc: the dense
     * vectors of the kept results, pw[] and isq2. A collection that finds the
     * weight table more than half full copies the live weights to new
     * indices, which those do not follow, and the checks would then read
     * stale weights, report wrong vectors or crash. The 2^22 entries here are
     * far from half full, which a larger test would reach, so the copy is off. */
    setenv("LIMDD_NO_WGT_GC", "1", 1);
    for (unsigned r = 0; r < repeat; r++) {
        for (int fam = BQD_FAMILY_SCALAR; fam <= BQD_FAMILY_PAULI; fam++) {
            if (only >= 0 && fam != only) continue;
            printf("== %s, exact weights, %u workers ==\n", bqd_family_name((bqd_family_t)fam), workers);
            session_begin(workers, 20);
            const int res = RUN(run_family, fam);
            session_end();
            printf("  %s\n", res ? "FAILED" : "ok");
            bad |= res;
        }
        for (int fam = BQD_FAMILY_SCALAR; fam <= BQD_FAMILY_PAULI; fam++) {
            if (only >= 0 && fam != only) continue;
            /* a small memo, which a collection clears, and which overwrites
             * entries more often in between */
            printf("== %s, collections between operations, %u workers ==\n",
                   bqd_family_name((bqd_family_t)fam), workers);
            session_begin(workers, 16);
            const int res = fam == BQD_FAMILY_SCALAR ? RUN(run_collections)
                                                     : RUN(run_collections_xp, fam);
            session_end();
            printf("  %s\n", res ? "FAILED" : "ok");
            bad |= res;
        }
    }
    if (only < 0 || only == 3) {
        printf("== two sessions in one Sylvan package, contention and the trigger, %u workers ==\n",
               workers);
        session_begin(workers, 20);
        const int res = RUN(run_sessions) | RUN(run_contention) | RUN(run_pauli_trigger);
        session_end();
        printf("  %s\n", res ? "FAILED" : "ok");
        bad |= res;
    }
    free(walked); free(checked);
    return bad;
}
