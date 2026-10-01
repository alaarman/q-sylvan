/*
 * Rule SM of the BQD's scalar family, Shannon where misaligned (skip:sec:sm
 * of the note "Binary Quotient Diagrams with Level Skipping", bqd-skip.tex,
 * whose labels are cited as skip:...), checked so that a node with the wrong
 * tag or flag, a level kept or skipped wrongly, or a wrong result cannot pass:
 *
 *   build         decode(build(f)) = f by bqd_to_vector and by a decoder
 *                 written here from skip:def:sm, bqd_eval agrees at every
 *                 index, building twice gives one edge, the full-support test
 *                 is right, and the cofactor at every level down to an edge's
 *                 node is the edge of the dense cofactor: on every vector over
 *                 {0, 1, -1, i} on 1 to 3 qubits, and on 4 to 8 qubits on
 *                 states with zeros that ignore random qubits, Clifford+T
 *                 circuit states, phase, product and coset states; a basis
 *                 state is bqd_basis_state's, a node per qubit
 *   structure     every node of every diagram here, built or computed: no
 *                 level whose two cofactors are equal (skip:def:sm); the tag
 *                 S exactly where the low cofactor is zero or the high
 *                 support is not inside the low one; a Q node with a nonzero
 *                 low edge and a ratio that is zero off the low support; an S
 *                 node with a nonzero high edge, which is the identity where
 *                 the low edge is zero; a low edge the identity or zero; the
 *                 full flag exactly where the node's function has no zero, no
 *                 other flag bit; the node the one bqd_from_vector builds for
 *                 its function; Cof and Compose give the node back from its
 *                 cofactors, and mk_nested or the Shannon rebuild from its two
 *                 stored edges; Ind is the node of the indicator; and at every
 *                 level an edge skips, the edge is both of its cofactors and
 *                 composes with itself to itself (skip:lem:virtual); Dep at
 *                 the root against the variables of the nodes it reaches
 *   operations    the product (bqd_multiply and bqd_product), the sum, the
 *                 ratio where the support allows it, Restr, Sub, and Compose
 *                 of two cofactors at level 0: on every pair of vectors on 1
 *                 and 2 qubits, a random partner for every vector on 3, and
 *                 random pairs on 4 to 8 qubits, among them pairs with one
 *                 ratio at the top and pairs whose supports nest; each result
 *                 the edge bqd_from_vector builds for the dense result, and a
 *                 diagram that passes the structure checks; the reference
 *                 stacks as they were after each; the shortcut cases of Apply
 *                 counted and met
 *   gates         every gate of the table (the projectors among them) and a
 *                 random exact 2x2, general, diagonal or a flip, at every
 *                 qubit; every set of controls above every target with X, Y,
 *                 Z, H, S, T, sqrt X and the random one; a phase with
 *                 controls anywhere; bqd_cgate_either with X, Y, Z and H both
 *                 ways, the Hadamard with the control below refused and the
 *                 state left alone; bqd_swap and bqd_perm of the three kinds
 *                 at every pair; bqd_x, restriction and projection at every
 *                 qubit; Pair at every qubit and both restrictions, with a
 *                 partner and with a zero operand; the diagonal gate on every
 *                 set of qubits by bqd_apply_diagonal, whose walk on full
 *                 support visits exactly the path of skip:prop:diag and the
 *                 same edge as the memoised Diag, and by bqd_diag_any; the
 *                 monomial; the local matvec on one and two qubits; MulOff on
 *                 a partner's support, on none and on all; PhaseMul and Exp
 *                 with a random exponent diagram, a root of unity or
 *                 1 + sqrt2: on every vector over {0, 1, -1, i} on 1 and 2
 *                 qubits, one vector in eight on 3 with one of each drawn, and
 *                 random states on 4 to 8 qubits; each result the edge
 *                 bqd_from_vector builds for the dense result and a diagram
 *                 that passes the structure checks, and the reference stacks
 *                 as they were. The cases of Gate at the root (the support on
 *                 one side, the ratio passed through, an S node above, a
 *                 ratio that depends on the target, a skipped target or
 *                 control), of Pair (both nested at the top, an S node, a
 *                 zero operand), of Diag and PhaseMul on S nodes and of
 *                 MulOff are counted and met; Side and iota of every node are
 *                 checked against its support among the structure checks
 *   circuits      random circuits of the qasm runner's gate set, called as
 *                 the runner calls them, on 2 to 8 qubits, the state checked
 *                 after every gate
 *   the two rules a session of bqd_init_rule with the copy rule runs it and
 *                 makes no node with a flag, and bqd_init gives the scalar
 *                 family SM and the other two the copy rule; on every function
 *                 of full support and every coset state the two rules have
 *                 the same number of nodes at every level (skip:cor:smd0); on
 *                 supports that misalign the counts are printed, not checked;
 *                 and bqd_init_rule with SM and the translation or Pauli
 *                 family exits with a message, in a child process
 *   collections   limdd_gc between operations in small tables, on protected
 *                 results on 8 and 6 qubits, each result checked as above, the
 *                 gates, selections and phases among the operations, and
 *                 after every collection every kept edge decoded and built
 *                 again
 *
 * Exact weights, in Q(w_8, sqrt2), throughout, on the number of workers in
 * BQD_SM_WORKERS (default 4).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

#include <sylvan.h>
#include <sylvan_int.h>

#include "qsylvan_bqd.h"
#include "qsylvan_bqd_ops.h"
#include "qsylvan_bqd_gates.h"
#include "qsylvan_bqd_sm.h"
#include "qsylvan_limdd_gc.h"
#include "qsylvan_limdd_inspect.h"
#include "qsylvan_simulator.h"
#include "sylvan_edge_weights_qisq2.h"

#define NQ 8
#define MAXV (1u << NQ)

/* Under a sanitizer the random draws are thinned by THIN; every kind of
 * check still runs, and every case it counts is still met. */
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
    printf("  %-60s %s (%u)\n", what, bad ? "FAILED" : "ok", tot);
}

/** A case the checks must have met, or they say nothing about it. */
static void
covered(const char *what, unsigned seen)
{
    expect(seen > 0, what, "never met");
    printf("  %-60s %u\n", what, seen);
}

/** Qubit q as a vector-index bit: qubit 0 is the most significant of n. */
static inline uint64_t
ibit(unsigned q, unsigned n)
{
    return UINT64_C(1) << (n - 1 - q);
}

/** A set of qubits, bit q for qubit q, as an index mask. */
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

static inline LIMDD_TARG
node_or_zero(BQD e)
{
    return limdd_edge_is_zero(e) ? 0 : limdd_target(e);
}

/** The level of a node: its variable, or n for the terminal of an n-qubit vector. */
static inline unsigned
level_of(LIMDD_TARG p, unsigned n)
{
    return p == LIMDD_TERMINAL ? n : limdd_node_var(p);
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

static bool
is_zero_vector(const EVBDD_WGT *f, uint64_t len)
{
    for (uint64_t x = 0; x < len; x++) if (f[x] != EVBDD_ZERO) return false;
    return true;
}

/** supp f inside supp g */
static bool
supp_inside(const EVBDD_WGT *f, const EVBDD_WGT *g, uint64_t len)
{
    for (uint64_t x = 0; x < len; x++) if (f[x] != EVBDD_ZERO && g[x] == EVBDD_ZERO) return false;
    return true;
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

/* --- the vector of a diagram, from skip:def:sm ----------------------------- */

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
 * The function of node p read at its own level v, 2^(n-v) values, by rule
 * SM with skipping edges: f_0 = l_0 . ext[N_0], and f_1 = l_1 . ext[N_0] .
 * ext[N_1] for a Q node and f_1 = l_1 . ext[N_1] for an S node, the plain
 * product. Written from the definition, not from the SM decoder, so that the
 * two check each other.
 */
static void
node_dense(LIMDD_TARG p, unsigned n, EVBDD_WGT *out)
{
    if (p == LIMDD_TERMINAL) { out[0] = EVBDD_ONE; return; }
    const unsigned v = limdd_node_var(p);
    const uint64_t h = UINT64_C(1) << (n - v - 1);
    const LIMDD low = limdd_node_low(p), high = limdd_node_high(p);
    const bool s = (limdd_node_flags(p) & BQD_SM_S) != 0;
    EVBDD_WGT *u1 = malloc(h * sizeof(EVBDD_WGT));
    if (limdd_edge_is_zero(low)) {
        zero_fill(out, h);
    } else {
        ext_dense(limdd_target(low), v + 1, n, out);
        const EVBDD_WGT c0 = limdd_lim_weight(limdd_label(low));
        for (uint64_t y = 0; y < h; y++) if (out[y] != EVBDD_ZERO) out[y] = wgt_mul(c0, out[y]);
    }
    if (limdd_edge_is_zero(high)) zero_fill(u1, h); else ext_dense(limdd_target(high), v + 1, n, u1);
    const EVBDD_WGT c1 = limdd_lim_weight(limdd_label(high));
    for (uint64_t y = 0; y < h; y++) {
        EVBDD_WGT w = EVBDD_ZERO;
        if (u1[y] != EVBDD_ZERO && (s || out[y] != EVBDD_ZERO))
            w = wgt_mul(c1, s ? u1[y] : wgt_mul(out[y], u1[y]));
        out[h + y] = w;
    }
    free(u1);
}

/** The function the edge e denotes read at level 0 of an n-qubit vector. */
static void
edge_dense(BQD e, unsigned n, EVBDD_WGT *out)
{
    const uint64_t len = UINT64_C(1) << n;
    if (limdd_edge_is_zero(e)) { zero_fill(out, len); return; }
    ext_dense(limdd_target(e), 0, n, out);
    const EVBDD_WGT c = limdd_lim_weight(limdd_label(e));
    for (uint64_t y = 0; y < len; y++) if (out[y] != EVBDD_ZERO) out[y] = wgt_mul(c, out[y]);
}

/* --- the structure of a diagram ------------------------------------------- */

typedef struct {
    unsigned nodes;         /* distinct nodes checked */
    unsigned redundant;     /* a node whose two cofactors are equal */
    unsigned tag;           /* the tag is not the one skip:def:sm gives */
    unsigned q_ratio;       /* a Q node whose ratio is not zero off its low support */
    unsigned s_edges;       /* an S node with a zero high edge, or a zero low and a labelled high */
    unsigned low_label;     /* a low edge neither the identity nor zero */
    unsigned full;          /* the full flag where the function has a zero, or not where it has none */
    unsigned flags;         /* a flag bit that is neither S nor full */
    unsigned not_canon;     /* a node that is not bqd_from_vector's for its function */
    unsigned cof;           /* Cof or Compose disagrees at a node or a skipped level */
    unsigned rebuild;       /* mk_nested or the Shannon rebuild does not give the node back */
    unsigned ind;           /* Ind is not the node of the indicator */
    unsigned iota;          /* iota is not the indicator as an exponent */
    unsigned side;          /* Side disagrees with the support at some qubit */
    unsigned dep;           /* Dep disagrees with the variables a root reaches */
    unsigned q_nodes, s_nodes, s_zero_low, q_zero_high, full_nodes, skipping;   /* met */
} structure_t;

static structure_t sf;

/** The checks on an edge read at level r, about the levels it skips. */
static void
check_edge(BQD e, unsigned r, unsigned n)
{
    if (limdd_edge_is_zero(e)) return;
    const unsigned v = level_of(limdd_target(e), n);
    if (v <= r) return;
    sf.skipping++;
    for (unsigned m = r; m < v; m++) {
        const BQD c0 = bqd_cofactor(e, m, 0), c1 = bqd_cofactor(e, m, 1);
        if (c0 != e || c1 != e || bqd_compose(m, c0, c1) != e) { sf.cof++; break; }
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
    const uint32_t flags = limdd_node_flags(p);
    const bool s = (flags & BQD_SM_S) != 0;
    const bool lz = limdd_edge_is_zero(low), hz = limdd_edge_is_zero(high);
    if (flags & ~(BQD_SM_S | BQD_SM_FULL)) sf.flags++;
    if (!lz && limdd_label(low) != LIMDD_LIM_IDENTITY) sf.low_label++;
    check_edge(low, v + 1, n);
    check_edge(high, v + 1, n);

    /* its function read at v, extended to all n qubits, and the two cofactors at v */
    const uint64_t len = UINT64_C(1) << n, hb = ibit(v, n);
    EVBDD_WGT *g = malloc(len * sizeof(EVBDD_WGT));
    EVBDD_WGT *g0 = malloc(len * sizeof(EVBDD_WGT)), *g1 = malloc(len * sizeof(EVBDD_WGT));
    ext_dense(p, 0, n, g);
    for (uint64_t x = 0; x < len; x++) { g0[x] = g[x & ~hb]; g1[x] = g[x | hb]; }
    if (same_vector(g0, g1, n)) sf.redundant++;

    const bool nested = !is_zero_vector(g0, len) && supp_inside(g1, g0, len);
    if (s == nested) sf.tag++;
    if (s) {
        sf.s_nodes++;
        if (hz) sf.s_edges++;
        if (lz) {
            sf.s_zero_low++;
            if (!hz && limdd_label(high) != LIMDD_LIM_IDENTITY) sf.s_edges++;
        }
    } else {
        sf.q_nodes++;
        if (lz) sf.tag++;
        if (hz) {
            sf.q_zero_high++;
        } else if (!lz) {
            /* the ratio, as a function of the levels below v, is zero off supp [A] */
            EVBDD_WGT *a = malloc(len * sizeof(EVBDD_WGT)), *r = malloc(len * sizeof(EVBDD_WGT));
            edge_dense(unit(limdd_target(low)), n, a);
            edge_dense(unit(limdd_target(high)), n, r);
            if (!supp_inside(r, a, len)) sf.q_ratio++;
            free(a); free(r);
        }
    }

    bool no_zero = true;
    for (uint64_t x = 0; x < len && no_zero; x++) no_zero = g[x] != EVBDD_ZERO;
    if (((flags & BQD_SM_FULL) != 0) != no_zero) sf.full++;
    if (no_zero) sf.full_nodes++;

    if (bqd_from_vector(g, n) != unit(p)) sf.not_canon++;

    /* its cofactors are the canonical edges of the two halves, and give it back */
    const BQD c0 = bqd_cofactor(unit(p), v, 0), c1 = bqd_cofactor(unit(p), v, 1);
    if (c0 != bqd_from_vector(g0, n) || c1 != bqd_from_vector(g1, n)
        || bqd_compose(v, c0, c1) != unit(p)) sf.cof++;

    /* and so do its two stored edges, by the constructor of its tag, and the
     * two scaled, which scales the node: the scale comes out onto the edge */
    const BQD back = s ? bqd_sm_rebuild_s(v, low, high) : bqd_sm_mk_nested(v, low, high);
    const EVBDD_WGT c = wgt_mul(isq2, pw[1 + p % 7]);
    const BQD cl = bqd_scale(low, c), ch = bqd_scale(high, c);
    const BQD scaled = s ? bqd_sm_rebuild_s(v, cl, ch) : bqd_sm_mk_nested(v, cl, high);
    if (back != unit(p) || scaled != bqd_scale(unit(p), c)) sf.rebuild++;

    /* Side at every qubit: where the support lies on one side of x_q; 2
     * above the node, where the function does not depend on x_q */
    for (unsigned q = 0; q < n; q++) {
        bool at0 = false, at1 = false;
        for (uint64_t x = 0; x < len; x++)
            if (g[x] != EVBDD_ZERO) { if (x & ibit(q, n)) at1 = true; else at0 = true; }
        const int want = (at0 && at1) ? 2 : (at1 ? 1 : 0);
        if (bqd_sm_side(p, q) != want) { sf.side++; break; }
    }

    /* Ind is the node of the indicator, and a node of its own kind, and iota
     * the indicator as an exponent; a point y of an exponent has qubit v at
     * bit v, where an index has it at bit n-1-v */
    for (uint64_t x = 0; x < len; x++) g[x] = (g[x] != EVBDD_ZERO) ? EVBDD_ONE : EVBDD_ZERO;
    if (unit(bqd_sm_ind(p)) != bqd_from_vector(g, n)) sf.ind++;
    const BQD_EXP io = mtbdd_refs_push(RUN(bqd_sm_iota, p));
    for (uint64_t x = 0; x < len; x++) {
        uint64_t y = 0;
        for (unsigned v = 0; v < n; v++) if (x & ibit(v, n)) y |= UINT64_C(1) << v;
        if (bqd_exp_eval(io, y) != (g[x] == EVBDD_ONE ? 1 : 0)) { sf.iota++; break; }
    }
    mtbdd_refs_pop(1);

    free(g); free(g0); free(g1);
    if (!lz) check_node(limdd_target(low), n);
    if (!hz) check_node(limdd_target(high), n);
}

/** The variables of the nodes reachable from p, bit q for qubit q, p included. */
static uint64_t
reach_vars(LIMDD_TARG p)
{
    if (p == LIMDD_TERMINAL || !set_add(walked, p)) return 0;
    uint64_t m = UINT64_C(1) << limdd_node_var(p);
    const LIMDD lo = limdd_node_low(p), hi = limdd_node_high(p);
    if (!limdd_edge_is_zero(lo)) m |= reach_vars(limdd_target(lo));
    if (!limdd_edge_is_zero(hi)) m |= reach_vars(limdd_target(hi));
    return m;
}

/** The structure of e, a diagram on n qubits. */
static void
check_diagram(BQD e, unsigned n)
{
    if (limdd_edge_is_zero(e)) return;
    check_edge(e, 0, n);
    check_node(limdd_target(e), n);
    set_clear(walked);
    const uint64_t m = reach_vars(limdd_target(e));
    for (unsigned q = 0; q < n; q++)
        if (bqd_sm_dep(limdd_target(e), q) != (int)((m >> q) & 1)) { sf.dep++; break; }
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
        { "a node whose two cofactors are equal", sf.redundant },
        { "a tag that is not skip:def:sm's", sf.tag },
        { "a Q node whose ratio is not zero off its low support", sf.q_ratio },
        { "an S node with a zero high edge or a labelled one over a zero low", sf.s_edges },
        { "a low edge neither the identity nor zero", sf.low_label },
        { "a full flag that is wrong", sf.full },
        { "a flag bit that is neither S nor full", sf.flags },
        { "a node that is not the one bqd_from_vector builds for its function", sf.not_canon },
        { "Cof or Compose disagrees at a node or a skipped level", sf.cof },
        { "mk_nested or the Shannon rebuild does not give a node back", sf.rebuild },
        { "Ind is not the node of the indicator", sf.ind },
        { "iota is not the indicator as an exponent", sf.iota },
        { "Side disagrees with the support at a qubit", sf.side },
        { "Dep disagrees with the variables a root reaches", sf.dep },
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
    printf("  %-60s %s (%u nodes: %u Q, %u S, %u full; %u skipping edges)\n", buf,
           bad ? "FAILED" : "ok", sf.nodes, sf.q_nodes, sf.s_nodes, sf.full_nodes, sf.skipping);
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
    for (uint64_t w = 1; w < (UINT64_C(1) << n) && nt < 96; w++)
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

/** A tensor product of random one-qubit states, |0>, |1>, |+> and |-> among them. */
static void
product_state(unsigned n, bool zeros, EVBDD_WGT *f)
{
    EVBDD_WGT a[NQ + 2][2];
    for (unsigned q = 0; q < n; q++) {
        switch (rnd_below(zeros ? 6 : 4)) {
        case 0:  a[q][0] = EVBDD_ONE;  a[q][1] = EVBDD_ONE;     break;
        case 1:  a[q][0] = EVBDD_ONE;  a[q][1] = EVBDD_MIN_ONE; break;
        case 2:  a[q][0] = EVBDD_ONE;  a[q][1] = pw[2];         break;
        case 3:  a[q][0] = rand_value(0.0, true); a[q][1] = rand_value(0.0, true); break;
        case 4:  a[q][0] = EVBDD_ONE;  a[q][1] = EVBDD_ZERO;    break;
        default: a[q][0] = EVBDD_ZERO; a[q][1] = EVBDD_ONE;     break;
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

/**
 * [x in a + C] . w_8^{P(x)} for a random linear code C of dimension k, a
 * random shift a and a random phase polynomial P of degree at most 3: a
 * function whose support is an affine subspace.
 */
static void
coset_state(unsigned n, EVBDD_WGT *f)
{
    const uint64_t len = UINT64_C(1) << n;
    uint8_t *in = calloc(len, 1);
    in[0] = 1;
    const unsigned k = (unsigned)rnd_below(n + 1);
    for (unsigned i = 0; i < k; i++) {
        const uint64_t b = rnd_below(len);
        for (uint64_t x = 0; x < len; x++) if (in[x] == 1 && !in[x ^ b]) in[x ^ b] = 2;
        for (uint64_t x = 0; x < len; x++) if (in[x] == 2) in[x] = 1;
    }
    const uint64_t a = rnd_below(len);
    phase_on(n, rand_qubits(n), f);
    for (uint64_t x = 0; x < len; x++) if (!in[x ^ a]) f[x] = EVBDD_ZERO;
    free(in);
}

enum { S_ZEROS, S_ALGEBRAIC, S_PHASE, S_CIRCUIT, S_PRODUCT, S_BASIS_PLUS, S_COSET, S_KINDS };

/**
 * A state of one of seven kinds: a function of random qubits with zeros, the
 * same over Z[sqrt2, i]/2 without zeros, a phase state of random qubits, a
 * Clifford+T circuit state, a product state, a basis state on random qubits
 * with |+> on the rest, and a coset state.
 */
static void
rand_state(unsigned n, unsigned kind, EVBDD_WGT *f)
{
    switch (kind) {
    case S_ZEROS:     rand_on(n, rand_qubits(n), 0.2 * (double)(1 + rnd_below(3)), rnd_below(2), f); break;
    case S_ALGEBRAIC: rand_on(n, rand_qubits(n), 0.0, true, f); break;
    case S_PHASE:     phase_on(n, rand_qubits(n), f); break;
    case S_CIRCUIT:   circuit_state(n, f); break;
    case S_PRODUCT:   product_state(n, true, f); break;
    case S_COSET:     coset_state(n, f); break;
    default: {
        const uint64_t D = imask(rand_qubits(n), n), x0 = rnd_below(UINT64_C(1) << n) & D;
        for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
            f[x] = ((x & D) == x0) ? EVBDD_ONE : EVBDD_ZERO;
    }
    }
}

/* --- building ------------------------------------------------------------- */

typedef struct {
    unsigned tot, decode, dense, eval, again, full, cof;
    unsigned skip, root_skip, zero_low_root;
} builds_t;

/** The cofactor of f at qubit v set to b, as a vector on all n qubits. */
static void
dense_cofactor(const EVBDD_WGT *f, unsigned n, unsigned v, int b, EVBDD_WGT *out)
{
    const uint64_t qb = ibit(v, n);
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) out[x] = f[b ? (x | qb) : (x & ~qb)];
}

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
    edge_dense(e, n, g);
    if (!same_vector(f, g, n)) b->dense++;
    for (uint64_t x = 0; x < len; x++) if (bqd_eval(e, n, x) != f[x]) { b->eval++; break; }
    if (bqd_from_vector(f, n) != e) b->again++;
    bool no_zero = true;
    for (uint64_t x = 0; x < len && no_zero; x++) no_zero = f[x] != EVBDD_ZERO;
    if (bqd_has_full_support(e) != no_zero) b->full++;

    /* the cofactor at every level down to the edge's node */
    if (!limdd_edge_is_zero(e)) {
        const unsigned t = level_of(limdd_target(e), n);
        for (unsigned v = 0; v <= t && v < n; v++) for (int c = 0; c < 2; c++) {
            dense_cofactor(f, n, v, c, g);
            if (bqd_cofactor(e, v, c) != bqd_from_vector(g, n)) { b->cof++; v = n; break; }
        }
        if (t > 0) b->root_skip++;
        if (t < n && limdd_edge_is_zero(limdd_node_low(limdd_target(e)))) b->zero_low_root++;
    }
    const unsigned before = sf.skipping;
    check_diagram(e, n);
    if (sf.skipping > before) b->skip++;
    free(g);
    return e;
}

static void
builds_report(const char *section, const builds_t *b)
{
    char what[96];
    snprintf(what, sizeof(what), "%s: decode(build(f)) = f", section);
    report(what, b->decode, b->tot);
    snprintf(what, sizeof(what), "%s: the decoder of skip:def:sm", section);
    report(what, b->dense, b->tot);
    snprintf(what, sizeof(what), "%s: bqd_eval at every index", section);
    report(what, b->eval, b->tot);
    snprintf(what, sizeof(what), "%s: built twice, one edge", section);
    report(what, b->again, b->tot);
    snprintf(what, sizeof(what), "%s: full support, from the flag", section);
    report(what, b->full, b->tot);
    snprintf(what, sizeof(what), "%s: the cofactor at every level to the node", section);
    report(what, b->cof, b->tot);
    snprintf(what, sizeof(what), "%s: diagrams with a skipping edge (new nodes)", section);
    covered(what, b->skip);
    snprintf(what, sizeof(what), "%s: root edges that skip", section);
    covered(what, b->root_skip);
    snprintf(what, sizeof(what), "%s: roots with a zero low cofactor", section);
    covered(what, b->zero_low_root);
}

/** Every vector over {0, 1, -1, i} on 1 to 3 qubits: 16, 256 and 65536 of them. */
static void
check_exhaustive_builds(void)
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
    covered("n <= 3, all: S nodes over a zero low cofactor", sf.s_zero_low);
    covered("n <= 3, all: S nodes over a misaligned pair", sf.s_nodes - sf.s_zero_low);
    covered("n <= 3, all: Q nodes with no high cofactor", sf.q_zero_high);
}

static void
check_random_builds(void)
{
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT));
    builds_t b;
    memset(&b, 0, sizeof(b));
    structure_begin();
    for (unsigned n = 4; n <= NQ; n++) for (int rep = 0; rep < 8; rep++)
        for (unsigned kind = 0; kind < S_KINDS; kind++) {
            rand_state(n, kind, f);
            build_checked(f, n, &b);
        }
    builds_report("n = 4..8, random", &b);
    structure_report("n = 4..8, random");

    /* a basis state has a node per qubit and is bqd_basis_state's */
    unsigned bad = 0, tot = 0;
    size_t counts[NQ];
    for (unsigned n = 1; n <= NQ; n++) for (int rep = 0; rep < 8; rep++) {
        const uint64_t y = rnd_below(UINT64_C(1) << n);
        for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) f[x] = (x == y) ? EVBDD_ONE : EVBDD_ZERO;
        const BQD e = bqd_from_vector(f, n);
        limdd_level_counts(e, counts, n);
        bool ok = e == bqd_basis_state(y, n);
        for (unsigned q = 0; q < n; q++) ok = ok && counts[q] == 1;
        tot++;
        if (!ok) bad++;
    }
    report("a basis state: a node per level, bqd_basis_state", bad, tot);
    free(f);
}

/* --- operations ------------------------------------------------------------ */

typedef struct { unsigned bad, tot; } tally_t;

/** Count got wrong unless it is the edge bqd_from_vector builds for h, and check its structure. */
static void
tally(tally_t *t, BQD got, const EVBDD_WGT *h, unsigned n)
{
    t->tot++;
    check_diagram(got, n);
    if (got != bqd_from_vector(h, n)) t->bad++;
}

typedef struct {
    tally_t mul, prod, add, ratio, restr, sub, compose, refs;
    unsigned nested_pair, s_pair, one_ratio, ratio_pairs, subset_true, subset_false;
} ops_t;

/** Whether [N] read at level v is nested there, as Apply decides it: Q, below v, or the terminal. */
static bool
nested_at(LIMDD_TARG N, unsigned v, unsigned n)
{
    return level_of(N, n) > v || !(limdd_node_flags(N) & BQD_SM_S);
}

/**
 * Every operation on the pair f, g on n qubits, with their edges ef and eg.
 * Compose is of the cofactor x_0 = 0 of f and x_0 = 1 of g, which do not
 * depend on qubit 0, and makes the function that agrees with f where x_0 is 0
 * and with g where it is 1.
 */
static void
ops_on_pair(const EVBDD_WGT *f, const EVBDD_WGT *g, BQD ef, BQD eg, unsigned n, ops_t *o)
{
    const uint64_t len = UINT64_C(1) << n;
    EVBDD_WGT h[MAXV];
    size_t d0[3], d1[3];
    limdd_refs_depths(d0);

    /* the shortcut cases of Apply, from the operands */
    const LIMDD_TARG Nf = node_or_zero(ef), Ng = node_or_zero(eg);
    if (Nf != 0 && Ng != 0 && Nf != LIMDD_TERMINAL && Ng != LIMDD_TERMINAL && Nf != Ng) {
        const unsigned lf = level_of(Nf, n), lg = level_of(Ng, n), v = lf < lg ? lf : lg;
        if (nested_at(Nf, v, n) && nested_at(Ng, v, n)) {
            o->nested_pair++;
            BQD a, ra, b, rb;
            RUN(bqd_sm_nest, Nf, v, &a, &ra);
            RUN(bqd_sm_nest, Ng, v, &b, &rb);
            if (ra == rb) o->one_ratio++;
        } else {
            o->s_pair++;
        }
    }

    for (uint64_t x = 0; x < len; x++) h[x] = wgt_mul(f[x], g[x]);
    tally(&o->mul, bqd_multiply(ef, eg), h, n);
    tally(&o->prod, bqd_product(ef, eg), h, n);

    for (uint64_t x = 0; x < len; x++) h[x] = wgt_add(f[x], g[x]);
    tally(&o->add, bqd_add(ef, eg), h, n);

    /* the ratio, where supp g lies inside supp f */
    if (supp_inside(g, f, len)) {
        o->ratio_pairs++;
        for (uint64_t x = 0; x < len; x++)
            h[x] = (f[x] == EVBDD_ZERO || g[x] == EVBDD_ZERO) ? EVBDD_ZERO : wgt_div(g[x], f[x]);
        tally(&o->ratio, bqd_sm_apply(BQD_SM_RATIO, ef, eg), h, n);
    }

    /* Sub, and Restr where g is not zero */
    const bool in = supp_inside(f, g, len);
    if (in) o->subset_true++; else o->subset_false++;
    o->sub.tot++;
    if (bqd_sm_sub(Nf, Ng) != (int)in) o->sub.bad++;
    if (!limdd_edge_is_zero(eg)) {
        for (uint64_t x = 0; x < len; x++) h[x] = (g[x] == EVBDD_ZERO) ? EVBDD_ZERO : f[x];
        tally(&o->restr, bqd_sm_restr(ef, eg), h, n);
    }

    /* Compose at level 0 of a cofactor of each */
    const uint64_t b0 = ibit(0, n);
    for (uint64_t x = 0; x < len; x++) h[x] = (x & b0) ? g[x] : f[x];
    tally(&o->compose, bqd_compose(0, bqd_cofactor(ef, 0, 0), bqd_cofactor(eg, 0, 1)), h, n);

    limdd_refs_depths(d1);
    o->refs.tot++;
    if (memcmp(d0, d1, sizeof(d0)) != 0) o->refs.bad++;
}

static void
ops_report(const char *section, const ops_t *o)
{
    char what[112];
    const struct { const char *name; const tally_t *t; } k[] = {
        { "bqd_multiply", &o->mul }, { "bqd_product", &o->prod }, { "bqd_add", &o->add },
        { "the ratio, supp g inside supp f", &o->ratio }, { "Restr", &o->restr },
        { "Sub against the dense supports", &o->sub },
        { "Compose of two cofactors at level 0", &o->compose },
        { "the reference stacks as they were", &o->refs },
    };
    for (size_t i = 0; i < sizeof(k) / sizeof(k[0]); i++) {
        snprintf(what, sizeof(what), "%s: %s", section, k[i].name);
        report(what, k[i].t->bad, k[i].t->tot);
    }
    snprintf(what, sizeof(what), "%s: pairs nested at the top (Apply on stored edges)", section);
    covered(what, o->nested_pair);
    snprintf(what, sizeof(what), "%s: pairs with one ratio at the top", section);
    covered(what, o->one_ratio);
    snprintf(what, sizeof(what), "%s: pairs with an S node at the top", section);
    covered(what, o->s_pair);
    snprintf(what, sizeof(what), "%s: supports that nest, and that do not", section);
    covered(what, o->subset_true < o->subset_false ? o->subset_true : o->subset_false);
}

/** Every pair of vectors over {0, 1, -1, i} on 1 and 2 qubits, and a random partner on 3. */
static void
check_exhaustive_ops(void)
{
    const EVBDD_WGT vals[4] = { EVBDD_ZERO, EVBDD_ONE, EVBDD_MIN_ONE, pw[2] };
    static EVBDD_WGT vec[256][4];
    static BQD e[256];
    EVBDD_WGT f[8], g[8];
    ops_t o;
    memset(&o, 0, sizeof(o));
    structure_begin();
    for (unsigned n = 1; n <= 2; n++) {
        const uint64_t len = UINT64_C(1) << n, count = UINT64_C(1) << (2 * len);
        for (uint64_t i = 0; i < count; i++) {
            for (uint64_t x = 0; x < len; x++) vec[i][x] = vals[(i >> (2 * x)) & 3];
            e[i] = bqd_from_vector(vec[i], n);
        }
        for (uint64_t i = 0; i < count; i++)
            for (uint64_t j = 0; j < count; j++)
                ops_on_pair(vec[i], vec[j], e[i], e[j], n, &o);
    }
    for (uint64_t i = 0; i < 65536; i += THIN) {
        for (uint64_t x = 0; x < 8; x++) f[x] = vals[(i >> (2 * x)) & 3];
        const uint64_t j = rnd_below(65536);
        for (uint64_t x = 0; x < 8; x++) g[x] = vals[(j >> (2 * x)) & 3];
        ops_on_pair(f, g, bqd_from_vector(f, 3), bqd_from_vector(g, 3), 3, &o);
    }
    ops_report("n <= 2 all pairs, n = 3", &o);
    structure_report("n <= 3, operations");
}

/**
 * Random pairs on 4 to 8 qubits: two random states; a state and its product
 * with another, whose support nests in the first; and two states with one
 * ratio at the top, a full-support low cofactor each times one shared factor.
 */
static void
check_random_ops(void)
{
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *g = malloc(MAXV * sizeof(EVBDD_WGT));
    EVBDD_WGT *r = malloc(MAXV * sizeof(EVBDD_WGT));
    ops_t o;
    memset(&o, 0, sizeof(o));
    structure_begin();
    for (unsigned n = 4; n <= NQ; n++) for (int rep = 0; rep < 60 / THIN; rep++) {
        const uint64_t len = UINT64_C(1) << n;
        rand_state(n, (unsigned)rnd_below(S_KINDS), f);
        rand_state(n, (unsigned)rnd_below(S_KINDS), g);
        switch (rnd_below(3)) {
        case 0:
            break;
        case 1:                                     /* supp g inside supp f */
            for (uint64_t x = 0; x < len; x++) g[x] = wgt_mul(f[x], g[x]);
            break;
        default: {                                  /* one ratio at qubit 0 */
            const uint64_t b0 = ibit(0, n);
            rand_state(n, S_ZEROS, r);
            rand_on(n, rand_qubits(n) & ~UINT64_C(1), 0.0, true, f);
            rand_on(n, rand_qubits(n) & ~UINT64_C(1), 0.0, true, g);
            for (uint64_t x = 0; x < len; x++) if (x & b0) {
                f[x] = wgt_mul(f[x & ~b0], r[x & ~b0]);
                g[x] = wgt_mul(g[x & ~b0], r[x & ~b0]);
            }
        }
        }
        const BQD ef = limdd_refs_push(bqd_from_vector(f, n));
        const BQD eg = limdd_refs_push(bqd_from_vector(g, n));
        ops_on_pair(f, g, ef, eg, n, &o);
        ops_on_pair(g, f, eg, ef, n, &o);
        limdd_refs_pop(2);
    }
    ops_report("n = 4..8, random", &o);
    structure_report("n = 4..8, operations");
    free(f); free(g); free(r);
}

/* --- gates, selections and phases ------------------------------------------ */

static void
dense_project(const EVBDD_WGT *f, unsigned n, unsigned q, int b, EVBDD_WGT *out)
{
    const uint64_t qb = ibit(q, n);
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
        out[x] = (((x & qb) != 0) == (b != 0)) ? f[x] : EVBDD_ZERO;
}

static void
dense_x(const EVBDD_WGT *f, unsigned n, unsigned q, EVBDD_WGT *out)
{
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) out[x] = f[x ^ ibit(q, n)];
}

/** (pi f)(x) = f(pi x), for the permutation `kind` of the qubits qa < qb. */
static void
dense_perm(const EVBDD_WGT *f, unsigned n, uint32_t kind, unsigned qa, unsigned qb, EVBDD_WGT *out)
{
    const uint64_t ba = ibit(qa, n), bb = ibit(qb, n);
    for (uint64_t y = 0; y < (UINT64_C(1) << n); y++) {
        const int xa = (y & ba) != 0, xb = (y & bb) != 0;
        int na = xa, nb = xb;
        if (kind == BQD_PERM_SWAP) { na = xb; nb = xa; }
        else if (kind == BQD_PERM_CX_DOWN) nb = xb ^ xa;
        else na = xa ^ xb;
        out[y] = f[(y & ~(ba | bb)) | (na ? ba : 0) | (nb ? bb : 0)];
    }
}

/** [x_b = 0] A|_{x_b = s} + [x_b = 1] B|_{x_b = t}. */
static void
dense_pair(const EVBDD_WGT *A, int s, const EVBDD_WGT *B, int t, unsigned b, unsigned n,
           EVBDD_WGT *out)
{
    const uint64_t bb = ibit(b, n);
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
        out[x] = (x & bb) ? B[t ? (x | bb) : (x & ~bb)] : A[s ? (x | bb) : (x & ~bb)];
}

/** f . phase^{x_A}, A an index mask. */
static void
dense_diag(const EVBDD_WGT *f, unsigned n, uint64_t A, EVBDD_WGT phase, EVBDD_WGT *out)
{
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
        out[x] = ((x & A) == A && f[x] != EVBDD_ZERO) ? wgt_mul(f[x], phase) : f[x];
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

/** A matrix entry: zero one time in four, so that diagonal and flip matrices come up. */
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

/** Put U in the table as the dynamic gate. The gate memo is keyed on the id, so the cache goes. */
static void
set_dynamic(const EVBDD_WGT *U)
{
    for (int i = 0; i < 4; i++) gates[GATEID_dynamic][i] = U[i];
    sylvan_clear_cache();
}

static bool
is_phase_gate(uint32_t gid)
{
    const uint64_t *U = gates[gid];
    return U[1] == EVBDD_ZERO && U[2] == EVBDD_ZERO && U[0] != EVBDD_ZERO && U[3] != EVBDD_ZERO;
}

static bool
is_flip_gate(uint32_t gid)
{
    const uint64_t *U = gates[gid];
    return U[0] == EVBDD_ZERO && U[3] == EVBDD_ZERO && U[1] != EVBDD_ZERO && U[2] != EVBDD_ZERO;
}

/** Whether bqd_cgate with these controls goes to Gate, and not to the monomial or the flip. */
static bool
reaches_gate(uint32_t gid, uint64_t cmask)
{
    if (cmask == 0) return !is_phase_gate(gid) && !is_flip_gate(gid);
    if (is_phase_gate(gid) && gates[gid][0] == EVBDD_ONE) return false;
    return !(is_flip_gate(gid) && (cmask & (cmask - 1)) == 0);
}

/*
 * The exponent diagram of the table v, v[x] at the vector index x, reduced
 * modulo r: the Shannon expansion from the top qubit down, protected as it is
 * made, since a node can make Sylvan collect.
 */
static BQD_EXP
exp_table(const int64_t *v, unsigned var, unsigned n, uint32_t r)
{
    if (var == n) return bqd_exp_const(v[0], r);
    const uint64_t h = UINT64_C(1) << (n - var - 1);
    const BQD_EXP lo = mtbdd_refs_push(exp_table(v, var + 1, n, r));
    const BQD_EXP hi = mtbdd_refs_push(exp_table(v + h, var + 1, n, r));
    const BQD_EXP e = bqd_exp_node(var, lo, hi);
    mtbdd_refs_pop(2);
    return e;
}

/** A random exponent table of the qubits in D (bit q for qubit q), its values in [lo, lo + span). */
static void
rand_exponents(unsigned n, uint64_t D, int64_t lo, unsigned span, int64_t *v)
{
    const uint64_t dm = imask(D, n);
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
        v[x] = (x & ~dm) ? v[x & dm] : lo + (int64_t)rnd_below(span);
}

/**
 * The nodes of e the diagonal walk visits, from the diagram alone
 * (skip:prop:diag): the path from the root that goes low at a node above the
 * top qubit of A left and high at a node of A, and ends after the last qubit
 * of A, or at the first edge that skips the top qubit of A left, whose node
 * is not counted. A is a set of qubits, bit q for qubit q.
 */
static uint32_t
diag_path(BQD e, uint64_t A, unsigned n)
{
    uint32_t visits = 0;
    LIMDD_TARG t = limdd_target(e);
    while (A != 0) {
        const unsigned a = (unsigned)__builtin_ctzll(A);
        const unsigned v = level_of(t, n);
        if (v > a) break;
        visits++;
        if (v < a) { t = limdd_target(limdd_node_low(t)); continue; }
        A &= A - 1;
        if (A != 0) t = limdd_target(limdd_node_high(t));
    }
    return visits;
}

enum {
    G_GATE, G_DYN, G_CGATE, G_CPHASE, G_EITHER, G_SWAP, G_PERM, G_X, G_RESTRICT, G_PROJECT,
    G_PAIR, G_DIAG, G_DIAG_ANY, G_DIAG_MEMO, G_VISITS, G_MONO, G_MATVEC1, G_MATVEC2, G_MULOFF,
    G_PHASEMUL, G_EXP, G_REFS, G_COUNT
};
static const char *g_name[G_COUNT] = {
    "bqd_gate, table gates", "bqd_gate, a random exact 2x2", "bqd_cgate, controls above",
    "bqd_cgate, a phase with controls anywhere", "bqd_cgate_either", "bqd_swap",
    "bqd_perm, the three kinds", "bqd_x", "bqd_restrict", "bqd_project", "bqd_pair",
    "bqd_apply_diagonal", "bqd_diag_any", "Diag memoised = the walk, full support",
    "the walk's visits are skip:prop:diag's path", "bqd_monomial", "bqd_local_matvec, one qubit",
    "bqd_local_matvec, two qubits", "bqd_mul_off", "bqd_phase_mul", "bqd_exp_state",
    "the reference stacks as they were",
};

enum {
    V_SIDE, V_PASS, V_S_ABOVE, V_Q_ABOVE, V_SKIP_TARGET, V_SKIP_CONTROL, V_ZERO_LOW, V_CANCEL,
    V_EITHER_BELOW, V_REFUSED, V_PAIR_NESTED, V_PAIR_S, V_PAIR_ZERO, V_DIAG_S, V_DIAG_WALK,
    V_DIAG_SKIP, V_PHASE_S, V_MULOFF_PHASE, V_COUNT
};
static const char *v_name[V_COUNT] = {
    "Gate at the root: the support on one side of the target",
    "Gate at the root: a Q node whose ratio passes the gate",
    "Gate at the root: an S node above the target",
    "Gate at the root: a Q node whose ratio depends on the target",
    "Gate on a skipped target", "Gate under a skipped control",
    "gates on a qubit whose low cofactor is zero", "gates that cancel an amplitude",
    "cgate_either with the control below", "cgate_either refused (H, control below)",
    "Pair: both nested at the top", "Pair: an S node at the top", "Pair: a zero operand",
    "Diag: an S node above the top of A", "Diag: the walk, on full support",
    "Diag: the state skips the top of A", "PhaseMul: an S node at the top",
    "MulOff that multiplies, not a subset",
};

typedef struct {
    tally_t t[G_COUNT];
    unsigned v[V_COUNT];
} gops_t;

static size_t own_depth[6];

static void
own_mark(void)
{
    limdd_refs_depths(own_depth);
    mtbdd_refs_depths(own_depth + 3);
}

/** e, with the calling worker's six reference stacks as own_mark found them, counted in t */
static BQD
own_check(tally_t *t, BQD e)
{
    size_t d[6];
    limdd_refs_depths(d);
    mtbdd_refs_depths(d + 3);
    t->tot++;
    if (memcmp(d, own_depth, sizeof(d)) != 0) t->bad++;
    return e;
}

#define OWN(o, e) (own_mark(), own_check(&(o)->t[G_REFS], (e)))

/** Which case of Gate the root of F is in, for gid on q under the controls cmask. */
static void
gate_root_case(BQD F, uint32_t gid, uint64_t cmask, unsigned q, unsigned n, unsigned *v)
{
    if (limdd_edge_is_zero(F) || !reaches_gate(gid, cmask)) return;
    const LIMDD_TARG N = limdd_target(F);
    const unsigned t = level_of(N, n);
    const unsigned p = (unsigned)__builtin_ctzll(cmask | (UINT64_C(1) << q));
    if (t > p) { v[p == q ? V_SKIP_TARGET : V_SKIP_CONTROL]++; return; }
    if (t == p) return;
    const bool s_node = (limdd_node_flags(N) & BQD_SM_S) != 0;
    const int side = (p == q) ? bqd_sm_side(N, q) : 2;
    if (side != 2 && gates[gid][side] != EVBDD_ZERO && gates[gid][2 + side] != EVBDD_ZERO)
        v[V_SIDE]++;
    else if (!s_node && !bqd_sm_dep(node_or_zero(limdd_node_high(N)), q)) v[V_PASS]++;
    else v[s_node ? V_S_ABOVE : V_Q_ABOVE]++;
}

/** A zero low cofactor at q under a nonzero high one, and a cancellation, from f to h. */
static void
gate_cases(const EVBDD_WGT *f, const EVBDD_WGT *h, unsigned n, unsigned q, unsigned *v)
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
    if (low_zero && nonzero) v[V_ZERO_LOW]++;
    if (cancel) v[V_CANCEL]++;
}

/** Whether the top of N read at l is an S node there, l at or above N. */
static bool
s_at(LIMDD_TARG N, unsigned l, unsigned n)
{
    return N != LIMDD_TERMINAL && level_of(N, n) == l && (limdd_node_flags(N) & BQD_SM_S);
}

static const uint32_t table_gates[] = { GATEID_H, GATEID_X, GATEID_Y, GATEID_Z, GATEID_S,
                                        GATEID_Sdag, GATEID_T, GATEID_Tdag, GATEID_sqrtX,
                                        GATEID_sqrtY, GATEID_proj0, GATEID_proj1 };
static const uint32_t ctl_gates[] = { GATEID_X, GATEID_Y, GATEID_Z, GATEID_H, GATEID_S,
                                      GATEID_T, GATEID_sqrtX, GATEID_dynamic };
static const uint32_t phase_gates[] = { GATEID_Z, GATEID_S, GATEID_Sdag, GATEID_T, GATEID_Tdag };
static const uint32_t either_gates[] = { GATEID_X, GATEID_Y, GATEID_Z, GATEID_H };

#define NUM(a) (sizeof(a) / sizeof((a)[0]))

/** gid on q under cmask, against dense_gate, its root case and its cases counted. */
static void
one_gate(const EVBDD_WGT *f, BQD F, unsigned n, uint32_t gid, uint64_t cmask, unsigned q,
         tally_t *t, gops_t *o)
{
    EVBDD_WGT h[MAXV];
    memcpy(h, f, (UINT64_C(1) << n) * sizeof(EVBDD_WGT));
    dense_gate(h, n, gates[gid], cmask, q);
    gate_root_case(F, gid, cmask, q, n, o->v);
    if (cmask == 0) gate_cases(f, h, n, q, o->v);
    tally(t, OWN(o, cmask ? bqd_cgate(F, gid, cmask, q, n) : bqd_gate(F, gid, q, n)), h, n);
}

/**
 * The gates and selections on f, F its edge, each against its dense result:
 * with `every`, at every qubit, every pair and every set of controls above a
 * target, with every gate, and otherwise one of each, drawn. g and G are a
 * partner, for Pair and MulOff.
 */
static void
gate_ops(const EVBDD_WGT *f, BQD F, const EVBDD_WGT *g, BQD G, unsigned n, bool every, gops_t *o)
{
    const uint64_t len = UINT64_C(1) << n;
    EVBDD_WGT h[MAXV], u[MAXV];
    const bool full = bqd_has_full_support(F);
    /* the draws where not every case is taken */
    const unsigned q1 = (unsigned)rnd_below(n), qc = (n >= 2) ? 1 + (unsigned)rnd_below(n - 1) : 0;
    const uint64_t cm1 = (n >= 2) ? 1 + rnd_below((UINT64_C(1) << qc) - 1) : 0;
    const unsigned a1 = (unsigned)rnd_below(n);
    const unsigned b1 = (n >= 2) ? (a1 + 1 + (unsigned)rnd_below(n - 1)) % n : 0;
    const size_t i1 = rnd_below(NUM(table_gates)), i2 = rnd_below(NUM(ctl_gates));
    const size_t i3 = rnd_below(NUM(either_gates));
    const unsigned pb = (unsigned)rnd_below(n), ps = (unsigned)rnd_below(4);
    const uint64_t A1 = 1 + rnd_below(len - 1);

    /* one qubit: the table gates and the dynamic one, restriction, projection, X */
    for (unsigned q = 0; q < n; q++) {
        if (!every && q != q1) continue;
        for (size_t i = 0; i < NUM(table_gates); i++) {
            if (!every && i != i1) continue;
            one_gate(f, F, n, table_gates[i], 0, q, &o->t[G_GATE], o);
        }
        one_gate(f, F, n, GATEID_dynamic, 0, q, &o->t[G_DYN], o);
        for (int b = 0; b < 2; b++) {
            dense_cofactor(f, n, q, b, h);
            tally(&o->t[G_RESTRICT], OWN(o, bqd_restrict(F, q, b)), h, n);
            dense_project(f, n, q, b, h);
            tally(&o->t[G_PROJECT], OWN(o, bqd_project(F, q, b)), h, n);
        }
        dense_x(f, n, q, h);
        tally(&o->t[G_X], OWN(o, bqd_x(F, q)), h, n);
        EVBDD_WGT M[4];
        for (int i = 0; i < 4; i++) M[i] = rand_entry();
        const uint32_t qq = q;
        dense_matvec(f, n, M, &qq, 1, h);
        tally(&o->t[G_MATVEC1], OWN(o, bqd_local_matvec(F, M, &qq, 1, n)), h, n);
    }

    /* controls above a target, and phases with controls anywhere */
    for (unsigned q = 1; q < n; q++) for (uint64_t cm = 1; cm < (UINT64_C(1) << q); cm++) {
        if (!every && (q != qc || cm != cm1)) continue;
        for (size_t i = 0; i < NUM(ctl_gates); i++) {
            if (!every && i != i2) continue;
            one_gate(f, F, n, ctl_gates[i], cm, q, &o->t[G_CGATE], o);
        }
    }
    for (unsigned q = 0; q < n && n >= 2; q++) {
        if (!every && q != q1) continue;
        uint64_t pm = 0;
        while (pm == 0) pm = rnd_below(len) & ~(UINT64_C(1) << q);
        one_gate(f, F, n, phase_gates[rnd_below(NUM(phase_gates))], pm, q, &o->t[G_CPHASE], o);
    }

    /* pairs of qubits: cgate_either both ways, the swap, Perm, the local matvec */
    for (unsigned a = 0; a < n; a++) for (unsigned b = 0; b < n; b++) {
        if (a == b || (!every && (a != a1 || b != b1))) continue;
        for (size_t i = 0; i < NUM(either_gates); i++) {
            if (!every && i != i3) continue;
            const uint32_t gid = either_gates[i];
            memcpy(h, f, len * sizeof(EVBDD_WGT));
            dense_gate(h, n, gates[gid], UINT64_C(1) << a, b);
            bool ok;
            const BQD r = OWN(o, bqd_cgate_either(F, gid, a, b, n, &ok));
            if (a > b) o->v[V_EITHER_BELOW]++;
            if (!ok) {
                /* only the Hadamard with the control below is refused, and leaves e */
                o->v[V_REFUSED]++;
                o->t[G_EITHER].tot++;
                if (gid != GATEID_H || a < b || r != F) o->t[G_EITHER].bad++;
                continue;
            }
            tally(&o->t[G_EITHER], r, h, n);
        }
        if (a > b) continue;
        dense_perm(f, n, BQD_PERM_SWAP, a, b, h);
        tally(&o->t[G_SWAP], OWN(o, bqd_swap(F, b, a, n)), h, n);
        for (uint32_t kind = BQD_PERM_SWAP; kind <= BQD_PERM_CX_UP; kind++) {
            dense_perm(f, n, kind, a, b, h);
            tally(&o->t[G_PERM], OWN(o, bqd_perm(F, kind, a, b)), h, n);
        }
        const bool flip = rnd_below(2) != 0;                 /* the two qubits either way */
        const uint32_t qs[2] = { flip ? b : a, flip ? a : b };
        EVBDD_WGT M[16];
        for (int i = 0; i < 16; i++) M[i] = rand_entry();
        dense_matvec(f, n, M, qs, 2, h);
        tally(&o->t[G_MATVEC2], OWN(o, bqd_local_matvec(F, M, qs, 2, n)), h, n);
    }

    /* Pair at every qubit, both restrictions, with the partner and with zero */
    for (unsigned b = 0; b < n; b++) for (int s = 0; s < 2; s++) for (int t = 0; t < 2; t++) {
        if (!every && (b != pb || (unsigned)(2 * s + t) != ps)) continue;
        const int z = (int)rnd_below(6);                 /* now and then a zero operand */
        const EVBDD_WGT *A = (z == 0) ? u : f, *B = (z == 1) ? u : g;
        const BQD EA = (z == 0) ? limdd_zero_edge() : F, EB = (z == 1) ? limdd_zero_edge() : G;
        if (z <= 1) zero_fill(u, len);
        if (z <= 1) o->v[V_PAIR_ZERO]++;
        const unsigned la = limdd_edge_is_zero(EA) ? n : level_of(limdd_target(EA), n);
        const unsigned lb = limdd_edge_is_zero(EB) ? n : level_of(limdd_target(EB), n);
        const unsigned l = la < lb ? la : lb;
        if (l < b) {
            const bool sa = !limdd_edge_is_zero(EA) && s_at(limdd_target(EA), l, n);
            const bool sb = !limdd_edge_is_zero(EB) && s_at(limdd_target(EB), l, n);
            o->v[(sa || sb) ? V_PAIR_S : V_PAIR_NESTED]++;
        }
        dense_pair(A, s, B, t, b, n, h);
        tally(&o->t[G_PAIR], OWN(o, bqd_pair(EA, s, EB, t, b)), h, n);
    }

    /* the diagonal gates on every set of qubits, or one, with a phase drawn */
    for (uint64_t Aq = 1; Aq < len; Aq++) {
        if (!every && Aq != A1) continue;
        const uint64_t A = imask(Aq, n);
        const EVBDD_WGT phase = pw[rnd_below(8)];
        dense_diag(f, n, A, phase, h);
        uint32_t visits = 7;
        const BQD d = OWN(o, bqd_apply_diagonal_counted(F, A, phase, n, &visits));
        tally(&o->t[G_DIAG], d, h, n);
        tally(&o->t[G_DIAG_ANY], OWN(o, bqd_diag_any(F, A, phase, n)), h, n);
        o->t[G_VISITS].tot++;
        if (visits != (full ? diag_path(F, Aq, n) : 0)) o->t[G_VISITS].bad++;
        if (!limdd_edge_is_zero(F)) {
            const unsigned a = (unsigned)__builtin_ctzll(Aq), t = level_of(limdd_target(F), n);
            if (full) o->v[V_DIAG_WALK]++;
            if (t > a) o->v[V_DIAG_SKIP]++;
            else if (t < a && (limdd_node_flags(limdd_target(F)) & BQD_SM_S)) o->v[V_DIAG_S]++;
        }
        if (full) {
            o->t[G_DIAG_MEMO].tot++;
            if (OWN(o, bqd_sm_diag(F, Aq, phase)) != d) o->t[G_DIAG_MEMO].bad++;
        }
        for (uint64_t x = 0; x < len; x++) u[x] = EVBDD_ONE;
        dense_diag(u, n, A, phase, h);
        tally(&o->t[G_MONO], OWN(o, bqd_monomial(A, phase, n)), h, n);
    }

    /* MulOff on the partner's support, on none, and on all */
    {
        const EVBDD_WGT c = pw[1 + rnd_below(7)];
        const int k = (int)rnd_below(5);
        const LIMDD_TARG U = (k == 0) ? 0 : (k == 1) ? LIMDD_TERMINAL : node_or_zero(G);
        for (uint64_t x = 0; x < len; x++) {
            const bool in = (k == 1) || (k > 1 && g[x] != EVBDD_ZERO);
            h[x] = (in || f[x] == EVBDD_ZERO) ? f[x] : wgt_mul(c, f[x]);
        }
        if (k > 1 && !limdd_edge_is_zero(F) && U != 0 && !bqd_sm_full(U)
            && !bqd_sm_sub(limdd_target(F), U)) o->v[V_MULOFF_PHASE]++;
        tally(&o->t[G_MULOFF], OWN(o, bqd_mul_off(F, c, U)), h, n);
    }

    /* PhaseMul and Exp with a random exponent, a root of unity or 1 + sqrt2 */
    {
        int64_t v[MAXV];
        const bool root = rnd_below(4) != 0;
        const EVBDD_WGT beta = root ? pw[1 + rnd_below(7)] : qisq2_lookup(1, 1, 1, 1, 0, 1, 0, 1);
        const uint32_t r = bqd_exp_order(beta);
        rand_exponents(n, rand_qubits(n), root ? 0 : -2, root ? 8 : 5, v);
        const BQD_EXP eps = mtbdd_refs_push(exp_table(v, 0, n, r));
        for (uint64_t x = 0; x < len; x++)
            h[x] = (f[x] == EVBDD_ZERO) ? f[x] : wgt_mul(f[x], bqd_exp_power(beta, v[x]));
        if (!limdd_edge_is_zero(F) && limdd_target(F) != LIMDD_TERMINAL && !bqd_exp_is_const(eps)) {
            const unsigned lN = level_of(limdd_target(F), n), le = bqd_exp_var(eps);
            if (s_at(limdd_target(F), lN < le ? lN : le, n)) o->v[V_PHASE_S]++;
        }
        tally(&o->t[G_PHASEMUL], OWN(o, bqd_phase_mul(F, eps, beta)), h, n);
        for (uint64_t x = 0; x < len; x++) h[x] = bqd_exp_power(beta, v[x]);
        tally(&o->t[G_EXP], OWN(o, bqd_exp_state(eps, beta, r)), h, n);
        mtbdd_refs_pop(1);
    }
}

static void
gops_report(const char *section, const gops_t *o)
{
    char what[112];
    for (int i = 0; i < G_COUNT; i++) {
        snprintf(what, sizeof(what), "%s: %s", section, g_name[i]);
        report(what, o->t[i].bad, o->t[i].tot);
    }
    for (int i = 0; i < V_COUNT; i++) {
        snprintf(what, sizeof(what), "%s: %s", section, v_name[i]);
        covered(what, o->v[i]);
    }
}

/** A random dynamic gate: general, diagonal or a flip, each drawn. */
static void
rand_dynamic(void)
{
    EVBDD_WGT U[4];
    for (int i = 0; i < 4; i++) U[i] = rand_entry();
    switch (rnd_below(4)) {
    case 0: U[1] = U[2] = EVBDD_ZERO; break;
    case 1: U[0] = U[3] = EVBDD_ZERO; break;
    default: break;
    }
    set_dynamic(U);
}

/**
 * Every vector over {0, 1, -1, i} on 1 and 2 qubits with every gate and
 * selection at every qubit, pair and set of controls; and one vector in
 * eight on 3 qubits with one of each, drawn. The partner of each is drawn
 * from the same set.
 */
static void
check_exhaustive_gates(void)
{
    const EVBDD_WGT vals[4] = { EVBDD_ZERO, EVBDD_ONE, EVBDD_MIN_ONE, pw[2] };
    EVBDD_WGT f[8], g[8];
    gops_t *o = calloc(1, sizeof(gops_t));
    structure_begin();
    rand_dynamic();
    for (unsigned n = 1; n <= 3; n++) {
        const uint64_t len = UINT64_C(1) << n, count = UINT64_C(1) << (2 * len);
        for (uint64_t i = (n == 3) ? rnd_below(8) : 0; i < count; i += (n == 3) ? 8 * THIN : 1) {
            for (uint64_t x = 0; x < len; x++) f[x] = vals[(i >> (2 * x)) & 3];
            const uint64_t j = rnd_below(count);
            for (uint64_t x = 0; x < len; x++) g[x] = vals[(j >> (2 * x)) & 3];
            if ((i & 1023) < 8) rand_dynamic();
            const BQD F = limdd_refs_push(bqd_from_vector(f, n));
            const BQD G = limdd_refs_push(bqd_from_vector(g, n));
            gate_ops(f, F, g, G, n, n <= 2, o);
            limdd_refs_pop(2);
        }
    }
    gops_report("n <= 3", o);
    structure_report("n <= 3, gates");
    free(o);
}

/** Random states on 4 to 8 qubits, of every kind, with every gate and selection, drawn. */
static void
check_random_gates(void)
{
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *g = malloc(MAXV * sizeof(EVBDD_WGT));
    gops_t *o = calloc(1, sizeof(gops_t));
    structure_begin();
    for (unsigned n = 4; n <= NQ; n++) for (int rep = 0; rep < 48 / THIN; rep++) {
        rand_dynamic();
        rand_state(n, (unsigned)rnd_below(S_KINDS), f);
        rand_state(n, (unsigned)rnd_below(S_KINDS), g);
        const BQD F = limdd_refs_push(bqd_from_vector(f, n));
        const BQD G = limdd_refs_push(bqd_from_vector(g, n));
        gate_ops(f, F, g, G, n, false, o);
        gate_ops(f, F, g, G, n, false, o);
        limdd_refs_pop(2);
    }
    gops_report("n = 4..8, random", o);
    structure_report("n = 4..8, gates");
    free(o); free(f); free(g);
}

/* --- circuits ---------------------------------------------------------------- */

/**
 * One gate of the qasm runner's set on state e, as the runner calls it, and
 * on the dense vector v: h, x, y, z, s, sdg, t, tdg; cx, cz and swap with the
 * qubits either way; cy and ch with the control above; cs, csdg and ccz.
 */
static BQD
runner_gate(BQD e, EVBDD_WGT *v, unsigned n)
{
    static const uint32_t one[] = { GATEID_H, GATEID_H, GATEID_X, GATEID_Y, GATEID_Z, GATEID_S,
                                    GATEID_Sdag, GATEID_T, GATEID_Tdag };
    const unsigned a = (unsigned)rnd_below(n);
    const unsigned k = (n >= 2) ? (unsigned)rnd_below(6) : 0;
    if (k <= 2) {
        const uint32_t gid = one[rnd_below(NUM(one))];
        dense_gate(v, n, gates[gid], 0, a);
        return bqd_gate(e, gid, a, n);
    }
    unsigned b = (unsigned)rnd_below(n - 1); if (b >= a) b++;
    if (k == 3) {                                        /* cx, cz, swap: either way */
        const unsigned w = (unsigned)rnd_below(3);
        if (w == 2) {
            EVBDD_WGT h[MAXV];
            dense_perm(v, n, BQD_PERM_SWAP, a < b ? a : b, a < b ? b : a, h);
            memcpy(v, h, (UINT64_C(1) << n) * sizeof(EVBDD_WGT));
            return bqd_swap(e, a, b, n);
        }
        const uint32_t gid = w ? GATEID_Z : GATEID_X;
        dense_gate(v, n, gates[gid], UINT64_C(1) << a, b);
        bool ok;
        const BQD r = bqd_cgate_either(e, gid, a, b, n, &ok);
        if (!ok) failures++;
        return r;
    }
    const unsigned lo = a < b ? a : b, hi = a < b ? b : a;
    if (k == 4) {                                        /* cy, ch: the control above */
        const uint32_t gid = rnd_below(2) ? GATEID_Y : GATEID_H;
        dense_gate(v, n, gates[gid], UINT64_C(1) << lo, hi);
        bool ok;
        return bqd_cgate_either(e, gid, lo, hi, n, &ok);
    }
    if (n < 3 || rnd_below(2)) {                         /* cs, csdg */
        const uint32_t gid = rnd_below(2) ? GATEID_S : GATEID_Sdag;
        dense_gate(v, n, gates[gid], UINT64_C(1) << lo, hi);
        return bqd_cgate(e, gid, UINT64_C(1) << lo, hi, n);
    }
    unsigned c = (unsigned)rnd_below(n);
    while (c == a || c == b) c = (unsigned)rnd_below(n);
    uint32_t q[3] = { a, b, c };                         /* ccz: the highest-numbered the target */
    for (int i = 0; i < 2; i++) for (int j = i + 1; j < 3; j++)
        if (q[j] < q[i]) { const uint32_t x = q[i]; q[i] = q[j]; q[j] = x; }
    const uint64_t cm = (UINT64_C(1) << q[0]) | (UINT64_C(1) << q[1]);
    dense_gate(v, n, gates[GATEID_Z], cm, q[2]);
    return bqd_cgate(e, GATEID_Z, cm, q[2], n);
}


/**
 * Random circuits of the runner's gate set from |0...0> on 2 to 8 qubits, as
 * the runner calls them, the state the edge bqd_from_vector builds for the
 * dense one after every gate, and its structure every fourth. A Hadamard
 * after a T gate makes supports that are not affine, where SM has S nodes.
 */
static void
check_circuits(void)
{
    EVBDD_WGT *v = malloc(MAXV * sizeof(EVBDD_WGT));
    tally_t t = { 0, 0 }, refs = { 0, 0 };
    unsigned with_s = 0;
    structure_begin();
    for (unsigned n = 2; n <= NQ; n++) for (int rep = 0; rep < 8 / THIN; rep++) {
        for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) v[x] = x ? EVBDD_ZERO : EVBDD_ONE;
        BQD e = bqd_basis_state(0, n);
        limdd_refs_pushptr(&e);
        const unsigned s_before = sf.s_nodes;
        for (unsigned i = 0; i < 10 * n; i++) {
            size_t d0[3], d1[3];
            limdd_refs_depths(d0);
            e = runner_gate(e, v, n);
            limdd_refs_depths(d1);
            refs.tot++;
            if (memcmp(d0, d1, sizeof(d0)) != 0) refs.bad++;
            t.tot++;
            if (e != bqd_from_vector(v, n)) t.bad++;
            if (i % 4 == 3) check_diagram(e, n);
        }
        if (sf.s_nodes > s_before) with_s++;
        limdd_refs_popptr(1);
    }
    report("circuits: the state after every gate", t.bad, t.tot);
    report("circuits: the reference stacks as they were", refs.bad, refs.tot);
    covered("circuits whose states had new S nodes", with_s);
    structure_report("circuits");
    free(v);
}

/* --- the two rules ---------------------------------------------------------- */

/** Whether a node reachable from p has a flag set. */
static bool
any_flag(LIMDD_TARG p)
{
    if (p == LIMDD_TERMINAL || !set_add(walked, p)) return false;
    if (limdd_node_flags(p) != 0) return true;
    const LIMDD lo = limdd_node_low(p), hi = limdd_node_high(p);
    return (!limdd_edge_is_zero(lo) && any_flag(limdd_target(lo)))
        || (!limdd_edge_is_zero(hi) && any_flag(limdd_target(hi)));
}

#define RULE_VECTORS 240

typedef struct {
    unsigned  n[RULE_VECTORS];
    bool      aligned[RULE_VECTORS];   /* full support or an affine one */
    size_t    counts[RULE_VECTORS][NQ];
    EVBDD_WGT v[RULE_VECTORS][MAXV];
} rule_set_t;

/** Build every vector of the set in the session there is, and count its nodes per level. */
static void
rule_counts(rule_set_t *rs, size_t counts[RULE_VECTORS][NQ], unsigned *flagged)
{
    *flagged = 0;
    for (unsigned i = 0; i < RULE_VECTORS; i++) {
        const BQD e = bqd_from_vector(rs->v[i], rs->n[i]);
        limdd_level_counts(e, counts[i], rs->n[i]);
        set_clear(walked);
        if (!limdd_edge_is_zero(e) && any_flag(limdd_target(e))) (*flagged)++;
    }
}

/*
 * The same vectors in a session of each rule: the copy rule's, which
 * bqd_init_rule gives on request and which makes no node with a flag, and
 * SM's, which bqd_init gives the scalar family. On full support and on affine
 * supports the two diagrams are one up to the tags (skip:cor:smd0), so the
 * counts per level agree; elsewhere they differ, and SM's may be the larger
 * (skip:sec:sm:checks), so those are printed. bqd_init gives the other two
 * families the copy rule.
 */
TASK_0(int, run_rules)
{
    const int before = failures;
    rule_set_t *rs = malloc(sizeof(rule_set_t));
    static size_t sm_counts[RULE_VECTORS][NQ];
    rng_state = UINT64_C(0x5EED5A1E0000);

    bqd_init_rule(BQD_FAMILY_SCALAR, BQD_ZERO_COPY, NQ, 1LL << 20, 1LL << 20, 1LL << 22,
                  1LL << 20);
    exact_constants();
    expect(bqd_zero_rule() == BQD_ZERO_COPY && !bqd_sm()
           && strcmp(bqd_zero_rule_name(bqd_zero_rule()), "copy") == 0,
           "the two rules", "bqd_init_rule with the copy rule does not give it");
    static const unsigned kinds[] = { S_ALGEBRAIC, S_PHASE, S_COSET, S_COSET, S_CIRCUIT, S_ZEROS };
    for (unsigned i = 0; i < RULE_VECTORS; i++) {
        const unsigned n = 2 + (unsigned)rnd_below(NQ - 1);
        const uint64_t len = UINT64_C(1) << n;
        rs->n[i] = n;
        rand_state(n, kinds[i % (sizeof(kinds) / sizeof(kinds[0]))], rs->v[i]);
        if (i % 7 == 6) product_state(n, false, rs->v[i]);
        bool full = true;
        for (uint64_t x = 0; x < len && full; x++) full = rs->v[i][x] != EVBDD_ZERO;
        /* a coset state, or a full one; the circuit states and those with
         * zeros are aligned only by chance, and not counted on */
        const unsigned k = kinds[i % (sizeof(kinds) / sizeof(kinds[0]))];
        rs->aligned[i] = full || (k == S_COSET && i % 7 != 6);
    }
    unsigned flagged;
    rule_counts(rs, rs->counts, &flagged);
    report("the copy rule: nodes with a flag", flagged, RULE_VECTORS);
    bqd_quit();

    bqd_init(BQD_FAMILY_SCALAR, NQ, 1LL << 20, 1LL << 20, 1LL << 22, 1LL << 20);
    exact_constants();
    expect(bqd_zero_rule() == BQD_ZERO_SM && bqd_sm()
           && strcmp(bqd_zero_rule_name(bqd_zero_rule()), "sm") == 0,
           "the two rules", "bqd_init does not give the scalar family SM");
    rule_counts(rs, sm_counts, &flagged);
    covered("SM: diagrams with a flag", flagged);
    unsigned bad = 0, tot = 0, other = 0;
    size_t copy_nodes = 0, sm_nodes = 0;
    for (unsigned i = 0; i < RULE_VECTORS; i++) {
        size_t c = 0, s = 0;
        for (unsigned q = 0; q < rs->n[i]; q++) { c += rs->counts[i][q]; s += sm_counts[i][q]; }
        if (rs->aligned[i]) {
            tot++;
            if (memcmp(rs->counts[i], sm_counts[i], rs->n[i] * sizeof(size_t)) != 0) bad++;
        } else {
            other++;
            copy_nodes += c;
            sm_nodes += s;
        }
    }
    report("full and affine supports: SM's level counts are the copy rule's", bad, tot);
    printf("  %-60s %zu SM, %zu copy (%u vectors)\n", "other supports: nodes", sm_nodes,
           copy_nodes, other);
    bqd_quit();
    for (int fam = BQD_FAMILY_X; fam <= BQD_FAMILY_PAULI; fam++) {
        bqd_init((bqd_family_t)fam, NQ, 1LL << 16, 1LL << 16, 1LL << 18, 1LL << 16);
        expect(bqd_zero_rule() == BQD_ZERO_COPY && !bqd_sm(), "the two rules",
               "bqd_init does not give the translation or Pauli family the copy rule");
        bqd_quit();
    }
    free(rs);
    return failures != before;
}

/* --- collections ------------------------------------------------------------ */

#define GC_WIDTH (NQ + 2)
#define POOL     8

/** Kept results on n qubits: protected edges beside their dense vectors. */
typedef struct {
    unsigned  n;
    BQD       e[POOL];
    EVBDD_WGT v[POOL][MAXV];
} pool_t;

static void
pool_fresh(pool_t *p, unsigned i)
{
    rand_state(p->n, (unsigned)rnd_below(S_KINDS), p->v[i]);
    p->e[i] = bqd_from_vector(p->v[i], p->n);
}

enum { C_MUL, C_ADD, C_RATIO, C_COMPOSE, C_COFACTOR, C_PRODUCT, C_GATE, C_SELECT, C_PERM,
       C_PAIR, C_DIAG, C_PHASE, C_KINDS };

/**
 * Results kept across collections, in two pools that share the tables, on NQ
 * and on NQ - 2 qubits. A random operation takes its operands from one pool
 * and puts its result back, and every one to four operations a collection
 * runs. After it every kept edge must still decode to its vector and still be
 * the edge bqd_from_vector builds, and the set of checked nodes starts empty,
 * since the indices it holds may have been handed out again.
 */
TASK_0(int, run_collections)
{
    const int before = failures;
    bqd_init_rule(BQD_FAMILY_SCALAR, BQD_ZERO_SM, GC_WIDTH, 1LL << 16, 1LL << 16, 1LL << 17,
                  1LL << 16);
    exact_constants();
    rng_state = UINT64_C(0x5EED6C5A0000);

    pool_t *pools = malloc(2 * sizeof(pool_t));
    EVBDD_WGT *h = malloc(MAXV * sizeof(EVBDD_WGT)), *g = malloc(MAXV * sizeof(EVBDD_WGT));
    tally_t t = { 0, 0 }, kept = { 0, 0 };
    unsigned collections = 0, freed = 0;
    pools[0].n = NQ;
    pools[1].n = NQ - 2;
    for (int k = 0; k < 2; k++) for (unsigned i = 0; i < POOL; i++) {
        pool_fresh(&pools[k], i);
        limdd_protect(&pools[k].e[i]);
    }
    structure_begin();

    unsigned next = 1 + (unsigned)rnd_below(4);
    for (int it = 0; it < 400 / THIN; it++) {
        pool_t *p = &pools[rnd_below(3) == 0];
        const unsigned n = p->n;
        const uint64_t len = UINT64_C(1) << n, b0 = ibit(0, n);
        const unsigned a = (unsigned)rnd_below(POOL), b = (unsigned)rnd_below(POOL);
        const EVBDD_WGT *fa = p->v[a], *fb = p->v[b];
        BQD r;
        switch (rnd_below(C_KINDS)) {
        case C_MUL:
            for (uint64_t x = 0; x < len; x++) h[x] = wgt_mul(fa[x], fb[x]);
            r = bqd_multiply(p->e[a], p->e[b]);
            break;
        case C_ADD:
            for (uint64_t x = 0; x < len; x++) h[x] = wgt_add(fa[x], fb[x]);
            r = bqd_add(p->e[a], p->e[b]);
            break;
        case C_RATIO: {                              /* f . g over f, which is g on supp f */
            for (uint64_t x = 0; x < len; x++) h[x] = (fa[x] == EVBDD_ZERO) ? EVBDD_ZERO : fb[x];
            const BQD m = limdd_refs_push(bqd_multiply(p->e[a], p->e[b]));
            r = bqd_sm_apply(BQD_SM_RATIO, p->e[a], m);
            limdd_refs_pop(1);
            break;
        }
        case C_COMPOSE: {
            for (uint64_t x = 0; x < len; x++) h[x] = (x & b0) ? fb[x] : fa[x];
            const BQD lo = limdd_refs_push(bqd_cofactor(p->e[a], 0, 0));
            const BQD hi = limdd_refs_push(bqd_cofactor(p->e[b], 0, 1));
            r = bqd_compose(0, lo, hi);
            limdd_refs_pop(2);
            break;
        }
        case C_COFACTOR: {
            const BQD e = p->e[a];
            const unsigned top = limdd_edge_is_zero(e) ? 0 : level_of(limdd_target(e), n);
            const unsigned v = (unsigned)rnd_below((top < n ? top : n - 1) + 1);
            const int c = (int)rnd_below(2);
            dense_cofactor(fa, n, v, c, h);
            r = bqd_cofactor(e, v, c);
            break;
        }
        case C_GATE:                                 /* a gate of the runner's set */
            memcpy(h, fa, len * sizeof(EVBDD_WGT));
            r = runner_gate(p->e[a], h, n);
            break;
        case C_SELECT: {                             /* restriction, projection or X */
            const unsigned q = (unsigned)rnd_below(n);
            const int c = (int)rnd_below(2);
            switch (rnd_below(3)) {
            case 0:  dense_cofactor(fa, n, q, c, h); r = bqd_restrict(p->e[a], q, c); break;
            case 1:  dense_project(fa, n, q, c, h);  r = bqd_project(p->e[a], q, c);  break;
            default: dense_x(fa, n, q, h);           r = bqd_x(p->e[a], q);           break;
            }
            break;
        }
        case C_PERM: {
            const unsigned qa = (unsigned)rnd_below(n - 1);
            const unsigned qb = qa + 1 + (unsigned)rnd_below(n - 1 - qa);
            const uint32_t kind = (uint32_t)rnd_below(3);
            dense_perm(fa, n, kind, qa, qb, h);
            r = bqd_perm(p->e[a], kind, qa, qb);
            break;
        }
        case C_PAIR: {
            const unsigned q = (unsigned)rnd_below(n);
            const int s0 = (int)rnd_below(2), t0 = (int)rnd_below(2);
            dense_pair(fa, s0, fb, t0, q, n, h);
            r = bqd_pair(p->e[a], s0, p->e[b], t0, q);
            break;
        }
        case C_DIAG: {
            const uint64_t A = imask(1 + rnd_below((UINT64_C(1) << n) - 1), n);
            const EVBDD_WGT phase = pw[1 + rnd_below(7)];
            dense_diag(fa, n, A, phase, h);
            r = bqd_apply_diagonal(p->e[a], A, phase, n);
            break;
        }
        case C_PHASE: {                              /* PhaseMul, or MulOff on b's support */
            const EVBDD_WGT beta = pw[1 + rnd_below(7)];
            if (rnd_below(2)) {
                const LIMDD_TARG U = node_or_zero(p->e[b]);
                for (uint64_t x = 0; x < len; x++)
                    h[x] = (fb[x] != EVBDD_ZERO || fa[x] == EVBDD_ZERO) ? fa[x] : wgt_mul(beta, fa[x]);
                r = bqd_mul_off(p->e[a], beta, U);
                break;
            }
            int64_t v[MAXV];
            rand_exponents(n, rand_qubits(n), 0, 8, v);
            const BQD_EXP eps = mtbdd_refs_push(exp_table(v, 0, n, bqd_exp_order(beta)));
            for (uint64_t x = 0; x < len; x++)
                h[x] = (fa[x] == EVBDD_ZERO) ? fa[x] : wgt_mul(fa[x], bqd_exp_power(beta, v[x]));
            r = bqd_phase_mul(p->e[a], eps, beta);
            mtbdd_refs_pop(1);
            break;
        }
        default:
            for (uint64_t x = 0; x < len; x++) h[x] = wgt_mul(fa[x], fb[x]);
            r = bqd_product(p->e[a], p->e[b]);
            break;
        }
        tally(&t, r, h, n);
        const unsigned d = (unsigned)rnd_below(POOL);
        p->e[d] = r;
        memcpy(p->v[d], h, len * sizeof(EVBDD_WGT));
        /* products and cofactors thin a pool out, and gates on gates grow the
         * exact weights, so now and then a fresh state */
        if (rnd_below(3) == 0) pool_fresh(p, (unsigned)rnd_below(POOL));

        if (--next > 0) continue;
        next = 1 + (unsigned)rnd_below(4);
        const size_t was = limdd_node_table_count();
        CALL(limdd_gc);
        collections++;
        if (limdd_node_table_count() < was) freed++;
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
    structure_report("collections");
    covered("collections", collections);
    covered("collections that freed nodes", freed);
    free(pools); free(h); free(g);
    bqd_quit();
    return failures != before;
}

/* --- harness -------------------------------------------------------------- */

TASK_0(int, run_sm)
{
    const int before = failures;
    bqd_init_rule(BQD_FAMILY_SCALAR, BQD_ZERO_SM, NQ, 1LL << 21, 1LL << 20, 1LL << 22, 1LL << 20);
    exact_constants();
    rng_state = UINT64_C(0x5EED5A000000);
    check_exhaustive_builds();
    check_random_builds();
    check_exhaustive_ops();
    check_random_ops();
    check_exhaustive_gates();
    check_random_gates();
    check_circuits();
    printf("  %zu nodes in the table\n", limdd_node_table_count());
    bqd_quit();
    return failures != before;
}

/**
 * bqd_init_rule with SM and the translation or the Pauli family exits with
 * a message before it makes a table, so it runs in a child, outside Lace and
 * Sylvan, and its exit status and stderr are what is checked.
 */
static int
sm_refused(void)
{
    const int before = failures;
    for (int fam = BQD_FAMILY_X; fam <= BQD_FAMILY_PAULI; fam++) {
        int fd[2];
        if (pipe(fd) != 0) { perror("pipe"); return 1; }
        fflush(stdout);
        const pid_t pid = fork();
        if (pid == 0) {
            dup2(fd[1], 2);
            close(fd[0]);
            bqd_init_rule((bqd_family_t)fam, BQD_ZERO_SM, 4, 1 << 10, 1 << 10, 1 << 10, 1 << 10);
            _exit(7);                                /* it went on */
        }
        close(fd[1]);
        char msg[256] = { 0 };
        size_t got = 0;
        ssize_t k;
        while (got < sizeof(msg) - 1 && (k = read(fd[0], msg + got, sizeof(msg) - 1 - got)) > 0)
            got += (size_t)k;
        close(fd[0]);
        int status = 0;
        waitpid(pid, &status, 0);
        const bool ok = WIFEXITED(status) && WEXITSTATUS(status) == 1
                     && strstr(msg, "scalar family") != NULL;
        char what[96];
        snprintf(what, sizeof(what), "bqd_init_rule(%s, SM) exits with a message",
                 bqd_family_name((bqd_family_t)fam));
        report(what, !ok, 1);
    }
    return failures != before;
}

/** A fresh Sylvan package, with 2^cache memo entries. */
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
    const unsigned workers = getenv("BQD_SM_WORKERS") ? (unsigned)atoi(getenv("BQD_SM_WORKERS")) : 4;
    walked = calloc(1, sizeof(nodeset_t));
    checked = calloc(1, sizeof(nodeset_t));
    if (walked == NULL || checked == NULL) { fprintf(stderr, "out of memory\n"); return 1; }
    int bad = 0, res;

    printf("== BQD, rule SM refused to the translation and Pauli families ==\n");
    res = sm_refused();
    printf("  %s\n", res ? "FAILED" : "ok");
    bad |= res;

    printf("== BQD, rule SM, exact weights, %u workers ==\n", workers);
    session_begin(workers, 20);
    res = RUN(run_sm);
    session_end();
    printf("  %s\n", res ? "FAILED" : "ok");
    bad |= res;

    printf("== BQD, the copy rule and SM on the same vectors, %u workers ==\n", workers);
    session_begin(workers, 20);
    res = RUN(run_rules);
    session_end();
    printf("  %s\n", res ? "FAILED" : "ok");
    bad |= res;

    /* a small memo, which a collection clears, and which overwrites entries
     * more often in between */
    printf("== BQD, rule SM, collections between operations, %u workers ==\n", workers);
    session_begin(workers, 16);
    res = RUN(run_collections);
    session_end();
    printf("  %s\n", res ? "FAILED" : "ok");
    bad |= res;

    free(walked); free(checked);
    return bad;
}
