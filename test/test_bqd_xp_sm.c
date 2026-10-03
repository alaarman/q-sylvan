/*
 * Rule SM of the BQD's translation and Pauli families (skip:sec:sm:xp of the
 * note "Binary Quotient Diagrams with Level Skipping", bqd-skip.tex, whose
 * labels are cited as skip:...), the core of qsylvan_bqd_xp_sm.h, checked in
 * both families so that a node with the wrong tag, flag or label, a level
 * kept or skipped wrongly, or a wrong result cannot pass:
 *
 *   build         decode(build(f)) = f by bqd_to_vector and by a decoder
 *                 written here from skip:def:smxp, bqd_eval agrees at every
 *                 index, building twice gives one edge, the full-support test
 *                 is right, and the cofactor at every level down to an edge's
 *                 top is the edge of the dense cofactor: on every vector over
 *                 {0, 1, -1, i} on 1 to 3 qubits, and on 4 to 8 qubits on
 *                 states with zeros that ignore random qubits, Clifford+T
 *                 circuit states, phase, product, coset and near-coset states;
 *                 a basis state is bqd_basis_state's and the monomials are
 *                 bqd_monomial's
 *   structure     every node of every diagram here, built or computed: no
 *                 level whose two cofactors are equal, nor in the Pauli family
 *                 opposite; the tag S exactly where the high cofactor, moved
 *                 by its least point, has a point outside the low support;
 *                 a Q node's ratio zero off the low support; an S node with a
 *                 nonzero high edge; a low edge the identity on a node, never
 *                 zero; a high label normal, with s.t = 0, and with no Z in
 *                 the translation family; the full flag exactly where the
 *                 node's function has no zero, no other flag bit; the node the
 *                 one bqd_from_vector builds for its function, which makes a
 *                 Pauli node's tag that of the translation family's node of
 *                 its function; Cof and Compose give the node back from its
 *                 cofactors; Ind is the node of the indicator; Side at every
 *                 qubit against the support; the levels an edge skips; and
 *                 Dep at the root against what it reaches
 *   tags          the tag at the root of every vector on 1 to 3 qubits, in a
 *                 session of each family: the two agree where the two roots
 *                 are at one level, and where the Pauli root is lower, its
 *                 edge has a Z above it (f_1 = -f_0, which that family skips)
 *   operations    the product (bqd_multiply and bqd_product), the sum, the
 *                 ratio where the supports allow it, Restr, Sub with a
 *                 translation, Compose of two cofactors at level 0, the public
 *                 cofactor, Canon of a labelled edge, the product, sum and
 *                 ratio of labelled operands, and the view at every level
 *                 against its defining identity: on every pair of vectors on 1
 *                 and 2 qubits, a random partner for every vector on 3, and
 *                 random pairs on 4 to 8 qubits, among them pairs whose
 *                 supports nest and pairs with one ratio at the top; each
 *                 result the edge bqd_from_vector builds for the dense result,
 *                 and a diagram that passes the structure checks; the
 *                 reference stacks as they were after each; the cases of
 *                 Apply (both operands viewed at the top with one
 *                 translation, a nonzero one among them, one ratio for the
 *                 sum, an S node) and of Canon (a relabelled Q node and S
 *                 node, and a relabelling that changes the tag, both ways)
 *                 counted and met, the change of tag on 1 to 3 qubits, where
 *                 random states rarely show it; the smallest S node that a
 *                 translation turns into a Q node (skip:lem:smxptag)
 *   selections, gates and phases
 *                 bqd_restrict, bqd_project and bqd_x on every qubit,
 *                 bqd_perm of each kind, bqd_swap, and CX, CY, CZ, CH and CS
 *                 with the control above and below (bqd_cgate_either) on
 *                 every pair, H and a random one-qubit gate (sqrt X and sqrt Y
 *                 among them) on every qubit, a gate with two controls above,
 *                 a phase on a random set of qubits (bqd_apply_diagonal),
 *                 bqd_phase_mul of a random exponent, bqd_exp_state, bqd_pair
 *                 with a random partner, bqd_mul_off in the translation family
 *                 (with the base 1 + sqrt2, which is no root of unity, among
 *                 them) and bqd_local_matvec on one or two qubits; Restrict,
 *                 Project, Pair, Gate and PhaseMul of labelled operands: on
 *                 every vector on 1 and 2 qubits and every third on 3, and on
 *                 random states on 4 to 8 qubits; each result the edge
 *                 bqd_from_vector builds for the dense one, with the structure
 *                 checks, the reference stacks as they were after each, and
 *                 the cases that take a node's view, an S node, the gate's
 *                 pass-through, a permutation above a Q and an S root, a
 *                 nested Pair, a phase on an S root, and an X or a
 *                 permutation that changes the tag of the root, counted and
 *                 met; H and a random gate on the state projected to one
 *                 side of the qubit, which the gate's one-sided case (Side)
 *                 takes
 *   Side          H on a support on one side of the qubit, from a cleared
 *                 memo, on 4 to 8 qubits: on x_q = 0 no high cofactor and
 *                 no Apply miss, which the cofactor recursion made before
 *                 Side was carried over; and bqd_xpsm_diagonal of the zero
 *                 edge
 *   circuits      random Clifford+T circuits on 2 to 8 qubits through the
 *                 runner's entry points, CS, CH, CCZ and CCX among them, the
 *                 state checked after every gate
 *   collections   limdd_gc between operations in small tables, on protected
 *                 results on 8 and 6 qubits, each result checked as above,
 *                 gates, permutations, selections, phases and Pair among the
 *                 operations, and after every collection every kept edge
 *                 decoded and built again
 *
 * Exact weights, in Q(w_8, sqrt2), throughout, on the number of workers in
 * BQD_XP_SM_WORKERS (default 4).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sylvan.h>
#include <sylvan_int.h>

#include "qsylvan_bqd.h"
#include "qsylvan_bqd_ops.h"
#include "qsylvan_bqd_gates.h"
#include "qsylvan_bqd_sm.h"
#include "qsylvan_bqd_xp.h"
#include "qsylvan_bqd_xp_sm.h"
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
    printf("  %-64s %s (%u)\n", what, bad ? "FAILED" : "ok", tot);
}

/** A case the checks must have met, or they say nothing about it. */
static void
covered(const char *what, unsigned seen)
{
    expect(seen > 0, what, "never met");
    printf("  %-64s %u\n", what, seen);
}

static bool pauli_session(void) { return bqd_family() == BQD_FAMILY_PAULI; }
static const char *fam_name(void) { return pauli_session() ? "P" : "X"; }

static inline bool parity(uint64_t v) { return __builtin_parityll(v) != 0; }

/** Qubit q as a vector-index bit: qubit 0 is the most significant of n. */
static inline uint64_t
ibit(unsigned q, unsigned n)
{
    return UINT64_C(1) << (n - 1 - q);
}

/** A set of qubits, bit q for qubit q (a LIM mask), as an index mask. */
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

static inline bool
is_s(LIMDD_TARG p)
{
    return p != LIMDD_TERMINAL && (limdd_node_flags(p) & BQD_SM_S) != 0;
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

/** supp f inside supp g */
static bool
supp_inside(const EVBDD_WGT *f, const EVBDD_WGT *g, uint64_t len)
{
    for (uint64_t x = 0; x < len; x++) if (f[x] != EVBDD_ZERO && g[x] == EVBDD_ZERO) return false;
    return true;
}

static int64_t
lexmin(const EVBDD_WGT *g, uint64_t len)
{
    for (uint64_t y = 0; y < len; y++) if (g[y] != EVBDD_ZERO) return (int64_t)y;
    return -1;
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

static nodeset_t *checked;

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

/* --- labels, as vector-index masks of an n-qubit vector ------------------- */

typedef struct {
    EVBDD_WGT c;
    uint64_t  s, t;
} lab_t;

static lab_t
lab_of(LIMDD_LIM lim, unsigned n)
{
    lab_t l;
    bqd_lim_masks(lim, n, &l.c, &l.s, &l.t);
    return l;
}

/** a . b, which acts as b and then a: moving b's Z past a's X costs (-1)^{s_b . t_a}. */
static lab_t
lab_mul(lab_t a, lab_t b)
{
    lab_t r = { wgt_mul(a.c, b.c), a.s ^ b.s, a.t ^ b.t };
    if (parity(b.s & a.t)) r.c = wgt_neg(r.c);
    return r;
}

/** (c Z^s X^t . g)(y) = c (-1)^{s.y} g(y ^ t), on a vector of length len. */
static void
lab_apply(lab_t l, const EVBDD_WGT *g, uint64_t len, EVBDD_WGT *out)
{
    for (uint64_t y = 0; y < len; y++) {
        EVBDD_WGT v = g[y ^ l.t];
        if (v != EVBDD_ZERO) {
            if (parity(l.s & y)) v = wgt_neg(v);
            v = wgt_mul(l.c, v);
        }
        out[y] = v;
    }
}

/**
 * The labelled edge l . e, normal (no X above its node), as the operations
 * take them, and not canonical in general: the zero edge for e zero.
 */
static BQD
lab_edge(lab_t l, BQD e, unsigned n)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG N = limdd_target(e);
    lab_t p = lab_mul(l, lab_of(limdd_label(e), n));
    p.t &= (UINT64_C(1) << (n - level_of(N, n))) - 1;
    return limdd_bundle(bqd_lim_make(p.c, p.s, p.t, n), N);
}

/** A random label of the family on n qubits, with a scalar, Z (Pauli) and X at every level. */
static lab_t
rand_lab(unsigned n)
{
    const EVBDD_WGT cs[6] = { EVBDD_ONE, EVBDD_MIN_ONE, pw[2], pw[1], pw[7], isq2 };
    lab_t l = { cs[rnd_below(6)], 0, rnd_below(UINT64_C(1) << n) };
    if (pauli_session()) l.s = rnd_below(UINT64_C(1) << n);
    return l;
}

/* --- the vector of a diagram, from skip:def:smxp ---------------------------- */

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
 * The function of node p read at its own level v, 2^(n-v) values, by rule SM
 * with labels: f_0 = ext[A], and f_1 = l_1 . (ext[A] . ext[R]) for a Q node
 * and l_1 . ext[R] for an S node, the plain product, the label acting on the
 * levels below v. Written from the definition, not from the decoder, so that
 * the two check each other.
 */
static void
node_dense(LIMDD_TARG p, unsigned n, EVBDD_WGT *out)
{
    if (p == LIMDD_TERMINAL) { out[0] = EVBDD_ONE; return; }
    const unsigned v = limdd_node_var(p);
    const uint64_t h = UINT64_C(1) << (n - v - 1);
    const LIMDD low = limdd_node_low(p), high = limdd_node_high(p);
    EVBDD_WGT *u = malloc(h * sizeof(EVBDD_WGT));
    if (limdd_edge_is_zero(low)) zero_fill(out, h); else ext_dense(limdd_target(low), v + 1, n, out);
    if (limdd_edge_is_zero(high)) {
        zero_fill(out + h, h);
    } else {
        ext_dense(limdd_target(high), v + 1, n, u);
        if (!is_s(p))
            for (uint64_t y = 0; y < h; y++)
                u[y] = (u[y] == EVBDD_ZERO || out[y] == EVBDD_ZERO) ? EVBDD_ZERO : wgt_mul(out[y], u[y]);
        lab_apply(lab_of(limdd_label(high), n), u, h, out + h);
    }
    free(u);
}

/** The function the edge e denotes read at level 0 of an n-qubit vector. */
static void
edge_dense(BQD e, unsigned n, EVBDD_WGT *out)
{
    const uint64_t len = UINT64_C(1) << n;
    if (limdd_edge_is_zero(e)) { zero_fill(out, len); return; }
    EVBDD_WGT *g = malloc(len * sizeof(EVBDD_WGT));
    ext_dense(limdd_target(e), 0, n, g);
    lab_apply(lab_of(limdd_label(e), n), g, len, out);
    free(g);
}

/* --- the structure of a diagram --------------------------------------------- */

typedef struct {
    unsigned nodes;         /* distinct nodes checked */
    unsigned redundant;     /* a node whose two cofactors are equal, or opposite in the Pauli family */
    unsigned tag;           /* the tag is not the one skip:def:smxp gives */
    unsigned q_ratio;       /* a Q node whose ratio is not zero off its low support */
    unsigned edges;         /* a zero or labelled low edge, or an S node with a zero high edge */
    unsigned label;         /* a high label with an X above its node, s.t odd, or a Z in X */
    unsigned full;          /* the full flag where the function has a zero, or not where it has none */
    unsigned flags;         /* a flag bit that is neither S nor full */
    unsigned not_canon;     /* a node that is not bqd_from_vector's for its function */
    unsigned cof;           /* Cof or Compose disagrees at a node or a skipped level */
    unsigned ind;           /* Ind is not the node of the indicator */
    unsigned dep;           /* Dep disagrees with what a root reaches */
    unsigned side;          /* Side disagrees with the support at some qubit */
    unsigned q_nodes, s_nodes, q_zero_high, full_nodes, skipping, z_skips, t_high;   /* met */
} structure_t;

static structure_t sf;

/**
 * The levels r.. an edge read at r skips, down to its top: there it is both
 * of its cofactors, or at a Z of the Pauli family two opposite ones, and
 * Compose gives it back. Only for an edge that is canonical: a low edge, an
 * S node's high edge and a root.
 */
static void
check_edge(BQD e, unsigned r, unsigned n)
{
    if (limdd_edge_is_zero(e)) return;
    const unsigned v = level_of(limdd_target(e), n);
    if (v <= r) return;
    sf.skipping++;
    const unsigned top = bqd_xp_top(e) < v ? bqd_xp_top(e) : v;
    for (unsigned m = r; m <= top && m < v; m++) {
        const BQD c0 = bqd_cofactor(e, m, 0), c1 = bqd_cofactor(e, m, 1);
        if (m < top && (c0 != e || c1 != e)) { sf.cof++; break; }
        if (m == top) {
            sf.z_skips++;
            if (c1 != bqd_scale(c0, EVBDD_MIN_ONE)) { sf.cof++; break; }
        }
        if (bqd_compose(m, c0, c1) != e) { sf.cof++; break; }
    }
}

static void check_node(LIMDD_TARG p, unsigned n);

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
    if (lz || limdd_label(low) != LIMDD_LIM_IDENTITY || (s && hz)) sf.edges++;
    if (lz) return;                                     /* nothing below makes sense */
    if (!hz) {
        const lab_t l1 = lab_of(limdd_label(high), n);
        const uint64_t blk = UINT64_C(1) << (n - level_of(limdd_target(high), n));
        if (l1.t >= blk || parity(l1.s & l1.t) || (!pauli_session() && l1.s != 0)) sf.label++;
        if (l1.t != 0) sf.t_high++;
    }
    check_edge(low, v + 1, n);
    if (s) check_edge(high, v + 1, n);

    /* its function at its own level, and the two halves */
    const uint64_t len = UINT64_C(1) << n, blk = UINT64_C(1) << (n - v), h = blk >> 1;
    EVBDD_WGT *b = malloc(blk * sizeof(EVBDD_WGT));
    node_dense(p, n, b);
    const EVBDD_WGT *f0 = b, *f1 = b + h;
    bool eq = true, opp = true;
    for (uint64_t y = 0; y < h; y++) {
        if (f1[y] != f0[y]) eq = false;
        if (f1[y] != wgt_neg(f0[y])) opp = false;
    }
    if (eq || (pauli_session() && opp)) sf.redundant++;

    /* the tag: S where X^{t_1} f_1 has a point where f_0 is zero */
    const int64_t m = lexmin(f1, h);
    bool want_s = false;
    if (m >= 0)
        for (uint64_t y = 0; y < h; y++)
            if (f1[y ^ (uint64_t)m] != EVBDD_ZERO && f0[y] == EVBDD_ZERO) want_s = true;
    if (s != want_s) sf.tag++;
    if (s) {
        sf.s_nodes++;
    } else {
        sf.q_nodes++;
        if (hz) {
            sf.q_zero_high++;
        } else {
            EVBDD_WGT *a = malloc(h * sizeof(EVBDD_WGT)), *r = malloc(h * sizeof(EVBDD_WGT));
            ext_dense(limdd_target(low), v + 1, n, a);
            ext_dense(limdd_target(high), v + 1, n, r);
            if (!supp_inside(r, a, h)) sf.q_ratio++;
            free(a); free(r);
        }
    }

    bool no_zero = true;
    for (uint64_t y = 0; y < blk && no_zero; y++) no_zero = b[y] != EVBDD_ZERO;
    if (((flags & BQD_SM_FULL) != 0) != no_zero) sf.full++;
    if (no_zero) sf.full_nodes++;

    /* the node is bqd_from_vector's for its function, extended to n qubits */
    EVBDD_WGT *g = malloc(len * sizeof(EVBDD_WGT));
    EVBDD_WGT *g0 = malloc(len * sizeof(EVBDD_WGT)), *g1 = malloc(len * sizeof(EVBDD_WGT));
    for (uint64_t x = 0; x < len; x++) g[x] = b[x & (blk - 1)];
    if (bqd_from_vector(g, n) != unit(p)) sf.not_canon++;

    /* Side at every qubit: whether the support lies in {x_q = 0}, the only
     * side it can lie on, since 0 is in it; not above the node, where the
     * function does not depend on x_q */
    for (unsigned q = 0; q < n; q++) {
        bool at1 = false;
        for (uint64_t x = 0; x < len; x++) if (g[x] != EVBDD_ZERO && (x & ibit(q, n))) at1 = true;
        if (bqd_xpsm_side(p, q) != !at1) { sf.side++; break; }
    }

    /* its cofactors are the canonical edges of the two halves, and give it back */
    const uint64_t hb = ibit(v, n);
    for (uint64_t x = 0; x < len; x++) { g0[x] = g[x & ~hb]; g1[x] = g[x | hb]; }
    const BQD c0 = bqd_cofactor(unit(p), v, 0), c1 = bqd_cofactor(unit(p), v, 1);
    if (c0 != bqd_from_vector(g0, n) || c1 != bqd_from_vector(g1, n)
        || bqd_compose(v, c0, c1) != unit(p)) sf.cof++;

    /* Ind is the node of the indicator */
    for (uint64_t x = 0; x < len; x++) g[x] = (g[x] != EVBDD_ZERO) ? EVBDD_ONE : EVBDD_ZERO;
    if (unit(bqd_xpsm_ind(p)) != bqd_from_vector(g, n)) sf.ind++;

    free(b); free(g); free(g0); free(g1);
    check_node(limdd_target(low), n);
    if (!hz) check_node(limdd_target(high), n);
}

/**
 * Dep from the diagram alone: a node of variable q, or a Z_q on a high label,
 * on a path from p through nodes above q.
 */
static bool
reaches(LIMDD_TARG p, unsigned q, unsigned n)
{
    if (p == 0 || p == LIMDD_TERMINAL || limdd_node_var(p) > q) return false;
    if (limdd_node_var(p) == q) return true;
    const LIMDD lo = limdd_node_low(p), hi = limdd_node_high(p);
    if (!limdd_edge_is_zero(hi) && (lab_of(limdd_label(hi), n).s & ibit(q, n))) return true;
    if (!limdd_edge_is_zero(lo) && reaches(limdd_target(lo), q, n)) return true;
    return !limdd_edge_is_zero(hi) && reaches(limdd_target(hi), q, n);
}

/** The structure of e, a diagram on n qubits. */
static void
check_diagram(BQD e, unsigned n)
{
    if (limdd_edge_is_zero(e)) return;
    check_edge(e, 0, n);
    check_node(limdd_target(e), n);
    for (unsigned q = 0; q < n; q++)
        if (bqd_xpsm_dep(limdd_target(e), q) != (int)reaches(limdd_target(e), q, n)) { sf.dep++; break; }
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
        { "a node whose two cofactors are equal or opposite", sf.redundant },
        { "a tag that is not skip:def:smxp's", sf.tag },
        { "a Q node whose ratio is not zero off its low support", sf.q_ratio },
        { "a zero or labelled low edge, or an S node with no high edge", sf.edges },
        { "a high label not normal, with s.t odd, or with a Z in X", sf.label },
        { "a full flag that is wrong", sf.full },
        { "a flag bit that is neither S nor full", sf.flags },
        { "a node that is not the one bqd_from_vector builds for its function", sf.not_canon },
        { "Cof or Compose disagrees at a node or a skipped level", sf.cof },
        { "Ind is not the node of the indicator", sf.ind },
        { "Dep disagrees with what a root reaches", sf.dep },
        { "Side disagrees with the support at a qubit", sf.side },
    };
    unsigned bad = 0;
    char buf[160];
    for (size_t i = 0; i < sizeof(k) / sizeof(k[0]); i++) {
        if (k[i].bad == 0) continue;
        snprintf(buf, sizeof(buf), "%u times %s", k[i].bad, k[i].what);
        expect(false, section, buf);
        bad += k[i].bad;
    }
    snprintf(buf, sizeof(buf), "%s %s: structure", fam_name(), section);
    printf("  %-64s %s (%u nodes: %u Q, %u S, %u full; %u skipping edges)\n", buf,
           bad ? "FAILED" : "ok", sf.nodes, sf.q_nodes, sf.s_nodes, sf.full_nodes, sf.skipping);
}

/* --- vectors: qubit q is bit n-1-q of an index ------------------------------- */

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
 * [x in a + C] . w_8^{P(x)} for a random linear code C, a random shift a and
 * a random phase polynomial P: a function whose support is an affine
 * subspace. Near a coset (skip:sec:sm:xp:checks): with `near`, a few points
 * of the coset are dropped or a few outside it added, where the two rules,
 * and the two kinds of node, part.
 */
static void
coset_state(unsigned n, bool near, EVBDD_WGT *f)
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
    if (near) {
        const unsigned flips = 1 + (unsigned)rnd_below(3);
        for (unsigned i = 0; i < flips; i++) in[rnd_below(len)] ^= 1;
    }
    const uint64_t a = rnd_below(len);
    phase_on(n, rand_qubits(n), f);
    for (uint64_t x = 0; x < len; x++) if (!in[x ^ a]) f[x] = EVBDD_ZERO;
    free(in);
}

enum { S_ZEROS, S_ALGEBRAIC, S_PHASE, S_CIRCUIT, S_PRODUCT, S_BASIS_PLUS, S_COSET, S_NEAR, S_KINDS };

/**
 * A state of one of eight kinds: a function of random qubits with zeros, the
 * same over Z[sqrt2, i]/2 without zeros, a phase state of random qubits, a
 * Clifford+T circuit state, a product state, a basis state on random qubits
 * with |+> on the rest, a coset state and a state near a coset.
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
    case S_COSET:     coset_state(n, false, f); break;
    case S_NEAR:      coset_state(n, true, f); break;
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
    unsigned skip, root_skip, root_z, root_t;
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

    /* the cofactor at every level down to the edge's top */
    if (!limdd_edge_is_zero(e)) {
        const unsigned t = bqd_xp_top(e) < n ? bqd_xp_top(e) : n;
        for (unsigned v = 0; v <= t && v < n; v++) for (int c = 0; c < 2; c++) {
            dense_cofactor(f, n, v, c, g);
            if (bqd_cofactor(e, v, c) != bqd_from_vector(g, n)) { b->cof++; v = n; break; }
        }
        if (level_of(limdd_target(e), n) > 0) b->root_skip++;
        if (t < level_of(limdd_target(e), n)) b->root_z++;
        if (lab_of(limdd_label(e), n).t != 0) b->root_t++;
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
    char what[112];
    const struct { const char *name; unsigned bad; } k[] = {
        { "decode(build(f)) = f", b->decode }, { "the decoder of skip:def:smxp", b->dense },
        { "bqd_eval at every index", b->eval }, { "built twice, one edge", b->again },
        { "full support, from the flag", b->full },
        { "the cofactor at every level to the top", b->cof },
    };
    for (size_t i = 0; i < sizeof(k) / sizeof(k[0]); i++) {
        snprintf(what, sizeof(what), "%s %s: %s", fam_name(), section, k[i].name);
        report(what, k[i].bad, b->tot);
    }
    snprintf(what, sizeof(what), "%s %s: diagrams with a skipping edge (new nodes)", fam_name(), section);
    covered(what, b->skip);
    snprintf(what, sizeof(what), "%s %s: root edges that skip", fam_name(), section);
    covered(what, b->root_skip);
    snprintf(what, sizeof(what), "%s %s: root edges with a translation", fam_name(), section);
    covered(what, b->root_t);
    if (pauli_session()) {
        snprintf(what, sizeof(what), "%s %s: root edges with a Z above their node", fam_name(), section);
        covered(what, b->root_z);
    }
}

/*
 * The root of every vector on 1 to 3 qubits, as the X session leaves it for
 * the Pauli session to compare: the level of its node, n for zero or a
 * terminal, and whether the node is an S node.
 */
#define EXHAUSTIVE (16u + 256u + 65536u)
static uint8_t root_level[EXHAUSTIVE], root_s[EXHAUSTIVE];

static unsigned
root_level_of(BQD e, unsigned n)
{
    return limdd_edge_is_zero(e) ? n : level_of(limdd_target(e), n);
}

/**
 * The Pauli root of a vector against the translation family's: one tag where
 * the two nodes are at one level, and where the Pauli node is lower, a Z on
 * the Pauli root above it, the f_1 = -f_0 that family skips; counted in *z.
 */
static bool
root_agrees(BQD e, unsigned n, unsigned i, unsigned *z)
{
    const unsigned lv = root_level_of(e, n);
    if (lv == root_level[i]) return lv == n || (is_s(limdd_target(e)) ? 1 : 0) == root_s[i];
    if (lv < root_level[i] || limdd_edge_is_zero(e) || bqd_xp_top(e) >= lv) return false;
    (*z)++;
    return true;
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
    unsigned i = 0, tags_bad = 0, roots_s = 0, roots_z = 0;
    for (unsigned n = 1; n <= 3; n++) {
        const uint64_t len = UINT64_C(1) << n;
        for (uint64_t idx = 0; idx < (UINT64_C(1) << (2 * len)); idx++, i++) {
            for (uint64_t x = 0; x < len; x++) f[x] = vals[(idx >> (2 * x)) & 3];
            const BQD e = build_checked(f, n, &b);
            const bool s = !limdd_edge_is_zero(e) && is_s(limdd_target(e));
            if (s) roots_s++;
            if (pauli_session()) {
                if (!root_agrees(e, n, i, &roots_z)) tags_bad++;
            } else {
                root_level[i] = (uint8_t)root_level_of(e, n);
                root_s[i] = s;
            }
        }
    }
    builds_report("n <= 3, all", &b);
    structure_report("n <= 3, all");
    char what[112];
    snprintf(what, sizeof(what), "%s n <= 3, all: S nodes", fam_name());
    covered(what, sf.s_nodes);
    snprintf(what, sizeof(what), "%s n <= 3, all: Q nodes with no high cofactor", fam_name());
    covered(what, sf.q_zero_high);
    snprintf(what, sizeof(what), "%s n <= 3, all: high labels with a translation", fam_name());
    covered(what, sf.t_high);
    if (pauli_session()) {
        covered("P n <= 3, all: skipped levels with a Z", sf.z_skips);
        report("P n <= 3, all: the root's tag is the translation family's", tags_bad, EXHAUSTIVE);
        covered("P n <= 3, all: roots with tag S", roots_s);
        covered("P n <= 3, all: roots below the translation family's, under a Z", roots_z);
    }
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

    /* a basis state is bqd_basis_state's, a node per level */
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
    char what[112];
    snprintf(what, sizeof(what), "%s a basis state: a node per level, bqd_basis_state", fam_name());
    report(what, bad, tot);

    /* the monomials, phase^{x_A}, a phase of -1 among them */
    bad = tot = 0;
    const EVBDD_WGT phases[4] = { EVBDD_MIN_ONE, pw[2], pw[1], pw[5] };
    for (unsigned n = 1; n <= NQ; n++) for (int rep = 0; rep < 8; rep++) {
        const uint64_t A = imask(1 + rnd_below((UINT64_C(1) << n) - 1), n);
        const EVBDD_WGT phase = phases[rnd_below(4)];
        for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) f[x] = ((x & A) == A) ? phase : EVBDD_ONE;
        const BQD e = bqd_monomial(A, phase, n);
        check_diagram(e, n);
        tot++;
        if (e != bqd_from_vector(f, n)) bad++;
    }
    snprintf(what, sizeof(what), "%s bqd_monomial", fam_name());
    report(what, bad, tot);
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
    tally_t mul, prod, add, ratio, restr, sub, compose, cof, canon, lmul, ladd, lratio, view, refs;
    unsigned view_pair, view_pair_t, add_one_ratio, s_pair, subset_true, subset_false;
    unsigned canon_q, canon_s, s_to_q, q_to_s;
} ops_t;

/**
 * Whether the views of a and b at the higher of their tops are the shortcut
 * Apply takes for a product or a ratio (one translation) and for a sum (one
 * translation and one ratio), and whether either has an S node there.
 */
static void
apply_cases(BQD a, BQD b, ops_t *o)
{
    if (limdd_edge_is_zero(a) || limdd_edge_is_zero(b)) return;
    const uint32_t ta = bqd_xp_top(a), tb = bqd_xp_top(b), v = ta < tb ? ta : tb;
    if (v >= limdd_level(LIMDD_TERMINAL)) return;
    BQD a0, ra, b0, rb;
    uint64_t sa, sb;
    if (!RUN(bqd_xpsm_view, a, v, &a0, &sa, &ra) || !RUN(bqd_xpsm_view, b, v, &b0, &sb, &rb)) {
        o->s_pair++;
        return;
    }
    if (limdd_edge_is_zero(ra) || limdd_edge_is_zero(rb) || sa == sb) {
        o->view_pair++;
        if (sa != 0 || sb != 0) o->view_pair_t++;
    }
    if (sa == sb && ra == rb && !limdd_edge_is_zero(ra)) o->add_one_ratio++;
}

/**
 * The view of the labelled edge a, whose function is fa, at every level j at
 * or above its top: where it is defined, F_0 = F0, F_1 = X^t (F0 . rho) and
 * supp rho inside supp F0, as functions of all n qubits.
 */
static void
check_views(BQD a, const EVBDD_WGT *fa, unsigned n, ops_t *o)
{
    if (limdd_edge_is_zero(a)) return;
    const uint64_t len = UINT64_C(1) << n;
    const uint32_t top = bqd_xp_top(a);
    EVBDD_WGT *d0 = malloc(len * sizeof(EVBDD_WGT)), *dr = malloc(len * sizeof(EVBDD_WGT));
    for (unsigned j = 0; j < n && j <= top; j++) {
        BQD F0, rho;
        uint64_t t;
        if (!RUN(bqd_xpsm_view, a, j, &F0, &t, &rho)) continue;
        edge_dense(F0, n, d0);
        edge_dense(rho, n, dr);
        const uint64_t T = imask(t, n), hb = ibit(j, n);
        bool ok = (T & ~((hb << 1) - 1)) == 0 && (T & hb) == 0;   /* t is of the levels below j */
        for (uint64_t x = 0; x < len && ok; x++) {
            const EVBDD_WGT want1 = (d0[x ^ T] == EVBDD_ZERO || dr[x ^ T] == EVBDD_ZERO)
                                  ? EVBDD_ZERO : wgt_mul(d0[x ^ T], dr[x ^ T]);
            ok = fa[x & ~hb] == d0[x] && fa[x | hb] == want1 && (dr[x] == EVBDD_ZERO || d0[x] != EVBDD_ZERO);
        }
        o->view.tot++;
        if (!ok) o->view.bad++;
    }
    free(d0); free(dr);
}

/** Canon of a labelled edge whose node is at its top with no X there, and the tag it ends with. */
static void
canon_cases(BQD a, BQD c, unsigned n, ops_t *o)
{
    if (limdd_edge_is_zero(a) || limdd_edge_is_zero(c)) return;
    const LIMDD_TARG N = limdd_target(a), C = limdd_target(c);
    if (N == LIMDD_TERMINAL || bqd_xp_top(a) != limdd_node_var(N)) return;
    if (lab_of(limdd_label(a), n).t & ibit(limdd_node_var(N), n)) return;
    if (is_s(N)) o->canon_s++; else o->canon_q++;
    if (C == LIMDD_TERMINAL || limdd_node_var(C) != limdd_node_var(N)) return;
    if (is_s(N) && !is_s(C)) o->s_to_q++;
    if (!is_s(N) && is_s(C)) o->q_to_s++;
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
    EVBDD_WGT h[MAXV], fa[MAXV], gb[MAXV], ef_d[MAXV], eg_d[MAXV];
    size_t d0[3], d1[3];
    limdd_refs_depths(d0);
    apply_cases(ef, eg, o);

    for (uint64_t x = 0; x < len; x++) h[x] = wgt_mul(f[x], g[x]);
    tally(&o->mul, bqd_multiply(ef, eg), h, n);
    tally(&o->prod, bqd_product(ef, eg), h, n);

    for (uint64_t x = 0; x < len; x++) h[x] = wgt_add(f[x], g[x]);
    tally(&o->add, bqd_add(ef, eg), h, n);

    /* the ratio, where supp g lies inside supp f */
    if (supp_inside(g, f, len)) {
        for (uint64_t x = 0; x < len; x++)
            h[x] = (f[x] == EVBDD_ZERO || g[x] == EVBDD_ZERO) ? EVBDD_ZERO : wgt_div(g[x], f[x]);
        tally(&o->ratio, bqd_xpsm_apply(BQD_SM_RATIO, ef, eg), h, n);
    }

    /* Sub with a random translation, and Restr where g is not zero */
    if (!limdd_edge_is_zero(ef) && !limdd_edge_is_zero(eg)) {
        const LIMDD_TARG X = limdd_target(ef), Y = limdd_target(eg);
        const uint64_t tau = rnd_below(UINT64_C(1) << n);           /* bit q for qubit q */
        const uint64_t T = imask(tau, n);
        edge_dense(unit(X), n, ef_d);
        edge_dense(unit(Y), n, eg_d);
        bool in = true;
        for (uint64_t x = 0; x < len && in; x++) in = ef_d[x] == EVBDD_ZERO || eg_d[x ^ T] != EVBDD_ZERO;
        if (in) o->subset_true++; else o->subset_false++;
        o->sub.tot++;
        if (bqd_xpsm_sub(X, Y, tau) != (int)in) o->sub.bad++;
    }
    if (!limdd_edge_is_zero(eg)) {
        for (uint64_t x = 0; x < len; x++) h[x] = (g[x] == EVBDD_ZERO) ? EVBDD_ZERO : f[x];
        tally(&o->restr, bqd_xpsm_restr(ef, eg), h, n);
    }

    /* Compose at level 0 of a cofactor of each, and the public cofactor at a level to the top */
    const uint64_t b0 = ibit(0, n);
    for (uint64_t x = 0; x < len; x++) h[x] = (x & b0) ? g[x] : f[x];
    tally(&o->compose, bqd_compose(0, bqd_cofactor(ef, 0, 0), bqd_cofactor(eg, 0, 1)), h, n);
    if (!limdd_edge_is_zero(ef)) {
        const unsigned top = bqd_xp_top(ef) < n - 1 ? bqd_xp_top(ef) : n - 1;
        const unsigned v = (unsigned)rnd_below(top + 1);
        const int c = (int)rnd_below(2);
        dense_cofactor(f, n, v, c, h);
        tally(&o->cof, bqd_cofactor(ef, v, c), h, n);
    }

    /* labelled operands: Canon, the views, the product, the sum, the ratio */
    const lab_t la = rand_lab(n), lb = rand_lab(n);
    const BQD A1 = limdd_refs_push(lab_edge(la, ef, n));
    const BQD B1 = limdd_refs_push(lab_edge(lb, eg, n));
    lab_apply(la, f, len, fa);
    lab_apply(lb, g, len, gb);
    const BQD ca = bqd_xpsm_canon(A1);
    tally(&o->canon, ca, fa, n);
    canon_cases(A1, ca, n, o);
    check_views(A1, fa, n, o);
    apply_cases(A1, B1, o);
    for (uint64_t x = 0; x < len; x++) h[x] = wgt_mul(fa[x], gb[x]);
    tally(&o->lmul, bqd_xpsm_apply(BQD_SM_MUL, A1, B1), h, n);
    for (uint64_t x = 0; x < len; x++) h[x] = wgt_add(fa[x], gb[x]);
    tally(&o->ladd, bqd_xpsm_apply(BQD_SM_ADD, A1, B1), h, n);
    for (uint64_t x = 0; x < len; x++) gb[x] = (fa[x] == EVBDD_ZERO) ? EVBDD_ZERO : gb[x];
    const BQD G1 = limdd_refs_push(bqd_from_vector(gb, n));
    for (uint64_t x = 0; x < len; x++)
        h[x] = (fa[x] == EVBDD_ZERO || gb[x] == EVBDD_ZERO) ? EVBDD_ZERO : wgt_div(gb[x], fa[x]);
    tally(&o->lratio, bqd_xpsm_apply(BQD_SM_RATIO, A1, G1), h, n);
    limdd_refs_pop(3);

    limdd_refs_depths(d1);
    o->refs.tot++;
    if (memcmp(d0, d1, sizeof(d0)) != 0) o->refs.bad++;
}

/**
 * The tallies of a section, and the cases it must have met. A relabelling
 * that changes a tag is rare on random states, so it is required of the
 * exhaustive section and only counted elsewhere.
 */
static void
ops_report(const char *section, const ops_t *o, bool exhaustive)
{
    char what[112];
    const struct { const char *name; const tally_t *t; } k[] = {
        { "bqd_multiply", &o->mul }, { "bqd_product", &o->prod }, { "bqd_add", &o->add },
        { "the ratio, supp g inside supp f", &o->ratio }, { "Restr", &o->restr },
        { "Sub with a translation, against the dense supports", &o->sub },
        { "Compose of two cofactors at level 0", &o->compose },
        { "bqd_cofactor at a level to the top", &o->cof },
        { "Canon of a labelled edge", &o->canon },
        { "the product of labelled operands", &o->lmul },
        { "the sum of labelled operands", &o->ladd },
        { "the ratio of labelled operands", &o->lratio },
        { "the view at every level, its identity", &o->view },
        { "the reference stacks as they were", &o->refs },
    };
    for (size_t i = 0; i < sizeof(k) / sizeof(k[0]); i++) {
        snprintf(what, sizeof(what), "%s %s: %s", fam_name(), section, k[i].name);
        report(what, k[i].t->bad, k[i].t->tot);
    }
    const struct { const char *name; unsigned seen; } c[] = {
        { "pairs viewed at the top with one translation", o->view_pair },
        { "  with a nonzero translation", o->view_pair_t },
        { "pairs with one ratio at the top (the sum's shortcut)", o->add_one_ratio },
        { "pairs with an S node or an X at the top", o->s_pair },
        { "supports that nest, and that do not", o->subset_true < o->subset_false
                                                 ? o->subset_true : o->subset_false },
        { "Canon of a relabelled Q node", o->canon_q },
        { "Canon of a relabelled S node", o->canon_s },
        { "a relabelling that turns an S node into a Q node", o->s_to_q },
        { "a relabelling that turns a Q node into an S node", o->q_to_s },
    };
    const size_t required = exhaustive ? sizeof(c) / sizeof(c[0]) : sizeof(c) / sizeof(c[0]) - 2;
    for (size_t i = 0; i < sizeof(c) / sizeof(c[0]); i++) {
        snprintf(what, sizeof(what), "%s %s: %s", fam_name(), section, c[i].name);
        if (i < required) covered(what, c[i].seen);
        else printf("  %-64s %u\n", what, c[i].seen);
    }
}

/**
 * The smallest S node a translation turns into a Q node (skip:lem:smxptag):
 * f = (0,1,1,1,0,0,1,1) has an S root, and X on qubit 1 of it, whose least
 * point is 0, a Q root, which Canon of the label product has to make.
 */
static void
check_tag_change(void)
{
    static const int fv[8] = { 0, 1, 1, 1, 0, 0, 1, 1 };
    EVBDD_WGT f[8], g[8];
    for (int x = 0; x < 8; x++) f[x] = fv[x] ? EVBDD_ONE : EVBDD_ZERO;
    for (int x = 0; x < 8; x++) g[x] = f[x ^ (int)ibit(1, 3)];
    const BQD e = bqd_from_vector(f, 3);
    const lab_t xq = { EVBDD_ONE, 0, ibit(1, 3) };
    const BQD c = bqd_xpsm_canon(lab_edge(xq, e, 3));
    const bool ok = is_s(limdd_target(e)) && limdd_node_var(limdd_target(e)) == 0
                 && !is_s(limdd_target(c)) && limdd_node_var(limdd_target(c)) == 0
                 && c == bqd_from_vector(g, 3);
    char what[112];
    snprintf(what, sizeof(what), "%s an S root that X on qubit 1 turns into a Q root", fam_name());
    report(what, !ok, 1);
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
            for (uint64_t j = 0; j < count; j += (n == 2 ? THIN : 1))
                ops_on_pair(vec[i], vec[j], e[i], e[j], n, &o);
    }
    for (uint64_t i = 0; i < 65536; i += THIN) {
        for (uint64_t x = 0; x < 8; x++) f[x] = vals[(i >> (2 * x)) & 3];
        const uint64_t j = rnd_below(65536);
        for (uint64_t x = 0; x < 8; x++) g[x] = vals[(j >> (2 * x)) & 3];
        ops_on_pair(f, g, bqd_from_vector(f, 3), bqd_from_vector(g, 3), 3, &o);
    }
    ops_report("n <= 2 all pairs, n = 3", &o, true);
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
    ops_report("n = 4..8, random", &o, false);
    structure_report("n = 4..8, operations");
    free(f); free(g); free(r);
}

/* --- selections, gates and phases ------------------------------------------- */

/** pi(x) for the permutation `kind` of qubits qa < qb, on a vector index. */
static uint64_t
perm_index(uint32_t kind, unsigned qa, unsigned qb, uint64_t x, unsigned n)
{
    const uint64_t ba = ibit(qa, n), bb = ibit(qb, n);
    const bool a = (x & ba) != 0, b = (x & bb) != 0;
    bool na = a, nb = b;
    switch (kind) {
    case BQD_PERM_SWAP:    na = b; nb = a; break;
    case BQD_PERM_CX_DOWN: nb = b ^ a;     break;     /* control qa, target qb */
    default:               na = a ^ b;     break;     /* control qb, target qa */
    }
    x &= ~(ba | bb);
    return x | (na ? ba : 0) | (nb ? bb : 0);
}

/** A vector index as a LIM mask, bit q for qubit q, the point an exponent is read at. */
static uint64_t
lim_point(uint64_t x, unsigned n)
{
    uint64_t y = 0;
    for (unsigned q = 0; q < n; q++) if ((x >> (n - 1 - q)) & 1) y |= UINT64_C(1) << q;
    return y;
}

/** The exponent diagram of the table v over LIM points of n qubits, reduced modulo r. */
static BQD_EXP
exp_of(const int64_t *v, unsigned var, unsigned n, uint64_t y, uint32_t r)
{
    if (var == n) return bqd_exp_const(v[y], r);
    const BQD_EXP lo = mtbdd_refs_push(exp_of(v, var + 1, n, y, r));
    const BQD_EXP hi = mtbdd_refs_push(exp_of(v, var + 1, n, y | (UINT64_C(1) << var), r));
    const BQD_EXP e = bqd_exp_node(var, lo, hi);
    mtbdd_refs_pop(2);
    return e;
}

/** A random exponent of n qubits: a function of a random set of them, small values. */
static BQD_EXP
rand_exp(unsigned n, uint32_t r, int64_t *v)
{
    const uint64_t D = rand_qubits(n), len = UINT64_C(1) << n;
    for (uint64_t y = 0; y < len; y++)
        v[y] = (y & ~D) ? v[y & D] : (int64_t)rnd_below(r ? r : 5) - (r ? 0 : 2);
    return exp_of(v, 0, n, 0, r);
}

/** beta^v for an integer v, a negative one included. */
static EVBDD_WGT
wpow(EVBDD_WGT beta, int64_t v)
{
    return bqd_exp_power(beta, v);
}

typedef struct {
    tally_t restr, proj, x, perm, pair, gate, sgate, cgate, either, swap, diag, phase, expst, muloff,
            matvec, lrestr, lproj, lpair, lgate, lphase, refs;
    unsigned nest_sel, s_sel, pass, side, s_gate, perm_q, perm_s, perm_tag, x_tag, phase_s, pair_nest;
} sel_t;

/**
 * Whether e's root is a stored node at its top v with a view there whose
 * translation has no X at any qubit of `avoid`: the case where the
 * selections and the gate take the node's view rather than its cofactors.
 */
static bool
root_viewed(BQD e, uint64_t avoid, BQD *rho)
{
    if (limdd_edge_is_zero(e) || limdd_target(e) == LIMDD_TERMINAL) return false;
    const uint32_t v = bqd_xp_top(e);
    if (v != limdd_node_var(limdd_target(e))) return false;
    BQD F0;
    uint64_t t;
    if (!RUN(bqd_xpsm_view, e, v, &F0, &t, rho)) return false;
    return (t & avoid) == 0;
}

/** Whether the ratio rho does not depend on x_q as the diagram stores it: Gate's pass-through test. */
static bool
free_of(BQD rho, unsigned q, unsigned n)
{
    if (limdd_edge_is_zero(rho)) return true;
    if (lab_of(limdd_label(rho), n).s & ibit(q, n)) return false;
    return !bqd_xpsm_dep(limdd_target(rho), q);
}

static bool
root_is_s(BQD e)
{
    return !limdd_edge_is_zero(e) && limdd_target(e) != LIMDD_TERMINAL && is_s(limdd_target(e));
}

/** The tag of the root against that of e's root at the same variable, where both are nodes there. */
static bool
tag_changed(BQD e, BQD r)
{
    if (limdd_edge_is_zero(e) || limdd_edge_is_zero(r)) return false;
    const LIMDD_TARG a = limdd_target(e), b = limdd_target(r);
    if (a == LIMDD_TERMINAL || b == LIMDD_TERMINAL || limdd_node_var(a) != limdd_node_var(b)) return false;
    return is_s(a) != is_s(b);
}

/** The random gates of the runner's set and two more, as (gate, controls, target) on n qubits. */
static const uint32_t one_gates[] = { GATEID_X, GATEID_Y, GATEID_Z, GATEID_H, GATEID_S, GATEID_Sdag,
                                      GATEID_T, GATEID_Tdag, GATEID_sqrtX, GATEID_sqrtY };
#define N_ONE_GATES (sizeof(one_gates) / sizeof(one_gates[0]))

/**
 * Every selection, gate and phase of rule SM on the vector f, with edge e, on
 * n qubits, against the dense result; g, with edge eg, is a partner for Pair.
 * With `all`, every qubit and every pair of qubits; otherwise one random each.
 */
static void
sel_on(const EVBDD_WGT *f, BQD e, const EVBDD_WGT *g, BQD eg, unsigned n, bool all, sel_t *o)
{
    const uint64_t len = UINT64_C(1) << n;
    EVBDD_WGT *h = malloc(len * sizeof(EVBDD_WGT)), *fa = malloc(len * sizeof(EVBDD_WGT));
    EVBDD_WGT *gb = malloc(len * sizeof(EVBDD_WGT));
    int64_t *ev = malloc(len * sizeof(int64_t));
    size_t d0[3], d1[3];
    limdd_refs_depths(d0);
    limdd_refs_push(e);
    limdd_refs_push(eg);
    const unsigned q0 = all ? 0 : (unsigned)rnd_below(n), q1 = all ? n : q0 + 1;

    /* Restrict, Project and X on every qubit, and the labelled operand of each */
    for (unsigned q = q0; q < q1; q++) {
        BQD rho;
        if (root_viewed(e, UINT64_C(1) << q, &rho) && bqd_xp_top(e) < q) o->nest_sel++;
        if (root_is_s(e) && bqd_xp_top(e) < q) o->s_sel++;
        for (int b = 0; b < 2; b++) {
            dense_cofactor(f, n, q, b, h);
            tally(&o->restr, bqd_restrict(e, q, b), h, n);
            for (uint64_t x = 0; x < len; x++) h[x] = (((x & ibit(q, n)) != 0) == b) ? f[x] : EVBDD_ZERO;
            tally(&o->proj, bqd_project(e, q, b), h, n);
        }
        for (uint64_t x = 0; x < len; x++) h[x] = f[x ^ ibit(q, n)];
        const BQD r = bqd_x(e, q);
        tally(&o->x, r, h, n);
        if (tag_changed(e, r)) o->x_tag++;

        const lab_t la = rand_lab(n);
        const BQD A1 = limdd_refs_push(lab_edge(la, e, n));
        lab_apply(la, f, len, fa);
        const int b = (int)rnd_below(2);
        dense_cofactor(fa, n, q, b, h);
        tally(&o->lrestr, bqd_xpsm_restrict(A1, q, b), h, n);
        for (uint64_t x = 0; x < len; x++) h[x] = (((x & ibit(q, n)) != 0) == b) ? fa[x] : EVBDD_ZERO;
        tally(&o->lproj, bqd_xpsm_project(A1, q, b), h, n);
        limdd_refs_pop(1);
    }

    /* the permutations, the swap and the two-qubit gates on every pair */
    unsigned pa = 0, pb = 0;
    if (!all && n >= 2) {
        pa = (unsigned)rnd_below(n);
        pb = (unsigned)rnd_below(n - 1); if (pb >= pa) pb++;
        if (pb < pa) { const unsigned x = pa; pa = pb; pb = x; }
    }
    for (unsigned qa = 0; qa < n; qa++) for (unsigned qb = qa + 1; qb < n; qb++) {
        if (!all && (qa != pa || qb != pb)) continue;
        if (!limdd_edge_is_zero(e) && bqd_xp_top(e) < qa) {
            if (root_is_s(e)) o->perm_s++; else o->perm_q++;
        }
        for (uint32_t kind = BQD_PERM_SWAP; kind <= BQD_PERM_CX_UP; kind++) {
            for (uint64_t x = 0; x < len; x++) h[x] = f[perm_index(kind, qa, qb, x, n)];
            const BQD r = bqd_perm(e, kind, qa, qb);
            tally(&o->perm, r, h, n);
            if (tag_changed(e, r) && bqd_xp_top(e) < qa) o->perm_tag++;
        }
        for (uint64_t x = 0; x < len; x++) h[x] = f[perm_index(BQD_PERM_SWAP, qa, qb, x, n)];
        tally(&o->swap, (rnd_below(2) ? bqd_swap(e, qa, qb, n) : bqd_swap(e, qb, qa, n)), h, n);

        /* CX, CY, CZ and CH, the control above and below the target, and CS */
        static const uint32_t cg[] = { GATEID_X, GATEID_Y, GATEID_Z, GATEID_H, GATEID_S };
        for (size_t i = 0; i < sizeof(cg) / sizeof(cg[0]); i++) for (int up = 0; up < 2; up++) {
            const unsigned c = up ? qa : qb, t = up ? qb : qa;
            memcpy(h, f, len * sizeof(EVBDD_WGT));
            dense_gate(h, n, gates[cg[i]], UINT64_C(1) << c, t);
            bool ok;
            const BQD r = bqd_cgate_either(e, cg[i], c, t, n, &ok);
            if (cg[i] == GATEID_H && !up) {
                o->either.tot++;
                if (ok) o->either.bad++;             /* no reordering identity for CH */
            } else {
                tally(up ? &o->cgate : &o->either, r, h, n);
            }
        }
    }

    /* one-qubit gates, and a gate with two controls above its target */
    for (unsigned q = q0; q < q1; q++) {
        const uint32_t gid = one_gates[rnd_below(N_ONE_GATES)];
        for (int k = 0; k < 2; k++) {
            const uint32_t id = k ? gid : GATEID_H;
            memcpy(h, f, len * sizeof(EVBDD_WGT));
            dense_gate(h, n, gates[id], 0, q);
            tally(&o->gate, bqd_gate(e, id, q, n), h, n);
        }
        BQD rho;
        if (root_viewed(e, UINT64_C(1) << q, &rho) && bqd_xp_top(e) < q && free_of(rho, q, n))
            o->pass++;
        if (root_is_s(e) && bqd_xp_top(e) < q) o->s_gate++;
        if (q >= 2) {
            const unsigned c1 = (unsigned)rnd_below(q), c2 = (unsigned)rnd_below(q);
            const uint64_t cm = (UINT64_C(1) << c1) | (UINT64_C(1) << c2);
            static const uint32_t ccg[] = { GATEID_H, GATEID_X, GATEID_Z, GATEID_sqrtY, GATEID_T };
            const uint32_t id = ccg[rnd_below(sizeof(ccg) / sizeof(ccg[0]))];
            memcpy(h, f, len * sizeof(EVBDD_WGT));
            dense_gate(h, n, gates[id], cm, q);
            tally(&o->cgate, bqd_cgate(e, id, cm, q, n), h, n);
        }
        /* the gate on a labelled operand, with a control above where there is one */
        const lab_t la = rand_lab(n);
        const BQD A1 = limdd_refs_push(lab_edge(la, e, n));
        lab_apply(la, f, len, fa);
        uint64_t cm = 0;
        if (q > 0 && rnd_below(2)) cm = UINT64_C(1) << rnd_below(q);
        memcpy(h, fa, len * sizeof(EVBDD_WGT));
        dense_gate(h, n, gates[gid], cm, q);
        tally(&o->lgate, bqd_xpsm_gate(A1, gid, cm, q), h, n);
        limdd_refs_pop(1);

        /* H and the random gate on the state projected to one side of q,
         * which the gate's Side takes where the root is above q */
        const int b = (int)rnd_below(2);
        for (uint64_t x = 0; x < len; x++) fa[x] = (((x & ibit(q, n)) != 0) == b) ? f[x] : EVBDD_ZERO;
        const BQD P = limdd_refs_push(bqd_project(e, q, b));
        if (!limdd_edge_is_zero(P) && bqd_xp_top(P) < q && bqd_xpsm_side(limdd_target(P), q))
            o->side++;
        for (int k = 0; k < 2; k++) {
            const uint32_t id = k ? gid : GATEID_H;
            memcpy(h, fa, len * sizeof(EVBDD_WGT));
            dense_gate(h, n, gates[id], 0, q);
            tally(&o->sgate, bqd_gate(P, id, q, n), h, n);
        }
        limdd_refs_pop(1);
    }

    /* diagonal gates: a phase on a random set of qubits, by the walk or PhaseMul */
    {
        const EVBDD_WGT phases[5] = { EVBDD_MIN_ONE, pw[2], pw[1], pw[7], pw[3] };
        const uint64_t A = imask(1 + rnd_below(len - 1), n);
        const EVBDD_WGT phase = phases[rnd_below(5)];
        for (uint64_t x = 0; x < len; x++) h[x] = ((x & A) == A) ? wgt_mul(phase, f[x]) : f[x];
        tally(&o->diag, bqd_apply_diagonal(e, A, phase, n), h, n);
        if (root_is_s(e)) o->phase_s++;
    }

    /* PhaseMul of a random exponent, on the edge and on a labelled one */
    {
        const bool pauli = pauli_session();
        const EVBDD_WGT bases[3] = { pw[1], pw[2], qisq2_lookup(1, 1, 1, 1, 0, 1, 0, 1) };
        const EVBDD_WGT beta = pauli ? pw[1 + rnd_below(7)] : bases[rnd_below(3)];
        const uint32_t r = pauli ? 8 : bqd_exp_order(beta);
        const BQD_EXP eps = mtbdd_refs_push(rand_exp(n, r, ev));
        for (uint64_t x = 0; x < len; x++) h[x] = f[x] == EVBDD_ZERO ? f[x]
                                              : wgt_mul(f[x], wpow(beta, ev[lim_point(x, n)]));
        tally(&o->phase, bqd_phase_mul(e, eps, beta), h, n);
        const lab_t la = rand_lab(n);
        const BQD A1 = limdd_refs_push(lab_edge(la, e, n));
        const BQD C1 = limdd_refs_push(bqd_xpsm_canon(A1));
        lab_apply(la, f, len, fa);
        for (uint64_t x = 0; x < len; x++) h[x] = fa[x] == EVBDD_ZERO ? fa[x]
                                              : wgt_mul(fa[x], wpow(beta, ev[lim_point(x, n)]));
        tally(&o->lphase, bqd_phase_mul(C1, eps, beta), h, n);
        limdd_refs_pop(2);
        if (rnd_below(4) == 0) {                         /* beta^eps itself */
            for (uint64_t x = 0; x < len; x++) h[x] = wpow(beta, ev[lim_point(x, n)]);
            tally(&o->expst, bqd_exp_state(eps, beta, r), h, n);
        }
        mtbdd_refs_pop(1);
    }

    /* MulOff in the translation family: k on supp [U] and c k off it */
    if (!pauli_session()) {
        const EVBDD_WGT cs[4] = { EVBDD_MIN_ONE, pw[1], pw[2], qisq2_lookup(1, 1, 1, 1, 0, 1, 0, 1) };
        const EVBDD_WGT c = cs[rnd_below(4)];
        const LIMDD_TARG U = limdd_edge_is_zero(eg) ? 0 : limdd_target(eg);
        edge_dense(U ? unit(U) : eg, n, gb);
        for (uint64_t x = 0; x < len; x++) h[x] = (gb[x] != EVBDD_ZERO) ? f[x] : wgt_mul(c, f[x]);
        tally(&o->muloff, bqd_mul_off(e, c, U), h, n);
    }

    /* Pair with the partner, at a random qubit, and on labelled operands */
    {
        const unsigned b = (unsigned)rnd_below(n);
        const int s = (int)rnd_below(2), t = (int)rnd_below(2);
        const uint64_t bb = ibit(b, n);
        for (uint64_t x = 0; x < len; x++)
            h[x] = (x & bb) ? g[t ? (x | bb) : (x & ~bb)] : f[s ? (x | bb) : (x & ~bb)];
        tally(&o->pair, bqd_pair(e, s, eg, t, b), h, n);
        const lab_t la = rand_lab(n), lb = rand_lab(n);
        const BQD A1 = limdd_refs_push(lab_edge(la, e, n));
        const BQD B1 = limdd_refs_push(lab_edge(lb, eg, n));
        lab_apply(la, f, len, fa);
        lab_apply(lb, g, len, gb);
        for (uint64_t x = 0; x < len; x++)
            h[x] = (x & bb) ? gb[t ? (x | bb) : (x & ~bb)] : fa[s ? (x | bb) : (x & ~bb)];
        tally(&o->lpair, bqd_xpsm_pair(A1, s, B1, t, b), h, n);
        BQD r1, r2;
        if (!limdd_edge_is_zero(A1) && !limdd_edge_is_zero(B1) && bqd_xp_top(A1) == bqd_xp_top(B1)
            && bqd_xp_top(A1) < b && root_viewed(A1, bb, &r1) && root_viewed(B1, bb, &r2))
            o->pair_nest++;
        limdd_refs_pop(2);
    }

    /* a local matrix on one or two qubits */
    {
        const uint32_t k = (n >= 2 && rnd_below(2)) ? 2 : 1;
        uint32_t qs[2];
        qs[0] = (uint32_t)rnd_below(n);
        if (k == 2) { qs[1] = (uint32_t)rnd_below(n - 1); if (qs[1] >= qs[0]) qs[1]++; }
        const EVBDD_WGT vals[6] = { EVBDD_ZERO, EVBDD_ONE, EVBDD_MIN_ONE, pw[2], pw[1], isq2 };
        EVBDD_WGT M[16];
        const uint32_t dim = 1u << k;
        for (uint32_t i = 0; i < dim * dim; i++) M[i] = vals[rnd_below(6)];
        for (uint64_t x = 0; x < len; x++) {
            uint32_t row = 0;
            for (uint32_t i = 0; i < k; i++) row |= (uint32_t)((x & ibit(qs[i], n)) != 0) << (k - 1 - i);
            EVBDD_WGT v = EVBDD_ZERO;
            for (uint32_t col = 0; col < dim; col++) {
                uint64_t y = x;
                for (uint32_t i = 0; i < k; i++) {
                    y &= ~ibit(qs[i], n);
                    if ((col >> (k - 1 - i)) & 1) y |= ibit(qs[i], n);
                }
                if (M[row * dim + col] != EVBDD_ZERO && f[y] != EVBDD_ZERO)
                    v = wgt_add(v, wgt_mul(M[row * dim + col], f[y]));
            }
            h[x] = v;
        }
        tally(&o->matvec, bqd_local_matvec(e, M, qs, k, n), h, n);
    }

    limdd_refs_pop(2);
    limdd_refs_depths(d1);
    o->refs.tot++;
    if (memcmp(d0, d1, sizeof(d0)) != 0) o->refs.bad++;
    free(h); free(fa); free(gb); free(ev);
}

static void
sel_report(const char *section, const sel_t *o, bool exhaustive)
{
    char what[112];
    const struct { const char *name; const tally_t *t; } k[] = {
        { "bqd_restrict", &o->restr }, { "bqd_project", &o->proj }, { "bqd_x", &o->x },
        { "bqd_perm, the three kinds", &o->perm }, { "bqd_swap", &o->swap },
        { "bqd_pair", &o->pair }, { "bqd_gate, H and a random one", &o->gate },
        { "bqd_gate on one side of the qubit, H and a random one", &o->sgate },
        { "bqd_cgate, controls above", &o->cgate },
        { "bqd_cgate_either, the control below", &o->either },
        { "bqd_apply_diagonal", &o->diag }, { "bqd_phase_mul", &o->phase },
        { "bqd_exp_state", &o->expst }, { "bqd_mul_off", &o->muloff },
        { "bqd_local_matvec", &o->matvec },
        { "Restrict of a labelled edge", &o->lrestr }, { "Project of a labelled edge", &o->lproj },
        { "Pair of labelled edges", &o->lpair }, { "Gate of a labelled edge", &o->lgate },
        { "PhaseMul of a relabelled state", &o->lphase },
        { "the reference stacks as they were", &o->refs },
    };
    for (size_t i = 0; i < sizeof(k) / sizeof(k[0]); i++) {
        if (k[i].t == &o->muloff && pauli_session()) continue;
        snprintf(what, sizeof(what), "%s %s: %s", fam_name(), section, k[i].name);
        report(what, k[i].t->bad, k[i].t->tot);
    }
    const struct { const char *name; unsigned seen; } c[] = {
        { "selections through a Q node's view above the qubit", o->nest_sel },
        { "selections through an S node above the qubit", o->s_sel },
        { "gates whose ratio passes through", o->pass },
        { "gates on a support on one side of the qubit (Side)", o->side },
        { "gates through an S node above the qubit", o->s_gate },
        { "permutations of a Q root above both qubits", o->perm_q },
        { "permutations of an S root above both qubits", o->perm_s },
        { "Pairs of two nested roots above b", o->pair_nest },
        { "diagonal gates on a state with an S root", o->phase_s },
        { "X that changes the tag of the root", o->x_tag },
        { "permutations that change the tag of the root", o->perm_tag },
    };
    const size_t required = exhaustive ? sizeof(c) / sizeof(c[0]) : sizeof(c) / sizeof(c[0]) - 2;
    for (size_t i = 0; i < sizeof(c) / sizeof(c[0]); i++) {
        snprintf(what, sizeof(what), "%s %s: %s", fam_name(), section, c[i].name);
        if (i < required) covered(what, c[i].seen);
        else printf("  %-64s %u\n", what, c[i].seen);
    }
}

/**
 * Every vector over {0, 1, -1, i} on 1 and 2 qubits, and on 3 every THIN'th
 * of them with a stride of SEL_STRIDE, each with a random partner.
 */
#define SEL_STRIDE 3
static void
check_exhaustive_sel(void)
{
    const EVBDD_WGT vals[4] = { EVBDD_ZERO, EVBDD_ONE, EVBDD_MIN_ONE, pw[2] };
    EVBDD_WGT f[8], g[8];
    sel_t o;
    memset(&o, 0, sizeof(o));
    structure_begin();
    for (unsigned n = 1; n <= 3; n++) {
        const uint64_t len = UINT64_C(1) << n, count = UINT64_C(1) << (2 * len);
        const uint64_t step = n == 3 ? SEL_STRIDE * THIN : 1;
        for (uint64_t i = rnd_below(step); i < count; i += step) {
            for (uint64_t x = 0; x < len; x++) f[x] = vals[(i >> (2 * x)) & 3];
            const uint64_t j = rnd_below(count);
            for (uint64_t x = 0; x < len; x++) g[x] = vals[(j >> (2 * x)) & 3];
            const BQD e = limdd_refs_push(bqd_from_vector(f, n));
            const BQD eg = limdd_refs_push(bqd_from_vector(g, n));
            sel_on(f, e, g, eg, n, true, &o);
            limdd_refs_pop(2);
        }
    }
    sel_report("n <= 3, selections, gates, phases", &o, true);
    structure_report("n <= 3, selections, gates, phases");
}

static void
check_random_sel(void)
{
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *g = malloc(MAXV * sizeof(EVBDD_WGT));
    sel_t o;
    memset(&o, 0, sizeof(o));
    structure_begin();
    for (unsigned n = 4; n <= NQ; n++) for (int rep = 0; rep < 40 / THIN; rep++) {
        rand_state(n, (unsigned)rnd_below(S_KINDS), f);
        rand_state(n, (unsigned)rnd_below(S_KINDS), g);
        const BQD e = limdd_refs_push(bqd_from_vector(f, n));
        const BQD eg = limdd_refs_push(bqd_from_vector(g, n));
        sel_on(f, e, g, eg, n, n <= 5, &o);
        limdd_refs_pop(2);
    }
    sel_report("n = 4..8, random, selections, gates, phases", &o, false);
    structure_report("n = 4..8, selections, gates, phases");
    free(f); free(g);
}

/**
 * The work of Side, H on a qubit q where the support lies on one side of it,
 * on random states on 4 to 8 qubits with x_q = 0 or x_q = 1, each from a
 * cleared memo. On x_q = 0 the gate is the restriction and a scale, with no
 * high cofactor (Cof1) and no Apply miss. Without Side the gate recursed on
 * the cofactors wherever a Q node's ratio has a node of q, which a one-sided
 * support puts in every ratio above q, and made Cof1, a product, at each:
 * gate 96 of ising_n16_s3, h q[15] on 28,612 nodes, took 25.8 s in the
 * translation family and 0.010 s in the scalar one. On x_q = 1 it is the
 * restriction and the phase -1 on it; there the root's translation has X_q,
 * and the restriction of the translation family may move a least point and
 * take a product in NEST, so that side is checked for its result only and its
 * work printed.
 */
static void
check_side_work(void)
{
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT)), *h = malloc(MAXV * sizeof(EVBDD_WGT));
    tally_t t = { 0, 0 }, work = { 0, 0 };
    unsigned side = 0;
    uint64_t most[2][2] = { { 0, 0 }, { 0, 0 } };   /* [x_q][Cof1, Apply] */
    structure_begin();
    for (unsigned n = 4; n <= NQ; n++) for (int rep = 0; rep < 24 / THIN; rep++) {
        const uint64_t len = UINT64_C(1) << n;
        const unsigned q = (unsigned)rnd_below(n);
        const int b = (int)rnd_below(2);
        rand_state(n, (unsigned)rnd_below(S_KINDS), f);
        for (uint64_t x = 0; x < len; x++) if (((x & ibit(q, n)) != 0) != b) f[x] = EVBDD_ZERO;
        const BQD e = limdd_refs_push(bqd_from_vector(f, n));
        if (!limdd_edge_is_zero(e) && bqd_xp_top(e) < q) side++;
        memcpy(h, f, len * sizeof(EVBDD_WGT));
        dense_gate(h, n, gates[GATEID_H], 0, q);
        sylvan_clear_cache();
        bqd_counts_reset();
        const BQD r = limdd_refs_push(bqd_gate(e, GATEID_H, q, n));
        uint64_t c[BQD_COUNTS];
        bqd_counts_read(c);
        tally(&t, r, h, n);
        if (b == 0) {
            work.tot++;
            if (c[BQD_COUNT_COF1] != 0 || c[BQD_COUNT_APPLY] != 0) work.bad++;
        }
        if (c[BQD_COUNT_COF1] > most[b][0]) most[b][0] = c[BQD_COUNT_COF1];
        if (c[BQD_COUNT_APPLY] > most[b][1]) most[b][1] = c[BQD_COUNT_APPLY];
        limdd_refs_pop(2);
    }
    char what[112];
    snprintf(what, sizeof(what), "%s Side: H on a support on one side of the qubit", fam_name());
    report(what, t.bad, t.tot);
    snprintf(what, sizeof(what), "%s Side: on x_q = 0 no Cof1 and no Apply miss", fam_name());
    report(what, work.bad, work.tot);
    printf("    most Cof1 / Apply misses: x_q = 0 %llu / %llu, x_q = 1 %llu / %llu\n",
           (unsigned long long)most[0][0], (unsigned long long)most[0][1],
           (unsigned long long)most[1][0], (unsigned long long)most[1][1]);
    snprintf(what, sizeof(what), "%s Side: states with the root above the qubit", fam_name());
    covered(what, side);
    structure_report("Side");
    free(f); free(h);
}

/**
 * bqd_xpsm_diagonal of the zero edge, in the Pauli family: zero, with no node
 * visited. Its entry assertion asked for full support, which the zero edge
 * does not have, so a debug build aborted on the branch that handles it.
 */
static void
check_diagonal_zero(void)
{
    if (!pauli_session()) return;
    uint32_t visits = 1;
    const BQD r = bqd_xpsm_diagonal(limdd_zero_edge(), ibit(0, NQ) | ibit(2, NQ), pw[1], NQ, &visits);
    report("P bqd_xpsm_diagonal of the zero edge: zero, no node visited",
           !limdd_edge_is_zero(r) || visits != 0, 1);
}

/**
 * Random Clifford+T circuits from |0...0>, gate by gate through the runner's
 * entry points (bqd_gate, bqd_cgate with one and two controls, bqd_cgate_either
 * with the control on either side, bqd_swap, CS and CCZ as the runner sends
 * them), the state against the dense one and the edge bqd_from_vector builds
 * after every gate.
 */
static void
check_circuits(void)
{
    EVBDD_WGT *f = malloc(MAXV * sizeof(EVBDD_WGT));
    tally_t t = { 0, 0 };
    unsigned full_gates = 0, s_states = 0;
    structure_begin();
    static const uint32_t one[] = { GATEID_H, GATEID_H, GATEID_S, GATEID_Sdag, GATEID_T, GATEID_Tdag,
                                    GATEID_X, GATEID_Y, GATEID_Z };
    for (unsigned n = 2; n <= NQ; n++) for (int rep = 0; rep < 12 / THIN + 1; rep++) {
        const uint64_t len = UINT64_C(1) << n;
        for (uint64_t x = 0; x < len; x++) f[x] = x ? EVBDD_ZERO : EVBDD_ONE;
        BQD e = limdd_refs_push(bqd_basis_state(0, n));
        const unsigned gates_n = 10 + (unsigned)rnd_below(8 * n);
        for (unsigned i = 0; i < gates_n; i++) {
            if (bqd_has_full_support(e)) full_gates++;
            if (root_is_s(e)) s_states++;
            const unsigned a = (unsigned)rnd_below(n);
            unsigned b = (unsigned)rnd_below(n - 1); if (b >= a) b++;
            BQD r;
            switch (rnd_below(6)) {
            case 0: case 1: case 2: {
                const uint32_t id = one[rnd_below(sizeof(one) / sizeof(one[0]))];
                dense_gate(f, n, gates[id], 0, a);
                r = bqd_gate(e, id, a, n);
                break;
            }
            case 3: {                                    /* cx, cy, cz with the control a */
                static const uint32_t cg[] = { GATEID_X, GATEID_Y, GATEID_Z };
                const uint32_t id = cg[rnd_below(3)];
                dense_gate(f, n, gates[id], UINT64_C(1) << a, b);
                bool ok;
                r = bqd_cgate_either(e, id, a, b, n, &ok);
                break;
            }
            case 4: {                                    /* cs, ch, ccz, ccx */
                const unsigned lo = a < b ? a : b, hi = a < b ? b : a;
                const uint64_t k = rnd_below(4);
                if (k == 0 || n < 3) {
                    dense_gate(f, n, gates[GATEID_S], UINT64_C(1) << lo, hi);
                    r = bqd_cgate(e, GATEID_S, UINT64_C(1) << lo, hi, n);
                } else if (k == 1) {
                    dense_gate(f, n, gates[GATEID_H], UINT64_C(1) << lo, hi);
                    r = bqd_cgate(e, GATEID_H, UINT64_C(1) << lo, hi, n);
                } else {
                    unsigned c = (unsigned)rnd_below(n);
                    while (c == lo || c == hi) c = (unsigned)rnd_below(n);
                    unsigned q[3] = { lo, hi, c };
                    for (int x = 0; x < 2; x++) for (int y = x + 1; y < 3; y++)
                        if (q[y] < q[x]) { const unsigned s = q[x]; q[x] = q[y]; q[y] = s; }
                    const uint64_t cm = (UINT64_C(1) << q[0]) | (UINT64_C(1) << q[1]);
                    const uint32_t id = k == 2 ? GATEID_Z : GATEID_X;
                    dense_gate(f, n, gates[id], cm, q[2]);
                    r = bqd_cgate(e, id, cm, q[2], n);
                }
                break;
            }
            default:
                for (uint64_t x = 0; x < len; x++) {
                    const uint64_t y = perm_index(BQD_PERM_SWAP, a < b ? a : b, a < b ? b : a, x, n);
                    if (y > x) { const EVBDD_WGT v = f[x]; f[x] = f[y]; f[y] = v; }
                }
                r = bqd_swap(e, a, b, n);
                break;
            }
            limdd_refs_pop(1);
            e = limdd_refs_push(r);
            tally(&t, e, f, n);
        }
        limdd_refs_pop(1);
    }
    char what[112];
    snprintf(what, sizeof(what), "%s Clifford+T circuits, the state after every gate", fam_name());
    report(what, t.bad, t.tot);
    snprintf(what, sizeof(what), "%s Clifford+T circuits: gates on full support", fam_name());
    covered(what, full_gates);
    snprintf(what, sizeof(what), "%s Clifford+T circuits: gates on a state with an S root", fam_name());
    covered(what, s_states);
    structure_report("circuits");
    free(f);
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

enum { C_MUL, C_ADD, C_RATIO, C_COMPOSE, C_COFACTOR, C_PRODUCT, C_CANON, C_GATE, C_PERM,
       C_SELECT, C_PHASE, C_PAIR, C_KINDS };

/**
 * Results kept across collections, in two pools that share the tables, on NQ
 * and on NQ - 2 qubits. A random operation takes its operands from one pool
 * and puts its result back, and every one to four operations a collection
 * runs. After it every kept edge must still decode to its vector and still be
 * the edge bqd_from_vector builds, and the set of checked nodes starts empty,
 * since the indices it holds may have been handed out again.
 */
static int
collections(bqd_family_t fam)
{
    const int before = failures;
    bqd_init_rule(fam, BQD_ZERO_SM, GC_WIDTH, 1LL << 16, 1LL << 16, 1LL << 17, 1LL << 16);
    exact_constants();
    rng_state = UINT64_C(0x5EED6C5A0000) + (uint64_t)fam;

    pool_t *pools = malloc(2 * sizeof(pool_t));
    EVBDD_WGT *h = malloc(MAXV * sizeof(EVBDD_WGT)), *g = malloc(MAXV * sizeof(EVBDD_WGT));
    tally_t t = { 0, 0 }, kept = { 0, 0 };
    unsigned collected = 0, freed = 0;
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
            r = bqd_xpsm_apply(BQD_SM_RATIO, p->e[a], m);
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
            const unsigned top = limdd_edge_is_zero(e) ? 0 : bqd_xp_top(e);
            const unsigned v = (unsigned)rnd_below((top < n ? top : n - 1) + 1);
            const int c = (int)rnd_below(2);
            dense_cofactor(fa, n, v, c, h);
            r = bqd_cofactor(e, v, c);
            break;
        }
        case C_CANON: {                              /* a random label times a kept edge */
            const lab_t l = rand_lab(n);
            lab_apply(l, fa, len, h);
            const BQD m = limdd_refs_push(lab_edge(l, p->e[a], n));
            r = bqd_xpsm_canon(m);
            limdd_refs_pop(1);
            break;
        }
        case C_GATE: {                               /* a one-qubit gate, or one with a control */
            const unsigned q = (unsigned)rnd_below(n);
            memcpy(h, fa, len * sizeof(EVBDD_WGT));
            if (q > 0 && rnd_below(2)) {
                const unsigned c = (unsigned)rnd_below(q);
                const uint32_t id = rnd_below(2) ? GATEID_H : GATEID_X;
                dense_gate(h, n, gates[id], UINT64_C(1) << c, q);
                r = bqd_cgate(p->e[a], id, UINT64_C(1) << c, q, n);
            } else {
                const uint32_t id = one_gates[rnd_below(N_ONE_GATES)];
                dense_gate(h, n, gates[id], 0, q);
                r = bqd_gate(p->e[a], id, q, n);
            }
            break;
        }
        case C_PERM: {
            const unsigned qa = (unsigned)rnd_below(n - 1);
            const unsigned qb = qa + 1 + (unsigned)rnd_below(n - 1 - qa);
            const uint32_t kind = (uint32_t)rnd_below(3);
            for (uint64_t x = 0; x < len; x++) h[x] = fa[perm_index(kind, qa, qb, x, n)];
            r = bqd_perm(p->e[a], kind, qa, qb);
            break;
        }
        case C_SELECT: {                             /* a restriction or a projection */
            const unsigned q = (unsigned)rnd_below(n);
            const int c = (int)rnd_below(2);
            if (rnd_below(2)) {
                dense_cofactor(fa, n, q, c, h);
                r = bqd_restrict(p->e[a], q, c);
            } else {
                for (uint64_t x = 0; x < len; x++)
                    h[x] = (((x & ibit(q, n)) != 0) == c) ? fa[x] : EVBDD_ZERO;
                r = bqd_project(p->e[a], q, c);
            }
            break;
        }
        case C_PHASE: {                              /* a phase on one or two qubits */
            const uint64_t A = ibit((unsigned)rnd_below(n), n) | ibit((unsigned)rnd_below(n), n);
            const EVBDD_WGT phase = pw[1 + rnd_below(7)];
            for (uint64_t x = 0; x < len; x++) h[x] = ((x & A) == A) ? wgt_mul(phase, fa[x]) : fa[x];
            r = bqd_apply_diagonal(p->e[a], A, phase, n);
            break;
        }
        case C_PAIR: {
            const unsigned q = (unsigned)rnd_below(n);
            const int s = (int)rnd_below(2), u = (int)rnd_below(2);
            const uint64_t qb = ibit(q, n);
            for (uint64_t x = 0; x < len; x++)
                h[x] = (x & qb) ? fb[u ? (x | qb) : (x & ~qb)] : fa[s ? (x | qb) : (x & ~qb)];
            r = bqd_pair(p->e[a], s, p->e[b], u, q);
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
        /* products and cofactors thin a pool out, so now and then a fresh state */
        if (rnd_below(3) == 0) pool_fresh(p, (unsigned)rnd_below(POOL));

        if (--next > 0) continue;
        next = 1 + (unsigned)rnd_below(4);
        const size_t was = limdd_node_table_count();
        RUN(limdd_gc);
        collected++;
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
    char what[112];
    snprintf(what, sizeof(what), "%s operations on kept results, collected every 1 to 4", fam_name());
    report(what, t.bad, t.tot);
    snprintf(what, sizeof(what), "%s kept results after a collection, decoded and rebuilt", fam_name());
    report(what, kept.bad, kept.tot);
    structure_report("collections");
    snprintf(what, sizeof(what), "%s collections", fam_name());
    covered(what, collected);
    snprintf(what, sizeof(what), "%s collections that freed nodes", fam_name());
    covered(what, freed);
    free(pools); free(h); free(g);
    bqd_quit();
    return failures != before;
}

/* --- harness -------------------------------------------------------------- */

/** A session of fam under SM, which bqd_init_rule gives it, and the checks above. */
static int
run_family(bqd_family_t fam)
{
    const int before = failures;
    bqd_init_rule(fam, BQD_ZERO_SM, NQ, 1LL << 21, 1LL << 20, 1LL << 22, 1LL << 20);
    exact_constants();
    expect(bqd_zero_rule() == BQD_ZERO_SM && bqd_xp_sm() && !bqd_sm(), bqd_family_name(fam),
           "bqd_init_rule with SM does not give this family's SM session");
    rng_state = UINT64_C(0x5EED5A000000) + (uint64_t)fam;
    check_exhaustive_builds();
    check_random_builds();
    check_tag_change();
    check_exhaustive_ops();
    check_random_ops();
    check_exhaustive_sel();
    check_random_sel();
    check_side_work();
    check_diagonal_zero();
    check_circuits();
    printf("  %zu nodes in the table\n", limdd_node_table_count());
    bqd_quit();
    return failures != before;
}

TASK_0(int, run_all)
{
    int bad = 0;
    for (int fam = BQD_FAMILY_X; fam <= BQD_FAMILY_PAULI; fam++) {
        printf("== %s, rule SM, exact weights ==\n", bqd_family_name((bqd_family_t)fam));
        const int res = run_family((bqd_family_t)fam);
        printf("  %s\n", res ? "FAILED" : "ok");
        bad |= res;
    }
    return bad;
}

TASK_0(int, run_collections)
{
    int bad = 0;
    for (int fam = BQD_FAMILY_X; fam <= BQD_FAMILY_PAULI; fam++) {
        printf("== %s, rule SM, collections between operations ==\n",
               bqd_family_name((bqd_family_t)fam));
        const int res = collections((bqd_family_t)fam);
        printf("  %s\n", res ? "FAILED" : "ok");
        bad |= res;
    }
    return bad;
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
    const unsigned workers = getenv("BQD_XP_SM_WORKERS") ? (unsigned)atoi(getenv("BQD_XP_SM_WORKERS"))
                                                          : 4;
    checked = calloc(1, sizeof(nodeset_t));
    if (checked == NULL) { fprintf(stderr, "out of memory\n"); return 1; }
    /* The collection tests keep weight indices across limdd_gc: the dense
     * vectors of the kept results, pw[] and isq2. A collection that finds the
     * weight table more than half full copies the live weights to new
     * indices, which those do not follow, and the checks would then read
     * stale weights, report wrong vectors or crash. The 2^22 entries here are
     * far from half full, which a larger test would reach, so the copy is off. */
    setenv("LIMDD_NO_WGT_GC", "1", 1);
    int bad = 0;

    printf("== BQD, rule SM in the translation and Pauli families, %u workers ==\n", workers);
    session_begin(workers, 20);
    bad |= RUN(run_all);
    session_end();

    /* a small memo, which a collection clears, and which overwrites entries
     * more often in between */
    session_begin(workers, 16);
    bad |= RUN(run_collections);
    session_end();

    free(checked);
    printf("%s\n", bad ? "FAILED" : "ok");
    return bad;
}
