/*
 * The gates without the high cofactor, skip:sec:ratio of the note "Binary
 * Quotient Diagrams with Level Skipping" (bqd-skip.tex), in every family they
 * apply to, each result the edge bqd_from_vector builds for the dense result,
 * so a result right in value and wrong in shape fails too:
 *
 *   exhaustive    every vector over {0, 1, -1, i} on 1 to 3 qubits: bqd_perm
 *                 with the swap and the CX with the control on either side at
 *                 every pair of qubits, bqd_x at every qubit, restriction and
 *                 projection at every qubit (RestrictR and ProjectR in the
 *                 scalar family), bqd_diag_any with every monomial and every
 *                 power of w_8 but 1, Pair at every qubit and every choice of
 *                 the two restrictions with a partner drawn at random, on
 *                 labelled operands in the translation and Pauli families, the
 *                 indicator node Ind, MulOff in the scalar and translation
 *                 families, PhaseMul with a random exponent diagram, and in the
 *                 translation family Canon of (w X^t, N) for every t, which is
 *                 CanonT's wherever t is not the least point
 *   random        the same, drawn, on random states on 4 to 8 qubits: with
 *                 zeros, on random qubits so that levels are skipped, phase,
 *                 circuit and product states, translated and signed in the
 *                 families that can; PhaseMul also with a base that is no root
 *                 of unity, whose exponents are integers
 *   entry points  bqd_gate with X and Y, bqd_cgate with X and Y under one
 *                 control above, bqd_cgate_either with X and Y and the control
 *                 below, which the Hadamard conjugation served for X only,
 *                 bqd_swap, and bqd_apply_diagonal on states with zeros, the
 *                 same edge as the product with the monomial
 *   collections   gate sequences of those gates, the Hadamard and random
 *                 phase multiplications, with limdd_gc between gates in small
 *                 tables, the state checked after every gate; and Sylvan
 *                 collections in the middle of the phase multiplications,
 *                 forced by bqd_exp_collect_every, which the operations
 *                 survive only with their exponents on the mtbdd_refs stack
 *   any worker    phase multiplications on 8 workers in a Sylvan table of
 *                 2^12 buckets, without bqd_exp_collect_every, so that the
 *                 collections in them start from whichever worker fills a
 *                 region, the result checked after every one
 *   stacks        after every operation above, the six reference stacks of
 *                 the calling worker (the LIMDD's and Sylvan's, each with
 *                 values, pointers and tasks) as the operation found them, and
 *                 every worker's at checkpoints, since a leak stays
 *   any thread    bqd_cgate_either with a flip and the control below the
 *                 target, bqd_diag_any, bqd_swap, bqd_gate and bqd_cgate
 *                 called from the main thread, outside any task, as the gates
 *                 of the header may be
 *   regression    the first 593 gates of clifford_T_circuit_20_700, 20
 *                 qubits, whose swaps at gates 573 and 592 took 20 s and more
 *                 than 540 s by the Hadamard conjugation: each swap from a
 *                 cold operation cache has to stay within 20,000 high
 *                 cofactors and 400,000 misses of the pointwise operations,
 *                 where the old route made 1.39 and 35.7 million at gate 573
 *                 with a cache of 2^20, and 7.9 and 88 million at gate 592
 *                 with one of 2^24; the node counts are the canonical ones
 *                 the old route gave, and the swap undone gives the state
 *                 back; and each of the 360 diagonal gates of the prefix,
 *                 all on states without full support, makes no high
 *                 cofactor and no pointwise operation (DiagR)
 *
 * Exact weights, in Q(w_8, sqrt2), throughout, on BQD_RATIO_WORKERS workers
 * (default 4).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sylvan.h>
#include <sylvan_int.h>

#include "qsylvan_bqd.h"
#include "qsylvan_bqd_exp.h"
#include "qsylvan_bqd_ops.h"
#include "qsylvan_bqd_gates.h"
#include "qsylvan_bqd_xp.h"
#include "qsylvan_limdd_gc.h"
#include "qsylvan_limdd_inspect.h"
#include "qsylvan_simulator.h"
#include "sylvan_edge_weights_qisq2.h"

#define NQ 8
#define MAXV (1u << NQ)

/* Under a sanitizer every check runs some twenty times slower: the draws on
 * three qubits are thinned by THIN there, and every kind still runs. */
#if defined(__SANITIZE_ADDRESS__) || defined(__SANITIZE_THREAD__)
#define THIN 8
#elif defined(__has_feature)
#if __has_feature(address_sanitizer) || __has_feature(thread_sanitizer) || \
    __has_feature(memory_sanitizer)
#define THIN 8
#endif
#endif
#ifndef THIN
#define THIN 1
#endif

static int failures = 0;
static EVBDD_WGT pw[8];                     /* w_8^e */
static EVBDD_WGT isq2;                      /* 1/sqrt2 */

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
    expect(bad == 0 && tot > 0, what, buf);
    printf("  %-60s %s (%u)\n", what, bad || !tot ? "FAILED" : "ok", tot);
}

static void
covered(const char *what, unsigned seen)
{
    expect(seen > 0, what, "never met");
    printf("  %-60s %u\n", what, seen);
}

static inline bool parity(uint64_t v) { return __builtin_parityll(v) != 0; }

/** Qubit q as a vector-index bit: qubit 0 is the most significant of n. */
static inline uint64_t
ibit(unsigned q, unsigned n)
{
    return UINT64_C(1) << (n - 1 - q);
}

static inline BQD
unit(LIMDD_TARG t)
{
    return limdd_bundle(LIMDD_LIM_IDENTITY, t);
}

static bool
same_vector(const EVBDD_WGT *a, const EVBDD_WGT *b, unsigned n)
{
    return memcmp(a, b, ((size_t)1 << n) * sizeof(EVBDD_WGT)) == 0;
}

/* --- states: qubit q is bit n-1-q of an index -------------------------------- */

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

/** A random function of the qubits in D (bit q for qubit q) alone, repeated over the others. */
static void
rand_on(unsigned n, uint64_t D, double pzero, bool alg, EVBDD_WGT *f)
{
    uint64_t dm = 0;
    for (unsigned q = 0; q < n; q++) if ((D >> q) & 1) dm |= ibit(q, n);
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
        f[x] = (x & ~dm) ? f[x & dm] : rand_value(pzero, alg);
}

/** U on qubit q where the qubits of cmask (bit c for qubit c) are 1. */
static void
dense_gate(EVBDD_WGT *v, unsigned n, const EVBDD_WGT *U, uint64_t cmask, unsigned q)
{
    const uint64_t qb = ibit(q, n);
    uint64_t cm = 0;
    for (unsigned c = 0; c < n; c++) if ((cmask >> c) & 1) cm |= ibit(c, n);
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++) {
        if ((x & qb) || (x & cm) != cm) continue;
        const EVBDD_WGT a = v[x], b = v[x | qb];
        v[x]      = wgt_add(wgt_mul(U[0], a), wgt_mul(U[1], b));
        v[x | qb] = wgt_add(wgt_mul(U[2], a), wgt_mul(U[3], b));
    }
}

/** A random Clifford+T circuit from |0...0>. */
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

/** c . w_8^{P(x)}, P a random polynomial over Z_8 of degree at most 3 in the qubits of D. */
static void
phase_on(unsigned n, uint64_t D, EVBDD_WGT *f)
{
    uint64_t dm = 0;
    for (unsigned q = 0; q < n; q++) if ((D >> q) & 1) dm |= ibit(q, n);
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

/** A tensor product of random one-qubit states, basis states among them. */
static void
product_state(unsigned n, EVBDD_WGT *f)
{
    EVBDD_WGT a[NQ][2];
    for (unsigned q = 0; q < n; q++) {
        switch (rnd_below(5)) {
        case 0:  a[q][0] = EVBDD_ONE;  a[q][1] = EVBDD_ZERO;    break;
        case 1:  a[q][0] = EVBDD_ZERO; a[q][1] = EVBDD_ONE;     break;
        case 2:  a[q][0] = EVBDD_ONE;  a[q][1] = EVBDD_MIN_ONE; break;
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

enum { S_ZEROS, S_ALGEBRAIC, S_PHASE, S_CIRCUIT, S_PRODUCT, S_KINDS };

static void
rand_state(unsigned n, unsigned kind, EVBDD_WGT *f)
{
    const uint64_t D = rnd_below(UINT64_C(1) << n);
    switch (kind) {
    case S_ZEROS:     rand_on(n, D, 0.2 * (double)(1 + rnd_below(3)), rnd_below(2), f); break;
    case S_ALGEBRAIC: rand_on(n, D, 0.0, true, f); break;
    case S_PHASE:     phase_on(n, D, f); break;
    case S_CIRCUIT:   circuit_state(n, f); break;
    default:          product_state(n, f); break;
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
    EVBDD_WGT g[MAXV];
    memcpy(g, f, len * sizeof(EVBDD_WGT));
    for (uint64_t x = 0; x < len; x++) {
        EVBDD_WGT v = g[x ^ t];
        if (v != EVBDD_ZERO && parity(s & x)) v = wgt_neg(v);
        f[x] = v;
    }
}

/* --- dense references ------------------------------------------------------ */

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

static void
dense_x(const EVBDD_WGT *f, unsigned n, unsigned q, EVBDD_WGT *out)
{
    for (uint64_t y = 0; y < (UINT64_C(1) << n); y++) out[y] = f[y ^ ibit(q, n)];
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

/** [x_b = 0] A|_{x_b = s} + [x_b = 1] B|_{x_b = t}. */
static void
dense_pair(const EVBDD_WGT *A, int s, const EVBDD_WGT *B, int t, unsigned b, unsigned n,
           EVBDD_WGT *out)
{
    const uint64_t bb = ibit(b, n);
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
        out[x] = (x & bb) ? B[t ? (x | bb) : (x & ~bb)] : A[s ? (x | bb) : (x & ~bb)];
}

static void
dense_diag(const EVBDD_WGT *f, unsigned n, uint64_t A, EVBDD_WGT phase, EVBDD_WGT *out)
{
    for (uint64_t x = 0; x < (UINT64_C(1) << n); x++)
        out[x] = ((x & A) == A) ? wgt_mul(f[x], phase) : f[x];
}

/** The function of node N read on n qubits. */
static void
node_vector(LIMDD_TARG N, unsigned n, EVBDD_WGT *out)
{
    bqd_to_vector(unit(N), n, out);
}

/**
 * The label c Z^s X^t (vector-index masks) of the session's family on the
 * node of e, and in out the function it denotes.
 */
static BQD
labelled(BQD e, EVBDD_WGT c, uint64_t s, uint64_t t, unsigned n, EVBDD_WGT *out)
{
    const uint64_t len = UINT64_C(1) << n;
    if (limdd_edge_is_zero(e)) { for (uint64_t y = 0; y < len; y++) out[y] = EVBDD_ZERO; return e; }
    EVBDD_WGT u[MAXV];
    node_vector(limdd_target(e), n, u);
    for (uint64_t y = 0; y < len; y++) {
        EVBDD_WGT v = u[y ^ t];
        if (v != EVBDD_ZERO) {
            if (parity(s & y)) v = wgt_neg(v);
            v = wgt_mul(c, v);
        }
        out[y] = v;
    }
    return limdd_bundle(bqd_lim_make(c, s, t, n), limdd_target(e));
}

/** A random label of the family on e's node, or e itself in the scalar family. */
static BQD
rand_labelled(BQD e, unsigned n, const EVBDD_WGT *f, EVBDD_WGT *out)
{
    const bqd_family_t fam = bqd_family();
    const uint64_t len = UINT64_C(1) << n;
    if (fam == BQD_FAMILY_SCALAR || rnd_below(2)) { memcpy(out, f, len * sizeof(EVBDD_WGT)); return e; }
    const uint64_t s = fam == BQD_FAMILY_PAULI ? rnd_below(len) : 0;
    return labelled(e, pw[rnd_below(8)], s, rnd_below(len), n, out);
}

/*
 * The exponent diagram of the table v, v[x] at the vector index x, reduced
 * modulo r: the Shannon expansion from the top qubit down. Protects what it
 * holds, since a node can make Sylvan collect.
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

/* --- the reference stacks ------------------------------------------------------ */

/*
 * An operation leaves every reference stack as it found it. A leaked push only
 * deepens a stack, and a missing one shows only when a collection lands where
 * it was due, so the depths are compared. The calling worker's six, the
 * LIMDD's and Sylvan's each with values, pointers and tasks, after every
 * operation with what they were before it (OWN). The other workers' run the
 * stolen halves of the operations, and reading them takes every worker
 * (TOGETHER), a millisecond or more on a loaded machine, so they are compared
 * with a baseline at checkpoints: after every random state, gate sequence and
 * round, and every 64 gates of the replay. A leak stays, so a checkpoint
 * finds every one since the baseline.
 */
enum { STACKS = 6, MAXW = 64 };

typedef struct {
    size_t d[STACKS];
    char pad[SYLVAN_SHARING_PAD - (STACKS * sizeof(size_t)) % SYLVAN_SHARING_PAD];
} stack_line_t;

static _Alignas(SYLVAN_SHARING_PAD) stack_line_t stack_base[MAXW];
static _Alignas(SYLVAN_SHARING_PAD) stack_line_t stack_now[MAXW];
static size_t own_before[STACKS];
static unsigned stack_checks, stack_bad;

static void
stack_depths(size_t d[STACKS])
{
    limdd_refs_depths(d);
    mtbdd_refs_depths(d + 3);
}

VOID_TASK_1(stack_read, stack_line_t *, into)
{
    stack_depths(into[LACE_WORKER_ID].d);
}

static bool
stacks_equal(const char *what, unsigned w, const size_t *was, const size_t *is)
{
    static const char *name[STACKS] = { "LIMDD values", "LIMDD pointers", "LIMDD tasks",
                                        "MTBDD values", "MTBDD pointers", "MTBDD tasks" };
    bool ok = true;
    for (int k = 0; k < STACKS; k++) {
        if (was[k] == is[k]) continue;
        if (ok && stack_bad < 20)
            fprintf(stderr, "FAIL reference stacks after %s: worker %u, %s at %zu, was %zu\n",
                    what, w, name[k], is[k], was[k]);
        ok = false;
    }
    return ok;
}

/** Every worker's stacks now, as the baseline of the checkpoints to come. */
static void
stacks_baseline(void)
{
    stack_checks = stack_bad = 0;
    TOGETHER(stack_read, stack_base);
}

/** A checkpoint: every worker's stacks as at the baseline, which a leak moves. */
static void
stacks_checkpoint(const char *what)
{
    TOGETHER(stack_read, stack_now);
    bool ok = true;
    for (unsigned w = 0; w < lace_workers(); w++)
        if (!stacks_equal(what, w, stack_base[w].d, stack_now[w].d)) ok = false;
    stack_checks++;
    if (!ok) {
        stack_bad++;
        memcpy(stack_base, stack_now, sizeof(stack_base));
    }
}

static void
own_mark(void)
{
    stack_depths(own_before);
}

static uint64_t
own_check(uint64_t result)
{
    size_t d[STACKS];
    stack_depths(d);
    stack_checks++;
    if (!stacks_equal("an operation", 0, own_before, d)) stack_bad++;
    return result;
}

/** e, an operation's result, with the calling worker's stacks as it found them */
#define OWN(e) (own_mark(), own_check(e))

static void
report_stacks(void)
{
    report("reference stacks: every operation and checkpoint balanced", stack_bad, stack_checks);
}

/* --- the checks -------------------------------------------------------------- */

typedef struct { unsigned bad, tot; } tally_t;

static void
tally(tally_t *t, BQD got, const EVBDD_WGT *h, unsigned n)
{
    t->tot++;
    if (got != bqd_from_vector(h, n)) t->bad++;
}

enum {
    T_PERM, T_X, T_RESTRICT, T_PROJECT, T_DIAG, T_PAIR, T_IND, T_MULOFF, T_PHASEMUL,
    T_PHASEINT, T_CANONT, T_GATE_X, T_GATE_Y, T_CX, T_CY, T_EITHER, T_SWAP, T_APPLY_DIAG,
    T_EXP, T_COUNT
};
static const char *tally_name[T_COUNT] = {
    "bqd_perm, SW, CX down and CX up", "bqd_x", "bqd_restrict (RestrictR)",
    "bqd_project (ProjectR)", "bqd_diag_any (DiagR)", "bqd_pair (Pair)", "bqd_ind (Ind)",
    "bqd_mul_off (MulOff)", "bqd_phase_mul, root-of-unity base", "bqd_phase_mul, integer exponents",
    "bqd_xp_canon of (w X^t, N) (CanonT)", "bqd_gate X", "bqd_gate Y",
    "bqd_cgate X, one control above", "bqd_cgate Y, one control above",
    "bqd_cgate_either X and Y, control below", "bqd_swap",
    "bqd_apply_diagonal with zeros, as the product", "bqd_exp_state (Exp)",
};

enum { C_SKIP_A, C_SKIP_B, C_ZEROS, C_MULOFF_WORK, C_CANONT, C_ROOT_SKIP, C_COUNT };
static const char *cover_name[C_COUNT] = {
    "perms whose state skips the upper qubit", "perms whose state skips the lower qubit",
    "diagonal gates on states with zeros", "MulOff that multiplies off a support",
    "translated edges that are not canonical (CanonT)", "states whose root edge skips",
};

/**
 * Every operation on F, the diagram of f on n qubits: with `all` every
 * choice of its parameters, otherwise one drawn. g is a partner function,
 * with G its diagram.
 */
static void
check_state(BQD F, const EVBDD_WGT *f, BQD G, const EVBDD_WGT *g, unsigned n, bool all,
            tally_t *t, unsigned *cov)
{
    const bqd_family_t fam = bqd_family();
    const uint64_t len = UINT64_C(1) << n;
    EVBDD_WGT h[MAXV], a[MAXV], b[MAXV], u[MAXV];
    const unsigned rootvar = limdd_edge_is_zero(F) ? n : limdd_level(limdd_target(F));
    if (!limdd_edge_is_zero(F) && rootvar > 0 && rootvar < n) cov[C_ROOT_SKIP]++;
    bool zeros = false;
    for (uint64_t x = 0; x < len; x++) if (f[x] == EVBDD_ZERO) zeros = true;

    /* permutations of two qubits */
    if (n >= 2) {
        const unsigned reps = all ? n * (n - 1) / 2 : 2;
        unsigned k = 0;
        for (unsigned qa = 0; qa < n; qa++) for (unsigned qb = qa + 1; qb < n; qb++) {
            unsigned pa = qa, pb = qb;
            if (!all) {
                if (k++ >= reps) break;
                pa = (unsigned)rnd_below(n - 1);
                pb = pa + 1 + (unsigned)rnd_below(n - 1 - pa);
            }
            for (uint32_t kind = 0; kind <= BQD_PERM_CX_UP; kind++) {
                dense_perm(f, n, kind, pa, pb, h);
                tally(&t[T_PERM], OWN(bqd_perm(F, kind, pa, pb)), h, n);
            }
            if (!limdd_edge_is_zero(F)) {
                bool dep_a = false, dep_b = false;
                for (uint64_t x = 0; x < len; x++) {
                    if (f[x] != f[x ^ ibit(pa, n)]) dep_a = true;
                    if (f[x] != f[x ^ ibit(pb, n)]) dep_b = true;
                }
                if (!dep_a) cov[C_SKIP_A]++;
                if (!dep_b) cov[C_SKIP_B]++;
            }
        }
    }

    /* X, restriction and projection */
    for (unsigned q = 0; q < n; q++) {
        if (!all && rnd_below(n) != 0) continue;
        dense_x(f, n, q, h);
        tally(&t[T_X], OWN(bqd_x(F, q)), h, n);
        for (int bb = 0; bb < 2; bb++) {
            dense_restrict(f, n, q, bb, h);
            tally(&t[T_RESTRICT], OWN(bqd_restrict(F, q, bb)), h, n);
            dense_project(f, n, q, bb, h);
            tally(&t[T_PROJECT], OWN(bqd_project(F, q, bb)), h, n);
        }
    }

    /* DiagR: every monomial with every phase, or one drawn */
    for (uint64_t A = 1; A < len; A++) for (int k = 1; k < 8; k++) {
        if (!all && rnd_below(len * 7 / 3 + 1) != 0) continue;
        dense_diag(f, n, A, pw[k], h);
        tally(&t[T_DIAG], OWN(bqd_diag_any(F, A, pw[k], n)), h, n);
        if (zeros) cov[C_ZEROS]++;
    }

    /* Pair, on labelled operands in the other families */
    for (unsigned bq = 0; bq < n; bq++) for (int st = 0; st < 4; st++) {
        if (!all && rnd_below(2 * n) != 0) continue;
        const BQD LA = rand_labelled(F, n, f, a), LB = rand_labelled(G, n, g, b);
        dense_pair(a, st >> 1, b, st & 1, bq, n, h);
        tally(&t[T_PAIR], OWN(bqd_pair(LA, st >> 1, LB, st & 1, bq)), h, n);
    }

    if (limdd_edge_is_zero(F)) return;
    const LIMDD_TARG N = limdd_target(F);

    /* Ind */
    node_vector(N, n, u);
    for (uint64_t x = 0; x < len; x++) h[x] = u[x] == EVBDD_ZERO ? EVBDD_ZERO : EVBDD_ONE;
    tally(&t[T_IND], unit(OWN(bqd_ind(N))), h, n);

    /* MulOff, k = F and the support of G's node */
    if (fam != BQD_FAMILY_PAULI && !limdd_edge_is_zero(G)) {
        const EVBDD_WGT c = rnd_below(2) ? pw[1 + rnd_below(7)] : rand_alg();
        node_vector(limdd_target(G), n, u);
        bool off = false;
        for (uint64_t x = 0; x < len; x++) {
            h[x] = (u[x] != EVBDD_ZERO || f[x] == EVBDD_ZERO) ? f[x] : wgt_mul(c, f[x]);
            if (h[x] != f[x]) off = true;
        }
        if (off) cov[C_MULOFF_WORK]++;
        tally(&t[T_MULOFF], OWN(bqd_mul_off(F, c, limdd_target(G))), h, n);
    }

    /* PhaseMul with a random exponent: a power of w_8 in every family, and
     * a base without finite order in the scalar and translation families */
    int64_t ev[MAXV];
    const int m = 1 + (int)rnd_below(7);
    const EVBDD_WGT beta = pw[m];
    const uint32_t r = bqd_exp_order(beta);
    for (uint64_t x = 0; x < len; x++) ev[x] = (int64_t)rnd_below(8);
    const BQD_EXP eps = mtbdd_refs_push(exp_table(ev, 0, n, r));
    for (uint64_t x = 0; x < len; x++) h[x] = wgt_mul(f[x], bqd_exp_power(beta, ev[x]));
    tally(&t[T_PHASEMUL], OWN(bqd_phase_mul(F, eps, beta)), h, n);
    mtbdd_refs_pop(1);
    if (fam != BQD_FAMILY_PAULI) {
        const EVBDD_WGT gamma = qisq2_lookup(1, 1, 1, 1, 0, 1, 0, 1);      /* 1 + sqrt2 */
        for (uint64_t x = 0; x < len; x++) ev[x] = (int64_t)rnd_below(5) - 2;
        const BQD_EXP ei = mtbdd_refs_push(exp_table(ev, 0, n, 0));
        for (uint64_t x = 0; x < len; x++) h[x] = wgt_mul(f[x], bqd_exp_power(gamma, ev[x]));
        tally(&t[T_PHASEINT], OWN(bqd_phase_mul(F, ei, gamma)), h, n);
        /* and Exp of the same exponent, a function of full support */
        for (uint64_t x = 0; x < len; x++) h[x] = bqd_exp_power(gamma, ev[x]);
        const BQD ex = OWN(bqd_exp_state(ei, gamma, 0));
        t[T_EXP].tot++;
        EVBDD_WGT d[MAXV];
        bqd_to_vector(ex, n, d);
        if (!same_vector(d, h, n) || ex != bqd_from_vector(h, n)) t[T_EXP].bad++;
        mtbdd_refs_pop(1);
    }

    /* CanonT: (w X^t, N) for every t, in the translation family */
    if (fam == BQD_FAMILY_X) {
        for (uint64_t tt = 0; tt < len; tt++) {
            if (!all && rnd_below(len / 4 + 1) != 0) continue;
            const BQD L = labelled(F, pw[1], 0, tt, n, h);
            const BQD c = OWN(bqd_xp_canon(L));
            if (c != L) cov[C_CANONT]++;
            tally(&t[T_CANONT], c, h, n);
        }
    }
}

/* The gates of the runner through their entry points. */
static void
check_entry_points(BQD F, const EVBDD_WGT *f, unsigned n, tally_t *t)
{
    const uint64_t len = UINT64_C(1) << n;
    EVBDD_WGT h[MAXV];
    for (unsigned q = 0; q < n; q++) {
        memcpy(h, f, len * sizeof(EVBDD_WGT));
        dense_gate(h, n, gates[GATEID_X], 0, q);
        tally(&t[T_GATE_X], OWN(bqd_gate(F, GATEID_X, q, n)), h, n);
        memcpy(h, f, len * sizeof(EVBDD_WGT));
        dense_gate(h, n, gates[GATEID_Y], 0, q);
        tally(&t[T_GATE_Y], OWN(bqd_gate(F, GATEID_Y, q, n)), h, n);
    }
    if (n < 2) return;
    for (int k = 0; k < 3; k++) {
        const unsigned c = (unsigned)rnd_below(n - 1);
        const unsigned tq = c + 1 + (unsigned)rnd_below(n - 1 - c);
        memcpy(h, f, len * sizeof(EVBDD_WGT));
        dense_gate(h, n, gates[GATEID_X], UINT64_C(1) << c, tq);
        tally(&t[T_CX], OWN(bqd_cgate(F, GATEID_X, UINT64_C(1) << c, tq, n)), h, n);
        memcpy(h, f, len * sizeof(EVBDD_WGT));
        dense_gate(h, n, gates[GATEID_Y], UINT64_C(1) << c, tq);
        tally(&t[T_CY], OWN(bqd_cgate(F, GATEID_Y, UINT64_C(1) << c, tq, n)), h, n);
        /* the control below the target, which the Hadamard conjugation served for X alone */
        for (int y = 0; y < 2; y++) {
            const uint32_t gid = y ? GATEID_Y : GATEID_X;
            memcpy(h, f, len * sizeof(EVBDD_WGT));
            dense_gate(h, n, gates[gid], UINT64_C(1) << tq, c);
            bool ok;
            const BQD r = OWN(bqd_cgate_either(F, gid, tq, c, n, &ok));
            tally(&t[T_EITHER], r, h, n);
            if (!ok) t[T_EITHER].bad++;
        }
        const unsigned sa = (unsigned)rnd_below(n);
        unsigned sb = (unsigned)rnd_below(n - 1); if (sb >= sa) sb++;
        dense_perm(f, n, BQD_PERM_SWAP, sa < sb ? sa : sb, sa < sb ? sb : sa, h);
        tally(&t[T_SWAP], OWN(bqd_swap(F, sa, sb, n)), h, n);
    }
    /* a diagonal gate on a state with a zero: DiagR, the product's edge */
    uint64_t A = 0;
    while (A == 0) A = rnd_below(len);
    const EVBDD_WGT phase = pw[1 + rnd_below(7)];
    dense_diag(f, n, A, phase, h);
    const BQD d = OWN(bqd_apply_diagonal(F, A, phase, n));
    tally(&t[T_APPLY_DIAG], d, h, n);
    if (!bqd_has_full_support(F) && d != bqd_multiply(F, bqd_monomial(A, phase, n)))
        t[T_APPLY_DIAG].bad++;
}

static void
check_family(void)
{
    const EVBDD_WGT vals[4] = { EVBDD_ZERO, EVBDD_ONE, EVBDD_MIN_ONE, pw[2] };
    EVBDD_WGT f[MAXV], g[MAXV];
    tally_t t[T_COUNT];
    unsigned cov[C_COUNT];
    memset(t, 0, sizeof(t));
    memset(cov, 0, sizeof(cov));

    stacks_baseline();
    for (unsigned n = 1; n <= 3; n++) {
        const uint64_t len = UINT64_C(1) << n, count = UINT64_C(1) << (2 * len);
        for (uint64_t idx = 0; idx < count; idx++) {
            if (n == 3 && THIN > 1 && rnd_below(THIN) != 0) continue;
            for (uint64_t x = 0; x < len; x++) f[x] = vals[(idx >> (2 * x)) & 3];
            const uint64_t j = rnd_below(count);
            for (uint64_t x = 0; x < len; x++) g[x] = vals[(j >> (2 * x)) & 3];
            const BQD F = bqd_from_vector(f, n), G = bqd_from_vector(g, n);
            check_state(F, f, G, g, n, true, t, cov);
            if (n <= 2 || rnd_below(16) == 0) check_entry_points(F, f, n, t);
        }
    }
    stacks_checkpoint("the exhaustive checks");
    printf("  exhaustive on 1 to 3 qubits done\n");
    for (unsigned n = 4; n <= NQ; n++) for (int rep = 0; rep < 24; rep++) {
        rand_state(n, (unsigned)rnd_below(S_KINDS), f);
        if (rnd_below(2)) twist(n, f);
        rand_state(n, (unsigned)rnd_below(S_KINDS), g);
        if (rnd_below(2)) twist(n, g);
        const BQD F = bqd_from_vector(f, n), G = bqd_from_vector(g, n);
        check_state(F, f, G, g, n, false, t, cov);
        check_entry_points(F, f, n, t);
        stacks_checkpoint("a random state");
    }

    const bqd_family_t fam = bqd_family();
    for (int i = 0; i < T_COUNT; i++) {
        if (fam == BQD_FAMILY_PAULI && (i == T_MULOFF || i == T_PHASEINT || i == T_EXP)) continue;
        if (fam != BQD_FAMILY_X && i == T_CANONT) continue;
        report(tally_name[i], t[i].bad, t[i].tot);
    }
    for (int i = 0; i < C_COUNT; i++) {
        if (fam == BQD_FAMILY_PAULI && i == C_MULOFF_WORK) continue;
        if (fam != BQD_FAMILY_X && i == C_CANONT) continue;
        covered(cover_name[i], cov[i]);
    }
    report_stacks();
}

/* --- collections ------------------------------------------------------------- */

/* Sylvan's collections, and those that ran while a phase multiplication was
 * under way: the hook runs inside the collection's frame, which every worker
 * has joined, so it sees the flag the caller set before the call. */
static unsigned sylvan_collections, sylvan_collections_inside;
static volatile int in_phase_mul;

VOID_TASK_0(count_sylvan_collection)
{
    sylvan_collections++;
    if (in_phase_mul) sylvan_collections_inside++;
}

/**
 * A gate sequence of the gates above, the Hadamard and random phase
 * multiplications, from a basis state, the state checked against the dense
 * one after every gate, with limdd_gc about every fifth gate. The state is
 * protected, and nothing else is live between gates. The phase
 * multiplications, of diagonal gates on any support and of random exponents,
 * run with a Sylvan collection at every 64th exponent node this worker makes
 * (bqd_exp_collect_every), which is started by that worker alone: Sylvan's
 * table is large enough that none other starts one.
 */
static void
gate_sequences(tally_t *t, unsigned *collected, int runs)
{
    EVBDD_WGT f[MAXV];
    BQD state = limdd_zero_edge();
    limdd_protect(&state);
    for (int run = 0; run < runs; run++) {
        const unsigned n = 2 + (unsigned)rnd_below(NQ - 1);
        const uint64_t len = UINT64_C(1) << n;
        const uint64_t x0 = rnd_below(len);
        for (uint64_t x = 0; x < len; x++) f[x] = (x == x0) ? EVBDD_ONE : EVBDD_ZERO;
        state = bqd_basis_state(x0, n);
        for (int k = 0; k < 40; k++) {
            const unsigned q = (unsigned)rnd_below(n);
            unsigned c = (unsigned)rnd_below(n - 1); if (c >= q) c++;
            switch (rnd_below(9)) {
            case 0: {
                const uint32_t gid = rnd_below(3) ? GATEID_H : GATEID_Y;
                dense_gate(f, n, gates[gid], 0, q);
                state = OWN(bqd_gate(state, gid, q, n));
                break;
            }
            case 1: {
                const uint32_t gid = rnd_below(2) ? GATEID_X : GATEID_Y;
                bool ok;
                dense_gate(f, n, gates[gid], UINT64_C(1) << c, q);
                state = OWN(bqd_cgate_either(state, gid, c, q, n, &ok));
                if (!ok) t->bad++;
                break;
            }
            case 2: {
                EVBDD_WGT h[MAXV];
                dense_perm(f, n, BQD_PERM_SWAP, c < q ? c : q, c < q ? q : c, h);
                memcpy(f, h, len * sizeof(EVBDD_WGT));
                state = OWN(bqd_swap(state, c, q, n));
                break;
            }
            case 3:
                dense_gate(f, n, gates[GATEID_X], 0, q);
                state = OWN(bqd_gate(state, GATEID_X, q, n));
                break;
            case 4: case 5: {
                /* a random exponent: many exponent nodes, which fill a small
                 * Sylvan table, so that it collects inside the multiplication */
                int64_t ev[MAXV];
                const EVBDD_WGT beta = pw[1 + rnd_below(7)];
                for (uint64_t x = 0; x < len; x++) ev[x] = (int64_t)rnd_below(8);
                const BQD_EXP eps = mtbdd_refs_push(exp_table(ev, 0, n, bqd_exp_order(beta)));
                for (uint64_t x = 0; x < len; x++) f[x] = wgt_mul(f[x], bqd_exp_power(beta, ev[x]));
                in_phase_mul = 1;
                state = OWN(bqd_phase_mul(state, eps, beta));
                in_phase_mul = 0;
                mtbdd_refs_pop(1);
                break;
            }
            default: {
                uint64_t A = 0;
                while (A == 0) A = rnd_below(len);
                const EVBDD_WGT phase = pw[1 + rnd_below(7)];
                dense_diag(f, n, A, phase, f);
                in_phase_mul = 1;
                state = OWN(bqd_apply_diagonal(state, A, phase, n));
                in_phase_mul = 0;
                break;
            }
            }
            t->tot++;
            if (state != bqd_from_vector(f, n)) t->bad++;
            if (rnd_below(5) == 0) { RUN(limdd_gc); (*collected)++; }
        }
        stacks_checkpoint("a gate sequence");
    }
    limdd_unprotect(&state);
}

TASK_1(int, run_collections, int, fam)
{
    bqd_init((bqd_family_t)fam, NQ + 2, 1LL << 16, 1LL << 16, 1LL << 18, 1LL << 16);
    exact_constants();
    rng_state = UINT64_C(0xC011EC7) + (uint64_t)fam;
    const int before = failures;
    tally_t t = { 0, 0 };
    unsigned collected = 0;
    sylvan_collections = sylvan_collections_inside = 0;
    stacks_baseline();
    bqd_exp_collect_every(64);
    gate_sequences(&t, &collected, 30);
    bqd_exp_collect_every(0);
    report("gate sequences, limdd_gc every ~5 gates", t.bad, t.tot);
    covered("collections between gates", collected);
    covered("Sylvan collections", sylvan_collections);
    covered("Sylvan collections inside a phase multiplication", sylvan_collections_inside);
    report_stacks();
    bqd_quit();
    return failures != before;
}

/* --- regression: the swaps that stalled clifford_T_circuit_20_700 ------------ */

/*
 * Its first 593 gates, each a name and one or two qubits, "cz4.13" for cz
 * q[4],q[13]. The circuit comes from the diagnosis of the BQD's time-outs;
 * gates 573 and 592 are the two swaps whose Hadamard conjugation stalled.
 */
static const char *const clifford_T_prefix =
    "t11 sdg4 t1 cz4.13 y18 t10 cz9.17 cx4.16 s3 s19 cx13.18 h17 cx16.2 cz16.1 t7 z18 sdg3 "
    "z19 t5 y8 z0 t4 cx16.17 cz5.4 x18 t17 s1 sdg2 sdg19 z7 t6 cz6.1 t7 s11 s3 swap18.13 "
    "t0 t12 t16 cx5.2 cx11.3 h3 h0 swap0.11 t2 s13 x9 x11 t5 y1 x16 swap9.5 cz19.11 y16 t2 "
    "sdg7 y8 h12 t2 x11 swap0.2 t13 t10 h19 h0 t17 x14 t4 cx14.19 cx3.6 cx2.19 cz8.4 x2 "
    "sdg19 cx9.11 cx4.11 x7 h3 swap7.2 z14 s9 s1 swap19.15 x16 cz19.1 t2 x8 t18 y19 x15 "
    "swap11.0 t7 y6 t7 cz8.3 t0 y6 cz19.11 t1 cx17.8 t19 cx0.6 t17 cx19.2 t16 z8 cz19.5 h3 "
    "y1 t19 cz14.5 x9 s13 cz2.0 cz0.18 x9 t3 swap0.12 cx3.14 z12 cz17.7 t14 swap0.1 t6 "
    "cx14.7 swap14.17 h16 h11 cz18.9 cx1.11 t8 cx3.19 swap2.15 cz6.17 sdg4 y18 t1 t8 t3 t8 "
    "y10 x19 cx1.3 cx4.18 h14 s2 sdg14 t2 swap4.17 t10 z13 t1 z6 t5 y13 x6 z8 t18 cx1.4 "
    "t12 t1 t1 s2 h3 cx14.7 t7 s15 t3 t18 swap15.2 t10 x11 h16 s18 t2 s11 t13 t10 t9 "
    "swap11.3 cx7.10 y2 s10 y4 cz4.15 cx2.13 x3 swap12.9 swap3.2 swap17.4 t10 t10 t0 sdg3 "
    "y0 y0 cz9.16 z6 t18 cx14.3 y7 cz9.18 t12 cx3.9 swap16.4 t13 cz11.8 x9 t17 swap12.16 "
    "y12 sdg14 t19 y5 swap5.1 cz0.4 cx13.3 swap17.8 swap5.8 t9 sdg11 t14 y17 t17 t10 y18 "
    "sdg3 cz0.3 y3 cz7.4 t3 x14 cx0.12 x12 swap8.4 cx16.12 t13 t9 cx15.8 h4 s3 sdg13 "
    "cx15.8 x12 z19 t18 cz14.19 cz3.11 z0 swap16.2 cz9.5 swap8.18 swap0.13 t6 x2 z5 t8 t6 "
    "x13 t0 cz2.1 x14 cz19.14 t2 cx3.17 s2 t8 cz19.17 cx8.6 cx9.4 t1 t4 y5 t3 cz7.18 t5 "
    "cz5.3 cz16.1 t19 cz12.4 t3 t17 t14 x9 swap11.1 t3 swap5.18 t0 cx18.13 sdg12 t17 t16 "
    "swap16.18 cx10.14 t14 t15 cx17.0 t5 cx4.19 t16 t9 cz2.19 t10 z8 z5 s6 z9 y1 cx5.3 t16 "
    "t14 sdg6 h16 swap1.10 cz9.11 cx11.6 sdg15 sdg9 t13 t17 s11 t10 t7 t17 x15 x18 cx6.3 "
    "z6 t18 swap16.1 t1 h14 swap17.1 t5 cz15.0 t4 t19 t0 t2 x6 cx4.11 cx17.6 y10 x13 t12 "
    "t18 t10 s19 t11 cz11.3 cx9.7 t10 t6 cx10.2 cz19.1 h3 t18 t6 z12 t14 h4 t1 t10 cz5.12 "
    "h15 swap13.14 t9 cz16.6 t12 t7 t15 swap16.5 cx10.11 cz6.1 swap19.1 sdg5 y12 cx4.7 t5 "
    "cx5.13 h18 cz10.14 t13 cx2.6 t14 swap7.6 t15 t4 y16 h3 h16 swap3.17 t0 t12 sdg18 t16 "
    "cx4.14 t12 t6 cz8.19 t14 cz12.3 cx13.14 y13 cz16.8 cz7.10 y9 cx13.10 cx3.17 y13 sdg18 "
    "y13 sdg16 s8 t10 swap16.19 t16 t0 swap19.11 z16 cx5.2 t11 cz18.6 cx5.10 t13 y2 s18 "
    "swap0.17 x11 swap14.10 cz9.1 cz6.4 t17 swap3.14 h14 s16 y6 t5 t7 t3 t12 t15 t1 cz1.0 "
    "cz14.16 x19 t17 t0 t17 cz11.14 x5 swap18.11 cx9.12 cz1.3 s1 t9 cx10.18 sdg16 t16 t8 "
    "t14 x3 swap5.19 y3 y10 z12 cz14.0 s15 x4 cx5.18 s15 t0 z10 sdg11 cx3.1 s12 sdg18 t9 "
    "t9 sdg18 swap5.1 z0 swap1.5 cz13.11 t7 sdg6 t3 t9 z11 t17 sdg3 s9 t11 swap19.9 t16 "
    "t13 t10 cx16.6 t14 t12 t4 sdg15 cx3.13 t18 t11 cx6.10 cx11.6 sdg14 swap10.17 cx8.16 "
    "s10 s18 s15 sdg19 cz8.9 cz3.7 t18 t4 t12 cx7.17 x5 t11 cz5.1 y1 t7 cx10.16 t19 "
    "cz10.16 t10 s18 y18 s8 h9 t13 s8 cx13.14 y1 z10 cx7.6 t4 cx17.8 t3 t11 swap2.1 h13 "
    "cz19.1 t6 cx9.19 sdg19 t5 t17 t3 h16 t0 h12 t4 sdg15 cz13.15 t4 t18 t14 sdg2 cz16.18 "
    "y10 t15 t19 swap6.3 x5 cx11.3 cx14.9 cz13.2 t12 cx0.15 swap0.7 t3 x15 swap7.14 t0 "
    "sdg10 t5 cz0.12 h7 t18 z7 z6 t1 t4 y4 t12 x10 y17 t19 swap16.2 ";

static BQD
replay_gate(BQD e, const char *name, int a, int b, uint32_t n)
{
    bool ok = true;
    BQD r;
    if (!strcmp(name, "cx"))        r = bqd_cgate_either(e, GATEID_X, (uint32_t)a, (uint32_t)b, n, &ok);
    else if (!strcmp(name, "cz"))   r = bqd_cgate_either(e, GATEID_Z, (uint32_t)a, (uint32_t)b, n, &ok);
    else if (!strcmp(name, "swap")) r = bqd_swap(e, (uint32_t)a, (uint32_t)b, n);
    else {
        static const struct { const char *name; uint32_t id; } one[] = {
            { "x", GATEID_X }, { "y", GATEID_Y }, { "z", GATEID_Z }, { "h", GATEID_H },
            { "s", GATEID_S }, { "sdg", GATEID_Sdag }, { "t", GATEID_T }, { "tdg", GATEID_Tdag },
        };
        uint32_t id = GATEID_I;
        for (size_t i = 0; i < sizeof(one) / sizeof(one[0]); i++)
            if (!strcmp(name, one[i].name)) id = one[i].id;
        if (id == GATEID_I) { fprintf(stderr, "unknown gate %s\n", name); exit(1); }
        r = bqd_gate(e, id, (uint32_t)a, n);
    }
    if (!ok) { fprintf(stderr, "gate %s %d,%d refused\n", name, a, b); exit(1); }
    return r;
}

TASK_1(int, run_regression, int, fam)
{
    enum { RN = 20 };
    /* the canonical node counts after gate 573, which the old route gave too,
     * and after gate 592, per family */
    static const size_t after573[3] = { 1089, 151, 71 };
    static const size_t after592[3] = { 218, 146, 63 };
    bqd_init((bqd_family_t)fam, RN, 1LL << 20, 1LL << 20, 1LL << 22, 1LL << 20);
    const int before = failures;
    BQD state = bqd_basis_state(0, RN), prev = state;
    limdd_protect(&state);
    limdd_protect(&prev);
    unsigned bad = 0, gate = 0;
    unsigned diag_gates = 0, diag_partial = 0, diag_costly = 0;
    uint64_t diag_phasemul = 0;
    stacks_baseline();
    const char *p = clifford_T_prefix;
    char what[128];
    while (*p) {
        while (*p == ' ') p++;
        if (!*p) break;
        char name[8];
        int len = 0, a = 0, b = -1;
        while (*p >= 'a' && *p <= 'z' && len < 7) name[len++] = *p++;
        name[len] = 0;
        a = (int)strtol(p, (char **)&p, 10);
        if (*p == '.') b = (int)strtol(p + 1, (char **)&p, 10);
        const bool watched = gate == 573 || gate == 592;
        if (watched) {
            sylvan_clear_cache();
            bqd_counts_reset();
            prev = state;
        }
        /* a diagonal gate, on any support, is DiagR: no high cofactor and no
         * pointwise operation at all, where the product with the monomial
         * made 186,000 and 578,000 on one T gate of hidden-shift_n30 */
        const bool diag = !strcmp(name, "t") || !strcmp(name, "tdg") || !strcmp(name, "s") ||
                          !strcmp(name, "sdg") || !strcmp(name, "z") || !strcmp(name, "cz");
        if (diag) {
            diag_gates++;
            if (!bqd_has_full_support(state)) diag_partial++;
            bqd_counts_reset();
        }
        state = OWN(replay_gate(state, name, a, b, RN));
        if (diag) {
            uint64_t c[BQD_COUNTS];
            bqd_counts_read(c);
            if (c[BQD_COUNT_COF1] != 0 || c[BQD_COUNT_APPLY] != 0) diag_costly++;
            diag_phasemul += c[BQD_COUNT_PHASEMUL];
        }
        if (watched) {
            uint64_t c[BQD_COUNTS];
            bqd_counts_read(c);
            const size_t nodes = limdd_nodecount(state, RN);
            const size_t want = gate == 573 ? after573[fam] : after592[fam];
            printf("  gate %u, %s q%d,q%d: %zu nodes, %llu high cofactors, %llu Apply misses\n",
                   gate, name, a, b, nodes, (unsigned long long)c[BQD_COUNT_COF1],
                   (unsigned long long)c[BQD_COUNT_APPLY]);
            snprintf(what, sizeof(what), "gate %u: within 20,000 Cof1 and 400,000 Apply misses", gate);
            report(what, c[BQD_COUNT_COF1] > 20000 || c[BQD_COUNT_APPLY] > 400000, 1);
            snprintf(what, sizeof(what), "gate %u: %zu nodes, the canonical count", gate, want);
            report(what, nodes != want, 1);
            snprintf(what, sizeof(what), "gate %u: the swap undone gives the state back", gate);
            report(what, OWN(replay_gate(state, name, a, b, RN)) != prev, 1);
        }
        gate++;
        if (gate % 64 == 0) stacks_checkpoint("64 gates of the prefix");
        if (limdd_gc_wanted()) RUN(limdd_gc);
    }
    stacks_checkpoint("the prefix");
    if (gate != 593) bad++;
    report("the prefix: 593 gates", bad, 1);
    printf("  %u diagonal gates, %u of them on a state without full support, "
           "%llu PhaseMul misses\n", diag_gates, diag_partial, (unsigned long long)diag_phasemul);
    report("diagonal gates: no Cof1 and no Apply miss", diag_costly, diag_gates);
    covered("diagonal gates on a state without full support", diag_partial);
    report_stacks();
    limdd_unprotect(&prev);
    limdd_unprotect(&state);
    bqd_quit();
    return failures != before;
}

/* --- Sylvan collections started by any worker --------------------------------- */

/*
 * The exponents of the phase multiplications are Sylvan MTBDDs, so a full
 * Sylvan table starts a collection in the middle of one, from whichever worker
 * fails to claim a region, and every other worker joins it wherever it is.
 * The collections section starts them from one worker, through
 * bqd_exp_collect_every; here nothing does, and a table of 2^12 buckets, 8
 * regions of 512, on 8 workers runs out of regions after a few
 * multiplications: some ten collections a family, started by one to six
 * different workers. The worker that starts one runs the collection's root,
 * and so its hooks.
 */
typedef struct {
    unsigned n;
    char pad[SYLVAN_SHARING_PAD - sizeof(unsigned)];
} count_line_t;

static _Alignas(SYLVAN_SHARING_PAD) count_line_t started_by[MAXW];

VOID_TASK_0(note_collection_start)
{
    started_by[LACE_WORKER_ID].n++;
}

/** f times beta^ev, ev a random exponent, as bqd_phase_mul of its diagram; F protected by the caller */
static BQD
random_phase_mul(BQD F, EVBDD_WGT *f, unsigned n, EVBDD_WGT beta, uint32_t r, int64_t lo, int64_t span)
{
    const uint64_t len = UINT64_C(1) << n;
    int64_t ev[MAXV];
    for (uint64_t x = 0; x < len; x++) ev[x] = lo + (int64_t)rnd_below((uint64_t)span);
    const BQD_EXP eps = mtbdd_refs_push(exp_table(ev, 0, n, r));
    for (uint64_t x = 0; x < len; x++) f[x] = wgt_mul(f[x], bqd_exp_power(beta, ev[x]));
    const BQD res = OWN(bqd_phase_mul(F, eps, beta));
    mtbdd_refs_pop(1);
    return res;
}

TASK_1(int, run_any_worker, int, fam)
{
    enum { ROUNDS = 32, STEPS = 6 };
    bqd_init((bqd_family_t)fam, NQ, 1LL << 20, 1LL << 20, 1LL << 22, 1LL << 20);
    exact_constants();
    rng_state = UINT64_C(0xA11C011EC7) + (uint64_t)fam;
    const int before = failures;
    memset(started_by, 0, sizeof(started_by));
    const EVBDD_WGT gamma = qisq2_lookup(1, 1, 1, 1, 0, 1, 0, 1);      /* 1 + sqrt2 */
    EVBDD_WGT f[MAXV];
    BQD state = limdd_zero_edge();
    limdd_protect(&state);
    tally_t t = { 0, 0 };
    stacks_baseline();
    for (int round = 0; round < ROUNDS; round++) {
        const unsigned n = NQ - (unsigned)rnd_below(3);
        rand_state(n, (unsigned)rnd_below(S_KINDS), f);
        if (rnd_below(2)) twist(n, f);
        state = bqd_from_vector(f, n);
        for (int k = 0; k < STEPS; k++) {
            switch (rnd_below(fam == BQD_FAMILY_PAULI ? 2 : 3)) {
            case 0: {
                const EVBDD_WGT beta = pw[1 + rnd_below(7)];
                state = random_phase_mul(state, f, n, beta, bqd_exp_order(beta), 0, 8);
                break;
            }
            case 1: {
                uint64_t A = 0;
                while (A == 0) A = rnd_below(UINT64_C(1) << n);
                const EVBDD_WGT phase = pw[1 + rnd_below(7)];
                dense_diag(f, n, A, phase, f);
                state = OWN(bqd_apply_diagonal(state, A, phase, n));
                break;
            }
            default:
                state = random_phase_mul(state, f, n, gamma, 0, -2, 5);
                break;
            }
            t.tot++;
            if (state != bqd_from_vector(f, n)) t.bad++;
        }
        stacks_checkpoint("a round of phase multiplications");
    }
    limdd_unprotect(&state);
    unsigned total = 0, starters = 0;
    for (unsigned w = 0; w < lace_workers(); w++) {
        total += started_by[w].n;
        if (started_by[w].n > 0) starters++;
    }
    report("phase multiplications, collections from any worker", t.bad, t.tot);
    covered("Sylvan collections", total);
    printf("  %-60s %u\n", "workers that started one", starters);
    report_stacks();
    bqd_quit();
    return failures != before;
}

/* --- from a thread that is not a worker ----------------------------------------- */

/*
 * The gates of qsylvan_bqd_gates.h are RUN macros, or plain functions that
 * reach the diagram only through RUN macros, so a program calls them from its
 * own thread, outside any task. A plain function that pushes on a worker's
 * reference stack itself reads a stack that thread does not have: the
 * controlled flip with the control below the target and DiagR did, and
 * crashed there. (bqd_apply_diagonal of qsylvan_bqd_ops.h is not among them:
 * it is from a Lace worker, and the Pauli family's walk pushes.) Every family,
 * on five qubits, the results checked as above.
 */
VOID_TASK_1(any_thread_init, int, fam)
{
    bqd_init((bqd_family_t)fam, NQ, 1LL << 18, 1LL << 18, 1LL << 20, 1LL << 18);
    exact_constants();
}

VOID_TASK_0(any_thread_quit)
{
    bqd_quit();
}

static int
run_any_thread(int fam)
{
    enum { N = 5, LEN = 1 << N };
    RUN(any_thread_init, fam);
    rng_state = UINT64_C(0x3A1D7) + (uint64_t)fam;
    const int before = failures;
    tally_t t = { 0, 0 };
    EVBDD_WGT f[LEN], h[LEN];
    for (int rep = 0; rep < 16; rep++) {
        rand_state(N, (unsigned)rnd_below(S_KINDS), f);
        if (rnd_below(2)) twist(N, f);
        const BQD F = bqd_from_vector(f, N);
        const unsigned c = 1 + (unsigned)rnd_below(N - 1), q = (unsigned)rnd_below(c);
        for (int y = 0; y < 2; y++) {                /* the control below the target */
            const uint32_t gid = y ? GATEID_Y : GATEID_X;
            memcpy(h, f, sizeof(h));
            dense_gate(h, N, gates[gid], UINT64_C(1) << c, q);
            bool ok;
            tally(&t, bqd_cgate_either(F, gid, c, q, N, &ok), h, N);
            if (!ok) t.bad++;
        }
        uint64_t A = 0;
        while (A == 0) A = rnd_below(LEN);
        const EVBDD_WGT phase = pw[1 + rnd_below(7)];
        dense_diag(f, N, A, phase, h);
        tally(&t, bqd_diag_any(F, A, phase, N), h, N);
        dense_perm(f, N, BQD_PERM_SWAP, q, c, h);
        tally(&t, bqd_swap(F, q, c, N), h, N);
        memcpy(h, f, sizeof(h));
        dense_gate(h, N, gates[GATEID_Y], 0, q);
        tally(&t, bqd_gate(F, GATEID_Y, q, N), h, N);
        memcpy(h, f, sizeof(h));
        dense_gate(h, N, gates[GATEID_X], UINT64_C(1) << q, c);
        tally(&t, bqd_cgate(F, GATEID_X, UINT64_C(1) << q, c, N), h, N);
    }
    report("gates called from the main thread", t.bad, t.tot);
    RUN(any_thread_quit);
    return failures != before;
}

/* --- harness ------------------------------------------------------------------ */

TASK_1(int, run_family, int, fam)
{
    bqd_init((bqd_family_t)fam, NQ, 1LL << 21, 1LL << 20, 1LL << 22, 1LL << 20);
    exact_constants();
    rng_state = UINT64_C(0x5EED2A7100) + (uint64_t)fam;
    const int before = failures;
    check_family();
    bqd_quit();
    return failures != before;
}

static void
session_begin(unsigned workers, int table, int cache)
{
    lace_start(workers, 0);
    sylvan_set_sizes(1LL << table, 1LL << table, 1LL << cache, 1LL << cache);
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
    setvbuf(stdout, NULL, _IOLBF, 0);
    const unsigned workers = getenv("BQD_RATIO_WORKERS") ? (unsigned)atoi(getenv("BQD_RATIO_WORKERS")) : 4;
    int bad = 0;
    for (int fam = BQD_FAMILY_SCALAR; fam <= BQD_FAMILY_PAULI; fam++) {
        printf("== %s, exact weights, %u workers ==\n", bqd_family_name((bqd_family_t)fam), workers);
        session_begin(workers, 20, 20);
        const int res = RUN(run_family, fam);
        session_end();
        printf("  %s\n", res ? "FAILED" : "ok");
        bad |= res;
    }
    for (int fam = BQD_FAMILY_SCALAR; fam <= BQD_FAMILY_PAULI; fam++) {
        printf("== %s, collections, %u workers ==\n", bqd_family_name((bqd_family_t)fam), workers);
        session_begin(workers, 16, 16);
        sylvan_gc_hook_postgc(count_sylvan_collection_CALL);
        const int res = RUN(run_collections, fam);
        session_end();
        printf("  %s\n", res ? "FAILED" : "ok");
        bad |= res;
    }
    for (int fam = BQD_FAMILY_SCALAR; fam <= BQD_FAMILY_PAULI; fam++) {
        printf("== %s, the swaps of clifford_T_circuit_20_700, %u workers ==\n",
               bqd_family_name((bqd_family_t)fam), workers);
        session_begin(workers, 20, 20);
        const int res = RUN(run_regression, fam);
        session_end();
        printf("  %s\n", res ? "FAILED" : "ok");
        bad |= res;
    }
    for (int fam = BQD_FAMILY_SCALAR; fam <= BQD_FAMILY_PAULI; fam++) {
        printf("== %s, Sylvan collections from any worker, 8 workers ==\n",
               bqd_family_name((bqd_family_t)fam));
        session_begin(8, 12, 16);
        sylvan_gc_hook_pregc(note_collection_start_CALL);
        const int res = RUN(run_any_worker, fam);
        session_end();
        printf("  %s\n", res ? "FAILED" : "ok");
        bad |= res;
    }
    for (int fam = BQD_FAMILY_SCALAR; fam <= BQD_FAMILY_PAULI; fam++) {
        printf("== %s, from the main thread, %u workers ==\n",
               bqd_family_name((bqd_family_t)fam), workers);
        session_begin(workers, 20, 20);
        const int res = run_any_thread(fam);
        session_end();
        printf("  %s\n", res ? "FAILED" : "ok");
        bad |= res;
    }
    return bad;
}
