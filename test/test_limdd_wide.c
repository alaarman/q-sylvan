/*
 * Copyright 2026 System Verification Lab, LIACS, Leiden University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * Tests for LIMDD gates on every qubit a Pauli word can name.
 *
 * Built with LIMDD_PAULI_WORDS=2 the library takes 128 qubits, but the gates
 * kept their controls, and the qubits whose Pauli a gate keeps on the edge it
 * caches, in one 64-bit mask, and limdd_eval read a basis state as one word.
 * A gate with its target or a control at qubit 64 or above then computed a
 * wrong state, with no assertion to catch it.
 *
 * The oracle needs no dense vector of 2^n amplitudes. A random Clifford+T
 * circuit on K qubits is run on qubits 0..K-1 of an n-qubit diagram and again
 * on qubits off..off+K-1 of the same diagram, the other qubits idle in |0>,
 * and the two states must have the same 2^K amplitudes. The weights are
 * exact, so one amplitude is one weight index. Both diagrams are read in one
 * session with no gate in between, so no collection renumbers the weights
 * between the two reads.
 *
 * That catches what changes with the placement, but not a mistake every
 * placement makes alike, such as a control list that drops its third
 * control. So the run on qubits 0..K-1 is also checked against a dense
 * vector of 2^K amplitudes that the same circuit is applied to.
 *
 * At one word the circuit moves onto qubits 56..63, where qubit 63 is the top
 * bit of a mask. At more words it also moves across every boundary between
 * two words, wholly above qubit 64, and onto the last K qubits. Wherever the
 * controls sit below qubit 64 the mask API, limdd_cgate, is checked against
 * the list API on the same circuits, the target on either side of qubit 64.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lace.h>
#include <sylvan.h>

#include "qsylvan.h"
#include "qsylvan_limdd_ops.h"
#include "qsylvan_limdd_gc.h"
#include "test_assert.h"

#define K      8            /* qubits a circuit acts on */
#define NAMP   (1u << K)
#define NOPS   24           /* gates per circuit, the first K of them H */
#define NCIRC  12           /* circuits per placement */
#define NPLACE 8            /* placements per diagram, qubits 0..K-1 first */

static uint64_t rng_state;
static uint64_t
rnd(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

/* --- random circuits ----------------------------------------------------- */

typedef enum {
    OP_GATE,      /* a 1-qubit gate on q[0] */
    OP_EITHER,    /* CX or CZ, control q[0], target q[1], in either order */
    OP_SWAP,      /* q[0] and q[1] */
    OP_CTRL,      /* `gate` on the highest of q[0..nctl], the rest controls */
} op_kind_t;

typedef struct {
    op_kind_t kind;
    uint32_t  gate;
    uint32_t  nctl;
    uint32_t  q[5];         /* distinct, in 0..K-1, in no particular order */
} op_t;

static op_t circuit[NOPS];

static void
random_circuit(uint64_t seed)
{
    static const uint32_t one[] = { GATEID_H, GATEID_S, GATEID_Sdag, GATEID_T,
                                    GATEID_Tdag, GATEID_X, GATEID_Y, GATEID_Z };
    /*
     * The controlled gates of the simulator's gate set and three more: CS,
     * CCZ, CCX, a Z with three controls and an X with four, the most a
     * control list holds. CS and CCZ are symmetric, so taking the highest
     * qubit as the target is what the simulator does with whichever qubit
     * the circuit called the target.
     */
    static const struct { uint32_t gate, nctl; } ctrl[] = {
        { GATEID_S, 1 }, { GATEID_Z, 2 }, { GATEID_X, 2 }, { GATEID_Z, 3 },
        { GATEID_X, 4 },
    };

    rng_state = seed;
    for (int i = 0; i < NOPS; i++) {
        op_t o;
        memset(&o, 0, sizeof(o));
        for (int j = 0; j < 5; j++) {
            bool fresh;
            do {
                o.q[j] = (uint32_t)(rnd() % K);
                fresh = true;
                for (int l = 0; l < j; l++) fresh &= (o.q[l] != o.q[j]);
            } while (!fresh);
        }
        const uint64_t r = rnd() % 17;
        if (r < 9) {
            o.kind = OP_GATE;
            o.gate = one[rnd() % 8];
        } else if (r < 11) {
            o.kind = OP_EITHER;
            o.gate = (rnd() & 1) ? GATEID_Z : GATEID_X;
        } else if (r < 12) {
            o.kind = OP_SWAP;
        } else {
            o.kind = OP_CTRL;
            o.gate = ctrl[r - 12].gate;
            o.nctl = ctrl[r - 12].nctl;
        }
        if (i < K) {
            /* H on every qubit first. A wrong amplitude only shows where the
             * right one is not zero, and after a few H gates most are zero. */
            o.kind = OP_GATE;
            o.gate = GATEID_H;
            o.q[0] = (uint32_t)i;
        }
        circuit[i] = o;
    }
}

/**
 * The target of an OP_CTRL on qubits `q`, the highest of them, with its
 * controls in `ctl` in the circuit's order. Returns the number of controls.
 */
static uint32_t
ctrl_split(const op_t *o, const uint32_t *q, uint32_t *target, uint32_t *ctl)
{
    uint32_t top = 0, k = 0;
    for (uint32_t j = 1; j <= o->nctl; j++) if (q[j] > q[top]) top = j;
    for (uint32_t j = 0; j <= o->nctl; j++) if (j != top) ctl[k++] = q[j];
    *target = q[top];
    return k;
}

/* Gates sent through limdd_cgate with the target at qubit 64 or above. */
static unsigned mask_wide;

/**
 * The circuit on qubits off..off+K-1 of the all-zero state on n qubits. With
 * `mask`, every gate whose controls sit above its target and below qubit 64
 * goes through limdd_cgate and a control mask. The others go through
 * limdd_cgate_list, or limdd_cgate_either and limdd_swap, which use the list.
 */
static LIMDD
run(uint32_t n, uint32_t off, bool mask)
{
    LIMDD e = limdd_all_zero_state(n);
    limdd_protect(&e);
    for (int i = 0; i < NOPS; i++) {
        const op_t *o = &circuit[i];
        uint32_t q[5];
        for (int j = 0; j < 5; j++) q[j] = o->q[j] + off;

        bool ok = true;
        switch (o->kind) {
        case OP_GATE:
            e = limdd_gate(e, o->gate, q[0], n);
            break;
        case OP_EITHER:
            if (mask && q[0] < q[1] && q[0] < 64) {
                mask_wide += (q[1] >= 64);
                e = limdd_cgate(e, o->gate, UINT64_C(1) << q[0], q[1], n);
            } else {
                e = limdd_cgate_either(e, o->gate, q[0], q[1], n, &ok);
            }
            break;
        case OP_SWAP:
            e = limdd_swap(e, q[0], q[1], n);
            break;
        case OP_CTRL: {
            /* The controls go to limdd_control_list in the circuit's order,
             * which it sorts. */
            uint32_t t, ctl[4];
            const uint32_t k = ctrl_split(o, q, &t, ctl);
            bool low = true;
            for (uint32_t j = 0; j < k; j++) low &= (ctl[j] < 64);
            if (mask && low) {
                uint64_t m = 0;
                for (uint32_t j = 0; j < k; j++) m |= UINT64_C(1) << ctl[j];
                mask_wide += (t >= 64);
                e = limdd_cgate(e, o->gate, m, t, n);
            } else {
                e = limdd_cgate_list(e, o->gate, limdd_control_list(ctl, k), t, n);
            }
            break;
        }
        }
        if (!ok) {
            fprintf(stderr, "limdd_cgate_either refused gate %u on %u, %u\n",
                    o->gate, q[0], q[1]);
            exit(1);
        }
    }
    limdd_unprotect(&e);
    return e;
}

/* --- the dense reference ------------------------------------------------- */

/*
 * The gate matrices come out of Q-Sylvan's gate table, as in test_limdd_ops,
 * so the dense vector and the diagram cannot disagree about what gate was
 * meant. Bit j of an index is qubit j of the circuit.
 */
typedef struct { double re, im; } cx;

static cx
cx_mul_add(cx a, cx b, cx c, cx d)       /* a*b + c*d */
{
    cx r = { a.re*b.re - a.im*b.im + c.re*d.re - c.im*d.im,
             a.re*b.im + a.im*b.re + c.re*d.im + c.im*d.re };
    return r;
}

static cx
wgt_cx(EVBDD_WGT w)
{
    const complex_t c = weight_as_complex(w);
    cx r = { c.r, c.i };
    return r;
}

static void
dense_gate(cx *v, uint32_t gateid, uint32_t t, uint32_t controls)
{
    const cx u00 = wgt_cx(gates[gateid][0]), u01 = wgt_cx(gates[gateid][1]);
    const cx u10 = wgt_cx(gates[gateid][2]), u11 = wgt_cx(gates[gateid][3]);
    const uint32_t bit = 1u << t;
    for (uint32_t i = 0; i < NAMP; i++) {
        if ((i & bit) || (i & controls) != controls) continue;
        const cx a = v[i], b = v[i | bit];
        v[i]       = cx_mul_add(u00, a, u01, b);
        v[i | bit] = cx_mul_add(u10, a, u11, b);
    }
}

static void
dense_run(cx *v)
{
    memset(v, 0, NAMP * sizeof(cx));
    v[0].re = 1.0;
    for (int i = 0; i < NOPS; i++) {
        const op_t *o = &circuit[i];
        switch (o->kind) {
        case OP_GATE:
            dense_gate(v, o->gate, o->q[0], 0);
            break;
        case OP_EITHER:
            dense_gate(v, o->gate, o->q[1], 1u << o->q[0]);
            break;
        case OP_SWAP: {
            const uint32_t a = 1u << o->q[0], b = 1u << o->q[1];
            for (uint32_t x = 0; x < NAMP; x++) {
                if ((x & a) && !(x & b)) {
                    const cx s = v[x]; v[x] = v[x ^ a ^ b]; v[x ^ a ^ b] = s;
                }
            }
            break;
        }
        case OP_CTRL: {
            uint32_t t, ctl[4], m = 0;
            const uint32_t k = ctrl_split(o, o->q, &t, ctl);
            for (uint32_t j = 0; j < k; j++) m |= 1u << ctl[j];
            dense_gate(v, o->gate, t, m);
            break;
        }
        }
    }
}

/* --- comparing amplitudes ------------------------------------------------- */

/* The 2^K amplitudes with qubits off..off+K-1 running through every value,
 * bit j of the index being qubit off+j, and every other qubit 0. */
static void
amplitudes(LIMDD e, uint32_t n, uint32_t off, EVBDD_WGT *out)
{
    bool bits[LIMDD_MAX_QUBITS];
    for (uint32_t y = 0; y < NAMP; y++) {
        for (uint32_t q = 0; q < n; q++) bits[q] = false;
        for (uint32_t j = 0; j < K; j++) bits[off + j] = (y >> j) & 1;
        out[y] = limdd_eval(e, bits, n);
    }
}

static unsigned
nonzero(LIMDD e, uint32_t n)
{
    EVBDD_WGT x[NAMP];
    amplitudes(e, n, 0, x);
    unsigned k = 0;
    for (uint32_t i = 0; i < NAMP; i++) k += (x[i] != EVBDD_ZERO);
    return k;
}

/** Counts the amplitudes of `a` at a_off that differ from those of `b` at b_off. */
static unsigned
differ(LIMDD a, uint32_t a_off, LIMDD b, uint32_t b_off, uint32_t n, int circ)
{
    EVBDD_WGT x[NAMP], y[NAMP];
    amplitudes(a, n, a_off, x);
    amplitudes(b, n, b_off, y);
    unsigned bad = 0;
    for (uint32_t i = 0; i < NAMP; i++) {
        if (x[i] == y[i]) continue;
        if (bad++ == 0) {
            const complex_t cx = weight_as_complex(x[i]), cy = weight_as_complex(y[i]);
            fprintf(stderr, "  circuit %d, amplitude %u: (%g,%g) on qubits %u.., "
                            "(%g,%g) on qubits %u..\n",
                    circ, i, cx.r, cx.i, a_off, cy.r, cy.i, b_off);
        }
    }
    return bad;
}

/** Counts the amplitudes of `e` at `off` that are not those of `want`. */
static unsigned
differ_dense(LIMDD e, uint32_t n, uint32_t off, const cx *want, int circ,
             const char *what)
{
    EVBDD_WGT x[NAMP];
    amplitudes(e, n, off, x);
    unsigned bad = 0;
    for (uint32_t i = 0; i < NAMP; i++) {
        const cx got = wgt_cx(x[i]);
        if (fabs(got.re - want[i].re) < 1e-9 && fabs(got.im - want[i].im) < 1e-9) continue;
        if (bad++ == 0) {
            fprintf(stderr, "  circuit %d, %s on qubits %u.., amplitude %u: (%g,%g), "
                            "dense (%g,%g)\n",
                    circ, what, off, i, got.re, got.im, want[i].re, want[i].im);
        }
    }
    return bad;
}

/*
 * Controlled Z gates on one target with one first control, which differ in
 * the controls after the first. Each is applied to the same state, so a gate
 * cache whose key lost the controls after the first would hand the result of
 * one to the next. The qubits are those of the circuit, in any order.
 */
#define PROBE_TARGET 7
static const struct { uint32_t k, c[4]; } probes[] = {
    { 1, { 3 } }, { 2, { 3, 4 } }, { 2, { 5, 3 } }, { 3, { 6, 3, 4 } },
    { 4, { 3, 4, 5, 6 } },
};
#define NPROBES (sizeof(probes) / sizeof(probes[0]))

/** The probes on `e` at `off`, against the same gates on `want`. */
static unsigned
test_probes(LIMDD e, uint32_t n, uint32_t off, const cx *want, int circ)
{
    unsigned bad = 0;
    for (size_t p = 0; p < NPROBES; p++) {
        uint32_t c[4], m = 0;
        for (uint32_t j = 0; j < probes[p].k; j++) {
            c[j] = probes[p].c[j] + off;
            m |= 1u << probes[p].c[j];
        }
        const LIMDD r = limdd_cgate_list(e, GATEID_Z, limdd_control_list(c, probes[p].k),
                                         PROBE_TARGET + off, n);
        cx w[NAMP];
        memcpy(w, want, sizeof(w));
        dense_gate(w, GATEID_Z, PROBE_TARGET, m);
        char what[32];
        snprintf(what, sizeof(what), "probe %zu", p);
        bad += differ_dense(r, n, off, w, circ, what);
    }
    return bad;
}

/* --- the tests ------------------------------------------------------------ */

/**
 * Every circuit on qubits 0..K-1 of one n-qubit diagram, and on qubits
 * off..off+K-1 of it for each offset. At every placement below qubit 64 the
 * circuit also goes through the mask API where its controls let it.
 */
static int
test_placements(uint32_t n, const uint32_t *offs, int noffs)
{
    test_assert(noffs < NPLACE);
    limdd_nodes_init(n, 1LL << 20, 1LL << 20, 1LL << 20, 1LL << 18);

    /* Placement 0 is qubits 0..K-1, placement i is offs[i-1] onwards. */
    unsigned bad[NPLACE] = { 0 }, bad_mask[NPLACE] = { 0 }, bad_probe[NPLACE] = { 0 };
    unsigned bad_dense = 0, live = 0;
    mask_wide = 0;
    for (int c = 0; c < NCIRC; c++) {
        random_circuit(UINT64_C(0x9e3779b97f4a7c15) + (uint64_t)c);
        cx want[NAMP];
        dense_run(want);
        LIMDD ref = run(n, 0, false);
        limdd_protect(&ref);
        live += nonzero(ref, n);
        bad_dense += differ_dense(ref, n, 0, want, c, "circuit");
        for (int i = 0; i <= noffs; i++) {
            const uint32_t off = (i == 0) ? 0 : offs[i - 1];
            LIMDD e = (i == 0) ? ref : run(n, off, false);
            limdd_protect(&e);
            if (i > 0) bad[i] += differ(ref, 0, e, off, n, c);
            bad_probe[i] += test_probes(e, n, off, want, c);
            if (off < 64) {
                const LIMDD m = run(n, off, true);
                bad_mask[i] += differ(e, off, m, off, n, c);
            }
            limdd_unprotect(&e);
        }
        limdd_unprotect(&ref);
    }
    limdd_nodes_quit();

    /* Two placements agree trivially where both are zero, so the circuits
     * have to leave enough amplitudes that are not. */
    const unsigned all = NCIRC * NAMP;
    printf("  %3u qubits, qubits   0..%3u: %4u of %u amplitudes not zero\n",
           n, K - 1, live, all);
    printf("  %3u qubits, qubits   0..%3u against a dense vector:  %4u of %u differ\n",
           n, K - 1, bad_dense, all);
    int res = (live < all / 2) | (bad_dense != 0);
    for (int i = 0; i <= noffs; i++) {
        const uint32_t off = (i == 0) ? 0 : offs[i - 1];
        if (i > 0)
            printf("  %3u qubits, qubits %3u..%3u against qubits 0..%u: %4u of %u differ\n",
                   n, off, off + K - 1, K - 1, bad[i], all);
        printf("  %3u qubits, qubits %3u..%3u, %zu gates on one state:  %4u of %zu differ\n",
               n, off, off + K - 1, NPROBES, bad_probe[i], NPROBES * all);
        if (off < 64)
            printf("  %3u qubits, qubits %3u..%3u, mask against list:  %4u of %u differ\n",
                   n, off, off + K - 1, bad_mask[i], all);
        res |= (bad[i] != 0) | (bad_mask[i] != 0) | (bad_probe[i] != 0);
    }
    /* A placement across qubit 64 must have sent some gate with its controls
     * below 64 and its target above through the mask. */
    if (n > 64) {
        printf("  %3u qubits, gates through a mask with the target at 64 or above: %u\n",
               n, mask_wide);
        res |= (mask_wide == 0);
    }
    return res;
}

/** The list does not depend on the order of the qubits, nor on a repeat. */
static int
test_control_list(void)
{
    const uint32_t hi = LIMDD_MAX_QUBITS - 1;
    const uint32_t a[] = { 2, 5 }, b[] = { 5, 2 }, c[] = { 5, 2, 5, 2, 2 };
    const uint32_t d[] = { 2 }, e[] = { 5, 2, 6 }, f[] = { 7, 6, 5, 2 };
    const uint32_t g[] = { hi, 60, hi }, h[] = { 60, hi };
    test_assert(limdd_control_list(NULL, 0) == 0);
    test_assert(limdd_control_list(a, 2) == limdd_control_list(b, 2));
    test_assert(limdd_control_list(a, 2) == limdd_control_list(c, 5));
    test_assert(limdd_control_list(a, 2) != limdd_control_list(d, 1));
    test_assert(limdd_control_list(a, 2) != limdd_control_list(e, 3));
    test_assert(limdd_control_list(e, 3) != limdd_control_list(f, 4));
    test_assert(limdd_control_list(g, 3) == limdd_control_list(h, 2));
    test_assert(limdd_control_list(g, 1) != limdd_control_list(h, 1));
    printf("  control lists: in any order, a repeat counts once: ok\n");
    return 0;
}

TASK_0(int, runtests)
{
    if (test_control_list()) return 1;

    /* 56..63: the top bit of a mask, as far as one word goes. */
    const uint32_t one_word[] = { 56 };
    if (test_placements(64, one_word, 1)) return 1;

    if (LIMDD_PAULI_WORDS >= 2) {
        /* Across every boundary between two words, wholly above the first
         * word, and on the last K qubits. */
        const uint32_t n = LIMDD_MAX_QUBITS;
        uint32_t offs[NPLACE - 1];
        int k = 0;
        for (uint32_t w = 1; w < LIMDD_PAULI_WORDS && k < NPLACE - 3; w++) offs[k++] = 64 * w - 4;
        offs[k++] = 70;
        offs[k++] = n - K;
        if (test_placements(n, offs, k)) return 1;
    }
    return 0;
}

int
main(void)
{
    printf("== LIMDD gates on every qubit, %d Pauli word%s, exact weights ==\n",
           LIMDD_PAULI_WORDS, LIMDD_PAULI_WORDS == 1 ? "" : "s");

    lace_start(1, 0);
    sylvan_set_sizes(1LL << 20, 1LL << 20, 1LL << 20, 1LL << 20);
    sylvan_init_package();
    /* Through the simulator's initialiser, so the gate table is populated.
     * Exact weights: see the note at the top. */
    qsylvan_init_simulator(1LL << 20, 1LL << 20, -1, QISQ2_MAP, NORM_LOW);

    const int res = RUN(runtests);

    sylvan_quit();
    lace_stop();
    return res;
}
