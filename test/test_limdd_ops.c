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
 * Tests for LIMDD addition and gate application.
 *
 * Everything here is checked against a dense state vector carried alongside
 * the diagram: the same gate is applied to both and all 2^n amplitudes are
 * compared. A diagram operation can be wrong in ways no structural check would
 * notice -- a branch swapped, a phase dropped, a control read off the wrong
 * qubit -- and all of them change an amplitude, so that is what is compared.
 *
 * The gate matrices come out of Q-Sylvan's own gate table rather than being
 * written out here, so the dense simulator and the diagram cannot disagree
 * about what gate was meant.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lace.h>
#include <sylvan.h>

#include "qsylvan.h"
#include "qsylvan_limdd_ops.h"
#include "sylvan_edge_weights_complex.h"
#include "sylvan_edge_weights_qisq2.h"
#include "test_assert.h"

#define NQUBITS 5
#define NBASIS  (1u << NQUBITS)

static uint64_t rng_state = UINT64_C(0x1234567890ABCDEF);
static uint64_t
rnd(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

typedef struct { double re, im; } cx;

static cx cx_add(cx a, cx b) { cx r = { a.re+b.re, a.im+b.im }; return r; }
static cx cx_mul(cx a, cx b)
{
    cx r = { a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re };
    return r;
}
static bool cx_eq(cx a, cx b)
{
    return fabs(a.re-b.re) < 1e-9 && fabs(a.im-b.im) < 1e-9;
}

static cx
wgt_cx(EVBDD_WGT w)
{
    const complex_t c = weight_as_complex(w);
    cx r = { c.r, c.i };
    return r;
}

/* --- the dense reference ------------------------------------------------- */

static void
dense_gate(cx *v, uint32_t gateid, uint32_t q, uint64_t controls)
{
    const cx u00 = wgt_cx(gates[gateid][0]), u01 = wgt_cx(gates[gateid][1]);
    const cx u10 = wgt_cx(gates[gateid][2]), u11 = wgt_cx(gates[gateid][3]);
    const unsigned bit = 1u << q;

    for (unsigned i = 0; i < NBASIS; i++) {
        if (i & bit) continue;
        if ((i & controls) != controls) continue;   /* controls must all be 1 */
        const unsigned j = i | bit;
        const cx a = v[i], b = v[j];
        v[i] = cx_add(cx_mul(u00, a), cx_mul(u01, b));
        v[j] = cx_add(cx_mul(u10, a), cx_mul(u11, b));
    }
}

static void
limdd_to_vector(LIMDD e, cx *out)
{
    bool bits[NQUBITS];
    for (unsigned i = 0; i < NBASIS; i++) {
        for (int k = 0; k < NQUBITS; k++) bits[k] = (i >> k) & 1;
        out[i] = wgt_cx(limdd_eval(e, bits, NQUBITS));
    }
}

static int
compare(LIMDD e, const cx *want, const char *what, int step)
{
    cx got[NBASIS];
    limdd_to_vector(e, got);
    for (unsigned i = 0; i < NBASIS; i++) {
        if (!cx_eq(got[i], want[i])) {
            fprintf(stderr, "%s step %d: amplitude %u is (%g,%g), should be (%g,%g)\n",
                    what, step, i, got[i].re, got[i].im, want[i].re, want[i].im);
            return 1;
        }
    }
    return 0;
}

/* --- tests ---------------------------------------------------------------- */

int
test_basis_state(void)
{
    LIMDD e = limdd_all_zero_state(NQUBITS);
    cx want[NBASIS];
    memset(want, 0, sizeof(want));
    want[0].re = 1.0;
    return compare(e, want, "all-zero state", 0);
}

int
test_single_gates(void)
{
    /* Each gate on each qubit, from |0..0>, against the dense reference. */
    const uint32_t ids[] = { GATEID_X, GATEID_Y, GATEID_Z, GATEID_H,
                             GATEID_S, GATEID_T, GATEID_Sdag, GATEID_Tdag };
    for (unsigned g = 0; g < sizeof(ids)/sizeof(ids[0]); g++) {
        for (uint32_t q = 0; q < NQUBITS; q++) {
            LIMDD e = limdd_all_zero_state(NQUBITS);
            cx want[NBASIS];
            memset(want, 0, sizeof(want));
            want[0].re = 1.0;

            e = limdd_gate(e, ids[g], q, NQUBITS);
            dense_gate(want, ids[g], q, 0);

            char label[64];
            snprintf(label, sizeof(label), "gate %u on qubit %u", ids[g], q);
            if (compare(e, want, label, 0)) return 1;
        }
    }
    return 0;
}

int
test_controlled_gates(void)
{
    /* Controls only bite once the control qubit is in superposition, so each
     * case starts by putting H on the controls. */
    /* Controls above targets only; see limdd_cgate. */
    for (uint32_t c = 0; c < NQUBITS; c++) {
        for (uint32_t t = c + 1; t < NQUBITS; t++) {
            LIMDD e = limdd_all_zero_state(NQUBITS);
            cx want[NBASIS];
            memset(want, 0, sizeof(want));
            want[0].re = 1.0;

            e = limdd_gate(e, GATEID_H, c, NQUBITS);
            dense_gate(want, GATEID_H, c, 0);

            e = limdd_cgate(e, GATEID_X, UINT64_C(1) << c, t, NQUBITS);
            dense_gate(want, GATEID_X, t, UINT64_C(1) << c);

            char label[64];
            snprintf(label, sizeof(label), "cnot %u->%u", c, t);
            if (compare(e, want, label, 0)) return 1;
        }
    }

    /* Two controls, which is what Grover's diffusion needs. */
    LIMDD e = limdd_all_zero_state(NQUBITS);
    cx want[NBASIS];
    memset(want, 0, sizeof(want));
    want[0].re = 1.0;
    for (uint32_t q = 0; q < 3; q++) {
        e = limdd_gate(e, GATEID_H, q, NQUBITS);
        dense_gate(want, GATEID_H, q, 0);
    }
    const uint64_t ctrl = 0x3;   /* qubits 0 and 1 */
    e = limdd_cgate(e, GATEID_Z, ctrl, 2, NQUBITS);
    dense_gate(want, GATEID_Z, 2, ctrl);
    return compare(e, want, "ccz", 0);
}

/**
 * SWAP, and controlled gates with the control below the target.
 *
 * These go through identities rather than the recursion, so they need
 * checking against the dense reference in their own right -- a wrong
 * Hadamard conjugation would still produce a valid state, just not this one.
 */
int
test_reversed_and_swap(void)
{
    for (uint32_t a = 0; a < NQUBITS; a++) {
        for (uint32_t b = 0; b < NQUBITS; b++) {
            if (a == b) continue;

            LIMDD e = limdd_all_zero_state(NQUBITS);
            cx want[NBASIS];
            memset(want, 0, sizeof(want));
            want[0].re = 1.0;

            /* Spread the amplitude first, or every gate below is a no-op. */
            for (uint32_t q = 0; q < NQUBITS; q++) {
                const uint32_t g = (q % 2) ? GATEID_H : GATEID_T;
                e = limdd_gate(e, g, q, NQUBITS);
                dense_gate(want, g, q, 0);
            }

            bool ok;
            e = limdd_cgate_either(e, GATEID_X, a, b, NQUBITS, &ok);
            test_assert(ok);
            dense_gate(want, GATEID_X, b, UINT64_C(1) << a);
            char l1[48]; snprintf(l1, sizeof(l1), "cx %u->%u", a, b);
            if (compare(e, want, l1, 0)) return 1;

            e = limdd_cgate_either(e, GATEID_Z, a, b, NQUBITS, &ok);
            test_assert(ok);
            dense_gate(want, GATEID_Z, b, UINT64_C(1) << a);
            char l2[48]; snprintf(l2, sizeof(l2), "cz %u->%u", a, b);
            if (compare(e, want, l2, 0)) return 1;

            e = limdd_swap(e, a, b, NQUBITS);
            /* SWAP as three CNOTs, applied to the reference the same way. */
            dense_gate(want, GATEID_X, b, UINT64_C(1) << a);
            dense_gate(want, GATEID_X, a, UINT64_C(1) << b);
            dense_gate(want, GATEID_X, b, UINT64_C(1) << a);
            char l3[48]; snprintf(l3, sizeof(l3), "swap %u<->%u", a, b);
            if (compare(e, want, l3, 0)) return 1;
        }
    }
    return 0;
}

int
test_random_circuits(void)
{
    /*
     * Random Clifford+T circuits, compared after every gate. Checking only at
     * the end would find the first mistake but say nothing about where, and
     * errors that cancel would go unnoticed.
     */
    const uint32_t ids[] = { GATEID_X, GATEID_Y, GATEID_Z, GATEID_H,
                             GATEID_S, GATEID_T };

    for (int trial = 0; trial < 40; trial++) {
        LIMDD e = limdd_all_zero_state(NQUBITS);
        cx want[NBASIS];
        memset(want, 0, sizeof(want));
        want[0].re = 1.0;

        for (int step = 0; step < 40; step++) {
            char label[48];
            snprintf(label, sizeof(label), "random trial %d", trial);

            if ((rnd() % 3) == 0) {
                const uint32_t c = rnd() % (NQUBITS - 1);
                const uint32_t t = c + 1 + (rnd() % (NQUBITS - 1 - c));
                const uint32_t g = ids[rnd() % 4];   /* X, Y, Z, H */
                e = limdd_cgate(e, g, UINT64_C(1) << c, t, NQUBITS);
                dense_gate(want, g, t, UINT64_C(1) << c);
            } else {
                const uint32_t q = rnd() % NQUBITS;
                const uint32_t g = ids[rnd() % 6];
                e = limdd_gate(e, g, q, NQUBITS);
                dense_gate(want, g, q, 0);
            }
            if (compare(e, want, label, step)) return 1;
        }
    }
    return 0;
}

int
test_addition(void)
{
    /* Sums of random reachable states, against the dense sum. */
    for (int trial = 0; trial < 200; trial++) {
        LIMDD a = limdd_all_zero_state(NQUBITS);
        LIMDD b = limdd_all_zero_state(NQUBITS);
        cx va[NBASIS], vb[NBASIS];
        memset(va, 0, sizeof(va)); va[0].re = 1.0;
        memset(vb, 0, sizeof(vb)); vb[0].re = 1.0;

        const uint32_t ids[] = { GATEID_X, GATEID_Z, GATEID_H, GATEID_S };
        for (int step = 0; step < 6; step++) {
            uint32_t q = rnd() % NQUBITS, g = ids[rnd() % 4];
            a = limdd_gate(a, g, q, NQUBITS);  dense_gate(va, g, q, 0);
            if (compare(a, va, "addition/build a", step)) {
                fprintf(stderr, "  (trial %d, gate %u on qubit %u)\n", trial, g, q);
                return 1;
            }
            q = rnd() % NQUBITS; g = ids[rnd() % 4];
            b = limdd_gate(b, g, q, NQUBITS);  dense_gate(vb, g, q, 0);
            if (compare(b, vb, "addition/build b", step)) {
                fprintf(stderr, "  (trial %d, gate %u on qubit %u)\n", trial, g, q);
                return 1;
            }
        }

        const LIMDD s = limdd_plus(a, b, 0);
        cx want[NBASIS];
        for (unsigned i = 0; i < NBASIS; i++) want[i] = cx_add(va[i], vb[i]);
        if (compare(s, want, "addition", trial)) {
            cx ca[NBASIS], cb[NBASIS], cs[NBASIS];
            limdd_to_vector(a, ca); limdd_to_vector(b, cb); limdd_to_vector(s, cs);
            fprintf(stderr, "  a == b ? %s   a=%llu b=%llu s=%llu\n",
                    a == b ? "yes" : "no",
                    (unsigned long long)a, (unsigned long long)b, (unsigned long long)s);
            for (unsigned i = 0; i < NBASIS; i++) {
                fprintf(stderr, "  [%2u] a=(%6.3f,%6.3f) dense_a=(%6.3f,%6.3f) | "
                                "b=(%6.3f,%6.3f) dense_b=(%6.3f,%6.3f) | sum=(%6.3f,%6.3f) want=(%6.3f,%6.3f)\n",
                        i, ca[i].re, ca[i].im, va[i].re, va[i].im,
                        cb[i].re, cb[i].im, vb[i].re, vb[i].im,
                        cs[i].re, cs[i].im, want[i].re, want[i].im);
            }
            return 1;
        }
    }
    return 0;
}

int
test_probabilities(void)
{
    for (int trial = 0; trial < 100; trial++) {
        LIMDD e = limdd_all_zero_state(NQUBITS);
        cx want[NBASIS];
        memset(want, 0, sizeof(want)); want[0].re = 1.0;

        const uint32_t ids[] = { GATEID_X, GATEID_Z, GATEID_H, GATEID_S, GATEID_T };
        for (int step = 0; step < 12; step++) {
            if ((rnd() % 3) == 0) {
                const uint32_t c = rnd() % (NQUBITS - 1);
                const uint32_t t = c + 1 + (rnd() % (NQUBITS - 1 - c));
                e = limdd_cgate(e, GATEID_X, UINT64_C(1) << c, t, NQUBITS);
                dense_gate(want, GATEID_X, t, UINT64_C(1) << c);
            } else {
                const uint32_t q = rnd() % NQUBITS, g = ids[rnd() % 5];
                e = limdd_gate(e, g, q, NQUBITS);
                dense_gate(want, g, q, 0);
            }
        }

        /* A unitary circuit from a basis state stays normalised, and the
         * diagram must agree without expanding it. */
        double n2 = 0.0;
        for (unsigned i = 0; i < NBASIS; i++) n2 += want[i].re*want[i].re + want[i].im*want[i].im;
        if (fabs(limdd_norm_squared(e, NQUBITS) - n2) > 1e-9) {
            fprintf(stderr, "trial %d: norm^2 is %g, dense says %g\n",
                    trial, limdd_norm_squared(e, NQUBITS), n2);
            return 1;
        }

        for (uint32_t q = 0; q < NQUBITS; q++) {
            double p1 = 0.0;
            for (unsigned i = 0; i < NBASIS; i++) {
                if (i & (1u << q)) p1 += want[i].re*want[i].re + want[i].im*want[i].im;
            }
            const double got = limdd_prob_qubit_one(e, q, NQUBITS);
            if (fabs(got - p1) > 1e-9) {
                fprintf(stderr, "trial %d qubit %u: P(1) is %g, dense says %g\n",
                        trial, q, got, p1);
                return 1;
            }
        }
    }
    return 0;
}

/* --- operations on diagrams that skip levels ----------------------------- */

/** A random scalar, the same five values under either backend. */
static EVBDD_WGT
random_scalar(void)
{
    const unsigned w = rnd() % 5;
    if (sylvan_get_edge_weight_type() == WGT_QISQ2) {
        switch (w) {
        case 0: return qisq2_lookup( 1,1, 0,1,  0,1, 0,1);   //  1
        case 1: return qisq2_lookup( 0,1, 0,1,  1,1, 0,1);   //  i
        case 2: return qisq2_lookup(-1,1, 0,1,  0,1, 0,1);   // -1
        case 3: return qisq2_lookup( 1,2, 0,1,  1,2, 0,1);   //  (1+i)/2
        default:return qisq2_lookup( 0,1, 0,1, -1,4, 0,1);   // -i/4
        }
    }
    static const double re[5] = { 1.0, 0.0, -1.0, 0.5, 0.0 };
    static const double im[5] = { 0.0, 1.0,  0.0, 0.5, -0.25 };
    return complex_lookup(re[w], im[w]);
}

/** A random label acting on qubits `level`..n-1 only. */
static LIMDD_LIM
random_label_above(uint32_t level)
{
    const uint64_t mask = ((UINT64_C(1) << NQUBITS) - 1) & ~((UINT64_C(1) << level) - 1);
    const limdd_pauli_t p = { rnd() & mask, rnd() & mask };
    return limdd_lim_make(p, random_scalar());
}

/**
 * A random edge to be read at `level`, built with makenode so that it may
 * skip levels: it lands on a node anywhere from `level` down, or on the
 * terminal, and every edge in it carries a random Pauli including on the
 * levels it skips.
 */
static LIMDD
random_skipping_edge(uint32_t level)
{
    if (level == NQUBITS || (rnd() % 4) == 0) {
        return limdd_bundle(random_label_above(level), LIMDD_TERMINAL);
    }
    const uint32_t k = level + (uint32_t)(rnd() % (NQUBITS - level));
    LIMDD lo = random_skipping_edge(k + 1);
    lo = limdd_bundle(LIMDD_LIM_IDENTITY, limdd_target(lo));   /* low carries no label */
    LIMDD hi = (rnd() % 6) == 0 ? limdd_zero_edge() : random_skipping_edge(k + 1);
    const LIMDD_TARG t = limdd_makenode(k, lo, hi);
    return limdd_bundle(random_label_above(level), t);
}

int
test_ops_on_skipping_diagrams(void)
{
    /*
     * Nothing builds a skipping diagram yet, so they are built by hand; the
     * operations must then treat them as what they denote. Each result goes
     * through eval, which test_limdd_node checks against the semantics.
     */
    const uint32_t ids[] = { GATEID_X, GATEID_Y, GATEID_Z, GATEID_H, GATEID_S, GATEID_T };

    for (int trial = 0; trial < 300; trial++) {
        const LIMDD e = random_skipping_edge(0);
        cx v[NBASIS];
        limdd_to_vector(e, v);

        /* every single-qubit gate on every qubit, skipped or not */
        for (uint32_t q = 0; q < NQUBITS; q++) {
            const uint32_t g = ids[rnd() % 6];
            cx want[NBASIS];
            memcpy(want, v, sizeof(want));
            dense_gate(want, g, q, 0);
            if (compare(limdd_gate(e, g, q, NQUBITS), want, "skipping/gate", trial)) {
                fprintf(stderr, "  (gate %u on qubit %u, root level %u)\n",
                        g, q, limdd_level(limdd_target(e)));
                return 1;
            }
        }

        /* a CNOT with random control above target */
        {
            const uint32_t c = rnd() % (NQUBITS - 1);
            const uint32_t t = c + 1 + (rnd() % (NQUBITS - 1 - c));
            cx want[NBASIS];
            memcpy(want, v, sizeof(want));
            dense_gate(want, GATEID_X, t, UINT64_C(1) << c);
            if (compare(limdd_cgate(e, GATEID_X, UINT64_C(1) << c, t, NQUBITS),
                        want, "skipping/cnot", trial)) {
                fprintf(stderr, "  (control %u, target %u)\n", c, t);
                return 1;
            }
        }

        /* swap, which routes through cx_reversed and so through cgate twice */
        {
            const uint32_t a = rnd() % NQUBITS;
            const uint32_t b = rnd() % NQUBITS;
            cx want[NBASIS];
            for (unsigned i = 0; i < NBASIS; i++) {
                const unsigned ba = (i >> a) & 1, bb = (i >> b) & 1;
                unsigned j = i & ~((1u << a) | (1u << b));
                j |= bb << a; j |= ba << b;
                want[j] = v[i];
            }
            if (compare(limdd_swap(e, a, b, NQUBITS), want, "skipping/swap", trial)) {
                fprintf(stderr, "  (swap %u <-> %u, root level %u)\n",
                        a, b, limdd_level(limdd_target(e)));
                return 1;
            }
        }

        /* addition with another skipping diagram */
        {
            const LIMDD f = random_skipping_edge(0);
            cx vf[NBASIS], want[NBASIS];
            limdd_to_vector(f, vf);
            for (unsigned i = 0; i < NBASIS; i++) want[i] = cx_add(v[i], vf[i]);
            if (compare(limdd_plus(e, f, 0), want, "skipping/plus", trial)) return 1;
        }

        /* norm and per-qubit probabilities, both unnormalised */
        {
            double n2 = 0.0;
            for (unsigned i = 0; i < NBASIS; i++) n2 += v[i].re*v[i].re + v[i].im*v[i].im;
            const double got = limdd_norm_squared(e, NQUBITS);
            if (fabs(got - n2) > 1e-9 * (1.0 + n2)) {
                fprintf(stderr, "trial %d: norm^2 is %g, dense says %g\n", trial, got, n2);
                return 1;
            }
            for (uint32_t q = 0; q < NQUBITS; q++) {
                double p1 = 0.0;
                for (unsigned i = 0; i < NBASIS; i++) {
                    if (i & (1u << q)) p1 += v[i].re*v[i].re + v[i].im*v[i].im;
                }
                const double gp = limdd_prob_qubit_one(e, q, NQUBITS);
                if (fabs(gp - p1) > 1e-9 * (1.0 + n2)) {
                    fprintf(stderr, "trial %d qubit %u: P(1) mass is %g, dense says %g\n",
                            trial, q, gp, p1);
                    return 1;
                }
            }
        }
    }
    return 0;
}

/**
 * Deferring the canonical form and applying it later must land in the same
 * place as never deferring it.
 *
 * The check is edge equality, not amplitude equality: two diagrams can denote
 * the same state and still differ, and what is claimed here is the stronger
 * thing -- that limdd_canonize reconstructs exactly the diagram the eager path
 * would have built. If that holds, batching changes only when the work is
 * done, never the answer.
 */
int
test_deferred_canonization(void)
{
    const uint32_t ids[] = { GATEID_X, GATEID_Y, GATEID_Z, GATEID_H,
                             GATEID_S, GATEID_T };

    for (int trial = 0; trial < 30; trial++) {
        const uint64_t seed = rng_state;

        /* Eager: canonical after every gate. */
        limdd_set_canon_policy(LIMDD_CANON_ALWAYS, 0);
        rng_state = seed;
        LIMDD eager = limdd_all_zero_state(NQUBITS);
        for (int step = 0; step < 25; step++) {
            if ((rnd() % 3) == 0) {
                const uint32_t c = rnd() % (NQUBITS - 1);
                const uint32_t t = c + 1 + (rnd() % (NQUBITS - 1 - c));
                eager = limdd_cgate(eager, ids[rnd() % 4], UINT64_C(1) << c, t, NQUBITS);
            } else {
                eager = limdd_gate(eager, ids[rnd() % 6], rnd() % NQUBITS, NQUBITS);
            }
        }

        /* Deferred: the same gates with the orbit search off, then rebuilt. */
        limdd_set_canon_policy(LIMDD_CANON_MANUAL, 0);
        rng_state = seed;
        LIMDD lazy = limdd_all_zero_state(NQUBITS);
        for (int step = 0; step < 25; step++) {
            if ((rnd() % 3) == 0) {
                const uint32_t c = rnd() % (NQUBITS - 1);
                const uint32_t t = c + 1 + (rnd() % (NQUBITS - 1 - c));
                lazy = limdd_cgate(lazy, ids[rnd() % 4], UINT64_C(1) << c, t, NQUBITS);
            } else {
                lazy = limdd_gate(lazy, ids[rnd() % 6], rnd() % NQUBITS, NQUBITS);
            }
        }
        lazy = limdd_canonize(lazy);

        if (eager != lazy) {
            cx ve[NBASIS], vl[NBASIS];
            limdd_to_vector(eager, ve);
            limdd_to_vector(lazy, vl);
            bool same_state = true;
            for (unsigned i = 0; i < NBASIS; i++)
                if (!cx_eq(ve[i], vl[i])) same_state = false;
            fprintf(stderr, "trial %d: eager edge %llu, deferred+canonized %llu "
                            "(states %s)\n", trial,
                    (unsigned long long)eager, (unsigned long long)lazy,
                    same_state ? "agree" : "DIFFER");
            limdd_set_canon_policy(LIMDD_CANON_ALWAYS, 0);
            return 1;
        }
    }
    limdd_set_canon_policy(LIMDD_CANON_ALWAYS, 0);
    return 0;
}

TASK_0(int, runtests)
{
    limdd_nodes_init(NQUBITS, 1LL << 18, 1LL << 18, 1LL << 18, 1LL << 18);

    if (test_basis_state()) return 1;
    printf("limdd all-zero state:                    ok\n");
    if (test_single_gates()) return 1;
    printf("every 1-qubit gate matches a dense sim:  ok\n");
    if (test_controlled_gates()) return 1;
    printf("controlled gates match a dense sim:      ok\n");
    if (test_addition()) return 1;
    printf("addition matches the dense sum:          ok\n");
    if (test_reversed_and_swap()) return 1;
    printf("reversed controls and swap:              ok\n");
    if (test_random_circuits()) return 1;
    printf("40 random Clifford+T circuits, per gate: ok\n");
    if (test_deferred_canonization()) return 1;
    printf("deferred canonization equals eager:      ok\n");
    if (test_probabilities()) return 1;
    printf("norm and measurement probabilities:      ok\n");
    if (test_ops_on_skipping_diagrams()) return 1;
    printf("every op on hand-built skipping diagrams: ok\n");

    printf("(%zu nodes, %zu LIMs)\n", limdd_node_table_count(), limdd_lim_table_count());
    limdd_nodes_quit();
    return 0;
}

static int
run_with(int backend, const char *name)
{
    printf("== LIMDD operations with %s edge weights ==\n", name);
    lace_start(1, 0);
    sylvan_set_sizes(1LL << 18, 1LL << 18, 1LL << 18, 1LL << 18);
    sylvan_init_package();
    /* Through the simulator's initialiser, so the gate table is populated. */
    qsylvan_init_simulator(1LL << 18, 1LL << 18, -1, backend, NORM_LOW);

    const int res = RUN(runtests);

    sylvan_quit();
    lace_stop();
    return res;
}

int
main(void)
{
    if (run_with(COMP_HASHMAP, "complex")) return 1;
    /* The same circuits with exact coefficients: T is a pi/4 rotation, so
     * Clifford+T stays inside Q[i,sqrt2]. */
    if (run_with(QISQ2_MAP, "exact (Q[i,sqrt2])")) return 1;
    return 0;
}
