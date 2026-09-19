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

#include <assert.h>
#include <stdlib.h>

#include <sylvan_int.h>

#include "qsylvan_limdd_ops.h"
#include "qsylvan_gates.h"

LIMDD
limdd_scale(LIMDD e, EVBDD_WGT w)
{
    if (w == EVBDD_ZERO || limdd_edge_is_zero(e)) return limdd_zero_edge();
    if (w == EVBDD_ONE) return e;

    const LIMDD_LIM l = limdd_label(e);
    return limdd_bundle(limdd_lim_make(limdd_lim_pauli(l),
                                       wgt_mul(limdd_lim_weight(l), w)),
                        limdd_target(e));
}

/** `l` times the edge `e`, i.e. the same target under a composed map. */
static inline LIMDD
lim_times_edge(LIMDD_LIM l, LIMDD e)
{
    if (limdd_lim_is_identity(l)) return e;
    if (limdd_edge_is_zero(e) || limdd_lim_is_zero(l)) return limdd_zero_edge();
    return limdd_bundle(limdd_lim_mul(l, limdd_label(e)), limdd_target(e));
}

void
limdd_cofactors(LIMDD e, uint32_t var, LIMDD *low, LIMDD *high)
{
    if (limdd_edge_is_zero(e)) {
        *low = *high = limdd_zero_edge();
        return;
    }

    const LIMDD_TARG t = limdd_target(e);
    assert(t != LIMDD_TERMINAL && "no cofactors below the last qubit");
    assert(limdd_node_var(t) == var);

    const LIMDD_LIM l = limdd_label(e);
    limdd_pauli_t p = limdd_lim_pauli(l);

    const uint64_t bit = UINT64_C(1) << var;
    const bool px = (p.x & bit) != 0;
    const bool pz = (p.z & bit) != 0;

    /* What is left of the LIM once this qubit's Pauli is taken off it. It acts
     * on everything below and so multiplies into both branches. */
    p.x &= ~bit;
    p.z &= ~bit;
    const LIMDD_LIM rest = limdd_lim_make(p, limdd_lim_weight(l));

    const LIMDD x0 = lim_times_edge(rest, limdd_node_low(t));
    const LIMDD x1 = lim_times_edge(rest, limdd_node_high(t));

    if (!px && !pz) {            /* I */
        *low = x0;
        *high = x1;
    } else if (!px && pz) {      /* Z: |1> picks up a minus */
        *low = x0;
        *high = limdd_scale(x1, EVBDD_MIN_ONE);
    } else if (px && !pz) {      /* X: the branches swap */
        *low = x1;
        *high = x0;
    } else {                     /* Y = iXZ: swap, and the two halves differ
                                  * by more than a sign -- -i below, +i above */
        *low = limdd_scale(x1, limdd_wgt_i_pow(3));
        *high = limdd_scale(x0, limdd_wgt_i_pow(1));
    }
}

/* --- addition ------------------------------------------------------------ */

TASK_IMPL_3(LIMDD, limdd_plus, LIMDD, a, LIMDD, b, uint32_t, var)
{
    if (limdd_edge_is_zero(a)) return b;
    if (limdd_edge_is_zero(b)) return a;

    const uint32_t nqubits = (uint32_t)limdd_lims_nqubits();
    if (var == nqubits) {
        /* Both are scalars on the terminal; no Pauli is left to act. */
        assert(limdd_target(a) == LIMDD_TERMINAL);
        assert(limdd_target(b) == LIMDD_TERMINAL);
        const EVBDD_WGT w = wgt_add(limdd_lim_weight(limdd_label(a)),
                                    limdd_lim_weight(limdd_label(b)));
        if (w == EVBDD_ZERO) return limdd_zero_edge();
        return limdd_bundle(limdd_lim_make(limdd_pauli_identity(), w),
                            LIMDD_TERMINAL);
    }

    /*
     * Addition commutes, so fix an order; then divide the first label out, so
     * that A|u> + B|v> and (cA)|u> + (cB)|v> share a cache entry instead of
     * occupying two.
     */
    if (a > b) { const LIMDD t = a; a = b; b = t; }

    const LIMDD_LIM A = limdd_label(a);
    const LIMDD_LIM Ainv = limdd_lim_inverse(A);
    const LIMDD na = limdd_bundle(LIMDD_LIM_IDENTITY, limdd_target(a));
    const LIMDD nb = lim_times_edge(Ainv, b);

    /*
     * cache_get3 ORs its second argument into the operation id, so a LIMDD
     * edge cannot go there: its LIM index occupies bits 40..62, exactly where
     * the id lives, and operations would collide in each other's entries.
     * Everything therefore goes in the two key slots, which are compared at
     * full width. `var` is left out because it is implied -- it is the
     * variable of na's node.
     */
    LIMDD res;
    if (cache_get3(CACHE_LIMDD_PLUS, 0, na, nb, &res)) {
        return lim_times_edge(A, res);
    }

    LIMDD a0, a1, b0, b1;
    limdd_cofactors(na, var, &a0, &a1);
    limdd_cofactors(nb, var, &b0, &b1);

    SPAWN(limdd_plus, a1, b1, var + 1);
    const LIMDD lo = CALL(limdd_plus, a0, b0, var + 1);
    const LIMDD hi = SYNC(limdd_plus);

    res = limdd_makeedge(var, lo, hi);
    cache_put3(CACHE_LIMDD_PLUS, 0, na, nb, res);

    return lim_times_edge(A, res);
}

/* --- applying a gate ----------------------------------------------------- */

/**
 * The gate's four entries, as edge weights. Q-Sylvan keeps them in the same
 * table as amplitudes, which is what lets an exact backend stay exact: a gate
 * built from qisq2 literals composes with qisq2 amplitudes without rounding.
 */
static inline void
gate_entries(uint32_t gateid, EVBDD_WGT u[4])
{
    u[0] = gates[gateid][0];
    u[1] = gates[gateid][1];
    u[2] = gates[gateid][2];
    u[3] = gates[gateid][3];
}

TASK_IMPL_4(LIMDD, limdd_gate, LIMDD, e, uint32_t, gateid, uint32_t, target,
            uint32_t, nqubits)
{
    if (limdd_edge_is_zero(e)) return e;

    LIMDD_TARG t = limdd_target(e);
    const uint32_t var = (t == LIMDD_TERMINAL) ? nqubits : limdd_node_var(t);
    assert(var <= target && "gate applied below its target");

    /*
     * The label is NOT divided out here, unlike in limdd_plus. Addition is
     * linear, so A|u> + A|v> = A(|u> + |v>) and the label can be taken outside
     * and put back. A gate cannot: U(A|v>) is not A(U|v>), because a gate and
     * a Pauli on the same qubit do not commute -- H then X is not X then H.
     * Pulling the label out and reapplying it afterwards silently applies the
     * gate to the wrong state.
     *
     * Nothing is lost by leaving it in: limdd_cofactors already pushes the
     * label through the node as it descends, so the recursion sees the right
     * state at every level. The cost is a colder cache, since it is keyed on
     * the whole edge rather than just the node.
     *
     * The gate and its target are small, so they ride in the operation id --
     * inside this op's own 2^40 block, which keeps it clear of its
     * neighbours. See the note in limdd_plus about why they cannot go in the
     * second argument.
     */
    assert(gateid < (1u << 20) && target < 64);
    const uint64_t opid = CACHE_LIMDD_GATE | ((uint64_t)gateid << 20) | target;

    LIMDD res;
    if (cache_get3(opid, 0, e, 0, &res)) return res;

    LIMDD lo, hi;
    limdd_cofactors(e, var, &lo, &hi);

    if (var == target) {
        /*
         * At the target the two branches mix:
         *   new |0> component = u00*lo + u01*hi
         *   new |1> component = u10*lo + u11*hi
         */
        EVBDD_WGT u[4];
        gate_entries(gateid, u);

        const LIMDD a = limdd_scale(lo, u[0]);
        const LIMDD b = limdd_scale(hi, u[1]);
        const LIMDD c = limdd_scale(lo, u[2]);
        const LIMDD d = limdd_scale(hi, u[3]);

        SPAWN(limdd_plus, c, d, var + 1);
        const LIMDD new_lo = CALL(limdd_plus, a, b, var + 1);
        const LIMDD new_hi = SYNC(limdd_plus);
        res = limdd_makeedge(var, new_lo, new_hi);
    } else {
        SPAWN(limdd_gate, hi, gateid, target, nqubits);
        const LIMDD new_lo = CALL(limdd_gate, lo, gateid, target, nqubits);
        const LIMDD new_hi = SYNC(limdd_gate);
        res = limdd_makeedge(var, new_lo, new_hi);
    }

    cache_put3(opid, 0, e, 0, res);
    return res;
}

TASK_IMPL_5(LIMDD, limdd_cgate, LIMDD, e, uint32_t, gateid, uint64_t, controls,
            uint32_t, target, uint32_t, nqubits)
{
    if (controls == 0) return CALL(limdd_gate, e, gateid, target, nqubits);
    if (limdd_edge_is_zero(e)) return e;

    /*
     * Every control must sit above the target in the variable order. The
     * recursion fixes a control by descending into its |1> branch, so a
     * control BELOW the target is still unresolved when the target's two
     * branches have to be mixed -- and the mixing coefficient would then
     * depend on a qubit not yet examined. Q-Sylvan's QMDD has the same
     * restriction and falls back to building the gate as a matrix and doing a
     * matrix-vector product; LIMDD has no matrix representation yet.
     */
    assert(controls < (UINT64_C(1) << target)
           && "limdd_cgate needs every control above the target");

    LIMDD_TARG t = limdd_target(e);
    const uint32_t var = (t == LIMDD_TERMINAL) ? nqubits : limdd_node_var(t);

    /* The label stays on, for the same reason as in limdd_gate. */
    assert(gateid < (1u << 20) && target < 64);
    const uint64_t opid = CACHE_LIMDD_CGATE | ((uint64_t)gateid << 20) | target;

    LIMDD res;
    if (cache_get3(opid, 0, e, controls, &res)) return res;

    LIMDD lo, hi;
    limdd_cofactors(e, var, &lo, &hi);

    const uint64_t bit = UINT64_C(1) << var;
    if (controls & bit) {
        /* A control qubit: the |0> branch is untouched, and only the |1>
         * branch continues carrying the remaining controls. */
        const LIMDD new_hi = CALL(limdd_cgate, hi, gateid, controls & ~bit,
                                  target, nqubits);
        res = limdd_makeedge(var, lo, new_hi);
    } else {
        SPAWN(limdd_cgate, hi, gateid, controls, target, nqubits);
        const LIMDD new_lo = CALL(limdd_cgate, lo, gateid, controls, target,
                                  nqubits);
        const LIMDD new_hi = SYNC(limdd_cgate);
        res = limdd_makeedge(var, new_lo, new_hi);
    }

    cache_put3(opid, 0, e, controls, res);
    return res;
}

/* --- states and probabilities -------------------------------------------- */

LIMDD
limdd_all_zero_state(uint32_t nqubits)
{
    LIMDD e = limdd_one_edge();
    for (int q = (int)nqubits - 1; q >= 0; q--) {
        e = limdd_makeedge((uint32_t)q, e, limdd_zero_edge());
    }
    return e;
}

/**
 * Sum of |amplitude|^2 over what `e` denotes.
 *
 * The edge's LIM is pushed down by limdd_cofactors rather than being peeled
 * off and its scalar applied separately. The total norm would survive that
 * shortcut -- a Pauli permutes basis states and only the scalar changes the
 * magnitude -- but the per-qubit version below would not: X exchanges a
 * qubit's |0> and |1> branches, which is exactly what it is measuring.
 * Treating both the same way keeps them from drifting apart.
 */
static double
norm_sq(LIMDD e, uint32_t var, uint32_t nqubits)
{
    if (limdd_edge_is_zero(e)) return 0.0;

    if (limdd_target(e) == LIMDD_TERMINAL) {
        assert(var == nqubits);
        const complex_t w = weight_as_complex(limdd_lim_weight(limdd_label(e)));
        return w.r * w.r + w.i * w.i;
    }

    LIMDD lo, hi;
    limdd_cofactors(e, var, &lo, &hi);
    return norm_sq(lo, var + 1, nqubits) + norm_sq(hi, var + 1, nqubits);
}

double
limdd_norm_squared(LIMDD e, uint32_t nqubits)
{
    return norm_sq(e, 0, nqubits);
}

static double
prob_one(LIMDD e, uint32_t var, uint32_t qubit, uint32_t nqubits)
{
    if (limdd_edge_is_zero(e)) return 0.0;
    if (limdd_target(e) == LIMDD_TERMINAL) return 0.0;

    LIMDD lo, hi;
    limdd_cofactors(e, var, &lo, &hi);

    if (var == qubit) return norm_sq(hi, var + 1, nqubits);
    return prob_one(lo, var + 1, qubit, nqubits)
         + prob_one(hi, var + 1, qubit, nqubits);
}

double
limdd_prob_qubit_one(LIMDD e, uint32_t qubit, uint32_t nqubits)
{
    return prob_one(e, 0, qubit, nqubits);
}
