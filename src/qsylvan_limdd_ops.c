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
#include <math.h>
#include <stdlib.h>

#include <sylvan_int.h>

#include "qsylvan_limdd_ops.h"
#include "qsylvan_limdd_canon.h"
#include "qsylvan_limdd_gc.h"
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
    const uint32_t lev = limdd_level(t);
    assert(var < limdd_lims_nqubits() && "no cofactors below the last qubit");
    assert(lev >= var && "an edge may skip levels, never climb them");

    const LIMDD_LIM l = limdd_label(e);
    limdd_pauli_t p = limdd_lim_pauli(l);

    const limdd_pauli_op_t op = limdd_pauli_get(p, var);
    const bool px = (op & 2u) != 0;
    const bool pz = (op & 1u) != 0;

    /* What is left of the LIM once this qubit's Pauli is taken off it. It acts
     * on everything below and so multiplies into both branches. */
    limdd_pauli_set(&p, var, LIMDD_PAULI_I);
    const LIMDD_LIM rest = limdd_lim_make(p, limdd_lim_weight(l));

    /*
     * At a level the edge skips, both cofactors are the edge itself: the
     * level holds |0>+|1>, amplitude 1 on either branch. Feeding one edge
     * into both slots of the dispatch below then does exactly what this
     * level's Pauli does to |0>+|1> -- I and X give (V, V), Z gives (V, -V),
     * Y gives (-iV, iV) -- so the skipped case needs no code of its own.
     */
    const bool skip = lev > var;
    const LIMDD x0 = skip ? limdd_bundle(rest, t) : lim_times_edge(rest, limdd_node_low(t));
    const LIMDD x1 = skip ? x0                    : lim_times_edge(rest, limdd_node_high(t));

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

/**
 * Split `e`'s label at `level`: the part acting on levels below it comes back
 * in `*hoisted` with scalar 1, and the returned edge keeps the rest, scalar
 * included. The two parts have disjoint support, so neither the split nor
 * putting them back together picks up a phase.
 *
 * This is how a gate jumps over levels an edge skips. Levels below `level`
 * are not touched by the gate, so whatever the label does there commutes
 * with it and can wait outside the recursion. Leaving it in would hand
 * makeedge children whose labels act above their own level.
 */
static inline LIMDD
lim_split_above(LIMDD e, uint32_t level, LIMDD_LIM *hoisted)
{
    const LIMDD_LIM l = limdd_label(e);
    limdd_pauli_t p = limdd_lim_pauli(l);
    const limdd_pauli_t above = limdd_pauli_split_below(&p, level);
    if (limdd_pauli_is_identity(above)) {
        *hoisted = LIMDD_LIM_IDENTITY;
        return e;
    }
    *hoisted = limdd_lim_make(above, EVBDD_ONE);
    return limdd_bundle(limdd_lim_make(p, limdd_lim_weight(l)), limdd_target(e));
}

/* --- addition ------------------------------------------------------------ */

TASK_IMPL_3(LIMDD, limdd_plus, LIMDD, a, LIMDD, b, uint32_t, var)
{
    /* Join a collection another worker has opened, or open one if a table is
     * nearly full. Here rather than between gates: a single gate can mint
     * millions of intermediates, and this is the only point inside one where
     * stopping is safe -- everything live is on the reference stacks. */
    sylvan_gc_test();
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
     * edge cannot go there: its LIM index occupies bits 32..63, which covers
     * where the id lives, and operations would collide in each other's entries. The
     * two edges therefore go in the key slots, which are compared at full
     * width, and `var` -- small -- rides in the id.
     *
     * `var` is part of the key because the operands do not determine it:
     * (I->terminal) + (I->terminal) is 2 at level n and 2 |+>^(n-k) at level
     * k. Addition is the one operation whose level is genuinely context.
     */
    LIMDD res;
    if (cache_get3(CACHE_LIMDD_PLUS, var, na, nb, &res)) {
        return lim_times_edge(A, res);
    }

    LIMDD a0, a1, b0, b1;
    limdd_cofactors(na, var, &a0, &a1);
    limdd_cofactors(nb, var, &b0, &b1);

    /* The four cofactors are live across the recursion, and so is the branch
     * result until makeedge has consumed it. A collection may run inside
     * either call, and only these stacks tell it so. */
    limdd_refs_push(a0); limdd_refs_push(a1);
    limdd_refs_push(b0); limdd_refs_push(b1);
    limdd_refs_spawn(SPAWN(limdd_plus, a1, b1, var + 1));
    const LIMDD lo = CALL(limdd_plus, a0, b0, var + 1);
    limdd_refs_push(lo);
    const LIMDD hi = limdd_refs_sync(SYNC(limdd_plus));
    limdd_refs_push(hi);

    res = limdd_makeedge(var, lo, hi);
    limdd_refs_pop(6);
    cache_put3(CACHE_LIMDD_PLUS, var, na, nb, res);

    return lim_times_edge(A, res);
}

/* --- the low-level operations ------------------------------------------- */

TASK_IMPL_3(LIMDD, limdd_times, LIMDD, a, LIMDD, b, uint32_t, var)
{
    sylvan_gc_test();
    if (limdd_edge_is_zero(a) || limdd_edge_is_zero(b)) return limdd_zero_edge();

    const uint32_t nqubits = (uint32_t)limdd_lims_nqubits();
    if (var == nqubits) {
        assert(limdd_target(a) == LIMDD_TERMINAL && limdd_target(b) == LIMDD_TERMINAL);
        const EVBDD_WGT w = wgt_mul(limdd_lim_weight(limdd_label(a)),
                                    limdd_lim_weight(limdd_label(b)));
        if (w == EVBDD_ZERO) return limdd_zero_edge();
        return limdd_bundle(limdd_lim_make(limdd_pauli_identity(), w), LIMDD_TERMINAL);
    }

    /*
     * The product commutes, so one order. Unlike the sum, the labels cannot
     * be divided out: a Pauli's translation acts on each factor separately,
     * (P u)(Q v) is not a label times u v, so the memo is on the two edges
     * as they are, and on the level for the reason limdd_plus gives.
     */
    if (a > b) { const LIMDD t = a; a = b; b = t; }
    LIMDD res;
    if (cache_get3(CACHE_LIMDD_TIMES, var, a, b, &res)) return res;

    LIMDD a0, a1, b0, b1;
    limdd_cofactors(a, var, &a0, &a1);
    limdd_cofactors(b, var, &b0, &b1);
    limdd_refs_push(a0); limdd_refs_push(a1);
    limdd_refs_push(b0); limdd_refs_push(b1);
    limdd_refs_spawn(SPAWN(limdd_times, a1, b1, var + 1));
    const LIMDD lo = CALL(limdd_times, a0, b0, var + 1);
    limdd_refs_push(lo);
    const LIMDD hi = limdd_refs_sync(SYNC(limdd_times));
    limdd_refs_push(hi);
    res = limdd_makeedge(var, lo, hi);
    limdd_refs_pop(6);
    cache_put3(CACHE_LIMDD_TIMES, var, a, b, res);
    return res;
}

LIMDD
limdd_reduce_root(LIMDD e, uint32_t var)
{
    if (limdd_edge_is_zero(e) || !limdd_get_high_determinism()) return e;
    return limdd_edge_canonical(var, e);
}

LIMDD
limdd_negate(LIMDD e)
{
    return limdd_reduce_root(limdd_scale(e, EVBDD_MIN_ONE), 0);
}

/*
 * Restriction and projection take the cofactors at every level down to q,
 * since a label may act on x_q with a Z even where the edge skips q, and
 * rebuild above q with makeedge; at q the restriction puts the chosen
 * cofactor on both sides, and the projection puts zero on the other.
 */
static inline uint64_t
qb_key(uint32_t q, int b)
{
    return ((uint64_t)q << 1) | (uint64_t)(b & 1);
}

TASK_IMPL_4(LIMDD, limdd_restrict, LIMDD, e, uint32_t, q, int, b, uint32_t, var)
{
    sylvan_gc_test();
    if (limdd_edge_is_zero(e)) return e;
    assert(var <= q);
    LIMDD res;
    if (cache_get3(CACHE_LIMDD_RESTRICT, var, e, qb_key(q, b), &res)) return res;

    LIMDD e0, e1;
    limdd_cofactors(e, var, &e0, &e1);
    limdd_refs_push(e0); limdd_refs_push(e1);
    if (var == q) {
        /*
         * Reduced first: makeedge passes an edge with two equal children
         * through as it is, which is right for an edge out of makeedge and
         * not for a cofactor, whose label is the parent's pushed down with
         * no search for its class. Unreduced, the same state came back
         * under two labels, 11 times in 1080 (test_lowlevel_ops). Not done
         * in makeedge itself: on float weights reducing an edge that is
         * reduced already is not the identity, and it widened the LIMDD
         * of a 12-qubit Clifford+4T state past 2^t (test_limdd_structure).
         */
        const LIMDD g = limdd_reduce_root(b ? e1 : e0, var + 1);
        res = limdd_makeedge(var, g, g);
        limdd_refs_pop(2);
    } else {
        limdd_refs_spawn(SPAWN(limdd_restrict, e1, q, b, var + 1));
        const LIMDD lo = limdd_refs_push(CALL(limdd_restrict, e0, q, b, var + 1));
        const LIMDD hi = limdd_refs_push(limdd_refs_sync(SYNC(limdd_restrict)));
        res = limdd_makeedge(var, lo, hi);
        limdd_refs_pop(4);
    }
    cache_put3(CACHE_LIMDD_RESTRICT, var, e, qb_key(q, b), res);
    return res;
}

TASK_IMPL_4(LIMDD, limdd_project, LIMDD, e, uint32_t, q, int, b, uint32_t, var)
{
    sylvan_gc_test();
    if (limdd_edge_is_zero(e)) return e;
    assert(var <= q);
    LIMDD res;
    if (cache_get3(CACHE_LIMDD_PROJECT, var, e, qb_key(q, b), &res)) return res;

    LIMDD e0, e1;
    limdd_cofactors(e, var, &e0, &e1);
    limdd_refs_push(e0); limdd_refs_push(e1);
    if (var == q) {
        res = b ? limdd_makeedge(var, limdd_zero_edge(), e1)
                : limdd_makeedge(var, e0, limdd_zero_edge());
        limdd_refs_pop(2);
    } else {
        limdd_refs_spawn(SPAWN(limdd_project, e1, q, b, var + 1));
        const LIMDD lo = limdd_refs_push(CALL(limdd_project, e0, q, b, var + 1));
        const LIMDD hi = limdd_refs_push(limdd_refs_sync(SYNC(limdd_project)));
        res = limdd_makeedge(var, lo, hi);
        limdd_refs_pop(4);
    }
    cache_put3(CACHE_LIMDD_PROJECT, var, e, qb_key(q, b), res);
    return res;
}

LIMDD
limdd_local_matvec(LIMDD e, const EVBDD_WGT *M, const uint32_t *qubits, uint32_t k,
                   uint32_t nqubits)
{
    (void)nqubits;
    const uint32_t dim = 1u << k;
    LIMDD *vc = malloc(dim * sizeof(LIMDD));
    if (vc == NULL) { fprintf(stderr, "sylvan: out of memory in limdd_local_matvec\n"); exit(1); }

    /* every intermediate on the reference stack: each operation may join a
     * collection */
    for (uint32_t c = 0; c < dim; c++) {
        LIMDD r = e;
        for (uint32_t i = 0; i < k; i++)
            r = limdd_restrict(r, qubits[i], (int)((c >> (k - 1 - i)) & 1), 0);
        vc[c] = limdd_refs_push(r);
    }
    LIMDD out = limdd_refs_push(limdd_zero_edge());
    for (uint32_t r = 0; r < dim; r++) {
        LIMDD w = limdd_zero_edge();
        for (uint32_t c = 0; c < dim; c++) {
            const EVBDD_WGT m = M[(size_t)r * dim + c];
            if (m == EVBDD_ZERO) continue;
            limdd_refs_push(w);
            w = limdd_plus(w, limdd_scale(vc[c], m), 0);
            limdd_refs_pop(1);
        }
        for (uint32_t i = 0; i < k; i++) {
            limdd_refs_push(w);
            w = limdd_project(w, qubits[i], (int)((r >> (k - 1 - i)) & 1), 0);
            limdd_refs_pop(1);
        }
        limdd_refs_push(w);
        const LIMDD next = limdd_plus(out, w, 0);
        limdd_refs_pop(2);
        out = limdd_refs_push(next);
    }
    limdd_refs_pop(1 + (long)dim);
    free(vc);
    return out;
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
    sylvan_gc_test();
    if (limdd_edge_is_zero(e)) return e;

    /*
     * The level to work at is the first one that matters: the target's own
     * node, or the gate's qubit if the edge skips it. Levels above that are
     * skipped and untouched by the gate, and the label's entries there are
     * hoisted out of the recursion (see below). This is also what keeps the
     * cache key sound without a level in it: `var` is a function of the edge
     * and the target, so equal keys mean equal work whatever level the caller
     * was at.
     */
    const uint32_t lev = limdd_level(limdd_target(e));
    const uint32_t var = lev < target ? lev : target;
    assert(var < nqubits);

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
    assert(gateid < (1u << 20) && target < LIMDD_MAX_QUBITS);
    const uint64_t opid = CACHE_LIMDD_GATE | ((uint64_t)gateid << 20) | target;

    LIMDD res;
    if (cache_get3(opid, 0, e, 0, &res)) return res;

    LIMDD_LIM hoisted;
    const LIMDD inner = lim_split_above(e, var, &hoisted);
    if (inner != e) {
        res = lim_times_edge(hoisted, CALL(limdd_gate, inner, gateid, target, nqubits));
        cache_put3(opid, 0, e, 0, res);
        return res;
    }

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

        limdd_refs_push(a); limdd_refs_push(b);
        limdd_refs_push(c); limdd_refs_push(d);
        limdd_refs_spawn(SPAWN(limdd_plus, c, d, var + 1));
        const LIMDD new_lo = CALL(limdd_plus, a, b, var + 1);
        limdd_refs_push(new_lo);
        const LIMDD new_hi = limdd_refs_sync(SYNC(limdd_plus));
        limdd_refs_push(new_hi);
        res = limdd_makeedge(var, new_lo, new_hi);
        limdd_refs_pop(6);
    } else {
        limdd_refs_push(lo); limdd_refs_push(hi);
        limdd_refs_spawn(SPAWN(limdd_gate, hi, gateid, target, nqubits));
        const LIMDD new_lo = CALL(limdd_gate, lo, gateid, target, nqubits);
        limdd_refs_push(new_lo);
        const LIMDD new_hi = limdd_refs_sync(SYNC(limdd_gate));
        limdd_refs_push(new_hi);
        res = limdd_makeedge(var, new_lo, new_hi);
        limdd_refs_pop(4);
    }

    cache_put3(opid, 0, e, 0, res);
    return res;
}

TASK_IMPL_5(LIMDD, limdd_cgate, LIMDD, e, uint32_t, gateid, uint64_t, controls,
            uint32_t, target, uint32_t, nqubits)
{
    sylvan_gc_test();
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

    /*
     * Work at the first level that matters: the target's node, or the
     * highest control if the edge skips it. A skipped control is a control
     * in superposition, and cofactoring it is exactly what entangles, so the
     * recursion has to stop there even though there is no node. As in
     * limdd_gate, `var` is a function of the key -- edge and controls -- so
     * the key needs no level of its own.
     */
    const uint32_t lev = limdd_level(limdd_target(e));
    const uint32_t top = (uint32_t)__builtin_ctzll(controls);
    const uint32_t var = lev < top ? lev : top;
    assert(var < nqubits);

    /* The label stays on, for the same reason as in limdd_gate. */
    assert(gateid < (1u << 20) && target < LIMDD_MAX_QUBITS);
    const uint64_t opid = CACHE_LIMDD_CGATE | ((uint64_t)gateid << 20) | target;

    LIMDD res;
    if (cache_get3(opid, 0, e, controls, &res)) return res;

    LIMDD_LIM hoisted;
    const LIMDD inner = lim_split_above(e, var, &hoisted);
    if (inner != e) {
        res = lim_times_edge(hoisted, CALL(limdd_cgate, inner, gateid, controls, target, nqubits));
        cache_put3(opid, 0, e, controls, res);
        return res;
    }

    LIMDD lo, hi;
    limdd_cofactors(e, var, &lo, &hi);

    const uint64_t bit = UINT64_C(1) << var;
    if (controls & bit) {
        /* A control qubit: the |0> branch is untouched, and only the |1>
         * branch continues carrying the remaining controls. */
        limdd_refs_push(lo); limdd_refs_push(hi);
        const LIMDD new_hi = CALL(limdd_cgate, hi, gateid, controls & ~bit,
                                  target, nqubits);
        limdd_refs_push(new_hi);
        res = limdd_makeedge(var, lo, new_hi);
        limdd_refs_pop(3);
    } else {
        limdd_refs_push(lo); limdd_refs_push(hi);
        limdd_refs_spawn(SPAWN(limdd_cgate, hi, gateid, controls, target, nqubits));
        const LIMDD new_lo = CALL(limdd_cgate, lo, gateid, controls, target,
                                  nqubits);
        limdd_refs_push(new_lo);
        const LIMDD new_hi = limdd_refs_sync(SYNC(limdd_cgate));
        limdd_refs_push(new_hi);
        res = limdd_makeedge(var, new_lo, new_hi);
        limdd_refs_pop(4);
    }

    cache_put3(opid, 0, e, controls, res);
    return res;
}

/* --- controls on either side, and swap -------------------------------- */

/** CNOT with the control BELOW the target, via Hadamard conjugation. */
static LIMDD
cx_reversed(LIMDD e, uint32_t c, uint32_t t, uint32_t nq)
{
    assert(c > t);
    e = limdd_gate(e, GATEID_H, c, nq);
    e = limdd_gate(e, GATEID_H, t, nq);
    e = limdd_cgate(e, GATEID_X, UINT64_C(1) << t, c, nq);  /* now t < c */
    e = limdd_gate(e, GATEID_H, t, nq);
    e = limdd_gate(e, GATEID_H, c, nq);
    return e;
}

LIMDD
limdd_cgate_either(LIMDD e, uint32_t gateid, uint32_t control, uint32_t target,
                   uint32_t nqubits, bool *ok)
{
    *ok = true;
    if (control < target) return limdd_cgate(e, gateid, UINT64_C(1) << control,
                                             target, nqubits);

    if (gateid == GATEID_Z) {
        /* CZ = diag(1,1,1,-1) is symmetric in its two qubits. */
        return limdd_cgate(e, GATEID_Z, UINT64_C(1) << target, control, nqubits);
    }
    if (gateid == GATEID_X) return cx_reversed(e, control, target, nqubits);

    *ok = false;
    return e;
}

LIMDD
limdd_swap(LIMDD e, uint32_t a, uint32_t b, uint32_t nqubits)
{
    if (a == b) return e;
    if (a > b) { const uint32_t t = a; a = b; b = t; }

    e = limdd_cgate(e, GATEID_X, UINT64_C(1) << a, b, nqubits);
    e = cx_reversed(e, b, a, nqubits);
    e = limdd_cgate(e, GATEID_X, UINT64_C(1) << a, b, nqubits);
    return e;
}

/* --- states and probabilities -------------------------------------------- */

/* --- counting the live nodes --------------------------------------------- */

typedef struct {
    LIMDD_TARG *slot;
    size_t      cap;      /* power of two */
    size_t      used;
} seen_t;

static bool
seen_add(seen_t *s, LIMDD_TARG t)
{
    size_t i = (size_t)((t * UINT64_C(0x9E3779B97F4A7C15)) & (s->cap - 1));
    for (;;) {
        if (s->slot[i] == 0) { s->slot[i] = t; s->used++; return true; }
        if (s->slot[i] == t) return false;
        i = (i + 1) & (s->cap - 1);
    }
}

static void
count_rec(LIMDD e, seen_t *s)
{
    if (limdd_edge_is_zero(e)) return;
    const LIMDD_TARG t = limdd_target(e);
    if (t == LIMDD_TERMINAL) return;
    if (!seen_add(s, t)) return;
    count_rec(limdd_node_low(t), s);
    count_rec(limdd_node_high(t), s);
}

uint64_t
limdd_countnodes(LIMDD e)
{
    /* Sized from the table so the probe sequence never fills; the walk visits
     * at most that many distinct nodes. */
    seen_t s;
    s.cap = 1;
    while (s.cap < 4 * (limdd_node_table_count() + 8)) s.cap <<= 1;
    s.slot = calloc(s.cap, sizeof(LIMDD_TARG));
    if (s.slot == NULL) return 0;
    s.used = 0;

    count_rec(e, &s);
    const uint64_t n = s.used;
    free(s.slot);
    return n;
}

/* --- what an inline LIM encoding would have to hold ----------------------- */

typedef struct {
    /* Indexed by qubit, so they must be sized by the qubit limit and not by
     * a literal 64. Kept as members rather than pointers so that the
     * memset(&st, 0, sizeof st) below still clears the counters. */
    uint64_t nodes[LIMDD_MAX_QUBITS];
    uint64_t support_sum[LIMDD_MAX_QUBITS];
    uint64_t support_max[LIMDD_MAX_QUBITS];
    uint64_t by_support[LIMDD_MAX_QUBITS + 1];
    uint64_t stab_uncomputed, stab_trivial;
    /* A group on n qubits has at most n generators, so n+1 buckets and a
     * slot for the overflow clamp. (2n is STAB_SCRATCH_ROWS, a bound on rows
     * handed to the reducer, which is a different quantity.) */
    uint64_t by_ngens[LIMDD_MAX_QUBITS + 2];
    uint64_t *handles; size_t nhandles, caphandles;
    LIMDD_LIM *lims;           /* distinct high LIMs seen */
    size_t    nlims, caplims;
    EVBDD_WGT *wgts;           /* distinct weights among live labels */
    size_t    nwgts, capwgts;
} lim_stats_t;

static void
lim_stats_rec(LIMDD e, seen_t *s, lim_stats_t *st, uint32_t nqubits)
{
    if (limdd_edge_is_zero(e)) return;
    const LIMDD_TARG t = limdd_target(e);
    if (t == LIMDD_TERMINAL) return;
    if (!seen_add(s, t)) return;

    const uint32_t var = limdd_node_var(t);
    const LIMDD high = limdd_node_high(t);
    st->nodes[var]++;

    {   /* the cached stabiliser group, read raw so nothing is computed here */
        const uint64_t g = limdd_node_stab_raw(t);
        if (g == 0) st->stab_uncomputed++;
        else if (limdd_stab_is_trivial(g)) st->stab_trivial++;
        else {
            const size_t k = limdd_stab_ngens(g);
            st->by_ngens[k <= LIMDD_MAX_QUBITS ? k : LIMDD_MAX_QUBITS + 1]++;
            bool fresh = true;
            for (size_t i = 0; i < st->nhandles; i++) if (st->handles[i] == g) { fresh = false; break; }
            if (fresh) {
                if (st->nhandles == st->caphandles) {
                    st->caphandles = st->caphandles ? st->caphandles * 2 : 1024;
                    st->handles = realloc(st->handles, st->caphandles * sizeof(uint64_t));
                }
                st->handles[st->nhandles++] = g;
            }
        }
    }

    if (!limdd_edge_is_zero(high)) {
        const LIMDD_LIM l = limdd_label(high);
        const limdd_pauli_t p = limdd_lim_pauli(l);
        const uint64_t sup = (uint64_t)limdd_pauli_weight(p);
        st->support_sum[var] += sup;
        if (sup > st->support_max[var]) st->support_max[var] = sup;

        bool fresh = true;
        for (size_t i = 0; i < st->nlims; i++) if (st->lims[i] == l) { fresh = false; break; }
        if (fresh) {
            if (st->nlims == st->caplims) {
                st->caplims = st->caplims ? st->caplims * 2 : 1024;
                st->lims = realloc(st->lims, st->caplims * sizeof(LIMDD_LIM));
            }
            st->lims[st->nlims++] = l;
            const EVBDD_WGT wv = limdd_lim_weight(l);
            bool wfresh = true;
            for (size_t i = 0; i < st->nwgts; i++)
                if (st->wgts[i] == wv) { wfresh = false; break; }
            if (wfresh) {
                if (st->nwgts == st->capwgts) {
                    st->capwgts = st->capwgts ? st->capwgts * 2 : 1024;
                    st->wgts = realloc(st->wgts, st->capwgts * sizeof(EVBDD_WGT));
                }
                st->wgts[st->nwgts++] = wv;
            }
            st->by_support[sup <= LIMDD_MAX_QUBITS ? sup : LIMDD_MAX_QUBITS]++;
        }
    }

    lim_stats_rec(limdd_node_low(t), s, st, nqubits);
    lim_stats_rec(high, s, st, nqubits);
}

void
limdd_report_lim_stats(FILE *out, LIMDD e, uint32_t nqubits)
{
    seen_t s;
    s.cap = 1;
    while (s.cap < 4 * (limdd_node_table_count() + 8)) s.cap <<= 1;
    s.slot = calloc(s.cap, sizeof(LIMDD_TARG));
    if (s.slot == NULL) return;
    s.used = 0;

    lim_stats_t st;
    memset(&st, 0, sizeof st);
    lim_stats_rec(e, &s, &st, nqubits);

    fprintf(out, "# level nodes mean_support max_support dense_bits_needed\n");
    uint64_t tot = 0, fit23 = 0;
    for (uint32_t v = 0; v < nqubits; v++) {
        if (st.nodes[v] == 0) continue;
        /* A high edge is read at v+1, so its Pauli is confined to v+1..n-1. */
        const uint32_t dense = 2u * (nqubits - v - 1);
        fprintf(out, "LIMSTAT %u %llu %.2f %llu %u\n", v,
                (unsigned long long)st.nodes[v],
                (double)st.support_sum[v] / (double)st.nodes[v],
                (unsigned long long)st.support_max[v], dense);
        tot += st.nodes[v];
        if (dense <= 23) fit23 += st.nodes[v];
    }
    fprintf(out, "LIMSTAT-TABLES pauli=%zu lim=%zu nodes=%zu stabcells=%zu\n",
            limdd_pauli_table_count(), limdd_lim_table_count(),
            limdd_node_table_count(), limdd_stab_table_count());
    fprintf(out, "LIMSTAT-TOTAL nodes=%llu fit_dense_23bit=%llu (%.1f%%) distinct_high_lims=%zu\n",
            (unsigned long long)tot, (unsigned long long)fit23,
            tot ? 100.0 * (double)fit23 / (double)tot : 0.0, st.nlims);
    fprintf(out, "LIMSTAT-WGT live_distinct=%zu interned_estimate=%llu\n",
            st.nwgts, (unsigned long long)wgt_table_entries_estimate());
    fprintf(out, "LIMSTAT-STAB uncomputed=%llu trivial=%llu nontrivial_distinct=%zu\n",
            (unsigned long long)st.stab_uncomputed, (unsigned long long)st.stab_trivial,
            st.nhandles);
    for (unsigned k = 0; k <= LIMDD_MAX_QUBITS + 1; k++) {
        if (st.by_ngens[k]) fprintf(out, "LIMGENS %u %llu\n", k, (unsigned long long)st.by_ngens[k]);
    }
    fprintf(out, "# distinct high LIMs by Pauli support\n");
    for (unsigned k = 0; k <= LIMDD_MAX_QUBITS; k++) {
        if (st.by_support[k]) fprintf(out, "LIMSUP %u %llu\n", k, (unsigned long long)st.by_support[k]);
    }

    free(st.wgts);
    free(st.handles);
    free(st.lims);
    free(s.slot);
}

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

    /*
     * Every level the edge skips holds |0>+|1>, of squared norm 2, whatever
     * Pauli the label puts there -- all four send it to a vector of squared
     * norm 2. So the skipped levels contribute a power of two and nothing
     * else, and the walk resumes at the target's own level.
     */
    const LIMDD_TARG t = limdd_target(e);
    const uint32_t lev = limdd_level(t);
    assert(lev >= var && lev <= nqubits);
    const double skipped = ldexp(1.0, (int)(lev - var));

    if (t == LIMDD_TERMINAL) {
        const complex_t w = weight_as_complex(limdd_lim_weight(limdd_label(e)));
        return skipped * (w.r * w.r + w.i * w.i);
    }

    /*
     * Memoised on the edge and the level. Without this the walk is a tree
     * walk over a DAG, so a shared subdiagram is re-summed once per path
     * reaching it -- 25% of the run on a 20-qubit Clifford+T circuit, and
     * that is only the norm the simulator reports at the end, not the
     * simulation. The level is part of the key for the same reason it is in
     * limdd_plus's: an edge does not say what level it is read at, and the
     * skipped levels above it each double the result.
     */
    union { double d; uint64_t u; } conv;
    if (cache_get3(CACHE_LIMDD_NORMSQ, var, e, 0, &conv.u)) return conv.d;

    LIMDD lo, hi;
    limdd_cofactors(e, lev, &lo, &hi);
    const double res = skipped * (norm_sq(lo, lev + 1, nqubits)
                                + norm_sq(hi, lev + 1, nqubits));
    conv.d = res;
    cache_put3(CACHE_LIMDD_NORMSQ, var, e, 0, conv.u);
    return res;
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

    const LIMDD_TARG t = limdd_target(e);
    const uint32_t lev = limdd_level(t);
    assert(lev >= var);

    /*
     * A skipped qubit has exactly half its mass in |1>: each of I, X, Z and
     * Y sends (1, 1) to a pair of equal magnitude. This covers the terminal
     * too, whose level is nqubits > qubit; before levels could be skipped an
     * edge reaching the terminal had passed the qubit and contributed 0
     * here, which would now silently report a qubit in superposition as
     * certainly |0>.
     */
    if (qubit < lev) return 0.5 * norm_sq(e, var, nqubits);

    const double skipped = ldexp(1.0, (int)(lev - var));
    LIMDD lo, hi;
    limdd_cofactors(e, lev, &lo, &hi);

    if (lev == qubit) return skipped * norm_sq(hi, lev + 1, nqubits);
    return skipped * (prob_one(lo, lev + 1, qubit, nqubits)
                    + prob_one(hi, lev + 1, qubit, nqubits));
}

double
limdd_prob_qubit_one(LIMDD e, uint32_t qubit, uint32_t nqubits)
{
    return prob_one(e, 0, qubit, nqubits);
}
