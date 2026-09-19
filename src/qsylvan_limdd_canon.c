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

#include "qsylvan_limdd_canon.h"

LIMDD_STAB
limdd_node_stab(LIMDD_TARG p)
{
    if (p == LIMDD_TERMINAL) return LIMDD_STAB_TRIVIAL;

    const uint64_t cached = limdd_node_stab_raw(p);
    if (cached != 0) return cached;

    /*
     * Recomputing rather than requiring the caller to have filled this in:
     * makeedge builds bottom-up so the children are normally already cached,
     * and a worker that arrives at a node another worker is still finishing
     * would otherwise have to wait. Both compute the same group.
     */
    const LIMDD low = limdd_node_low(p);
    const LIMDD high = limdd_node_high(p);
    const LIMDD_STAB s0 = limdd_edge_is_zero(low)  ? LIMDD_STAB_TRIVIAL
                                                   : limdd_node_stab(limdd_target(low));
    const LIMDD_STAB s1 = limdd_edge_is_zero(high) ? LIMDD_STAB_TRIVIAL
                                                   : limdd_node_stab(limdd_target(high));

    const LIMDD_STAB s = limdd_stab_of_node(limdd_node_var(p), low, high, s0, s1);
    limdd_node_set_stab_raw(p, s);
    return s;
}

LIMDD_STAB
limdd_edge_stab(LIMDD e)
{
    if (limdd_edge_is_zero(e)) return LIMDD_STAB_TRIVIAL;

    const LIMDD_STAB s = limdd_node_stab(limdd_target(e));
    const limdd_pauli_t q = limdd_lim_pauli(limdd_label(e));
    if (limdd_pauli_is_identity(q)) return s;

    /* L v is stabilised by L G L^-1, which for Pauli words is G up to the
     * sign their commutation gives. */
    LIMDD_LIM gens[LIMDD_MAX_QUBITS];
    const size_t k = limdd_stab_ngens(s);
    for (size_t i = 0; i < k; i++) {
        const LIMDD_LIM g = limdd_stab_gen(s, i);
        gens[i] = (limdd_pauli_commutation_phase(limdd_lim_pauli(g), q) == 0)
                ? g
                : limdd_lim_make(limdd_lim_pauli(g), wgt_neg(limdd_lim_weight(g)));
    }
    return limdd_stab_make(gens, k);
}

LIMDD
limdd_edge_canonical(LIMDD e)
{
    if (limdd_edge_is_zero(e)) return limdd_zero_edge();

    const LIMDD_TARG t = limdd_target(e);
    if (t == LIMDD_TERMINAL) return e;

    const LIMDD_STAB s = limdd_node_stab(t);
    if (limdd_stab_is_trivial(s)) return e;

    /*
     * L and L G denote the same state for every G in Stab(t), so the label is
     * only determined up to that right coset. Minimising over it is what makes
     * two edges for one state identical.
     *
     * limdd_stab_min_coset also offers a sign, which this coset does not have:
     * -L G is a different state, since no stabiliser group contains -I. So the
     * flip is undone. The word it chose is unaffected either way, and for a
     * fixed word the coset holds exactly one element.
     */
    bool neg;
    LIMDD_LIM m = limdd_stab_min_coset(limdd_label(e), LIMDD_STAB_TRIVIAL, s, NULL, &neg);
    if (neg) m = limdd_lim_make(limdd_lim_pauli(m), wgt_neg(limdd_lim_weight(m)));

    return limdd_bundle(m, t);
}

/**
 * `a` with the Pauli (x, z) placed at qubit `var`.
 *
 * `a` acts only above `var`, so the supports are disjoint and the product
 * picks up no phase: the bits are simply set. Anything that does need a phase,
 * such as composing X with Z at the same qubit, applies it to the scalar
 * itself -- see the swapped branch of makeedge.
 */
static LIMDD_LIM
lim_with_pauli_at(LIMDD_LIM a, uint32_t var, bool x, bool z)
{
    const uint64_t bit = UINT64_C(1) << var;
    limdd_pauli_t w = limdd_lim_pauli(a);
    assert(((w.x | w.z) & bit) == 0 && "the label already acts at this level");
    if (x) w.x |= bit;
    if (z) w.z |= bit;
    return limdd_lim_make(w, limdd_lim_weight(a));
}

/**
 * Which of the two candidate nodes is kept.
 *
 * Any fixed total order gives a canonical form; this one compares the two
 * children then the label. What matters is that the two candidates form the
 * same unordered pair however the state was presented, so that taking the
 * smaller always lands on the same node.
 */
static bool
node_less(LIMDD_TARG lo_a, LIMDD_TARG hi_a, LIMDD_LIM lab_a,
          LIMDD_TARG lo_b, LIMDD_TARG hi_b, LIMDD_LIM lab_b)
{
    if (lo_a != lo_b) return lo_a < lo_b;
    if (hi_a != hi_b) return hi_a < hi_b;
    return lab_a < lab_b;
}

/** Intern (var, I->lo, lab->hi) and make sure its group is cached. */
static LIMDD_TARG
intern(uint32_t var, LIMDD_TARG lo, LIMDD_LIM lab, LIMDD_TARG hi,
       LIMDD_STAB s_lo, LIMDD_STAB s_hi)
{
    const LIMDD low = limdd_bundle(LIMDD_LIM_IDENTITY, lo);
    const LIMDD high = limdd_bundle(lab, hi);
    int created;
    const LIMDD_TARG t = limdd_makenode_ex(var, low, high, &created);

    /*
     * A fresh bucket may be one a collection has recycled, so whatever is in
     * its cache slot belongs to whoever had it before and must be replaced,
     * not read.
     */
    if (created || limdd_node_stab_raw(t) == 0) {
        limdd_node_set_stab_raw(t, limdd_stab_of_node(var, low, high, s_lo, s_hi));
    }
    return t;
}

LIMDD
limdd_makeedge(uint32_t var, LIMDD low, LIMDD high)
{
    const bool lz = limdd_edge_is_zero(low);
    const bool hz = limdd_edge_is_zero(high);
    if (lz && hz) return limdd_zero_edge();

    /*
     * One dead branch. |0>(x)v and |1>(x)v differ by X at this level, so they
     * are one node, and the live branch is always put on the low side. This is
     * what makes |0..0> and |10..0> share every node they can.
     */
    if (lz || hz) {
        const LIMDD live = lz ? high : low;
        const LIMDD_TARG v = limdd_target(live);
        const LIMDD_TARG t = intern(var, v, LIMDD_LIM_ZERO, LIMDD_TERMINAL,
                                    limdd_node_stab(v), LIMDD_STAB_TRIVIAL);
        return limdd_edge_canonical(
            limdd_bundle(lim_with_pauli_at(limdd_label(live), var, lz, false), t));
    }

    const LIMDD_TARG v0 = limdd_target(low);
    const LIMDD_TARG v1 = limdd_target(high);
    const LIMDD_STAB s0 = limdd_node_stab(v0);
    const LIMDD_STAB s1 = limdd_node_stab(v1);

    const LIMDD_LIM a = limdd_label(low);

    /* Divide the low label out; what is left is all the node keeps. */
    const LIMDD_LIM bp = limdd_lim_mul(limdd_lim_inverse(a), limdd_label(high));

    /* Candidate as given. */
    LIMDD_LIM g1; bool neg1;
    const LIMDD_LIM lab1 = limdd_stab_min_coset(bp, s0, s1, &g1, &neg1);

    /*
     * Candidate with the branches swapped. Applying X at this level gives
     * |0>(x)B'|v1> + |1>(x)|v0>; dividing B' out in turn leaves the node
     * (v1, B'^-1, v0), so the swap costs one inverse and a second minimisation
     * over the same two groups in the other order.
     */
    LIMDD_LIM g2; bool neg2;
    const LIMDD_LIM lab2 = limdd_stab_min_coset(limdd_lim_inverse(bp), s1, s0, &g2, &neg2);

    if (node_less(v0, v1, lab1, v1, v0, lab2)) {
        const LIMDD_TARG t = intern(var, v0, lab1, v1, s0, s1);
        /*
         * Undoing the minimisation costs G on the branches above and, when the
         * sign rule flipped the label, a Z at this level -- which is exactly
         * the operator that negates the high branch and leaves the low one.
         */
        const LIMDD_LIM r = limdd_lim_mul(a, g1);
        return limdd_edge_canonical(limdd_bundle(lim_with_pauli_at(r, var, false, neg1), t));
    } else {
        const LIMDD_TARG t = intern(var, v1, lab2, v0, s1, s0);

        LIMDD_LIM r = limdd_lim_mul(limdd_lim_mul(a, bp), g2);
        if (neg2) {
            /* X here from the swap, Z from the sign rule, and XZ = -iY: the
             * factor cannot be dropped, it is part of the map. */
            r = lim_with_pauli_at(r, var, true, true);
            r = limdd_lim_make(limdd_lim_pauli(r),
                               wgt_mul(limdd_lim_weight(r), limdd_wgt_i_pow(3)));
        } else {
            r = lim_with_pauli_at(r, var, true, false);
        }
        return limdd_edge_canonical(limdd_bundle(r, t));
    }
}
