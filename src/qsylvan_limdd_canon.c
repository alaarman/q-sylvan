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

/* Declared in qsylvan_limdd_ops.h, which includes this module's header, so it
 * cannot be included from here. */
uint64_t limdd_countnodes(LIMDD e);

/** `l` composed onto `e`'s label: the same target under a composed map. */
static inline LIMDD
lim_times_edge(LIMDD_LIM l, LIMDD e)
{
    if (limdd_lim_is_identity(l)) return e;
    if (limdd_edge_is_zero(e) || limdd_lim_is_zero(l)) return limdd_zero_edge();
    return limdd_bundle(limdd_lim_mul(l, limdd_label(e)), limdd_target(e));
}

/* Read-mostly: set once at configuration time, then only read. */
static bool high_determinism = true;

static limdd_canon_policy_t canon_policy = LIMDD_CANON_ALWAYS;
static uint64_t canon_interval = 0;      /* work units between rebuilds */
static uint64_t canon_before = 0, canon_after = 0;

/*
 * The work counter, one cache line per worker.
 *
 * A single counter would be simpler and wrong: the natural thing to count is
 * work done -- nodes made, operations applied -- which happens inside the
 * parallel recursions, so every worker would read-modify-write one line
 * millions of times and spend more on the counting than on the canonical form
 * it is scheduling. Each worker owning a line means no invalidation at all,
 * and the sum is taken only when the driver asks, once per gate.
 */
typedef struct {
    uint64_t n;
    char pad[SYLVAN_SHARING_PAD - sizeof(uint64_t)];
} canon_counter_t;

static canon_counter_t *canon_counters = NULL;
static unsigned canon_nworkers = 0;
/* A rebuild calls makeedge for every node it visits; counting those would
 * schedule the next rebuild immediately. */
static bool canon_in_progress = false;

static void
canon_counters_reset(void)
{
    for (unsigned i = 0; i < canon_nworkers; i++) canon_counters[i].n = 0;
}

void
limdd_canon_count(void)
{
    if (canon_counters == NULL || canon_in_progress) return;
    WorkerP *w = lace_get_worker();
    /* The thread that called lace_start is not a worker; its work is counted
     * in slot 0, which no worker uses for anything else while it runs. */
    canon_counters[w ? (unsigned)w->worker : 0].n++;
}

void
limdd_set_canon_policy(limdd_canon_policy_t policy, uint64_t interval)
{
    canon_policy = policy;
    canon_interval = interval ? interval : 1;

    if (canon_counters == NULL) {
        canon_nworkers = (unsigned)lace_workers();
        if (canon_nworkers == 0) canon_nworkers = 1;
        canon_counters = (canon_counter_t *)
            sylvan_alloc_padded(canon_nworkers * sizeof(canon_counter_t));
        if (canon_counters == NULL) {
            fprintf(stderr, "sylvan: out of memory for the LIMDD work counters\n");
            exit(1);
        }
    }
    canon_counters_reset();
    /* Deferring means operations run without the orbit search. */
    high_determinism = (policy == LIMDD_CANON_ALWAYS);
}

void
limdd_canon_last(uint64_t *before, uint64_t *after)
{
    if (before) *before = canon_before;
    if (after)  *after  = canon_after;
}

bool
limdd_canon_due(void)
{
    if (canon_policy == LIMDD_CANON_ALWAYS ||
        canon_policy == LIMDD_CANON_MANUAL) return false;

    uint64_t total = 0;
    for (unsigned i = 0; i < canon_nworkers; i++) total += canon_counters[i].n;
    return total >= canon_interval;
}

void
limdd_set_high_determinism(bool on)
{
    high_determinism = on;
}

bool
limdd_get_high_determinism(void)
{
    return high_determinism;
}

/**
 * The group of what `t` denotes when reached by an edge read at `level`: the
 * node's own group, extended by X_j for each level in [level, limdd_level(t))
 * that such an edge skips. This is the group BEFORE any label is applied,
 * which is what makeedge's coset and limdd_stab_of_node consume; the group of
 * an edge's actual state is limdd_edge_stab.
 *
 * `level` must not exceed the target's level. Reading a target from higher up
 * than its edge really starts would add generators for levels the label
 * already acts on, and nothing here can detect that.
 */
static LIMDD_STAB
stab_at_level(uint32_t level, LIMDD_TARG t)
{
    const uint32_t lev = limdd_level(t);
    assert(level <= lev && "an edge may skip levels, never climb them");
    return limdd_stab_extend_skipped(limdd_node_stab(t), level, lev);
}

LIMDD_STAB
limdd_node_stab(LIMDD_TARG p)
{
    if (p == LIMDD_TERMINAL) return LIMDD_STAB_TRIVIAL;

    /* The common case, answered from the node itself: no cache line outside
     * this node is touched, so it costs nothing and contends with nobody. */
    if (limdd_node_stab_is_trivial(p)) return LIMDD_STAB_TRIVIAL;

    const uint64_t cached = limdd_node_stab_raw(p);
    if (cached != 0) return cached;

    /*
     * Recomputing rather than requiring the caller to have filled this in:
     * makeedge builds bottom-up so the children are normally already cached,
     * and a worker that arrives at a node another worker is still finishing
     * would otherwise have to wait. Both compute the same group.
     */
    const uint32_t var = limdd_node_var(p);
    const LIMDD low = limdd_node_low(p);
    const LIMDD high = limdd_node_high(p);
    /* Each child as read at var+1, so over any levels its edge skips. The
     * cache stays keyed on the node alone: a node fixes its children and its
     * level, hence what they skip. */
    const LIMDD_STAB s0 = limdd_edge_is_zero(low)  ? LIMDD_STAB_TRIVIAL
                                                   : stab_at_level(var + 1, limdd_target(low));
    const LIMDD_STAB s1 = limdd_edge_is_zero(high) ? LIMDD_STAB_TRIVIAL
                                                   : stab_at_level(var + 1, limdd_target(high));

    const LIMDD_STAB s = limdd_stab_of_node(var, low, high, s0, s1);
    /*
     * Not cached while the canonical form is off. limdd_stab_of_node needs
     * canonical children, and a node built without them can have isomorphic
     * children stored apart, which makes the group come back too small. Too
     * small is harmless to use once, but cached it would survive into the
     * rebuild -- and there a missing X_var is a missed skip.
     */
    if (high_determinism) limdd_node_set_stab_raw(p, s);
    return s;
}

LIMDD_STAB
limdd_edge_stab(uint32_t level, LIMDD e)
{
    if (limdd_edge_is_zero(e)) return LIMDD_STAB_TRIVIAL;

    const LIMDD_STAB s = stab_at_level(level, limdd_target(e));
    const limdd_pauli_t q = limdd_lim_pauli(limdd_label(e));
    if (limdd_pauli_is_identity(q)) return s;

    /* L v is stabilised by L G L^-1, which for Pauli words is G up to the
     * sign their commutation gives. That covers the skipped levels too: X_j
     * anticommutes with a Z or Y at j and comes back as -X_j, which is right,
     * since that level then holds |0>-|1>. */
    /* One walk, not one per index: limdd_stab_gen restarts at the head. */
    LIMDD_LIM gens[LIMDD_MAX_QUBITS];
    size_t k = 0;
    for (LIMDD_STAB c = s; c != LIMDD_STAB_TRIVIAL; c = limdd_stab_tail(c)) {
        const LIMDD_LIM g = limdd_stab_head(c);
        gens[k++] = (limdd_pauli_commutation_phase(limdd_lim_pauli(g), q) == 0)
                  ? g
                  : limdd_lim_make(limdd_lim_pauli(g), wgt_neg(limdd_lim_weight(g)));
    }
    return limdd_stab_make(gens, k);
}

/**
 * Rebuild the node `t` and everything under it, canonically.
 *
 * Returns the edge makeedge produced, which carries whatever label the
 * normalisation moved upwards; the caller composes its own label onto it. The
 * result is memoised on the node, not the edge, because that label is the only
 * thing two edges into one node differ by.
 */
static LIMDD
canonize_node(LIMDD_TARG t)
{
    if (t == LIMDD_TERMINAL) return limdd_one_edge();

    LIMDD res;
    if (cache_get3(CACHE_LIMDD_CANONIZE, 0, t, 0, &res)) return res;

    const LIMDD low = limdd_node_low(t);
    const LIMDD high = limdd_node_high(t);

    /* Children first: makeedge needs canonical children, and this is what
     * makes that true rather than assumed. */
    const LIMDD lo = limdd_edge_is_zero(low) ? limdd_zero_edge()
        : lim_times_edge(limdd_label(low), canonize_node(limdd_target(low)));
    const LIMDD hi = limdd_edge_is_zero(high) ? limdd_zero_edge()
        : lim_times_edge(limdd_label(high), canonize_node(limdd_target(high)));

    res = limdd_makeedge(limdd_node_var(t), lo, hi);
    cache_put3(CACHE_LIMDD_CANONIZE, 0, t, 0, res);
    return res;
}

LIMDD
limdd_canonize(LIMDD e)
{
    if (limdd_edge_is_zero(e)) return e;

    canon_before = limdd_countnodes(e);

    /*
     * makeedge has to do the full job here whatever the policy says, so the
     * flag is flipped for the duration. The operation cache is cleared
     * afterwards: its entries map edges to edges built under the OLD regime,
     * and a hit would hand back a non-canonical result.
     */
    const bool saved = high_determinism;
    high_determinism = true;
    canon_in_progress = true;

    /*
     * The root's label goes back onto an edge makeedge had already reduced,
     * so it is reduced once more, at level 0: otherwise the root would be the
     * one edge in the diagram still free to carry an unreduced label. An edge
     * straight to the terminal is not exempt, since it may skip every level.
     */
    const LIMDD_TARG t = limdd_target(e);
    const LIMDD res = limdd_edge_canonical(0, t == LIMDD_TERMINAL ? e
        : lim_times_edge(limdd_label(e), canonize_node(t)));

    canon_in_progress = false;
    high_determinism = saved;
    cache_clear();

    canon_after = limdd_countnodes(res);
    canon_counters_reset();

    if (canon_policy == LIMDD_CANON_ADAPTIVE) {
        /*
         * Tune the interval by what the rebuild was worth. A big reduction
         * means the diagram was carrying a lot that the canonical form would
         * have merged, so go sooner; almost no reduction means the rebuild
         * barely paid for itself, so wait longer. Bounded, so neither runs
         * away.
         */
        const double gain = canon_after ? (double)canon_before / canon_after : 1.0;
        if (gain > 2.0 && canon_interval > 1) canon_interval /= 2;
        else if (gain < 1.2 && canon_interval < (1u << 20)) canon_interval *= 2;
    }
    return res;
}

LIMDD
limdd_edge_canonical(uint32_t level, LIMDD e)
{
    if (limdd_edge_is_zero(e)) return limdd_zero_edge();
    if (!high_determinism) return e;   // no groups are maintained

    /* No early-out for the terminal: read above level nqubits it denotes
     * |+>^k, whose group is anything but trivial. */
    const LIMDD_TARG t = limdd_target(e);
    const LIMDD_STAB s = stab_at_level(level, t);
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
    limdd_pauli_t w = limdd_lim_pauli(a);
    assert(limdd_pauli_get(w, var) == LIMDD_PAULI_I
           && "the label already acts at this level");
    limdd_pauli_set(&w, var, (limdd_pauli_op_t)((x ? 2u : 0u) | (z ? 1u : 0u)));
    return limdd_lim_make(w, limdd_lim_weight(a));
}

/**
 * A total order on LIMs, by VALUE.
 *
 * Comparing interned indices instead would be cheaper and wrong. An index
 * records when a label was first seen, not what it is, so the same two labels
 * can compare either way depending on what the program interned earlier -- and
 * a collection that sweeps one of them and hands back a different index when
 * it reappears changes the answer outright. The order below depends on nothing
 * but the label.
 */
static int
lim_cmp(LIMDD_LIM a, LIMDD_LIM b)
{
    if (a == b) return 0;

    const limdd_pauli_t pa = limdd_lim_pauli(a), pb = limdd_lim_pauli(b);
    const int pc = limdd_pauli_cmp(pa, pb);
    if (pc != 0) return pc;

    /* As doubles, which is enough to order them and works for any backend.
     * Exact weights lose precision here, but only the ORDER matters and it
     * stays a deterministic function of the value. */
    const complex_t ca = weight_as_complex(limdd_lim_weight(a));
    const complex_t cb = weight_as_complex(limdd_lim_weight(b));
    if (ca.r != cb.r) return ca.r < cb.r ? -1 : 1;
    if (ca.i != cb.i) return ca.i < cb.i ? -1 : 1;
    return 0;
}

/**
 * Which of the two candidate nodes is kept.
 *
 * Any fixed total order gives a canonical form; this one compares the two
 * children then the label. What matters is that the order depends only on the
 * two candidates, so that taking the smaller lands on the same node however
 * the state was presented and whatever the tables have seen before.
 *
 * Child indices are safe to compare: both nodes already exist when this runs,
 * so their indices are fixed, and a collection leaves them where they are. The
 * labels are not -- they are interned during this very call -- so they go
 * through lim_cmp.
 */
static bool
node_less(LIMDD_TARG lo_a, LIMDD_TARG hi_a, LIMDD_LIM lab_a,
          LIMDD_TARG lo_b, LIMDD_TARG hi_b, LIMDD_LIM lab_b)
{
    if (lo_a != lo_b) return lo_a < lo_b;
    if (hi_a != hi_b) return hi_a < hi_b;
    return lim_cmp(lab_a, lab_b) < 0;
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

    /* Without high determinism no group was needed to build this node, and
     * computing one here would cost exactly what that mode exists to avoid. */
    if (!high_determinism) return t;

    /*
     * Nothing to record when the node already says its group is trivial: the
     * bit is in the node, limdd_node_stab answers from it, and the cache is
     * never consulted for this node again. Skipping the store is the point of
     * that bit -- otherwise every node would still write one word into an
     * array shared by all the workers, and those writes are not confined to
     * the region a worker allocates from, since a worker computes the group
     * of whatever node it descends into.
     */
    if (limdd_node_stab_is_trivial(t)) return t;

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
    /* One unit of work. Counted here rather than per gate, so the trigger
     * measures what the circuit actually costs instead of how many gates it
     * happens to contain. This runs inside the parallel recursions, which is
     * why the counter is per worker. */
    limdd_canon_count();

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
                                    high_determinism ? stab_at_level(var + 1, v)
                                                     : LIMDD_STAB_TRIVIAL,
                                    LIMDD_STAB_TRIVIAL);
        const LIMDD e = limdd_bundle(
            lim_with_pauli_at(limdd_label(live), var, lz, false), t);
        return high_determinism ? limdd_edge_canonical(var, e) : e;
    }

    const LIMDD_TARG v0 = limdd_target(low);
    const LIMDD_TARG v1 = limdd_target(high);

    const LIMDD_LIM a = limdd_label(low);

    /*
     * THE SKIP RULE, the half that needs no group. Equal child edges mean the
     * node would denote |0>(x)w + |1>(x)w = (|0>+|1>)(x)w -- which is what
     * `low` already denotes when read at this level, so the level is simply
     * not stored. Nothing has to be put on the label: `a` has no entry at
     * `var`, and `low` is already reduced at `var`, since extending its group
     * by X_var only adds a pivot in a column no reduced label has set.
     */
    if (low == high) return low;

    /* Divide the low label out; what is left is all the node keeps. */
    const LIMDD_LIM bp = limdd_lim_mul(limdd_lim_inverse(a), limdd_label(high));

    if (!high_determinism) {
        /*
         * The other group-free case: |0>w - |1>w is (|0>-|1>)(x)w, the same
         * skip with Z at this level. Everything else keeps the label as it
         * came out of the division, with no search for its class.
         */
        if (v0 == v1 && limdd_pauli_is_identity(limdd_lim_pauli(bp))
                     && limdd_lim_weight(bp) == EVBDD_MIN_ONE) {
            return limdd_bundle(lim_with_pauli_at(a, var, false, true), v0);
        }
        const LIMDD_TARG t = intern(var, v0, bp, v1, LIMDD_STAB_TRIVIAL,
                                    LIMDD_STAB_TRIVIAL);
        return limdd_bundle(a, t);
    }

    /* The children as read at var+1, over whatever levels they skip. */
    const LIMDD_STAB s0 = stab_at_level(var + 1, v0);
    const LIMDD_STAB s1 = stab_at_level(var + 1, v1);

    /* Candidate as given. */
    LIMDD_LIM g1; bool neg1;
    const LIMDD_LIM lab1 = limdd_stab_min_coset(bp, s0, s1, &g1, &neg1);

    /*
     * THE SKIP RULE, the half that needs the groups. lab1 == I with one
     * target means bp lies in +-Stab(v0): the high branch is +- the low one,
     * the node is (|0> +- |1>)(x)w, and this level is not stored. The sign
     * rule has already resolved the +- to +, recording the flip in neg1, and
     * that flip is exactly a Z at this level -- so this is the unswapped
     * branch's own label formula with the intern() left out. |0>w + i|1>w,
     * a product state but no Pauli image of |0>+|1>, gives lab1 = iI and
     * correctly falls through.
     *
     * It must come BEFORE the swapped candidate: on this input lab1 and lab2
     * tie, node_less then takes the swapped branch, and that deposits XZ at
     * this level rather than Z. Same state, another label, and only the
     * reduction against X_var would bring the two together again.
     */
    if (v0 == v1 && lab1 == LIMDD_LIM_IDENTITY) {
        const LIMDD_LIM r = limdd_lim_mul(a, g1);
        return limdd_edge_canonical(var,
            limdd_bundle(lim_with_pauli_at(r, var, false, neg1), v0));
    }

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
        return limdd_edge_canonical(var, limdd_bundle(lim_with_pauli_at(r, var, false, neg1), t));
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
        return limdd_edge_canonical(var, limdd_bundle(r, t));
    }
}
