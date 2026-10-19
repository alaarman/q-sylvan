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
#include <string.h>

#include <sylvan_int.h>

#include "qsylvan_limdd_canon.h"
#include "qsylvan_limdd_profile.h"
#include "qsylvan_limdd_class.h"
#include "sylvan_edge_weights_complex.h"

/* Declared in qsylvan_limdd_ops.h, which includes this module's header, so it
 * cannot be included from here. */
uint64_t limdd_countnodes(LIMDD e);

/** `l` composed onto `e`'s label: the same target under a composed map. */
static inline LIMDD
lim_times_edge(LIMDD_LIM l, LIMDD e)
{
    if (limdd_lim_is_identity(l)) return e;
    if (limdd_edge_is_zero(e) || limdd_lim_is_zero(l)) return limdd_zero_edge();
    /* On floats a product within the tolerance of 0 is the zero map, and the
     * zero edge points at the terminal, not at e's target. */
    const LIMDD_LIM m = limdd_lim_mul(l, limdd_label(e));
    return limdd_lim_is_zero(m) ? limdd_zero_edge() : limdd_bundle(m, limdd_target(e));
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
    limdd_prof_init();
    /* Deferring means operations run without the orbit search; the lazy form
     * runs it where a comparison needs it, so it keeps it on. */
    high_determinism = (policy == LIMDD_CANON_ALWAYS || policy == LIMDD_CANON_LAZY);
    /* The lazy form's tables come with the node table, or now if that
     * exists; switching it off takes effect with the next session. */
    limdd_class_request(policy == LIMDD_CANON_LAZY);
    limdd_nodes_attach_class();
}

bool
limdd_canon_is_lazy(void)
{
    return limdd_class_on();
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
    if (canon_policy == LIMDD_CANON_ALWAYS || canon_policy == LIMDD_CANON_LAZY ||
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

#ifdef LIMDD_PROFILE

limdd_prof_line_t *limdd_prof = NULL;
static unsigned prof_nworkers = 0;

/*
 * Fingerprint sets for the lazy estimate: open addressing over 64-bit
 * fingerprints, 0 meaning empty. Kept only with one worker (prof_sets_on), so
 * nothing here is shared; when one fills past 3/4 the estimate stops and the
 * report says so.
 */
#define PROF_SET_LOG 23
typedef struct { uint64_t *slot; uint64_t n; } prof_set_t;
static prof_set_t prof_tuples, prof_classes;
static bool prof_sets_on = false, prof_sets_full = false;

static uint64_t
prof_mix(uint64_t x)
{
    x ^= x >> 30; x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27; x *= UINT64_C(0x94d049bb133111eb);
    x ^= x >> 31;
    return x ? x : 1;
}

/** Inserts fp; returns whether it was there already. */
static bool
prof_set_seen(prof_set_t *s, uint64_t fp)
{
    const uint64_t mask = (UINT64_C(1) << PROF_SET_LOG) - 1;
    for (uint64_t i = fp & mask;; i = (i + 1) & mask) {
        if (s->slot[i] == fp) return true;
        if (s->slot[i] == 0) {
            if (s->n >= (mask + 1) / 4 * 3) { prof_sets_full = true; return false; }
            s->slot[i] = fp; s->n++;
            return false;
        }
    }
}

void
limdd_prof_init(void)
{
    if (limdd_prof != NULL) return;
    prof_nworkers = (unsigned)lace_workers();
    if (prof_nworkers == 0) prof_nworkers = 1;
    limdd_prof = (limdd_prof_line_t *)sylvan_alloc_padded(prof_nworkers * sizeof(limdd_prof_line_t));
    if (limdd_prof == NULL) { fprintf(stderr, "sylvan: out of memory for the LIMDD profile\n"); exit(1); }
    memset(limdd_prof, 0, prof_nworkers * sizeof(limdd_prof_line_t));
    if (prof_nworkers == 1) {
        prof_tuples.slot  = (uint64_t *)calloc(UINT64_C(1) << PROF_SET_LOG, sizeof(uint64_t));
        prof_classes.slot = (uint64_t *)calloc(UINT64_C(1) << PROF_SET_LOG, sizeof(uint64_t));
        prof_sets_on = prof_tuples.slot != NULL && prof_classes.slot != NULL;
    }
}

void
limdd_prof_report(FILE *out)
{
    if (limdd_prof == NULL) return;
    static const char *name[LP_N] = {
        "calls", "zero", "dead", "eq", "same", "gen", "gen_triv", "gen_created",
        "gen_exact", "gen_merged", "lz_tseen", "lz_classnew", "lz_collide",
        "lz_collide_found", "mincoset", "mincoset_hit", "zassenhaus", "stab_of_node",
        "lazy_known", "lazy_resolve", "lazy_member", "lazy_filtered", "lazy_alias"
    };
    for (int k = 0; k < LP_N; k++) {
        uint64_t t = 0;
        for (unsigned i = 0; i < prof_nworkers; i++) t += limdd_prof[i].c[k];
        fprintf(out, "LIMPROF %s %llu\n", name[k], (unsigned long long)t);
    }
    fprintf(out, "LIMPROF lazy_estimate %s\n", !prof_sets_on ? "off (more than one worker)"
                                             : prof_sets_full ? "stopped (a set filled)" : "complete");
}

/*
 * The lazy estimate for a node with two different children: whether this
 * factored node (var, v0, bp, v1) was asked for before (a lazy makeedge would
 * find it in the node table), else whether any node with its class key
 * (var, {v0, v1}) was (it would have to be compared), else neither (it would be
 * stored without its canonical form). Returns whether it was a collision.
 */
static bool
prof_lazy(uint32_t var, LIMDD_TARG v0, LIMDD_LIM bp, LIMDD_TARG v1)
{
    if (!prof_sets_on || prof_sets_full) return false;
    const uint64_t lo = v0 < v1 ? v0 : v1, hi = v0 < v1 ? v1 : v0;
    const uint64_t cls = prof_mix(((uint64_t)var << 48) ^ prof_mix(lo << 32 | hi));
    const uint64_t tup = prof_mix(cls ^ prof_mix((uint64_t)bp << 32 | v0));
    const bool tseen = prof_set_seen(&prof_tuples, tup);
    const bool cseen = prof_set_seen(&prof_classes, cls);
    if (tseen) { LIMDD_PROF(LP_LZ_TSEEN); return false; }
    if (!cseen) { LIMDD_PROF(LP_LZ_CLASSNEW); return false; }
    LIMDD_PROF(LP_LZ_COLLIDE);
    return true;
}

#else

void limdd_prof_init(void) {}
void limdd_prof_report(FILE *out) { (void)out; }

#endif

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
    if (limdd_node_stab_is_trivial(p)) {
#ifndef NDEBUG
        /* The bit must say what the group computation says. */
        const uint32_t v = limdd_node_var(p);
        const LIMDD lo = limdd_node_low(p), hi = limdd_node_high(p);
        assert(limdd_stab_of_node(v, lo, hi, stab_at_level(v + 1, limdd_target(lo)),
                                  stab_at_level(v + 1, limdd_target(hi))) == LIMDD_STAB_TRIVIAL
               && "a node's trivial-group bit disagrees with its group");
#endif
        return LIMDD_STAB_TRIVIAL;
    }

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
    LIMDD lo = limdd_edge_is_zero(low) ? limdd_zero_edge()
        : lim_times_edge(limdd_label(low), canonize_node(limdd_target(low)));
    LIMDD hi = limdd_edge_is_zero(high) ? limdd_zero_edge()
        : lim_times_edge(limdd_label(high), canonize_node(limdd_target(high)));
    if (limdd_l2_on) {
        /* the node is read as its children scaled to norm 1, as limdd_cofactors
         * reads it; makeedge puts the factor back on the edge it returns */
        const double f = limdd_l2_node_factor(t);
        lo = limdd_scale_real(lo, f);
        hi = limdd_scale_real(hi, f);
    }

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
 * An order on LIMs, by VALUE.
 *
 * Comparing interned indices instead would be cheaper and wrong. An index
 * records when a label was first seen, not what it is, so the same two labels
 * can compare either way depending on what the program interned earlier -- and
 * a collection that sweeps one of them and hands back a different index when
 * it reappears changes the answer outright. The order below depends on nothing
 * but the label.
 *
 * On float weights a component decides only when it differs by more than the
 * weight table's tolerance; within it the next one does. Rounding leaves
 * residues such as 1e-17 where the exact value is 0, and compared exactly
 * their sign would pick the candidate: one class then gets two canonical
 * nodes, and the weight table cannot merge the two labels either, since its
 * grid has a cell boundary at 0. Near such ties the order is not transitive,
 * but makeedge only ever compares its two candidates. Only the lazy form does
 * this (with exact weights it compares exactly too), so the eager form keeps
 * every choice it made before.
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
    double tol = 0.0;
    if (limdd_class_on() && sylvan_get_edge_weight_type() == WGT_COMPLEX_128) {
        /* relative to the larger modulus, see labels_match */
        tol = sylvan_edge_weights_tolerance() * fmax(hypot((double)ca.r, (double)ca.i),
                                                     hypot((double)cb.r, (double)cb.i));
    }
    if (fabs((double)ca.r - (double)cb.r) > tol) return ca.r < cb.r ? -1 : 1;
    if (fabs((double)ca.i - (double)cb.i) > tol) return ca.i < cb.i ? -1 : 1;
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
       LIMDD_STAB s_lo, LIMDD_STAB s_hi, int *created_out)
{
    const LIMDD low = limdd_bundle(LIMDD_LIM_IDENTITY, lo);
    const LIMDD high = limdd_bundle(lab, hi);
    int created;
    const LIMDD_TARG t = limdd_makenode_ex(var, low, high, &created);
    if (created_out != NULL) *created_out = created;

    /* Without high determinism no group was needed to build this node, and
     * computing one here would cost exactly what that mode exists to avoid.
     * The lazy form computes groups only where a comparison asks for them
     * (limdd_node_stab, on demand). */
    if (!high_determinism || limdd_class_on()) return t;

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

/**
 * In the lazy form on float weights, a canonical label with the components
 * that are within the weight tolerance of 0 (relative to its modulus) set to
 * +0: rounding leaves residues such as 1e-17 or -0 where the exact value is 0,
 * and the weight table, whose grid has a cell boundary at 0, would intern
 * them apart, so the canonical node of one class would come out twice. The
 * label moves by less than the tolerance, which is what the table merges by.
 * The eager form keeps its labels as they were.
 */
static LIMDD_LIM
lim_clean(LIMDD_LIM l)
{
    if (!limdd_class_on() || sylvan_get_edge_weight_type() != WGT_COMPLEX_128) return l;
    const complex_t c = weight_as_complex(limdd_lim_weight(l));
    const double tol = sylvan_edge_weights_tolerance() * hypot((double)c.r, (double)c.i);
    const bool zr = fabs((double)c.r) <= tol, zi = fabs((double)c.i) <= tol;
    if ((!zr || ((double)c.r == 0.0 && !signbit((double)c.r))) &&
        (!zi || ((double)c.i == 0.0 && !signbit((double)c.i)))) return l;
    return limdd_lim_make(limdd_lim_pauli(l),
                          complex_lookup(zr ? 0.0 : (double)c.r, zi ? 0.0 : (double)c.i));
}

/*
 * THE CANONICAL NODE of the factored node (var, I->v0, bp->v1): which child goes
 * low, the minimal label, and what undoes them, as canon_root assembles it.
 * Shared by the eager form, which interns the result, and the lazy one, which
 * only compares it; both need the children's groups.
 */
typedef struct {
    LIMDD_TARG lo, hi;      /* the canonical node (var, I->lo, lab->hi) */
    LIMDD_LIM  lab;
    LIMDD_STAB s_lo, s_hi;  /* the groups of lo and hi as read at var+1 */
    LIMDD_LIM  g;           /* the minimisation's witness */
    bool       neg;         /* the sign rule flipped the label */
    bool       swapped;     /* lo is v1: the branches were swapped */
    bool       skip;        /* v0 == v1 and the level is not stored at all */
} canon_key_t;

static void
canon_key(uint32_t var, LIMDD_TARG v0, LIMDD_LIM bp, LIMDD_TARG v1, canon_key_t *k)
{
    /* The children as read at var+1, over whatever levels they skip. */
    const LIMDD_STAB s0 = stab_at_level(var + 1, v0);
    const LIMDD_STAB s1 = stab_at_level(var + 1, v1);

    /* Candidate as given. */
    LIMDD_LIM g1; bool neg1;
    const LIMDD_LIM lab1 = lim_clean(limdd_stab_min_coset(bp, s0, s1, &g1, &neg1));

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
    k->skip = v0 == v1 && lab1 == LIMDD_LIM_IDENTITY;
    if (k->skip) {
        k->lo = k->hi = v0; k->lab = lab1; k->s_lo = s0; k->s_hi = s1;
        k->g = g1; k->neg = neg1; k->swapped = false;
        return;
    }

    /*
     * Candidate with the branches swapped. Applying X at this level gives
     * |0>(x)B'|v1> + |1>(x)|v0>; dividing B' out in turn leaves the node
     * (v1, B'^-1, v0), so the swap costs one inverse and a second minimisation
     * over the same two groups in the other order.
     */
    /*
     * Computed only when it can matter. node_less compares the two targets
     * first and reaches the label only when they are equal, so with v0 < v1
     * the unswapped candidate wins on the targets alone and this second
     * minimisation is never read. Measured over five million calls, the two
     * always agree on the Pauli WORD anyway -- a Pauli is its own inverse and
     * both reduce by the same pivot set -- and differ only in the scalar,
     * which is exactly what the label comparison looks at when v0 == v1.
     */
    LIMDD_LIM g2 = LIMDD_LIM_IDENTITY; bool neg2 = false;
    LIMDD_LIM lab2 = LIMDD_LIM_IDENTITY;
    /*
     * On floats the inverse of a label above 1/tolerance is merged with 0, and
     * the swapped node would then have the zero map on a live branch, which
     * the table cannot tell from a dead one. Such a node is kept in the order
     * it came in. Dropping the low branch instead would be wrong: the children
     * are not unit vectors under low normalisation, so a large label does not
     * make the low branch small.
     */
    bool can_swap = true;
    if (v0 >= v1) {
        const LIMDD_LIM bpinv = limdd_lim_inverse(bp);
        if (limdd_lim_is_zero(bpinv)) can_swap = false;
        else lab2 = lim_clean(limdd_stab_min_coset(bpinv, s1, s0, &g2, &neg2));
    }

    if (!can_swap || node_less(v0, v1, lab1, v1, v0, lab2)) {
        k->lo = v0; k->hi = v1; k->lab = lab1; k->s_lo = s0; k->s_hi = s1;
        k->g = g1; k->neg = neg1; k->swapped = false;
    } else {
        k->lo = v1; k->hi = v0; k->lab = lab2; k->s_lo = s1; k->s_hi = s0;
        k->g = g2; k->neg = neg2; k->swapped = true;
    }
}

/**
 * The root label r with  a|node> = r |canonical node>, where `node` is the
 * factored node (var, I->v0, bp->v1) that `k` came from and `a` the low label
 * that was divided out of it.
 */
static LIMDD_LIM
canon_root(uint32_t var, LIMDD_LIM a, LIMDD_LIM bp, const canon_key_t *k)
{
    if (!k->swapped) {
        /*
         * Undoing the minimisation costs G on the branches above and, when the
         * sign rule flipped the label, a Z at this level -- which is exactly
         * the operator that negates the high branch and leaves the low one.
         */
        const LIMDD_LIM r = limdd_lim_mul(a, k->g);
        return lim_with_pauli_at(r, var, false, k->neg);
    }
    LIMDD_LIM r = limdd_lim_mul(limdd_lim_mul(a, bp), k->g);
    if (k->neg) {
        /* X here from the swap, Z from the sign rule, and XZ = -iY: the
         * factor cannot be dropped, it is part of the map. */
        r = lim_with_pauli_at(r, var, true, true);
        r = limdd_lim_make(limdd_lim_pauli(r),
                           wgt_mul(limdd_lim_weight(r), limdd_wgt_i_pow(3)));
    } else {
        r = lim_with_pauli_at(r, var, true, false);
    }
    return r;
}

/*
 * THE LAZY FORM (LIMDD_CANON_LAZY; qsylvan_limdd_class.h has the design).
 */

/** The canonical node of stored node m, memoised on m (cleared by a collection). */
static void
member_key(uint32_t var, LIMDD_TARG m, canon_key_t *k)
{
    const LIMDD low = limdd_node_low(m), high = limdd_node_high(m);
    const LIMDD_TARG v0 = limdd_target(low), v1 = limdd_target(high);
    uint64_t r1, r2;
    if (cache_get6(CACHE_LIMDD_CLASSKEY, m, 0, 0, 0, 0, &r1, &r2)) {
        k->swapped = (r1 >> 32) & 1;
        k->lab = (LIMDD_LIM)(r1 & 0xffffffff);
        k->g = (LIMDD_LIM)(r2 & 0xffffffff);
        k->neg = (r2 >> 32) & 1;
        k->lo = k->swapped ? v1 : v0;
        k->hi = k->swapped ? v0 : v1;
        k->skip = false;
        k->s_lo = k->s_hi = LIMDD_STAB_TRIVIAL;   /* not read by the comparison */
        return;
    }
    canon_key(var, v0, limdd_label(high), v1, k);
    cache_put6(CACHE_LIMDD_CLASSKEY, m, 0, 0, 0, 0,
               (uint64_t)k->lab | (uint64_t)k->swapped << 32,
               (uint64_t)k->g | (uint64_t)k->neg << 32);
}

static double
lim_modulus(LIMDD_LIM l)
{
    const complex_t c = weight_as_complex(limdd_lim_weight(l));
    return hypot((double)c.r, (double)c.i);
}

/**
 * Whether two canonical labels name one LIM. Interned labels are equal when
 * their indices are, and on exact weights only then. On floats the canonical
 * labels of two members of a class are computed from different inputs and can
 * round a few ulps apart; the weight table merges such values, unless they
 * fall either side of a cell of its grid. So floats compare the Pauli word
 * exactly and the scalar within the table's tolerance: the lazy form then
 * merges what the table would, boundary or not.
 */
static bool
labels_match(LIMDD_LIM x, LIMDD_LIM y)
{
    if (x == y) return true;
    if (sylvan_get_edge_weight_type() != WGT_COMPLEX_128) return false;
    if (limdd_lim_pauli_ref(x) != limdd_lim_pauli_ref(y)) return false;
    const complex_t a = weight_as_complex(limdd_lim_weight(x));
    const complex_t b = weight_as_complex(limdd_lim_weight(y));
    /* Relative: under low normalisation labels span many decades, and an
     * absolute tolerance would call any two small ones equal. */
    const double scale = fmax(hypot((double)a.r, (double)a.i), hypot((double)b.r, (double)b.i));
    return hypot((double)a.r - (double)b.r, (double)a.i - (double)b.i)
           <= sylvan_edge_weights_tolerance() * scale;
}

/**
 * The class of the stored factored node t = (var, I->v0, bp->v1), v0 != v1:
 * (I, t) if t is its class's member, else (map, m) with |t> = map |m>. Every
 * worker that asks gets the same answer, so it is recorded with a plain store.
 */
static LIMDD
class_resolve(uint32_t var, LIMDD_TARG t, LIMDD_TARG v0, LIMDD_LIM bp, LIMDD_TARG v1)
{
    LIMDD_PROF(LP_LZ_RESOLVE);
    const LIMDD_TARG lo = v0 < v1 ? v0 : v1, hi = v0 < v1 ? v1 : v0;
    /* The filter needs the scalars as doubles, which on exact weights is a
     * conversion of rationals each time; there the canonical keys compare in
     * one step anyway, so it is used on floats only. */
    const bool filter = sylvan_get_edge_weight_type() == WGT_COMPLEX_128;
    const double mt = filter ? lim_modulus(bp) : 0.0;
    /* On exact weights the scalar up to units: where the scan starts, and an
     * exact filter, since every member of the class has the same one. */
    const uint64_t ut = limdd_class_scalar_key(bp);
    canon_key_t kt;
    bool have_kt = false;
    LIMDD res = 0;
    for (uint32_t j = 0;; j++) {
        int claimed;
        const LIMDD_TARG m = limdd_class_slot(var, lo, hi, ut, j, t, &claimed);
        if (m == t) { res = limdd_bundle(LIMDD_LIM_IDENTITY, t); break; }
        LIMDD_PROF(LP_LZ_MEMBER);

        /*
         * The modulus of the high scalar is shared by the class, inverted when
         * the member has its children the other way round; a LIM image only
         * multiplies it by scalars of modulus 1. The tolerance is far wider
         * than any rounding, so the filter never parts two members of a class.
         */
        const LIMDD mhigh = limdd_node_high(m);
        if (ut != 0 && limdd_target(limdd_node_low(m)) == v0 &&
            limdd_class_scalar_key(limdd_label(mhigh)) != ut) {
            LIMDD_PROF(LP_LZ_FILTERED); continue;
        }
        if (filter) {
            const double mm = lim_modulus(limdd_label(mhigh));
            const double want = limdd_target(limdd_node_low(m)) == v0 ? mm : 1.0 / mm;
            if (fabs(mt - want) > 1e-6 * fmax(mt, want)) { LIMDD_PROF(LP_LZ_FILTERED); continue; }
        }

        if (!have_kt) { canon_key(var, v0, bp, v1, &kt); have_kt = true; }
        canon_key_t km;
        member_key(var, m, &km);
        if (kt.lo != km.lo || kt.hi != km.hi || !labels_match(kt.lab, km.lab)) continue;

        /* |t> = rt |C> and |m> = rm |C>, so |t> = rt rm^-1 |m>. */
        const LIMDD_LIM rt = canon_root(var, LIMDD_LIM_IDENTITY, bp, &kt);
        const LIMDD_LIM rm = canon_root(var, LIMDD_LIM_IDENTITY, limdd_label(mhigh), &km);
        const LIMDD_LIM rminv = limdd_lim_inverse(rm);   /* see bp in makeedge_canonical */
        res = limdd_bundle(limdd_lim_is_zero(rminv) ? limdd_lim_rdiv(rt, rm)
                                                    : limdd_lim_mul(rt, rminv), m);
        LIMDD_PROF(LP_LZ_ALIAS);
        break;
    }
    limdd_node_set_class_raw(t, (uint64_t)res);
    return res;
}

/**
 * The factored node (var, I->v0, bp->v1), v0 != v1, in the form it is stored in
 * lazily: canon_key with both children's groups taken as trivial, so with no
 * group computed. The children go in index order and the sign rule fixes the
 * label's sign; canon_root undoes both, as for the canonical node. When both
 * groups are trivial this IS the canonical node: a LIM between two such nodes
 * with the same children in the same order is then +-1 at the top and the
 * identity below (qsylvan_limdd_class.h), so the class is {(v0, +-b, v1),
 * (v1, +-b^-1, v0)} and this picks one of the four. Isomorphic presentations
 * then meet in the node table, without a comparison and without a node for
 * each of them. Otherwise it only fixes the order and the sign, and the class
 * table does the rest; members are stored this way, which is what lets the
 * class key use the scalar.
 */
static void
lazy_key(LIMDD_TARG v0, LIMDD_LIM bp, LIMDD_TARG v1, canon_key_t *k)
{
    assert(v0 != v1);
    k->skip = false;
    k->s_lo = k->s_hi = LIMDD_STAB_TRIVIAL;
    k->swapped = v0 > v1;
    k->lo = k->swapped ? v1 : v0;
    k->hi = k->swapped ? v0 : v1;
    const LIMDD_LIM b = k->swapped ? limdd_lim_inverse(bp) : bp;
    k->lab = lim_clean(limdd_stab_min_coset(b, LIMDD_STAB_TRIVIAL, LIMDD_STAB_TRIVIAL,
                                            &k->g, &k->neg));
}

/**
 * makeedge for two different children, lazily: the factored node is interned
 * in the order and sign of lazy_key, and handed out if it is its class's
 * member. Only when another node may share its class is anything canonical
 * computed, and then only to compare.
 *
 * On floats the factored node is interned as it comes. Ordering it costs the
 * swapped half an inverse, and its edge a product with the label, each of
 * which rounds; on the 79 amendment circuits that drift kept more members of
 * a class apart than the normal form brought together (60_950: a peak of
 * 903672 nodes against 17077, grover_n30 filling its table), and the scalar
 * key that the order serves is exact-only anyway.
 */
static LIMDD
makeedge_lazy(uint32_t var, LIMDD_LIM a, LIMDD_TARG v0, LIMDD_LIM bp, LIMDD_TARG v1)
{
    if (sylvan_get_edge_weight_type() != WGT_QISQ2) {
        int created;
        const LIMDD_TARG t = intern(var, v0, bp, v1, LIMDD_STAB_TRIVIAL, LIMDD_STAB_TRIVIAL, &created);
        uint64_t st = created ? 0 : limdd_node_class_raw(t);
        if (st == 0) st = (uint64_t)class_resolve(var, t, v0, bp, v1);
        else LIMDD_PROF(LP_LZ_KNOWN);
        return lim_times_edge(a, (LIMDD)st);
    }
    canon_key_t k;
    lazy_key(v0, bp, v1, &k);
    int created;
    const LIMDD_TARG t = intern(var, k.lo, k.lab, k.hi, LIMDD_STAB_TRIVIAL, LIMDD_STAB_TRIVIAL, &created);
    uint64_t st = created ? 0 : limdd_node_class_raw(t);
    if (st == 0) st = (uint64_t)class_resolve(var, t, k.lo, k.lab, k.hi);
    else LIMDD_PROF(LP_LZ_KNOWN);
    return lim_times_edge(canon_root(var, a, bp, &k), (LIMDD)st);
}

static LIMDD makeedge_canonical(uint32_t var, LIMDD low, LIMDD high);

LIMDD
limdd_makeedge(uint32_t var, LIMDD low, LIMDD high)
{
    const LIMDD e = makeedge_canonical(var, low, high);
    if (!limdd_l2_on || limdd_edge_is_zero(e)) return e;

    /*
     * The node denotes its vector at norm 1 under L2 (see qsylvan_limdd_node.h),
     * so the edge carries the norm that was divided out: sqrt(1 + |b|^2) for
     * the node made at `var`, b the scalar of its high label, or sqrt(2) when
     * the level was skipped and the edge points below `var`. The canonical
     * form chose the node and the label's Pauli word, and a positive factor
     * changes neither, so the edge is still reduced.
     */
    const LIMDD_TARG t = limdd_target(e);
    const double f = limdd_level(t) > var ? M_SQRT2 : 1.0 / limdd_l2_node_factor(t);
    return limdd_scale_real(e, f);
}

/** limdd_makeedge without the L2 factor: the canonical node, and the label
 *  that makes up the difference with the low-factored node. */
static LIMDD
makeedge_canonical(uint32_t var, LIMDD low, LIMDD high)
{
    /* One unit of work. Counted here rather than per gate, so the trigger
     * measures what the circuit actually costs instead of how many gates it
     * happens to contain. This runs inside the parallel recursions, which is
     * why the counter is per worker. */
    limdd_canon_count();
    LIMDD_PROF(LP_CALLS);
    const bool lazy = limdd_class_on();

    const bool lz = limdd_edge_is_zero(low);
    const bool hz = limdd_edge_is_zero(high);
    if (lz && hz) { LIMDD_PROF(LP_ZERO); return limdd_zero_edge(); }

    /*
     * One dead branch. |0>(x)v and |1>(x)v differ by X at this level, so they
     * are one node, and the live branch is always put on the low side. This is
     * what makes |0..0> and |10..0> share every node they can.
     */
    if (lz || hz) {
        LIMDD_PROF(LP_DEAD);
        const LIMDD live = lz ? high : low;
        const LIMDD_TARG v = limdd_target(live);
        const LIMDD_TARG t = intern(var, v, LIMDD_LIM_ZERO, LIMDD_TERMINAL,
                                    high_determinism && !lazy ? stab_at_level(var + 1, v)
                                                              : LIMDD_STAB_TRIVIAL,
                                    LIMDD_STAB_TRIVIAL, NULL);
        const LIMDD e = limdd_bundle(
            lim_with_pauli_at(limdd_label(live), var, lz, false), t);
        /* Lazily the root label is not reduced; the node is unique anyway. */
        return high_determinism && !lazy ? limdd_edge_canonical(var, e) : e;
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
    if (low == high) { LIMDD_PROF(LP_EQ); return low; }

    /*
     * Divide the low label out; what is left is all the node keeps. Through the
     * inverse, as every other quotient of the canonical form is taken: on
     * floats a single division rounds differently, and that alone kept
     * isomorphic nodes apart (on 40_700 eager float 18383 nodes against 593).
     * Only when the inverse is merged with 0 -- a label above 1/tolerance --
     * is the quotient taken in one division, which is then the only way to
     * get it.
     */
    const LIMDD_LIM ainv = limdd_lim_inverse(a);
    const LIMDD_LIM bp = limdd_lim_is_zero(ainv) ? limdd_lim_ldiv(a, limdd_label(high))
                                                 : limdd_lim_mul(ainv, limdd_label(high));

    /*
     * On floats the quotient is the zero map when it is within the weight
     * table's tolerance of 0: the table has merged the high branch's scalar,
     * relative to the low one, with 0, as an EVDD's low normalisation does, so
     * the node already has a dead high branch. Built as one, it is the
     * canonical dead-branch node. With the zero label kept, the high edge would
     * point at a node rather than the terminal, a second form of one node,
     * and the coset search would be handed the zero map.
     */
    if (limdd_lim_is_zero(bp)) return makeedge_canonical(var, low, limdd_zero_edge());

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
                                    LIMDD_STAB_TRIVIAL, NULL);
        return limdd_bundle(a, t);
    }

    if (lazy && v0 != v1 && limdd_class_wanted(var, v0, v1)) return makeedge_lazy(var, a, v0, bp, v1);

    canon_key_t k;
    canon_key(var, v0, bp, v1, &k);

#ifdef LIMDD_PROFILE
    bool prof_collide = false;
    if (v0 == v1) {
        LIMDD_PROF(LP_SAME);
    } else {
        LIMDD_PROF(LP_GEN);
        if (k.s_lo == LIMDD_STAB_TRIVIAL && k.s_hi == LIMDD_STAB_TRIVIAL) LIMDD_PROF(LP_GEN_TRIV);
        prof_collide = prof_lazy(var, v0, bp, v1);
    }
#endif

    /* The skip rule (see canon_key): the level is not stored. */
    if (k.skip) {
        const LIMDD_LIM r = limdd_lim_mul(a, k.g);
        const LIMDD e = limdd_bundle(lim_with_pauli_at(r, var, false, k.neg), v0);
        return lazy ? e : limdd_edge_canonical(var, e);
    }

    int created;
    const LIMDD_TARG t = intern(var, k.lo, k.lab, k.hi, k.s_lo, k.s_hi, &created);
#ifdef LIMDD_PROFILE
    if (v0 != v1) {
        LIMDD_PROF(created ? LP_GEN_CREATED : (!k.swapped && k.lab == bp) ? LP_GEN_EXACT : LP_GEN_MERGED);
        if (prof_collide && !created) LIMDD_PROF(LP_LZ_COLLIDE_FOUND);
    }
#else
    (void)created;
#endif
    const LIMDD e = limdd_bundle(canon_root(var, a, bp, &k), t);
    return lazy ? e : limdd_edge_canonical(var, e);
}
