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
#include <stdatomic.h>
#include <stdlib.h>

#include <sylvan_int.h>
#include <sylvan_platform.h>

#include "qsylvan_limdd_node.h"
#include "qsylvan_limdd_stab.h"

/*
 * Written once by limdd_nodes_init and read-only afterwards, so there is no
 * line here for two workers to bounce between caches.
 */
static llmsset_t limdd_nodes   = NULL;
static size_t    limdd_nqubits = 0;

/* One slot per node bucket; see limdd_node_stab_raw in the header. */
static _Atomic(uint64_t) *node_stab = NULL;

static LIMDD zero_edge = 0;
static LIMDD one_edge  = 0;

/*
 * Field masks, for the layout in the header (NODES). The variable field,
 * bits 47..62 of the low word, sits where EVBDD has its own, but the two
 * layouts do not agree. Bits 32..45 of the low word are the flags of the
 * diagram that reads the node, which the LIMDD leaves at zero and the BQD's
 * rule SM sets (bit 45 for its tag S, bit 44 for full; qsylvan_bqd_sm.h),
 * and bit 45 is EVBDD's wgt_val flag. The high word is a LIM and a target,
 * all 64 bits, and EVBDD keeps its mark in bit 63 of it. A walker that read
 * LIMDD buckets as EVBDD nodes would misread both; none does, since only
 * this file reads them, and collection marks in llmsset's own bitmap, not in
 * the bucket.
 *
 * Bit 63 of the low word is set iff this node's stabiliser group is provably
 * trivial by the cheap test in stab_trivially_trivial() below.
 */
static const uint64_t limdd_stab_triv_mask = UINT64_C(0x8000000000000000);
static const uint64_t limdd_var_mask      = UINT64_C(0x7fff800000000000); // bits 47..62
static const uint64_t limdd_low_zero_mask = UINT64_C(0x0000400000000000); // bit 46
static const uint64_t limdd_flags_mask    = UINT64_C(0x00003fff00000000); // bits 32..45
static const uint64_t limdd_lim_mask      = UINT64_C(0xffffffff00000000); // bits 32..63
static const uint64_t limdd_targ_mask     = UINT64_C(0x00000000ffffffff); // bits 0..31

#define LIMDD_VAR_SHIFT 47
#define LIMDD_FLAGS_SHIFT 32

/** A node as it sits in the table: exactly one 16-byte llmsset bucket. */
typedef struct __attribute__((packed)) limddnode {
    uint64_t low, high;
} *limddnode_t;

static inline limddnode_t
limdd_getnode(LIMDD_TARG p)
{
    assert(limdd_nodes != NULL);
    assert(p != LIMDD_TERMINAL);
    return (limddnode_t) llmsset_index_to_ptr(limdd_nodes, p);
}

bool
limdd_node_stab_is_trivial(LIMDD_TARG p)
{
    if (p == LIMDD_TERMINAL) return true;   /* the empty product fixes nothing */
    return (limdd_getnode(p)->low & limdd_stab_triv_mask) != 0;
}

/**
 * Is this node's stabiliser group trivial, by a test that costs nothing?
 *
 * Reading limdd_stab_of_node: with neither branch dead, the only generators
 * are the diagonal part -- the elements of s0 whose word lies in s1's span --
 * and one anti-diagonal representative, which exists only when the two child
 * edges point at the same node. If EITHER child group is trivial the diagonal
 * part is empty: with s0 trivial there is nothing to draw from, and with s1
 * trivial its span is {0}, which admits only the identity word, and no
 * stabiliser group holds a non-identity element with the identity word. So
 * the group is exactly trivial when neither branch is dead, the targets
 * differ, and one child contributes nothing.
 *
 * A child contributes nothing when its own group is trivial AND its edge
 * skips no level: a skipped level adds an X_j generator (see
 * limdd_stab_extend_skipped), which would survive the intersection.
 *
 * Conservative in one direction only -- when both children have non-trivial
 * groups this answers false and the caller computes properly -- and it is a
 * pure function of the node's contents, so every caller that builds the same
 * node derives the same bit and hash consing is unaffected.
 */
static bool
stab_trivially_trivial(uint32_t var, LIMDD low, LIMDD high)
{
    if (limdd_lim_is_zero(limdd_label(low)) ||
        limdd_lim_is_zero(limdd_label(high))) return false;   /* a dead branch adds Z */

    const LIMDD_TARG lt = limdd_target(low), ht = limdd_target(high);
    if (lt == ht) return false;                               /* anti-diagonal may fire */

    return (limdd_node_stab_is_trivial(lt) && limdd_level(lt) == var + 1)
        || (limdd_node_stab_is_trivial(ht) && limdd_level(ht) == var + 1);
}

uint64_t
limdd_node_stab_raw(LIMDD_TARG p)
{
    assert(node_stab != NULL);
    return atomic_load_explicit(&node_stab[p], memory_order_relaxed);
}

void
limdd_node_set_stab_raw(LIMDD_TARG p, uint64_t v)
{
    assert(node_stab != NULL);
    atomic_store_explicit(&node_stab[p], v, memory_order_relaxed);
}

LIMDD
limdd_zero_edge(void)
{
    return zero_edge;
}

LIMDD
limdd_one_edge(void)
{
    return one_edge;
}

uint32_t
limdd_node_var(LIMDD_TARG p)
{
    return (uint32_t)((limdd_getnode(p)->low & limdd_var_mask) >> LIMDD_VAR_SHIFT);
}

uint32_t
limdd_level(LIMDD_TARG p)
{
    return p == LIMDD_TERMINAL ? (uint32_t)limdd_nqubits : limdd_node_var(p);
}

uint32_t
limdd_node_flags(LIMDD_TARG p)
{
    return (uint32_t)((limdd_getnode(p)->low & limdd_flags_mask) >> LIMDD_FLAGS_SHIFT);
}

LIMDD
limdd_node_low(LIMDD_TARG p)
{
    const limddnode_t n = limdd_getnode(p);
    /*
     * The low label is never stored: normalisation leaves it either the
     * identity or zero, and one bit says which. Rebuilding the edge here is
     * what lets the rest of the code treat both children uniformly.
     */
    const LIMDD_LIM lim = (n->low & limdd_low_zero_mask) ? LIMDD_LIM_ZERO
                                                         : LIMDD_LIM_IDENTITY;
    return limdd_bundle(lim, n->low & limdd_targ_mask);
}

LIMDD
limdd_node_high(LIMDD_TARG p)
{
    const limddnode_t n = limdd_getnode(p);
    return limdd_bundle((n->high & limdd_lim_mask) >> LIMDD_TARG_BITS,
                        n->high & limdd_targ_mask);
}

LIMDD_TARG
limdd_makenode(uint32_t var, LIMDD low, LIMDD high)
{
    int created;
    return limdd_makenode_ex(var, low, high, &created);
}

LIMDD_TARG
limdd_makenode_ex(uint32_t var, LIMDD low, LIMDD high, int *created)
{
    return limdd_makenode_flags_ex(var, low, high, 0, created);
}

LIMDD_TARG
limdd_makenode_flags(uint32_t var, LIMDD low, LIMDD high, uint32_t flags)
{
    int created;
    return limdd_makenode_flags_ex(var, low, high, flags, &created);
}

LIMDD_TARG
limdd_makenode_flags_ex(uint32_t var, LIMDD low, LIMDD high, uint32_t flags, int *created)
{
    assert(limdd_nodes != NULL);
    assert(var < limdd_nqubits);
    assert(flags < (UINT32_C(1) << LIMDD_NODE_FLAG_BITS));

    const LIMDD_LIM low_lim  = limdd_label(low);
    const LIMDD_LIM high_lim = limdd_label(high);

    /*
     * The low edge carries no label of its own. A caller that has not factored
     * it out to the parent would silently lose it here, so this is an error
     * rather than something to normalise away: only the caller knows what to
     * multiply the parent edge by.
     */
    assert(low_lim == LIMDD_LIM_IDENTITY || low_lim == LIMDD_LIM_ZERO);
    assert(high_lim < LIMDD_LIM_MAX && "LIM table exceeds the edge's LIM field");

    /*
     * A zero edge denotes the zero vector whatever it points at, so the target
     * is pinned to the terminal. Without this the same node could be written
     * as many ways as there are nodes to dangle a dead branch off, and lookup
     * would stop finding nodes that already exist.
     */
    LIMDD_TARG low_targ  = limdd_lim_is_zero(low_lim)  ? LIMDD_TERMINAL : limdd_target(low);
    LIMDD_TARG high_targ = limdd_lim_is_zero(high_lim) ? LIMDD_TERMINAL : limdd_target(high);

    assert(low_targ <= LIMDD_TARG_MAX && high_targ <= LIMDD_TARG_MAX);

    /* A child may skip levels but must lie strictly below this node. */
    assert(limdd_lim_is_zero(low_lim)  || var < limdd_level(low_targ));
    assert(limdd_lim_is_zero(high_lim) || var < limdd_level(high_targ));

    /*
     * No `if (low == high) return low` here, and there never will be: see the
     * header. Skipping a level is limdd_makeedge's decision, made once the low
     * label is factored out and the high label is canonical, and it is the
     * parent edge that carries the skipped level's Pauli. makenode stores what
     * it is handed.
     */

    struct limddnode n;
    n.low  = ((uint64_t)var << LIMDD_VAR_SHIFT)
           | (limdd_lim_is_zero(low_lim) ? limdd_low_zero_mask : 0)
           | (stab_trivially_trivial(var, low, high) ? limdd_stab_triv_mask : 0)
           | (((uint64_t)flags << LIMDD_FLAGS_SHIFT) & limdd_flags_mask)
           | low_targ;
    n.high = ((uint64_t)high_lim << LIMDD_TARG_BITS) | high_targ;

    const uint64_t res = llmsset_lookup(limdd_nodes, n.low, n.high, created);
    if (res == 0) {
        fprintf(stderr, "sylvan: LIMDD node table is full\n");
        exit(1);
    }
    return res;
}

size_t
limdd_node_table_size(void)
{
    return limdd_nodes == NULL ? 0 : llmsset_get_size(limdd_nodes);
}

size_t
limdd_node_table_count(void)
{
    if (limdd_nodes == NULL) return 0;
    /* buckets 0 and 1 are reserved at creation; 1 is the terminal */
    const size_t occupied = (size_t)llmsset_count_marked(limdd_nodes);
    return occupied < 2 ? 0 : occupied - 2;
}

/**
 * Amplitude of basis state `b` under edge `e`, where `e` lives at `level` and
 * bit k % 64 of word k / 64 of `b` is qubit k.
 *
 * Applying the edge's LIM `w * P` to the vector v that its target denotes:
 * P maps the basis state |c> to a phase times |c XOR x>, so the only c that
 * contributes to <b| is c = b XOR x. Per qubit the phase is 1 for I and X,
 * (-1)^c_k for Z, and i*(-1)^c_k for Y, giving
 *
 *     amp(e, b) = w * i^(#Y) * (-1)^popcount(z & c) * amp(target, c)
 *
 * and the two phase factors combine into the single power of i below.
 */
static EVBDD_WGT
eval_edge(LIMDD e, const uint64_t b[LIMDD_PAULI_WORDS], uint32_t level)
{
    if (limdd_edge_is_zero(e)) return EVBDD_ZERO;

    const LIMDD_LIM lim = limdd_label(e);
    const limdd_pauli_t p = limdd_lim_pauli(lim);

    /*
     * An edge at `level` spans qubits level..n-1, so its LIM may not act on
     * anything above it. A violation means a label was attached at the wrong
     * depth, which would otherwise show up only as a wrong amplitude much
     * later.
     */
    assert(!limdd_pauli_acts_below(p, level));

    /* i per Y, and -1 == i^2 per Z or Y meeting a set bit. This also advances
     * `b` to the basis state that actually reaches the target. */
    uint64_t c[LIMDD_PAULI_WORDS];
    for (unsigned i = 0; i < LIMDD_PAULI_WORDS; i++) c[i] = b[i];
    const unsigned k = limdd_pauli_apply_basis(p, c);

    EVBDD_WGT w = limdd_lim_weight(lim);
    if (k != 0) w = wgt_mul(w, limdd_wgt_i_pow(k));

    const LIMDD_TARG t = limdd_target(e);
    const uint32_t lev = limdd_level(t);
    assert(lev >= level && "an edge may skip levels, never climb them");
    (void)level;   /* only the assertions consult it */
    if (t == LIMDD_TERMINAL) return w;

    /*
     * Levels level..lev-1 are skipped. Each holds |0>+|1>, whose amplitude is
     * 1 on either branch, so beyond the phase already folded into k above --
     * (-1)^c_j for a Z there, -i(-1)^c_j for a Y, nothing for I or X -- they
     * contribute nothing, and the walk resumes at the target's own level.
     */
    const LIMDD child = (c[LIMDD_PAULI_LANE(lev)] & LIMDD_PAULI_BIT(lev))
                      ? limdd_node_high(t) : limdd_node_low(t);
    const EVBDD_WGT sub = eval_edge(child, c, lev + 1);
    if (sub == EVBDD_ZERO) return EVBDD_ZERO;

    return wgt_mul(w, sub);
}

EVBDD_WGT
limdd_eval(LIMDD e, const bool *bits, size_t nqubits)
{
    assert(nqubits == limdd_nqubits);

    uint64_t b[LIMDD_PAULI_WORDS] = { 0 };
    for (size_t k = 0; k < nqubits; k++) {
        if (bits[k]) b[LIMDD_PAULI_LANE(k)] |= LIMDD_PAULI_BIT(k);
    }
    return eval_edge(e, b, 0);
}

int
limdd_gc_mark_node(LIMDD_TARG p)
{
    assert(limdd_nodes != NULL);
    if (p == LIMDD_TERMINAL) return 0;   /* reserved, never swept */
    return llmsset_mark(limdd_nodes, p);
}

void
limdd_gc_clear_nodes(void)
{
    llmsset_clear_data(limdd_nodes);
}

void
limdd_gc_rehash_nodes(void)
{
    llmsset_clear_hashes(limdd_nodes);
    if (llmsset_rehash(limdd_nodes) != 0) {
        fprintf(stderr, "sylvan: LIMDD nodes could not all be rehashed\n");
        exit(1);
    }
}

void
limdd_gc_purge_stab_cache(void)
{
    for (uint64_t b = 2; b < limdd_nodes->table_size; b++) {
        const uint64_t word = limdd_nodes->bitmap2[b / 64];
        const uint64_t mask = UINT64_C(0x8000000000000000) >> (b & 63);
        if (!(word & mask)) atomic_store_explicit(&node_stab[b], 0, memory_order_relaxed);
    }
}

void
limdd_nodes_init(size_t nqubits, size_t node_tablesize,
                 size_t pauli_tablesize, size_t lim_tablesize,
                 size_t stab_tablesize)
{
    limdd_lims_init(nqubits, pauli_tablesize, lim_tablesize);
    limdd_stab_init(stab_tablesize);

    limdd_nqubits = nqubits;
    limdd_nodes = llmsset_create(node_tablesize, node_tablesize);
    if (limdd_nodes == NULL) {
        fprintf(stderr, "sylvan: could not create the LIMDD node table\n");
        exit(1);
    }

    node_stab = calloc(node_tablesize, sizeof(_Atomic(uint64_t)));
    if (node_stab == NULL) {
        fprintf(stderr, "sylvan: could not allocate the LIMDD stabiliser cache\n");
        exit(1);
    }

    /*
     * The edge weight table is collected by copying: survivors move to a new
     * table and every EVBDD is rewritten with their new indices. A LIMDD
     * cannot be rewritten that way -- a LIM is interned under (Pauli, weight),
     * so new weight indices mean new LIM indices, new node contents and new
     * node indices, i.e. the whole forest. Left on, an EVBDD operation
     * anywhere in the program would silently invalidate every live LIMDD.
     * See qsylvan_limdd_gc.h.
     */
    evbdd_set_auto_gc_wgt_table(false);

    /* One reference stack per worker, so operations can be interrupted by a
     * collection without losing what they are holding. Declared here rather
     * than included: qsylvan_limdd_gc.h reaches this header through canon.h,
     * so including it back would be circular. */
    {
        extern void limdd_refs_init(void);
        limdd_refs_init();
    }

    /* The operation cache may hold an earlier session's entries, keyed on the
     * indices these new tables hand out again (qsylvan_limdd_gc.c). */
    {
        extern void limdd_gc_new_session(void);
        limdd_gc_new_session();
    }

    zero_edge = limdd_bundle(LIMDD_LIM_ZERO, LIMDD_TERMINAL);
    one_edge  = limdd_bundle(LIMDD_LIM_IDENTITY, LIMDD_TERMINAL);
}

void
limdd_nodes_quit(void)
{
    {
        extern void limdd_refs_quit(void);     /* declared here for the reason above */
        limdd_refs_quit();
    }
    if (limdd_nodes != NULL) {
        llmsset_free(limdd_nodes);
        limdd_nodes = NULL;
    }
    free(node_stab);
    node_stab = NULL;
    limdd_nqubits = 0;
    zero_edge = 0;
    one_edge = 0;
    limdd_stab_quit();
    limdd_lims_quit();
}
