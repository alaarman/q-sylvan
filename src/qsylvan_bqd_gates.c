/*
 * Copyright 2025 Q-Sylvan contributors
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

#include "qsylvan_bqd_gates.h"
#include "qsylvan_bqd_exp.h"
#include "qsylvan_bqd_sm.h"
#include "qsylvan_bqd_xp.h"
#include "qsylvan_gates.h"
#include "qsylvan_limdd_gc.h"

/* --- edges ---------------------------------------------------------------- */

static inline EVBDD_WGT
scalar_of(BQD e)
{
    return limdd_lim_weight(limdd_label(e));
}

/**
 * c . [t], the zero edge for c = 0. The label is the one bqd_build writes,
 * interned through the per-worker memo of bqd_lim_word.
 */
static inline BQD
edge(EVBDD_WGT c, LIMDD_TARG t)
{
    if (c == EVBDD_ZERO) return limdd_zero_edge();
    return limdd_bundle(bqd_lim_word(c, 0, 0), t);
}

static inline BQD
unit(LIMDD_TARG t)
{
    return limdd_bundle(LIMDD_LIM_IDENTITY, t);
}

static inline BQD
scale(EVBDD_WGT c, BQD e)
{
    if (c == EVBDD_ZERO || limdd_edge_is_zero(e)) return limdd_zero_edge();
    if (c == EVBDD_ONE) return e;
    return edge(wgt_mul(c, scalar_of(e)), limdd_target(e));
}

/* --- pointwise operations -------------------------------------------------- */

/*
 * Four operations, the two a state needs and the two the quotient rule is
 * made of. xprod(a, b) is a.b where a is nonzero and b where it is zero, the
 * (.) of def:bqd, and xquot(b, a) its inverse, b/a where a is nonzero and b
 * where it is zero. A node's high cofactor is xprod(low, high), and the high
 * child of the node of (f_0, f_1) is xquot(f_1, f_0), both up to the scale
 * the representative divides out.
 */
enum { OP_MUL = 0, OP_ADD = 1, OP_XPROD = 2, OP_XQUOT = 3 };

static uint64_t
op_cache_id(int op)
{
    switch (op) {
    case OP_MUL:   return CACHE_BQD_MUL;
    case OP_ADD:   return CACHE_BQD_ADD;
    case OP_XPROD: return CACHE_BQD_XPROD;
    default:       return CACHE_BQD_XQUOT;
    }
}

static EVBDD_WGT
op_scalar(int op, EVBDD_WGT a, EVBDD_WGT b)
{
    switch (op) {
    case OP_MUL:   return wgt_mul(a, b);
    case OP_ADD:   return wgt_add(a, b);
    case OP_XPROD: return (a == EVBDD_ZERO) ? b : wgt_mul(a, b);
    default:       return (b == EVBDD_ZERO) ? a : wgt_div(a, b);
    }
}

TASK_DECL_3(BQD, bqd_apply_op, int, BQD, BQD);
TASK_DECL_1(BQD, bqd_high_cofactor, LIMDD_TARG);
TASK_DECL_3(BQD, bqd_compose_scalar, uint32_t, BQD, BQD);
TASK_DECL_3(BQD, bqd_cofactor_scalar, BQD, uint32_t, int);
TASK_DECL_3(BQD, bqd_restrict_scalar, BQD, uint32_t, int);
TASK_DECL_3(BQD, bqd_project_scalar, BQD, uint32_t, int);
TASK_DECL_3(BQD, bqd_rebuild, uint32_t, BQD, BQD);

/**
 * The function with cofactors (a, b), at level var, in canonical form.
 *
 * Equal cofactors are a function that does not depend on x_var, and its
 * canonical edge is theirs, read at var: the level is skipped (Compose of
 * skip:alg:constructors). Both are canonical, so equal functions are equal
 * edges and that is the whole test, O(1). With a zero, the node is (0, b)
 * and its representative is b's, so the node stores b's node on a high edge
 * with the identity and the scale goes on top. Otherwise the node is f / c_a,
 * whose low half is a's node, a representative already, and whose ratio is
 * xquot(b / c_a, a's node): the scale has to come off b before the quotient
 * and not after, because where a is zero the quotient copies b and does not
 * divide it.
 */
TASK_IMPL_3(BQD, bqd_compose_scalar, uint32_t, var, BQD, a, BQD, b)
{
    if (a == b) return a;
    if (limdd_edge_is_zero(a)) {
        if (limdd_edge_is_zero(b)) return limdd_zero_edge();
        const LIMDD_TARG node = limdd_makenode(var, limdd_zero_edge(), unit(limdd_target(b)));
        return edge(scalar_of(b), node);
    }
    const EVBDD_WGT ca = scalar_of(a);
    const BQD b1 = scale(wgt_div(EVBDD_ONE, ca), b);
    const BQD r = CALL(bqd_apply_op, OP_XQUOT, b1, unit(limdd_target(a)));
    const LIMDD_TARG node = limdd_makenode(var, unit(limdd_target(a)), r);
    return edge(ca, node);
}

/**
 * The high cofactor of a node's function: xprod(low, high), both read at the
 * level below the node, which either may skip. Memoised.
 */
TASK_IMPL_1(BQD, bqd_high_cofactor, LIMDD_TARG, t)
{
    const BQD lo = limdd_node_low(t), hi = limdd_node_high(t);
    if (limdd_edge_is_zero(lo)) return hi;                  /* the copy everywhere */
    if (limdd_edge_is_zero(hi)) return limdd_zero_edge();
    uint64_t hit;
    if (cache_get3(CACHE_BQD_COF1, t, 0, 0, &hit)) return (BQD)hit;
    bqd_count(LACE_WORKER_ID, BQD_COUNT_COF1);
    const BQD r = CALL(bqd_apply_op, OP_XPROD, lo, hi);
    cache_put3(CACHE_BQD_COF1, t, 0, 0, (uint64_t)r);
    return r;
}

/**
 * op(f, g) for two edges read at one level. The scalars are taken out first,
 * as far as the operation allows, so the memo is on the two nodes and what is
 * left of the scalars:
 *
 *     mul    c.F  d.G  =  cd . (F G)
 *     add    c.F  d.G  =  c . (F + (d/c) G)
 *     xprod  c.F  d.G  =  d . xprod(c F, G)        (the copy keeps c)
 *     xquot  c.F  d.G  =  c . xquot(F, d G)
 *
 * The recursion is at the higher of the two nodes' levels, the smaller
 * variable. An operand whose node is lower, or is the terminal, skips that
 * level, and is then its own cofactor on both sides (skip:lem:virtual), so
 * the two reach the terminal together only when both are there. The key
 * holds no level: an edge denotes the same function at every level it is
 * read at, extended to the levels it skips (skip:lem:ext).
 */
TASK_IMPL_3(BQD, bqd_apply_op, int, op, BQD, f, BQD, g)
{
    const bool fz = limdd_edge_is_zero(f), gz = limdd_edge_is_zero(g);
    switch (op) {
    case OP_MUL:   if (fz || gz) return limdd_zero_edge(); break;
    case OP_ADD:   if (fz) return g; if (gz) return f; break;
    case OP_XPROD: if (fz) return g; if (gz) return limdd_zero_edge(); break;
    default:       if (gz) return f; if (fz) return limdd_zero_edge(); break;
    }

    const LIMDD_TARG tf = limdd_target(f), tg = limdd_target(g);
    const EVBDD_WGT cf = scalar_of(f), cg = scalar_of(g);
    if (tf == LIMDD_TERMINAL && tg == LIMDD_TERMINAL)
        return edge(op_scalar(op, cf, cg), LIMDD_TERMINAL);

    EVBDD_WGT outer, k;
    BQD F, G;
    switch (op) {
    case OP_MUL:   outer = wgt_mul(cf, cg); k = EVBDD_ONE;       F = unit(tf);    G = unit(tg);    break;
    case OP_ADD:   outer = cf;  k = wgt_div(cg, cf);             F = unit(tf);    G = edge(k, tg); break;
    case OP_XPROD: outer = cg;  k = cf;                          F = edge(k, tf); G = unit(tg);    break;
    default:       outer = cf;  k = cg;                          F = unit(tf);    G = edge(k, tg); break;
    }

    uint64_t hit;
    if (cache_get3(op_cache_id(op), tf, tg, k, &hit)) return scale(outer, (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_APPLY);

    const uint32_t lf = limdd_level(tf), lg = limdd_level(tg);
    const uint32_t var = lf < lg ? lf : lg;
    const BQD f0 = CALL(bqd_cofactor_scalar, F, var, 0);
    const BQD g0 = CALL(bqd_cofactor_scalar, G, var, 0);
    const BQD f1 = CALL(bqd_cofactor_scalar, F, var, 1);
    const BQD g1 = CALL(bqd_cofactor_scalar, G, var, 1);

    /* the two cofactors of the result are independent */
    limdd_refs_spawn(SPAWN(bqd_apply_op, op, f0, g0));
    const BQD hi = limdd_refs_push(CALL(bqd_apply_op, op, f1, g1));
    const BQD lo = limdd_refs_sync(SYNC(bqd_apply_op));
    limdd_refs_pop(1);

    const BQD r = CALL(bqd_compose_scalar, var, lo, hi);
    cache_put3(op_cache_id(op), tf, tg, k, (uint64_t)r);
    return scale(outer, r);
}

/*
 * The entry points below are where the translation and Pauli families part:
 * their labels are not scalars, so they go to qsylvan_bqd_xp.h. The scalar
 * recursion calls the scalar tasks directly and never tests the family. A
 * session of rule SM, which is the scalar family's, goes to qsylvan_bqd_sm.h
 * before either, so the copy rule's code never meets an SM node.
 */
TASK_IMPL_2(BQD, bqd_multiply, BQD, f, BQD, g)
{
    if (bqd_sm()) return CALL(bqd_sm_apply, BQD_SM_MUL, f, g);
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_apply_op, BQD_OP_MUL, f, g);
    return CALL(bqd_apply_op, OP_MUL, f, g);
}

TASK_IMPL_2(BQD, bqd_add, BQD, f, BQD, g)
{
    if (bqd_sm()) return CALL(bqd_sm_apply, BQD_SM_ADD, f, g);
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_apply_op, BQD_OP_ADD, f, g);
    return CALL(bqd_apply_op, OP_ADD, f, g);
}

TASK_IMPL_3(BQD, bqd_compose, uint32_t, var, BQD, a, BQD, b)
{
    if (bqd_sm()) return CALL(bqd_sm_compose, var, a, b);
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_compose, var, a, b);
    return CALL(bqd_compose_scalar, var, a, b);
}

TASK_IMPL_3(BQD, bqd_cofactor_scalar, BQD, e, uint32_t, var, int, b)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG t = limdd_target(e);
    assert(var <= limdd_level(t) && "an edge is read at or above its node");
    /* an edge that skips var denotes a function that does not depend on
     * x_var, so it is both of its cofactors there (skip:lem:virtual) */
    if (limdd_level(t) > var) return e;
    const BQD c = b ? CALL(bqd_high_cofactor, t) : limdd_node_low(t);
    return scale(scalar_of(e), c);
}

/* In the other two families the cofactor is a label product, which Canon
 * makes the canonical edge (skip:lem:xcof, skip:prop:xcanon). */
TASK_IMPL_3(BQD, bqd_cofactor, BQD, e, uint32_t, var, int, b)
{
    if (bqd_sm()) return CALL(bqd_sm_cof, e, var, b);
    if (bqd_family() == BQD_FAMILY_SCALAR) return CALL(bqd_cofactor_scalar, e, var, b);
    const BQD c = limdd_refs_push(CALL(bqd_xp_cofactor, e, var, b));
    const BQD r = CALL(bqd_xp_canon, c);
    limdd_refs_pop(1);
    return r;
}

BQD
bqd_scale(BQD e, EVBDD_WGT c)
{
    if (bqd_family() != BQD_FAMILY_SCALAR) return bqd_xp_scale(e, c);
    return scale(c, e);
}

BQD
bqd_negate(BQD e)
{
    return bqd_scale(e, EVBDD_MIN_ONE);
}

/*
 * Restriction and projection, both linear, so memoised on the node. At q the
 * restriction is the chosen cofactor itself: it does not depend on x_q, so as
 * an edge read at q it skips q and needs no node there. The projection
 * composes the chosen cofactor with zero. An edge whose node is below q, or is
 * the terminal, skips q and denotes a function that does not depend on x_q
 * (skip:alg:gates): it is its own restriction, and its projection is the node
 * at q with the edge on one side and zero on the other. Both are answered
 * before the lookup, in O(1).
 *
 * Above q both take each value of their result from one value of their
 * operand, so they commute with the pair of a node's low cofactor and ratio
 * (skip:lem:select): the recursion goes down the node's two stored edges and
 * Rebuild makes the node of the two results, RestrictR and ProjectR of
 * skip:alg:perm. No high cofactor is made above q, and none but the one at q.
 */
TASK_IMPL_3(BQD, bqd_restrict_scalar, BQD, e, uint32_t, q, int, b)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG t = limdd_target(e);
    const uint32_t var = limdd_level(t);
    if (var > q) return e;
    const uint64_t key = ((uint64_t)q << 1) | (uint64_t)(b & 1);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_RESTRICT, t, key, 0, &hit)) return scale(scalar_of(e), (BQD)hit);

    const BQD f0 = limdd_node_low(t);
    BQD r;
    if (var == q) {
        r = b ? CALL(bqd_high_cofactor, t) : f0;
    } else {
        limdd_refs_spawn(SPAWN(bqd_restrict_scalar, limdd_node_high(t), q, b));
        const BQD lo = limdd_refs_push(CALL(bqd_restrict_scalar, f0, q, b));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_restrict_scalar)));
        r = CALL(bqd_rebuild, var, lo, hi);
        limdd_refs_pop(2);
    }
    cache_put3(CACHE_BQD_RESTRICT, t, key, 0, (uint64_t)r);
    return scale(scalar_of(e), r);
}

TASK_IMPL_3(BQD, bqd_project_scalar, BQD, e, uint32_t, q, int, b)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG t = limdd_target(e);
    const uint32_t var = limdd_level(t);
    if (var > q) return b ? CALL(bqd_compose_scalar, q, limdd_zero_edge(), e)
                          : CALL(bqd_compose_scalar, q, e, limdd_zero_edge());
    const uint64_t key = ((uint64_t)q << 1) | (uint64_t)(b & 1);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_PROJECT, t, key, 0, &hit)) return scale(scalar_of(e), (BQD)hit);

    const BQD f0 = limdd_node_low(t);
    BQD r;
    if (var == q) {
        r = b ? CALL(bqd_compose_scalar, var, limdd_zero_edge(), CALL(bqd_high_cofactor, t))
              : CALL(bqd_compose_scalar, var, f0, limdd_zero_edge());
    } else {
        limdd_refs_spawn(SPAWN(bqd_project_scalar, limdd_node_high(t), q, b));
        const BQD lo = limdd_refs_push(CALL(bqd_project_scalar, f0, q, b));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_project_scalar)));
        r = CALL(bqd_rebuild, var, lo, hi);
        limdd_refs_pop(2);
    }
    cache_put3(CACHE_BQD_PROJECT, t, key, 0, (uint64_t)r);
    return scale(scalar_of(e), r);
}

/* In the other two families a cofactor is a labelled edge, so the chosen one
 * goes through Canon where the scalar family returns it (skip:alg:xrestrict). */
TASK_IMPL_3(BQD, bqd_restrict, BQD, e, uint32_t, q, int, b)
{
    if (bqd_sm()) return CALL(bqd_sm_restrict, e, q, b);
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_restrict, e, q, b);
    return CALL(bqd_restrict_scalar, e, q, b);
}

TASK_IMPL_3(BQD, bqd_project, BQD, e, uint32_t, q, int, b)
{
    if (bqd_sm()) return CALL(bqd_sm_project, e, q, b);
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_project, e, q, b);
    return CALL(bqd_project_scalar, e, q, b);
}

BQD
bqd_local_matvec(BQD e, const EVBDD_WGT *M, const uint32_t *qubits, uint32_t k, uint32_t n)
{
    (void)n;
    const uint32_t dim = 1u << k;
    BQD *vc = malloc(dim * sizeof(BQD));
    if (vc == NULL) { fprintf(stderr, "sylvan: out of memory in bqd_local_matvec\n"); exit(1); }
    for (uint32_t c = 0; c < dim; c++) {
        BQD r = e;
        for (uint32_t i = 0; i < k; i++) r = bqd_restrict(r, qubits[i], (int)((c >> (k - 1 - i)) & 1));
        vc[c] = limdd_refs_push(r);
    }
    BQD out = limdd_refs_push(limdd_zero_edge());
    for (uint32_t r = 0; r < dim; r++) {
        BQD w = limdd_zero_edge();
        for (uint32_t c = 0; c < dim; c++) {
            const EVBDD_WGT m = M[(size_t)r * dim + c];
            if (m == EVBDD_ZERO) continue;
            limdd_refs_push(w);
            w = bqd_add(w, bqd_scale(vc[c], m));
            limdd_refs_pop(1);
        }
        for (uint32_t i = 0; i < k; i++) {
            limdd_refs_push(w);
            w = bqd_project(w, qubits[i], (int)((r >> (k - 1 - i)) & 1));
            limdd_refs_pop(1);
        }
        limdd_refs_push(w);
        const BQD next = bqd_add(out, w);
        limdd_refs_pop(2);
        out = limdd_refs_push(next);
    }
    limdd_refs_pop(1 + (long)dim);
    free(vc);
    return out;
}

/* --- states --------------------------------------------------------------- */

/* In the other two families every basis state is one chain under a label. */
BQD
bqd_basis_state(uint64_t x, uint32_t n)
{
    if (bqd_sm()) return bqd_sm_basis_state(x, n);
    if (bqd_family() != BQD_FAMILY_SCALAR) return bqd_xp_basis_state(x, n);
    LIMDD_TARG t = LIMDD_TERMINAL;
    for (uint32_t v = n; v-- > 0; ) {
        const bool one = (x >> (n - 1 - v)) & 1;
        t = one ? limdd_makenode(v, limdd_zero_edge(), unit(t))
                : limdd_makenode(v, unit(t), limdd_zero_edge());
    }
    return unit(t);
}

/*
 * A zero edge anywhere means a zero amplitude: a zero low edge makes f_0 zero,
 * and a zero high edge makes f_1 = f_0 (.) 0 zero. Without one, every low
 * cofactor and every ratio is nonzero everywhere, and so is f_1 = f_0 r.
 */
TASK_1(int, bqd_full_rec, LIMDD_TARG, t)
{
    if (t == LIMDD_TERMINAL) return 1;
    const BQD lo = limdd_node_low(t), hi = limdd_node_high(t);
    if (limdd_edge_is_zero(lo) || limdd_edge_is_zero(hi)) return 0;
    uint64_t hit;
    if (cache_get3(CACHE_BQD_FULL, t, 0, 0, &hit)) return (int)hit - 1;
    SPAWN(bqd_full_rec, limdd_target(lo));
    const int h = CALL(bqd_full_rec, limdd_target(hi));
    const int l = SYNC(bqd_full_rec);
    const int r = h && l;
    cache_put3(CACHE_BQD_FULL, t, 0, 0, (uint64_t)r + 1);   /* 0 would read as a miss */
    return r;
}

bool
bqd_has_full_support(BQD e)
{
    if (limdd_edge_is_zero(e)) return false;
    if (bqd_sm()) return bqd_sm_full(limdd_target(e));          /* the node's flag, O(1) */
    return RUN(bqd_full_rec, limdd_target(e)) != 0;
}

/* --- gates without the high cofactor (skip:sec:ratio) ---------------------- */

/*
 * A node stores its low cofactor and the ratio of the two, not its high
 * cofactor, so a recursion on cofactors makes one (Cof1, a product) at every
 * node above the level it acts on and forms each new ratio by a quotient in
 * Compose. On the circuits that stalled that was 100 to 1000 times the work of
 * the gate's own recursion. An operation that takes each value of its result
 * from one value of an operand commutes with the pair of a node's low
 * cofactor and ratio (skip:lem:select): the permutations of two qubits, X,
 * restriction and projection. Those recurse on a node's two stored edges above
 * the qubits they act on, and Rebuild makes the node of the two results, with
 * no Cof1 and no Apply there (skip:prop:perm). A diagonal gate is not one, but
 * it multiplies by beta^eps for an exponent diagram eps (qsylvan_bqd_exp.h),
 * and PhaseMul keeps it at the stored edges on any support
 * (skip:prop:phasemul). Only the level a gate acts on, and Pair's level b,
 * still make a high cofactor and a quotient.
 *
 * A node below is written 0 where a function is zero: 0 is no bucket llmsset
 * hands out.
 */

/** The node of an edge, or 0 for the zero edge. */
static inline LIMDD_TARG
node_or_zero(BQD e)
{
    return limdd_edge_is_zero(e) ? 0 : limdd_target(e);
}

/**
 * Ind of skip:alg:constructors, the node of the indicator of supp [t], by the
 * support shadow (skip:lem:ind): the support of a node is {0} x supp [N_0]
 * and {1} x supp [N_1], so the indicator's children are the indicators of
 * t's, with the identity, and a node whose two are one node skips its level.
 * The terminal on full support. Memoised on t.
 */
TASK_1(LIMDD_TARG, bqd_ind_rec, LIMDD_TARG, t)
{
    if (t == LIMDD_TERMINAL) return t;
    uint64_t hit;
    if (cache_get3(CACHE_BQD_IND, t, 0, 0, &hit)) return (LIMDD_TARG)hit;
    const LIMDD_TARG n0 = node_or_zero(limdd_node_low(t)), n1 = node_or_zero(limdd_node_high(t));
    if (n0 != 0) SPAWN(bqd_ind_rec, n0);
    const LIMDD_TARG i1 = (n1 != 0) ? CALL(bqd_ind_rec, n1) : 0;
    const LIMDD_TARG i0 = (n0 != 0) ? SYNC(bqd_ind_rec) : 0;
    const LIMDD_TARG r = (i0 != 0 && i0 == i1)
        ? i0
        : limdd_makenode(limdd_node_var(t), i0 ? unit(i0) : limdd_zero_edge(),
                         i1 ? unit(i1) : limdd_zero_edge());
    cache_put3(CACHE_BQD_IND, t, 0, 0, r);
    return r;
}

/**
 * iota(t) of skip:alg:phasemul, the indicator of supp [t] as a 0/1 exponent:
 * the support shadow again, as an exponent node. Memoised on t; t is 0 for
 * the zero function.
 */
TASK_1(BQD_EXP, bqd_iota_rec, LIMDD_TARG, t)
{
    if (t == 0) return mtbdd_int64(0);
    if (t == LIMDD_TERMINAL) return mtbdd_int64(1);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_IOTA, t, 0, 0, &hit)) return (BQD_EXP)hit;
    mtbdd_refs_spawn(SPAWN(bqd_iota_rec, node_or_zero(limdd_node_low(t))));
    const BQD_EXP i1 = mtbdd_refs_push(CALL(bqd_iota_rec, node_or_zero(limdd_node_high(t))));
    const BQD_EXP i0 = mtbdd_refs_sync(SYNC(bqd_iota_rec));
    mtbdd_refs_pop(1);
    const BQD_EXP r = bqd_exp_node(limdd_node_var(t), i0, i1);
    cache_put3(CACHE_BQD_IOTA, t, 0, 0, r);
    return r;
}

/**
 * Whether supp [a] lies inside supp [u], on the two support shadows alone:
 * the optional test of MulOff (skip:alg:rebuild), which needs no arithmetic.
 * Memoised on the pair; 0 is the zero function.
 */
TASK_2(int, bqd_supp_subset, LIMDD_TARG, a, LIMDD_TARG, u)
{
    if (a == 0 || a == u) return 1;
    if (u == 0) return 0;
    if (CALL(bqd_full_rec, u)) return 1;
    uint64_t hit;
    if (cache_get3(CACHE_BQD_SUBSET, a, u, 0, &hit)) return (int)hit - 1;
    const uint32_t la = limdd_level(a), lu = limdd_level(u);
    const uint32_t v = la < lu ? la : lu;
    /* a node below v is both of its cofactors there */
    const LIMDD_TARG a0 = la > v ? a : node_or_zero(limdd_node_low(a));
    const LIMDD_TARG a1 = la > v ? a : node_or_zero(limdd_node_high(a));
    const LIMDD_TARG u0 = lu > v ? u : node_or_zero(limdd_node_low(u));
    const LIMDD_TARG u1 = lu > v ? u : node_or_zero(limdd_node_high(u));
    SPAWN(bqd_supp_subset, a0, u0);
    const int s1 = CALL(bqd_supp_subset, a1, u1);
    const int s0 = SYNC(bqd_supp_subset);
    const int r = s0 && s1;
    cache_put3(CACHE_BQD_SUBSET, a, u, 0, (uint64_t)r + 1);   /* 0 would read as a miss */
    return r;
}

TASK_DECL_4(BQD, bqd_phase_mul_rec, BQD, BQD_EXP, EVBDD_WGT, uint32_t);

/**
 * MulOff of skip:alg:rebuild: the canonical edge of [k] . c^{[x not in supp [u]]},
 * which is k on the support of u and c k off it. It is k when c is 1, when u
 * has full support and when supp [k] is inside supp [u]; otherwise it is the
 * phase multiplication by c^{1 - iota(u)}.
 */
TASK_3(BQD, bqd_mul_off_rec, BQD, k, EVBDD_WGT, c, LIMDD_TARG, u)
{
    if (limdd_edge_is_zero(k) || c == EVBDD_ONE) return k;
    if (u == 0) return scale(c, k);
    if (CALL(bqd_full_rec, u) || CALL(bqd_supp_subset, limdd_target(k), u)) return k;
    const uint32_t r = bqd_exp_order(c);
    const BQD_EXP one = mtbdd_refs_push(bqd_exp_const(1, r));
    const BQD_EXP io = mtbdd_refs_push(CALL(bqd_iota_rec, u));
    const BQD_EXP eps = mtbdd_refs_push(CALL(bqd_exp_sub, one, io, r));
    const BQD res = CALL(bqd_phase_mul_rec, k, eps, c, r);
    mtbdd_refs_pop(3);
    return res;
}

/**
 * Split of skip:alg:rebuild: the canonical edges of the low cofactor and of
 * the ratio at var of the function e denotes read at var (skip:lem:ratio).
 * Where e skips var both cofactors are e, and the ratio is the indicator of
 * its support; a zero low edge makes the ratio c times the high edge
 * everywhere; otherwise the ratio is the high edge on the support of the low
 * child and c times it off that support, which MulOff makes, and which is the
 * high edge itself for c = 1.
 */
VOID_TASK_4(bqd_split, BQD, e, uint32_t, var, BQD *, lo, BQD *, ratio)
{
    if (limdd_edge_is_zero(e)) { *lo = *ratio = e; return; }
    const LIMDD_TARG t = limdd_target(e);
    if (limdd_level(t) > var) {
        *lo = e;
        *ratio = unit(CALL(bqd_ind_rec, t));
        return;
    }
    const EVBDD_WGT c = scalar_of(e);
    const BQD low = limdd_node_low(t), high = limdd_node_high(t);
    if (limdd_edge_is_zero(low)) {
        *lo = low;
        *ratio = scale(c, high);
        return;
    }
    *lo = scale(c, low);
    *ratio = CALL(bqd_mul_off_rec, high, c, limdd_target(low));
}

/**
 * Rebuild of skip:alg:rebuild: the canonical edge at var of the function g
 * with the low cofactor lo = E(g_0) and the ratio r = E(rho(g)). A zero low
 * cofactor is Compose's case. Otherwise g_0 = a [U], and the node of g / a has
 * the low child U and, as its high edge, the ratio on supp [U] and r / a off
 * it, MulOff(r, 1/a, U); the level is skipped where that is the indicator of
 * supp [U], the terminal on full support (skip:prop:rebuild).
 */
TASK_IMPL_3(BQD, bqd_rebuild, uint32_t, var, BQD, lo, BQD, r)
{
    if (limdd_edge_is_zero(lo)) return CALL(bqd_compose_scalar, var, lo, r);
    const EVBDD_WGT a = scalar_of(lo);
    const LIMDD_TARG U = limdd_target(lo);
    if (limdd_edge_is_zero(r)) return edge(a, limdd_makenode(var, unit(U), r));
    const BQD k = (a == EVBDD_ONE) ? r : CALL(bqd_mul_off_rec, r, wgt_div(EVBDD_ONE, a), U);
    if (limdd_label(k) == LIMDD_LIM_IDENTITY && limdd_target(k) == CALL(bqd_ind_rec, U)) return lo;
    return edge(a, limdd_makenode(var, unit(U), k));
}

/**
 * XQ of skip:alg:perm, X on qubit q. Memoised on the node and q, the scalar
 * outside. Above q it is a selection, so the node's two stored edges are
 * treated and rebuilt; at q the two cofactors trade places, which needs the
 * high one and Compose's quotient. An edge that skips q is its own image.
 */
TASK_2(BQD, bqd_xq_rec, BQD, e, uint32_t, q)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG t = limdd_target(e);
    const uint32_t var = limdd_level(t);
    if (var > q) return e;
    uint64_t hit;
    if (cache_get3(CACHE_BQD_XQ, t, q, 0, &hit)) return scale(scalar_of(e), (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    BQD r;
    if (var == q) {
        const BQD f1 = limdd_refs_push(CALL(bqd_high_cofactor, t));
        r = CALL(bqd_compose_scalar, q, f1, limdd_node_low(t));
        limdd_refs_pop(1);
    } else {
        limdd_refs_spawn(SPAWN(bqd_xq_rec, limdd_node_high(t), q));
        const BQD lo = limdd_refs_push(CALL(bqd_xq_rec, limdd_node_low(t), q));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_xq_rec)));
        r = CALL(bqd_rebuild, var, lo, hi);
        limdd_refs_pop(2);
    }
    cache_put3(CACHE_BQD_XQ, t, q, 0, (uint64_t)r);
    return scale(scalar_of(e), r);
}

/**
 * Pair of skip:alg:perm, [x_b = 0] A|_{x_b = s} + [x_b = 1] B|_{x_b = t}: the
 * function that takes its values from A where x_b is 0 and from B where it is
 * 1, each restricted. Linear in the pair, so the scale g of the first nonzero
 * operand comes out, and the memo is on the two nodes, B's scale over g, s, t
 * and b. Above b it is a selection of the two operands, whose low cofactors
 * and ratios Split gives at the higher of their levels, and Rebuild puts back
 * together; at b it composes the two chosen cofactors, and below b the two
 * operands themselves.
 */
TASK_5(BQD, bqd_pair_rec, BQD, A, int, s, BQD, B, int, t, uint32_t, b)
{
    const bool az = limdd_edge_is_zero(A), bz = limdd_edge_is_zero(B);
    if (az && bz) return A;
    const EVBDD_WGT g = az ? scalar_of(B) : scalar_of(A);
    const BQD A1 = az ? A : unit(limdd_target(A));
    const BQD B1 = bz ? B : az ? unit(limdd_target(B)) : scale(wgt_div(EVBDD_ONE, g), B);
    const uint64_t kappa = bz ? 0 : (uint64_t)scalar_of(B1);
    const uint64_t key = ((uint64_t)node_or_zero(B1) << 16) | ((uint64_t)b << 2)
                       | ((uint64_t)(s & 1) << 1) | (uint64_t)(t & 1);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_PAIR, node_or_zero(A1), key, kappa, &hit)) return scale(g, (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    const uint32_t la = limdd_level(limdd_target(A1)), lb = limdd_level(limdd_target(B1));
    const uint32_t l = la < lb ? la : lb;
    BQD r;
    if (l > b) {
        r = CALL(bqd_compose_scalar, b, A1, B1);
    } else if (l == b) {
        const BQD a0 = limdd_refs_push(CALL(bqd_cofactor_scalar, A1, b, s));
        const BQD b0 = limdd_refs_push(CALL(bqd_cofactor_scalar, B1, b, t));
        r = CALL(bqd_compose_scalar, b, a0, b0);
        limdd_refs_pop(2);
    } else {
        BQD A0, Ar, B0, Br;
        CALL(bqd_split, A1, l, &A0, &Ar);
        limdd_refs_push(A0);
        limdd_refs_push(Ar);
        CALL(bqd_split, B1, l, &B0, &Br);
        limdd_refs_push(B0);
        limdd_refs_push(Br);
        limdd_refs_spawn(SPAWN(bqd_pair_rec, Ar, s, Br, t, b));
        const BQD lo = limdd_refs_push(CALL(bqd_pair_rec, A0, s, B0, t, b));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_pair_rec)));
        r = CALL(bqd_rebuild, l, lo, hi);
        limdd_refs_pop(6);
    }
    cache_put3(CACHE_BQD_PAIR, node_or_zero(A1), key, kappa, (uint64_t)r);
    return scale(g, r);
}

/**
 * Perm of skip:alg:perm: the permutation `kind` of the qubits qa < qb, as
 * (pi f)(x) = f(pi x), memoised on the node, the kind and the two qubits, the
 * scalar outside. Above qa the node's two stored edges are permuted and
 * Rebuild makes the node, with no high cofactor; at qa the two cofactors are
 * f_0 and Cof1, and the new ones are Pairs at qb, or for the CX with the
 * control above f_0 and X_qb f_1; an edge that skips qa is a function that
 * does not depend on x_qa, whose image is itself for the CX with the control
 * below, and otherwise has the cofactors (f, X_qb f) or the two restrictions
 * at qb.
 */
TASK_4(BQD, bqd_perm_rec, BQD, e, uint32_t, kind, uint32_t, qa, uint32_t, qb)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG N = limdd_target(e);
    const uint32_t var = limdd_level(N);
    if (var > qa && kind == BQD_PERM_CX_UP) return e;
    const uint64_t key = (uint64_t)kind | ((uint64_t)qa << 8) | ((uint64_t)qb << 16);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_PERM, N, key, 0, &hit)) return scale(scalar_of(e), (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PERM);

    const BQD E = unit(N);
    BQD r;
    if (var > qa) {
        BQD g0, g1;
        if (kind == BQD_PERM_CX_DOWN) {
            g0 = limdd_refs_push(E);
            g1 = limdd_refs_push(CALL(bqd_xq_rec, E, qb));
        } else {
            limdd_refs_spawn(SPAWN(bqd_restrict_scalar, E, qb, 1));
            g0 = limdd_refs_push(CALL(bqd_restrict_scalar, E, qb, 0));
            g1 = limdd_refs_push(limdd_refs_sync(SYNC(bqd_restrict_scalar)));
        }
        r = CALL(bqd_compose_scalar, qa, g0, g1);
        limdd_refs_pop(2);
    } else if (var < qa) {
        limdd_refs_spawn(SPAWN(bqd_perm_rec, limdd_node_high(N), kind, qa, qb));
        const BQD lo = limdd_refs_push(CALL(bqd_perm_rec, limdd_node_low(N), kind, qa, qb));
        const BQD hi = limdd_refs_push(limdd_refs_sync(SYNC(bqd_perm_rec)));
        r = CALL(bqd_rebuild, var, lo, hi);
        limdd_refs_pop(2);
    } else {
        const BQD f0 = limdd_node_low(N);
        const BQD f1 = limdd_refs_push(CALL(bqd_high_cofactor, N));
        BQD g0, g1;
        if (kind == BQD_PERM_CX_DOWN) {
            g0 = f0;
            g1 = limdd_refs_push(CALL(bqd_xq_rec, f1, qb));
            limdd_refs_push(g0);
        } else {
            /* SW: (Pair(f0,0; f1,0), Pair(f0,1; f1,1)); CX up: (Pair(f0,0; f1,1), Pair(f1,0; f0,1)) */
            const bool sw = kind == BQD_PERM_SWAP;
            limdd_refs_spawn(SPAWN(bqd_pair_rec, sw ? f0 : f1, sw ? 1 : 0, sw ? f1 : f0, 1, qb));
            g0 = limdd_refs_push(CALL(bqd_pair_rec, f0, 0, f1, sw ? 0 : 1, qb));
            g1 = limdd_refs_push(limdd_refs_sync(SYNC(bqd_pair_rec)));
        }
        r = CALL(bqd_compose_scalar, qa, g0, g1);
        limdd_refs_pop(3);
    }
    cache_put3(CACHE_BQD_PERM, N, key, 0, (uint64_t)r);
    return scale(scalar_of(e), r);
}

/**
 * The least point of supp [t], a mask of qubits, down one path: the low edge
 * where it is not zero, and the high edge otherwise, whose support the high
 * cofactor has (skip:lem:minpoint with no translation).
 */
static uint64_t
minpoint_scalar(LIMDD_TARG t)
{
    uint64_t y = 0;
    while (t != LIMDD_TERMINAL) {
        const BQD lo = limdd_node_low(t);
        if (!limdd_edge_is_zero(lo)) { t = limdd_target(lo); continue; }
        y |= UINT64_C(1) << limdd_node_var(t);
        t = limdd_target(limdd_node_high(t));
    }
    return y;
}

/**
 * Exp of skip:alg:phasemul, the canonical edge of beta^eps, a function of full
 * support whose ratio at the root's variable is beta^{eps_1 - eps_0}. That is
 * not the constant 1, since eps is reduced modulo the order of beta; on floats
 * the order test can miss and a ratio can round to 1, and the level is then
 * skipped as Rebuild would. Memoised on beta and eps; the two calls overlap.
 */
TASK_IMPL_3(BQD, bqd_exp_state, BQD_EXP, eps, EVBDD_WGT, beta, uint32_t, r)
{
    if (bqd_sm()) return CALL(bqd_sm_exp, eps, beta, r);
    if (bqd_exp_is_const(eps)) return edge(bqd_exp_power(beta, bqd_exp_value(eps)), LIMDD_TERMINAL);
    uint64_t hit;
    if (cache_get3(CACHE_BQD_EXP, beta, eps, r, &hit)) return (BQD)hit;
    const uint32_t v = bqd_exp_var(eps);
    const BQD_EXP e0 = bqd_exp_cof(eps, v, 0), e1 = bqd_exp_cof(eps, v, 1);
    const BQD_EXP d = mtbdd_refs_push(CALL(bqd_exp_sub, e1, e0, r));
    limdd_refs_spawn(SPAWN(bqd_exp_state, e0, beta, r));
    const BQD h = limdd_refs_push(CALL(bqd_exp_state, d, beta, r));
    const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_exp_state)));
    const BQD res = (h == limdd_one_edge())
                  ? lo : edge(scalar_of(lo), limdd_makenode(v, unit(limdd_target(lo)), h));
    limdd_refs_pop(2);
    mtbdd_refs_pop(1);
    cache_put3(CACHE_BQD_EXP, beta, eps, r, (uint64_t)res);
    return res;
}

/**
 * PhaseMul of skip:alg:phasemul: the canonical edge of [k] . beta^eps on any
 * support, r the order of beta. Memoised on the node, eps and beta, the scalar
 * outside, so the node is split with the scale 1 and no MulOff. With g_0 =
 * k_0 beta^{eps_0} = alpha [U] and alpha = beta^v, v = eps_0 at the least
 * point of supp k_0, the high edge Rebuild would store is rho(k) times
 * beta^{eps_1 - eps_0} on supp k_0 and beta^{eps_1 - v} off it, so one call on
 * the ratio with that exponent makes it, the multiplication off the support
 * included. The two calls depend on the input alone, and overlap.
 */
TASK_IMPL_4(BQD, bqd_phase_mul_rec, BQD, k, BQD_EXP, eps, EVBDD_WGT, beta, uint32_t, r)
{
    if (limdd_edge_is_zero(k)) return k;
    if (bqd_exp_is_const(eps)) return scale(bqd_exp_power(beta, bqd_exp_value(eps)), k);
    const LIMDD_TARG N = limdd_target(k);
    if (N == LIMDD_TERMINAL) return scale(scalar_of(k), CALL(bqd_exp_state, eps, beta, r));
    uint64_t hit;
    if (cache_get3(CACHE_BQD_PHASEMUL, N, eps, beta, &hit)) return scale(scalar_of(k), (BQD)hit);
    bqd_count(LACE_WORKER_ID, BQD_COUNT_PHASEMUL);

    const uint32_t lN = limdd_level(N), le = bqd_exp_var(eps);
    const uint32_t l = lN < le ? lN : le;
    /* the low cofactor and the ratio of [N] at l: (N, ind N) where N skips l */
    const BQD k0 = lN > l ? unit(N) : limdd_node_low(N);
    const BQD kr = lN > l ? unit(CALL(bqd_ind_rec, N)) : limdd_node_high(N);
    const BQD_EXP e0 = bqd_exp_cof(eps, l, 0), e1 = bqd_exp_cof(eps, l, 1);
    BQD res;
    if (limdd_edge_is_zero(k0)) {
        const BQD h = limdd_refs_push(CALL(bqd_phase_mul_rec, kr, e1, beta, r));
        res = CALL(bqd_compose_scalar, l, k0, h);
        limdd_refs_pop(1);
    } else if (limdd_edge_is_zero(kr)) {
        const BQD lo = CALL(bqd_phase_mul_rec, k0, e0, beta, r);
        res = edge(scalar_of(lo), limdd_makenode(l, unit(limdd_target(lo)), kr));
    } else {
        const LIMDD_TARG U0 = limdd_target(k0);
        const int64_t v = bqd_exp_eval(e0, minpoint_scalar(U0));
        const BQD_EXP io = mtbdd_refs_push(CALL(bqd_iota_rec, U0));
        const BQD_EXP sel = mtbdd_refs_push(CALL(bqd_exp_sel, io, e0, v));
        const BQD_EXP e2 = mtbdd_refs_push(CALL(bqd_exp_sub, e1, sel, r));
        limdd_refs_spawn(SPAWN(bqd_phase_mul_rec, k0, e0, beta, r));
        const BQD h = limdd_refs_push(CALL(bqd_phase_mul_rec, kr, e2, beta, r));
        const BQD lo = limdd_refs_push(limdd_refs_sync(SYNC(bqd_phase_mul_rec)));
        const LIMDD_TARG U = limdd_target(lo);
        res = (limdd_label(h) == LIMDD_LIM_IDENTITY && limdd_target(h) == CALL(bqd_ind_rec, U))
            ? lo : edge(scalar_of(lo), limdd_makenode(l, unit(U), h));
        limdd_refs_pop(2);
        mtbdd_refs_pop(3);
    }
    cache_put3(CACHE_BQD_PHASEMUL, N, eps, beta, (uint64_t)res);
    return scale(scalar_of(k), res);
}

/* --- gates ---------------------------------------------------------------- */

/**
 * U on qubit q, conditioned on the controls in `cmask` (all above q) that are
 * still pending, which is Gate of skip:alg:gates. Linear, so the memo is on
 * the node and the scalar is multiplied back on. The work is at the first
 * level that matters, p: the highest pending control, or q when none is left.
 * Where the node is above p any variable passes the gate to both cofactors;
 * at p a control keeps its low cofactor and passes the gate, without that
 * control, to its high one, and at q the new cofactors are the rows of U
 * applied to the old ones. Where the node is below p, or is the terminal, the
 * edge skips p and is both of its cofactors there (skip:lem:virtual): a
 * skipped control gives (e, gate(e)), and a skipped target the row sums of U
 * times e.
 *
 * The key's third word is cmask with bit q set, which is unambiguous since
 * every control is below bit q. It holds no level: a control is dropped from
 * cmask once passed, so the key names the controls still to come, and with
 * the node it determines the result whatever level the edge is read at.
 */
TASK_4(BQD, bqd_cgate_rec, BQD, e, uint32_t, gid, uint64_t, cmask, uint32_t, q)
{
    if (limdd_edge_is_zero(e)) return e;
    const LIMDD_TARG t = limdd_target(e);
    const uint32_t var = limdd_level(t);
    const uint64_t key = cmask | (UINT64_C(1) << q);
    const uint32_t p = (uint32_t)__builtin_ctzll(key);

    uint64_t hit;
    if (cache_get3(CACHE_BQD_GATE, t, gid, key, &hit)) return scale(scalar_of(e), (BQD)hit);

    const EVBDD_WGT u00 = gates[gid][0], u01 = gates[gid][1];
    const EVBDD_WGT u10 = gates[gid][2], u11 = gates[gid][3];
    BQD r;
    if (var > p) {
        const BQD u = unit(t);
        if (p < q) {
            const BQD hi = CALL(bqd_cgate_rec, u, gid, cmask & ~(UINT64_C(1) << p), q);
            r = CALL(bqd_compose_scalar, p, u, hi);
        } else {
            r = CALL(bqd_compose_scalar, q, scale(wgt_add(u00, u01), u), scale(wgt_add(u10, u11), u));
        }
    } else {
        const BQD f0 = limdd_node_low(t);
        const BQD f1 = CALL(bqd_high_cofactor, t);
        BQD lo, hi;
        if (var < p) {
            limdd_refs_spawn(SPAWN(bqd_cgate_rec, f0, gid, cmask, q));
            hi = limdd_refs_push(CALL(bqd_cgate_rec, f1, gid, cmask, q));
            lo = limdd_refs_sync(SYNC(bqd_cgate_rec));
            limdd_refs_pop(1);
        } else if (var < q) {
            lo = f0;
            hi = CALL(bqd_cgate_rec, f1, gid, cmask & ~(UINT64_C(1) << var), q);
        } else {
            lo = CALL(bqd_apply_op, OP_ADD, scale(u00, f0), scale(u01, f1));
            hi = CALL(bqd_apply_op, OP_ADD, scale(u10, f0), scale(u11, f1));
        }
        r = CALL(bqd_compose_scalar, var, lo, hi);
    }
    cache_put3(CACHE_BQD_GATE, t, gid, key, (uint64_t)r);
    return scale(scalar_of(e), r);
}

/** diag(u00, u11) with both nonzero: a scalar times a monomial phase */
static inline bool
is_phase(uint32_t gid)
{
    return gates[gid][1] == EVBDD_ZERO && gates[gid][2] == EVBDD_ZERO &&
           gates[gid][0] != EVBDD_ZERO && gates[gid][3] != EVBDD_ZERO;
}

/** qubit q as a vector-index bit, the convention of bqd_apply_diagonal */
static inline uint64_t
vbit(uint32_t q, uint32_t n)
{
    return UINT64_C(1) << (n - 1 - q);
}

/** e . phase^{x_A}; bqd_apply_diagonal routes by support itself. */
static inline BQD
monomial_phase(BQD e, uint64_t A, EVBDD_WGT phase, uint32_t n)
{
    return bqd_apply_diagonal(e, A, phase, n);
}

/** [[0, u01], [u10, 0]] with both nonzero: u10 . X . diag(1, u01/u10), X and Y among them */
static inline bool
is_flip(uint32_t gid)
{
    return gates[gid][0] == EVBDD_ZERO && gates[gid][3] == EVBDD_ZERO &&
           gates[gid][1] != EVBDD_ZERO && gates[gid][2] != EVBDD_ZERO;
}

/* --- the entry points of skip:sec:ratio ------------------------------------ */

/*
 * Each dispatches on the family, as the other entry points do: the recursion
 * above for the scalar family, and qsylvan_bqd_xp.h for the other two.
 */
TASK_IMPL_4(BQD, bqd_perm, BQD, e, uint32_t, kind, uint32_t, qa, uint32_t, qb)
{
    assert(qa < qb && kind <= BQD_PERM_CX_UP);
    if (bqd_sm()) return CALL(bqd_sm_perm, e, kind, qa, qb);
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_perm, e, kind, qa, qb);
    return CALL(bqd_perm_rec, e, kind, qa, qb);
}

TASK_IMPL_2(BQD, bqd_x, BQD, e, uint32_t, q)
{
    if (bqd_sm()) return CALL(bqd_sm_xq, e, q);
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_x, e, q);
    return CALL(bqd_xq_rec, e, q);
}

TASK_IMPL_5(BQD, bqd_pair, BQD, A, int, s, BQD, B, int, t, uint32_t, b)
{
    if (bqd_sm()) return CALL(bqd_sm_pair, A, s, B, t, b);
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_pair, A, s, B, t, b);
    return CALL(bqd_pair_rec, A, s, B, t, b);
}

TASK_IMPL_1(LIMDD_TARG, bqd_ind, LIMDD_TARG, t)
{
    if (bqd_sm()) return CALL(bqd_sm_ind, t);
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_ind, t);
    return CALL(bqd_ind_rec, t);
}

TASK_IMPL_3(BQD, bqd_mul_off, BQD, k, EVBDD_WGT, c, LIMDD_TARG, u)
{
    if (bqd_sm()) return CALL(bqd_sm_mul_off, k, c, u);
    switch (bqd_family()) {
    case BQD_FAMILY_SCALAR: return CALL(bqd_mul_off_rec, k, c, u);
    case BQD_FAMILY_X:      return CALL(bqd_xp_mul_off, k, c, u);
    default:
        fprintf(stderr, "sylvan: MulOff is the scalar and the translation family's\n");
        exit(1);
    }
}

TASK_IMPL_3(BQD, bqd_phase_mul, BQD, k, BQD_EXP, eps, EVBDD_WGT, beta)
{
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_phase_mul, k, eps, beta);
    if (beta == EVBDD_ONE) return k;
    if (bqd_sm()) return CALL(bqd_sm_phase_mul, k, eps, beta, bqd_exp_order(beta));
    return CALL(bqd_phase_mul_rec, k, eps, beta, bqd_exp_order(beta));
}

/*
 * DiagR of skip:alg:phasemul and its counterparts in the other families: the
 * phase multiplication by the monomial, on any support. In the Pauli family a
 * phase that is not a power of w_8 has no exponent in its base, and takes the
 * product with the monomial as before. A task, so that the RUN in the header
 * makes it callable from any thread, as the other entry points are.
 */
TASK_IMPL_4(BQD, bqd_diag_any, BQD, e, uint64_t, A, EVBDD_WGT, phase, uint32_t, n)
{
    if (limdd_edge_is_zero(e) || phase == EVBDD_ONE) return e;
    if (A == 0) return bqd_scale(e, phase);
    if (bqd_sm()) return CALL(bqd_sm_diag, e, bqd_sm_var_mask(A, n), phase);
    if (bqd_family() == BQD_FAMILY_PAULI && bqd_xp_w8_log(phase) < 0) {
        limdd_refs_push(e);
        const BQD m = limdd_refs_push(bqd_monomial(A, phase, n));
        const BQD res = CALL(bqd_multiply, e, m);
        limdd_refs_pop(2);
        return res;
    }
    uint64_t vars = 0;                               /* qubit q is bit n-1-q of A */
    for (uint64_t rest = A; rest != 0; rest &= rest - 1)
        vars |= UINT64_C(1) << (n - 1 - (uint32_t)__builtin_ctzll(rest));
    const uint32_t r = bqd_family() == BQD_FAMILY_PAULI ? 8 : bqd_exp_order(phase);
    limdd_refs_push(e);
    const BQD_EXP eps = mtbdd_refs_push(bqd_exp_monomial(vars, 1, r));
    const BQD res = CALL(bqd_phase_mul, e, eps, phase);
    mtbdd_refs_pop(1);
    limdd_refs_pop(1);
    return res;
}

/*
 * The two entry points serve every family. The monomial goes through
 * bqd_apply_diagonal, which picks the product for the family; X, Y and every
 * other flip through X of skip:alg:perm after the diagonal gate that makes it
 * one; and any other gate through the family's recursion: bqd_cgate_rec here,
 * and Gate of qsylvan_bqd_xp.h, on the top of a labelled edge, in the other
 * two. The CX with one control is Perm, and a controlled flip is Perm between
 * two monomials: C(u10 X D) = u10^{x_c} . CX . p^{x_c x_t} for D = diag(1, p).
 */
TASK_IMPL_4(BQD, bqd_gate, BQD, e, uint32_t, gid, uint32_t, target, uint32_t, n)
{
    if (is_phase(gid)) {
        /* diag(u00, u11) = u00 . diag(1, u11/u00), a scalar and a monomial */
        const EVBDD_WGT u00 = gates[gid][0];
        const EVBDD_WGT p = wgt_div(gates[gid][3], u00);
        const BQD d = (p == EVBDD_ONE) ? e : monomial_phase(e, vbit(target, n), p, n);
        return bqd_scale(d, u00);
    }
    if (is_flip(gid)) {
        /* [[0, u01], [u10, 0]] = u10 . X . diag(1, u01/u10) */
        const EVBDD_WGT u10 = gates[gid][2];
        const EVBDD_WGT p = wgt_div(gates[gid][1], u10);
        const BQD d = limdd_refs_push(p == EVBDD_ONE ? e : monomial_phase(e, vbit(target, n), p, n));
        const BQD x = CALL(bqd_x, d, target);
        limdd_refs_pop(1);
        return bqd_scale(x, u10);
    }
    if (bqd_sm()) return CALL(bqd_sm_gate, e, gid, 0, target);
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_cgate_rec, e, gid, 0, target);
    return CALL(bqd_cgate_rec, e, gid, 0, target);
}

/**
 * The flip gid on `target` controlled by the one qubit `control`, on either
 * side: u10^{x_c} . CX . p^{x_c x_t}, the CX being Perm with the control above
 * or below. A task, since bqd_cgate_either, which is not one, calls it with RUN.
 */
TASK_5(BQD, bqd_controlled_flip, BQD, e, uint32_t, gid, uint32_t, control, uint32_t, target,
       uint32_t, n)
{
    const EVBDD_WGT u10 = gates[gid][2];
    const EVBDD_WGT p = wgt_div(gates[gid][1], u10);
    BQD d = limdd_refs_push(p == EVBDD_ONE ? e
                            : monomial_phase(e, vbit(control, n) | vbit(target, n), p, n));
    const BQD x = (control < target) ? CALL(bqd_perm, d, BQD_PERM_CX_DOWN, control, target)
                                     : CALL(bqd_perm, d, BQD_PERM_CX_UP, target, control);
    limdd_refs_pop(1);
    if (u10 == EVBDD_ONE) return x;
    limdd_refs_push(x);
    d = monomial_phase(x, vbit(control, n), u10, n);
    limdd_refs_pop(1);
    return d;
}

TASK_IMPL_5(BQD, bqd_cgate, BQD, e, uint32_t, gid, uint64_t, cmask, uint32_t, target,
            uint32_t, n)
{
    if (cmask == 0) return CALL(bqd_gate, e, gid, target, n);
    if (is_phase(gid) && gates[gid][0] == EVBDD_ONE) {
        /* phase u11 where every control and the target are 1: one monomial */
        uint64_t A = vbit(target, n);
        for (uint32_t c = 0; c < n; c++) if ((cmask >> c) & 1) A |= vbit(c, n);
        return monomial_phase(e, A, gates[gid][3], n);
    }
    if ((cmask >> target) != 0) {
        fprintf(stderr, "sylvan: bqd_cgate needs every control above the target "
                        "for a gate that is not a monomial phase\n");
        exit(1);
    }
    if (is_flip(gid) && (cmask & (cmask - 1)) == 0)
        return CALL(bqd_controlled_flip, e, gid, (uint32_t)__builtin_ctzll(cmask), target, n);
    if (bqd_sm()) return CALL(bqd_sm_gate, e, gid, cmask, target);
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_cgate_rec, e, gid, cmask, target);
    return CALL(bqd_cgate_rec, e, gid, cmask, target);
}

/*
 * A control below the target: a monomial phase has its controls anywhere, and
 * a flip is the controlled flip above, whose CX is Perm with the control
 * below, CX_ba of skip:alg:perm; no Hadamard conjugation.
 */
BQD
bqd_cgate_either(BQD e, uint32_t gid, uint32_t control, uint32_t target, uint32_t n, bool *ok)
{
    *ok = true;
    if (control < target) return bqd_cgate(e, gid, UINT64_C(1) << control, target, n);
    if (is_phase(gid) && gates[gid][0] == EVBDD_ONE)
        return bqd_cgate(e, gid, UINT64_C(1) << control, target, n);
    if (is_flip(gid)) return RUN(bqd_controlled_flip, e, gid, control, target, n);
    *ok = false;
    return e;
}

/* The swap is Perm, SW of skip:alg:perm, and no longer three CNOTs. */
BQD
bqd_swap(BQD e, uint32_t a, uint32_t b, uint32_t n)
{
    (void)n;
    if (a == b) return e;
    if (a > b) { const uint32_t t = a; a = b; b = t; }
    return bqd_perm(e, BQD_PERM_SWAP, a, b);
}
