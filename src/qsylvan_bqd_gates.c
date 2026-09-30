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
 * recursion calls the scalar tasks directly and never tests the family.
 */
TASK_IMPL_2(BQD, bqd_multiply, BQD, f, BQD, g)
{
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_apply_op, BQD_OP_MUL, f, g);
    return CALL(bqd_apply_op, OP_MUL, f, g);
}

TASK_IMPL_2(BQD, bqd_add, BQD, f, BQD, g)
{
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_apply_op, BQD_OP_ADD, f, g);
    return CALL(bqd_apply_op, OP_ADD, f, g);
}

TASK_IMPL_3(BQD, bqd_compose, uint32_t, var, BQD, a, BQD, b)
{
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
 * Restriction and projection, both linear, so memoised on the node. Above q
 * the two cofactors are treated and composed again. At q the restriction is
 * the chosen cofactor itself: it does not depend on x_q, so as an edge read
 * at q it skips q and needs no node there. The projection composes the
 * chosen cofactor with zero. An edge whose node is below q, or is the
 * terminal, skips q and denotes a function that does not depend on x_q
 * (skip:alg:gates): it is its own restriction, and its projection is the
 * node at q with the edge on one side and zero on the other. Both are
 * answered before the lookup, in O(1).
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
        const BQD f1 = CALL(bqd_high_cofactor, t);
        limdd_refs_spawn(SPAWN(bqd_restrict_scalar, f1, q, b));
        const BQD lo = limdd_refs_push(CALL(bqd_restrict_scalar, f0, q, b));
        const BQD hi = limdd_refs_sync(SYNC(bqd_restrict_scalar));
        limdd_refs_pop(1);
        r = CALL(bqd_compose_scalar, var, lo, hi);
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
        const BQD f1 = CALL(bqd_high_cofactor, t);
        limdd_refs_spawn(SPAWN(bqd_project_scalar, f1, q, b));
        const BQD lo = limdd_refs_push(CALL(bqd_project_scalar, f0, q, b));
        const BQD hi = limdd_refs_sync(SYNC(bqd_project_scalar));
        limdd_refs_pop(1);
        r = CALL(bqd_compose_scalar, var, lo, hi);
    }
    cache_put3(CACHE_BQD_PROJECT, t, key, 0, (uint64_t)r);
    return scale(scalar_of(e), r);
}

/* In the other two families a cofactor is a labelled edge, so the chosen one
 * goes through Canon where the scalar family returns it (skip:alg:xrestrict). */
TASK_IMPL_3(BQD, bqd_restrict, BQD, e, uint32_t, q, int, b)
{
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_restrict, e, q, b);
    return CALL(bqd_restrict_scalar, e, q, b);
}

TASK_IMPL_3(BQD, bqd_project, BQD, e, uint32_t, q, int, b)
{
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
    return RUN(bqd_full_rec, limdd_target(e)) != 0;
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

/*
 * The two entry points serve every family. The monomial goes through
 * bqd_apply_diagonal, which picks the product for the family, and any other
 * gate through the family's recursion: bqd_cgate_rec here, and Gate of
 * qsylvan_bqd_xp.h, on the top of a labelled edge, in the other two.
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
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_cgate_rec, e, gid, 0, target);
    return CALL(bqd_cgate_rec, e, gid, 0, target);
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
    if (bqd_family() != BQD_FAMILY_SCALAR) return CALL(bqd_xp_cgate_rec, e, gid, cmask, target);
    return CALL(bqd_cgate_rec, e, gid, cmask, target);
}

static BQD
cx_reversed(BQD e, uint32_t control, uint32_t target, uint32_t n)
{
    /* CNOT_{control -> target} with control > target, as
     * (H (x) H) CNOT_{target -> control} (H (x) H) */
    e = bqd_gate(e, GATEID_H, control, n);
    e = bqd_gate(e, GATEID_H, target, n);
    e = bqd_cgate(e, GATEID_X, UINT64_C(1) << target, control, n);
    e = bqd_gate(e, GATEID_H, target, n);
    e = bqd_gate(e, GATEID_H, control, n);
    return e;
}

BQD
bqd_cgate_either(BQD e, uint32_t gid, uint32_t control, uint32_t target, uint32_t n, bool *ok)
{
    *ok = true;
    if (control < target) return bqd_cgate(e, gid, UINT64_C(1) << control, target, n);
    if (is_phase(gid) && gates[gid][0] == EVBDD_ONE)
        return bqd_cgate(e, gid, UINT64_C(1) << control, target, n);
    if (gid == GATEID_X) return cx_reversed(e, control, target, n);
    *ok = false;
    return e;
}

BQD
bqd_swap(BQD e, uint32_t a, uint32_t b, uint32_t n)
{
    if (a == b) return e;
    if (a > b) { const uint32_t t = a; a = b; b = t; }
    e = bqd_cgate(e, GATEID_X, UINT64_C(1) << a, b, n);
    e = cx_reversed(e, b, a, n);
    e = bqd_cgate(e, GATEID_X, UINT64_C(1) << a, b, n);
    return e;
}
