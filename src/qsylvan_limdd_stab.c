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
#include <string.h>

#include <sylvan_int.h>
#include <sylvan_platform.h>

#include "qsylvan_limdd_stab.h"

/* Written once at init, read-only afterwards. */
static llmsset_t stab_table = NULL;

/**
 * A generator during reduction: a Pauli word plus a sign.
 *
 * Reduction works on this unpacked form rather than on interned LIMs, because
 * row operations touch every row repeatedly and interning each intermediate
 * would put a table lookup inside the inner loop of an O(n^3) algorithm. Only
 * the final rows are interned.
 */
typedef struct {
    limdd_pauli_t p;
    bool neg;
} gen_row_t;

/*
 * Reduction scratch, per worker.
 *
 * Thread-local, so it is private by construction: two threads' TLS blocks are
 * separate allocations and cannot share a cache line, which is the reason this
 * is not a pool indexed by worker id. Sized for a redundant input -- callers
 * building a group from two children hand over both children's generators at
 * once, so the input can exceed the n generators the result can have.
 */
#define STAB_SCRATCH_ROWS (2 * LIMDD_MAX_QUBITS)
static SYLVAN_TLS gen_row_t scratch[STAB_SCRATCH_ROWS];

/**
 * A row of the subgroup-intersection elimination: a Pauli word together with
 * the set of Stab(v0) generators that produced it.
 *
 * Only the word is eliminated; the mask rides along, so when a row's word
 * cancels to the identity the mask names a product of Stab(v0) generators that
 * also lies in Stab(v1). That is one generator of the intersection.
 */
/*
 * Which generators went into a row.
 *
 * One bit per generator, and a group on n qubits has up to n of them, so a
 * single word only reaches 64 qubits. It used to be a single word, guarded by
 * an assertion -- which compiles out, leaving `1 << i` with i >= 64 as
 * undefined behaviour and a silently wrong answer above 64 qubits, exactly
 * where a wide Pauli word is the point. It is now as wide as the Pauli word.
 */
typedef struct { uint64_t w[LIMDD_PAULI_WORDS]; } genmask_t;

static inline void genmask_clear(genmask_t *m)
{ for (unsigned i = 0; i < LIMDD_PAULI_WORDS; i++) m->w[i] = 0; }

static inline void genmask_set(genmask_t *m, size_t b)
{ genmask_clear(m); m->w[LIMDD_PAULI_LANE(b)] = LIMDD_PAULI_BIT(b); }

static inline void genmask_xor(genmask_t *a, const genmask_t *b)
{ for (unsigned i = 0; i < LIMDD_PAULI_WORDS; i++) a->w[i] ^= b->w[i]; }

static inline bool genmask_is_zero(const genmask_t *m)
{ for (unsigned i = 0; i < LIMDD_PAULI_WORDS; i++) if (m->w[i]) return false; return true; }

static inline bool genmask_test(const genmask_t *m, size_t b)
{ return (m->w[LIMDD_PAULI_LANE(b)] & LIMDD_PAULI_BIT(b)) != 0; }

typedef struct {
    limdd_pauli_t w;
    genmask_t mask0;  /* which generators of s0 went into this row */
    genmask_t mask1;  /* and which of s1 */
    size_t   pivot;   /* leading column, for rows 0..top-1 of the RREF */
} isect_row_t;

static SYLVAN_TLS isect_row_t isect[STAB_SCRATCH_ROWS];

static void
die(void)
{
    fprintf(stderr, "sylvan: LIMDD stabiliser table is full\n");
    exit(1);
}

size_t
limdd_stab_table_size(void)
{
    return stab_table == NULL ? 0 : llmsset_get_size(stab_table);
}

size_t
limdd_stab_table_count(void)
{
    if (stab_table == NULL) return 0;
    const size_t occupied = (size_t)llmsset_count_marked(stab_table);
    return occupied < 2 ? 0 : occupied - 2;
}

/* --- the packed list ----------------------------------------------------- */

static inline LIMDD_STAB
stab_cons(LIMDD_LIM head, LIMDD_STAB tail)
{
    int created;
    const uint64_t s = llmsset_lookup(stab_table, head, tail, &created);
    if (s == 0) die();
    return s;
}

static inline LIMDD_LIM
stab_head(LIMDD_STAB s)
{
    assert(s != LIMDD_STAB_TRIVIAL);
    return ((const uint64_t *)llmsset_index_to_ptr(stab_table, s))[0];
}

static inline LIMDD_STAB
stab_tail(LIMDD_STAB s)
{
    assert(s != LIMDD_STAB_TRIVIAL);
    return ((const uint64_t *)llmsset_index_to_ptr(stab_table, s))[1];
}

size_t
limdd_stab_ngens(LIMDD_STAB s)
{
    size_t n = 0;
    while (s != LIMDD_STAB_TRIVIAL) { n++; s = stab_tail(s); }
    return n;
}

LIMDD_LIM
limdd_stab_head(LIMDD_STAB s)
{
    return stab_head(s);
}

LIMDD_STAB
limdd_stab_tail(LIMDD_STAB s)
{
    return stab_tail(s);
}

LIMDD_LIM
limdd_stab_gen(LIMDD_STAB s, size_t i)
{
    while (i-- > 0) s = stab_tail(s);
    return stab_head(s);
}

/*
 * A group is a linked list, so limdd_stab_gen(s, i) walks from the head every
 * time and reading all k generators by index costs O(k^2) pointer chases.
 * That showed up as the second hottest function in the profile, under
 * stack_generators. Everything below walks the list once instead. The ORDER
 * is preserved exactly: Pauli products do not commute, and the callers that
 * multiply tagged generators back together recover the phase from it.
 */

/* --- reduction ----------------------------------------------------------- */

/**
 * Column `c` of the symplectic matrix: columns 0..n-1 are the X bits of qubits
 * 0..n-1, columns n..2n-1 the Z bits. Any fixed order gives a canonical RREF;
 * this is the usual (X|Z) one.
 */
static inline bool
row_bit(const gen_row_t *r, size_t c, size_t nqubits)
{
    return limdd_pauli_column(r->p, c, nqubits);
}

/**
 * dst <- src * dst, as Pauli operators.
 *
 * The product of two Pauli words carries a factor i^k. Here k is always 0 or
 * 2: elements of a stabiliser group commute, and for commuting Pauli words the
 * two orderings agree, which forces the factor to be real. A k of 1 or 3 means
 * the caller passed anticommuting generators, so the set does not generate a
 * stabiliser group and the assertion is the right response.
 */
static void
add_row(gen_row_t *dst, const gen_row_t *src, size_t nqubits)
{
    limdd_pauli_t prod = src->p;
    const unsigned k = limdd_pauli_rightmul(&prod, dst->p, nqubits);
    assert((k & 1) == 0 && "stabiliser generators must commute");

    dst->p = prod;
    dst->neg = dst->neg ^ src->neg ^ (k == 2);
}

/**
 * Reduce `rows[0..n)` to row reduced echelon form in place, returning the
 * number of independent rows, which are left at the front.
 *
 * Note that elimination runs over EVERY other row, not only those below the
 * pivot. Eliminating downwards alone gives row echelon form, which is not
 * canonical: the same group then has several echelon forms and interning them
 * would hand out several handles for one group. (This is exactly the bug found
 * in the reference implementation and fixed there.)
 */
static size_t
rref(gen_row_t *rows, size_t n, size_t nqubits)
{
    size_t top = 0;
    for (size_t c = 0; c < 2 * nqubits && top < n; c++) {
        size_t pivot = top;
        while (pivot < n && !row_bit(&rows[pivot], c, nqubits)) pivot++;
        if (pivot == n) continue;

        if (pivot != top) {
            const gen_row_t t = rows[top];
            rows[top] = rows[pivot];
            rows[pivot] = t;
        }

        for (size_t r = 0; r < n; r++) {
            if (r != top && row_bit(&rows[r], c, nqubits)) {
                add_row(&rows[r], &rows[top], nqubits);
            }
        }
        top++;
    }

    /*
     * Dependent rows reduce to the identity word. A -I among them would mean
     * the generators multiply to -I, which no state is stabilised by; such a
     * set is not a stabiliser group and the caller has a bug.
     */
    for (size_t r = top; r < n; r++) {
        assert(limdd_pauli_is_identity(rows[r].p));
        assert(!rows[r].neg && "generators multiply to -I; not a stabiliser group");
    }
    return top;
}

/* --- construction -------------------------------------------------------- */

static gen_row_t
unpack(LIMDD_LIM g)
{
    const EVBDD_WGT w = limdd_lim_weight(g);
    assert((w == EVBDD_ONE || w == EVBDD_MIN_ONE)
           && "a stabiliser generator's scalar must be +1 or -1");

    gen_row_t r;
    r.p = limdd_lim_pauli(g);
    r.neg = (w == EVBDD_MIN_ONE);
    return r;
}

static LIMDD_LIM
pack(gen_row_t r)
{
    return limdd_lim_make(r.p, r.neg ? EVBDD_MIN_ONE : EVBDD_ONE);
}

LIMDD_STAB
limdd_stab_make(const LIMDD_LIM *gens, size_t ngens)
{
    assert(stab_table != NULL);
    assert(ngens <= STAB_SCRATCH_ROWS);

    const size_t nqubits = limdd_lims_nqubits();

    for (size_t i = 0; i < ngens; i++) scratch[i] = unpack(gens[i]);

    const size_t k = rref(scratch, ngens, nqubits);

    /* Built back to front, so the list reads in RREF order. */
    LIMDD_STAB s = LIMDD_STAB_TRIVIAL;
    for (size_t i = k; i-- > 0; ) s = stab_cons(pack(scratch[i]), s);
    return s;
}

LIMDD_STAB
limdd_stab_extend_skipped(LIMDD_STAB s, uint32_t from, uint32_t to)
{
    assert(from <= to && to <= limdd_lims_nqubits());
#ifndef NDEBUG
    /* The prepend is only RREF if nothing in `s` acts below `to`. */
    const uint64_t below = (to >= 64) ? ~UINT64_C(0) : ((UINT64_C(1) << to) - 1);
    for (LIMDD_STAB c = s; c != LIMDD_STAB_TRIVIAL; c = stab_tail(c)) {
        const limdd_pauli_t p = limdd_lim_pauli(stab_head(c));
        assert(!limdd_pauli_acts_below(p, to)
               && "group acts below the level it is extended to");
    }
#endif
    for (uint32_t j = to; j-- > from; ) {
        const limdd_pauli_t xj = limdd_pauli_single(j, LIMDD_PAULI_X);
        s = stab_cons(limdd_lim_make(xj, EVBDD_ONE), s);
    }
    return s;
}

bool
limdd_stab_contains(LIMDD_STAB s, LIMDD_LIM g)
{
    const size_t nqubits = limdd_lims_nqubits();
    gen_row_t r = unpack(g);

    /*
     * Reduce `r` by each row in turn. The rows are in RREF, so each owns a
     * pivot column no other row touches: cancelling that column can never undo
     * an earlier cancellation, and one pass over the rows suffices.
     */
    for (; s != LIMDD_STAB_TRIVIAL; s = stab_tail(s)) {
        const gen_row_t row = unpack(stab_head(s));

        size_t c = 0;
        while (c < 2 * nqubits && !row_bit(&row, c, nqubits)) c++;
        assert(c < 2 * nqubits && "an RREF row is never the identity");

        if (row_bit(&r, c, nqubits)) add_row(&r, &row, nqubits);
    }

    return limdd_pauli_is_identity(r.p) && !r.neg;
}

LIMDD_LIM
limdd_stab_element(LIMDD_STAB s, uint64_t k)
{
    LIMDD_LIM acc = LIMDD_LIM_IDENTITY;
    /* Consume k rather than index it. The loop runs once per generator, not
     * once per set bit, so `k >> i` was evaluated for every i below the
     * generator count -- undefined once a group can have 64 or more, which a
     * wider Pauli word allows. Shifting k down is defined at any width and
     * needs no bound on the caller. */
    for (; s != LIMDD_STAB_TRIVIAL; s = stab_tail(s), k >>= 1) {
        if (k & 1) acc = limdd_lim_mul(acc, stab_head(s));
    }
    return acc;
}

void
limdd_stab_fprint(FILE *out, LIMDD_STAB s, const char *indent)
{
    if (s == LIMDD_STAB_TRIVIAL) {
        fprintf(out, "%s<trivial>\n", indent);
        return;
    }
    for (; s != LIMDD_STAB_TRIVIAL; s = stab_tail(s)) {
        fprintf(out, "%s", indent);
        limdd_lim_fprint(out, stab_head(s));
        fprintf(out, "\n");
    }
}

/* --- the stabiliser group of a node -------------------------------------- */

/** The element of `s` whose Pauli word is `w`, or 0 if the group has none. */
static LIMDD_LIM
stab_find(LIMDD_STAB s, limdd_pauli_t w)
{
    const size_t nqubits = limdd_lims_nqubits();
    LIMDD_LIM acc = LIMDD_LIM_IDENTITY;

    for (; s != LIMDD_STAB_TRIVIAL; s = stab_tail(s)) {
        const LIMDD_LIM g = stab_head(s);
        const gen_row_t row = unpack(g);

        size_t c = 0;
        while (c < 2 * nqubits && !row_bit(&row, c, nqubits)) c++;
        assert(c < 2 * nqubits && "an RREF row is never the identity");

        gen_row_t cur = { w, false };
        if (row_bit(&cur, c, nqubits)) {
            limdd_pauli_t prod = w;
            limdd_pauli_rightmul(&prod, row.p, nqubits);
            w = prod;
            acc = limdd_lim_mul(acc, g);
        }
    }

    /*
     * `acc` is the product of the generators that were used, so its word is
     * the original `w` exactly when what remains is the identity. The signs
     * take care of themselves: the generators commute, so limdd_lim_mul never
     * produces an imaginary factor here.
     */
    return limdd_pauli_is_identity(w) ? acc : 0;
}

/** `a` with its sign flipped. `a` must have scalar +-1. */
static LIMDD_LIM
lim_negate(LIMDD_LIM a)
{
    return limdd_lim_make(limdd_lim_pauli(a), wgt_neg(limdd_lim_weight(a)));
}

/** Q a Q, i.e. `a` negated iff its word anticommutes with `q`. */
static LIMDD_LIM
conjugate(LIMDD_LIM a, limdd_pauli_t q)
{
    return (limdd_pauli_commutation_phase(limdd_lim_pauli(a), q) == 0)
         ? a : lim_negate(a);
}

/** True iff `w` and `q` commute. */
static inline bool
commutes(limdd_pauli_t w, limdd_pauli_t q)
{
    return limdd_pauli_commutation_phase(w, q) == 0;
}

/**
 * Generators of { A in s0 : word(A) lies in the span of s1's words }, written
 * into `out`; returns how many.
 *
 * Standard Zassenhaus: stack s0's rows carrying a tag and s1's rows carrying
 * none, eliminate the words, and read the tags off the rows that cancelled.
 */
static size_t
stack_generators(LIMDD_STAB s0, LIMDD_STAB s1, size_t *nrows)
{
    const size_t nqubits = limdd_lims_nqubits();
    const size_t k0 = limdd_stab_ngens(s0);
    const size_t k1 = limdd_stab_ngens(s1);
    assert(k0 + k1 <= STAB_SCRATCH_ROWS);
    assert(k0 <= LIMDD_MAX_QUBITS && k1 <= LIMDD_MAX_QUBITS);

    size_t n = 0;
    {
        LIMDD_STAB c = s0;
        for (size_t i = 0; i < k0; i++, c = stab_tail(c)) {
            isect[n].w = limdd_lim_pauli(stab_head(c));
            genmask_set(&isect[n].mask0, i);
            genmask_clear(&isect[n].mask1);
            n++;
        }
        c = s1;
        for (size_t j = 0; j < k1; j++, c = stab_tail(c)) {
            isect[n].w = limdd_lim_pauli(stab_head(c));
            genmask_clear(&isect[n].mask0);
            genmask_set(&isect[n].mask1, j);
            n++;
        }
    }

    /*
     * Eliminate the words only; the tags ride along. Signs are deliberately
     * not tracked -- both callers recover them afterwards by multiplying the
     * tagged generators back together, which gets the phases right for free.
     */
    size_t top = 0;
    for (size_t c = 0; c < 2 * nqubits && top < n; c++) {
        size_t pivot = top;
        while (pivot < n) {
            const gen_row_t r = { isect[pivot].w, false };
            if (row_bit(&r, c, nqubits)) break;
            pivot++;
        }
        if (pivot == n) continue;

        if (pivot != top) {
            const isect_row_t t = isect[top];
            isect[top] = isect[pivot];
            isect[pivot] = t;
        }

        for (size_t r = 0; r < n; r++) {
            if (r == top) continue;
            const gen_row_t g = { isect[r].w, false };
            if (!row_bit(&g, c, nqubits)) continue;
            limdd_pauli_xor(&isect[r].w, isect[top].w);
            genmask_xor(&isect[r].mask0, &isect[top].mask0);
            genmask_xor(&isect[r].mask1, &isect[top].mask1);
        }
        /* Recorded here, where it is already known. limdd_stab_min_coset
         * would otherwise rescan each row from column 0 to find it again --
         * O(top * 2n) bit tests on every coset minimisation, and this is the
         * hottest function in the profile. */
        isect[top].pivot = c;
        top++;
    }

    *nrows = n;
    return top;
}

/** The product of the generators of `s` named by `mask`, in index order. */
static LIMDD_LIM
product_of(LIMDD_STAB s, const genmask_t *mask)
{
    LIMDD_LIM acc = LIMDD_LIM_IDENTITY;
    for (size_t i = 0; s != LIMDD_STAB_TRIVIAL; i++, s = stab_tail(s)) {
        if (genmask_test(mask, i)) acc = limdd_lim_mul(acc, stab_head(s));
    }
    return acc;
}

/**
 * Generators of { A in s0 : word(A) lies in the span of s1's words }, written
 * into `out`; returns how many.
 *
 * Standard Zassenhaus: the rows that cancel to the identity are combinations
 * of s0 and s1 rows that agree, and their s0 tag names the element.
 */
static size_t
intersect_gens(LIMDD_STAB s0, LIMDD_STAB s1, LIMDD_LIM *out)
{
    size_t n;
    const size_t top = stack_generators(s0, s1, &n);

    size_t nout = 0;
    for (size_t r = top; r < n; r++) {
        assert(limdd_pauli_is_identity(isect[r].w));
        if (genmask_is_zero(&isect[r].mask0)) continue;  /* inside s1 alone */
        out[nout++] = product_of(s0, &isect[r].mask0);
    }
    return nout;
}

/**
 * True iff `w` is the sign the rule prefers: nonnegative imaginary part, and
 * among the two purely real options the nonnegative one.
 */
static bool
sign_is_canonical(EVBDD_WGT w)
{
    const complex_t c = weight_as_complex(w);
    if (c.i > 1e-14) return true;
    if (c.i < -1e-14) return false;
    return c.r >= 0.0;
}

LIMDD_LIM
limdd_stab_min_coset(LIMDD_LIM b, LIMDD_STAB s0, LIMDD_STAB s1,
                     LIMDD_LIM *witness, bool *negated)
{
    assert(stab_table != NULL);
    assert(!limdd_lim_is_zero(b) && "the zero map has no coset to minimise");

    /*
     * Memoised on its three arguments, all of which are interned handles, so
     * the key names the problem exactly. Worth it because the expensive part
     * -- stack_generators, which runs a full elimination over both groups --
     * depends only on s0 and s1, and makeedge asks for the same pair of child
     * groups over and over as a diagram is built. Two results, so this uses
     * the six-word cache: the representative, and the witness with the sign
     * flip in its top bit.
     */
    {
        uint64_t r1, r2;
        if (cache_get6(CACHE_LIMDD_MINCOSET, b, s0, s1, 0, 0, &r1, &r2)) {
            if (witness != NULL) *witness = (LIMDD_LIM)(r2 & ~(UINT64_C(1) << 63));
            if (negated != NULL) *negated = (r2 >> 63) != 0;
            return (LIMDD_LIM)r1;
        }
    }

    const size_t nqubits = limdd_lims_nqubits();

    size_t n;
    const size_t top = stack_generators(s0, s1, &n);

    /*
     * Reduce b's word by the basis. The basis is in RREF, so clearing every
     * pivot column leaves the one element of the coset with zeros there, which
     * is the least under the column order -- any other element differs by some
     * combination of basis rows and so carries a 1 at that combination's
     * leading pivot.
     */
    gen_row_t cur = { limdd_lim_pauli(b), false };
    genmask_t mask0, mask1;
    genmask_clear(&mask0); genmask_clear(&mask1);

    for (size_t r = 0; r < top; r++) {
        const size_t c = isect[r].pivot;     /* recorded by stack_generators */
        assert(c < 2 * nqubits && "an eliminated row is never the identity");
#ifndef NDEBUG
        {   /* the recorded pivot must still be the row's leading column */
            const gen_row_t row = { isect[r].w, false };
            size_t chk = 0;
            while (chk < 2 * nqubits && !row_bit(&row, chk, nqubits)) chk++;
            assert(chk == c);
        }
#endif
        if (!row_bit(&cur, c, nqubits)) continue;
        limdd_pauli_xor(&cur.p, isect[r].w);
        genmask_xor(&mask0, &isect[r].mask0);
        genmask_xor(&mask1, &isect[r].mask1);
    }

    /*
     * The tags name the G and H that produce that word. Multiplying them back
     * out in order is what recovers the phase: the words alone cannot, since
     * XOR forgets every factor of i the products pick up.
     */
    const LIMDD_LIM g = product_of(s0, &mask0);
    const LIMDD_LIM h = product_of(s1, &mask1);
    LIMDD_LIM e = limdd_lim_mul(limdd_lim_mul(g, b), h);
    assert(limdd_pauli_equals(limdd_lim_pauli(e), cur.p));

    const bool flip = !sign_is_canonical(limdd_lim_weight(e));
    if (flip) e = limdd_lim_make(limdd_lim_pauli(e), wgt_neg(limdd_lim_weight(e)));

    cache_put6(CACHE_LIMDD_MINCOSET, b, s0, s1, 0, 0,
               (uint64_t)e, (uint64_t)g | ((uint64_t)flip << 63));

    if (witness != NULL) *witness = g;
    if (negated != NULL) *negated = flip;
    return e;
}

LIMDD_STAB
limdd_stab_of_node(uint32_t var, LIMDD low, LIMDD high,
                   LIMDD_STAB s0, LIMDD_STAB s1)
{
    assert(stab_table != NULL);
    assert(var < limdd_lims_nqubits());

    const bool low_zero  = limdd_lim_is_zero(limdd_label(low));
    const bool high_zero = limdd_lim_is_zero(limdd_label(high));
    assert(!(low_zero && high_zero) && "the zero vector has no stabiliser group");
    assert((low_zero || limdd_lim_is_identity(limdd_label(low)))
           && "the low edge must be normalised before its group is taken");

    const limdd_pauli_t z_here = limdd_pauli_single(var, LIMDD_PAULI_Z);

    /* At most one generator per qubit, plus the anti-diagonal coset rep. */
    LIMDD_LIM out[LIMDD_MAX_QUBITS + 1];
    size_t nout = 0;

    if (high_zero) {
        /*
         * |0>(x)|v0>. Both I and Z fix |0>, so Z at this level joins Stab(v0)
         * untouched -- this is the one case that gains a generator rather than
         * losing some to the intersection.
         */
        out[nout++] = limdd_lim_make(z_here, EVBDD_ONE);
        for (LIMDD_STAB c = s0; c != LIMDD_STAB_TRIVIAL; c = stab_tail(c))
            out[nout++] = stab_head(c);
        return limdd_stab_make(out, nout);
    }

    const LIMDD_LIM b = limdd_label(high);
    const limdd_pauli_t q = limdd_lim_pauli(b);

    if (low_zero) {
        /*
         * |1>(x)B|v1>. Z fixes |1> only up to a sign, so it is -Z that belongs
         * here, and Stab(v1) arrives conjugated by B.
         */
        out[nout++] = limdd_lim_make(z_here, EVBDD_MIN_ONE);
        for (LIMDD_STAB c = s1; c != LIMDD_STAB_TRIVIAL; c = stab_tail(c))
            out[nout++] = conjugate(stab_head(c), q);
        return limdd_stab_make(out, nout);
    }

    /*
     * --- diagonal part: the subgroup intersection ---
     *
     * Skipped outright when either side is trivial, and that is exact rather
     * than a heuristic. The diagonal part is the set of elements of s0 whose
     * word lies in the span of s1. With s0 trivial there is nothing to draw
     * from; with s1 trivial its span is {0}, which admits only the identity
     * word, and no stabiliser group carries a non-identity element with the
     * identity word. Either way the intersection is empty.
     *
     * Worth a branch because it is the common case: on the circuits measured
     * 75% to 100% of nodes have a trivial group, so their parents take this
     * path and never enter the O(n^3) elimination. It also matches the bit
     * that limdd_makenode_ex records in the node, which is derived from the
     * same argument.
     */
    LIMDD_LIM common[LIMDD_MAX_QUBITS];
    const size_t ncommon = (limdd_stab_is_trivial(s0) || limdd_stab_is_trivial(s1))
                         ? 0 : intersect_gens(s0, s1, common);

    for (size_t i = 0; i < ncommon; i++) {
        const LIMDD_LIM a = common[i];
        const limdd_pauli_t wa = limdd_lim_pauli(a);

        const LIMDD_LIM a1 = stab_find(s1, wa);
        assert(a1 != 0 && "the intersection produced a word s1 does not have");

        /*
         * Conjugating by B flips A's sign when they anticommute, and the two
         * branches must then agree: Z here supplies the extra -1 when they do
         * not already match. So the level's Pauli is forced, never chosen.
         */
        const bool flip = !commutes(wa, q);
        const bool differ = (limdd_lim_weight(a) != limdd_lim_weight(a1));
        const bool need_z = (flip != differ);

        limdd_pauli_t w = wa;
        if (need_z) limdd_pauli_set(&w, var, LIMDD_PAULI_Z);
        out[nout++] = limdd_lim_make(w, limdd_lim_weight(a));
    }

    /* --- anti-diagonal coset representative --- */
    if (limdd_target(low) == limdd_target(high)) {
        const EVBDD_WGT beta = limdd_lim_weight(b);
        const EVBDD_WGT beta_sq = wgt_mul(beta, beta);

        /*
         * Swapping the branches sends |0>|v0> to |1>B|v1> and back, and the
         * round trip multiplies by B^2 times the Pauli's own i's. It closes
         * only when B's scalar squares to +-1: to +1 with X here, to -1 with Y,
         * whose extra factor of i is what absorbs the difference.
         */
        limdd_pauli_t w = q;
        EVBDD_WGT scalar;
        bool found = true;

        if (beta_sq == EVBDD_ONE) {
            limdd_pauli_set(&w, var, LIMDD_PAULI_X);
            scalar = beta;
        } else if (beta_sq == EVBDD_MIN_ONE) {
            limdd_pauli_set(&w, var, LIMDD_PAULI_Y);
            scalar = wgt_mul(limdd_wgt_i_pow(3), beta);   /* -i * beta */
        } else {
            found = false;
            scalar = EVBDD_ZERO;
        }

        if (found) {
            assert((scalar == EVBDD_ONE || scalar == EVBDD_MIN_ONE)
                   && "an anti-diagonal generator must have scalar +-1");
            out[nout++] = limdd_lim_make(w, scalar);
        }
    }

    return limdd_stab_make(out, nout);
}

void
limdd_gc_mark_stab(LIMDD_STAB s)
{
    assert(stab_table != NULL);
    for (; s != LIMDD_STAB_TRIVIAL; s = stab_tail(s)) {
        limdd_gc_mark_lim(stab_head(s));
        /* Lists share tails, so a cell already marked means the rest is too. */
        if (llmsset_mark(stab_table, s) == 0) return;
    }
}

void
limdd_gc_clear_stabs(void)
{
    llmsset_clear_data(stab_table);
}

void
limdd_gc_rehash_stabs(void)
{
    llmsset_clear_hashes(stab_table);
    if (llmsset_rehash(stab_table) != 0) {
        fprintf(stderr, "sylvan: LIMDD generator sets could not all be rehashed\n");
        exit(1);
    }
}

void
limdd_stab_init(size_t stab_tablesize)
{
    stab_table = llmsset_create(stab_tablesize, stab_tablesize);
    if (stab_table == NULL) {
        fprintf(stderr, "sylvan: could not create the LIMDD stabiliser table\n");
        exit(1);
    }
}

void
limdd_stab_quit(void)
{
    if (stab_table != NULL) {
        llmsset_free(stab_table);
        stab_table = NULL;
    }
}
