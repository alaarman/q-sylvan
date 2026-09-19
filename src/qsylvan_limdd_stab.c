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

static void
die(void)
{
    fprintf(stderr, "sylvan: LIMDD stabiliser table is full\n");
    exit(1);
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
limdd_stab_gen(LIMDD_STAB s, size_t i)
{
    while (i-- > 0) s = stab_tail(s);
    return stab_head(s);
}

/* --- reduction ----------------------------------------------------------- */

/**
 * Column `c` of the symplectic matrix: columns 0..n-1 are the X bits of qubits
 * 0..n-1, columns n..2n-1 the Z bits. Any fixed order gives a canonical RREF;
 * this is the usual (X|Z) one.
 */
static inline bool
row_bit(const gen_row_t *r, size_t c, size_t nqubits)
{
    const uint64_t w = (c < nqubits) ? r->p.x : r->p.z;
    const size_t k = (c < nqubits) ? c : c - nqubits;
    return (w >> k) & 1;
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
    for (size_t i = 0; s != LIMDD_STAB_TRIVIAL; i++, s = stab_tail(s)) {
        if ((k >> i) & 1) acc = limdd_lim_mul(acc, stab_head(s));
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
