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

#include "qsylvan_limdd_lim.h"
#include "sylvan_edge_weights_complex.h"

LIMDD_LIM LIMDD_LIM_IDENTITY = 0;
LIMDD_LIM LIMDD_LIM_ZERO     = 0;

/*
 * Written once by limdd_lims_init and only read afterwards. Read-only shared
 * data cannot be false-shared: a line nobody writes is never invalidated.
 */
static llmsset_t pauli_table = NULL;
static llmsset_t lim_table   = NULL;
static size_t    lim_nqubits = 0;

/**
 * i^0, i^1, i^2, i^3 as edge weights.
 *
 * Multiplying two Pauli words produces a factor of i^k, and that factor is
 * folded into the scalar rather than stored alongside the word. Interning the
 * four powers once turns that fold into a single table multiply.
 */
static EVBDD_WGT i_pow[4];

static void
die(const char *what)
{
    fprintf(stderr, "sylvan: LIMDD %s table is full\n", what);
    exit(1);
}

size_t
limdd_lims_nqubits(void)
{
    return lim_nqubits;
}

LIMDD_PAULI_REF
limdd_pauli_intern(limdd_pauli_t p)
{
    assert(pauli_table != NULL);
    assert(limdd_pauli_is_canonical(p, lim_nqubits));

    int created;
    uint64_t ref = llmsset_lookup(pauli_table, p.x, p.z, &created);
    if (ref == 0) die("Pauli");
    return ref;
}

limdd_pauli_t
limdd_pauli_deref(LIMDD_PAULI_REF ref)
{
    assert(pauli_table != NULL);
    const uint64_t *bucket = (const uint64_t *)llmsset_index_to_ptr(pauli_table, ref);
    limdd_pauli_t p = { bucket[0], bucket[1] };
    return p;
}

/**
 * llmsset counts occupied buckets, and buckets 0 and 1 are reserved and marked
 * occupied when the table is created, so they have to come off the total.
 */
static size_t
table_count(llmsset_t dbs)
{
    if (dbs == NULL) return 0;
    const size_t occupied = (size_t)llmsset_count_marked(dbs);
    return occupied < 2 ? 0 : occupied - 2;
}

size_t
limdd_pauli_table_count(void)
{
    return table_count(pauli_table);
}

size_t
limdd_lim_table_count(void)
{
    return table_count(lim_table);
}

LIMDD_LIM
limdd_lim_make(limdd_pauli_t p, EVBDD_WGT w)
{
    assert(lim_table != NULL);

    /*
     * The zero map does not depend on the Pauli word, so every zero-scalar LIM
     * must land on the same index. Skipping this would hand out a different
     * index per Pauli word for one and the same map, and node lookup would
     * then fail to find nodes that already exist.
     */
    if (w == EVBDD_ZERO) return LIMDD_LIM_ZERO;

    const LIMDD_PAULI_REF pref = limdd_pauli_intern(p);

    int created;
    uint64_t lim = llmsset_lookup(lim_table, pref, w, &created);
    if (lim == 0) die("LIM");
    return lim;
}

limdd_pauli_t
limdd_lim_pauli(LIMDD_LIM lim)
{
    assert(lim_table != NULL);
    const uint64_t *bucket = (const uint64_t *)llmsset_index_to_ptr(lim_table, lim);
    return limdd_pauli_deref(bucket[0]);
}

EVBDD_WGT
limdd_lim_weight(LIMDD_LIM lim)
{
    assert(lim_table != NULL);
    const uint64_t *bucket = (const uint64_t *)llmsset_index_to_ptr(lim_table, lim);
    return (EVBDD_WGT)bucket[1];
}

LIMDD_LIM
limdd_lim_mul(LIMDD_LIM a, LIMDD_LIM b)
{
    if (limdd_lim_is_zero(a) || limdd_lim_is_zero(b)) return LIMDD_LIM_ZERO;
    if (limdd_lim_is_identity(a)) return b;
    if (limdd_lim_is_identity(b)) return a;

    limdd_pauli_t pa = limdd_lim_pauli(a);
    const limdd_pauli_t pb = limdd_lim_pauli(b);

    // pa becomes the product; the returned k is the power of i left over
    const unsigned k = limdd_pauli_rightmul(&pa, pb, lim_nqubits);

    EVBDD_WGT w = wgt_mul(limdd_lim_weight(a), limdd_lim_weight(b));
    if (k != 0) w = wgt_mul(w, i_pow[k]);

    return limdd_lim_make(pa, w);
}

LIMDD_LIM
limdd_lim_inverse(LIMDD_LIM lim)
{
    if (limdd_lim_is_zero(lim)) return LIMDD_LIM_ZERO;
    if (limdd_lim_is_identity(lim)) return LIMDD_LIM_IDENTITY;

    // A Pauli word squares to the identity with no leftover phase, so the word
    // is its own inverse and only the scalar has to be inverted.
    const EVBDD_WGT w = wgt_div(EVBDD_ONE, limdd_lim_weight(lim));
    return limdd_lim_make(limdd_lim_pauli(lim), w);
}

void
limdd_lim_fprint(FILE *out, LIMDD_LIM lim)
{
    char buf[LIMDD_MAX_QUBITS + 1];
    limdd_pauli_to_string(limdd_lim_pauli(lim), lim_nqubits, buf);
    wgt_fprint(out, limdd_lim_weight(lim));
    fprintf(out, " * %s", buf);
}

void
limdd_lims_init(size_t nqubits, size_t pauli_tablesize, size_t lim_tablesize)
{
    if (nqubits > LIMDD_MAX_QUBITS) {
        fprintf(stderr, "sylvan: LIMDD supports at most %d qubits, got %zu\n",
                LIMDD_MAX_QUBITS, nqubits);
        exit(1);
    }

    lim_nqubits  = nqubits;
    pauli_table  = llmsset_create(pauli_tablesize, pauli_tablesize);
    lim_table    = llmsset_create(lim_tablesize, lim_tablesize);

    /*
     * The four powers of i. A Pauli product can produce any of them, so all
     * four must exist in the weight table; that in turn requires the weights
     * to be complex.
     */
    i_pow[0] = EVBDD_ONE;
    i_pow[1] = complex_lookup(0.0, 1.0);
    i_pow[2] = EVBDD_MIN_ONE;
    i_pow[3] = complex_lookup(0.0, -1.0);

    if (i_pow[1] == i_pow[0] || i_pow[1] == i_pow[2]) {
        fprintf(stderr, "sylvan: LIMDD needs complex edge weights; "
                        "i is not representable in the current weight type\n");
        exit(1);
    }

    LIMDD_LIM_IDENTITY = limdd_lim_make(limdd_pauli_identity(), EVBDD_ONE);

    /*
     * LIMDD_LIM_ZERO must be interned by hand rather than through
     * limdd_lim_make, because that function maps every zero-scalar LIM onto
     * this value and so cannot be used to create it.
     */
    int created;
    const LIMDD_PAULI_REF id_ref = limdd_pauli_intern(limdd_pauli_identity());
    LIMDD_LIM_ZERO = llmsset_lookup(lim_table, id_ref, EVBDD_ZERO, &created);
    if (LIMDD_LIM_ZERO == 0) die("LIM");
}

void
limdd_lims_quit(void)
{
    if (lim_table != NULL) {
        llmsset_free(lim_table);
        lim_table = NULL;
    }
    if (pauli_table != NULL) {
        llmsset_free(pauli_table);
        pauli_table = NULL;
    }
    lim_nqubits = 0;
    LIMDD_LIM_IDENTITY = 0;
    LIMDD_LIM_ZERO = 0;
}
