#include "cmap.h"

#include <assert.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "atomics.h"
#include "cmap.h"
#include "fast_hash.h"
#include "util.h"

#undef CACHE_LINE
#undef CACHE_LINE_SIZE

#define CACHE_LINE 8
#define CACHE_LINE_SIZE 256

// how many "blocks" of 64 bits for a single table entry
#define entry_size (2*sizeof(fl_t)/8)

typedef union {
    complex_t       c;
    uint64_t        d[entry_size];
} bucket_t;

// float "equality" tolerance
static long double TOLERANCE = 1e-14l;

/* Merging rule.
 * CMAP_TOL_ABS    (default, historical): two weights are equal when their real
 *                 and imaginary parts each differ by less than TOLERANCE.
 * CMAP_TOL_HYBRID (see Brand et al., "Numerical Errors ... With Edge-Weighted
 *                 Decision Diagrams"): two non-zero weights are equal when they
 *                 differ by less than TOL_REL *relative* to their magnitude,
 *                 and a weight of magnitude at most TOL_ZERO is collapsed to
 *                 exactly 0.  Relative merging bounds the merging error by the
 *                 flow through a weight rather than by its sensitivity, which
 *                 makes it insensitive to the magnitude of the weight and hence
 *                 to the normalisation strategy.
 */
static int         TOL_MODE = CMAP_TOL_ABS;
static long double TOL_REL  = 1e-14l;   // relative threshold  (delta_rel)
static long double TOL_ZERO = 0.0l;     // zero-collapse threshold (delta_0)

void
cmap_set_hybrid_tolerance(double rel, double zero)
{
    TOL_MODE = CMAP_TOL_HYBRID;
    TOL_REL  = (long double) rel;
    TOL_ZERO = (long double) zero;
}

void
cmap_set_absolute_tolerance(double tol)
{
    TOL_MODE  = CMAP_TOL_ABS;
    TOLERANCE = (long double) tol;
}

int
cmap_get_tolerance_mode()
{
    return TOL_MODE;
}

static inline long double
cmag(const complex_t *v)
{
    return hypotl((long double) v->r, (long double) v->i);
}
static const uint64_t EMPTY = 14738995463583502973ull;
static const uint64_t LOCK  = 14738995463583502974ull;
static const uint64_t CL_MASK = -(1ULL << CACHE_LINE);

/**
\typedef Lockless hastable database.
*/
typedef struct cmap_s cmap_t;
struct cmap_s {
    size_t              size;
    size_t              mask;
    size_t              threshold;
    int                 seen_0;
    bucket_t  __attribute__(( __aligned__(32)))       *table;
    // Q: should this 32 change to 16 now that we use doubles instead of
    // long doubles for the real and imaginary components?
};

static void __attribute__((unused))
print_bucket_floats(bucket_t *b)
{
    printf("%.60Lf, %.60Lf\n", (long double) b->c.r, (long double) b->c.i);
}

static void __attribute__((unused))
print_bucket_bits(bucket_t* b)
{
    printf("%016" PRIu64, b->d[0]);
    for (unsigned int k = 1; k < entry_size; k++) {
        printf(" %016" PRIu64, b->d[k]);
    }
    printf("\n");
}

double
cmap_get_tolerance()
{
    return TOLERANCE;
}

static bool
complex_close(complex_t *in_table, const complex_t* to_insert)
{
    if (TOL_MODE == CMAP_TOL_HYBRID) {
        long double ma = cmag(in_table), mb = cmag(to_insert);
        // zero-collapse: a weight of magnitude <= TOL_ZERO *is* zero
        bool za = (ma <= TOL_ZERO), zb = (mb <= TOL_ZERO);
        if (za || zb) return (za && zb);
        if (TOL_REL == 0.0l) {
            return ((in_table->r == to_insert->r) &&
                    (in_table->i == to_insert->i));
        }
        // relative: |a - b| <= delta_rel * max(|a|,|b|)   (symmetric)
        long double dr = (long double) in_table->r - (long double) to_insert->r;
        long double di = (long double) in_table->i - (long double) to_insert->i;
        long double m  = (ma > mb) ? ma : mb;
        return (hypotl(dr, di) <= TOL_REL * m);
    }

    if (TOLERANCE == 0.0) {
         return ((in_table->r == to_insert->r) && 
                 (in_table->i == to_insert->i));
    }
    else {
        return ((flt_abs(in_table->r - to_insert->r) < TOLERANCE) && 
                (flt_abs(in_table->i - to_insert->i) < TOLERANCE));
    }
    
}

int
cmap_find_or_put(const void *dbs, const void *_v, uint64_t *ret)
{
    complex_t *v = (complex_t *) _v;
    cmap_t *cmap = (cmap_t *) dbs;
    bucket_t *val  = (bucket_t *) v;

    // Round the value to compute the hash with, but store the actual value v
    bucket_t round_v;
    if (TOL_MODE == CMAP_TOL_HYBRID) {
        // Log-polar quantisation: cells of *relative* width TOL_REL, so that
        // relatively-close weights land in the same (or an adjacent) cell,
        // exactly as absolutely-close weights do on the uniform grid below.
        long double m = cmag(v);
        if (m <= TOL_ZERO || TOL_REL == 0.0l) {
            round_v.c.r = (m <= TOL_ZERO) ? 0.0 : v->r;
            round_v.c.i = (m <= TOL_ZERO) ? 0.0 : v->i;
        }
        else {
            long double kr = floorl(logl(m) / log1pl(TOL_REL) + 0.5l);
            // angular buckets of width TOL_REL, wrapped so that -pi and +pi agree
            long double nbuckets = floorl(2.0l * 3.14159265358979323846264338327950288l / TOL_REL) + 1.0l;
            long double ka = floorl(atan2l((long double) v->i,
                                           (long double) v->r) / TOL_REL + 0.5l);
            ka = fmodl(ka, nbuckets);
            if (ka < 0.0l) ka += nbuckets;
            round_v.c.r = (fl_t) kr;
            round_v.c.i = (fl_t) ka;
        }
    }
    else if (TOLERANCE == 0.0) {
        round_v.c.r = v->r;
        round_v.c.i = v->i;
    }
    else {
        round_v.c.r = flt_round(v->r / TOLERANCE) * TOLERANCE;
        round_v.c.i = flt_round(v->i / TOLERANCE) * TOLERANCE;
    }

    // fix 0 possibly having a sign
    if(round_v.c.r == 0.0) round_v.c.r = 0.0;
    if(round_v.c.i == 0.0) round_v.c.i = 0.0;

    //printf("(%.3f,%.3f) ",(float)round_v.c.r,(float)round_v.c.i);
    //print_bucket_bits(&round_v); 
    
    uint32_t hash  = SuperFastHash(&round_v, sizeof(complex_t), 0);
    uint32_t prime = odd_primes[hash & PRIME_MASK];

    assert (val->d[0] != LOCK);
    assert (val->d[0] != EMPTY);

    // Insert/lookup `v`
    for (unsigned int c = 0; c < cmap->threshold; c++) {
        uint64_t            ref = hash & cmap->mask;
        uint64_t            line_end = (ref & CL_MASK) + CACHE_LINE_SIZE;
        for (size_t i = 0; i < CACHE_LINE_SIZE; i++) {
            
            // 1. Get bucket
            bucket_t *bucket = &cmap->table[ref];

            // 2. If bucket empty, insert new value here
            if (bucket->d[0] == EMPTY) {
                if (cas(&bucket->d[0], EMPTY, LOCK)) {
                    *ret = ref;
                    // write backwards (overwrite bucket->d[0] last)
                    for (int k = entry_size-1; k >= 0; k--) {
                        atomic_write (&bucket->d[k], val->d[k]);
                    }
                    return 0;
                }
            }

            // 3. Bucket not empty, wait for lock
            while (atomic_read(&bucket->d[0]) == LOCK) {}

            // 4. Bucket contains some complex value, check if close to `v`
            complex_t *in_table = (complex_t *)bucket;
            if (complex_close(in_table, v)) {
                *ret = ref;
                return 1;
            }

            // If unsuccessful, try next
            ref += 1;
            ref = ref == line_end ? line_end - CACHE_LINE_SIZE : ref;
        }
        hash += prime << CACHE_LINE;
    }
    // amplitude table full, unable to add
    return -1;
}

void *
cmap_get(const void *dbs, const uint64_t ref)
{
    cmap_t *cmap = (cmap_t *) dbs;
    return &(cmap->table[ref].c);
}

uint64_t
cmap_count_entries(const void *dbs)
{
    cmap_t *cmap = (cmap_t *) dbs;
    uint64_t entries = 0;
    for (unsigned int c = 0; c < cmap->size; c++) {
        if (cmap->table[c].d[0] != EMPTY)
            entries++;
    }
    return entries;
}

void
print_bitvalues(const void *dbs, const uint64_t ref)
{
    cmap_t *cmap = (cmap_t *) dbs;
    bucket_t* b = cmap_get(cmap, ref);
    printf("%016" PRIu64, b->d[0]);
    for (unsigned int k = 1; k < entry_size; k++) {
        printf(" %016" PRIu64, b->d[k]);
    }
}

void *
cmap_create(uint64_t size, double tolerance)
{
    TOLERANCE = tolerance;
    cmap_t  *cmap = calloc (1, sizeof(cmap_t));
    cmap->size = size;
    cmap->mask = cmap->size - 1;
    cmap->table = calloc (cmap->size, sizeof(bucket_t));
    for (unsigned int c = 0; c < cmap->size; c++) {
        cmap->table[c].d[0] = EMPTY;
    }
    cmap->threshold = cmap->size / 100;
    cmap->threshold = min(cmap->threshold, 1ULL << 16);
    cmap->seen_0 = 0;
    return (void *) cmap;
}

void
cmap_free(void *dbs)
{
    cmap_t * cmap = (cmap_t *) dbs;
    free (cmap->table);
    free (cmap);
}
