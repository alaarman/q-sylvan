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
    bool                neighbours;     /* hybrid rule: probe the near cells too */
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

double
cmap_get_merge_tolerance()
{
    return TOL_MODE == CMAP_TOL_HYBRID ? (double) TOL_REL : (double) TOLERANCE;
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

/*
 * One probe sequence from `hash`: find a value close to `v`, or, if `insert`,
 * put `val` at the first empty bucket. Returns 1 if found, 0 if inserted, 2 if
 * not found (an empty bucket ends the sequence, since nothing is ever removed
 * between collections), and -1 if the table is full.
 */
static int
probe_sequence(cmap_t *cmap, uint32_t hash, const bucket_t *val, const complex_t *v,
               bool insert, uint64_t *ret)
{
    const uint32_t prime = odd_primes[hash & PRIME_MASK];
    for (unsigned int c = 0; c < cmap->threshold; c++) {
        uint64_t            ref = hash & cmap->mask;
        uint64_t            line_end = (ref & CL_MASK) + CACHE_LINE_SIZE;
        for (size_t i = 0; i < CACHE_LINE_SIZE; i++) {

            // 1. Get bucket
            bucket_t *bucket = &cmap->table[ref];

            // 2. If bucket empty, insert new value here
            if (bucket->d[0] == EMPTY) {
                if (!insert) return 2;
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

static inline uint32_t
hash_cell(fl_t r, fl_t i)
{
    bucket_t round_v;
    round_v.c.r = r;
    round_v.c.i = i;
    // fix 0 possibly having a sign
    if(round_v.c.r == 0.0) round_v.c.r = 0.0;
    if(round_v.c.i == 0.0) round_v.c.i = 0.0;
    return SuperFastHash(&round_v, sizeof(complex_t), 0);
}

int
cmap_find_or_put(const void *dbs, const void *_v, uint64_t *ret)
{
    complex_t *v = (complex_t *) _v;
    cmap_t *cmap = (cmap_t *) dbs;
    bucket_t *val  = (bucket_t *) v;

    assert (val->d[0] != LOCK);
    assert (val->d[0] != EMPTY);

    // Round the value to compute the hash with, but store the actual value v
    if (TOL_MODE == CMAP_TOL_HYBRID) {
        long double m = cmag(v);
        if (m <= TOL_ZERO || TOL_REL == 0.0l) {
            return probe_sequence(cmap, hash_cell((m <= TOL_ZERO) ? 0.0 : v->r,
                                                  (m <= TOL_ZERO) ? 0.0 : v->i),
                                  val, v, true, ret);
        }
        /*
         * Log-polar quantisation, in cells of relative width 2*TOL_REL in
         * magnitude and 2*TOL_REL in angle. Two weights within TOL_REL of
         * each other are less than half a cell apart in both coordinates, so
         * the other one lies in this weight's own cell or in the neighbour
         * on the side of the nearer edge -- four cells in all. A lookup that
         * misses its own cell asks those neighbours before inserting:
         * without that, two weights that the rule calls equal but that fall
         * either side of a cell edge were never merged, and a LIMDD's
         * canonical form, which compares interned labels by index, kept two
         * nodes for one class (grover_iter_n8_it12 under L2 at 1e-9: labels
         * 3e-10 apart, interned twice).
         *
         * Two limits. The probes are not atomic with the insert, so two
         * workers inserting close values in neighbouring cells at once can
         * both insert; only a one-cell race is excluded, by the CAS. And the
         * cell coordinates are computed in long double, which is double on
         * arm64: past |log m| of about 45 their rounding is a sizeable part
         * of a cell at 1e-14, and close pairs can again be missed.
         */
        const long double w = 2.0l * TOL_REL;
        const long double fr = logl(m) / log1pl(w);
        const long double kr = floorl(fr + 0.5l);
        // angular cells of width w, wrapped so that -pi and +pi agree
        const long double nbuckets = floorl(2.0l * 3.14159265358979323846264338327950288l / w) + 1.0l;
        const long double fa = atan2l((long double) v->i, (long double) v->r) / w;
        const long double ka = floorl(fa + 0.5l);
        const long double dr = (fr >= kr) ? 1.0l : -1.0l;
        const long double da = (fa >= ka) ? 1.0l : -1.0l;
        long double cr[4] = { kr, kr + dr, kr,      kr + dr };
        long double ca[4] = { ka, ka,      ka + da, ka + da };
        uint32_t h[4];
        for (int k = 0; k < 4; k++) {
            long double a = fmodl(ca[k], nbuckets);
            if (a < 0.0l) a += nbuckets;
            h[k] = hash_cell((fl_t) cr[k], (fl_t) a);
        }
        if (!cmap->neighbours) return probe_sequence(cmap, h[0], val, v, true, ret);
        for (int k = 0; k < 4; k++) {
            const int r = probe_sequence(cmap, h[k], val, v, false, ret);
            if (r != 2) return r;          // found, or full
        }
        return probe_sequence(cmap, h[0], val, v, true, ret);
    }
    else if (TOLERANCE == 0.0) {
        return probe_sequence(cmap, hash_cell(v->r, v->i), val, v, true, ret);
    }
    else {
        return probe_sequence(cmap, hash_cell(flt_round(v->r / TOLERANCE) * TOLERANCE,
                                              flt_round(v->i / TOLERANCE) * TOLERANCE),
                              val, v, true, ret);
    }
}

void
cmap_set_neighbour_probe(void *dbs, bool on)
{
    ((cmap_t *) dbs)->neighbours = on;
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
    for (uint64_t c = 0; c < cmap->size; c++) {
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
    for (uint64_t c = 0; c < cmap->size; c++) {
        cmap->table[c].d[0] = EMPTY;
    }
    cmap->threshold = cmap->size / 100;
    cmap->threshold = min(cmap->threshold, 1ULL << 16);
    cmap->seen_0 = 0;
    cmap->neighbours = true;
    return (void *) cmap;
}

void
cmap_free(void *dbs)
{
    cmap_t * cmap = (cmap_t *) dbs;
    free (cmap->table);
    free (cmap);
}
