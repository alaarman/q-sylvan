#include "qisq2_map.h"

#include <assert.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "atomics.h"
#include "fast_hash.h"
#include "util.h"
#include "../sylvan_edge_weights_qisq2.h"

#undef CACHE_LINE
#undef CACHE_LINE_SIZE

#define CACHE_LINE 8
#define CACHE_LINE_SIZE 256

// how many "blocks" of 64 bits for a single table entry
#define entry_size (4*sizeof(mpq_t)/8)

typedef union {
    qisq2_t         c;
    uint64_t        d[entry_size];
} bucket_t;

// float "equality" tolerance
static long double TOLERANCE = 1e-14l;
/*
 * A bucket's first word is the (_mp_alloc, _mp_size) pair of its real part's
 * numerator. EMPTY is 0, so that the table can come from calloc and its pages
 * stay untouched until a weight lands in them: with a non-zero marker written
 * into every 128-byte bucket, a 2^24 table was 2 GB resident before the first
 * gate, whatever the run then used. A stored numerator always has storage
 * (qisq2_map_insert sees to it, since GMP leaves a fresh zero unallocated), so
 * its word is never 0; and an allocation of -1 never occurs, so neither is it
 * LOCK.
 */
static const uint64_t EMPTY = 0;
static const uint64_t LOCK  = UINT64_MAX;
static const uint64_t CL_MASK = -(1ULL << CACHE_LINE);

/**
\typedef Lockless hastable database.
*/
typedef struct qisq2_map_s qisq2_map_t;
struct qisq2_map_s {
    size_t              size;
    size_t              mask;
    size_t              threshold;
    int                 seen_0;
    bucket_t  __attribute__(( __aligned__(32)))       *table;
    // Q: should this 32 change to 16 now that we use doubles instead of
    // long doubles for the real and imaginary components?
};

static void __attribute__((unused))
print_bucket_bits(bucket_t* b)
{
    printf("%016" PRIu64, b->d[0]);
    for (unsigned int k = 1; k < entry_size; k++) {
        printf(" %016" PRIu64, b->d[k]);
    }
    printf("\n");
}

/*
Hash the number num
This function assumes that num is reduced already

It reads the limbs where they are: copying the eight numerators and
denominators out first cost eight allocations per lookup. 0 and 1 still hash to
0 and 1, so that they take buckets 0 and 1 of a fresh table.
*/
static inline uint64_t
qisq2_mix(uint64_t h)
{
    h ^= h >> 33; h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33; h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;
    return h;
}

static inline uint64_t
qisq2_mpz_hash(mpz_srcptr z, uint64_t h)
{
    const mp_limb_t *l = mpz_limbs_read(z);
    const size_t n = mpz_size(z);
    h = qisq2_mix(h ^ (uint64_t)(int64_t)z->_mp_size);
    for (size_t i = 0; i < n; i++) h = qisq2_mix(h ^ (uint64_t)l[i]);
    return h;
}

uint32_t
qisq2_hash(qisq2_t *num) {
    const bool a0 = mpq_sgn(num->a) == 0, b0 = mpq_sgn(num->b) == 0,
               c0 = mpq_sgn(num->c) == 0, d0 = mpq_sgn(num->d) == 0;
    if (a0 && b0 && c0 && d0) return 0;
    if (b0 && c0 && d0 && mpq_cmp_ui(num->a, 1, 1) == 0) return 1;
    uint64_t h = 0x9e3779b97f4a7c15ULL;
    h = qisq2_mpz_hash(mpq_numref(num->a), h); h = qisq2_mpz_hash(mpq_denref(num->a), h);
    h = qisq2_mpz_hash(mpq_numref(num->b), h); h = qisq2_mpz_hash(mpq_denref(num->b), h);
    h = qisq2_mpz_hash(mpq_numref(num->c), h); h = qisq2_mpz_hash(mpq_denref(num->c), h);
    h = qisq2_mpz_hash(mpq_numref(num->d), h); h = qisq2_mpz_hash(mpq_denref(num->d), h);
    return (uint32_t)(h ^ (h >> 32));
}

void qisq2_reduce(qisq2_t *num) {
    mpq_canonicalize(num->a);
    mpq_canonicalize(num->b);
    mpq_canonicalize(num->c);
    mpq_canonicalize(num->d);
}

/* Same sign and length in all eight integers: read from the entry itself, so a
 * probe that meets another value rarely has to follow its limb pointers. */
static inline bool
qisq2_same_shape(const qisq2_t *x, const qisq2_t *y)
{
    return mpq_numref(x->a)->_mp_size == mpq_numref(y->a)->_mp_size &&
           mpq_denref(x->a)->_mp_size == mpq_denref(y->a)->_mp_size &&
           mpq_numref(x->b)->_mp_size == mpq_numref(y->b)->_mp_size &&
           mpq_denref(x->b)->_mp_size == mpq_denref(y->b)->_mp_size &&
           mpq_numref(x->c)->_mp_size == mpq_numref(y->c)->_mp_size &&
           mpq_denref(x->c)->_mp_size == mpq_denref(y->c)->_mp_size &&
           mpq_numref(x->d)->_mp_size == mpq_numref(y->d)->_mp_size &&
           mpq_denref(x->d)->_mp_size == mpq_denref(y->d)->_mp_size;
}

static bool
qisq2_equal(qisq2_t *in_table, const qisq2_t* to_insert)
{
    if (!qisq2_same_shape(in_table, to_insert)) return false;
    if (mpq_equal(in_table->a,to_insert->a) == 0){
        return false;
    }
    if (mpq_equal(in_table->b,to_insert->b) == 0){
        return false;
    }
    if (mpq_equal(in_table->c,to_insert->c) == 0){
        return false;
    }
    if (mpq_equal(in_table->d,to_insert->d) == 0){
        return false;
    }
    return true;
}


void qisq2_clear(qisq2_t *num) {
    mpq_clear(num->a);
    mpq_clear(num->b);
    mpq_clear(num->c);
    mpq_clear(num->d);
    //qisq2_reduce(num);
}


static int qisq2_map_insert(qisq2_map_t *qisq2_map, qisq2_t *v, uint64_t *ret);

int
qisq2_map_find_or_put(const void *dbs, const void *_v, uint64_t *ret)
{
    qisq2_t *v = (qisq2_t *)_v;
    qisq2_reduce(v);
    return qisq2_map_insert((qisq2_map_t *) dbs, v, ret);
}

/*
 * The same for a value that is already reduced, as every result of GMP's mpq
 * arithmetic is: reducing it again costs a gcd per component and changes
 * nothing. `v` is only read.
 */
int
qisq2_map_find_or_put_reduced(const void *dbs, const qisq2_t *v, uint64_t *ret)
{
#ifndef NDEBUG
    mpq_t chk;
    mpq_init(chk);
    mpq_set(chk, v->a); mpq_canonicalize(chk); assert(mpq_equal(chk, v->a));
    mpq_set(chk, v->b); mpq_canonicalize(chk); assert(mpq_equal(chk, v->b));
    mpq_set(chk, v->c); mpq_canonicalize(chk); assert(mpq_equal(chk, v->c));
    mpq_set(chk, v->d); mpq_canonicalize(chk); assert(mpq_equal(chk, v->d));
    mpq_clear(chk);
#endif
    return qisq2_map_insert((qisq2_map_t *) dbs, (qisq2_t *) v, ret);
}

static int
qisq2_map_insert(qisq2_map_t *qisq2_map, qisq2_t *v, uint64_t *ret)
{
    bucket_t *val  = (bucket_t *) v;

    uint32_t hash = qisq2_hash(v);
    uint32_t prime = odd_primes[hash & PRIME_MASK];

    (void)val;

    // Insert/lookup `v`
    for (unsigned int c = 0; c < qisq2_map->threshold; c++) {
        uint64_t            ref = hash & qisq2_map->mask;
        uint64_t            line_end = (ref & CL_MASK) + CACHE_LINE_SIZE;
        for (size_t i = 0; i < CACHE_LINE_SIZE; i++) {
            
            // 1. Get bucket
            bucket_t *bucket = &qisq2_map->table[ref];

            // 2. If bucket empty, insert new value here
            if (bucket->d[0] == EMPTY) {
                if (cas(&bucket->d[0], EMPTY, LOCK)) {
                    *ret = ref;
                    /*
                     * The table gets its OWN limbs. Publishing `v`'s words
                     * directly would hand the caller's GMP storage to the
                     * table while the caller still holds it, so whichever of
                     * the two cleared first would corrupt the other. With a
                     * copy, find_or_put never consumes its argument and the
                     * caller can always release it. The copy is made here,
                     * after winning the CAS, so a thread that loses does not
                     * make one.
                     */
                    qisq2_t owned;
                    qisq2_init(&owned);
                    mpq_set(owned.a, v->a);
                    mpq_set(owned.b, v->b);
                    mpq_set(owned.c, v->c);
                    mpq_set(owned.d, v->d);
                    /* A numerator of 0 may have no storage, which would make
                     * the bucket's first word 0, i.e. EMPTY. */
                    if (mpq_numref(owned.a)->_mp_alloc == 0)
                        mpz_realloc2(mpq_numref(owned.a), 64);

                    bucket_t *src = (bucket_t *)&owned;
                    assert(src->d[0] != EMPTY && src->d[0] != LOCK);
                    // write backwards (overwrite bucket->d[0] last)
                    for (int k = entry_size-1; k >= 0; k--) {
                        atomic_write (&bucket->d[k], src->d[k]);
                    }
                    return 0;
                }
            }

            // 3. Bucket not empty, wait for lock
            while (atomic_read(&bucket->d[0]) == LOCK) {}

            // 4. Bucket contains some complex value, check if equal to `v`
            qisq2_t *in_table = (qisq2_t *)bucket;
            if (qisq2_equal(in_table, v)) {
                *ret = ref;
                /* `v` belongs to the caller; clearing it here is what made a
                 * lookup consume its argument. */
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
qisq2_map_get(const void *dbs, const uint64_t ref)
{
    qisq2_map_t *qisq2_map = (qisq2_map_t *) dbs;
    return &(qisq2_map->table[ref].c);
}

uint64_t
qisq2_map_count_entries(const void *dbs)
{
    qisq2_map_t *qisq2_map = (qisq2_map_t *) dbs;
    uint64_t entries = 0;
    for (uint64_t c = 0; c < qisq2_map->size; c++) {
        if (qisq2_map->table[c].d[0] != EMPTY)
            entries++;
    }
    return entries;
}

void
qisq2_map_print_bitvalues(const void *dbs, const uint64_t ref)
{
    qisq2_map_t *qisq2_map = (qisq2_map_t *) dbs;
    bucket_t *b = (bucket_t *) qisq2_map_get(qisq2_map, ref);
    printf("%016" PRIu64, b->d[0]);
    for (unsigned int k = 1; k < entry_size; k++) {
        printf(" %016" PRIu64, b->d[k]);
    }
}

void *
qisq2_map_create(uint64_t size, double tolerance)
{
    qisq2_map_t  *qisq2_map = calloc (1, sizeof(qisq2_map_t));
    qisq2_map->size = size;
    qisq2_map->mask = qisq2_map->size - 1;
    /* calloc, and EMPTY is 0: pages are only touched when used. */
    qisq2_map->table = calloc (qisq2_map->size, sizeof(bucket_t));
    qisq2_map->threshold = qisq2_map->size / 100;
    qisq2_map->threshold = min(qisq2_map->threshold, 1ULL << 16);
    qisq2_map->seen_0 = 0;
    return (void *) qisq2_map;
}

void 
qisq2_map_free_elt(qisq2_t *num){
    qisq2_clear(num);
}

void
qisq2_print(qisq2_t *x)
{
    printf("\t%f + %f sqrt(2) + %f i + %f i sqrt(2)\n",
               mpq_get_d(x->a), mpq_get_d(x->b), mpq_get_d(x->c), mpq_get_d(x->d));
    gmp_printf("\t%#4Qx + \n\t%#4Qx sqrt(2) +\n\t%#4Qx i +\n\t%#4Qx i sqrt(2)\n\n", x->a, x->b, x->c, x->d);
}

void
qisq2_map_free(void *dbs)
{
    qisq2_map_t * qisq2_map = (qisq2_map_t *) dbs;
    for (uint64_t c = 0; c < qisq2_map->size; c++) {
        if (qisq2_map->table[c].d[0] != EMPTY){
            /*
            // sanity check for debugging
            if (c!=(qisq2_hash(&(qisq2_map->table[c].c))& qisq2_map->mask)){
                uint32_t hash = qisq2_hash(&(qisq2_map->table[c].c));
                uint64_t ref = hash & qisq2_map->mask;
                uint64_t line_end = (ref & CL_MASK) + CACHE_LINE_SIZE;
                printf("ERROR: hash table failure \n");
                uint64_t d = (uint64_t)c;
                printf("Hash value of entry is: %ld \n", d);
                printf("Hash value of entry should be: %ld \n", qisq2_hash(&(qisq2_map->table[c].c))& qisq2_map->mask);
                printf("or equivalently: %d \n", hash);
                qisq2_print(&(qisq2_map->table[c].c));
                printf("line_end: %ld \n",line_end);
                uint32_t prime = odd_primes[hash & PRIME_MASK];
                printf("prime: %d \n",prime);
                hash += prime << CACHE_LINE;
                printf("new hash: %d \n", hash);
                ref = hash & qisq2_map->mask;
                printf("new ref: %ld \n", ref);
                //exit(0);
            }
            */
            qisq2_map_free_elt(&(qisq2_map->table[c].c));
        }
    }
    free (qisq2_map->table);
    free (qisq2_map);
}

uint64_t
qisq2_size(const void *dbs, const uint64_t ref)
{
    // return biggest number of bits used for some integer in qisq number ref
    //
    // More concrete: for qisq number anum/aden + ... + dnum/dden i sqrt(2)
    // return maxsize(anum,aden,bnum,...,dden) 
    qisq2_map_t * qisq2_map = (qisq2_map_t *) dbs;
    qisq2_t *num = &(qisq2_map->table[ref].c);

    mpz_t anum, aden, bnum, bden, cnum, cden, dnum, dden;
    mpz_init(anum); mpz_init(aden); 
    mpz_init(bnum); mpz_init(bden);
    mpz_init(cnum); mpz_init(cden); 
    mpz_init(dnum); mpz_init(dden);
    uint64_t max_size = 0;

    mpq_get_num(anum, num->a);
    mpq_get_den(aden, num->a);
    mpq_get_num(bnum, num->b);
    mpq_get_den(bden, num->b);
    mpq_get_num(cnum, num->c);
    mpq_get_den(cden, num->c);
    mpq_get_num(dnum, num->d);
    mpq_get_den(dden, num->d);
    max_size = max(mpz_sizeinbase(anum,2),max_size);
    max_size = max(mpz_sizeinbase(aden,2),max_size);
    max_size = max(mpz_sizeinbase(bnum,2),max_size);
    max_size = max(mpz_sizeinbase(bden,2),max_size);
    max_size = max(mpz_sizeinbase(cnum,2),max_size);
    max_size = max(mpz_sizeinbase(cden,2),max_size);
    max_size = max(mpz_sizeinbase(bnum,2),max_size);
    max_size = max(mpz_sizeinbase(dden,2),max_size);

    mpz_clear(anum); mpz_clear(aden); 
    mpz_clear(bnum); mpz_clear(bden); 
    mpz_clear(cnum); mpz_clear(cden); 
    mpz_clear(dnum); mpz_clear(dden);

    return max_size;
}

double
qisq2_map_get_tolerance()
{
    return 0.0;
}
