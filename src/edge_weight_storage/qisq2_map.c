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
static const uint64_t EMPTY = 14738995463583502973ull;
static const uint64_t LOCK  = 14738995463583502974ull;
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
*/
uint32_t 
qisq2_hash(qisq2_t *num) {
    mpz_t anum, aden, bnum, bden, cnum, cden, dnum, dden;
    mpz_init(anum); mpz_init(aden); mpz_init(bnum); mpz_init(bden);
    mpz_init(cnum); mpz_init(cden); mpz_init(dnum); mpz_init(dden);
    
    mpz_t temp;
    mpz_init(temp);
    uint32_t result;

    mpq_get_num(anum, num->a);
    mpq_get_den(aden, num->a);
    mpq_get_num(bnum, num->b);
    mpq_get_den(bden, num->b);
    mpq_get_num(cnum, num->c);
    mpq_get_den(cden, num->c);
    mpq_get_num(dnum, num->d);
    mpq_get_den(dden, num->d);

    // hash 0 and 1 to 0 resp. 1
    if (!(mpz_cmp_si(anum,0) || mpz_cmp_si(bnum,0) || mpz_cmp_si(cnum,0) || mpz_cmp_si(dnum,0))){
        result = 0;
    }
    else if (!(mpz_cmp_si(anum,1) || mpz_cmp_si(aden,1) || mpz_cmp_si(bnum,0) || mpz_cmp_si(cnum,0) || mpz_cmp_si(dnum,0))){
        result = 1;
    }
    else {
        mpz_addmul_ui(temp, anum, 3443);
        mpz_addmul_ui(temp, aden, 17);
        mpz_addmul_ui(temp, bnum, 8093);
        mpz_addmul_ui(temp, bden, 9);
        mpz_addmul_ui(temp, cnum, 7314);
        mpz_addmul_ui(temp, cden, 71);
        mpz_addmul_ui(temp, dnum, 2817);
        mpz_addmul_ui(temp, dden, 7);
        unsigned long int res = mpz_get_ui(temp);
        result = (uint32_t) res;
    }
    
    mpz_clear(anum); mpz_clear(aden); 
    mpz_clear(bnum); mpz_clear(bden); 
    mpz_clear(cnum); mpz_clear(cden); 
    mpz_clear(dnum); mpz_clear(dden);

    mpz_clear(temp);

    return result;
}

void qisq2_reduce(qisq2_t *num) {
    mpq_canonicalize(num->a);
    mpq_canonicalize(num->b);
    mpq_canonicalize(num->c);
    mpq_canonicalize(num->d);
}

static bool
qisq2_equal(qisq2_t *in_table, const qisq2_t* to_insert)
{
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


int
qisq2_map_find_or_put(const void *dbs, const void *_v, uint64_t *ret)
{
    qisq2_t *v = (qisq2_t *)_v;
    qisq2_map_t *qisq2_map = (qisq2_map_t *) dbs;
    bucket_t *val  = (bucket_t *) v;

    qisq2_reduce(v);
    uint32_t hash = qisq2_hash(v);
    uint32_t prime = odd_primes[hash & PRIME_MASK];

    assert (val->d[0] != LOCK);
    assert (val->d[0] != EMPTY);

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
                    // write backwards (overwrite bucket->d[0] last)
                    for (int k = entry_size-1; k >= 0; k--) {
                        atomic_write (&bucket->d[k], val->d[k]);
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
                qisq2_clear(v);
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
    for (unsigned int c = 0; c < qisq2_map->size; c++) {
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
    qisq2_map->table = calloc (qisq2_map->size, sizeof(bucket_t));
    for (unsigned int c = 0; c < qisq2_map->size; c++) {
        qisq2_map->table[c].d[0] = EMPTY;
    }
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
    for (unsigned int c = 0; c < qisq2_map->size; c++) {
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
