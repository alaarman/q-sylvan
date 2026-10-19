#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sylvan_edge_weights_qisq2.h>
#include <sylvan_platform.h>
#include <assert.h>


/**********************<Some static utility functions>*************************/


/**********************</Some static utility functions>************************/




/*****************<Implementation of edge_weights interface>*******************/

void qisq2_init(qisq2_t *num) {
    mpq_init(num->a);
    mpq_init(num->b);
    mpq_init(num->c);
    mpq_init(num->d);
}

/* qisq2_clear and qisq2_reduce live in edge_weight_storage/qisq2_map.c, which
 * needs them itself; defining them here as well is a duplicate symbol. */

void weight_qisq2_copy(qisq2_t *a, qisq2_t *a_copy){
    qisq2_init(a_copy);
    mpq_set(a_copy->a,a->a);
    mpq_set(a_copy->b,a->b);
    mpq_set(a_copy->c,a->c);
    mpq_set(a_copy->d,a->d);
}

qisq2_t *
/*
 * Raw memory: the four mpq_t are NOT initialised. Either qisq2_init it, or
 * fill it with weight_value, which initialises as it copies. Release it with
 * weight_qisq2_free once it holds GMP storage.
 */
weight_qisq2_malloc()
{
    qisq2_t *res = malloc(sizeof(qisq2_t));
    return res;
}

void weight_qisq2_set(qisq2_t *x, qisq2_t *y);

void
_weight_qisq2_value(void *wgt_store, EVBDD_WGT a, qisq2_t *res)
{
    /*
     * A deep copy, not a struct assignment. Each mpq_t holds a pointer to its
     * limbs, so copying the struct would leave `res` sharing storage with the
     * table entry -- and the very next qisq2_reduce would canonicalise in
     * place and rewrite the stored weight. That is what made wgt_neg change
     * the value it was reading.
     *
     * `res` is initialised here rather than by the caller, because callers
     * pass both weight_malloc'd buffers and plain stack variables. It owns its
     * limbs afterwards and must be released with weight_qisq2_free.
     */
    qisq2_init(res);
    weight_qisq2_set(res, (qisq2_t *)(wgt_store_get(wgt_store, a)));
}

void
weight_qisq2_free(qisq2_t *a)
{
    qisq2_clear(a);
    free(a);
}

EVBDD_WGT
_weight_qisq2_lookup_ptr(qisq2_t *a, void *wgt_store)
{
    // TODO: catch czero() / cone() here?
    uint64_t res;
    bool success;

    int present = wgt_store_find_or_put(wgt_store, a, &res);
    if (present == -1) {
        success = false;
    } else if (present == 0) { 
        success = true;
        wgt_table_gc_inc_entries_estimate();
    } else {
        success = true;
    }

    if (!success) {
        fprintf(stderr, "Amplitude table full!\n");
        exit(1);
    }
    return (EVBDD_WGT) res; 
}

EVBDD_WGT
weight_qisq2_lookup(qisq2_t *a)
{
    return _weight_qisq2_lookup_ptr(a, wgt_storage);
}

void
init_qisq2_one_zero(void *wgt_store)
{
    /* each qisq2_make initialises limbs of its own, and the table keeps a
     * copy, so each is released after its lookup */
    qisq2_t a;
    a = qisq2_one();     EVBDD_ONE     = _weight_qisq2_lookup_ptr(&a, wgt_store); qisq2_clear(&a);
    a = qisq2_zero();    EVBDD_ZERO    = _weight_qisq2_lookup_ptr(&a, wgt_store); qisq2_clear(&a);
    a = qisq2_mone();    EVBDD_MIN_ONE = _weight_qisq2_lookup_ptr(&a, wgt_store); qisq2_clear(&a);

}

void
weight_qisq2_abs_sqr(qisq2_t *a)
{
    qisq2_t temp;
    qisq2_init(&temp);
    weight_qisq2_complexConjugate(&temp,a);
    weight_qisq2_mul(a,&temp);
    qisq2_clear(&temp);
}

void
weight_qisq2_abs(qisq2_t *a)
{
    qisq2_t temp;
    qisq2_init(&temp);
    weight_qisq2_complexConjugate(&temp,a);
    weight_qisq2_mul(a,&temp);
    qisq2_clear(&temp);
    if (mpq_get_d(a->b)==0.0){
        mpq_set_d(a->a,sqrt(mpq_get_d(a->a)));
    } 
    else if (mpq_get_d(a->a)==0.0){
        mpq_set_d(a->b,sqrt(mpq_get_d(a->b)));
    }
    else{
        //TODO: implement abs for general case
        printf("NotImplementedError \n");
        printf("absolute value of \n");
        weight_qisq2_print(a);
        printf("cannot be calculated in qisq2\n");
        exit(0);
    }
}

void
weight_qisq2_neg(qisq2_t *num)
{
    mpq_neg(num->a, num->a);
    mpq_neg(num->b, num->b);
    mpq_neg(num->c, num->c);
    mpq_neg(num->d, num->d);
    qisq2_reduce(num);
}

void
weight_qisq2_conj(qisq2_t *x)
{
    mpq_neg(x->c, x->c);
    mpq_neg(x->d, x->d);
    qisq2_reduce(x);
}

void
weight_qisq2_sqr(qisq2_t *a)
{
    weight_qisq2_mul(a, a);
}

void
weight_qisq2_add(qisq2_t *x, qisq2_t *y)
{
    qisq2_t result;
    qisq2_init(&result);
    mpq_add(result.a, x->a, y->a);
    mpq_add(result.b, x->b, y->b);
    mpq_add(result.c, x->c, y->c);
    mpq_add(result.d, x->d, y->d);
    qisq2_reduce(&result);

    qisq2_clear(x);             /* see weight_qisq2_mul */
    *x = result;
}

void
weight_qisq2_sub(qisq2_t *x, qisq2_t *y)
{
    qisq2_t result;
    qisq2_init(&result);
    mpq_sub(result.a, x->a, y->a);
    mpq_sub(result.b, x->b, y->b);
    mpq_sub(result.c, x->c, y->c);
    mpq_sub(result.d, x->d, y->d);
    qisq2_reduce(&result);

    qisq2_clear(x);             /* see weight_qisq2_mul */
    *x = result;
}

void
weight_qisq2_mul(qisq2_t *x, qisq2_t *y)
{
    qisq2_t result;
    qisq2_init(&result);
    mpq_t temp;
    mpq_init(temp);
    mpq_t two;
    mpq_init(two);
    mpq_set_ui(two, 2, 1);

    mpq_mul(result.a, x->a, y->a);
    mpq_mul(temp, x->b, y->b);
    mpq_mul(temp, temp, two);
    mpq_add(result.a, result.a, temp);

    mpq_mul(temp, x->c, y->c);
    mpq_sub(result.a, result.a, temp);

    mpq_mul(temp, x->d, y->d);
    mpq_mul(temp, temp, two);
    mpq_sub(result.a, result.a, temp);

    mpq_mul(result.b, x->a, y->b);

    mpq_mul(temp, x->c, y->d);
    mpq_sub(result.b, result.b, temp);

    mpq_mul(temp, x->b, y->a);
    mpq_add(result.b, result.b, temp);

    mpq_mul(temp, x->d, y->c);
    mpq_sub(result.b, result.b, temp);

    mpq_mul(result.c, x->a, y->c);

    mpq_mul(temp, x->b, y->d);
    mpq_mul(temp, temp, two);
    mpq_add(result.c, result.c, temp);

    mpq_mul(temp, x->c, y->a);
    mpq_add(result.c, result.c, temp);

    mpq_mul(temp, x->d, y->b);
    mpq_mul(temp, temp, two);
    mpq_add(result.c, result.c, temp);

    mpq_mul(result.d, x->a, y->d);

    mpq_mul(temp, x->d, y->a);
    mpq_add(result.d, result.d, temp);

    mpq_mul(temp, x->b, y->c);
    mpq_add(result.d, result.d, temp);

    mpq_mul(temp, x->c, y->b);
    mpq_add(result.d, result.d, temp);
    
    mpq_clear(temp);
    mpq_clear(two);

    /*
     * x owns its limbs, and the struct assignment below replaces the
     * pointers to them, so they are released first or never. They were
     * never: every wgt_mul, wgt_div, wgt_add and wgt_sub that missed the
     * cache leaked the limbs of its first operand, and a divide leaked those
     * of its intermediate products as well -- 2.9 million blocks, 70 MB, over
     * one run of test_bqd. The clear was commented out
     * when a "copy" here could still alias a table entry. Nothing can now:
     * weight_value deep-copies (6ff7afe), and every caller passes a value it
     * made that way or initialised itself. y may be x (weight_qisq2_sqr),
     * which is safe, since result is complete before x is cleared.
     */
    qisq2_clear(x);
    *x = result;
    qisq2_reduce(x);
}

void
weight_qisq2_complexConjugate(qisq2_t *result, qisq2_t *x) {    
    mpq_set(result->a, x->a);
    mpq_set(result->b, x->b);
    mpq_neg(result->c, x->c);
    mpq_neg(result->d, x->d);
    
    qisq2_reduce(result);
}

void
weight_qisq2_sqrttwoConjugate(qisq2_t *result, qisq2_t *x) {
    mpq_set(result->a, x->a);
    mpq_neg(result->b, x->b);
    mpq_set(result->c, x->c);
    mpq_neg(result->d, x->d);

    qisq2_reduce(result);
}

void
qisq2_mul_clear_input(qisq2_t *x, qisq2_t *y) {
    /* A struct assignment here would alias x's limbs, and the qisq2_clear
     * below would then free them out from under x. */
    qisq2_t *tmp;
    tmp = weight_qisq2_malloc();
    qisq2_init(tmp);
    weight_qisq2_set(tmp, x);
    weight_qisq2_mul(x,y);
    qisq2_clear(tmp);
    free(tmp);
}

void
weight_qisq2_set(qisq2_t *x, qisq2_t *y){
    mpq_set(x->a,y->a);
    mpq_set(x->b,y->b);
    mpq_set(x->c,y->c);
    mpq_set(x->d,y->d);
}

void
weight_qisq2_div(qisq2_t *x, qisq2_t *y)
{
    //check if y is non-zero
    //assert(!(mpq_sgn(y->a)==0 || mpq_sgn(y->b)==0 || mpq_sgn(y->c)==0 || mpq_sgn(y->d)==0));

    qisq2_t *temp;
    qisq2_t *cc;
    qisq2_t *sc;
    temp = weight_qisq2_malloc();
    qisq2_init(temp); 
    cc = weight_qisq2_malloc();
    sc = weight_qisq2_malloc();
    qisq2_init(cc); 
    qisq2_init(sc);

    weight_qisq2_set(temp,y);//*temp = *y; // ---> gaat dit goed??? Er moet een mpq_copy functie aangeroepen worden.

    weight_qisq2_complexConjugate(cc, temp);
    qisq2_mul_clear_input(temp, cc); 

    weight_qisq2_sqrttwoConjugate(sc, temp);
    qisq2_mul_clear_input(temp, sc); 

    // y is real fraction now
    // this can be checked with the following assertion (for debugging)
    //assert(mpq_cmp_ui(temp.b, 0, 1) == 0 && mpq_cmp_ui(temp.c, 0, 1) == 0 && mpq_cmp_ui(temp.d, 0, 1) == 0);

    // multiply numerator x with conjugates

    weight_qisq2_mul(x, cc);
    qisq2_mul_clear_input(x, sc);

    qisq2_clear(cc);
    qisq2_clear(sc);
    free(cc);
    free(sc);

    // calculate 1/y
    mpq_inv(temp->a, temp->a);

    qisq2_mul_clear_input(x, temp);
    qisq2_clear(temp);
    free(temp);
}

bool
weight_qisq2_eq(qisq2_t *x, qisq2_t *y)
{
    if (mpq_equal(x->a,y->a) == 0){
        return false;
    }
    if (mpq_equal(x->b,y->b) == 0){
        return false;
    }
    if (mpq_equal(x->c,y->c) == 0){
        return false;
    }
    if (mpq_equal(x->d,y->d) == 0){
        return false;
    }
    return true;
}

bool
weight_qisq2_eps_close(qisq2_t *x, qisq2_t *y, double eps)
{
    if (eps == 0.0){
        return weight_qisq2_eq(x, y);
    }
    else {
        // TODO
        printf("weight_qisq2_eps_close() is not implemented yet");
        return weight_qisq2_eq(x, y);
    }    
}

bool
weight_qisq2_greater(qisq2_t *x, qisq2_t *y)
{
    if (mpq_cmp(x->a,y->a) > 0){
        return false;
    }
    if (mpq_cmp(x->b,y->b) > 0){
        return false;
    }
    if (mpq_cmp(x->c,y->c) > 0){
        return false;
    }
    if (mpq_cmp(x->c,y->c) > 0){
        return false;
    }
    return true;
}

void
weight_qisq2_fprint(FILE *stream, qisq2_t *x)
{
    fprintf(stream, "\t%f + %f sqrt(2) + %f i + %f i sqrt(2)\n",
               mpq_get_d(x->a), mpq_get_d(x->b), mpq_get_d(x->c), mpq_get_d(x->d));
}

void
weight_qisq2_print(qisq2_t *x)
{
    printf("\t%f + %f sqrt(2) + %f i + %f i sqrt(2)\n",
               mpq_get_d(x->a), mpq_get_d(x->b), mpq_get_d(x->c), mpq_get_d(x->d));
    gmp_printf("\t%#4Qx + \n\t%#4Qx sqrt(2) +\n\t%#4Qx i +\n\t%#4Qx i sqrt(2)\n\n", x->a, x->b, x->c, x->d);
}

/*****************</Implementation of edge_weights interface>******************/


/*******************<Arithmetic on interned weights, by index>*****************/

/*
 * The generic wgt_* path copies both operands out of the table (eight mpq per
 * operand, each an allocation), works in fresh temporaries and frees all of it
 * again: on exact runs malloc and free took a fifth of the time. These read the
 * stored operands in place -- table entries are never written after they are
 * published -- and compute into temporaries that belong to the calling thread
 * and keep their limbs from one call to the next. Every result of GMP's mpq
 * arithmetic is reduced, so it goes into the table without reducing it again.
 * The values are exactly those of weight_qisq2_mul and friends.
 */
typedef struct {
    bool    ready;
    qisq2_t r;          /* the result */
    qisq2_t inv;        /* 1/y for a division */
    mpq_t   t, u;       /* products and partial sums */
    mpq_t   p, q, n;    /* the inverse's p, q and norm */
} qisq2_scratch_t;

static SYLVAN_TLS qisq2_scratch_t qisq2_scratch;

static inline qisq2_scratch_t *
scratch(void)
{
    qisq2_scratch_t *s = &qisq2_scratch;
    if (!s->ready) {
        qisq2_init(&s->r);
        qisq2_init(&s->inv);
        mpq_init(s->t); mpq_init(s->u);
        mpq_init(s->p); mpq_init(s->q); mpq_init(s->n);
        s->ready = true;
    }
    return s;
}

static inline const qisq2_t *
stored(EVBDD_WGT a)
{
    return (const qisq2_t *) wgt_store_get(wgt_storage, a);
}

static EVBDD_WGT
intern_reduced(const qisq2_t *v)
{
    uint64_t res;
    const int present = qisq2_map_find_or_put_reduced(wgt_storage, v, &res);
    if (present == -1) {
        fprintf(stderr, "Amplitude table full!\n");
        exit(1);
    }
    if (present == 0) wgt_table_gc_inc_entries_estimate();
    return (EVBDD_WGT) res;
}

/* +i or -i, read off the stored value: 1 for i, -1 for -i, 0 otherwise. */
static inline int
unit_i(const qisq2_t *x)
{
    if (mpq_sgn(x->a) != 0 || mpq_sgn(x->b) != 0 || mpq_sgn(x->d) != 0) return 0;
    if (mpq_cmp_si(x->c, 1, 1) == 0) return 1;
    if (mpq_cmp_si(x->c, -1, 1) == 0) return -1;
    return 0;
}

/* r = x * y; r is neither x nor y. */
static void
mul_into(qisq2_t *r, const qisq2_t *x, const qisq2_t *y, mpq_ptr t)
{
    /* a: xa ya + 2 xb yb - xc yc - 2 xd yd */
    mpq_mul(r->a, x->a, y->a);
    mpq_mul(t, x->b, y->b); mpq_mul_2exp(t, t, 1); mpq_add(r->a, r->a, t);
    mpq_mul(t, x->c, y->c); mpq_sub(r->a, r->a, t);
    mpq_mul(t, x->d, y->d); mpq_mul_2exp(t, t, 1); mpq_sub(r->a, r->a, t);
    /* b: xa yb + xb ya - xc yd - xd yc */
    mpq_mul(r->b, x->a, y->b);
    mpq_mul(t, x->b, y->a); mpq_add(r->b, r->b, t);
    mpq_mul(t, x->c, y->d); mpq_sub(r->b, r->b, t);
    mpq_mul(t, x->d, y->c); mpq_sub(r->b, r->b, t);
    /* c: xa yc + xc ya + 2 xb yd + 2 xd yb */
    mpq_mul(r->c, x->a, y->c);
    mpq_mul(t, x->c, y->a); mpq_add(r->c, r->c, t);
    mpq_mul(t, x->b, y->d); mpq_mul_2exp(t, t, 1); mpq_add(r->c, r->c, t);
    mpq_mul(t, x->d, y->b); mpq_mul_2exp(t, t, 1); mpq_add(r->c, r->c, t);
    /* d: xa yd + xd ya + xb yc + xc yb */
    mpq_mul(r->d, x->a, y->d);
    mpq_mul(t, x->d, y->a); mpq_add(r->d, r->d, t);
    mpq_mul(t, x->b, y->c); mpq_add(r->d, r->d, t);
    mpq_mul(t, x->c, y->b); mpq_add(r->d, r->d, t);
}

/* r = i^k x for k = 1 (i), 2 (-1) or 3 (-i); i (a + b sqrt2 + i(c + d sqrt2))
 * = -(c + d sqrt2) + i (a + b sqrt2). */
static void
unit_into(qisq2_t *r, const qisq2_t *x, int k)
{
    switch (k) {
    case 1:  mpq_neg(r->a, x->c); mpq_neg(r->b, x->d); mpq_set(r->c, x->a); mpq_set(r->d, x->b); break;
    case 2:  mpq_neg(r->a, x->a); mpq_neg(r->b, x->b); mpq_neg(r->c, x->c); mpq_neg(r->d, x->d); break;
    default: mpq_set(r->a, x->c); mpq_set(r->b, x->d); mpq_neg(r->c, x->a); mpq_neg(r->d, x->b); break;
    }
}

/*
 * r = 1/y, y != 0. With y = e + f sqrt2 + i(g + h sqrt2), |y|^2 = p + q sqrt2
 * for p = e^2 + 2f^2 + g^2 + 2h^2 and q = 2(ef + gh), and multiplying by the
 * conjugate p - q sqrt2 leaves the rational n = p^2 - 2q^2, which is |y|^2
 * times |y'|^2 for y' the image of y under sqrt2 -> -sqrt2, so it is not 0:
 *     1/y = [(e p - 2 f q) + (f p - e q) sqrt2 + i((2 h q - g p) + (g q - h p) sqrt2)] / n.
 * About twenty rational operations, where dividing through two conjugates
 * took some two hundred.
 */
static void
inv_into(qisq2_t *r, const qisq2_t *y, qisq2_scratch_t *s)
{
    mpq_ptr t = s->t, u = s->u, p = s->p, q = s->q, n = s->n;
    mpq_mul(p, y->a, y->a);
    mpq_mul(t, y->b, y->b); mpq_mul_2exp(t, t, 1); mpq_add(p, p, t);
    mpq_mul(t, y->c, y->c); mpq_add(p, p, t);
    mpq_mul(t, y->d, y->d); mpq_mul_2exp(t, t, 1); mpq_add(p, p, t);
    mpq_mul(q, y->a, y->b);
    mpq_mul(t, y->c, y->d); mpq_add(q, q, t); mpq_mul_2exp(q, q, 1);
    mpq_mul(n, p, p);
    mpq_mul(t, q, q); mpq_mul_2exp(t, t, 1); mpq_sub(n, n, t);
    assert(mpq_sgn(n) > 0);
    /* a = (e p - 2 f q) / n */
    mpq_mul(r->a, y->a, p); mpq_mul(t, y->b, q); mpq_mul_2exp(t, t, 1); mpq_sub(r->a, r->a, t);
    mpq_div(r->a, r->a, n);
    /* b = (f p - e q) / n */
    mpq_mul(r->b, y->b, p); mpq_mul(t, y->a, q); mpq_sub(r->b, r->b, t);
    mpq_div(r->b, r->b, n);
    /* c = (2 h q - g p) / n */
    mpq_mul(u, y->d, q); mpq_mul_2exp(u, u, 1); mpq_mul(t, y->c, p); mpq_sub(r->c, u, t);
    mpq_div(r->c, r->c, n);
    /* d = (g q - h p) / n */
    mpq_mul(r->d, y->c, q); mpq_mul(t, y->d, p); mpq_sub(r->d, r->d, t);
    mpq_div(r->d, r->d, n);
}

/*
 * Products and quotients over one common denominator.
 *
 * mul_into and inv_into work on the four components as separate rationals:
 * sixteen products and twelve sums, each of which GMP reduces with a gcd of
 * its own. Once weights grow to thousands of bits that is the run: GMP's gcd
 * was 62% of an exact BQD run on clifford_T_circuit_20_700 and 90% of an EVDD
 * run on synth_qft_n6_e3. Here the operands are lifted to four integers over
 * one denominator, multiplied (or inverted) in integers, and only the result
 * is reduced, four gcds in all. The value is the same, and so, since the
 * reduced form is unique, are the stored rationals.
 *
 * For operands whose eight integers all fit a limb the separate rationals are
 * cheaper (no lifting, gcds of single limbs), so those keep the old path.
 */
typedef struct { mpz_t a, b, c, d, den; } qisq2_cd_t;

typedef struct {
    bool        ready;
    qisq2_cd_t  x, y, inv, r;
    mpz_t       t, p, q, n;
} qisq2_cd_scratch_t;

static SYLVAN_TLS qisq2_cd_scratch_t qisq2_cd_scratch;

static inline void
cd_init(qisq2_cd_t *x)
{
    mpz_init(x->a); mpz_init(x->b); mpz_init(x->c); mpz_init(x->d); mpz_init(x->den);
}

static inline qisq2_cd_scratch_t *
cd_scratch(void)
{
    qisq2_cd_scratch_t *s = &qisq2_cd_scratch;
    if (!s->ready) {
        cd_init(&s->x); cd_init(&s->y); cd_init(&s->inv); cd_init(&s->r);
        mpz_init(s->t); mpz_init(s->p); mpz_init(s->q); mpz_init(s->n);
        s->ready = true;
    }
    return s;
}

#ifndef QISQ2_CD_MIN_LIMBS
#define QISQ2_CD_MIN_LIMBS 2
#endif
/* True iff every numerator and denominator of x has fewer than
 * QISQ2_CD_MIN_LIMBS limbs. */
static inline bool
qisq2_small(const qisq2_t *x)
{
    return mpz_size(mpq_numref(x->a)) < QISQ2_CD_MIN_LIMBS && mpz_size(mpq_denref(x->a)) < QISQ2_CD_MIN_LIMBS &&
           mpz_size(mpq_numref(x->b)) < QISQ2_CD_MIN_LIMBS && mpz_size(mpq_denref(x->b)) < QISQ2_CD_MIN_LIMBS &&
           mpz_size(mpq_numref(x->c)) < QISQ2_CD_MIN_LIMBS && mpz_size(mpq_denref(x->c)) < QISQ2_CD_MIN_LIMBS &&
           mpz_size(mpq_numref(x->d)) < QISQ2_CD_MIN_LIMBS && mpz_size(mpq_denref(x->d)) < QISQ2_CD_MIN_LIMBS;
}

/* r = x over the least common denominator of its components. */
static void
to_cd(qisq2_cd_t *r, const qisq2_t *x)
{
    mpz_lcm(r->den, mpq_denref(x->a), mpq_denref(x->b));
    mpz_lcm(r->den, r->den, mpq_denref(x->c));
    mpz_lcm(r->den, r->den, mpq_denref(x->d));
    mpq_srcptr q[4] = { x->a, x->b, x->c, x->d };
    mpz_ptr z[4] = { r->a, r->b, r->c, r->d };
    for (int i = 0; i < 4; i++) {
        if (mpz_cmp(mpq_denref(q[i]), r->den) == 0) {
            mpz_set(z[i], mpq_numref(q[i]));
        } else {
            mpz_divexact(z[i], r->den, mpq_denref(q[i]));
            mpz_mul(z[i], z[i], mpq_numref(q[i]));
        }
    }
}

/* x = r, each component reduced on its own: the stored form. */
static void
from_cd(qisq2_t *x, const qisq2_cd_t *r)
{
    mpq_ptr q[4] = { x->a, x->b, x->c, x->d };
    mpz_srcptr z[4] = { r->a, r->b, r->c, r->d };
    for (int i = 0; i < 4; i++) {
        mpq_set_num(q[i], z[i]);
        mpq_set_den(q[i], r->den);
        mpq_canonicalize(q[i]);
    }
}

/* r = x y in integers; r is neither x nor y. Both denominators are positive. */
static void
cd_mul(qisq2_cd_t *r, const qisq2_cd_t *x, const qisq2_cd_t *y, mpz_ptr t)
{
    /* a: xa ya + 2 xb yb - xc yc - 2 xd yd */
    mpz_mul(r->a, x->a, y->a);
    mpz_mul(t, x->b, y->b); mpz_addmul_ui(r->a, t, 2);
    mpz_submul(r->a, x->c, y->c);
    mpz_mul(t, x->d, y->d); mpz_submul_ui(r->a, t, 2);
    /* b: xa yb + xb ya - xc yd - xd yc */
    mpz_mul(r->b, x->a, y->b); mpz_addmul(r->b, x->b, y->a);
    mpz_submul(r->b, x->c, y->d); mpz_submul(r->b, x->d, y->c);
    /* c: xa yc + xc ya + 2 xb yd + 2 xd yb */
    mpz_mul(r->c, x->a, y->c); mpz_addmul(r->c, x->c, y->a);
    mpz_mul(t, x->b, y->d); mpz_addmul_ui(r->c, t, 2);
    mpz_mul(t, x->d, y->b); mpz_addmul_ui(r->c, t, 2);
    /* d: xa yd + xd ya + xb yc + xc yb */
    mpz_mul(r->d, x->a, y->d); mpz_addmul(r->d, x->d, y->a);
    mpz_addmul(r->d, x->b, y->c); mpz_addmul(r->d, x->c, y->b);
    mpz_mul(r->den, x->den, y->den);
}

/*
 * r = 1/y in integers, as inv_into: with y = (E + F sqrt2 + i(G + H sqrt2))/D,
 * P = E^2 + 2F^2 + G^2 + 2H^2, Q = 2(EF + GH) and N = P^2 - 2Q^2 > 0,
 *     1/y = D [(EP - 2FQ) + (FP - EQ) sqrt2 + i((2HQ - GP) + (GQ - HP) sqrt2)] / N.
 */
static void
cd_inv(qisq2_cd_t *r, const qisq2_cd_t *y, qisq2_cd_scratch_t *s)
{
    mpz_ptr t = s->t, P = s->p, Q = s->q, N = s->n;
    mpz_mul(P, y->a, y->a);
    mpz_mul(t, y->b, y->b); mpz_addmul_ui(P, t, 2);
    mpz_addmul(P, y->c, y->c);
    mpz_mul(t, y->d, y->d); mpz_addmul_ui(P, t, 2);
    mpz_mul(Q, y->a, y->b); mpz_addmul(Q, y->c, y->d); mpz_mul_2exp(Q, Q, 1);
    mpz_mul(N, P, P);
    mpz_mul(t, Q, Q); mpz_submul_ui(N, t, 2);
    assert(mpz_sgn(N) > 0);
    mpz_mul(r->a, y->a, P); mpz_mul(t, y->b, Q); mpz_submul_ui(r->a, t, 2); mpz_mul(r->a, r->a, y->den);
    mpz_mul(r->b, y->b, P); mpz_submul(r->b, y->a, Q); mpz_mul(r->b, r->b, y->den);
    mpz_mul(r->c, y->d, Q); mpz_mul_2exp(r->c, r->c, 1); mpz_submul(r->c, y->c, P); mpz_mul(r->c, r->c, y->den);
    mpz_mul(r->d, y->c, Q); mpz_submul(r->d, y->d, P); mpz_mul(r->d, r->d, y->den);
    mpz_set(r->den, N);
}

EVBDD_WGT
qisq2_wgt_mul(EVBDD_WGT a, EVBDD_WGT b)
{
    qisq2_scratch_t *s = scratch();
    const qisq2_t *x = stored(a), *y = stored(b);
    if (b == EVBDD_MIN_ONE) { unit_into(&s->r, x, 2); return intern_reduced(&s->r); }
    if (a == EVBDD_MIN_ONE) { unit_into(&s->r, y, 2); return intern_reduced(&s->r); }
    int k;
    if ((k = unit_i(y)) != 0) { unit_into(&s->r, x, k == 1 ? 1 : 3); return intern_reduced(&s->r); }
    if ((k = unit_i(x)) != 0) { unit_into(&s->r, y, k == 1 ? 1 : 3); return intern_reduced(&s->r); }
    if (qisq2_small(x) && qisq2_small(y)) {
        mul_into(&s->r, x, y, s->t);
    } else {
        qisq2_cd_scratch_t *c = cd_scratch();
        to_cd(&c->x, x); to_cd(&c->y, y);
        cd_mul(&c->r, &c->x, &c->y, c->t);
        from_cd(&s->r, &c->r);
    }
    return intern_reduced(&s->r);
}

EVBDD_WGT
qisq2_wgt_inv(EVBDD_WGT b)
{
    qisq2_scratch_t *s = scratch();
    const qisq2_t *y = stored(b);
    if (b == EVBDD_MIN_ONE) return b;
    const int k = unit_i(y);
    if (k != 0) { unit_into(&s->r, stored(EVBDD_ONE), k == 1 ? 3 : 1); return intern_reduced(&s->r); }
    if (qisq2_small(y)) {
        inv_into(&s->r, y, s);
    } else {
        qisq2_cd_scratch_t *c = cd_scratch();
        to_cd(&c->y, y);
        cd_inv(&c->inv, &c->y, c);
        from_cd(&s->r, &c->inv);
    }
    return intern_reduced(&s->r);
}

EVBDD_WGT
qisq2_wgt_div(EVBDD_WGT a, EVBDD_WGT b)
{
    if (a == EVBDD_ONE) return qisq2_wgt_inv(b);
    qisq2_scratch_t *s = scratch();
    const qisq2_t *x = stored(a), *y = stored(b);
    if (b == EVBDD_MIN_ONE) { unit_into(&s->r, x, 2); return intern_reduced(&s->r); }
    const int k = unit_i(y);       /* x / (+-i) = x * (-+i) */
    if (k != 0) { unit_into(&s->r, x, k == 1 ? 3 : 1); return intern_reduced(&s->r); }
    if (qisq2_small(x) && qisq2_small(y)) {
        inv_into(&s->inv, y, s);
        mul_into(&s->r, x, &s->inv, s->t);
    } else {
        /* x * (1/y) without reducing 1/y in between */
        qisq2_cd_scratch_t *c = cd_scratch();
        to_cd(&c->x, x); to_cd(&c->y, y);
        cd_inv(&c->inv, &c->y, c);
        cd_mul(&c->r, &c->x, &c->inv, c->t);
        from_cd(&s->r, &c->r);
    }
    return intern_reduced(&s->r);
}

EVBDD_WGT
qisq2_wgt_add(EVBDD_WGT a, EVBDD_WGT b, bool sub)
{
    qisq2_scratch_t *s = scratch();
    const qisq2_t *x = stored(a), *y = stored(b);
    if (sub) {
        mpq_sub(s->r.a, x->a, y->a); mpq_sub(s->r.b, x->b, y->b);
        mpq_sub(s->r.c, x->c, y->c); mpq_sub(s->r.d, x->d, y->d);
    } else {
        mpq_add(s->r.a, x->a, y->a); mpq_add(s->r.b, x->b, y->b);
        mpq_add(s->r.c, x->c, y->c); mpq_add(s->r.d, x->d, y->d);
    }
    return intern_reduced(&s->r);
}

EVBDD_WGT
qisq2_wgt_neg(EVBDD_WGT a)
{
    qisq2_scratch_t *s = scratch();
    unit_into(&s->r, stored(a), 2);
    return intern_reduced(&s->r);
}

EVBDD_WGT
qisq2_wgt_conj(EVBDD_WGT a)
{
    qisq2_scratch_t *s = scratch();
    const qisq2_t *x = stored(a);
    mpq_set(s->r.a, x->a); mpq_set(s->r.b, x->b);
    mpq_neg(s->r.c, x->c); mpq_neg(s->r.d, x->d);
    return intern_reduced(&s->r);
}

void
qisq2_wgt_as_doubles(EVBDD_WGT a, double *re, double *im)
{
    const qisq2_t *q = stored(a);
    *re = mpq_get_d(q->a) + mpq_get_d(q->b) * sqrt(2.0);
    *im = mpq_get_d(q->c) + mpq_get_d(q->d) * sqrt(2.0);
}

/******************</Arithmetic on interned weights, by index>*****************/
