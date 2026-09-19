#include <stdio.h>
#include <stdlib.h>
#include <sylvan_edge_weights_qisq2.h>
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
    qisq2_t a;
    a = qisq2_one();     EVBDD_ONE     = _weight_qisq2_lookup_ptr(&a, wgt_store);
    a = qisq2_zero();    EVBDD_ZERO    = _weight_qisq2_lookup_ptr(&a, wgt_store);
    a = qisq2_mone();    EVBDD_MIN_ONE = _weight_qisq2_lookup_ptr(&a, wgt_store);

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

    //qisq2_clear(x); // deallocate previous memory of x

    // copy result of multiplication to x
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
