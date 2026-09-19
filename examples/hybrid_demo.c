/* Demonstration for the hybrid (relative + zero-collapse) edge-weight merging
 * rule.  Applies H to every qubit of |0..0>, which makes the root weight of the
 * state 2^{-n/2} -- an exponentially small *scale factor*.  An absolute merging
 * tolerance cannot distinguish that legitimate value from noise once
 * 2^{-n/2} < delta, and the state collapses; a relative tolerance can.
 * This is the regime of mqt-core issue #575.
 */
#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sylvan.h>
#include <qsylvan.h>

typedef struct {
    int n;
    int layers;
    int strat;
    int hybrid;
    double tol;
    double rel;
    double zero;
} demo_opts_t;

// Sylvan operations must run inside a Lace worker (Lace >= 1.6); the thread
// that calls lace_start() is not one.
TASK_1(int, run_demo, demo_opts_t*, opts)
{
    int n      = opts->n;
    int layers = opts->layers;
    int strat  = opts->strat;
    int hybrid = opts->hybrid;
    double tol = opts->tol, rel = opts->rel, zero = opts->zero;

    QMDD s = qmdd_create_all_zero_state(n);
    for (int l = 0; l < layers; l++)
        for (int k = 0; k < n; k++)
            s = qmdd_gate(s, GATEID_H, k);

    double norm = qmdd_get_magnitude(s, n);
    bool *zerostate = calloc(n, sizeof(bool));
    complex_t amp = qmdd_get_amplitude(s, zerostate, n);
    /* after an even number of H layers the state is exactly |0..0> */
    double want = (layers % 2 == 0) ? 1.0 : pow(2.0, -0.5*n);
    printf("%d,%d,%d,%s,%.3e,%.3e,%.17g,%.17g,%.3e,%" PRIu64 "\n",
           n, layers, strat, hybrid ? "hybrid" : "absolute",
           hybrid ? rel : tol, hybrid ? zero : 0.0,
           norm, (double) amp.r, fabs((double) amp.r - want),
           (uint64_t) evbdd_countnodes(s));
    free(zerostate);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 5) {
        fprintf(stderr, "usage: %s <nqubits> <layers> <norm-strat 0|1|2|3> "
                        "<abs TOL | rel REL ZERO>\n", argv[0]);
        return 1;
    }
    demo_opts_t opts;
    opts.n       = atoi(argv[1]);
    opts.layers  = atoi(argv[2]);
    opts.strat   = atoi(argv[3]);
    opts.hybrid  = (strcmp(argv[4], "rel") == 0);
    opts.tol = 1e-14; opts.rel = 1e-14; opts.zero = 0.0;
    if (opts.hybrid) { opts.rel = atof(argv[5]); opts.zero = (argc > 6) ? atof(argv[6]) : 0.0; }
    else             { opts.tol = atof(argv[5]); }

    lace_start(1, 0);
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<20, 1LL<<20);
    sylvan_init_package();
    if (opts.hybrid) sylvan_edge_weights_set_hybrid_tolerance(opts.rel, opts.zero);
    qsylvan_init_simulator(1LL<<23, 1LL<<23, opts.tol, COMP_HASHMAP, opts.strat);

    int res = RUN(run_demo, &opts);

    sylvan_quit();
    lace_stop();
    return res;
}
