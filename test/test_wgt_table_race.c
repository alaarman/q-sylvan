/*
 * The edge weight tables under contention: every worker looks up the same
 * fresh values at the same time, and each value must come back with one
 * index, whichever worker inserted it.
 *
 * A table publishes an entry by writing its words and then its first word
 * last, over a LOCK, and a reader that finds the LOCK waits and compares.
 * With volatile accesses in place of acquire and release that is only right
 * on a CPU that keeps stores in order: on ARM a reader could compare against
 * an entry, or GMP limbs, not yet visible, find "not equal" and insert the
 * value again, and the exact LIMDD and BQD then came out with two edges for
 * one state now and then. This checks the tables directly, for the exact
 * backend and the float one.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sylvan.h>
#include <sylvan_int.h>

#include "qsylvan_simulator.h"
#include "sylvan_edge_weights_complex.h"
#include "sylvan_edge_weights_qisq2.h"

#define NWORKERS 8
#define NVALUES  4096
#define ROUNDS   40

static EVBDD_WGT got[NWORKERS][NVALUES];
static int backend;
static unsigned round_no;

/* value i of this round: fresh in every round, so the lookups insert */
VOID_TASK_0(lookup_all)
{
    const unsigned w = (unsigned)lace_get_worker()->worker;
    for (unsigned i = 0; i < NVALUES; i++) {
        const long v = (long)(round_no * NVALUES + i + 1);
        got[w][i] = (backend == QISQ2_MAP)
            ? qisq2_lookup(v, 3, -(v % 7), 5, (v * 13) % 11, 7, v % 3, 2)
            : complex_lookup(1.0 + (double)v * 1e-3, (double)(v % 97) * 0.5);
    }
}

TASK_0(int, run)
{
    unsigned dup = 0;
    for (round_no = 0; round_no < ROUNDS; round_no++) {
        TOGETHER(lookup_all);
        for (unsigned i = 0; i < NVALUES; i++)
            for (unsigned w = 1; w < NWORKERS; w++)
                if (got[w][i] != got[0][i]) { dup++; break; }
    }
    return (int)dup;
}

static int
session(int b, const char *name)
{
    backend = b;
    lace_start(NWORKERS, 0);
    sylvan_set_sizes(1LL << 20, 1LL << 20, 1LL << 20, 1LL << 20);
    sylvan_init_package();
    qsylvan_init_simulator(1LL << 22, 1LL << 22, -1, b, NORM_LOW);
    const int dup = RUN(run);
    sylvan_quit();
    lace_stop();
    printf("  %-6s %u values, %d workers: %d stored under two indices %s\n",
           name, NVALUES * ROUNDS, NWORKERS, dup, dup ? "FAILED" : "ok");
    return dup != 0;
}

int
main(void)
{
    int bad = 0;
    bad |= session(QISQ2_MAP, "qisq2");
    bad |= session(COMP_HASHMAP, "float");
    return bad;
}
