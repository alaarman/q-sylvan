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

/*
 * Counters for where makeedge spends its canonical form, compiled in only with
 * -DLIMDD_PROFILE (they cost nothing otherwise). They answer whether skipping
 * the canonical form for a node that no stored node can be isomorphic to would
 * pay: how often makeedge takes each path, how often the children's groups are
 * trivial, how often the coset search runs its elimination, and -- by keeping
 * fingerprints of what it has seen -- how often a node's isomorphism class key
 * (var, {low target, high target}) was new when it was made.
 *
 * One cache line of counters per worker, summed only when the report is
 * printed. The fingerprint sets are kept only with one worker, which is how the
 * benchmarks run; with more they are not kept and the report says so.
 */

#ifndef QSYLVAN_LIMDD_PROFILE_H
#define QSYLVAN_LIMDD_PROFILE_H

#include <stdio.h>

#ifdef LIMDD_PROFILE

#include <sylvan_int.h>

typedef enum {
    LP_CALLS,          /* makeedge_canonical calls */
    LP_ZERO,           /* both branches zero */
    LP_DEAD,           /* one branch zero */
    LP_EQ,             /* low == high: the skip rule without a group */
    LP_SAME,           /* one child node, two edges (the swap and the skip rule need groups) */
    LP_GEN,            /* two different child nodes */
    LP_GEN_TRIV,       /* ... with both child groups trivial */
    LP_GEN_CREATED,    /* ... and the node was new */
    LP_GEN_EXACT,      /* ... found, and the factored node was already the canonical one */
    LP_GEN_MERGED,     /* ... found, but only after the canonical form changed the node */
    LP_LZ_TSEEN,       /* lazy estimate: this factored node was asked for before */
    LP_LZ_CLASSNEW,    /* lazy estimate: first node ever with this class key */
    LP_LZ_COLLIDE,     /* lazy estimate: class key seen before, factored node new */
    LP_LZ_COLLIDE_FOUND, /* ... and the canonical form found an existing node */
    LP_MINCOSET,       /* limdd_stab_min_coset calls */
    LP_MINCOSET_HIT,   /* ... answered from the operation cache */
    LP_ZASSENHAUS,     /* stack_generators runs (the O(n^3) elimination) */
    LP_STAB_OF_NODE,   /* limdd_stab_of_node runs */
    LP_LZ_KNOWN,       /* lazy form: a stored factored node whose class was known */
    LP_LZ_RESOLVE,     /* lazy form: class resolutions */
    LP_LZ_MEMBER,      /* ... other members met in the class's slots */
    LP_LZ_FILTERED,    /* ... of which the modulus filter set aside */
    LP_LZ_ALIAS,       /* ... resolutions that found an image already stored */
    LP_N
} limdd_prof_t;

typedef struct {
    uint64_t c[LP_N];
    char pad[SYLVAN_SHARING_PAD - (LP_N * sizeof(uint64_t)) % SYLVAN_SHARING_PAD];
} limdd_prof_line_t;

extern limdd_prof_line_t *limdd_prof;

static inline void
limdd_prof_count(limdd_prof_t k)
{
    if (limdd_prof == NULL) return;
    WorkerP *w = lace_get_worker();
    limdd_prof[w ? (unsigned)w->worker : 0].c[k]++;
}

#define LIMDD_PROF(k) limdd_prof_count(k)

#else

#define LIMDD_PROF(k) ((void)0)

#endif

/** Allocates the counters (once Lace runs); a no-op without LIMDD_PROFILE. */
void limdd_prof_init(void);

/** Prints the summed counters as LIMPROF lines; a no-op without LIMDD_PROFILE. */
void limdd_prof_report(FILE *out);

#endif
