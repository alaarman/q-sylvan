/*
 * Copyright 2025 Q-Sylvan contributors
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

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "qsylvan_limdd_inspect.h"
#include "qsylvan_limdd_lim.h"
#include "sylvan_edge_weights.h"
#include "edge_weight_storage/qisq2_map.h"

/*
 * A visited set of its own, rather than llmsset's mark bitmap.
 *
 * The collector's bitmap is the only marking the node table has, and borrowing
 * it here would make looking at a diagram destroy the information a concurrent
 * collection depends on. This is a plain open-addressed set, local to the
 * call, so inspecting a diagram cannot disturb anything.
 */
typedef struct {
    LIMDD_TARG *slot;
    size_t      mask;
    size_t      fill;
} seen_t;

static bool
seen_init(seen_t *s, size_t hint)
{
    size_t cap = 64;
    while (cap < hint * 2) cap <<= 1;
    s->slot = calloc(cap, sizeof(LIMDD_TARG));
    if (s->slot == NULL) return false;
    s->mask = cap - 1;
    s->fill = 0;
    return true;
}

static void seen_free(seen_t *s) { free(s->slot); s->slot = NULL; }

static bool seen_grow(seen_t *s);

/** True if `t` was already present; inserts it otherwise. */
static bool
seen_add(seen_t *s, LIMDD_TARG t)
{
    assert(t != 0 && "0 is the empty slot, and never a real target");
    for (;;) {
        size_t i = (size_t)(t * UINT64_C(0x9E3779B97F4A7C15)) & s->mask;
        for (size_t probe = 0; probe <= s->mask; probe++, i = (i + 1) & s->mask) {
            if (s->slot[i] == t) return true;
            if (s->slot[i] == 0) {
                if ((s->fill + 1) * 4 > (s->mask + 1) * 3) break;  /* over 3/4: grow */
                s->slot[i] = t;
                s->fill++;
                return false;
            }
        }
        if (!seen_grow(s)) {
            fprintf(stderr, "sylvan: out of memory while inspecting a LIMDD\n");
            exit(1);
        }
    }
}

static bool
seen_grow(seen_t *s)
{
    const size_t oldcap = s->mask + 1;
    LIMDD_TARG *old = s->slot;
    LIMDD_TARG *neu = calloc(oldcap * 2, sizeof(LIMDD_TARG));
    if (neu == NULL) return false;
    s->slot = neu;
    s->mask = oldcap * 2 - 1;
    s->fill = 0;
    for (size_t i = 0; i < oldcap; i++)
        if (old[i] != 0) seen_add(s, old[i]);
    free(old);
    return true;
}

/* --- the walk ------------------------------------------------------------ */

static void
walk(LIMDD e, seen_t *seen, size_t *counts, size_t nqubits)
{
    if (limdd_edge_is_zero(e)) return;
    const LIMDD_TARG t = limdd_target(e);
    if (t == LIMDD_TERMINAL) return;
    if (seen_add(seen, t)) return;

    const uint32_t var = limdd_node_var(t);
    if (counts != NULL && var < nqubits) counts[var]++;

    walk(limdd_node_low(t), seen, counts, nqubits);
    walk(limdd_node_high(t), seen, counts, nqubits);
}

size_t
limdd_level_counts(LIMDD e, size_t *counts, size_t nqubits)
{
    if (counts != NULL) memset(counts, 0, nqubits * sizeof(*counts));
    seen_t seen;
    if (!seen_init(&seen, 1024)) {
        fprintf(stderr, "sylvan: out of memory while inspecting a LIMDD\n");
        exit(1);
    }
    walk(e, &seen, counts, nqubits);
    const size_t n = seen.fill;
    seen_free(&seen);
    return n;
}

size_t
limdd_width(LIMDD e, size_t nqubits)
{
    size_t *counts = calloc(nqubits ? nqubits : 1, sizeof(size_t));
    if (counts == NULL) {
        fprintf(stderr, "sylvan: out of memory while inspecting a LIMDD\n");
        exit(1);
    }
    (void)limdd_level_counts(e, counts, nqubits);
    size_t w = 0;
    for (size_t i = 0; i < nqubits; i++) if (counts[i] > w) w = counts[i];
    free(counts);
    return w;
}

size_t
limdd_nodecount(LIMDD e, size_t nqubits)
{
    return limdd_level_counts(e, NULL, nqubits);
}

/* --- algebraic bit size -------------------------------------------------- */

static uint64_t
bits_of(LIMDD_LIM lim)
{
    if (limdd_lim_is_zero(lim)) return 0;
    return qisq2_size(wgt_storage, (uint64_t)limdd_lim_weight(lim));
}

static uint64_t
bits_walk(LIMDD e, seen_t *seen)
{
    if (limdd_edge_is_zero(e)) return 0;
    uint64_t best = bits_of(limdd_label(e));
    const LIMDD_TARG t = limdd_target(e);
    if (t == LIMDD_TERMINAL) return best;
    if (seen_add(seen, t)) return best;
    const uint64_t lo = bits_walk(limdd_node_low(t), seen);
    const uint64_t hi = bits_walk(limdd_node_high(t), seen);
    if (lo > best) best = lo;
    if (hi > best) best = hi;
    return best;
}

uint64_t
limdd_max_wgt_bits(LIMDD e, size_t nqubits)
{
    (void)nqubits;
    if (sylvan_get_edge_weight_type() != WGT_QISQ2) return 0;
    seen_t seen;
    if (!seen_init(&seen, 1024)) {
        fprintf(stderr, "sylvan: out of memory while inspecting a LIMDD\n");
        exit(1);
    }
    const uint64_t res = bits_walk(e, &seen);
    seen_free(&seen);
    return res;
}

/* --- dot ----------------------------------------------------------------- */

/** One arrow, with its LIM as the label. A zero edge is simply not drawn. */
static void
dot_arrow(FILE *out, LIMDD_TARG from, LIMDD e, bool high)
{
    if (limdd_edge_is_zero(e)) return;
    const LIMDD_TARG t = limdd_target(e);
    fprintf(out, "  n%llu -> ", (unsigned long long)from);
    if (t == LIMDD_TERMINAL) fprintf(out, "term");
    else                     fprintf(out, "n%llu", (unsigned long long)t);
    fprintf(out, " [style=%s, label=\"", high ? "solid" : "dashed");
    limdd_lim_fprint(out, limdd_label(e));
    fprintf(out, "\"];\n");
}

static void
dot_walk(FILE *out, LIMDD e, seen_t *seen, size_t nqubits)
{
    if (limdd_edge_is_zero(e)) return;
    const LIMDD_TARG t = limdd_target(e);
    if (t == LIMDD_TERMINAL) return;
    if (seen_add(seen, t)) return;

    fprintf(out, "  n%llu [shape=circle, label=\"q%u\"];\n",
            (unsigned long long)t, limdd_node_var(t));

    const LIMDD lo = limdd_node_low(t), hi = limdd_node_high(t);
    dot_arrow(out, t, lo, false);
    dot_arrow(out, t, hi, true);
    dot_walk(out, lo, seen, nqubits);
    dot_walk(out, hi, seen, nqubits);
}

void
limdd_fprintdot(FILE *out, LIMDD e, size_t nqubits)
{
    fprintf(out, "digraph limdd {\n");
    fprintf(out, "  rankdir=TB;\n  node [fontsize=10];\n  edge [fontsize=8];\n");

    if (limdd_edge_is_zero(e)) {
        fprintf(out, "  root [shape=point];\n  zero [shape=box, label=\"0\"];\n");
        fprintf(out, "  root -> zero;\n}\n");
        return;
    }

    fprintf(out, "  root [shape=point];\n");
    fprintf(out, "  root -> %s", limdd_target(e) == LIMDD_TERMINAL ? "term" : "n");
    if (limdd_target(e) != LIMDD_TERMINAL)
        fprintf(out, "%llu", (unsigned long long)limdd_target(e));
    fprintf(out, " [label=\"");
    limdd_lim_fprint(out, limdd_label(e));           /* the root label */
    fprintf(out, "\"];\n");

    seen_t seen;
    if (!seen_init(&seen, 1024)) {
        fprintf(stderr, "sylvan: out of memory while inspecting a LIMDD\n");
        exit(1);
    }
    dot_walk(out, e, &seen, nqubits);
    seen_free(&seen);

    fprintf(out, "  term [shape=box, label=\"1\"];\n");
    fprintf(out, "}\n");
}

int
limdd_writedot(const char *filename, LIMDD e, size_t nqubits)
{
    FILE *f = fopen(filename, "w");
    if (f == NULL) return -1;
    limdd_fprintdot(f, e, nqubits);
    return fclose(f);
}
