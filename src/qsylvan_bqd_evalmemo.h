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

/**
 * The memo of one bqd_eval call, private to the three eval_node walks
 * (qsylvan_bqd.c, qsylvan_bqd_sm.c, qsylvan_bqd_xp_sm.c).
 *
 * A quotient node's high cofactor at a point is the product of both its
 * children there, so the walk forks at every such node, and without a memo it
 * visits a node once per path to it: 2^n calls on a state of full support.
 * Keyed on (node, point), each pair is evaluated once, so a call costs at most
 * a visit per node and point; in the scalar family the point at a node is the
 * suffix of x below its level, one per node. The memo lives for one call
 * only, in the caller's frame and its own heap block, so it shares nothing.
 */

#ifndef QSYLVAN_BQD_EVALMEMO_H
#define QSYLVAN_BQD_EVALMEMO_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "sylvan_edge_weights.h"

typedef struct {
    uint64_t  *node;    /* 0 marks an empty slot: node indices start at 2 */
    uint64_t  *point;
    EVBDD_WGT *value;
    size_t     cap;     /* a power of two */
    size_t     used;
} bqd_eval_memo_t;

static inline void
bqd_evm_alloc(bqd_eval_memo_t *m, size_t cap)
{
    m->node = (uint64_t *)calloc(cap, sizeof(uint64_t));
    m->point = (uint64_t *)malloc(cap * sizeof(uint64_t));
    m->value = (EVBDD_WGT *)malloc(cap * sizeof(EVBDD_WGT));
    if (m->node == NULL || m->point == NULL || m->value == NULL) {
        fprintf(stderr, "sylvan: out of memory in bqd_eval\n");
        exit(1);
    }
    m->cap = cap;
    m->used = 0;
}

static inline void
bqd_evm_init(bqd_eval_memo_t *m)
{
    bqd_evm_alloc(m, 256);
}

static inline void
bqd_evm_free(bqd_eval_memo_t *m)
{
    free(m->node); free(m->point); free(m->value);
    m->node = m->point = NULL; m->value = NULL;
}

static inline size_t
bqd_evm_slot(const bqd_eval_memo_t *m, uint64_t node, uint64_t point)
{
    uint64_t h = node * UINT64_C(0x9E3779B97F4A7C15) ^ (point + UINT64_C(0x632BE59BD9B4E019)) * UINT64_C(0xC2B2AE3D27D4EB4F);
    h ^= h >> 31;
    size_t i = (size_t)h & (m->cap - 1);
    while (m->node[i] != 0 && (m->node[i] != node || m->point[i] != point)) i = (i + 1) & (m->cap - 1);
    return i;
}

static inline bool
bqd_evm_get(const bqd_eval_memo_t *m, uint64_t node, uint64_t point, EVBDD_WGT *v)
{
    const size_t i = bqd_evm_slot(m, node, point);
    if (m->node[i] == 0) return false;
    *v = m->value[i];
    return true;
}

static inline void
bqd_evm_put(bqd_eval_memo_t *m, uint64_t node, uint64_t point, EVBDD_WGT v)
{
    if (2 * (m->used + 1) > m->cap) {
        bqd_eval_memo_t g;
        bqd_evm_alloc(&g, 2 * m->cap);
        for (size_t i = 0; i < m->cap; i++) {
            if (m->node[i] == 0) continue;
            const size_t j = bqd_evm_slot(&g, m->node[i], m->point[i]);
            g.node[j] = m->node[i]; g.point[j] = m->point[i]; g.value[j] = m->value[i];
            g.used++;
        }
        bqd_evm_free(m);
        *m = g;
    }
    const size_t i = bqd_evm_slot(m, node, point);
    if (m->node[i] == 0) m->used++;
    m->node[i] = node; m->point[i] = point; m->value[i] = v;
}

#endif
