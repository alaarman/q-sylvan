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

/**
 * Looking at a diagram rather than computing with one.
 *
 * Width and per-level counts exist because the theory is stated in those
 * terms -- a stabiliser state is a Tower-LIMDD of width one, and a state with
 * m independent Pauli stabilisers has width at most 2^(n-m) -- and a node
 * TOTAL cannot express either. Nothing here is on an operation's path.
 */

#ifndef QSYLVAN_LIMDD_INSPECT_H
#define QSYLVAN_LIMDD_INSPECT_H

#include <stdio.h>
#include "qsylvan_limdd_node.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Reachable nodes per level, written into counts[0..nqubits-1].
 *
 * A node sits at the level of its variable, so counts[k] is how many distinct
 * nodes the diagram has at qubit k. Skipped levels are genuinely empty: a
 * reduced LIMDD has no node there, and that is the point of the reduction, so
 * a zero is a fact about the diagram and not a gap in the walk. The terminal
 * is not counted.
 *
 * Returns the total, which is the sum of the counts.
 */
size_t limdd_level_counts(LIMDD e, size_t *counts, size_t nqubits);

/** The largest per-level count: the diagram's width. */
size_t limdd_width(LIMDD e, size_t nqubits);

/** Distinct reachable nodes below `e`, excluding the terminal. */
size_t limdd_nodecount(LIMDD e, size_t nqubits);

/**
 * Write the diagram as Graphviz dot.
 *
 * Edges carry their LIM, low edges dashed and high edges solid, and each node
 * is labelled with its qubit. Reading a tower off this is the quickest way to
 * see whether a state came out in the form the theory says it should.
 */
void limdd_fprintdot(FILE *out, LIMDD e, size_t nqubits);

/** limdd_fprintdot to a file, by name. Returns 0 on success. */
int limdd_writedot(const char *filename, LIMDD e, size_t nqubits);

#ifdef __cplusplus
}
#endif

#endif
