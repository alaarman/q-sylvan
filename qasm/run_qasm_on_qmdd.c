/**
 * Copyright 2024 System Verification Lab, LIACS, Leiden University
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
 * 
 */

#include <inttypes.h>
#include <argp.h>
#include <sys/time.h>

#include "qsylvan.h"
#include "qsylvan_limdd_canon.h"
#include "qsylvan_limdd_ops.h"
#include "qsylvan_limdd_gc.h"
#include "qsylvan_limdd_inspect.h"
#include "qsylvan_bqd.h"
#include "qsylvan_bqd_ops.h"
#include "qsylvan_bqd_gates.h"
#include "qsylvan_qasm_parser.h"

/**********************<Arguments (configured via argp)>***********************/

static int workers = 1;
static int rseed = 0;
static bool count_nodes = false;
static bool count_qisq2_size = false;
static bool calc_measurement_prob = false;
static bool output_vector = false;
static size_t min_tablesize = 1LL<<25;
static size_t max_tablesize = 1LL<<25;
static size_t min_cachesize = 1LL<<16;
static size_t max_cachesize = 1LL<<16;
static bool cache_size_set = false;
static size_t min_wgt_tab_size = 1LL<<23;
static size_t max_wgt_tab_size = 1LL<<23;
static double tolerance = 1e-14;
static int wgt_table_type = COMP_HASHMAP;
static int wgt_norm_strat = NORM_L2;
static bool wgt_inv_caching = true;
/*
 * The merging rule: when are two edge weights the same table entry. Below
 * zero keeps the historical absolute rule; a relative tolerance selects the
 * hybrid rule, which compares relative to magnitude and collapses anything
 * smaller than zero_tolerance. See cmap.c.
 */
static int canon_policy = -1;        /* -1: follow --dd */
static uint64_t canon_interval = 100000;
static const char *canon_name = "always";
static double rel_tolerance = -1;
static double zero_tolerance = 1e-14;
static bool zero_tolerance_set = false;   /* did the user ask for one? */
static bool force_absolute = false;       /* --merging=abs: keep the historical rule */
static bool node_tab_size_set = false;    /* was --node-tab-size given? */
static int  lim_tab_size_log2  = 0;       /* --lim-tab-size; 0 = derive from the node table */
static bool lim_stats = false;
static char *trace_file = NULL;   /* --trace: per-gate telemetry */
static FILE *trace_out = NULL;
typedef enum { DD_QMDD, DD_LIMDD, DD_LIMDD_HEUR, DD_BQD } dd_kind_t;
static dd_kind_t dd_kind = DD_QMDD;
static const char *dd_kind_name = "qmdd";
static int reorder_qubits = 0;
static char* qasm_inputfile = NULL;
static char* json_outputfile = NULL;
/*
 * -v on the LIMDD and BQD arms: the amplitudes, in the order fprint_stats
 * prints them, taken before the diagram's tables are freed (up to
 * DECODE_MAX_QUBITS). The QMDD arm reads its state when it prints.
 */
#define DECODE_MAX_QUBITS 24     /* 2^24 weights, 128 MB */
static EVBDD_WGT *final_vector = NULL;


static struct argp_option options[] =
{
    {"workers", 'w', "<workers>", 0, "Number of workers/threads (default=1)", 0},
    {"rseed", 'r', "<random-seed>", 0, "Set random seed", 0},
    {"norm-strat", 's', "<low|max|min|l2>", 0, "Edge weight normalization strategy", 0},
    {"edge-weight-type", 'e', "<float|qisq2>", 0, "Edge weight type (default float)", 0},
    {"tol", 't', "<tolerance>", 0, "Tolerance for deciding edge weights equal (default=1e-14)", 0},
    {"json", 'j', "<filename>", 0, "Write stats to given filename as json", 0},
    {"count-nodes", 'c', 0, 0, "Track maximum number of nodes", 0},
    {"cache-size", 1010, "<size>", 0, "log2 of the operation cache size (default 16)", 0},
    {"trace", 1013, "<filename>", 0, "Write per-gate telemetry as CSV: after every gate, the running T count, the node count, the LIMDD width and the largest algebraic weight in bits. This is what relates the theorems to a run -- width against 2^t, and coefficient size against the bound of Section 5.", 0},
    {"lim-stats", 1009, 0, 0, "For limdd: report the Pauli support of the high-edge LIM per level, to size an inline encoding", 0},
    {"count-qisq-size", 'q', 0, 0, "Count the number of bits of the largest qisq value", 0},
    {"calc-measurement-prob", 'm', 0, 0, "Calculate the probability on a specific outcome of the final state", 0},
    {"state-vector", 'v', 0, 0, "Also output the complete state vector", 0},
    {"node-tab-size", 1000, "<size>", 0, "log2 of max node table size (max 40)", 0},
    {"lim-tab-size", 1012, "<size>", 0, "LIMDD only: log2 of the LIM table size. Default is four times the node table. A LIMDD mints far more labels than it keeps, so this is what fills first on wide circuits; it has no QMDD counterpart, so raising it does not change how a QMDD is resourced.", 0},
    {"wgt-tab-size", 1001, "<size>", 0, "log2 of max edge weigth table size (max 30 (23 if node table >2^30))", 0},
    {"reorder", 1002, 0, 0, "Reorders the qubits once such that (most) controls occur before targets in the variable order.", 0},
    {"reorder-swaps", 1003, 0, 0, "Reorders the qubits such that all controls occur before targets (requires inserting SWAP gates).", 0},
    {"disable-inv-caching", 1004, 0, 0, "Disable storing inverse of MUL and DIV in cache.", 0},
    {"canon", 1008, "<always|never|ops:K|adaptive>", 0, "LIMDD only: when to apply the canonical form. always (default) canonicalises inside every operation; never leaves it off; ops:K rebuilds after every K node-building operations; adaptive tunes the interval by how much the last rebuild helped.", 0},
    {"rel-tol", 1006, "<tolerance>", 0, "Relative tolerance; selects the hybrid merging rule (default: off for qmdd, --tol for limdd, which needs it)", 0},
    {"zero-tol", 1007, "<tolerance>", 0, "Zero-collapse tolerance for the hybrid merging rule (default: 1e-14 for qmdd, 0 for limdd, whose weights are legitimately tiny)", 0},
    {"merging", 1011, "<abs|hybrid>", 0, "Which merging rule to use, overriding the per-diagram default. abs is the historical single absolute tolerance (--tol); hybrid is relative plus zero-collapse. A LIMDD defaults to hybrid and needs --merging=abs to be held to the absolute rule; exact (qisq2) weights ignore both.", 0},
    {"dd", 'd', "<qmdd|limdd|limdd-heur|bqd>", 0, "Decision diagram to simulate with (default qmdd). limdd applies the full canonical form; limdd-heur skips the search for a canonical high-edge label and only divides the low label out, which is cheaper per node but stops nodes that are the same state up to a LIM from merging. bqd is the binary quotient diagram with scalar labels, on the LIMDD's gate set: the diagonal gates by the paper's O(n) algorithm while the state has full support, everything else by a recursion that is exact and has no size bound. It needs -e qisq2 once a gate can cancel amplitudes. It has no norm or measurement algorithm, so -m and -v decode the state when it has at most 24 qubits and report -1 otherwise.", 0},
    {0, 0, 0, 0, 0, 0}
};

static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
    switch (key) {
    case 'w':
        workers = atoi(arg);
        break;
    case 'r':
        rseed = atoi(arg);
        break;
    case 's':
        if (strcmp(arg, "low")==0) wgt_norm_strat = NORM_LOW;
        else if (strcmp(arg, "max")==0) wgt_norm_strat = NORM_MAX;
        else if (strcmp(arg, "min")==0) wgt_norm_strat = NORM_MIN;
        else if (strcasecmp(arg, "l2")==0) wgt_norm_strat = NORM_L2;
        else argp_usage(state);
        break;
    case 'e':
        if (strcmp(arg, "float")==0) wgt_table_type = COMP_HASHMAP;
        else if (strcasecmp(arg, "qisq2")==0) wgt_table_type = QISQ2_MAP;
        else argp_usage(state);
        break;
    case 't':
        tolerance = atof(arg);
        break;
    case 'j':
        json_outputfile = arg;
        break;
    case 'c':
        count_nodes = true;
        break;
    case 1009:
        lim_stats = true;
        break;
    case 1013:
        trace_file = arg;
        break;
    case 1012:
        if (atoi(arg) > 40) argp_usage(state);
        lim_tab_size_log2 = atoi(arg);
        break;
    case 1011:
        if (strcasecmp(arg, "abs")==0 || strcasecmp(arg, "absolute")==0) force_absolute = true;
        else if (strcasecmp(arg, "hybrid")==0) force_absolute = false;
        else argp_usage(state);
        break;
    case 1010:
        if (atoi(arg) > 30) argp_usage(state);
        min_cachesize = max_cachesize = 1LL<<(atoi(arg));
        cache_size_set = true;
        break;
    case 'q':
        count_qisq2_size = true;
        break;
    case 'm':
        calc_measurement_prob = true;
        break;
    case 'v':
        output_vector = true;
        break;
    case 1000:
        if (atoi(arg) > 40) argp_usage(state);
        max_tablesize = 1LL<<(atoi(arg));
        node_tab_size_set = true;
        break;
    case 1001:
        if (atoi(arg) > 30) argp_usage(state);
        max_wgt_tab_size = 1LL<<(atoi(arg));
        break;
    case 1002:
        reorder_qubits = 1;
        break;
    case 1003:
        reorder_qubits = 2;
        break;
    case 1004:
        wgt_inv_caching = false;
        break;
    case 1008:
        canon_name = arg;
        if (strcmp(arg, "always") == 0) canon_policy = LIMDD_CANON_ALWAYS;
        else if (strcmp(arg, "never") == 0) canon_policy = LIMDD_CANON_MANUAL;
        else if (strcmp(arg, "adaptive") == 0) canon_policy = LIMDD_CANON_ADAPTIVE;
        else if (strncmp(arg, "ops:", 4) == 0) {
            canon_policy = LIMDD_CANON_OPS;
            canon_interval = (uint64_t) atoll(arg + 4);
            if (canon_interval == 0) argp_error(state, "--canon=ops:K needs K > 0");
        }
        else argp_error(state, "unknown canonization policy '%s'", arg);
        break;
    case 1006:
        rel_tolerance = atof(arg);
        break;
    case 1007:
        zero_tolerance = atof(arg);
        zero_tolerance_set = true;
        break;
    case 'd':
        dd_kind_name = arg;
        if (strcmp(arg, "qmdd") == 0)             dd_kind = DD_QMDD;
        else if (strcmp(arg, "limdd") == 0)       dd_kind = DD_LIMDD;
        else if (strcmp(arg, "limdd-heur") == 0)  dd_kind = DD_LIMDD_HEUR;
        else if (strcmp(arg, "bqd") == 0)         dd_kind = DD_BQD;
        else argp_error(state, "unknown dd type '%s'", arg);
        break;
    case ARGP_KEY_ARG:
        if (state->arg_num >= 1) argp_usage(state);
        qasm_inputfile = arg;
        break;
    case ARGP_KEY_END:
        if (state->arg_num < 1) argp_usage(state);
        break;
    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}
static struct argp argp = { options, parse_opt, "<qasm_file>", 0, 0, 0, 0 };

/*********************</Arguments (configured via argp)>***********************/


// NOTE: If there are only measurements at the end of the circuit, 'final_nodes' 
// and 'norm' will contain the node count and the norm of the state QMDD before
// the measurements.
typedef struct stats_s {
    uint64_t applied_gates;
    uint64_t final_nodes;
    uint64_t max_nodes;
    uint64_t t_count;        /* T and Tdg gates applied */
    uint64_t final_width;    /* LIMDD only */
    uint64_t max_width;
    uint64_t final_qisq_size;
    uint64_t max_qisq_size;
    uint64_t shots;
    double simulation_time;
    double norm;
    double normed_prob;
    double unnormed_prob;
    double first_qubit_prob;
    QMDD final_state;
} stats_t;
stats_t stats;


void fprint_stats(FILE *stream, quantum_circuit_t* circuit)
{
    fprintf(stream, "{\n");
    fprintf(stream, "  \"measurement_results\": {\n");
    fprintf(stream, "    \""); fprint_creg(stream, circuit); fprintf(stream, "\": 1\n");
    fprintf(stream, "  },\n");
    if (output_vector && (dd_kind == DD_QMDD || final_vector != NULL))
    {
        fprintf(stream, "  \"state_vector\": [\n");
        for (int k = 0; k < (1<<(circuit->qreg_size)); k++) {
            bool *x = int_to_bitarray(k, circuit->qreg_size, !(circuit->reversed_qubit_order));
            const complex_t c = (dd_kind == DD_QMDD)
                ? qmdd_get_amplitude(stats.final_state, x, circuit->qreg_size)
                : weight_as_complex(final_vector[k]);
            fprintf(stream, "    [\n");
            fprintf(stream, "      %.16lf,\n", c.r);
            fprintf(stream, "      %.16lf\n", c.i);
            if (k == (1<<(circuit->qreg_size))-1)
                fprintf(stream, "    ]\n");
            else
                fprintf(stream, "    ],\n");
            free(x);
        }
        fprintf(stream, "  ],\n");
    }
    fprintf(stream, "  \"statistics\": {\n");
    fprintf(stream, "    \"applied_gates\": %" PRIu64 ",\n", stats.applied_gates);
    fprintf(stream, "    \"benchmark\": \"%s\",\n", circuit->name);
    fprintf(stream, "    \"final_nodes\": %" PRIu64 ",\n", stats.final_nodes);
    fprintf(stream, "    \"max_nodes\": %" PRIu64 ",\n", stats.max_nodes);
    fprintf(stream, "    \"t_count\": %" PRIu64 ",\n", stats.t_count);
    fprintf(stream, "    \"final_width\": %" PRIu64 ",\n", stats.final_width);
    fprintf(stream, "    \"max_width\": %" PRIu64 ",\n", stats.max_width);
    fprintf(stream, "    \"final_qisq_size\": %" PRIu64 ",\n", stats.final_qisq_size);
    fprintf(stream, "    \"max_qisq_size\": %" PRIu64 ",\n", stats.max_qisq_size);
    fprintf(stream, "    \"n_qubits\": %d,\n", circuit->qreg_size);
    fprintf(stream, "    \"norm\": %.5e,\n", stats.norm);
    fprintf(stream, "    \"unnormed_measurement_prob\": %.5e,\n", stats.unnormed_prob);
    fprintf(stream, "    \"normed_measurement_prob\": %.5e,\n", stats.normed_prob);
    fprintf(stream, "    \"first_qubit_measurement_prob\": %.5e,\n", stats.first_qubit_prob);
    fprintf(stream, "    \"reorder\": %d,\n", reorder_qubits);
    fprintf(stream, "    \"seed\": %d,\n", rseed);
    fprintf(stream, "    \"shots\": %" PRIu64 ",\n", stats.shots);
    fprintf(stream, "    \"simulation_time\": %lf,\n", stats.simulation_time);
    fprintf(stream, "    \"tolerance\": %.5e,\n", tolerance);
    fprintf(stream, "    \"wgt_inv_caching\": %d,\n", wgt_inv_caching);
    fprintf(stream, "    \"dd\": \"%s\",\n", dd_kind_name);
    fprintf(stream, "    \"canon\": \"%s\",\n", canon_name);
    fprintf(stream, "    \"merging_rule\": \"%s\",\n",
            rel_tolerance >= 0 ? "hybrid" : "absolute");
    fprintf(stream, "    \"rel_tol\": %.5e,\n", rel_tolerance);
    fprintf(stream, "    \"zero_tol\": %.5e,\n", zero_tolerance);
    fprintf(stream, "    \"wgt_norm_strat\": %d,\n", wgt_norm_strat);
    fprintf(stream, "    \"wgt_type\": %d,\n", wgt_table_type);
    fprintf(stream, "    \"min_node_tab_size\": %" PRId64 ",\n", min_tablesize);
    fprintf(stream, "    \"max_node_tab_size\": %" PRId64 ",\n", max_tablesize);
    fprintf(stream, "    \"min_wgt_tab_size\": %" PRId64 ",\n", min_wgt_tab_size);
    fprintf(stream, "    \"max_wgt_tab_size\": %" PRId64 ",\n", max_wgt_tab_size);
    fprintf(stream, "    \"workers\": %d\n", workers);
    fprintf(stream, "  }\n");
    fprintf(stream, "}\n");
}

static double
wctime()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec + 1E-6 * tv.tv_usec);
}

/**
 * Here we match the name of a gate in QASM to
 * the GATEID 
 */
QMDD apply_gate(QMDD state, quantum_op_t* gate, BDDVAR nqubits)
{
    // TODO: move this relation between parsed quantum_op and internal gate
    // somewhere else?
    stats.applied_gates++;

  if (strcmp(gate->name, "id") == 0) {
        stats.applied_gates--;
        return state;
    }
    else if (strcmp(gate->name, "x") == 0) {
        return qmdd_gate(state, GATEID_X, gate->targets[0]);
    }
    else if (strcmp(gate->name, "y") == 0) {
        return qmdd_gate(state, GATEID_Y, gate->targets[0]);
    }
    else if (strcmp(gate->name, "z") == 0) {
        return qmdd_gate(state, GATEID_Z, gate->targets[0]);
    }
    else if (strcmp(gate->name, "h") == 0) {
        return qmdd_gate(state, GATEID_H, gate->targets[0]);
    }
    else if (strcmp(gate->name, "s") == 0) {
        return qmdd_gate(state, GATEID_S, gate->targets[0]);
    }
    else if (strcmp(gate->name, "sdg") == 0) {
        return qmdd_gate(state, GATEID_Sdag, gate->targets[0]);
    }
    else if (strcmp(gate->name, "t") == 0) {
        stats.t_count++;   /* counted on both paths: it is a property of the
                            * circuit, and the bounds are stated in terms of it */
        return qmdd_gate(state, GATEID_T, gate->targets[0]);
    }
    else if (strcmp(gate->name, "tdg") == 0) {
        stats.t_count++;
        return qmdd_gate(state, GATEID_Tdag, gate->targets[0]);
    }
    else if (strcmp(gate->name, "sx") == 0) {
        return qmdd_gate(state, GATEID_sqrtX, gate->targets[0]);
    }
    else if (strcmp(gate->name, "sxdg") == 0) {
        return qmdd_gate(state, GATEID_sqrtXdag, gate->targets[0]);
    }
    else if (strcmp(gate->name, "rx") == 0) {
        return qmdd_gate(state, GATEID_Rx(gate->angle[0]), gate->targets[0]);
    }
    else if (strcmp(gate->name, "ry") == 0) {
        return qmdd_gate(state, GATEID_Ry(gate->angle[0]), gate->targets[0]);
    }
    else if (strcmp(gate->name, "rz") == 0) {
        return qmdd_gate(state, GATEID_Rz(gate->angle[0]), gate->targets[0]);
    }
    else if (strcmp(gate->name, "p") == 0) {
        return qmdd_gate(state, GATEID_Phase(gate->angle[0]), gate->targets[0]);
    }
    else if (strcmp(gate->name, "u2") == 0) {
        fl_t pi_over_2 = flt_acos(0.0);
        return qmdd_gate(state, GATEID_U(pi_over_2, gate->angle[0], gate->angle[1]), gate->targets[0]);
    }
    else if (strcmp(gate->name, "u") == 0) {
        return qmdd_gate(state, GATEID_U(gate->angle[0], gate->angle[1], gate->angle[2]), gate->targets[0]);
    }
    else if (strcmp(gate->name, "cx") == 0) {
        return qmdd_cgate(state, GATEID_X, gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "cy") == 0) {
        return qmdd_cgate(state, GATEID_Y, gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "cz") == 0) {
        return qmdd_cgate(state, GATEID_Z, gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "ch") == 0) {
        return qmdd_cgate(state, GATEID_H, gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "csx") == 0) {
        return qmdd_cgate(state, GATEID_sqrtX, gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "crx") == 0) {
        return qmdd_cgate(state, GATEID_Rx(gate->angle[0]), gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "cry") == 0) {
        return qmdd_cgate(state, GATEID_Ry(gate->angle[0]), gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "crz") == 0) {
        return qmdd_cgate(state, GATEID_Rz(gate->angle[0]), gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "cp") == 0) {
        return qmdd_cgate(state, GATEID_Phase(gate->angle[0]), gate->ctrls[0], gate->targets[0], nqubits);
    }
    /*
     * cs, csdg and ccz are level-3 diagonal gates that every arm applies as
     * one gate. Their T-count is that of the ancilla-free Clifford+T
     * decomposition, 3 for a controlled S and 7 for a CCZ, so that t_count
     * stays the T-count of the circuit and the bounds stated in it apply.
     */
    else if (strcmp(gate->name, "cs") == 0 || strcmp(gate->name, "csdg") == 0) {
        stats.t_count += 3;
        return qmdd_cgate(state, gate->name[2] == 'd' ? GATEID_Sdag : GATEID_S,
                          gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "ccz") == 0) {
        stats.t_count += 7;
        return qmdd_cgate2(state, GATEID_Z, gate->ctrls[0], gate->ctrls[1], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "cu") == 0) {
        return qmdd_cgate(state, GATEID_U(gate->angle[0], gate->angle[1], gate->angle[2]), gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "ccx") == 0) {
        return qmdd_cgate2(state, GATEID_X, gate->ctrls[0], gate->ctrls[1], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "c3x") == 0) {
        return qmdd_cgate3(state, GATEID_X, gate->ctrls[0], gate->ctrls[1], gate->ctrls[2], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "c3sx") == 0) {
        return qmdd_cgate3(state, GATEID_sqrtX, gate->ctrls[0], gate->ctrls[1], gate->ctrls[2], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "swap") == 0) {
        // no native SWAP gates in Q-Sylvan
        stats.applied_gates += 4;
        return qmdd_circuit_swap(state, gate->targets[0], gate->targets[1]);
    }
    else if (strcmp(gate->name, "cswap") == 0) {
        // no native CSWAP gates in Q-Sylvan
        stats.applied_gates += 4;
        // CCNOT
        state = qmdd_cgate2(state, GATEID_X, gate->ctrls[0], gate->targets[0], gate->targets[1], nqubits);
        // upside down CCNOT (equivalent)
        state = qmdd_cgate(state, GATEID_H, gate->ctrls[0], gate->targets[0], nqubits);
        state = qmdd_cgate2(state, GATEID_Z, gate->ctrls[0], gate->targets[0], gate->targets[1], nqubits);
        state = qmdd_cgate(state, GATEID_H, gate->ctrls[0], gate->targets[0], nqubits);
        // CCNOT
        state = qmdd_cgate2(state, GATEID_X, gate->ctrls[0], gate->targets[0], gate->targets[1], nqubits);

        return state;
    }
    else if (strcmp(gate->name, "rccx") == 0) {
        // no native RCCX (simplified Toffoli) gates in Q-Sylvan
        stats.applied_gates += 3;
        state = qmdd_cgate2(state, GATEID_X, gate->ctrls[0], gate->ctrls[1], gate->targets[0], nqubits);
        state = qmdd_gate(state, GATEID_X, gate->ctrls[1]);
        state = qmdd_cgate2(state, GATEID_Z, gate->ctrls[0], gate->ctrls[1], gate->targets[0], nqubits);
        state = qmdd_gate(state, GATEID_X, gate->ctrls[1]);
        return state;
    }
    else if (strcmp(gate->name, "rzz") == 0 ) {
        // no native RZZ gates in Q-Sylvan
        stats.applied_gates += 2;
        state = qmdd_cgate(state, GATEID_X, gate->targets[0], gate->targets[1], nqubits);
        state = qmdd_gate(state, GATEID_Phase(gate->angle[0]), gate->targets[1]);
        state = qmdd_cgate(state, GATEID_X, gate->targets[0], gate->targets[1], nqubits);
        return state;
    }
    else if (strcmp(gate->name, "rxx") == 0) {
        // no native RXX gates in Q-Sylvan
        fl_t pi = flt_acos(0.0) * 2;
        stats.applied_gates += 6;
        state = qmdd_gate(state, GATEID_U(pi/2.0, gate->angle[0], 0), gate->targets[0]);
        state = qmdd_gate(state, GATEID_H, gate->targets[1]);
        state = qmdd_cgate(state, GATEID_X, gate->targets[0], gate->targets[1], nqubits);
        state = qmdd_gate(state, GATEID_Phase(-(gate->angle[0])), gate->targets[1]);
        state = qmdd_cgate(state, GATEID_X, gate->targets[0], gate->targets[1], nqubits);
        state = qmdd_gate(state, GATEID_H, gate->targets[1]);
        state = qmdd_gate(state, GATEID_U(pi/2.0, -pi, pi - gate->angle[0]), gate->targets[0]);
        return state;
    }
    else {
        fprintf(stderr, "Gate '%s' currently unsupported\n", gate->name);
        return state;
    }
}


QMDD measure(QMDD state, quantum_op_t *meas, quantum_circuit_t* circuit)
{
    double p;
    int m;
    printf("measure qubit %d, store result in creg[%d]\n", meas->targets[0], meas->meas_dest);
    qmdd_measure_qubit(state, meas->targets[0], circuit->qreg_size, &m, &p);
    circuit->creg[meas->meas_dest] = m;
    return state;
}



/* --- LIMDD and BQD simulation ------------------------------------------- */

/*
 * The two diagrams are called the same way: a BQD edge is a LIMDD edge, and
 * qsylvan_bqd_gates.h has the signatures of limdd_gate, limdd_cgate,
 * limdd_cgate_either and limdd_swap. So one dispatch serves both, through
 * this table, and the two arms differ only in how they start, what they do
 * between gates and what they can report at the end.
 */
typedef struct {
    const char *name;
    LIMDD (*gate)(LIMDD e, uint32_t gateid, uint32_t target, uint32_t n);
    LIMDD (*cgate)(LIMDD e, uint32_t gateid, uint64_t cmask, uint32_t target, uint32_t n);
    LIMDD (*cgate_either)(LIMDD e, uint32_t gateid, uint32_t c, uint32_t t, uint32_t n, bool *ok);
    LIMDD (*swap)(LIMDD e, uint32_t a, uint32_t b, uint32_t n);
} dd_gate_ops_t;

static LIMDD op_limdd_gate(LIMDD e, uint32_t g, uint32_t t, uint32_t n) { return limdd_gate(e, g, t, n); }
static LIMDD op_limdd_cgate(LIMDD e, uint32_t g, uint64_t m, uint32_t t, uint32_t n) { return limdd_cgate(e, g, m, t, n); }
static LIMDD op_bqd_gate(LIMDD e, uint32_t g, uint32_t t, uint32_t n) { return bqd_gate(e, g, t, n); }
static LIMDD op_bqd_cgate(LIMDD e, uint32_t g, uint64_t m, uint32_t t, uint32_t n) { return bqd_cgate(e, g, m, t, n); }

static const dd_gate_ops_t LIMDD_GATES = { "limdd", op_limdd_gate, op_limdd_cgate, limdd_cgate_either, limdd_swap };
static const dd_gate_ops_t BQD_GATES   = { "bqd",   op_bqd_gate,   op_bqd_cgate,   bqd_cgate_either,   bqd_swap };

/**
 * The Clifford+T gate set, which is what the LIMDD and BQD paths support.
 *
 * Arbitrary-angle rotations are deliberately absent. With floating point
 * weights they would work, but the point of running LIMDDs here is to pair
 * them with exact coefficients, and e^(i*theta) is in Q[i,sqrt2] only for
 * theta a multiple of pi/4. Refusing is better than quietly producing a
 * number whose exactness is a fiction.
 *
 * Returns false and leaves *state alone if the gate is not supported.
 */
static bool
dd_apply_gate(const dd_gate_ops_t *dd, LIMDD *state, quantum_op_t *gate, BDDVAR nqubits)
{
    const uint32_t t = gate->targets[0];
    uint32_t id = 0;
    stats.applied_gates++;

    if      (strcmp(gate->name, "id")   == 0) return true;
    else if (strcmp(gate->name, "x")    == 0) id = GATEID_X;
    else if (strcmp(gate->name, "y")    == 0) id = GATEID_Y;
    else if (strcmp(gate->name, "z")    == 0) id = GATEID_Z;
    else if (strcmp(gate->name, "h")    == 0) id = GATEID_H;
    else if (strcmp(gate->name, "s")    == 0) id = GATEID_S;
    else if (strcmp(gate->name, "sdg")  == 0) id = GATEID_Sdag;
    else if (strcmp(gate->name, "t")    == 0) id = GATEID_T;
    else if (strcmp(gate->name, "tdg")  == 0) id = GATEID_Tdag;
    else if (strcmp(gate->name, "swap") == 0) {
        *state = dd->swap(*state, gate->targets[0], gate->targets[1], nqubits);
        return true;
    }
    else if (strcmp(gate->name, "cx")   == 0 || strcmp(gate->name, "cy") == 0 ||
             strcmp(gate->name, "cz")   == 0 || strcmp(gate->name, "ch") == 0) {
        const uint32_t cid = (gate->name[1] == 'x') ? GATEID_X
                           : (gate->name[1] == 'y') ? GATEID_Y
                           : (gate->name[1] == 'z') ? GATEID_Z : GATEID_H;
        bool ok;
        *state = dd->cgate_either(*state, cid, gate->ctrls[0], t, nqubits, &ok);
        if (!ok) {
            fprintf(stderr, "%s: '%s' with control %u below target %u has no "
                            "exact reordering identity\n",
                    dd->name, gate->name, gate->ctrls[0], t);
            return false;
        }
        return true;
    }
    else if (strcmp(gate->name, "cs") == 0 || strcmp(gate->name, "csdg") == 0) {
        /* symmetric, so the lower qubit is the control: a cgate needs
         * every control above the target. T-count as in apply_gate. */
        const uint32_t a = gate->ctrls[0] < t ? gate->ctrls[0] : t;
        const uint32_t b = gate->ctrls[0] < t ? t : gate->ctrls[0];
        stats.t_count += 3;
        *state = dd->cgate(*state, gate->name[2] == 'd' ? GATEID_Sdag : GATEID_S,
                           UINT64_C(1) << a, b, nqubits);
        return true;
    }
    else if (strcmp(gate->name, "ccz") == 0) {
        /* symmetric in all three: the highest-numbered qubit is the target */
        uint32_t q[3] = { (uint32_t)gate->ctrls[0], (uint32_t)gate->ctrls[1], t };
        for (int i = 0; i < 2; i++) for (int j = i + 1; j < 3; j++)
            if (q[j] < q[i]) { const uint32_t x = q[i]; q[i] = q[j]; q[j] = x; }
        stats.t_count += 7;
        *state = dd->cgate(*state, GATEID_Z, (UINT64_C(1) << q[0]) | (UINT64_C(1) << q[1]),
                           q[2], nqubits);
        return true;
    }
    else {
        fprintf(stderr, "%s: gate '%s' is outside Clifford+T\n", dd->name, gate->name);
        return false;
    }

    if (id == GATEID_T || id == GATEID_Tdag) stats.t_count++;

    *state = dd->gate(*state, id, t, nqubits);
    return true;
}

/**
 * One line of per-gate telemetry.
 *
 * Width and coefficient size both need a walk of the whole diagram, so this
 * is only done when --trace asks for it: on a circuit with a large diagram it
 * costs more than the gates do.
 */
static void
trace_row(uint64_t gate_idx, const char *gate, LIMDD state, uint32_t n)
{
    if (trace_out == NULL) return;
    const size_t nodes = limdd_nodecount(state, n);
    const size_t width = limdd_width(state, n);
    const uint64_t bits = limdd_max_wgt_bits(state, n);
    if (width > stats.max_width) stats.max_width = width;
    fprintf(trace_out, "%" PRIu64 ",%s,%" PRIu64 ",%zu,%zu,%" PRIu64 "\n",
            gate_idx, gate, stats.t_count, nodes, width, bits);
}

TASK_1(int, limdd_simulate_circuit, quantum_circuit_t*, circuit)
{
    const BDDVAR n = circuit->qreg_size;
    if (n > LIMDD_MAX_QUBITS) {
        fprintf(stderr, "limdd: %u qubits, but a Pauli word holds one bit per "
                        "qubit per component so the limit is %d\n",
                n, LIMDD_MAX_QUBITS);
        return 1;
    }

    const double t_start = wctime();
    LIMDD state = limdd_all_zero_state(n);
    limdd_protect(&state);
    uint64_t gate_idx = 0;
    if (trace_out != NULL)
        fprintf(trace_out, "gate_index,gate,t_count,nodes,width,wgt_bits\n");

    for (quantum_op_t *op = circuit->operations; op != NULL; op = op->next) {
        if (op->type == op_gate) {
            if (!dd_apply_gate(&LIMDD_GATES, &state, op, n)) { limdd_unprotect(&state); return 1; }
        }
        else if (op->type == op_measurement) {
            break;   /* the probability is taken from the final state below */
        }
        /* Between operations, never inside one: limdd_canonize rebuilds
         * nodes bottom-up and needs a root that is not half-built. */
        if (limdd_canon_due()) state = limdd_canonize(state);

        if (count_nodes) {
            const uint64_t c = limdd_countnodes(state);
            if (c > stats.max_nodes) stats.max_nodes = c;
        }

        if (trace_out != NULL) {
            trace_row(gate_idx, op->type == op_gate ? op->name : "-", state, n);
        }
        gate_idx++;

        /*
         * Nothing is freed until it is collected, and a few hundred gates
         * intern far more LIMs than stay reachable. `state` is protected, so
         * a collection here keeps exactly what the circuit still needs.
         * Between gates is the only safe moment: the collector has explicit
         * roots and does not see a half-finished operation's temporaries.
         */
        /*
         * The thresholds, the every-16th-gate amortisation of the bitmap
         * scans, and the damping that stops a collection which freed nothing
         * from being repeated all live in limdd_gc_wanted. They were once
         * spelled out here instead, which is how the damping came to be
         * written and never used: this loop tested its own undamped copy.
         */
        if (limdd_gc_wanted()) CALL(limdd_gc);
    }

    stats.final_width = limdd_width(state, n);
    if (stats.final_width > stats.max_width) stats.max_width = stats.final_width;
    stats.final_qisq_size = limdd_max_wgt_bits(state, n);

    /* A final rebuild, so a deferred run is compared in its canonical form
     * rather than mid-batch. */
    if (!limdd_get_high_determinism()) state = limdd_canonize(state);

    stats.simulation_time = wctime() - t_start;
    stats.norm = limdd_norm_squared(state, n);

    /*
     * Computed from the final state, not inside the measurement branch: the
     * paper's benchmark circuits carry no measure statement at all, so that
     * branch never runs and the figure stayed at its initial zero -- which
     * looked exactly like a wrong state.
     *
     * The quantity is the QMDD path's: the probability of |0> on the TOP
     * qubit as a fraction of the norm. Reporting P(1), or P on some other
     * qubit, would make the two columns disagree where they do not.
     */
    {
        const double p1 = limdd_prob_qubit_one(state, 0, n);
        stats.first_qubit_prob = (stats.norm > 0.0)
                               ? (stats.norm - p1) / stats.norm : 1e10;
    }
    stats.final_nodes = limdd_countnodes(state);
    if (lim_stats) limdd_report_lim_stats(stderr, state, n);

    /* -v reads the amplitudes now, while the tables exist: fprint_stats
     * runs after limdd_nodes_quit. It used to read stats.final_state, which
     * only the QMDD arm sets, and loop forever on the null edge. */
    if (output_vector) {
        if (n > DECODE_MAX_QUBITS) {
            fprintf(stderr, "limdd: %u qubits is too many for -v\n", n);
        } else {
            const uint64_t len = UINT64_C(1) << n;
            final_vector = malloc(len * sizeof(EVBDD_WGT));
            bool *bits = malloc(n * sizeof(bool));
            for (uint64_t k = 0; k < len; k++) {
                bool *x = int_to_bitarray(k, n, !(circuit->reversed_qubit_order));
                for (uint32_t q = 0; q < n; q++) bits[q] = x[n - 1 - q];   /* as qmdd_get_amplitude */
                final_vector[k] = limdd_eval(state, bits, n);
                free(x);
            }
            free(bits);
        }
    }
    limdd_unprotect(&state);
    return 0;
}

/**
 * The BQD arm: the binary quotient diagram with scalar labels, on the same
 * Clifford+T gate set and through the same dispatch as the LIMDD arm.
 *
 * The paper proves two operations on it (tab:ops), the product on full
 * support (prop:prodscalar) and one monomial of a diagonal gate (prop:diag),
 * and those are what z, s, sdg, t, tdg, cz, cs, csdg and ccz take while the
 * state has full support: O(n) per gate. Every other gate, and a diagonal one
 * on a state with a zero, goes through the recursion of qsylvan_bqd_gates.h,
 * which is correct on every function and has no bound; the paper shows none
 * is possible for the Hadamard in general. An IQP circuit, a Hadamard on
 * every qubit and then diagonal gates, therefore stays on the proved path
 * after its first layer.
 *
 * Exact weights are what the arm is for. On floats the quotient rule is not
 * reliable once a gate can cancel: a residue where the value should be 0
 * turns a copy point into a ratio point and an amplitude of size 1 into 0,
 * and no tolerance avoids that (test_bqd_gates.c has the measurements). So
 * on floats only the gates that cannot cancel are taken: the diagonal ones,
 * and a Hadamard on a qubit no gate has touched, which is the IQP layer.
 */
static bool
is_phase_gate(const char *name)
{
    static const char *const phase[] = { "id", "z", "s", "sdg", "t", "tdg", "cz",
                                         "cs", "csdg", "ccz", NULL };
    for (int i = 0; phase[i] != NULL; i++) if (strcmp(name, phase[i]) == 0) return true;
    return false;
}

static bool
bqd_float_takes(const quantum_op_t *op, uint64_t *touched)
{
    uint64_t qs = UINT64_C(1) << op->targets[0];
    if (op->targets[1] >= 0) qs |= UINT64_C(1) << op->targets[1];
    for (int c = 0; c < 3; c++) if (op->ctrls[c] >= 0) qs |= UINT64_C(1) << op->ctrls[c];
    const bool ok = is_phase_gate(op->name)
                 || (strcmp(op->name, "h") == 0 && (*touched & qs) == 0);
    *touched |= qs;
    return ok;
}

TASK_1(int, bqd_simulate_circuit, quantum_circuit_t*, circuit)
{
    const uint32_t n = (uint32_t)circuit->qreg_size;
    if (n > 63) {
        fprintf(stderr, "bqd: %u qubits, but a vector index is one 64-bit word\n", n);
        return 1;
    }
    const bool exact = (sylvan_get_edge_weight_type() == WGT_QISQ2);
    uint64_t touched = 0;

    const double t_start = wctime();
    BQD state = bqd_basis_state(0, n);
    limdd_protect(&state);
    uint64_t gate_idx = 0;
    if (trace_out != NULL)
        fprintf(trace_out, "gate_index,gate,t_count,nodes,width,wgt_bits\n");

    for (quantum_op_t *op = circuit->operations; op != NULL; op = op->next) {
        if (op->type == op_measurement) break;   /* nothing to take from it */
        if (op->type != op_gate) continue;
        if (!exact && !bqd_float_takes(op, &touched)) {
            fprintf(stderr, "bqd: '%s' can cancel amplitudes, and on float weights "
                            "the quotient rule is not reliable then; use -e qisq2\n",
                    op->name);
            limdd_unprotect(&state);
            return 1;
        }
        if (!dd_apply_gate(&BQD_GATES, &state, op, n)) { limdd_unprotect(&state); return 1; }

        if (count_nodes) {
            const uint64_t c = limdd_countnodes(state);
            if (c > stats.max_nodes) stats.max_nodes = c;
        }
        if (trace_out != NULL) trace_row(gate_idx, op->name, state, n);
        gate_idx++;

        /* Between gates only, as on the LIMDD arm: the root is protected and
         * no operation's temporaries are live. */
        if (limdd_gc_wanted()) CALL(limdd_gc);
    }

    stats.final_width = limdd_width(state, n);
    if (stats.final_width > stats.max_width) stats.max_width = stats.final_width;
    stats.final_qisq_size = limdd_max_wgt_bits(state, n);
    stats.final_nodes = limdd_countnodes(state);
    stats.simulation_time = wctime() - t_start;

    /*
     * No norm and no measurement: the paper leaves every summation query
     * open, and computing one here by decoding would make the arm look as if
     * it had an algorithm it does not have. The two figures are taken from
     * the decoded vector when that fits, as a correctness check against the
     * other arms, and reported as -1 otherwise. The decode is not timed.
     */
    stats.norm = -1.0;
    stats.first_qubit_prob = -1.0;
    if (calc_measurement_prob || output_vector) {
        if (n > DECODE_MAX_QUBITS) {
            fprintf(stderr, "bqd: %u qubits is too many to decode; norm and "
                            "measurement probability are reported as -1\n", n);
        } else {
            const uint64_t len = UINT64_C(1) << n;
            EVBDD_WGT *vec = malloc(len * sizeof(EVBDD_WGT));
            bqd_to_vector(state, n, vec);
            double norm = 0.0, p_top_zero = 0.0;
            for (uint64_t x = 0; x < len; x++) {
                const complex_t c = weight_as_complex(vec[x]);
                const double m = c.r * c.r + c.i * c.i;
                norm += m;
                if ((x >> (n - 1)) == 0) p_top_zero += m;
            }
            stats.norm = norm;
            stats.first_qubit_prob = (norm > 0.0) ? p_top_zero / norm : 1e10;
            if (output_vector) {
                /* the vector has qubit q at bit n-1-q, and qubit q's value in
                 * print order k is x[n-1-q] (see qmdd_get_amplitude), so bit
                 * j of the index is x[j] */
                final_vector = malloc(len * sizeof(EVBDD_WGT));
                for (uint64_t k = 0; k < len; k++) {
                    bool *x = int_to_bitarray(k, n, !(circuit->reversed_qubit_order));
                    uint64_t idx = 0;
                    for (uint32_t j = 0; j < n; j++) if (x[j]) idx |= UINT64_C(1) << j;
                    final_vector[k] = vec[idx];
                    free(x);
                }
            }
            free(vec);
        }
    }
    limdd_unprotect(&state);
    return 0;
}

void simulate_circuit(quantum_circuit_t* circuit)
{
    double t_start = wctime();
    QMDD state = qmdd_create_all_zero_state(circuit->qreg_size);
    quantum_op_t *op = circuit->operations;
    while (op != NULL) {
        if (op->type == op_gate) {
            state = apply_gate(state, op, circuit->qreg_size);

        }
        else if (op->type == op_measurement) {
            if (circuit->has_intermediate_measurements) {
                state = measure(state, op, circuit);
            }
            else {
                double p;
                // don't set state = post measurement state
                qmdd_measure_all(state, circuit->qreg_size, circuit->creg, &p);
                if (circuit->reversed_qubit_order) {
                    reverse_bit_array(circuit->creg, circuit->qreg_size);
                }
                break;
            }
        }
        if (count_nodes) {
            uint64_t count = evbdd_countnodes(state);
            if (count > stats.max_nodes) stats.max_nodes = count;
        }
        if (count_qisq2_size) {
            uint64_t count = evbdd_qisqsize(state);
            if (count > stats.max_qisq_size) stats.max_qisq_size = count;
        }
        op = op->next;
    }
    if (count_qisq2_size) {
        uint64_t count = evbdd_qisqsize(state);
        stats.final_qisq_size = count;
    }
    stats.simulation_time = wctime() - t_start;
    stats.final_state = state;
    stats.shots = 1;
    stats.final_nodes = evbdd_countnodes(state);
    stats.norm = qmdd_get_norm(state, circuit->qreg_size);

    if (calc_measurement_prob){
        srand(12345);
        uint64_t k  = rand();
        bool *x = int_to_bitarray(k, circuit->qreg_size, !(circuit->reversed_qubit_order));
        complex_t c = qmdd_get_amplitude(stats.final_state, x, circuit->qreg_size);
        /*
         * c.i*c.i, not c.r*c.i. The second term was a typo, so this reported
         * sqrt(re^2 + re*im) rather than the modulus -- right only when the
         * amplitude is real, and not even monotone in |c| otherwise.
         *
         * Left alone: the name says probability but this is |c|, an
         * amplitude, and the probability is |c|^2. Changing that would move
         * published numbers, so it is flagged rather than silently altered.
         */
        stats.unnormed_prob = sqrt(c.r*c.r+c.i*c.i);
        if (stats.norm != 0.0){
            stats.normed_prob = stats.unnormed_prob/stats.norm;
        }
        else {
            stats.normed_prob = 1e10;
        }
        free(x);
    }
    if (calc_measurement_prob){
        /*
         * Through evbdd_get_topvar, not the root node's children: a QMDD
         * skips qubit 0 when both of its cofactors are the same, and then the
         * root node decides another qubit, or the root is the terminal and
         * its "node" is a reserved slot. The probability was of the wrong
         * qubit then (qreg q[2]; h q[0] gave 1.0), or an assertion on one
         * qubit. A skipped qubit 0 gets the same child twice, so 0.5.
         */
        EVBDD low_edge, high_edge;
        BDDVAR top;
        evbdd_get_topvar(state, 0, &top, &low_edge, &high_edge);
        double low_norm = qmdd_unnormed_prob(low_edge, 1, circuit->qreg_size);
        double high_norm = qmdd_unnormed_prob(high_edge, 1, circuit->qreg_size);
        if (low_norm != 0.0 || high_norm != 0.0){
            stats.first_qubit_prob = low_norm/(low_norm+high_norm);
        }
        else {
            stats.first_qubit_prob = 1e10;
        }
    }

}


/**
 * Sylvan operations must run inside a Lace worker (Lace >= 1.6); the thread
 * that calls lace_start() is not one. Both the simulation and writing the
 * stats (which computes amplitudes) touch the EVBDD/QMDD tables.
 */
VOID_TASK_1(run_simulation, quantum_circuit_t*, circuit)
{
    if (dd_kind != DD_QMDD) {
        /* Inside the task: llmsset claims buckets from a per-worker region,
         * and the thread that called lace_start is not a worker. */
        /* --canon wins; otherwise --dd picks: limdd canonicalises eagerly,
         * limdd-heur never does. */
        if (dd_kind == DD_BQD)
            ;                       /* no canonical form to schedule: every node is built canonical */
        else if (canon_policy >= 0)
            limdd_set_canon_policy((limdd_canon_policy_t) canon_policy, canon_interval);
        else
            limdd_set_canon_policy(dd_kind == DD_LIMDD ? LIMDD_CANON_ALWAYS
                                                       : LIMDD_CANON_MANUAL, 0);
        /*
         * Four tables plus a word per node bucket, so the node table's own
         * size is not a safe default here: at 2^25 that is over 2 GB before a
         * single node exists. Capped by default; --node-tab-size lifts the
         * cap, and is the only way to raise these tables, since unlike the
         * weight table they are allocated once and never grow.
         */
        size_t lt = node_tab_size_set ? max_tablesize : min_tablesize;
        if (!node_tab_size_set && lt > (1LL<<23)) lt = 1LL<<23;

        /*
         * The LIM table wants to be several times the node table, not equal
         * to it. Every live node holds one label, but the arithmetic mints
         * far more: each limdd_lim_mul as a label is pushed through a node,
         * each limdd_scale in the four-way gate mix, each inverse. Those are
         * dead as soon as the operation ends but occupy buckets until the
         * next collection, and measured between collections the table runs
         * about eight times the node count. At parity a 20-qubit, 700-gate
         * circuit filled the LIM table with its node table an eighth used.
         *
         * Four times, capped by what an edge can address -- beyond
         * LIMDD_LIM_MAX the index does not fit in an edge however much room
         * the table has.
         */
        size_t lim_t = lim_tab_size_log2 ? (1ULL << lim_tab_size_log2) : (lt << 2);
        if (lim_t > LIMDD_LIM_MAX) lim_t = LIMDD_LIM_MAX;

        /* The circuit's width, not LIMDD_MAX_QUBITS: the recursions stop when
         * they reach it, so a value larger than the diagram is deep sends
         * them past the terminal. */
        if (dd_kind == DD_BQD) {
            /* The same tables, read under the quotient rule; see qsylvan_bqd.h. */
            bqd_init(BQD_FAMILY_SCALAR, circuit->qreg_size, lt, lt, lim_t, lt);
            if (CALL(bqd_simulate_circuit, circuit) != 0) {
                bqd_quit();
                exit(1);
            }
            bqd_quit();
        } else {
            limdd_nodes_init(circuit->qreg_size, lt, lt, lim_t, lt);
            if (CALL(limdd_simulate_circuit, circuit) != 0) {
                limdd_nodes_quit();
                exit(1);
            }
            limdd_nodes_quit();
        }
    } else {
        simulate_circuit(circuit);
    }

    if (json_outputfile != NULL) {
        FILE *fp = fopen(json_outputfile, "w");
        fprint_stats(fp, circuit);
        fclose(fp);
    } else {
        fprint_stats(stdout, circuit);
    }
}


int main(int argc, char *argv[])
{
    argp_parse(&argp, argc, argv, 0, 0, 0);
    quantum_circuit_t* circuit = parse_qasm_file(qasm_inputfile);
    if (reorder_qubits)
        optimize_qubit_order(circuit, reorder_qubits == 2);

    /*
     * LIMDD needs the relative merging rule, and gets it unless asked
     * otherwise.
     *
     * A LIMDD stores no weight on a low edge: the whole scale of a subdiagram
     * rides on the parent, so weights range over many orders of magnitude
     * instead of sitting near 1 the way a QMDD's do. A single ABSOLUTE
     * threshold cannot serve that range -- every weight smaller than it merges
     * with every other, which silently rewrites the state. On a 10-qubit,
     * 300-gate Clifford+T circuit the absolute rule loses norm (0.94 rather
     * than 1) and the answer moves with the weight table's size, since which
     * entry a weight meets first depends on where it hashes. Relative merging
     * bounds the error by the flow through a weight rather than by its
     * magnitude, which is exactly the insensitivity a LIMDD needs; see
     * cmap.c and Brand et al.
     *
     * QMDD keeps the historical absolute default: its normalisation holds
     * weights near 1, and the paper's numbers were taken that way.
     */
    /*
     * Allocate the whole weight table up front.
     *
     * limdd_gc does collect it now, and wgt_table_gc_init_new doubles it
     * when it runs, so in principle --wgt-tab-size could be a maximum to
     * grow towards. In practice growth never gets the chance: a collection
     * can only happen between gates, and one gate on a large diagram mints
     * millions of weights, so the table fills mid-gate and dies before
     * anything can enlarge it. Measured -- starting at 2^23 and growing,
     * rand_n20_d700 fails at every maximum up to 2^28; starting at 2^24 it
     * finishes in 53s.
     */
    if (dd_kind != DD_QMDD) {
        min_wgt_tab_size = max_wgt_tab_size;
    }

    /*
     * A bigger operation cache for LIMDD, which needs one far more than a
     * QMDD does. Every node it builds runs a coset minimisation whose result
     * is memoised, and those entries are large in number and expensive to
     * recompute, where a QMDD's are neither. Measured on a 32-qubit graph
     * state: LIMDD goes from 32.0s at the 2^16 default to 11.2s at 2^22,
     * while the same circuit on a QMDD moves from 2.33s to 1.80s and two
     * other circuits get slightly slower. 2^20 is 33 MB and captures nearly
     * all of it -- 11.3s against 10.7s at four times the size.
     */
    if (dd_kind != DD_QMDD && !cache_size_set) {
        min_cachesize = max_cachesize = 1LL << 20;
    }

    if (dd_kind != DD_QMDD && rel_tolerance < 0 && !force_absolute) {
        /*
         * 1e-12, not the storage tolerance.
         *
         * A LIMDD's whole claim is that a stabiliser state costs O(n) nodes,
         * and that holds only if subdiagrams denoting the same state actually
         * merge. Their scalars are reached by different multiplication orders,
         * so they differ in the last bits, and at 1e-14 those differences
         * survive: a 20-qubit, 400-gate Clifford circuit -- no T gates at all
         * -- came out as 7816 nodes with 92% of the stabiliser groups trivial,
         * where exact arithmetic gives 19 nodes and none trivial. One decade
         * coarser absorbs the accumulated rounding and reproduces the exact
         * diagram node for node on 18 of 19 circuits measured, norms and
         * probabilities matching throughout.
         *
         * Erring coarse is the safe direction here, and measurably so: across
         * that suite 1e-12 never produced FEWER nodes than exact arithmetic,
         * which is what over-merging would look like. It is also why this does
         * not track -t, whose default of 1e-14 is a storage tolerance for
         * weights near 1 rather than a statement about accumulated error.
         */
        rel_tolerance = 1e-12;
        /*
         * And no zero-collapse, unless asked for. The hybrid rule normally
         * rounds anything tiny down to exactly zero, which suits a QMDD,
         * where a weight that small is numerical dust. In a LIMDD it can be
         * the entire scale of the state: the root weight of a 300-gate
         * circuit is a product of a few hundred factors and is legitimately
         * far below 1e-14. Collapsing it deletes the state rather than
         * cleaning it up, and costs about 1% of the norm on q_10 at one
         * weight table size out of six.
         */
        if (!zero_tolerance_set) zero_tolerance = 0.0;
    }

    /*
     * The normalisation strategy is a QMDD notion: which of a node's two
     * weights is divided out. A LIMDD has no such choice -- the low edge
     * carries the identity and all scale goes up -- so the strategy is unused,
     * and demanding one the backend supports would refuse exact weights for no
     * reason. qisq2 has no absolute value and so rejects max, min and L2.
     */
    if (dd_kind != DD_QMDD && !qsylvan_norm_supported(wgt_table_type, wgt_norm_strat)) {
        wgt_norm_strat = NORM_LOW;
    }

    if (trace_file != NULL) {
        trace_out = fopen(trace_file, "w");
        if (trace_out == NULL) {
            fprintf(stderr, "cannot write the trace to %s\n", trace_file);
            exit(1);
        }
    }

    if (dd_kind == DD_BQD) canon_name = "n/a";

    if (rseed == 0) rseed = time(NULL);
    srand(rseed);
    
    // Standard Lace initialization
    lace_start(workers, 0);

    // Simple Sylvan initialization
    sylvan_set_sizes(min_tablesize, max_tablesize, min_cachesize, max_cachesize);
    sylvan_init_package();
    /*
     * Before qsylvan_init_simulator, not after: the rule decides how a weight
     * is hashed, and the table is populated during initialisation (1, 0 and
     * -1, then the gate entries). Switching afterwards leaves those at
     * positions the new hash never probes, and the next lookup of 1 creates a
     * second entry for it -- two indices for one value, which breaks every
     * comparison that goes through EVBDD_ONE.
     */
    if (force_absolute) rel_tolerance = -1;
    if (rel_tolerance >= 0)
        sylvan_edge_weights_set_hybrid_tolerance(rel_tolerance, zero_tolerance);

    qsylvan_init_simulator(min_wgt_tab_size, max_wgt_tab_size, tolerance, wgt_table_type, wgt_norm_strat);
    wgt_set_inverse_chaching(wgt_inv_caching);

    RUN(run_simulation, circuit);

    if (trace_out != NULL) fclose(trace_out);

    sylvan_quit();
    lace_stop();
    free_quantum_circuit(circuit);

    return 0;
}
