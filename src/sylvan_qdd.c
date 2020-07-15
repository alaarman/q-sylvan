/*
 * Copyright 2011-2016 Formal Methods and Tools, University of Twente
 * Copyright 2016-2017 Tom van Dijk, Johannes Kepler University Linz
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

#include <sylvan_int.h>
#include <sylvan_qdd.h>

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "sylvan_qdd_int.h"

//static int granularity = 1; // default

// TODO: make a cleaner separation between the code in
// - sylvan_qdd.h
// - sylvan_qdd.c
// - sylvan_qdd_int.h
// - sylvan_qdd_int.c
// consistent with for example the equivalent mtbdd files.


/****************< (bit level) manipulation of QDD / qddnode_t >***************/
/**
 * QDD node structure (128 bits)
 * 
 * 64 bits low:
 *       4 bits: lower 4 bits of 7 bit variable/qubit number of this node
 *      20 bits: index of edge weight of low edge in Ctable (AMP)
 *      40 bits: low edge pointer to next node (PTR)
 * 64 bits high:
 *       1 bit:  marked/unmarked flag (TODO: use first bit of PTR for this? (reserved for sylvan_invalid))
 *       3 bits: upper 3 bits of 7 bit variable/qubit number of this node
 *      20 bits: index of edge weight of high edge in Ctable (AMP)
 *      40 bits: high edge pointer to next node (PTR)
 */
typedef struct __attribute__((packed)) qddnode {
    QDD low, high;
} *qddnode_t; // 16 bytes


/**
 * Gets only the AMP information of a QDD edge `q`.
 */
static inline AMP
QDD_AMP(QDD q)
{
    return (q & 0x0fffff0000000000) >> 40; // 20 bits
}

/**
 * Gets only the PTR information of a QDD edge `q`.
 */
static inline PTR
QDD_PTR(QDD q)
{
    return q & 0x000000ffffffffff; // 40 bits
}

/**
 * 2 BDDVARs (assumed max 8 bits each)
 */
static inline uint32_t
QDD_PARAM_PACK_16(BDDVAR a, BDDVAR b) 
{
    return b<<8 | a;
}

/**
 *  24 bit gateid, 2 possible qubit parameters (e.g. control/target)
 */
static inline uint64_t
GATE_OPID_40(uint32_t gateid, BDDVAR a, BDDVAR b)
{
    uint64_t res = ((uint64_t)b)<<32 | ((uint64_t)a)<<24 | gateid;
    return res;
}

/**
 * 24 bits gateid, 5 possible qubit parameters (e.g. target/control/range)
 */
static inline uint64_t
GATE_OPID_64(uint32_t gateid, BDDVAR a, BDDVAR b, BDDVAR c, BDDVAR d, BDDVAR e)
{
    uint64_t res = ((uint64_t)e)<<56 | 
                   ((uint64_t)d)<<48 | 
                   ((uint64_t)c)<<40 | 
                   ((uint64_t)b)<<32 | 
                   ((uint64_t)a)<<24 | 
                   gateid;
    return res;
}

/**
 * Gets the variable number of a given node `n`.
 */
static inline BDDVAR
qddnode_getvar(qddnode_t n)
{
    return (BDDVAR) ( (n->low >> 60) | ((n->high >> 56) & 0x70) ); // 7 bits
}

/**
 * Gets the low edge of `n` with the AMP and PTR information, but without the
 * (halved) variable information.
 */
static inline QDD
qddnode_getlow(qddnode_t n)
{
    return (QDD) n->low & 0x0fffffffffffffff; // 60 bits
}

/**
 * Gets the high edge of `n` with the AMP and PTR information, but without the
 * (halved) variable information.
 */
static inline QDD
qddnode_gethigh(qddnode_t n)
{
    return (QDD) n->high & 0x0fffffffffffffff; // 60 bits
}

/**
 * Gets only the PTR of the low edge of `n`.
 */
static inline PTR
qddnode_getptrlow(qddnode_t n)
{
    return (PTR) QDD_PTR(n->low);
}

/**
 * Gets only the PTR of the high edge of `n`.
 */
static inline PTR
qddnode_getptrhigh(qddnode_t n)
{
    return (PTR) QDD_PTR(n->high);
}

/**
 * Gets only the AMP of the low edge of `n`.
 */
static inline AMP
qddnode_getamplow(qddnode_t n)
{
    return (AMP) QDD_AMP(n->low);
}

/**
 * Gets only the AMP of the high edge of `n`.
 */
static inline AMP
qddnode_getamphigh(qddnode_t n)
{
    return (AMP) QDD_AMP(n->high);
}

/**
 * Gets the value of the "marked" flag.
 */
static inline bool
qddnode_getmark(qddnode_t n)
{
    return n->high & 0x8000000000000000 ? 1 : 0;
}

/**
 * Sets the value of the "marked" flag to `mark`.
 */
static inline void
qddnode_setmark(qddnode_t n, bool mark)
{
    if (mark) n->high |= 0x8000000000000000; // set 1th bit from left to 1
    else      n->high &= 0x7fffffffffffffff; // set 1th bit from left to 0
}

/**
 * Gets the node `p` is pointing to.
 * TODO (?) return special node for when p == QDD_TERMINAL
 */
static inline qddnode_t
QDD_GETNODE(PTR p)
{
    return (qddnode_t) llmsset_index_to_ptr(nodes, p);
}

/**
 * Packs a PTR and AMP into a single 64 bit QDD.
 */
static inline QDD
qdd_bundle_ptr_amp(PTR p, AMP a)
{
    assert (p <= 0x000000fffffffffe);   // avoid clash with sylvan_invalid
    assert (a <= 0x00000000000fffff);
    return (a << 40 | p);
}


static inline QDD
qdd_bundle_low(BDDVAR var, PTR p, AMP a)
{
    // on the low edge we store the bottom 4 bits of the 8 bit var
    assert (var <= 0x7f);
    QDD q = qdd_bundle_ptr_amp(p, a);
    q = ((uint64_t)var << 60) | q;
    return q;
}


static inline QDD
qdd_bundle_high(BDDVAR var, PTR p, AMP a)
{
    // on the high edge we store the top 4 bits of the 8 bit var
    assert (var <= 0x7f);
    QDD q = qdd_bundle_ptr_amp(p, a);
    return (((uint64_t)var & 0x70) << 56) | q;
}


static inline void 
qddnode_make(qddnode_t n, BDDVAR var, PTR low, PTR high, AMP a, AMP b)
{
    n->low  = qdd_bundle_low(var, low, a);
    n->high = qdd_bundle_high(var, high, b);
}


/***************</ (bit level) manipulation of QDD / qddnode_t >***************/


static PTR
_qdd_makenode(BDDVAR var, PTR low, PTR high, AMP a, AMP b)
{
    struct qddnode n;

    qddnode_make(&n, var, low, high, a, b);

    PTR result;
    int created;
    PTR index = llmsset_lookup(nodes, n.low, n.high, &created);
    if (index == 0) {
        LACE_ME;

        mtbdd_refs_push(n.low);//mtbdd_refs_push(low);
        mtbdd_refs_push(n.high);//mtbdd_refs_push(high);
        sylvan_gc();
        mtbdd_refs_pop(2);

        index = llmsset_lookup(nodes, n.low, n.high, &created);
        if (index == 0) {
            fprintf(stderr, "BDD Unique table full, %zu of %zu buckets filled!\n", llmsset_count_marked(nodes), llmsset_get_size(nodes));
            exit(1);
        }
    }

    if (created) sylvan_stats_count(BDD_NODES_CREATED);
    else sylvan_stats_count(BDD_NODES_REUSED);

    result = index;
    return result;
} 

static QDD // (PTR and AMP, but the amp is the norm weight from below)
qdd_makenode(BDDVAR var, QDD low_edge, QDD high_edge)
{ 
    PTR low_ptr  = QDD_PTR(low_edge);
    AMP low_amp  = QDD_AMP(low_edge);
    PTR high_ptr = QDD_PTR(high_edge);
    AMP high_amp = QDD_AMP(high_edge);

    PTR p;
    AMP norm;

    // Edges with weight 0 should point straight to terminal.
    if(low_amp  == C_ZERO) low_ptr  = QDD_TERMINAL;
    if(high_amp == C_ZERO) high_ptr = QDD_TERMINAL;

    // If both low and high are the same (both PTR and AMP) return low
    if(low_ptr == high_ptr && low_amp == high_amp) {
        return low_edge;
    }
    else{
        // If the edges are not the same
        if(low_amp != C_ZERO){
            // Normalize using low
            norm     = low_amp;
            low_amp  = C_ONE;
            high_amp = Cdiv(high_amp, norm);
        }
        else{
            // Unless low amp is 0, then norm using high
            norm     = high_amp;
            high_amp = C_ONE;
        }
        p = _qdd_makenode(var, low_ptr, high_ptr, low_amp, high_amp);
        return qdd_bundle_ptr_amp(p, norm);
    }
}


/**
 * Helper function for recursive unmarking
 */
static void
qdd_unmark_rec(QDD qdd)
{
    if (QDD_PTR(qdd) == QDD_TERMINAL) return;
    qddnode_t n = QDD_GETNODE(QDD_PTR(qdd));
    if (!qddnode_getmark(n)) return;
    qddnode_setmark(n, 0);
    qdd_unmark_rec(qddnode_getlow(n));
    qdd_unmark_rec(qddnode_gethigh(n));
}


/**
 * Counts nodes in the qdd by marking them.
 */
static uint64_t
qdd_nodecount_mark(QDD qdd)
{
    if (QDD_PTR(qdd) == QDD_TERMINAL) return 0; // don't (repeat) count terminal
    qddnode_t n = QDD_GETNODE(QDD_PTR(qdd));
    if (qddnode_getmark(n)) return 0;
    qddnode_setmark(n, 1);
    return 1 + qdd_nodecount_mark(qddnode_getlow(n)) + qdd_nodecount_mark(qddnode_gethigh(n));
}

uint64_t
qdd_countnodes(QDD qdd)
{
    uint64_t res = qdd_nodecount_mark(qdd) + 1; // (+ 1 for terminal "node")
    qdd_unmark_rec(qdd);
    return res;
}


/**************************<cleaning amplitude table>**************************/

void
clean_amplitude_table(QDD qdds[], int n_qdds)
{
    LACE_ME;
    // 1. Create new amp table
    init_new_empty_table();

    // 2. Fill new table with amps present in given QDDs
    for (int i = 0; i < n_qdds; i++) qdds[i] = _fill_new_amp_table(qdds[i]);

    // 3. Delete old amp table
    delete_old_table();

    // 4. Any cache we migh have is now invalid because the same amplitudes 
    //    might now have different indices in the amp table
    cache_clear();
}

TASK_IMPL_1(QDD, _fill_new_amp_table, QDD, qdd)
{
    // Check cache
    QDD res;
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_QDD_CLEAN_AMP_TABLE, 0LL, qdd, 0LL, &res)) {
            return res;
        }
    }

    // Move amp from old to new table, get new index
    AMP new_amp = move_from_old_to_new(QDD_AMP(qdd));
    qdd = qdd_bundle_ptr_amp(QDD_PTR(qdd), new_amp);

    // If terminal, return
    if (QDD_PTR(qdd) == QDD_TERMINAL) return qdd;
    
    // Recursive for children
    QDD low, high;
    qddnode_t n = QDD_GETNODE(QDD_PTR(qdd));
    bdd_refs_spawn(SPAWN(_fill_new_amp_table, qddnode_gethigh(n)));
    low = CALL(_fill_new_amp_table, qddnode_getlow(n));
    bdd_refs_push(low);
    high = bdd_refs_sync(SYNC(_fill_new_amp_table));
    bdd_refs_pop(1);

    // We don't need to use the 'qdd_makenode()' function which normalizes the 
    // amplitudes, because the QDD doesn't actually change, only the AMP indices,
    // but none of the actual values.
    PTR ptr = _qdd_makenode(qddnode_getvar(n), QDD_PTR(low), QDD_PTR(high), 
                                               QDD_AMP(low), QDD_AMP(high));

    // Put in cache, return
    res = qdd_bundle_ptr_amp(ptr, new_amp);
    if (cachenow) cache_put3(CACHE_QDD_CLEAN_AMP_TABLE, 0LL, qdd, 0LL, res);
    return res;
}

/*************************</cleaning amplitude table>**************************/



/*******************************<applying gates>*******************************/

TASK_IMPL_2(QDD, qdd_plus, QDD, a, QDD, b)
{
    // Trivial cases
    if(QDD_AMP(a) == C_ZERO) return b;
    if(QDD_AMP(b) == C_ZERO) return a;

    // Check cache
    QDD res;
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_QDD_PLUS, sylvan_false, a, b, &res)) {
            sylvan_stats_count(QDD_PLUS_CACHED);
            return res;
        }
    }

    // Get var(a) and var(b)
    QDD low_a, low_b, high_a, high_b;
    BDDVAR var_a = UINT32_MAX, var_b = UINT32_MAX, topvar;
    if (QDD_PTR(a) != QDD_TERMINAL) {
        qddnode_t node = QDD_GETNODE(QDD_PTR(a));
        var_a  = qddnode_getvar(node);
    }
    if (QDD_PTR(b) != QDD_TERMINAL) {
        qddnode_t node = QDD_GETNODE(QDD_PTR(b));
        var_b  = qddnode_getvar(node);
    }

    // Check which is topvar, insert node for QDD with skipped case
    low_a  = qdd_bundle_ptr_amp(QDD_PTR(a), C_ONE);
    high_a = qdd_bundle_ptr_amp(QDD_PTR(a), C_ONE);
    low_b  = qdd_bundle_ptr_amp(QDD_PTR(b), C_ONE);
    high_b = qdd_bundle_ptr_amp(QDD_PTR(b), C_ONE);
    if (var_a <= var_b) { // didn't skip in QDD a
        qddnode_t node = QDD_GETNODE(QDD_PTR(a));
        low_a  = qddnode_getlow(node);
        high_a = qddnode_gethigh(node);
        topvar = var_a;
    }
    if (var_b <= var_a) { // didn't skip in QDD b
        qddnode_t node = QDD_GETNODE(QDD_PTR(b));
        low_b  = qddnode_getlow(node);
        high_b = qddnode_gethigh(node);
        topvar = var_b;
    }

    // Base/terminal case: same target and same variable
    if(QDD_PTR(a) == QDD_PTR(b) && var_a == var_b){
        AMP sum = Cadd(QDD_AMP(a), QDD_AMP(b));
        res = qdd_bundle_ptr_amp(QDD_PTR(a), sum);
        return res;
    }

    // If not base/terminal case, pass edge weight of current edge down
    AMP amp_la, amp_ha, amp_lb, amp_hb;
    amp_la = Cmul(QDD_AMP(a), QDD_AMP(low_a));
    amp_ha = Cmul(QDD_AMP(a), QDD_AMP(high_a));
    amp_lb = Cmul(QDD_AMP(b), QDD_AMP(low_b));
    amp_hb = Cmul(QDD_AMP(b), QDD_AMP(high_b));
    low_a  = qdd_bundle_ptr_amp(QDD_PTR(low_a),  amp_la);
    high_a = qdd_bundle_ptr_amp(QDD_PTR(high_a), amp_ha);
    low_b  = qdd_bundle_ptr_amp(QDD_PTR(low_b),  amp_lb);
    high_b = qdd_bundle_ptr_amp(QDD_PTR(high_b), amp_hb);

    // Recursive calls down
    QDD low, high;
    bdd_refs_spawn(SPAWN(qdd_plus, high_a, high_b));
    low = CALL(qdd_plus, low_a, low_b);
    bdd_refs_push(low);
    high = bdd_refs_sync(SYNC(qdd_plus));
    bdd_refs_pop(1);

    // Put in cache, return
    res = qdd_makenode(topvar, low, high);
    if (cachenow) {
        if (cache_put3(CACHE_QDD_PLUS, sylvan_false, a, b, res)) 
            sylvan_stats_count(QDD_PLUS_CACHEDPUT);
    }
    return res;
}

TASK_IMPL_3(QDD, qdd_gate, QDD, q, uint32_t, gate, BDDVAR, target)
{
    // Check cache
    QDD res;
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_QDD_GATE, GATE_OPID_40(gate, target, 0), q, sylvan_false, &res)) {
            sylvan_stats_count(QDD_GATE_CACHED);
            return res;
        }
    }

    // TODO: make has_skipped(BDDVAR target) function
    bool skipped = false;
    bool at_targ = false;
    BDDVAR var;
    if(QDD_PTR(q) == QDD_TERMINAL) {
        skipped = true;
    }
    else {
        qddnode_t node = QDD_GETNODE(QDD_PTR(q));
        var = qddnode_getvar(node);
        if (var >  target) skipped = true;
        if (var == target) at_targ = true;
    }

    QDD low, high;
    if (skipped) {
        low  = qdd_bundle_ptr_amp(QDD_PTR(q), C_ONE);
        high = qdd_bundle_ptr_amp(QDD_PTR(q), C_ONE);
        var  = target;
    }
    else {
        qddnode_t node = QDD_GETNODE(QDD_PTR(q));
        low  = qddnode_getlow(node);
        high = qddnode_gethigh(node);
    } // has_skipped(BDDVAR target) function could do all this work

    if (skipped || at_targ) { // apply gate here
        AMP a_u00 = Cmul(QDD_AMP(low), gates[gate][0]);
        AMP a_u10 = Cmul(QDD_AMP(low), gates[gate][2]);
        AMP b_u01 = Cmul(QDD_AMP(high), gates[gate][1]);
        AMP b_u11 = Cmul(QDD_AMP(high), gates[gate][3]);
        QDD qdd1 = qdd_makenode(target, qdd_bundle_ptr_amp(QDD_PTR(low), a_u00), 
                                        qdd_bundle_ptr_amp(QDD_PTR(low), a_u10));
        QDD qdd2 = qdd_makenode(target, qdd_bundle_ptr_amp(QDD_PTR(high),b_u01),
                                        qdd_bundle_ptr_amp(QDD_PTR(high),b_u11));
        res = qdd_plus(qdd1, qdd2);
    }
    else { // not at target qubit yet, recursive calls down
        bdd_refs_spawn(SPAWN(qdd_gate, high, gate, target));
        low = CALL(qdd_gate, low, gate, target);
        bdd_refs_push(low);
        high = bdd_refs_sync(SYNC(qdd_gate));
        bdd_refs_pop(1);
        res  = qdd_makenode(var, low, high);
    }

    // Multiply root amp of sum with input root amp, add to cache, return
    AMP new_root_amp = Cmul(QDD_AMP(q), QDD_AMP(res));
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    if (cachenow) {
        if (cache_put3(CACHE_QDD_GATE, GATE_OPID_40(gate, target, 0), q, sylvan_false, res)) 
            sylvan_stats_count(QDD_GATE_CACHEDPUT);
    }
    return res;
}

TASK_IMPL_4(QDD, qdd_cgate, QDD, q, uint32_t, gate, BDDVAR, c, BDDVAR, t)
{
    assert (c < t);

    // Check cache
    QDD res;
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_QDD_CGATE, GATE_OPID_40(gate, c, t), q, sylvan_false, &res)) {
            sylvan_stats_count(QDD_CGATE_CACHED);
            return res;
        }
    }

    // TODO: make has_skipped(BDDVAR target) function
    bool skipped = false;
    bool at_targ = false;
    BDDVAR var;
    if(QDD_PTR(q) == QDD_TERMINAL) {
        skipped = true;
    }
    else {
        qddnode_t node = QDD_GETNODE(QDD_PTR(q));
        var = qddnode_getvar(node);
        if (var >  c) skipped = true;
        if (var == c) at_targ = true;
    }

    QDD low, high;
    if (skipped) {
        low  = qdd_bundle_ptr_amp(QDD_PTR(q), C_ONE);
        high = qdd_bundle_ptr_amp(QDD_PTR(q), C_ONE);
        var  = c;
    }
    else {
        qddnode_t node = QDD_GETNODE(QDD_PTR(q));
        low  = qddnode_getlow(node);
        high = qddnode_gethigh(node);
    }// has_skipped(BDDVAR target) function could do all this work

    if (skipped || at_targ) { // apply gate to high, leave low unchanged
        high = qdd_gate(high, gate, t);
    }
    else { // not at target qubit yet, recursive calls down
        bdd_refs_spawn(SPAWN(qdd_cgate, high, gate, c, t));
        low = CALL(qdd_cgate, low, gate, c, t);
        bdd_refs_push(low);
        high = bdd_refs_sync(SYNC(qdd_cgate));
        bdd_refs_pop(1);
    }
    res  = qdd_makenode(var, low, high);

    // Multiply root amp of sum with input root amp, add to cache, return
    AMP new_root_amp = Cmul(QDD_AMP(q), QDD_AMP(res));
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    if (cachenow) {
        if (cache_put3(CACHE_QDD_CGATE, GATE_OPID_40(gate, c, t), q, sylvan_false, res)) 
            sylvan_stats_count(QDD_CGATE_CACHEDPUT);
    }
    return res;
}

/******************************</applying gates>*******************************/



/*********************<applying (controlled) sub-circuits>*********************/

QDD
qdd_circuit_swap(QDD qdd, BDDVAR qubit1, BDDVAR qubit2)
{
    assert (qubit1 < qubit2);
    
    LACE_ME; // TODO: make this function a LACE function instead?
    QDD res;

    // CNOT
    res = qdd_cgate(qdd, GATEID_X, qubit1, qubit2);
    // upside down CNOT (equivalent)
    res = qdd_gate(res, GATEID_H, qubit1);
    res = qdd_cgate(res, GATEID_Z, qubit1, qubit2);
    res = qdd_gate(res, GATEID_H, qubit1);
    // CNOT
    res = qdd_cgate(res, GATEID_X, qubit1, qubit2);

    return res;
}

QDD
qdd_circuit_swap_range(QDD qdd, BDDVAR first, BDDVAR last)
{
    QDD res = qdd;
    BDDVAR a, b;
    int num_qubits = (last - first) + 1;
    for (int j = 0; j < (int)(num_qubits/2); j++) {
        a = first + j;
        b = last  - j;
        res = qdd_circuit_swap(res, a, b);
    }
    return res;
}

QDD
qdd_circuit_QFT(QDD qdd, BDDVAR first, BDDVAR last)
{
    LACE_ME;

    int k;
    QDD res = qdd;
    BDDVAR a, b;
    for (a = first; a <= last; a++) {
        
        // H gate on current qubit
        res = qdd_gate(res, GATEID_H, a);

        // Controlled phase gates on all qubits below
        for (b = a+1; b <= last; b++) {
            k = (b - a) + 1;
            res = qdd_cgate(res, GATEID_Rk(k), a, b);
        }
    }

    // Note that we're not swapping the qubit order in this function

    return res;
}

QDD
qdd_circuit_QFT_inv(QDD qdd, BDDVAR first, BDDVAR last)
{
    LACE_ME;
    
    int k;
    QDD res = qdd;
    BDDVAR a, b;

    // Note that we're not swapping the qubit order in this function
    
    // H gates and phase gates (but now backwards)
    for (a = last + 1; a-- > first; ) { // weird for loop because BDDVARs are unsigned

        // Controlled phase gates (negative angles this time)
        for (b = last; b >= (a+1); b--){
            k = (b - a) + 1;
            res = qdd_cgate(res, GATEID_Rk_dag(k), a, b);
        }

        // H on current qubit
        res = qdd_gate(res, GATEID_H, a);
    }

    return res;
}

QDD
qdd_circuit(QDD qdd, uint32_t circ_id, BDDVAR t1, BDDVAR t2)
{
    switch (circ_id) {  // don't judge me please
        case CIRCID_swap          : return qdd_circuit_swap(qdd, t1, t2);
        case CIRCID_swap_range    : return qdd_circuit_swap_range(qdd, t1, t2);
        case CIRCID_QFT           : return qdd_circuit_QFT(qdd, t1, t2);
        case CIRCID_QFT_inv       : return qdd_circuit_QFT_inv(qdd, t1, t2);
        case CIRCID_phi_add_a     : return qdd_phi_add(qdd, t1, t2, shor_bits_a);
        case CIRCID_phi_add_N     : return qdd_phi_add(qdd, t1, t2, shor_bits_N);
        case CIRCID_phi_add_a_inv : return qdd_phi_add_inv(qdd, t1, t2, shor_bits_a);
        case CIRCID_phi_add_N_inv : return qdd_phi_add_inv(qdd, t1, t2, shor_bits_N);
        default :
            assert ("Invalid circuit ID" && false);
    }
}

TASK_IMPL_6(QDD, qdd_ccircuit, QDD, qdd, uint32_t, circ_id, BDDVAR*, cs, uint32_t, ci, BDDVAR, t1, BDDVAR, t2)
{
    // Cache lookup
    QDD res;
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_QDD_SUBCIRC, sylvan_false, qdd, 
                       GATE_OPID_64(circ_id, cs[0], cs[1], cs[2], t1, t2), 
                       &res)) {
            return res;
        }
    }

    // Get current control qubit
    BDDVAR c = cs[ci];

    // If no more control qubits, apply sub-circ here
    if (c == UINT8_MAX || ci > MAX_CONTROLS) {
        res = qdd_circuit(qdd, circ_id, t1, t2);
    }
    else {
        // Check if skipped control node
        bool skipped = false;
        if (QDD_PTR(qdd) == QDD_TERMINAL) skipped = true;
        else {
            qddnode_t node = QDD_GETNODE(QDD_PTR(qdd));
            if (qddnode_getvar(node) > c) skipped = true;
        }

        // Check if we need to control here
        BDDVAR var;
        QDD low, high;
        bool control_here = false;
        if (skipped) { // var > control
            low  = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
            high = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
            var  = c;
            control_here = true;
        }
        else {
            // not skipped, either (var == control) or (var < control)
            qddnode_t node = QDD_GETNODE(QDD_PTR(qdd));
            var  = qddnode_getvar(node);
            low  = qddnode_getlow(node);
            high = qddnode_gethigh(node);
            if (var == c) control_here = true;
        }

        if (control_here) {
            ci++; // next control qubit
            high = CALL(qdd_ccircuit, high, circ_id, cs, ci, t1, t2);
        }
        else {
            // recursive call to both children
            bdd_refs_spawn(SPAWN(qdd_ccircuit, high, circ_id, cs, ci, t1, t2));
            low = CALL(qdd_ccircuit, low, circ_id, cs, ci, t1, t2);
            bdd_refs_push(low);
            high = bdd_refs_sync(SYNC(qdd_ccircuit));
            bdd_refs_pop(1);
        }
        res = qdd_makenode(var, low, high); 
    }

    // Multiply root amp of sum with input root amp, add to cache, return
    AMP new_root_amp = Cmul(QDD_AMP(qdd), QDD_AMP(res));
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    if (cachenow) {
        cache_put3(CACHE_QDD_SUBCIRC, sylvan_false, qdd, 
                   GATE_OPID_64(CIRCID_swap, cs[0], cs[1], cs[2], t1, t2), 
                   res);
    }
    return res;
}

TASK_IMPL_4(QDD, qdd_all_control_phase, QDD, qdd, BDDVAR, k, BDDVAR, n, bool*, x)
{
    // TODO: remove LACE, no branching in this function
    assert(k < n);
    
    bool skipped_k = false;
    qddnode_t node;
    if (QDD_PTR(qdd) == QDD_TERMINAL) {
        skipped_k = true;
    }
    else {
        node = QDD_GETNODE(QDD_PTR(qdd));
        BDDVAR var = qddnode_getvar(node);
        if(var > k) {
            skipped_k = true;
        }
    }

    QDD low, high;
    if (skipped_k) {
        // insert skipped node
        low  = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
        high = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
    }
    else {
        // case var == k (var < k shouldn't happen)
        low  = qddnode_getlow(node);
        high = qddnode_gethigh(node);
    }

    // terminal case, apply phase depending on x[k] (control k on 0 or 1)
    if (k == (n-1)) {
        if (x[k] == 1) {
            AMP new_amp = Cmul(QDD_AMP(high), Clookup(Cmake(-1.0, 0.0)));
            high = qdd_bundle_ptr_amp(QDD_PTR(high), new_amp);
        }
        else {
            AMP new_amp = Cmul(QDD_AMP(low), Clookup(Cmake(-1.0, 0.0)));
            low = qdd_bundle_ptr_amp(QDD_PTR(low), new_amp);
        }
    }
    // non terminal case, choose low/high depending on x[k] (control k on 0 or 1)
    else {
        if (x[k] == 1) {
            k++; // next level
            high = CALL(qdd_all_control_phase, high, k, n, x);
            k--;
        }
        else {
            k++;
            low = CALL(qdd_all_control_phase, low, k, n, x);
            k--;
        }
    }

    QDD res = qdd_makenode(k, low, high);

    // multiply by existing edge weight on qdd
    AMP new_root_amp = Cmul(QDD_AMP(qdd), QDD_AMP(res));
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    return res;
}



/********************</applying (controlled) sub-circuits>*********************/



/***********************************<Grover>***********************************/

QDD 
_qdd_grover_iteration(QDD qdd, BDDVAR n, bool* flag)
{
    LACE_ME; // TODO: change stuff so we don't need to LACE_ME in every function

    // "oracle" call  (apply -1 flag to desired amplitude)
    qdd = qdd_all_control_phase(qdd, n, flag);

    // H on all qubits
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    // Phase on all amplitudes except |000...0>
    bool x[n]; 
    for(BDDVAR k = 0; k < n; k++) x[k] = 0;
    qdd = qdd_all_control_phase(qdd, n, x);
    AMP new_root_amp = Cmul(QDD_AMP(qdd), Clookup(Cmake(-1.0, 0.0)));
    qdd = qdd_bundle_ptr_amp(QDD_PTR(qdd), new_root_amp);

    // H on all qubits
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    return qdd;
}

QDD
qdd_grover(BDDVAR n, bool* flag)
{   
    LACE_ME;

    // not entirely sure about this, book says R <= ceil(pi/4 * sqrt(N))
    uint32_t R = floor( 3.14159265359/4.0 * sqrt( pow(2,n) ) );

    // start with all zero state |000...0>
    QDD qdd = qdd_create_all_zero_state(n);

    // H on all qubits
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    // Grover iterations
    for (uint32_t i = 0; i < R; i++) qdd = _qdd_grover_iteration(qdd, n, flag);

    return qdd;
}

/***********************************</Grover>**********************************/


/*******************************<Shor components>******************************/

uint64_t 
inverse_mod(uint64_t a, uint64_t N) {
    int t = 0;
	int newt = 1;
    int r = N;
    int newr = a;
    int h;
    while(newr != 0) {
        int quotient = r / newr;
        h = t;
        t = newt;
        newt = h - quotient * newt;
        h = r;
        r = newr;
        newr = h - quotient * newr;
    }
    if(r > 1)
        printf("ERROR: a is not invertible\n");
    if(t < 0)
    	t = t + N;
    return t;
}

QDD 
qdd_phi_add(QDD qdd, BDDVAR first, BDDVAR last, bool* a) 
{
    LACE_ME;

    QDD res = qdd;

    int k;
    int num_qubits = (last - first) + 1;
    BDDVAR qubit;
    for (int i = 0; i < num_qubits; i++) {
        qubit = first + i;
        for (int j = i; j < num_qubits; j++){
            if (a[j] == 1) {
                k = (j - i) + 1;
                res = qdd_gate(res, GATEID_Rk(k), qubit);
            }
        }
    }
    return res;
}

QDD 
qdd_phi_add_inv(QDD qdd, BDDVAR first, BDDVAR last, bool* a) 
{
    // These are all phase gates, so they'll all commute, so this is the exact
    // same function als qdd_phi_add() but with inversed angles.
    LACE_ME;

    QDD res = qdd;

    int k;
    int num_qubits = (last - first) + 1;
    BDDVAR qubit;
    for (int i = 0; i < num_qubits; i++) {
        qubit = first + i;
        for (int j = i; j < num_qubits; j++){
            if (a[j] == 1) {
                k = (j - i) + 1;
                res = qdd_gate(res, GATEID_Rk_dag(k), qubit);
            }
        }
    }
    return res;
}

QDD
qdd_phi_add_mod(QDD qdd, BDDVAR* cs, uint64_t a, uint64_t N)
{
    LACE_ME;
    // clear cache (this function is called with different a, and cached results
    // are not parameterized on a)
    sylvan_clear_cache();
    shor_set_globals(a, N); // set bitvalues of a/N (N says the same though)

    // 1.  controlled(c1,c2) phi_add(a)
    qdd = qdd_ccircuit(qdd, CIRCID_phi_add_a, cs, shor_wires.targ_first, shor_wires.targ_last);
    // 2.  phi_add_inv(N)
    qdd = qdd_circuit(qdd, CIRCID_phi_add_N_inv, shor_wires.targ_first, shor_wires.targ_last);
    // 3.  QFT_inv
    qdd = qdd_circuit(qdd, CIRCID_QFT_inv, shor_wires.targ_first, shor_wires.targ_last);
    // 4.  CNOT (control = carry wire? = first of phi ADD, target = helper)
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    qdd = qdd_cgate(qdd, GATEID_Z, shor_wires.helper, shor_wires.targ_first)
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    // 5.  QFT
    qdd = qdd_circuit(qdd, CIRCID_QFT, shor_wires.targ_first, shor_wires.targ_last);
    // 6.  controlled phi_add(N) (control = helper)
    BDDVAR controls[] = {shor_wires.helper, QDD_INVALID_VAR, QDD_INVALID_VAR};
    qdd = qdd_ccircuit(qdd, CIRCID_phi_add_N, controls, shor_wires.targ_first, shor_wires.targ_last);
    // 7. controlled(c1, c2) phi_add_inv(a)
    qdd = qdd_ccircuit(qdd, CIRCID_phi_add_a_inv, cs, shor_wires.targ_first, shor_wires.targ_last);
    // 8.  QFT_inv
    qdd = qdd_circuit(qdd, CIRCID_QFT_inv, shor_wires.targ_first, shor_wires.targ_last);
    // 9.  X on same wire as control of CNOT in 4/10
    qdd = qdd_gate(qdd, GATEID_X, shor_wires.targ_first);
    // 10. CNOT
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    qdd = qdd_cgate(qdd, GATEID_Z, shor_wires.helper, shor_wires.targ_first)
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    // 11. X on same wire as control of CNOT in 4/10
    qdd = qdd_gate(qdd, GATEID_X, shor_wires.targ_first);
    // 12. QFT
    qdd = qdd_circuit(qdd, CIRCID_QFT, shor_wires.targ_first, shor_wires.targ_last); 
    // 13. controlled(c1,c2) phi_add(a)
    qdd = qdd_ccircuit(qdd, CIRCID_phi_add_a, cs, shor_wires.targ_first, shor_wires.targ_last); 

    return qdd;
}

QDD
qdd_phi_add_mod_inv(QDD qdd, BDDVAR* cs, uint64_t a, uint64_t N)
{
    // Inverse of function above
    LACE_ME;
    // clear cache (this function is called with different a, and cached results
    // are not parameterized on a)
    sylvan_clear_cache();
    shor_set_globals(a, N); // set bitvalues of a/N (N says the same though)

    // 13. controlled(c1,c2) phi_add_inv(a)
    qdd = qdd_ccircuit(qdd, CIRCID_phi_add_a_inv, cs, shor_wires.targ_first, shor_wires.targ_last);
    // 12. QFT^-1
    qdd = qdd_circuit(qdd, CIRCID_QFT_inv, shor_wires.targ_first, shor_wires.targ_last);
    // 11. X^-1 = X
    qdd = qdd_gate(qdd, GATEID_X, shor_wires.targ_first);
    // 10. CNOT^-1 = CNOT
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    qdd = qdd_cgate(qdd, GATEID_Z, shor_wires.helper, shor_wires.targ_first)
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    // 9.  X^-1 = X
    qdd = qdd_gate(qdd, GATEID_X, shor_wires.targ_first);
    // 8.  (QFT^-1)^-1 = QFT
    qdd = qdd_circuit(qdd, CIRCID_QFT, shor_wires.targ_first, shor_wires.targ_last);
    // 7.  controlled(c1, c2) phi_add(a)
    qdd = qdd_ccircuit(qdd, CIRCID_phi_add_a, cs, shor_wires.targ_first, shor_wires.targ_last);
    // 6.  controlled phi_add_inv(N) (control = helper)
    BDDVAR controls[] = {shor_wires.helper, QDD_INVALID_VAR, QDD_INVALID_VAR};
    qdd = qdd_ccircuit(qdd, CIRCID_phi_add_N_inv, controls, shor_wires.targ_first, shor_wires.targ_last);
    // 5.  QFT^-1
    qdd = qdd_circuit(qdd, CIRCID_QFT_inv, shor_wires.targ_first, shor_wires.targ_last);
    // 4. CNOT^-1 = CNOT (control = carry wire? = first of phi ADD, target = helper)
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    qdd = qdd_cgate(qdd, GATEID_Z, shor_wires.helper, shor_wires.targ_first)
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    // 3. (QFT^-1)^-1 = QFT
    qdd = qdd_circuit(qdd, CIRCID_QFT, shor_wires.targ_first, shor_wires.targ_last);
    // 2.  phi_add(N)
    qdd = qdd_circuit(qdd, CIRCID_phi_add_N, shor_wires.targ_first, shor_wires.targ_last);
    // 1.  controlled(c1,c2) phi_add_inv(a)
    qdd = qdd_ccircuit(qdd, CIRCID_phi_add_a_inv, cs, shor_wires.targ_first, shor_wires.targ_last);

    return qdd;
}


QDD
qdd_cmult(QDD qdd, uint64_t a, uint64_t N)
{
    // this implements the _controlled_ cmult operation
    // 1. QFT on bottom register
    qdd = qdd_circuit(qdd, CIRCID_QFT, shor_wires.targ_first, shor_wires.targ_last);

    // 2. loop over k = {0, n-1}
    uint64_t t = a;
    BDDVAR cs[] = {shor_wires.top, QDD_INVALID_VAR, QDD_INVALID_VAR};
    for (BDDVAR i = shor_wires.ctrl_first; i <= shor_wires.ctrl_last; i++) { // loop the other way in QMDD imp?
        // 2a. double controlled phi_add_mod(a* 2^k)
        cs[1] = i;
        qdd = qdd_phi_add_mod(qdd, cs, t, N);
        t = (2*t) % N;
    }

    // 3. QFT^-1 on bottom register
    qdd = qdd_circuit(qdd, CIRCID_QFT_inv, shor_wires.targ_first, shor_wires.targ_last);

    return qdd;
}

QDD
qdd_cmult_inv(QDD qdd, uint64_t a, uint64_t N)
{
    // not quite inverse of above
    // 1. QFT on bottom register
    qdd = qdd_circuit(qdd, CIRCID_QFT, shor_wires.targ_first, shor_wires.targ_last);
    
    // 2. same loop over k but with phi_add_mod_inv
    uint64_t t = a;
    BDDVAR cs[] = {shor_wires.top, QDD_INVALID_VAR, QDD_INVALID_VAR};
    for (BDDVAR i = shor_wires.ctrl_first; i <= shor_wires.ctrl_last; i++) { // loop the other way in QMDD imp?
        // 2a. double controlled phi_add_mod_inv(a* 2^k)
        cs[1] = i;
        qdd = qdd_phi_add_mod_inv(qdd, cs, t, N);
        t = (2*t) % N;
    }

    // 3. QFT^-1 on bottom register
    qdd = qdd_circuit(qdd, CIRCID_QFT_inv, shor_wires.targ_first, shor_wires.targ_last);

    return qdd;
}

QDD
qdd_shor_ua(QDD qdd,  uint64_t a, uint64_t N)
{
    LACE_ME;

    // 1. controlled Cmult(a)
    qdd = qdd_cmult(qdd, a, N);

    // 2. controlled swap top/bottom registers
    BDDVAR cs[] = {shor_wires.top, QDD_INVALID_VAR, QDD_INVALID_VAR};
    for (uint32_t i = shor_wires.ctrl_first; i <= shor_wires.ctrl_last; i++) {
        qdd = qdd_ccircuit(qdd, CIRCID_swap, cs, i, shor_wires.targ_first+i); // maybe +i-1
    }

    // 3. controlled Cmult_inv(a^-1)
    uint64_t a_inv = inverse_mod(a, N);
    qdd = qdd_cmult_inv(qdd, a_inv, N);

    return qdd;
}

uint64_t
shor_period_finding(uint64_t a, uint64_t N)
{
    // TODO: circuit (quantum period finding of f(x) = a^x mod N)
    // create QDD |0>|11..1>
    uint32_t num_qubits = 2*shor_n + 3;
    bool x[num_qubits];
    x[0] = 0;
    for (BDDVAR k = 1; k < num_qubits; k++) x[k] = 1;
    QDD qdd = qdd_create_basis_state(num_qubits, x);

    uint64_t as[2*shor_n];
	as[2*shor_n-1] = a;
	uint64_t new_a = a;
	for(int i = 2*shor_n-2; i >= 0; i--) {
		new_a = new_a * new_a;
		new_a = new_a % N;
		as[i] = new_a;
	}

    printf("as:[");
    for (uint32_t i = 0; i < 2*shor_n; i++) printf("%ld, ", as[i]);
    printf("]\n");

    LACE_ME;

    int m_outcomes[2*shor_n];
    int m_outcome;
    double m_prob;

    FILE *fp;

    for (uint32_t i = 0; i < 2*shor_n; i++) {
        printf("shor it = %d/%d\n", i+1, 2*shor_n);
        // H on top wire
        qdd = qdd_gate(qdd, GATEID_H, shor_wires.top);
        // controlled Ua^...
        qdd = qdd_shor_ua(qdd, as[i], N);

        // phase gates based on previous measurement
        int k = 2; // First gate needs to be R^dag(2) = S^dag
        for (int j = i-1; j >= 0; j--) {
            if (m_outcomes[j] == 1)
                qdd = qdd_gate(qdd, GATEID_Rk_dag(k), shor_wires.top);
            k = k << 1; // 2^(iteration)
        }

        // H on top wire
        qdd = qdd_gate(qdd, GATEID_H, shor_wires.top);

        fp = fopen("shor_it.dot", "w");
        qdd_fprintdot(fp, qdd, false);
        fclose(fp);

        // measure q0
        assert(qdd_is_unitvector(qdd, num_qubits));
        qdd = qdd_measure_qubit(qdd, shor_wires.top, num_qubits, &m_outcome, &m_prob);
        m_outcomes[i] = m_outcome;
        assert(qdd_is_unitvector(qdd, num_qubits));
        

        // make sure q0 is in the |0> state
        if (m_outcome == 1) qdd = qdd_gate(qdd, GATEID_X, shor_wires.top);
    }

    // turn measurement outcomes into an integer
    uint64_t res = 0;
    for (uint32_t i = 0; i < 2*shor_n; i++) {
        int index = 2*shor_n-1-i;
        int mmm = m_outcomes[index];
        printf("m_outcomes[%d] = %d\n", index, mmm);
        res = (res << 1) + m_outcomes[index];
    }
    return res;
}

void
shor_set_globals(uint64_t a, uint64_t N) 
{
    shor_n = ceil(log2(N)); // number of bits for N (not the number of qubits!)  
    uint64_t p2 = 1;
    for (uint32_t i = 0; i < shor_n; i++) { // LSB in bits[0], MSB in bits[63]
        shor_bits_a[i] = a & p2;
        shor_bits_N[i] = N & p2;
        p2 = p2 << 1;
    }
    // Set wire numbers
    shor_wires.top        = 0;
    shor_wires.ctrl_first = 1;
    shor_wires.ctrl_last  = shor_n;
    shor_wires.helper     = shor_n + 1; // easier to have this in the middle
    shor_wires.targ_first = shor_n + 2;
    shor_wires.targ_last  = 2*shor_n + 2;
}

uint64_t 
my_gcd (uint64_t a, uint64_t b) // clash with gcd in sylvan_mtbdd.c ...
{
  uint64_t c;
  while ( a != 0 ) { c = a; a = b%a;  b = c; }
  return b;
}

void
run_shor(uint64_t N)
{
    // The classical part
    srand(42);
    uint64_t a;
	do {
		a = rand() % N;
	} while (my_gcd(a, N) != 1 || a == 1);
    // for testing, fix a
    a = 7;  // for a = 11, m=0 or m=8 (?)
            // for a = 7,  m=0, 64, 128, 192

    shor_set_globals(a, N);
    
    printf("input N        = %ld [", N);
    for (uint32_t i=0; i<shor_n; i++) printf("%d,", shor_bits_N[i]);
    printf("]\n");
    printf("n (bits for N) = %d\n",  shor_n);
    printf("random a       = %ld [", a);
    for (uint32_t i=0; i<shor_n; i++) printf("%d,", shor_bits_a[i]);
    printf("]\n\n");

    printf("wires:\n");
    printf("top:        %d\n", shor_wires.top);
    printf("ctrl_first: %d\n", shor_wires.ctrl_first);
    printf("ctrl_last:  %d\n", shor_wires.ctrl_last);
    printf("helper:     %d\n", shor_wires.helper);
    printf("targ_first: %d\n", shor_wires.targ_first);
    printf("targ_last:  %d\n\n", shor_wires.targ_last);

    uint64_t b = shor_period_finding(a, N);
    printf("b = %ld\n", b);
    // for b the following is true: 
    // b/2^l = x/r, 
    // where l = the number of bits required for N, r = the period we want

    // TODO: post processing
}

/******************************</Shor components>******************************/



/***********************<measurements and probabilities>***********************/

// Container for disguising doubles as ints so they can go in Sylvan's cache
typedef union {
    double   as_double;
    uint64_t as_int;
} prob_container_t;

QDD
qdd_measure_q0(QDD qdd, BDDVAR nvars, int *m, double *p)
{  
    LACE_ME;

    // get probabilities for q0 = |0> and q0 = |1>
    double prob_low, prob_high, prob_root;
    qddnode_t node;
    bool skipped = false;
    if (QDD_PTR(qdd) == QDD_TERMINAL) {
        skipped = true;
    }
    else {
        node = QDD_GETNODE(QDD_PTR(qdd));
        if (qddnode_getvar(node) != 0) {
            skipped = true;
        }
    }
    QDD low, high;
    if (skipped) {
        // if skipped q0 is a don't care, treat separately?
        low  = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
        high = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
    }
    else {
        low  = qddnode_getlow(node);
        high = qddnode_gethigh(node);
    }
    prob_low  = qdd_unnormed_prob(low,  1, nvars);
    prob_high = qdd_unnormed_prob(high, 1, nvars);
    prob_root = _prob(QDD_AMP(qdd));
    prob_low  *= prob_root;
    prob_high *= prob_root;
    if (fabs(prob_low + prob_high - 1.0) > TOLERANCE) {
        printf("prob sum = %.5lf\n", prob_low + prob_high);
        assert("probabilities don't sum to 1" && false);
    }

    // flip a coin
    float rnd = ((float)rand())/RAND_MAX;
    *m = (rnd < prob_low) ? 0 : 1;
    *p = prob_low;

    // produce post-measurement state
    AMP norm;
    if (*m == 0) {
        high = QDD_TERMINAL;
        low  = qdd_bundle_ptr_amp(QDD_PTR(low), C_ONE);
        norm = Clookup(Cmake(sqrt(prob_low), 0.0));
    }
    else {
        low  = QDD_TERMINAL;
        high = qdd_bundle_ptr_amp(QDD_PTR(high), C_ONE);
        norm = Clookup(Cmake(sqrt(prob_high), 0.0));
    }

    QDD res = qdd_makenode(0, low, high);
    AMP normalized = Cdiv(QDD_AMP(qdd), norm);
    res = qdd_bundle_ptr_amp(QDD_PTR(res), normalized);
    return res;
}

QDD
qdd_measure_qubit(QDD qdd, BDDVAR k, BDDVAR nvars, int *m, double *p)
{
    if (k == 0) return qdd_measure_q0(qdd, nvars, m, p);
    qdd = qdd_circuit_swap(qdd, 0, k);
    qdd = qdd_measure_q0(qdd, nvars, m, p);
    qdd = qdd_circuit_swap(qdd, 0, k);
    return qdd;
}

QDD
qdd_measure_all(QDD qdd, BDDVAR n, bool* ms, double *p)
{
    LACE_ME;

    qddnode_t node;
    bool skipped;
    BDDVAR var;
    double prob_low, prob_high, prob_path = 1.0, prob_roots = 1.0;

    for (BDDVAR k=0; k < n; k++) {
        // find relevant node (assuming it should be the next one)
        skipped = false;
        if (QDD_PTR(qdd) == QDD_TERMINAL) {
            skipped = true;
        }
        else {
            node = QDD_GETNODE(QDD_PTR(qdd));
            var = qddnode_getvar(node);
            assert(var >= k);
            if (var > k) skipped = true;
        }
        QDD low, high;
        if (skipped) {
            // if skipped q0 is a don't care, treat separately?
            low  = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
            high = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
        }
        else {
            low  = qddnode_getlow(node);
            high = qddnode_gethigh(node);
        }

        prob_low  = qdd_unnormed_prob(low,  k+1, n);
        prob_high = qdd_unnormed_prob(high, k+1, n);
        prob_roots *= _prob(QDD_AMP(qdd));
        prob_high = prob_high * prob_roots / prob_path;
        prob_low  = prob_low  * prob_roots / prob_path;

        if (fabs(prob_low + prob_high - 1.0) > TOLERANCE) {
            printf("prob sum = %.10lf\n", prob_low + prob_high);
            assert("probabilities don't sum to 1" && false);
        }

        // flip a coin
        float rnd = ((float)rand())/RAND_MAX;
        ms[k] = (rnd < prob_low) ? 0 : 1;

        // Get next edge
        qdd        = (ms[k] == 0) ? low : high;
        prob_path *= (ms[k] == 0) ? prob_low : prob_high;
    }

    *p = prob_path;
    return qdd_create_basis_state(n, ms);
}

TASK_IMPL_3(double, qdd_unnormed_prob, QDD, qdd, BDDVAR, topvar, BDDVAR, nvars)
{
    assert(topvar <= nvars);

    if (topvar == nvars) {
        assert(QDD_PTR(qdd) == QDD_TERMINAL);
        return _prob(QDD_AMP(qdd));
    }

    // Look in cache
    bool cachenow = 1;
    if (cachenow) {
        uint64_t prob_bits;
        if (cache_get(CACHE_QDD_PROB, qdd, QDD_PARAM_PACK_16(topvar, nvars), &prob_bits)) {
            sylvan_stats_count(QDD_PROB_CACHED);
            prob_container_t container = (prob_container_t) prob_bits;
            return container.as_double;
        }
    }

    // Check if the node we want is being skipped
    bool skipped = false;
    qddnode_t node;
    if (QDD_PTR(qdd) == QDD_TERMINAL) skipped = true;
    else {
        node = QDD_GETNODE(QDD_PTR(qdd));
        BDDVAR var = qddnode_getvar(node);
        assert (var >= topvar);
        if (var > topvar) skipped = true;
    }

    // get low and high
    QDD low, high;
    if (skipped) {
        low  = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
        high = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
    }
    else {
        low  = qddnode_getlow(node);
        high = qddnode_gethigh(node);
    }

    double prob_low, prob_high, prob_root, prob_res; // "prob" = absolute amps squared
    BDDVAR nextvar = topvar + 1;

    // these are not q/bdds, so no need to protect from gc right?
    SPAWN(qdd_unnormed_prob, high, nextvar, nvars);
    prob_low = CALL(qdd_unnormed_prob, low, nextvar, nvars);
    prob_high = SYNC(qdd_unnormed_prob);
    prob_root = _prob(QDD_AMP(qdd));
    prob_res = prob_root * (prob_low + prob_high);

    // Put in cache and return
    if (cachenow) {
        prob_container_t container = (prob_container_t) prob_res;
        if (cache_put(CACHE_QDD_PROB, qdd, QDD_PARAM_PACK_16(topvar, nvars), container.as_int))
            sylvan_stats_count(QDD_PROB_CACHEDPUT);
    }
    return prob_res;
}

AMP
qdd_get_amplitude(QDD q, bool* basis_state)
{

    AMP a = C_ONE;
    for (;;) {
        a = Cmul(a, QDD_AMP(q));
        
        // if the current edge is pointing to the terminal, we're done.
        if (QDD_PTR(q) == QDD_TERMINAL) break;

        // now we need to choose low or high edge of next node
        qddnode_t node = QDD_GETNODE(QDD_PTR(q));
        BDDVAR var     = qddnode_getvar(node);

        // Condition low/high choice on basis state vector[var]
        if (basis_state[var] == 0)
            q = node->low;
        else
            q = node->high;
    }

    // TODO: return complex struct instead of the index?
    return a;
}

double
_prob(AMP a) // TODO: rename to "abs_squared"
{
    // move to qdd_int file?
    complex_t c = Cvalue(a);
    double abs = sqrt( (c.r*c.r) + (c.i*c.i) );
    return (abs*abs);
}

/***********************<measurements and probabilities>***********************/



/***************************<Initial state creation>***************************/

QDD
qdd_create_all_zero_state(BDDVAR n)
{
    bool x[n];
    for(BDDVAR k=0; k<n; k++) x[k] = 0;
    return qdd_create_basis_state(n, x);
}

QDD
qdd_create_basis_state(BDDVAR n, bool* x)
{
    // start at terminal, and build backwards
    QDD low, high, prev = QDD_TERMINAL;

    for(int k = n-1; k >=0; k--){
        if(x[k] == 0){
            low = qdd_bundle_ptr_amp(QDD_PTR(prev), C_ONE);
            high = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ZERO);
        }
        else if(x[k] == 1){
            low = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ZERO);
            high = qdd_bundle_ptr_amp(QDD_PTR(prev), C_ONE);
        }
        // add node to unique table
        prev = qdd_makenode(k, low, high);
    }
    return prev;
}

/**************************</Initial state creation>***************************/


bool
qdd_equivalent(QDD a, QDD b, int n, bool exact, bool verbose)
{
    bool has_next = true;
    AMP amp_a, amp_b;
    bool x[n];
    for(int k=0; k<n; k++) x[k] = 0;
    while(has_next){
        amp_a = qdd_get_amplitude(a, x);
        amp_b = qdd_get_amplitude(b, x);
        if(exact){
            if(!CexactEqual(Cvalue(amp_a), Cvalue(amp_b))){
                if(verbose){
                    _print_bitstring(x, n, true);
                    printf(", amp a ="); Cprint(Cvalue(amp_a));
                    printf(" != amp b ="); Cprint(Cvalue(amp_b));
                    printf("\n");
                }
                return false;
            }
        }
        else{
            if(!CapproxEqual(Cvalue(amp_a), Cvalue(amp_b))){
                if(verbose){
                    _print_bitstring(x, n, true);
                    printf(", amp a ="); Cprint(Cvalue(amp_a));
                    printf(" !~= amp b ="); Cprint(Cvalue(amp_b));
                    printf("\n");
                }
                return false;
            }
        }
        has_next = _next_bitstring(x, n);
    }
    return true;
}

bool
qdd_is_unitvector(QDD qdd, BDDVAR n)
{
    bool has_next = true;
    AMP a;
    bool x[n];
    for(BDDVAR k=0; k<n; k++) x[k] = 0;

    double sum_abs_squares = 0.0;
    while(has_next){
        a = qdd_get_amplitude(qdd, x);
        sum_abs_squares += _prob(a);
        has_next = _next_bitstring(x, n);
    }

    if (abs(sum_abs_squares - 1.0) < TOLERANCE)
        return true;
    else
        return false;
}

bool
_next_bitstring(bool *x, int n)
{
    // binary add 1
    bool success = false;
    for(int k=0; k<n; k++){
        if(x[k] == 1){
            x[k] = 0;
        }
        else if(x[k] == 0) {
            x[k] = 1;
            success = true;
            break;
        }
    }
    return success;
}

void
_print_bitstring(bool *x, int n, bool backwards)
{
    if (backwards)
        for(int k=n-1; k>=0; k--) printf("%d", x[k]);
    else 
        for(int k=0; k<n; k++) printf("%d", x[k]);
}



/**************************<printing & file writing>***************************/

/**
 * Pretty prints the information contained in `n`.
 */
static void qddnode_pprint(qddnode_t n)
{
    AMP amp_low  = qddnode_getamplow(n);
    AMP amp_high = qddnode_getamphigh(n);
    printf("[var=%d, low=%lx, high=%lx, ", 
             qddnode_getvar(n),
             qddnode_getptrlow(n),
             qddnode_getptrhigh(n));
    if(amp_low == C_ZERO)      printf("a=C_ZERO, ");
    else if(amp_low == C_ONE)  printf("a=C_ONE, ");
    else {
        printf("a=%lx, ",amp_low);
        printf("("); Cprint(Cvalue(amp_high)); printf(")");
    }                      
    if(amp_high == C_ZERO)     printf("b=C_ZERO ");
    else if(amp_high == C_ONE) printf("b=C_ONE, ");
    else {                     
        printf("b=%lx", amp_high);
        printf("("); Cprint(Cvalue(amp_high)); printf(")");
    }
    printf("]\n");
}

void
_print_qdd(QDD q)
{
    if(QDD_PTR(q) != QDD_TERMINAL){
        qddnode_t node = QDD_GETNODE(QDD_PTR(q));
        if(!qddnode_getmark(node)){
            qddnode_setmark(node, 1);
            printf("%lx\t", QDD_PTR(q));
            qddnode_pprint(node);
            _print_qdd(qddnode_getlow(node));
            _print_qdd(qddnode_gethigh(node));
        }
    }
}

void
qdd_printnodes(QDD q)
{
    printf("root edge: %lx, %lx = ",QDD_PTR(q), QDD_AMP(q));
    Cprint(Cvalue(QDD_AMP(q)));
    printf("\n");
    _print_qdd(q);
    qdd_unmark_rec(q);
}

static void
qdd_fprintdot_label(FILE *out, AMP a)
{
    fprintf(out, ", label=\"");
    if (a == C_ONE) {}
    else if (a == C_ZERO) { fprintf(out, "0"); }
    else {
        complex_t val = Cvalue(a);
        if (val.r != 0.0) fprintf(out, "%.3lf", val.r);
        if (val.i > 0.0) fprintf(out, "+%.3lfi", val.i);
        else if (val.i < 0.0) fprintf(out, "-%.3lfi", val.i);
    }
    fprintf(out, "\"");
}

static void
qdd_fprintdot_rec(FILE *out, QDD qdd, bool draw_zeros)
{
    // terminal node
    if(QDD_PTR(qdd) == QDD_TERMINAL) return;

    qddnode_t n = QDD_GETNODE(QDD_PTR(qdd));
    if (qddnode_getmark(n)) return;
    qddnode_setmark(n, 1);

    // add this node
    fprintf(out, "%" PRIu64 " [label=\"%" PRIu32 "\"];\n",
            QDD_PTR(qdd), qddnode_getvar(n));

    
    // children of this node
    qdd_fprintdot_rec(out, qddnode_getlow(n), draw_zeros);
    qdd_fprintdot_rec(out, qddnode_gethigh(n), draw_zeros);

    // add edge from this node to each child (unless weight 0)
    if (draw_zeros || qddnode_getamplow(n) != C_ZERO) {
        fprintf(out, "%" PRIu64 " -> %" PRIu64 " [style=dashed",
                    QDD_PTR(qdd), qddnode_getptrlow(n));
        qdd_fprintdot_label(out, qddnode_getamplow(n));
        fprintf(out, "];\n");
    }
    if (draw_zeros || qddnode_getamphigh(n) != C_ZERO) {
        fprintf(out, "%" PRIu64 " -> %" PRIu64 " [style=solid",
                    QDD_PTR(qdd), qddnode_getptrhigh(n));
        qdd_fprintdot_label(out, qddnode_getamphigh(n));
        fprintf(out, "];\n");
    }
}

void
qdd_fprintdot(FILE *out, QDD qdd, bool draw_zeros)
{
    fprintf(out, "digraph \"DD\" {\n");
    fprintf(out, "center = true;\n");
    fprintf(out, "edge [dir=forward];\n");
    fprintf(out, "root [style=invis];\n");
    fprintf(out, "root -> %" PRIu64 " [style=solid", QDD_PTR(qdd));
    qdd_fprintdot_label(out, QDD_AMP(qdd));
    fprintf(out, "];\n");

    // terminal node
    fprintf(out, "%lu [shape=box, label=\"T\"];\n", QDD_TERMINAL);

    // recursively add nodes
    qdd_fprintdot_rec(out, qdd, draw_zeros);
    qdd_unmark_rec(qdd);

    fprintf(out, "}\n");
}

/**************************</printing & file writing>**************************/
