#include "circuit_to_Sylvan.h"

#include "sylvan.h"
#include "sylvan_qdd_complex.h"

void final_measure(QDD qdd, bool* measurements, BDDVAR nvars, BDDVAR runs, bool show)
{
    // Initialise variables
    bool *ms = malloc(nvars * sizeof(bool));
    double *p = malloc(nvars * sizeof(double));
    int *hits = malloc(pow(2,nvars) * sizeof(int));
    for(int i = 0; i < pow(2,nvars); i++) { hits[i] = 0; }

    // Run the circuit <runs> times
    for(BDDVAR i = 0; i < runs; i++) { 
        qdd_measure_all(qdd, nvars, ms, p);
        hits[bitarray_to_int(ms, nvars, false)]++;
    }

    // Reformat circuit results based on what qubits were measured
    BDDVAR index, j, sum = 0;
    for (BDDVAR i = 0; i < nvars; i++) { if (measurements[i]) sum++; }
    int *probs = malloc(pow(2, sum) * sizeof(int));
    for (int k = 0; k < pow(2, sum); k++) probs[k] = 0;
    for (int k = 0; k < (1 << (nvars)); k++) {
        bool *x = int_to_bitarray(k, nvars, false);
        j = sum-1;
        index = 0;
        for (BDDVAR i = nvars; i > 0; i--)
            if (measurements[i-1]) { index += x[i-1] * pow(2, j--); }
        probs[index] += hits[k];
        free(x);
    }

    // Print reformatted circuit results if <show> is toggled
    if (show) {
        for(int k = 0; k < pow(2, sum); k++) {
            if (probs[k] != 0) {
                bool *x = int_to_bitarray(k, nvars, false);
                j = sum-1;
                for (BDDVAR i = nvars; i > 0 ; i--) {
                    if (measurements[i-1]) {
                        printf("%d", x[j]);
                        j--;
                    }
                    else
                        printf("_");
                }
                printf(": %d\n", probs[k]);
                free(x);
            }
        }
    }

    // Free variables
    free(ms);
    free(p);
    free(hits);
    free(probs);
}

TASK_IMPL_1(BDDVAR, get_gate_id, Gate, gate)
{
    // Initialise variables
    BDDVAR gate_id;

    // Switch case on the id of the gate_struct
    if (gate.id == gate_I.id)
        gate_id = GATEID_I;
    else if (gate.id == gate_X.id)
        gate_id = GATEID_X;
    else if (gate.id == gate_Y.id)
        gate_id = GATEID_Y;
    else if (gate.id == gate_Z.id)
        gate_id = GATEID_Z;
    else if (gate.id == gate_H.id)
        gate_id = GATEID_H;
    else if (gate.id == gate_sX.id)
        gate_id = GATEID_sqrtX;
    else if (gate.id == gate_sY.id)
        gate_id = GATEID_sqrtY;
    else if (gate.id == gate_S.id)
        gate_id = GATEID_S;
    else if (gate.id == gate_Sd.id)
        gate_id = GATEID_Sdag;
    else if (gate.id == gate_T.id)
        gate_id = GATEID_T;
    else if (gate.id == gate_Td.id)
        gate_id = GATEID_Tdag;
    // Rx Ry and Rz gates have arbitrary rotations, Sylvan GATEIDs are allocated dynamically
    else if (gate.id == gate_Rx.id)
        gate_id = GATEID_Rx(gate.rotation);
    else if (gate.id == gate_Ry.id)
        gate_id = GATEID_Ry(gate.rotation);
    else if (gate.id == gate_Rz.id)
        gate_id = GATEID_Rz(gate.rotation);
    else {
        fprintf(stderr, "Unknown gate: %d\n", gate.id);
        exit(EXIT_FAILURE);
    }
    return gate_id;
}

TASK_IMPL_4(QDD, apply_controlled_gate, QDD, qdd, Gate, gate, BDDVAR, k, BDDVAR, n)
{
    // Create controlled gate qdd
    QDD qdd_column = handle_control_matrix(gate, k, n);
    // multiply with vector
    qdd = qdd_matvec_mult(qdd_column, qdd, n);
    return qdd;
}

TASK_IMPL_4(QDD, apply_gate, QDD, qdd, Gate, gate, BDDVAR, i, BDDVAR, n)
{
    // Get the corresponding Sylvan GATEID of <gate>
    BDDVAR gate_id = get_gate_id(gate);
    // Apply <gate> using the corresponding controls (if any)
    if (gate.controlSize == 0)
        qdd = qdd_gate(qdd, gate_id, i);
    // else if (gate.controlSize > 0 && gate.controlSize <= 4) {
    else
        qdd = apply_controlled_gate(qdd, gate, i, n);
    // }
    // else {
    //     // More than 3 controls are not fully implemented yet, easier to throw on all cases
    //     fprintf(stderr, "Control-gates with more than 3 controls not implemented.\n");
    //     exit(EXIT_FAILURE);
    // }
    return qdd;
}

TASK_IMPL_3(QDD, handle_control_matrix, Gate, gate, BDDVAR, k, BDDVAR, n)
{
    // Initialise variables
    QDD qdd, qdd_next;
    bool sorted = true;
    BDDVAR temp, gate_id, control = k, index = 0;
    int *c_options = malloc(n * sizeof(int));
    for (BDDVAR i = 0; i < n; i++) c_options[i] = -1;
    // Check if target is below controls
    for (BDDVAR j = 0; j < gate.controlSize; j++) {
        if (gate.control[j] > k)
            sorted = false;
        if (gate.control[j] > control) {
            control = gate.control[j];
            index = j;
        }
    }
    qdd = qdd_create_all_identity_matrix(n);
    // If target is not below some controls, swap target with lowest control
    if (!sorted) {
        qdd = circuit_swap_matrix(k, control, n);
        temp = k;
        k = gate.control[index];
        gate.control[index] = temp;
    }
    // Get the corresponding Sylvan GATEID from <gate>
    gate_id = get_gate_id(gate);
    // Set the indexes of all controls of <gate> to 1
    for (BDDVAR i = 0; i < gate.controlSize; i++)
        c_options[gate.control[i]] = 1;
    // Set the index the target of <gate> to 2
    c_options[k] = 2;
    // Create the gate QDD
    qdd_next = qdd_create_multi_cgate_rec(n, c_options, gate_id, 0);
    qdd = qdd_matmat_mult(qdd_next, qdd, n);
    if (!sorted) {
        qdd_next = circuit_swap_matrix(temp, k, n);
        qdd = qdd_matmat_mult(qdd_next, qdd, n);
        temp = k;
        k = gate.control[index];
        gate.control[index] = temp;
    }
    return qdd;
}

TASK_IMPL_3(QDD, circuit_swap_matrix, BDDVAR, qubit1, BDDVAR, qubit2, BDDVAR, n)
{
    // Initialise variables
    QDD qdd, qdd_next;
    int *c_options = malloc(n * sizeof(int));
    for (BDDVAR i = 0; i < n; i++) c_options[i] = -1;
    // Set qubit1 to 1
    c_options[qubit1] = 1;
    // Set qubit2 to 2
    c_options[qubit2] = 2;

    // Perform swap
    qdd = qdd_create_multi_cgate_rec(n, c_options, GATEID_X, 0);
    qdd_next = qdd_create_single_qubit_gate(n, qubit1, GATEID_H);
    qdd = qdd_matmat_mult(qdd_next, qdd, n);
    qdd_next = qdd_create_multi_cgate_rec(n, c_options, GATEID_Z, 0);
    qdd = qdd_matmat_mult(qdd_next, qdd, n);
    qdd_next = qdd_create_single_qubit_gate(n, qubit1, GATEID_H);
    qdd = qdd_matmat_mult(qdd_next, qdd, n);
    qdd_next = qdd_create_multi_cgate_rec(n, c_options, GATEID_X, 0);
    qdd = qdd_matmat_mult(qdd_next, qdd, n);
    return qdd;
}

QDD run_circuit_matrix(C_struct c_s, BDDVAR runs, BDDVAR limit, bool show)
{
    // Inisialise variables
    LACE_ME;
    Gate gate;
    bool* measurements = malloc(c_s.nvars * sizeof(bool));
    for (BDDVAR i = 0; i < c_s.nvars; i++) measurements[i] = false;
    BDDVAR gate_id;
    BDDVAR *gateids = malloc(c_s.nvars * sizeof(BDDVAR));
    QDD qdd, qdd_column;
    QDD vec = qdd_create_all_zero_state(c_s.nvars);
    qdd = qdd_create_single_qubit_gates_same(c_s.nvars, GATEID_I);
    qdd_column = qdd;

    // Loop over all places in the circuit
    for (BDDVAR j = 0; j < c_s.depth; j++) {
        for (BDDVAR i = 0; i < c_s.nvars; i++) {
            gate = c_s.circuit[i][j];
            // Preset gateid (overwritten if needed)
            gateids[i] = GATEID_I;
            // Skip barrier (does not affect runs)
            // Skip control gates (controls are used by target gate)
            if (gate.id == gate_barrier.id || gate.id == gate_ctrl.id || gate.id == gate_ctrl_c.id || gate.id == gate_I.id)
                continue;
            // Set measurement flag
            else if (gate.id == gate_measure.id)
                measurements[i] = true;
            // Handle controlled gates
            else if (gate.controlSize != 0) {
                qdd_column = handle_control_matrix(gate, i, c_s.nvars);
                // Separately multiply with gate QDD (not together with column)
                qdd = qdd_matmat_mult(qdd_column, qdd, c_s.nvars);
            }
            // Set single gate in gateids
            else {
                gate_id = get_gate_id(gate);
                gateids[i] = gate_id;
            }
        }
        // Skip barrier (whole column has no effect)
        if (gate.id != gate_barrier.id) {
            // Create gate QDD out of current column (gateids)
            qdd_column = qdd_create_single_qubit_gates(c_s.nvars, gateids);
            // Multiply with previous columns
            qdd = qdd_matmat_mult(qdd_column, qdd, c_s.nvars);
        }
        // If nodecount of column QDD is over limit...
        if (qdd_countnodes(qdd) > limit) {
            // Multiply with statevector QDD and reset column QDD
            vec = qdd_matvec_mult(qdd, vec, c_s.nvars);
            qdd = qdd_create_single_qubit_gates_same(c_s.nvars, GATEID_I);
        }
    }
    // Final multiply with statevector QDD
    vec = qdd_matvec_mult(qdd, vec, c_s.nvars);
    // Measure all qubits
    final_measure(vec, measurements, c_s.nvars, runs, show);
    // Free variables
    free(measurements);
    free(gateids);
    return vec;
}

QDD run_c_struct(C_struct c_s, BDDVAR runs, bool show)
{
    // Inisialise variables
    LACE_ME;
    Gate gate;
    bool* measurements = malloc(c_s.nvars * sizeof(bool));
    for (BDDVAR i = 0; i < c_s.nvars; i++) measurements[i] = false;
    QDD qdd = qdd_create_all_zero_state(c_s.nvars);

    // Loop over all places in the circuit
    for (BDDVAR j = 0; j < c_s.depth; j++) {
        for (BDDVAR i = 0; i < c_s.nvars; i++) {
            gate = c_s.circuit[i][j];
            // Skip barrier (does not affect runs)
            // Skip control gates (controls are used by target gate)
            if (gate.id == gate_barrier.id || gate.id == gate_ctrl.id || gate.id == gate_ctrl_c.id || gate.id == gate_I.id)
                continue;
            // Set measurement flag
            else if (gate.id == gate_measure.id)
                measurements[i] = true;
            // Apply gate
            else
                qdd = apply_gate(qdd, gate, i, c_s.nvars);
        }
    }
    // Measure all qubits
    final_measure(qdd, measurements, c_s.nvars, runs, show);
    // Free variables
    free(measurements);
    return qdd;
}

/**
 * Runs QASM circuit given by <filename> and prints the results.
 * Using the -m flag activates gate-gate multiplication runs, gate-statevector runs are used otherwise.
 * 
 * PARAMETERS:
 * - filename: the path to the file containing the QASM circuit code
 * 
 * FLAGS:
 * [-r runs] (optional) the number of runs to perform
 * [-s seed] (optional) the randomness seed to be used
 * [-m matrix] (optional) the boundaray value of nodes in a tree before multiplying with the state vector
 * [-o optimize] (optional) optimize the circuit if true. This option will remove negating gates before running
 * 
 * NOTE:
 * Since multiplying a gate-QDD with a gate-QDD is more expensive than multiplying a gate-QDD with
 * a statevector-QDD, gate-gate multiplication is usually slower. However, circuits which use a lot of 
 * uncomputation can lead to smaller resulting gate QDDs and possibly lead to faster runs.
 */
int main(int argc, char *argv[])
{
    // Initialise flag parameters
    char *filename = argv[1];
    BDDVAR runs = 100;
    BDDVAR seed = 100;
    BDDVAR matrix = 0;
    bool optimize = false;
    int opt;

    // Read flags from cmd and set parameters
    while((opt = getopt(argc, argv, "s:r:m:o")) != -1) {
        switch(opt) {
            case 'r':
                runs = atoi(optarg);
                break;
            case 's':
                seed = atoi(optarg);
                break;
            case 'm':
                matrix = abs(atoi(optarg));
                break;
            case 'o':
                optimize = true;
                break;
            default:
                fprintf(stderr, "usage: %s file [-r runs][-s seed][-m matrix_node_limit][-o optimize]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    // Set randomness seed
    srand(seed);
    // Check if a file is given, if not, return an error
    if(access(filename, F_OK) != 0)
    {
        fprintf(stderr, "Invalid QASM file.\n");
        exit(EXIT_FAILURE);
    }

    // Standard Lace initialization
    int workers = 1;
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<16, -1, COMP_HASHMAP, NORM_LARGEST);
    qdd_set_testing_mode(true); // turn on internal sanity tests

    // Create a circuit struct representing the QASM circuit in the given file
    C_struct c_s = make_c_struct(filename, optimize);
    // Run the circuit based on method
    if (matrix != 0)
        run_circuit_matrix(c_s, runs, matrix, true);
    else
        run_c_struct(c_s, runs, true);
    // Free variables
    delete_c_struct(&c_s);
    sylvan_quit();
    lace_exit();

    return 0;
}