#include <qsylvan_gates.h>
#include <sylvan_int.h>
#include <sylvan_edge_weights_complex.h>
#include <sylvan_edge_weights_qisq2.h>


static long double Pi;    // set value of global Pi

uint64_t gates[n_predef_gates+256+256][4];

/********************** <dynamic custom rotation gates> ***********************/


// store complex values of dynamic gate to re-initialize gate after gc
complex_t dynamic_gate[4];

uint32_t
GATEID_Rz(fl_t theta)
{
    // clear cache to invalidate cached results for GATEID_dynamic
    sylvan_clear_cache();

    // initialize (and store for gc)
    dynamic_gate[0] = cmake_angle(-theta/2.0, 1);
    dynamic_gate[1] = czero();
    dynamic_gate[2] = czero();
    dynamic_gate[3] = cmake_angle(theta/2.0, 1);
    gates[GATEID_dynamic][0] = weight_lookup(&dynamic_gate[0]); // u00
    gates[GATEID_dynamic][1] = weight_lookup(&dynamic_gate[1]); // u01
    gates[GATEID_dynamic][2] = weight_lookup(&dynamic_gate[2]); // u10
    gates[GATEID_dynamic][3] = weight_lookup(&dynamic_gate[3]); // u11

    // return (temporary) gate_id for this gate
    return GATEID_dynamic;
}


uint32_t
GATEID_Rx(fl_t theta)
{
    // clear cache to invalidate cached results for GATEID_dynamic
    sylvan_clear_cache();

    // initialize (and store for gc)
    dynamic_gate[0] = cmake(flt_cos(theta/2.0), 0.0);
    dynamic_gate[1] = cmake(0.0, -flt_sin(theta/2.0));
    dynamic_gate[2] = cmake(0.0, -flt_sin(theta/2.0));
    dynamic_gate[3] = cmake(flt_cos(theta/2.0), 0.0);
    gates[GATEID_dynamic][0] = weight_lookup(&dynamic_gate[0]); // u00
    gates[GATEID_dynamic][1] = weight_lookup(&dynamic_gate[1]); // u01
    gates[GATEID_dynamic][2] = weight_lookup(&dynamic_gate[2]); // u10
    gates[GATEID_dynamic][3] = weight_lookup(&dynamic_gate[3]); // u11

    // return (temporary) gate_id for this gate
    return GATEID_dynamic;
}

uint32_t
GATEID_Ry(fl_t theta)
{
    // clear cache to invalidate cached results for GATEID_dynamic
    sylvan_clear_cache();

    // initialize (and store for gc)
    dynamic_gate[0] = cmake( flt_cos(theta/2.0), 0.0);
    dynamic_gate[1] = cmake(-flt_sin(theta/2.0), 0.0);
    dynamic_gate[2] = cmake( flt_sin(theta/2.0), 0.0);
    dynamic_gate[3] = cmake( flt_cos(theta/2.0), 0.0);
    gates[GATEID_dynamic][0] = weight_lookup(&dynamic_gate[0]); // u00
    gates[GATEID_dynamic][1] = weight_lookup(&dynamic_gate[1]); // u01
    gates[GATEID_dynamic][2] = weight_lookup(&dynamic_gate[2]); // u10
    gates[GATEID_dynamic][3] = weight_lookup(&dynamic_gate[3]); // u11

    // return (temporary) gate_id for this gate
    return GATEID_dynamic;
}

uint32_t
GATEID_Phase(fl_t theta)
{
    // clear cache to invalidate cached results for GATEID_dynamic
    sylvan_clear_cache();
    
    // initialize (and store for gc)
    dynamic_gate[0] = cmake(1.0, 0.0);
    dynamic_gate[1] = cmake(0.0, 0.0);
    dynamic_gate[2] = cmake(0.0, 0.0);
    dynamic_gate[3] = cmake_angle(theta, 1);
    gates[GATEID_dynamic][0] = weight_lookup(&dynamic_gate[0]); // u00
    gates[GATEID_dynamic][1] = weight_lookup(&dynamic_gate[1]); // u01
    gates[GATEID_dynamic][2] = weight_lookup(&dynamic_gate[2]); // u10
    gates[GATEID_dynamic][3] = weight_lookup(&dynamic_gate[3]); // u11

    // return (temporary) gate_id for this gate
    return GATEID_dynamic;
}

uint32_t
GATEID_U(fl_t theta, fl_t phi, fl_t lambda)
{
    // clear cache to invalidate cached results for GATEID_dynamic
    sylvan_clear_cache();

    // initialize (and store for gc)
    dynamic_gate[0] = cmake(flt_cos(theta/2.0), 0.0);
    dynamic_gate[1] = cmul(cmake_angle(lambda,1), cmake(-flt_sin(theta/2.0), 0));
    dynamic_gate[2] = cmul(cmake_angle(phi,1), cmake(flt_sin(theta/2.0), 0));
    dynamic_gate[3] = cmul(cmake_angle(phi+lambda,1), cmake(flt_cos(theta/2.0), 0));
    gates[GATEID_dynamic][0] = weight_lookup(&dynamic_gate[0]); // u00
    gates[GATEID_dynamic][1] = weight_lookup(&dynamic_gate[1]); // u01
    gates[GATEID_dynamic][2] = weight_lookup(&dynamic_gate[2]); // u10
    gates[GATEID_dynamic][3] = weight_lookup(&dynamic_gate[3]); // u11

    // return (temporary) gate_id for this gate
    return GATEID_dynamic;
}

/********************* </dynamic custom rotation gates> ***********************/


/*************************** <dynamic custom gates> ***************************/
void
qmdd_gates_init()
{
    Pi = 2.0 * flt_acos(0.0);

    // initialize 2x2 gates (complex values from gates currently stored in 
    // same table as complex amplitude values)
    uint32_t k;

    k = GATEID_I;
    gates[k][0] = EVBDD_ONE;  gates[k][1] = EVBDD_ZERO;
    gates[k][2] = EVBDD_ZERO; gates[k][3] = EVBDD_ONE;

    k = GATEID_proj0;
    gates[k][0] = EVBDD_ONE;  gates[k][1] = EVBDD_ZERO;
    gates[k][2] = EVBDD_ZERO; gates[k][3] = EVBDD_ZERO;

    k = GATEID_proj1;
    gates[k][0] = EVBDD_ZERO; gates[k][1] = EVBDD_ZERO;
    gates[k][2] = EVBDD_ZERO; gates[k][3] = EVBDD_ONE;

    k = GATEID_X;
    gates[k][0] = EVBDD_ZERO; gates[k][1] = EVBDD_ONE;
    gates[k][2] = EVBDD_ONE;  gates[k][3] = EVBDD_ZERO;

    k = GATEID_Y;
    gates[k][0] = EVBDD_ZERO; gates[k][1] = complex_lookup(0.0, -1.0);
    gates[k][2] = complex_lookup(0.0, 1.0);  gates[k][3] = EVBDD_ZERO;

    k = GATEID_Z;
    gates[k][0] = EVBDD_ONE;  gates[k][1] = EVBDD_ZERO;
    gates[k][2] = EVBDD_ZERO; gates[k][3] = EVBDD_MIN_ONE;

    k = GATEID_H;
    gates[k][0] = gates[k][1] = gates[k][2] = complex_lookup(1.0/flt_sqrt(2.0),0);
    gates[k][3] = complex_lookup(-1.0/flt_sqrt(2.0),0);

    k = GATEID_S;
    gates[k][0] = EVBDD_ONE;  gates[k][1] = EVBDD_ZERO;
    gates[k][2] = EVBDD_ZERO; gates[k][3] = complex_lookup(0.0, 1.0);

    k = GATEID_Sdag;
    gates[k][0] = EVBDD_ONE;  gates[k][1] = EVBDD_ZERO;
    gates[k][2] = EVBDD_ZERO; gates[k][3] = complex_lookup(0.0, -1.0);

    k = GATEID_T;
    gates[k][0] = EVBDD_ONE;  gates[k][1] = EVBDD_ZERO;
    gates[k][2] = EVBDD_ZERO; gates[k][3] = complex_lookup(1.0/flt_sqrt(2.0), 1.0/flt_sqrt(2.0));

    k = GATEID_Tdag;
    gates[k][0] = EVBDD_ONE;  gates[k][1] = EVBDD_ZERO;
    gates[k][2] = EVBDD_ZERO; gates[k][3] = complex_lookup(1.0/flt_sqrt(2.0), -1.0/flt_sqrt(2.0));

    k = GATEID_sqrtX;
    gates[k][0] = complex_lookup(0.5, 0.5); gates[k][1] = complex_lookup(0.5,-0.5);
    gates[k][2] = complex_lookup(0.5,-0.5); gates[k][3] = complex_lookup(0.5, 0.5);

    k = GATEID_sqrtXdag;
    gates[k][0] = complex_lookup(0.5,-0.5); gates[k][1] = complex_lookup(0.5, 0.5);
    gates[k][2] = complex_lookup(0.5, 0.5); gates[k][3] = complex_lookup(0.5,-0.5);

    k = GATEID_sqrtY;
    gates[k][0] = complex_lookup(0.5, 0.5); gates[k][1] = complex_lookup(-0.5,-0.5);
    gates[k][2] = complex_lookup(0.5, 0.5); gates[k][3] = complex_lookup(0.5, 0.5);

    k = GATEID_sqrtYdag;
    gates[k][0] = complex_lookup(0.5,-0.5); gates[k][1] = complex_lookup(0.5,-0.5);
    gates[k][2] = complex_lookup(-0.5,0.5); gates[k][3] = complex_lookup(0.5,-0.5);

    qmdd_phase_gates_init(255);

    // init dynamic gate 
    // (necessary when qmdd_gates_init() is called after gc to re-init all gates)
    k = GATEID_dynamic;
    gates[k][0] = weight_lookup(&dynamic_gate[0]);
    gates[k][1] = weight_lookup(&dynamic_gate[1]);
    gates[k][2] = weight_lookup(&dynamic_gate[2]);
    gates[k][3] = weight_lookup(&dynamic_gate[3]);
}

/**
 * The R_k gates for the qisq2 backend.
 *
 * R_k is diag(1, e^(2*pi*i/2^k)), and e^(2*pi*i/2^k) lies in Q[i,sqrt2] only
 * for k <= 3: 1, -1, i, and (1+i)/sqrt2. Beyond that the root of unity is not
 * in the field at all, so there is no qisq2 number to store and the gate
 * cannot be represented -- this is a property of the field, not a gap in the
 * implementation.
 *
 * The unrepresentable ones are left as the zero matrix rather than filled from
 * a complex_t, which is what qmdd_phase_gates_init would do: passing a
 * complex_t to the qisq2 table reinterprets two doubles as four GMP rationals
 * and corrupts memory. A zero gate gives a visibly wrong answer instead of an
 * unpredictable one.
 */
void
qmdd_phase_gates_qisq2_init(int n)
{
    for (int k = 0; k <= n; k++) {
        const uint32_t fwd = GATEID_Rk(k);
        const uint32_t bwd = GATEID_Rk_dag(k);

        gates[fwd][0] = EVBDD_ONE;  gates[fwd][1] = EVBDD_ZERO;
        gates[fwd][2] = EVBDD_ZERO;
        gates[bwd][0] = EVBDD_ONE;  gates[bwd][1] = EVBDD_ZERO;
        gates[bwd][2] = EVBDD_ZERO;

        switch (k) {
        case 0:  // e^(2 pi i) = 1
            gates[fwd][3] = EVBDD_ONE;
            gates[bwd][3] = EVBDD_ONE;
            break;
        case 1:  // e^(i pi) = -1
            gates[fwd][3] = EVBDD_MIN_ONE;
            gates[bwd][3] = EVBDD_MIN_ONE;
            break;
        case 2:  // e^(i pi/2) = i
            gates[fwd][3] = qisq2_lookup(0,1, 0,1,  1,1, 0,1);
            gates[bwd][3] = qisq2_lookup(0,1, 0,1, -1,1, 0,1);
            break;
        case 3:  // e^(i pi/4) = (1+i)/sqrt2
            gates[fwd][3] = qisq2_lookup(0,1, 1,2, 0,1,  1,2);
            gates[bwd][3] = qisq2_lookup(0,1, 1,2, 0,1, -1,2);
            break;
        default: // not in Q[i,sqrt2]
            gates[fwd][0] = EVBDD_ZERO;
            gates[bwd][0] = EVBDD_ZERO;
            gates[fwd][3] = EVBDD_ZERO;
            gates[bwd][3] = EVBDD_ZERO;
            break;
        }
    }
}

void
qmdd_phase_gates_init(int n)
{
    // add gate R_k to gates table
    // (note that R_0 = I, R_1 = Z, R_2 = S, R_4 = T)
    uint32_t gate_id;
    fl_t angle;
    complex_t cartesian;
    for (int k=0; k<=n; k++) {
        // forward rotation
        angle = 2*Pi / (fl_t)(1<<k);
        cartesian = cmake_angle(angle, 1);
        gate_id = GATEID_Rk(k);
        gates[gate_id][0] = EVBDD_ONE;  gates[gate_id][1] = EVBDD_ZERO;
        gates[gate_id][2] = EVBDD_ZERO; gates[gate_id][3] = weight_lookup(&cartesian);

        // backward rotation
        angle = -2*Pi / (fl_t)(1<<k);
        cartesian = cmake_angle(angle, 1);
        gate_id = GATEID_Rk_dag(k);
        gates[gate_id][0] = EVBDD_ONE;  gates[gate_id][1] = EVBDD_ZERO;
        gates[gate_id][2] = EVBDD_ZERO; gates[gate_id][3] = weight_lookup(&cartesian);
    }
}


// ---------------- < gate definitions for qisq2 > ----------------

void
qmdd_gates_qisq2_init()
{
        // initialize 2x2 gates (complex values from gates currently stored in 
    // same table as complex amplitude values)
    uint32_t k;

    k = GATEID_I;
    gates[k][0] = EVBDD_ONE;  gates[k][1] = EVBDD_ZERO;
    gates[k][2] = EVBDD_ZERO; gates[k][3] = EVBDD_ONE;

    k = GATEID_proj0;
    gates[k][0] = EVBDD_ONE;  gates[k][1] = EVBDD_ZERO;
    gates[k][2] = EVBDD_ZERO; gates[k][3] = EVBDD_ZERO;

    k = GATEID_proj1;
    gates[k][0] = EVBDD_ZERO; gates[k][1] = EVBDD_ZERO;
    gates[k][2] = EVBDD_ZERO; gates[k][3] = EVBDD_ONE;

    k = GATEID_X;
    gates[k][0] = EVBDD_ZERO; gates[k][1] = EVBDD_ONE;
    gates[k][2] = EVBDD_ONE;  gates[k][3] = EVBDD_ZERO;

    k = GATEID_Y;
    gates[k][0] = EVBDD_ZERO; gates[k][1] = qisq2_lookup(0,1,0,1,-1,1,0,1);
    gates[k][2] = qisq2_lookup(0,1,0,1, 1,1,0,1);  gates[k][3] = EVBDD_ZERO;

    k = GATEID_Z;
    gates[k][0] = EVBDD_ONE;  gates[k][1] = EVBDD_ZERO;
    gates[k][2] = EVBDD_ZERO; gates[k][3] = EVBDD_MIN_ONE;

    k = GATEID_H;
    gates[k][0] = gates[k][1] = gates[k][2] = qisq2_lookup(0,1,1,2,0,1,0,1);
    gates[k][3] = qisq2_lookup(0,1,-1,2,0,1,0,1);;

    k = GATEID_S;
    gates[k][0] = EVBDD_ONE;  gates[k][1] = EVBDD_ZERO;
    gates[k][2] = EVBDD_ZERO; gates[k][3] = qisq2_lookup(0,1,0,1,1,1,0,1);

    k = GATEID_Sdag;
    gates[k][0] = EVBDD_ONE;  gates[k][1] = EVBDD_ZERO;
    gates[k][2] = EVBDD_ZERO; gates[k][3] = qisq2_lookup(0,1,0,1,-1,1,0,1);

    k = GATEID_T;
    gates[k][0] = EVBDD_ONE;  gates[k][1] = EVBDD_ZERO;
    gates[k][2] = EVBDD_ZERO; gates[k][3] = qisq2_lookup(0,1,1,2,0,1,1,2);

    k = GATEID_Tdag;
    gates[k][0] = EVBDD_ONE;  gates[k][1] = EVBDD_ZERO;
    gates[k][2] = EVBDD_ZERO; gates[k][3] = qisq2_lookup(0,1,1,2,0,1,-1,2);

    qisq2_lookup(1,2,0,1,1,2,0,1);

    k = GATEID_sqrtX;
    gates[k][0] = qisq2_lookup(1,2,0,1, 1,2,0,1); gates[k][1] = qisq2_lookup(1,2,0,1,-1,2,0,1);
    gates[k][2] = qisq2_lookup(1,2,0,1,-1,2,0,1); gates[k][3] = qisq2_lookup(1,2,0,1,1,2,0,1);

    k = GATEID_sqrtXdag;
    gates[k][0] = qisq2_lookup(1,2,0,1,-1,2,0,1); gates[k][1] = qisq2_lookup(1,2,0,1,1,2,0,1);
    gates[k][2] = qisq2_lookup(1,2,0,1,1,2,0,1); gates[k][3] = qisq2_lookup(1,2,0,1,-1,2,0,1);

    k = GATEID_sqrtY;
    gates[k][0] = qisq2_lookup(1,2,0,1,1,2,0,1); gates[k][1] = qisq2_lookup(-1,2,0,1,-1,2,0,1);
    gates[k][2] = qisq2_lookup(1,2,0,1,1,2,0,1); gates[k][3] = qisq2_lookup( 1,2,0,1, 1,2,0,1);

    k = GATEID_sqrtYdag;
    gates[k][0] = qisq2_lookup(1,2,0,1,-1,2,0,1); gates[k][1] = qisq2_lookup(1,2,0,1,-1,2,0,1);
    gates[k][2] = qisq2_lookup(-1,2,0,1,1,2,0,1); gates[k][3] = qisq2_lookup(1,2,0,1,-1,2,0,1);

    qmdd_phase_gates_qisq2_init(255);
}

// ---------------- </ gate definitions for qisq2 > ----------------