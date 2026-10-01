"""
The LIMDD arm of run_qasm_on_qmdd on qubits 64 and above.

Built with LIMDD_PAULI_WORDS=2 the runner takes 128 qubits. Until gates took
their controls as a list of qubits, a gate with its target or a control at
qubit 64 or above computed a wrong state there, and no test ran the runner
that high. On a register that large the runner reports no amplitudes, only
the norm and P(q[0] = 0). So a circuit on 8 qubits runs once for each of its
qubits j, with qubit j on q[0] and the other seven on the last seven qubits
of the register, and P(q[0] = 0) must be the marginal of qubit j in the state
vector the QMDD arm computes for the circuit on 8 qubits. The register has
128 qubits, and 68, which puts the seven across the boundary at qubit 64.
Idle qubits below the circuit make each gate slower, so there are none.

The runner is RUN_QASM_ON_QMDD_WIDE, or else ./build-w2/qasm/run_qasm_on_qmdd,
the two-word build of the CI. Without it the test is skipped.
"""
import json
import os
import random
import subprocess
import pytest


SIM_QASM = os.environ.get('RUN_QASM_ON_QMDD_WIDE', './build-w2/qasm/run_qasm_on_qmdd')
SIZES = ['--node-tab-size', '18', '--wgt-tab-size', '18']
HEADER = ["OPENQASM 2.0;", 'include "qelib1.inc";']
K = 8           # qubits of a circuit


def circuit(seed, gates=32):
    """H on every qubit, random gates from what the LIMDD arm takes, and H on
    every qubit again, as (name, qubits). cs, csdg and ccz go through its
    control list, cx, cz and swap through limdd_cgate_either and limdd_swap.
    With no H in between, most marginals at the end are not 1/2, which a
    wrong state would leave alike far more often."""
    rng = random.Random(seed)
    ops = [('h', [j]) for j in range(K)]
    for _ in range(gates):
        q = rng.sample(range(K), 3)
        x = rng.random()
        if x < 0.3:
            ops.append((rng.choice(['t', 'tdg', 's', 'sdg', 'z', 'x']), q[:1]))
        elif x < 0.45:
            ops.append((rng.choice(['cx', 'cz']), q[:2]))
        elif x < 0.5:
            ops.append(('swap', q[:2]))
        elif x < 0.75:
            ops.append((rng.choice(['cs', 'csdg']), q[:2]))
        else:
            ops.append(('ccz', q))
    return ops + [('h', [j]) for j in range(K)]


def write(path, ops, n, place):
    """The circuit with its qubit j on q[place[j]] of an n-qubit register."""
    lines = HEADER + [f"qreg q[{n}];"]
    lines += [f"{g} {','.join(f'q[{place[j]}]' for j in qs)};" for g, qs in ops]
    path.write_text("\n".join(lines) + "\n")


def run(path, dd, vector=False):
    out = subprocess.run([SIM_QASM, str(path), '-d', dd, '-e', 'qisq2', '-s', 'low',
                          *SIZES, *(['--state-vector'] if vector else [])],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    assert out.returncode == 0, out.stderr.decode()
    return json.loads(out.stdout)


@pytest.mark.skipif(not os.path.exists(SIM_QASM), reason=f"no two-word runner at {SIM_QASM}")
@pytest.mark.parametrize("seed", [3, 4])
def test_limdd_above_qubit_64(seed, tmp_path):
    ops = circuit(seed)

    # The marginals of the circuit on its own 8 qubits, from the QMDD arm.
    # Bit j of an index is qubit j. The QMDD arm reports P(q[0] = 0) only
    # after a measurement, so the LIMDD arm's on the same 8 qubits checks it.
    ref = tmp_path / "ref.qasm"
    write(ref, ops, K, list(range(K)))
    prob = [abs(complex(*a)) ** 2 for a in run(ref, 'qmdd', vector=True)['state_vector']]
    p0 = [sum(p for i, p in enumerate(prob) if not (i >> j) & 1) for j in range(K)]
    assert abs(p0[0] - run(ref, 'limdd')['statistics']['first_qubit_measurement_prob']) < 1e-5

    bad = []
    for n in (128, 68):
        for j in range(K):
            rest = iter(range(n - K + 1, n))
            place = [0 if i == j else next(rest) for i in range(K)]
            path = tmp_path / f"wide_{n}_{j}.qasm"
            write(path, ops, n, place)
            stats = run(path, 'limdd')['statistics']
            assert abs(stats['norm'] - 1.0) < 1e-5, (n, j, stats['norm'])
            got = stats['first_qubit_measurement_prob']
            if abs(got - p0[j]) > 1e-5:
                bad.append(f"qubit {j} on q[0], the rest on q[{n - K + 1}..{n - 1}]: "
                           f"P(0) is {got}, should be {p0[j]}")
    assert not bad, "\n".join(bad)
