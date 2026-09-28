"""
The BQD and LIMDD arms of run_qasm_on_qmdd against the EVDD arm.

Neither a BQD nor a LIMDD prints its state the way a QMDD does; -v decodes it,
and the decoded vector must be the EVDD's, amplitude for amplitude, on the
Clifford+T gate set the two arms share, Hadamards anywhere. The BQD arm must
also refuse what it cannot do reliably: float weights once a gate can cancel.
"""
import os
import random
import subprocess
import json
import pytest


SIM_QASM = os.environ.get('RUN_QASM_ON_QMDD', './build/qasm/run_qasm_on_qmdd')
SIZES = ['--node-tab-size', '18', '--wgt-tab-size', '18']
HEADER = ["OPENQASM 2.0;", 'include "qelib1.inc";']


def clifford_t_circuit(n, gates, seed):
    """Random gates from the arms' set. cy and ch only with the control above
    the target, where both arms take them; cx, cz and swap either way."""
    rng = random.Random(seed)
    lines = HEADER + [f"qreg q[{n}];"]
    one = ['h', 'h', 'x', 'y', 'z', 's', 'sdg', 't', 'tdg']
    for _ in range(gates):
        k = rng.randrange(6) if n > 1 else 0
        if k <= 2 or n == 1:
            lines.append(f"{rng.choice(one)} q[{rng.randrange(n)}];")
        elif k == 3:
            a, b = rng.sample(range(n), 2)
            lines.append(f"{rng.choice(['cx', 'cz', 'cs', 'csdg', 'swap'])} q[{a}],q[{b}];")
        elif k == 4:
            a, b = sorted(rng.sample(range(n), 2))
            lines.append(f"{rng.choice(['cy', 'ch'])} q[{a}],q[{b}];")
        elif n > 2:
            a, b, c = rng.sample(range(n), 3)
            lines.append(f"ccz q[{a}],q[{b}],q[{c}];")
    return "\n".join(lines)


def clifford_circuit(n, gates, seed):
    """Clifford gates only, Hadamards anywhere: stabiliser states, whose exact
    amplitudes stay small however long the circuit runs."""
    rng = random.Random(seed)
    lines = HEADER + [f"qreg q[{n}];"]
    one = ['h', 'h', 'x', 'y', 'z', 's', 'sdg']
    for _ in range(gates):
        if rng.randrange(2):
            lines.append(f"{rng.choice(one)} q[{rng.randrange(n)}];")
        else:
            a, b = rng.sample(range(n), 2)
            lines.append(f"{rng.choice(['cx', 'cz', 'swap'])} q[{a}],q[{b}];")
    return "\n".join(lines)


def iqp_circuit(n, gates, seed):
    """A Hadamard on every qubit, then level-3 diagonal gates."""
    rng = random.Random(seed)
    lines = HEADER + [f"qreg q[{n}];"] + [f"h q[{q}];" for q in range(n)]
    for _ in range(gates):
        k = rng.randrange(4) if n > 2 else rng.randrange(2)
        if k == 0:
            lines.append(f"{rng.choice(['t', 'tdg', 'z', 's', 'sdg'])} q[{rng.randrange(n)}];")
        elif k == 1 and n > 1:
            a, b = rng.sample(range(n), 2)
            lines.append(f"{rng.choice(['cz', 'cs', 'csdg'])} q[{a}],q[{b}];")
        elif n > 2:
            a, b, c = rng.sample(range(n), 3)
            lines.append(f"ccz q[{a}],q[{b}],q[{c}];")
    return "\n".join(lines)


def run(path, dd, backend, check=True, sizes=SIZES):
    out = subprocess.run([SIM_QASM, path, '-d', dd, '-e', backend, '-s', 'low',
                          '--state-vector', '-m', *sizes],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if not check:
        return out
    assert out.returncode == 0, out.stderr.decode()
    data = json.loads(out.stdout)
    vec = [complex(*a) for a in data['state_vector']]
    return vec, data['statistics']


def agree(path, dd, backend, n):
    vd, sd = run(path, dd, backend)
    vq, sq = run(path, 'qmdd', backend)
    assert len(vd) == len(vq) == 2 ** n
    assert max(abs(a - b) for a, b in zip(vd, vq)) < 1e-9, (path, dd)
    assert abs(sd['norm'] - sq['norm']) < 1e-9
    assert abs(sd['first_qubit_measurement_prob'] - sq['first_qubit_measurement_prob']) < 1e-9
    assert sd['t_count'] == sq['t_count']


@pytest.mark.parametrize("dd", ['bqd', 'limdd'])
def test_clifford_t_exact(dd, tmp_path):
    for n in range(2, 8):
        for seed in range(3):
            path = tmp_path / f"ct_{n}_{seed}.qasm"
            path.write_text(clifford_t_circuit(n, 8 * n, 1000 * n + seed))
            agree(str(path), dd, 'qisq2', n)


@pytest.mark.parametrize("backend", ['qisq2', 'float'])
def test_iqp(backend, tmp_path):
    for n in range(2, 9):
        for seed in range(3):
            path = tmp_path / f"iqp_{n}_{seed}.qasm"
            path.write_text(iqp_circuit(n, 6 * n, 100 * n + seed))
            agree(str(path), 'bqd', backend, n)


@pytest.mark.parametrize("kind, n, gates, log_nodes", [
    ('clifford', 8, 1500, 14),
    ('iqp', 12, 3000, 12),
])
def test_bqd_collects_between_gates(kind, n, gates, log_nodes, tmp_path):
    """A node table small enough that the runner collects between gates, on
    diagrams that skip the qubits a state does not depend on: the sweep, the
    cleared memo and the buckets built on again must leave the EVDD's state.
    Each table is the smallest the runner takes, or at least twice the
    smallest in which its circuit finishes, and the collections are counted,
    so that the case cannot stop collecting unnoticed."""
    path = tmp_path / f"{kind}_{n}.qasm"
    make = clifford_circuit if kind == 'clifford' else iqp_circuit
    path.write_text(make(n, gates, 7))
    small = ['--node-tab-size', str(log_nodes), '--wgt-tab-size', '18']
    vd, sd = run(str(path), 'bqd', 'qisq2', sizes=small)
    vq, sq = run(str(path), 'qmdd', 'qisq2')
    assert sd['limdd_collections'] >= 2, sd['limdd_collections']
    assert len(vd) == len(vq) == 2 ** n
    assert max(abs(a - b) for a, b in zip(vd, vq)) < 1e-9


def test_first_qubit_when_qubit_zero_is_skipped(tmp_path):
    """qreg q[2]; h q[0]: the QMDD skips qubit 0's node, and -m once read
    the probability of the other qubit then."""
    for n in (1, 2, 3):
        path = tmp_path / f"h_{n}.qasm"
        path.write_text("\n".join(HEADER + [f"qreg q[{n}];", "h q[0];"]))
        for dd in ('qmdd', 'limdd', 'bqd'):
            _, st = run(str(path), dd, 'qisq2')
            assert abs(st['first_qubit_measurement_prob'] - 0.5) < 1e-12, (n, dd)


@pytest.mark.parametrize("body", [
    ["h q[0];", "h q[1];", "t q[0];", "h q[1];"],      # a Hadamard on a touched qubit
    ["h q[0];", "cx q[0],q[1];"],                       # not diagonal
])
def test_bqd_float_refuses_cancelling_gates(body, tmp_path):
    path = tmp_path / "bad.qasm"
    path.write_text("\n".join(HEADER + ["qreg q[2];"] + body))
    out = run(str(path), 'bqd', 'float', check=False)
    assert out.returncode != 0
    assert b'bqd:' in out.stderr


def test_rotation_refused_on_exact_weights(tmp_path):
    path = tmp_path / "rz.qasm"
    path.write_text("\n".join(HEADER + ["qreg q[2];", "h q[0];", "rz(0.3) q[0];"]))
    out = run(str(path), 'qmdd', 'qisq2', check=False)
    assert out.returncode == 1
    assert b'angle' in out.stderr
