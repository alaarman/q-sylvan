"""
The BQD arm of run_qasm_on_qmdd against the EVDD arm, on the IQP fragment.

A BQD has no amplitude algorithm of its own; -v decodes the whole state, and
that decoded vector must be the EVDD's, on both weight backends. The arm must
also refuse what it has no algorithm for.
"""
import os
import random
import subprocess
import json
import pytest


SIM_QASM = os.environ.get('RUN_QASM_ON_QMDD', './build/qasm/run_qasm_on_qmdd')
SIZES = ['--node-tab-size', '18', '--wgt-tab-size', '18']


def iqp_circuit(n, gates, seed):
    """A Hadamard layer, then level-3 diagonal gates (as gen_iqp.py draws them)."""
    rng = random.Random(seed)
    lines = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{n}];"]
    lines += [f"h q[{q}];" for q in range(n)]
    for _ in range(gates):
        kind = rng.choices(['t', 'c1', 'c2', 'cs', 'ccz'],
                           weights=[3, 1, 1 if n > 1 else 0, 2 if n > 1 else 0, 1 if n > 2 else 0])[0]
        if kind == 't':
            lines.append(f"t q[{rng.randrange(n)}];")
        elif kind == 'c1':
            lines.append(f"{rng.choice(['z', 's', 'sdg'])} q[{rng.randrange(n)}];")
        elif kind == 'c2':
            a, b = rng.sample(range(n), 2)
            lines.append(f"cz q[{a}],q[{b}];")
        elif kind == 'cs':
            a, b = rng.sample(range(n), 2)
            lines.append(f"{rng.choice(['cs', 'csdg'])} q[{a}],q[{b}];")
        else:
            a, b, c = rng.sample(range(n), 3)
            lines.append(f"ccz q[{a}],q[{b}],q[{c}];")
    return "\n".join(lines)


def run(path, dd, backend):
    out = subprocess.run([SIM_QASM, path, '-d', dd, '-e', backend, '-s', 'low',
                          '--state-vector', '-m', *SIZES],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    assert out.returncode == 0, out.stderr.decode()
    data = json.loads(out.stdout)
    vec = [complex(*a) for a in data['state_vector']]
    return vec, data['statistics']


@pytest.mark.parametrize("backend", ['qisq2', 'float'])
def test_bqd_matches_evdd(backend, tmp_path):
    # from two qubits: on one, the EVDD arm reports no first-qubit probability
    for n in range(2, 9):
        for seed in range(3):
            path = tmp_path / f"iqp_{n}_{seed}.qasm"
            path.write_text(iqp_circuit(n, 6 * n, 100 * n + seed))
            vb, sb = run(str(path), 'bqd', backend)
            vq, sq = run(str(path), 'qmdd', backend)
            assert len(vb) == len(vq) == 2 ** n
            assert max(abs(a - b) for a, b in zip(vb, vq)) < 1e-9, (n, seed)
            assert abs(sb['norm'] - sq['norm']) < 1e-9
            assert abs(sb['first_qubit_measurement_prob'] - sq['first_qubit_measurement_prob']) < 1e-9
            assert sb['t_count'] == sq['t_count']


@pytest.mark.parametrize("body", [
    ["h q[0];", "h q[1];", "t q[0];", "h q[1];"],      # a Hadamard after the layer
    ["h q[0];", "h q[1];", "cx q[0],q[1];"],           # not diagonal
    ["h q[0];", "t q[0];", "h q[1];"],                 # a gate before the layer is complete
    ["h q[0];"],                                       # the layer never completes
])
def test_bqd_refuses(body, tmp_path):
    path = tmp_path / "bad.qasm"
    path.write_text("\n".join(["OPENQASM 2.0;", 'include "qelib1.inc";', "qreg q[2];"] + body))
    out = subprocess.run([SIM_QASM, str(path), '-d', 'bqd', *SIZES],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    assert out.returncode != 0
    assert b'bqd:' in out.stderr
