"""
The LIMDD arm of run_qasm_on_qmdd under --canon=lazy.

The lazy canonical form keeps one stored node per class of states that are LIM
images of each other, as --canon=always does, but searches for the canonical
label only when another node may share the class. So the state must be the
QMDD arm's, the JSON must say which form ran, and on exact (qisq2) weights the
final diagram must have exactly as many nodes as under --canon=always.

The runner is RUN_QASM_ON_QMDD, or else ./build/qasm/run_qasm_on_qmdd.
"""
import json
import os
import random
import subprocess
import pytest


SIM_QASM = os.environ.get('RUN_QASM_ON_QMDD', './build/qasm/run_qasm_on_qmdd')
SIZES = ['--node-tab-size', '18', '--wgt-tab-size', '18']
HEADER = ["OPENQASM 2.0;", 'include "qelib1.inc";']
N = 7


def circuit(seed, gates=120):
    """Random Clifford+T gates from what the LIMDD arm takes, with H often
    enough that the state is far from a basis state."""
    rng = random.Random(seed)
    lines = HEADER + [f"qreg q[{N}];"]
    for _ in range(gates):
        q = rng.sample(range(N), 3)
        x = rng.random()
        if x < 0.25:
            lines.append(f"h q[{q[0]}];")
        elif x < 0.55:
            lines.append(f"{rng.choice(['t', 'tdg', 's', 'sdg', 'z', 'x', 'y'])} q[{q[0]}];")
        elif x < 0.75:
            lines.append(f"{rng.choice(['cx', 'cz'])} q[{q[0]}],q[{q[1]}];")
        elif x < 0.85:
            lines.append(f"{rng.choice(['cs', 'csdg'])} q[{q[0]}],q[{q[1]}];")
        else:
            lines.append(f"ccz q[{q[0]}],q[{q[1]}],q[{q[2]}];")
    return "\n".join(lines) + "\n"


def run(path, args):
    out = subprocess.run([SIM_QASM, str(path), *args, *SIZES, '-c', '--state-vector'],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    assert out.returncode == 0, out.stderr.decode()
    data = json.loads(out.stdout)
    return data, [complex(*a) for a in data['state_vector']]


@pytest.mark.skipif(not os.path.exists(SIM_QASM), reason=f"no runner at {SIM_QASM}")
@pytest.mark.parametrize("seed", range(1, 13))
def test_limdd_lazy_matches_qmdd(seed, tmp_path):
    path = tmp_path / f"c{seed}.qasm"
    path.write_text(circuit(seed))
    _, ref = run(path, ['-d', 'qmdd', '-e', 'float', '-s', 'l2'])

    for weights in (['-e', 'float', '-s', 'low'], ['-e', 'float', '-s', 'l2'], ['-e', 'qisq2']):
        eager, _ = run(path, ['-d', 'limdd', '--canon=always', *weights])
        lazy, vec = run(path, ['-d', 'limdd', '--canon=lazy', *weights])
        assert lazy.get('canon', lazy['statistics'].get('canon')) == 'lazy'
        assert abs(lazy['statistics']['norm'] - 1.0) < 1e-9, (weights, lazy['statistics']['norm'])
        diff = max(abs(a - b) for a, b in zip(vec, ref))
        assert diff < 1e-9, (weights, diff)
        if weights[1] == 'qisq2':
            assert lazy['statistics']['final_nodes'] == eager['statistics']['final_nodes'], weights
            assert lazy['statistics']['max_nodes'] == eager['statistics']['max_nodes'], weights
