"""
The LIMDD arm of run_qasm_on_qmdd under -s l2.

With complex weights -s l2 reads every LIMDD node at norm 1 (limdd_set_l2),
which changes what each stored label means but not the state the circuit
computes. So the state vector must be the QMDD arm's, the JSON must say which
normalisation ran (3 for L2, 0 for low), and a LIMDD run without -s keeps the
low normalisation it has always had, although the runner's default strategy
is L2. Exact (qisq2) weights cannot take the square root, so there -s l2 falls
back to low.

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


def circuit(seed, gates=60):
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
    out = subprocess.run([SIM_QASM, str(path), *args, *SIZES, '--state-vector'],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    assert out.returncode == 0, out.stderr.decode()
    data = json.loads(out.stdout)
    return data, [complex(*a) for a in data['state_vector']]


def strategy(data):
    return data.get('wgt_norm_strat', data['statistics'].get('wgt_norm_strat'))


@pytest.mark.skipif(not os.path.exists(SIM_QASM), reason=f"no runner at {SIM_QASM}")
@pytest.mark.parametrize("seed", [1, 2, 3, 4, 5])
def test_limdd_l2_matches_qmdd(seed, tmp_path):
    path = tmp_path / f"c{seed}.qasm"
    path.write_text(circuit(seed))
    _, ref = run(path, ['-d', 'qmdd', '-e', 'float', '-s', 'l2'])

    for args, strat in ((['-d', 'limdd', '-e', 'float', '-s', 'l2'], 3),
                        (['-d', 'limdd', '-e', 'float', '-s', 'low'], 0),
                        (['-d', 'limdd', '-e', 'float'], 0),
                        (['-d', 'limdd', '-e', 'qisq2', '-s', 'l2'], 0)):
        data, vec = run(path, args)
        assert strategy(data) == strat, (args, strategy(data))
        assert abs(data['statistics']['norm'] - 1.0) < 1e-9, (args, data['statistics']['norm'])
        diff = max(abs(a - b) for a, b in zip(vec, ref))
        assert diff < 1e-9, (args, diff)
