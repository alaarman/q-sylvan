"""
The BQD and LIMDD arms of run_qasm_on_qmdd against the EVDD arm.

Neither a BQD nor a LIMDD prints its state the way a QMDD does; -v decodes it,
and the decoded vector must be the EVDD's, amplitude for amplitude, on the
Clifford+T gate set the two arms share, Hadamards anywhere. The BQD arm runs
in each of its three label families (--bqd-family): the scalar one, the
X-BQD and the Pauli-BQD, each under both zero rules (--bqd-zero): rule SM,
the default, and the paper's copy rule. On float weights it runs every gate
with a warning; there the gates that cannot cancel, the diagonal ones and a
first Hadamard, must still give the EVDD's state, and the others must at
least run.
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


def run(path, dd, backend, check=True, sizes=SIZES, family=None, extra=(), zero=None):
    fam = ['--bqd-family', family] if family else []
    rule = ['--bqd-zero', zero] if zero else []
    out = subprocess.run([SIM_QASM, path, '-d', dd, '-e', backend, '-s', 'low',
                          '--state-vector', '-m', *sizes, *fam, *rule, *extra],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if not check:
        return out
    assert out.returncode == 0, out.stderr.decode()
    data = json.loads(out.stdout)
    vec = [complex(*a) for a in data['state_vector']]
    if family:
        assert data['statistics']['bqd_family'] == family
    if dd == 'bqd':
        # the rule asked for, or the default, SM in every family
        want = zero or 'sm'
        assert data['statistics']['bqd_zero'] == want
    else:
        assert data['statistics']['bqd_zero'] == 'n/a'
    return vec, data['statistics']


def agree(path, dd, backend, n, family=None, zero=None, sizes=SIZES):
    vd, sd = run(path, dd, backend, family=family, zero=zero, sizes=sizes)
    vq, sq = run(path, 'qmdd', backend)
    assert len(vd) == len(vq) == 2 ** n
    assert max(abs(a - b) for a, b in zip(vd, vq)) < 1e-9, (path, dd)
    assert abs(sd['norm'] - sq['norm']) < 1e-9
    assert abs(sd['first_qubit_measurement_prob'] - sq['first_qubit_measurement_prob']) < 1e-9
    assert sd['t_count'] == sq['t_count']


@pytest.mark.parametrize("dd, zero", [('bqd', 'sm'), ('bqd', 'copy'), ('bqd', None),
                                      ('limdd', None)])
def test_clifford_t_exact(dd, zero, tmp_path):
    for n in range(2, 8):
        for seed in range(3):
            path = tmp_path / f"ct_{n}_{seed}.qasm"
            path.write_text(clifford_t_circuit(n, 8 * n, 1000 * n + seed))
            agree(str(path), dd, 'qisq2', n, zero=zero)


@pytest.mark.parametrize("zero", ['sm', 'copy'])
def test_clifford_t_long_exact(zero, tmp_path):
    """Longer Clifford+T circuits on up to 10 qubits, whose states have zeros
    on supports that are not affine after the first T gates meet a Hadamard,
    which is where the two rules make different diagrams; on the third the
    two differ by 22 nodes. The second takes the copy rule 7 s against SM's
    1 s, and a fourth, 160 gates on 10 qubits with the seed 70001, took it 26 s
    and tables of 2^20 against SM's 2 s, so it is not here."""
    for n, gates, seed in ((8, 240, 56000), (8, 240, 56001), (10, 160, 70000)):
        path = tmp_path / f"ctl_{n}_{seed}.qasm"
        path.write_text(clifford_t_circuit(n, gates, seed))
        agree(str(path), 'bqd', 'qisq2', n, zero=zero)


@pytest.mark.parametrize("zero", ['sm', 'copy', None])
@pytest.mark.parametrize("family", ['x', 'pauli'])
def test_bqd_families_clifford_t_exact(family, zero, tmp_path):
    """The X-BQD and the Pauli-BQD on the circuits of the scalar arm above,
    under each zero rule and the default."""
    for n in range(2, 8):
        for seed in range(3):
            path = tmp_path / f"ct_{n}_{seed}.qasm"
            path.write_text(clifford_t_circuit(n, 8 * n, 1000 * n + seed))
            agree(str(path), 'bqd', 'qisq2', n, family, zero=zero)


@pytest.mark.parametrize("family", ['x', 'pauli'])
def test_bqd_families_clifford_t_long_exact(family, tmp_path):
    """The longer circuits of test_clifford_t_long_exact in the X-BQD and the
    Pauli-BQD under rule SM, where S nodes with a translation on their high
    edge, and the change of tag a translation can make, meet every gate."""
    for n, gates, seed in ((8, 240, 56000), (8, 240, 56001), (10, 160, 70000)):
        path = tmp_path / f"ctl_{n}_{seed}.qasm"
        path.write_text(clifford_t_circuit(n, gates, seed))
        agree(str(path), 'bqd', 'qisq2', n, family, zero='sm')


@pytest.mark.parametrize("zero", ['sm', 'copy'])
@pytest.mark.parametrize("backend", ['qisq2', 'float'])
def test_iqp(backend, zero, tmp_path):
    for n in range(2, 9):
        for seed in range(3):
            path = tmp_path / f"iqp_{n}_{seed}.qasm"
            path.write_text(iqp_circuit(n, 6 * n, 100 * n + seed))
            agree(str(path), 'bqd', backend, n, zero=zero)


@pytest.mark.parametrize("zero", ['sm', 'copy'])
@pytest.mark.parametrize("backend", ['qisq2', 'float'])
@pytest.mark.parametrize("family", ['x', 'pauli'])
def test_bqd_families_iqp(family, backend, zero, tmp_path):
    for n in range(2, 9):
        for seed in range(3):
            path = tmp_path / f"iqp_{n}_{seed}.qasm"
            path.write_text(iqp_circuit(n, 6 * n, 100 * n + seed))
            agree(str(path), 'bqd', backend, n, family, zero=zero)


@pytest.mark.parametrize("workers", [1, 4])
@pytest.mark.parametrize("kind, n, gates, log_nodes, family, zero", [
    ('clifford', 8, 4500, 14, None, 'sm'),
    ('iqp', 12, 3000, 12, None, 'sm'),
    ('clifford_t', 8, 400, 14, None, 'sm'),
    ('clifford', 8, 4500, 14, None, 'copy'),
    ('iqp', 12, 3000, 12, None, 'copy'),
    ('clifford_t', 7, 400, 15, None, 'copy'),
    ('clifford', 8, 4500, 14, 'x', 'sm'),
    ('iqp', 12, 3000, 12, 'x', 'sm'),
    ('clifford_t', 8, 400, 15, 'x', 'sm'),
    ('clifford', 8, 4500, 14, 'x', 'copy'),
    ('iqp', 12, 3000, 12, 'x', 'copy'),
    ('clifford', 8, 9000, 13, 'pauli', 'sm'),
    ('iqp', 12, 3000, 12, 'pauli', 'sm'),
    ('clifford_t', 8, 400, 15, 'pauli', 'sm'),
    ('clifford', 8, 9000, 13, 'pauli', 'copy'),
    ('iqp', 12, 3000, 12, 'pauli', 'copy'),
])
def test_bqd_collects_between_gates(kind, n, gates, log_nodes, family, zero, workers, tmp_path):
    """A node table small enough that the runner collects between gates, on
    diagrams that skip the qubits a state does not depend on: the sweep, the
    cleared memo and the buckets built on again must leave the EVDD's state,
    in every family, on one worker and on four. Each table is the smallest the
    runner takes, or at least twice the smallest in which its circuit
    finishes, and the collections are counted, so that the case cannot stop
    collecting unnoticed. A swap, a CX with the control below and a diagonal
    gate make no high cofactor above their qubits, so the Clifford circuit
    leaves little garbage per gate and runs 4500 gates to collect several
    times; the Pauli-BQD keeps so few nodes on it that it runs twice as many.
    Its smallest table is twice the one in which one worker finishes, since
    four workers each hold a region of 512 buckets and fill a table of 2^12
    early. The Clifford+T circuit has states with zeros on supports that are
    not affine, where rule SM makes S nodes and the copy rule copies; the copy
    rule makes so many more nodes within one of its gates that it runs on a
    qubit fewer, in a table twice as large. The X-BQD and the Pauli-BQD run it
    under SM in 2^15, the smallest table in which the Pauli-BQD finishes and
    twice the X-BQD's; one size more and either collects once. Under the copy
    rule those two take five to seven times as long on it, and are left to
    the scalar family."""
    path = tmp_path / f"{kind}_{n}.qasm"
    make = {'clifford': clifford_circuit, 'iqp': iqp_circuit,
            'clifford_t': clifford_t_circuit}[kind]
    path.write_text(make(n, gates, 7))
    small = ['--node-tab-size', str(log_nodes), '--wgt-tab-size', '18']
    vd, sd = run(str(path), 'bqd', 'qisq2', sizes=small, family=family,
                 extra=['-w', str(workers)], zero=zero)
    vq, sq = run(str(path), 'qmdd', 'qisq2')
    assert sd['limdd_collections'] >= 2, sd['limdd_collections']
    assert len(vd) == len(vq) == 2 ** n
    assert max(abs(a - b) for a, b in zip(vd, vq)) < 1e-9


@pytest.mark.parametrize("log_nodes", [12, 13])
def test_bqd_pauli_words_do_not_fill_first(log_nodes, tmp_path):
    """The Pauli-BQD makes a new Pauli word with many of its label products,
    far more of them than new nodes. With the Pauli-word table the size of
    the node table, as the runner once gave it, this circuit stopped with
    that table full while the node and LIM tables were nearly empty, and no
    collection was asked for. The runner now gives it the LIM table's size,
    and the collection trigger watches it."""
    path = tmp_path / "cl_14.qasm"
    path.write_text(clifford_circuit(14, 350, 7))
    sizes = ['--node-tab-size', str(log_nodes), '--lim-tab-size', '20', '--wgt-tab-size', '18']
    vd, sd = run(str(path), 'bqd', 'qisq2', sizes=sizes, family='pauli')
    vq, sq = run(str(path), 'qmdd', 'qisq2')
    assert len(vd) == len(vq) == 2 ** 14
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


@pytest.mark.parametrize("family, zero", [(None, 'sm'), (None, 'copy'), ('x', 'sm'), ('x', 'copy'),
                                          ('pauli', 'sm'), ('pauli', 'copy')])
def test_bqd_float_runs_with_a_warning(family, zero, tmp_path):
    """Float weights are taken for every gate in every family, with a warning
    on stderr, and gates that can cancel must not crash the run."""
    for n in range(2, 7):
        for seed in range(3):
            path = tmp_path / f"ct_{n}_{seed}.qasm"
            path.write_text(clifford_t_circuit(n, 8 * n, 1000 * n + seed))
            out = run(str(path), 'bqd', 'float', check=False, family=family, zero=zero)
            assert out.returncode == 0, (path, out.stderr)
            assert b'bqd: warning: float weights' in out.stderr


def test_bqd_family_needs_bqd(tmp_path):
    path = tmp_path / "h.qasm"
    path.write_text("\n".join(HEADER + ["qreg q[1];", "h q[0];"]))
    out = run(str(path), 'limdd', 'qisq2', check=False, family='x')
    assert out.returncode != 0
    assert b'bqd-family' in out.stderr


def test_bqd_zero_options(tmp_path):
    """--bqd-zero needs -d bqd, knows two rules, and takes either in every
    family; the JSON names the rule the run had, its default, SM, included."""
    path = tmp_path / "h.qasm"
    path.write_text("\n".join(HEADER + ["qreg q[2];", "h q[0];", "cx q[0],q[1];"]))
    out = run(str(path), 'limdd', 'qisq2', check=False, zero='copy')
    assert out.returncode != 0 and b'bqd-zero' in out.stderr
    out = run(str(path), 'bqd', 'qisq2', check=False, zero='shannon')
    assert out.returncode != 0 and b'zero rule' in out.stderr
    for family in ('x', 'pauli'):
        for zero in (None, 'sm', 'copy'):
            run(str(path), 'bqd', 'qisq2', family=family, zero=zero)
    for zero in (None, 'sm', 'copy'):
        run(str(path), 'bqd', 'qisq2', zero=zero)
        run(str(path), 'bqd', 'qisq2', family='scalar', zero=zero)
    _, st = run(str(path), 'qmdd', 'qisq2')
    assert st['bqd_zero'] == 'n/a'


def test_rotation_refused_on_exact_weights(tmp_path):
    path = tmp_path / "rz.qasm"
    path.write_text("\n".join(HEADER + ["qreg q[2];", "h q[0];", "rz(0.3) q[0];"]))
    out = run(str(path), 'qmdd', 'qisq2', check=False)
    assert out.returncode == 1
    assert b'angle' in out.stderr
