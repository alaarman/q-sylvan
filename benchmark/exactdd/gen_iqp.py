#!/usr/bin/env python3
"""Random IQP circuits in the level-3 diagonal gate set, as OPENQASM 2.0, dependency-free.

An IQP circuit is a layer of Hadamards, a diagonal unitary, and a layer of
Hadamards before measurement.  What is emitted here is the first two: the
Hadamard on every qubit, then the diagonal gates, so that every intermediate
state is a phase state.  The closing layer is left out on purpose: the binary
quotient diagram has no Hadamard algorithm, and the three arms (EVDD, LIMDD,
BQD) are compared on the same file, so the file stops where all three can go.

The diagonal gates are the third level of the Clifford hierarchy on at most
three qubits, which all three arms of run_qasm_on_qmdd apply as one gate and
exactly on the qisq2 backend:

    t                       one qubit, level 3
    z, s, sdg, cz           diagonal Clifford
    cs, csdg                two qubits, level 3
    ccz                     three qubits, level 3

so the phase polynomial is a level-3 one over Z_8, and the BQD's per-level
bound is sum_{i<3} C(v, i) + 1 (thm:size).  The cs and ccz gates matter: with
t and cz alone the high cofactor of every level is a Pauli times the low one,
the LIMDD keeps width 1, and the file tells nothing.  cp with an angle is not
used, since it is not exact on qisq2.

Four classes, drawn with weights --t-frac (default 1/3), --cs-frac (2/9),
--ccz-frac (1/9) and the diagonal Cliffords for the rest.  Operands are an
ordered tuple of distinct qubits drawn uniformly.  The count in the file name
is the number of DIAGONAL gates; the Hadamard layer adds n more.  The T-count
the runner reports counts cs as 3 and ccz as 7, their ancilla-free Clifford+T
decompositions.
"""
import argparse, os, random

CLIFF_1 = ['z', 's', 'sdg']
CLIFF_2 = ['cz']
CS = ['cs', 'csdg']


def _qubits(rng, n, k):
    return rng.sample(range(n), k)


def _draw(rng, n, t_frac, cs_frac, ccz_frac):
    kinds = ['t', 'cliff', 'cs', 'ccz']
    w = [t_frac, max(0.0, 1.0 - t_frac - cs_frac - ccz_frac), cs_frac, ccz_frac]
    if n < 3:
        w[3] = 0.0
    if n < 2:
        w[2] = 0.0
    kind = rng.choices(kinds, weights=w, k=1)[0]
    if kind == 't':
        return f"t q[{_qubits(rng, n, 1)[0]}];"
    if kind == 'cs':
        a, b = _qubits(rng, n, 2)
        return f"{rng.choice(CS)} q[{a}],q[{b}];"
    if kind == 'ccz':
        a, b, c = _qubits(rng, n, 3)
        return f"ccz q[{a}],q[{b}],q[{c}];"
    if n >= 2 and rng.random() < 0.5:
        a, b = _qubits(rng, n, 2)
        return f"{rng.choice(CLIFF_2)} q[{a}],q[{b}];"
    return f"{rng.choice(CLIFF_1)} q[{_qubits(rng, n, 1)[0]}];"


def generate(n, gates, seed, t_frac=1.0 / 3.0, cs_frac=2.0 / 9.0, ccz_frac=1.0 / 9.0):
    rng = random.Random(seed)
    head = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{n}];"]
    body = [f"h q[{q}];" for q in range(n)]
    body += [_draw(rng, n, t_frac, cs_frac, ccz_frac) for _ in range(gates)]
    # no creg, no measure, no trailing newline: the artifact's convention
    return "\n".join(head + body)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--qubits", type=int, nargs="+", required=True)
    ap.add_argument("--gates", type=int, nargs="+", required=True,
                    help="number of diagonal gates after the Hadamard layer")
    ap.add_argument("--t-frac", type=float, default=1.0 / 3.0)
    ap.add_argument("--cs-frac", type=float, default=2.0 / 9.0)
    ap.add_argument("--ccz-frac", type=float, default=1.0 / 9.0)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--name", default="iqp")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    print(f"{'file':32s}{'gates':>7}{'t':>6}{'cs':>6}{'ccz':>6}")
    for n in a.qubits:
        for g in a.gates:
            txt = generate(n, g, a.seed + 1000 * n + g, a.t_frac, a.cs_frac, a.ccz_frac)
            f = os.path.join(a.out, f"{a.name}_{n}_{g}.qasm")
            with open(f, "w") as fh:
                fh.write(txt)
            cnt = lambda p: sum(1 for l in txt.split("\n") if l.startswith(p))
            print(f"{os.path.basename(f):32s}{g:>7}{cnt('t '):>6}{cnt('cs'):>6}{cnt('ccz'):>6}")
