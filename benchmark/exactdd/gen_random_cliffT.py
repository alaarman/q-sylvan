#!/usr/bin/env python3
"""Random Clifford+T circuits as OPENQASM 2.0, dependency-free.

Reproduces the gate distribution of the surviving paper artifact
(git ref quist/qisq2_tests, clifford_T_circuit_<qubits>_<gates>.qasm), which is
NOT what benchmark/generate_qasm_random_clif_T_gate.py in the q-sylvan tree
emits today.  Measured over all 49 artifact files (41650 gates):

    t     33.81%     cz  11.18%   sdg 5.67%  y 5.57%  x 5.46%
                     swap 10.91%  h   5.62%  s 5.50%  z 5.43%
                     cx  10.85%

i.e. three equiprobable classes -- T / single-qubit Clifford / two-qubit
Clifford -- with the gate drawn uniformly inside its class:

    1q classes : h s sdg x y z        (six; NO sx/sxdg)
    2q classes : cx cz swap           (three; NO cy)

The in-tree generator differs on all three counts: it has sx/sxdg and cy, it
has swap commented out, and its globals set the T fraction to 85%.  sx/sxdg
are also refused by the LIMDD path, so the in-tree script cannot produce a
runnable LIMDD benchmark at all.

Two-qubit operand order: the artifact draws an ordered pair (c,t), c != t,
uniformly -- 6803 of its 13721 two-qubit gates have control < target and 6918
have control > target, i.e. 50/50.  The in-tree generator instead rejects
until target > control.  Default here is "random" to match the artifact;
--ctrl-order ascending reproduces the in-tree convention.

Two sizing modes:
  --gates N      emit exactly N gates          (the artifact's convention; the
                 number in the artifact's filename is a GATE COUNT, not depth)
  --depth D      emit gates until the circuit's ASAP depth is exactly D
                 (what the paper's prose literally says)

Depth is the standard DAG longest path: level[q] for each wire, a gate lands at
1 + max(level of its operands).  No estimate of gates-per-layer is used -- the
--depth mode grows the circuit and stops the moment the depth hits the target.
"""
import argparse, os, random

ONE_QUBIT = ['h', 's', 'sdg', 'x', 'y', 'z']
TWO_QUBIT = ['cx', 'cz', 'swap']

CLASSES = ['t', 'one', 'two']


def _draw(rng, n, t_frac, ctrl_order):
    """One gate: returns (qasm line, qubit tuple)."""
    p_rest = 0.5 * (1.0 - t_frac)
    kind = rng.choices(CLASSES, weights=[t_frac, p_rest, p_rest], k=1)[0]
    if kind == 'two' and n > 1:
        c = t = 0
        if ctrl_order == 'ascending':
            while t <= c:
                c = rng.randint(0, n - 1)
                t = rng.randint(0, n - 1)
        else:
            while t == c:
                c = rng.randint(0, n - 1)
                t = rng.randint(0, n - 1)
        return f"{rng.choice(TWO_QUBIT)} q[{c}],q[{t}];", (c, t)
    if kind == 't':
        q = rng.randint(0, n - 1)
        return f"t q[{q}];", (q,)
    q = rng.randint(0, n - 1)
    return f"{rng.choice(ONE_QUBIT)} q[{q}];", (q,)


def generate(n, t_frac, seed, ctrl_order='random', gates=None, depth=None):
    """Returns (qasm text, gate count, measured depth)."""
    assert (gates is None) != (depth is None), "give exactly one of gates/depth"
    rng = random.Random(seed)
    body, level, d, count = [], [0] * n, 0, 0
    while True:
        if gates is not None and count >= gates:
            break
        if depth is not None and d >= depth:
            break
        line, qs = _draw(rng, n, t_frac, ctrl_order)
        lvl = max(level[q] for q in qs) + 1
        if depth is not None and lvl > depth:
            continue          # would overshoot the target depth; redraw
        for q in qs:
            level[q] = lvl
        d = max(d, lvl)
        body.append(line)
        count += 1
    head = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{n}];"]
    # No creg / no measure: the artifact has none either.  The paper's
    # "single-qubit measurement on the top qubit" is done by the runner (-m).
    # No trailing newline: all 49 artifact files end on the last ';'.
    return "\n".join(head + body), count, d


def measure_depth(path):
    """Independent re-measurement of a written file (used for verification)."""
    import re
    n, level, d = None, None, 0
    for ln in open(path):
        ln = ln.strip()
        if ln.startswith('qreg'):
            n = int(re.search(r'\[(\d+)\]', ln).group(1))
            level = [0] * n
            continue
        if not ln or ln.startswith(('OPENQASM', 'include', 'creg')):
            continue
        qs = [int(x) for x in re.findall(r'q\[(\d+)\]', ln)]
        if not qs:
            continue
        lvl = max(level[q] for q in qs) + 1
        for q in qs:
            level[q] = lvl
        d = max(d, lvl)
    return d


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--qubits", type=int, nargs="+", required=True)
    ap.add_argument("--gates", type=int, nargs="+",
                    help="gate counts (artifact convention)")
    ap.add_argument("--depth", type=int, nargs="+",
                    help="target ASAP depths (paper prose)")
    ap.add_argument("--t-frac", type=float, default=1.0 / 3.0,
                    help="fraction of T gates (default 1/3, artifact measures 0.3381)")
    ap.add_argument("--ctrl-order", choices=['random', 'ascending'], default='random')
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--name", default="clifford_T_circuit",
                    help="filename stem; file is <stem>_<qubits>_<size>.qasm")
    a = ap.parse_args()
    if (a.gates is None) == (a.depth is None):
        ap.error("give exactly one of --gates / --depth")
    os.makedirs(a.out, exist_ok=True)
    sizes = a.gates if a.gates is not None else a.depth
    print(f"{'file':44s}{'gates':>7}{'depth':>7}{'T%':>7}")
    for n in a.qubits:
        for s in sizes:
            seed = a.seed + 1000 * n + s
            kw = {'gates': s} if a.gates is not None else {'depth': s}
            txt, cnt, d = generate(n, a.t_frac, seed, a.ctrl_order, **kw)
            f = os.path.join(a.out, f"{a.name}_{n}_{s}.qasm")
            with open(f, "w") as fh:
                fh.write(txt)
            tn = sum(1 for l in txt.split("\n") if l.startswith("t "))
            assert measure_depth(f) == d, "depth re-measure disagrees"
            print(f"{os.path.basename(f):44s}{cnt:>7}{d:>7}{100.0*tn/cnt:>7.2f}")
