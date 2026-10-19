#!/usr/bin/env python3
"""Mirror (echo) circuits over random Clifford+T, as OPENQASM 2.0, dependency-free.

A mirror circuit runs a circuit U forward, a middle layer M, and then U
backwards:

    [H on every qubit]  U  M  V            V = U^-1, or a quasi-inverse (pauli)

U is a random Clifford+T circuit with the sampler of gen_random_cliffT.py
(imported, not copied; same --seed convention, seed + 1000*n + g), so the
first half of mirror-<mode>_<n>_<g> is gate for gate the file
clifford_T_circuit_<n>_<g>.qasm that gen_random_cliffT.py writes with the same
--seed.  --base instead mirrors existing QASM files, e.g. the paper's
artifact (benchmark/exactdd/README.md, orig/qasm/).  V reverses U and inverts
each gate: t<->tdg, s<->sdg, cs<->csdg, the rest are self-inverse.

Middle layers (--middle):

  none   M = I.  V U = I, so the final state is |0...0> (|+...+> with --plus):
         the identity miter of DD-based equivalence checking.
  pauli  M = a uniformly random Pauli string P, and V is the quasi-inverse of
         Proctor et al., "Scalable randomized benchmarking of quantum
         computers using mirror circuits", PRL 129 (2022): P is carried
         through U^-1 gate by gate as a Pauli frame Q.  A Clifford gate is
         inverted as usual and conjugates the frame; a T (or Tdg) gate whose
         qubit carries X or Y in the frame is NOT inverted, since
         T Q T = e^{i pi/4} Q there, and Q passes through unchanged.  Then
         V P U is a Pauli, and the final state is a computational basis state
         up to phase (a product of |+>/|-> with --plus).  With plain U^-1 the
         T gates would leave U^-1 P U non-Pauli and the state nontrivial.
  x0     M = X on one qubit (--x-qubit, default the last, q[n-1], the bottom
         variable, far from q[0] that the runner's -m measures), V = U^-1:
         a Loschmidt echo with a butterfly operator, as in OTOC experiments.
         U^-1 X U is a scrambled operator, the final state is nontrivial.

Why this family.  In none and pauli modes the answer is a product state known
in advance (n+1 EVDD nodes with the terminal for |0...0>, 1 for |+...+>),
while the state in the middle is the random circuit's.  In none mode an exact
EVDD walks back through the very diagrams it built on the way out (the state
after the k-th gate of V is the state before the k-th last gate of U), so its
second half costs what the first did and it ends on the product state.  A
float EVDD unwinds approximations: rounding residues that should cancel keep
nodes apart, and its final size, norm and top-qubit probability show how far
it drifted.  x0 is the same echo with a nontrivial (from |0>: sparse) answer.
--plus makes the middle state dense with many distinct amplitudes, which is
where float stops merging and exact pays for long coefficients.  Measured on
mirror-none-plus_12_300 (-w 1, run_arm.sh): qisq2_low 6.5 s, 1 node at the
end, norm 1; float_l2 1.3 s, 3,549 nodes (of 4,096 possible); float_low
0.18 s, norm 0.0.

The count in the file name is the number of gates in U (one half); the file
has 2g gates, plus the middle layer and n Hadamards with --plus.  No creg, no
measure, no trailing newline: the artifact's convention.  Gates emitted:
those of U (h s sdg x y z t cx cz swap from the sampler; t tdg s sdg cs csdg
h x y z id cx cy cz ch swap ccz from a --base file), their inverses, and
x y z in the middle -- all inside the set the LIMDD path accepts.
"""
import argparse, os, random, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_random_cliffT as rct  # noqa: E402  (same directory)

INVERSE = {'t': 'tdg', 'tdg': 't', 's': 'sdg', 'sdg': 's', 'cs': 'csdg', 'csdg': 'cs',
           'h': 'h', 'x': 'x', 'y': 'y', 'z': 'z', 'id': 'id',
           'cx': 'cx', 'cy': 'cy', 'cz': 'cz', 'ch': 'ch', 'swap': 'swap', 'ccz': 'ccz'}

_GATE = re.compile(r'^([a-z]+)\s+(q\[\d+\](?:\s*,\s*q\[\d+\])*)\s*;$')


def parse(text):
    """(n, [(gate, (qubits...)), ...]) of a single-qreg 'q' QASM file."""
    n, ops = None, []
    for ln in text.split("\n"):
        ln = ln.strip()
        if not ln or ln.startswith(('OPENQASM', 'include', '//')):
            continue
        if ln.startswith('qreg'):
            if n is not None:
                raise ValueError("more than one qreg")
            m = re.fullmatch(r'qreg\s+q\[(\d+)\]\s*;', ln)
            if not m:
                raise ValueError(f"qreg must be q[n]: {ln}")
            n = int(m.group(1))
            continue
        m = _GATE.match(ln)
        if not m or m.group(1) not in INVERSE:
            raise ValueError(f"gate not mirrored: {ln}")
        ops.append((m.group(1), tuple(int(x) for x in re.findall(r'q\[(\d+)\]', m.group(2)))))
    if n is None:
        raise ValueError("no qreg")
    return n, ops


def random_u(n, gates, seed, t_frac, ctrl_order):
    """U exactly as gen_random_cliffT.generate draws it for these arguments."""
    txt, _, _ = rct.generate(n, t_frac, seed, ctrl_order, gates=gates)
    return parse(txt)[1]


def _conjugate(g, qs, x, z):
    """Pauli frame (x, z bits, phase dropped) through a Clifford gate g.
    g and g^-1 act alike on the frame for every gate here."""
    if g in ('x', 'y', 'z', 'id'):
        return
    if g == 'h':
        q, = qs
        x[q], z[q] = z[q], x[q]
    elif g in ('s', 'sdg'):
        q, = qs
        z[q] ^= x[q]
    elif g == 'cx':
        c, t = qs
        x[t] ^= x[c]
        z[c] ^= z[t]
    elif g == 'cz':
        a, b = qs
        z[a] ^= x[b]
        z[b] ^= x[a]
    elif g == 'swap':
        a, b = qs
        x[a], x[b] = x[b], x[a]
        z[a], z[b] = z[b], z[a]
    else:
        raise ValueError(f"--middle pauli: no Pauli-frame rule for {g}")


def mirror(n, u, middle, plus=False, mid_seed=0, x_qubit=-1):
    """Returns (H layer + U lines, M lines, V lines, expected final state).

    expected: a string over 0/1 (basis state, q[0] first), over +/- with
    --plus, or None where the state is nontrivial (x0)."""
    pre = [f"h q[{q}];" for q in range(n)] if plus else []
    if middle == 'none':
        mid, back = [], [(INVERSE[g], qs) for g, qs in reversed(u)]
        exp = ('+' if plus else '0') * n
    elif middle == 'x0':
        q = x_qubit % n
        mid, back = [('x', (q,))], [(INVERSE[g], qs) for g, qs in reversed(u)]
        exp = None
    elif middle == 'pauli':
        rng = random.Random(mid_seed)
        p = [rng.choice('IXYZ') for _ in range(n)]
        mid = [(c.lower(), (q,)) for q, c in enumerate(p) if c != 'I']
        x = [int(c in 'XY') for c in p]
        z = [int(c in 'YZ') for c in p]
        back = []
        for g, qs in reversed(u):
            if g in ('t', 'tdg'):
                # T Q T = phase * Q when Q has X or Y on T's qubit: keep T.
                back.append((g if x[qs[0]] else INVERSE[g], qs))
            else:
                back.append((INVERSE[g], qs))
                _conjugate(g, qs, x, z)
        # final = Q |0..0> = |x>, or Q |+..+> = (+/-) per qubit by z
        exp = ''.join('-+'[1 - b] for b in z) if plus else ''.join(map(str, x))
    else:
        raise ValueError(middle)
    fmt = lambda ops: [f"{g} " + ",".join(f"q[{q}]" for q in qs) + ";" for g, qs in ops]
    return pre + fmt(u), fmt(mid), fmt(back), exp


def qasm(n, pre_fwd, mid, back):
    head = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{n}];"]
    return "\n".join(head + pre_fwd + mid + back)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--middle", choices=['none', 'pauli', 'x0'], nargs="+", default=['none'],
                    help="middle layer(s); one file per mode")
    ap.add_argument("--plus", action="store_true",
                    help="prepend H on every qubit, so U acts on |+...+>")
    ap.add_argument("--x-qubit", type=int, default=-1,
                    help="qubit of the X in --middle x0 (negative counts from the end; default -1)")
    ap.add_argument("--qubits", type=int, nargs="+", help="random U: qubit counts")
    ap.add_argument("--gates", type=int, nargs="+", help="random U: gates per half")
    ap.add_argument("--pairs", nargs="+", metavar="N:G",
                    help="random U: explicit (qubits, gates) pairs instead of the grid")
    ap.add_argument("--base", nargs="+", metavar="FILE",
                    help="mirror these QASM files instead of a random U; the file "
                         "stem replaces <n>_<g> in the name")
    ap.add_argument("--t-frac", type=float, default=1.0 / 3.0)
    ap.add_argument("--ctrl-order", choices=['random', 'ascending'], default='random')
    ap.add_argument("--seed", type=int, default=1234,
                    help="as gen_random_cliffT.py: U uses seed + 1000*n + g")
    ap.add_argument("--name", default="mirror")
    a = ap.parse_args()

    jobs = []   # (name suffix, n, U, seed)
    if a.base:
        if a.qubits or a.gates or a.pairs:
            ap.error("--base excludes --qubits/--gates/--pairs")
        for f in a.base:
            n, u = parse(open(f).read())
            stem = os.path.splitext(os.path.basename(f))[0]
            jobs.append((stem, n, u, a.seed + sum(map(ord, stem))))
    else:
        pairs = [tuple(map(int, p.split(':'))) for p in a.pairs] if a.pairs else \
                [(n, g) for n in (a.qubits or []) for g in (a.gates or [])]
        if not pairs:
            ap.error("give --qubits and --gates, --pairs, or --base")
        for n, g in pairs:
            s = a.seed + 1000 * n + g
            jobs.append((f"{n}_{g}", n, random_u(n, g, s, a.t_frac, a.ctrl_order), s))

    os.makedirs(a.out, exist_ok=True)
    print(f"{'file':46s}{'n':>4}{'gates':>7}{'T':>6}  final state")
    for mode in a.middle:
        tag = f"{a.name}-{mode}" + ("-plus" if a.plus else "")
        for suffix, n, u, s in jobs:
            pre_fwd, mid, back, exp = mirror(n, u, mode, a.plus, s ^ 0x5EED, a.x_qubit)
            txt = qasm(n, pre_fwd, mid, back)
            f = os.path.join(a.out, f"{tag}_{suffix}.qasm")
            with open(f, "w") as fh:
                fh.write(txt)
            body = pre_fwd + mid + back
            tc = sum(1 for l in body if l.startswith(('t ', 'tdg ')))
            # the runner counts cs/csdg as 3 T and ccz as 7; only --base can bring them
            print(f"{os.path.basename(f):46s}{n:>4}{len(body):>7}{tc:>6}  "
                  f"{'|' + exp + '>' if exp else 'nontrivial'}")
