#!/usr/bin/env python3
"""Rotation algorithms compiled to Clifford+T by Qiskit's synthesis, as OPENQASM 2.0.

NEEDS QISKIT (>= 2.2, tested with 2.5.2) AND NUMPY, unlike the other gen_*.py
here.  Run it with the venv that has them, e.g. venv-qiskit/bin/python made
from amendment/requirements-qiskit.txt (amendment/REPRODUCE.md, section 8).

The family.  A textbook algorithm with arbitrary-angle rotations is built in
Qiskit, transpiled (optimization_level 1, fixed seed) to the basis
{h, s, sdg, x, y, z, cx, cz, swap, rz}, and every rz(theta) is then replaced
by a single-qubit Clifford+T word:

  --synth gridsynth --eps E   Ross-Selinger (qiskit.synthesis.gridsynth_rz),
                              operator-norm error <= E, about 3 log2(1/E) T
                              gates and 2.5x that many gates per rotation;
  --synth sk --depth d        Solovay-Kitaev (SolovayKitaevDecomposition,
                              basis h,t,tdg, basic-approximation depth
                              --sk-basic-depth), recursion degree d; roughly
                              5x more gates per extra degree, and much worse
                              error per gate than gridsynth.

An rz whose angle is a multiple of pi/4 is not approximated: it is t^k up to a
global phase and is emitted exactly.  Global phases are dropped (no rz is ever
controlled).  The words are joined and a peephole pass merges runs of diagonal
single-qubit gates on one wire into the shortest of t, s, s t, z, z t, sdg, tdg
(exponent mod 8) and cancels adjacent h h, x x, y y; it never touches a
two-qubit gate.  Only h s sdg t tdg x y z cx cz swap are emitted, so the EVDD
and LIMDD paths of run_qasm_on_qmdd both accept the files and qisq2 represents
every gate exactly.

Algorithms (--alg), qubit q[0] is the top variable of the diagram:

  qft         X gates writing a seeded basis state |x>, then QFT_n (Qiskit's
              QFTGate, with the final swaps).  Exactly, QFT|x> is a product
              state; the synthesised words are slightly off-diagonal, so the
              compiled state carries O(eps) entangled residue on every branch.
  qft_mirror  |x>, the synthesised QFT word W, then the gate-wise inverse of the
              SAME W (reversed, t<->tdg, s<->sdg).  W^-1 W = 1 exactly, so the
              final state is exactly |x>: a chain of n nodes for an exact
              backend, whatever eps.  A float backend has to cancel tens of
              thousands of rounded products to get there.
  tfim        transverse-field Ising chain H = -J sum Z_i Z_{i+1} - h sum X_i
              (open boundary), from |0..0>, --steps first-order Trotter steps
              of --dt: per step rzz(-2 J dt) on every bond, then rx(-2 h dt)
              on every qubit.  All rotations have the same generic angle, so
              the circuit repeats one Clifford+T word.
  qaoa        MaxCut QAOA on a seeded random 3-regular graph (pairing model),
              |+..+>, --layers p of exp(-i gamma C) (rzz(-gamma) per edge) and
              exp(-i beta sum X) (rx(2 beta) per qubit), fixed generic angles.

Why it is here: exact Clifford+T arithmetic is exact on every gate of these
files, so whatever cancellation the algorithm has (the mirror is the extreme
case) is performed exactly, while the float backend sees the same long words
and accumulates rounding.  The cost for exact arithmetic is circuit length: a
gridsynth word at eps 1e-6 is about 160 gates with 60 T gates, so coefficient
sizes in Z[i, 1/sqrt2] grow with every rotation that does not cancel.

First measurements (validation runs, 2026-10-05, EVDD -s low / -s l2, w1):
the exact arm is 50-300x slower than float on every family here, because the
approximate words make the state dense (2^n nodes for qft/tfim/qaoa at small
n) and the Z[i,1/sqrt2] coefficients long (qft n=6 eps 1e-3: 18.5 s exact,
0.06 s float).  On qft_mirror only the exact arm ends at the n+1 node chain:
float -s low loses the norm (0 at n=4,5,8,12), float -s l2 keeps the
probability but ends with 13/25/183/1067 nodes at n=4/5/8/12 against 5/6/9/13,
and the float LIMDD ends with norm 4e-7 (n=3), 1.54 (n=4), 10.4 (n=8).

Each instance gets a sidecar <name>.json: the algorithm and its parameters,
the synthesis settings, gate count, T count (t + tdg), per-rotation statistics
(gates, T, worst operator-norm error of one word), and for n <= --fid-max the
fidelity |<psi_exact|psi_synth>|^2 against the unapproximated circuit (for
qft_mirror the exact state is |x>), with both states computed by Qiskit's
Statevector, the synthesised one from the written QASM text.

File names: synth_<alg>_n<n>[_r<steps>][_p<layers>]_<tag>.qasm, with tag e3
for --eps 1e-3 (e<k> for 1e-k, else eps<value>) and sk<d> for Solovay-Kitaev.
As the other generators: no creg, no measure, no trailing newline.
"""
import argparse, json, math, os, random, sys, time
from collections import Counter

try:
    import numpy as np
    from qiskit import QuantumCircuit, transpile
    from qiskit.circuit.library import QFTGate, RZGate
    from qiskit.quantum_info import Operator, Statevector
    from qiskit.synthesis import SolovayKitaevDecomposition, gridsynth_rz
except ImportError as e:              # pragma: no cover
    sys.exit(f"gen_synth.py needs qiskit >= 2.2 and numpy ({e}); run it with the "
             "venv-qiskit python, see the module docstring")

BASIS = ['h', 's', 'sdg', 'x', 'y', 'z', 'cx', 'cz', 'swap', 'rz']
EMIT = {'h', 's', 'sdg', 't', 'tdg', 'x', 'y', 'z', 'cx', 'cz', 'swap'}
DIAG = {'t': 1, 's': 2, 'z': 4, 'sdg': 6, 'tdg': 7}          # T exponent mod 8
DIAG_WORD = {0: [], 1: ['t'], 2: ['s'], 3: ['s', 't'], 4: ['z'], 5: ['z', 't'],
             6: ['sdg'], 7: ['tdg']}
INVOLUTION = {'h', 'x', 'y'}
INVERSE = {'t': 'tdg', 'tdg': 't', 's': 'sdg', 'sdg': 's'}


# ---------------------------------------------------------------- algorithms

def build_qft(n, seed):
    rng = random.Random(seed)
    x = rng.randrange(1 << n)
    core = QuantumCircuit(n)
    core.append(QFTGate(n), range(n))
    prep = [('x', (k,)) for k in range(n) if (x >> k) & 1]
    return core, prep, {'x': x}


def build_tfim(n, steps, dt, J, hx):
    core = QuantumCircuit(n)
    for _ in range(steps):
        for i in range(n - 1):
            core.rzz(-2.0 * J * dt, i, i + 1)
        for i in range(n):
            core.rx(-2.0 * hx * dt, i)
    return core, [], {'steps': steps, 'dt': dt, 'J': J, 'h': hx}


def random_regular_graph(n, d, rng):
    """Pairing model with rejection: a uniform simple d-regular graph."""
    if n * d % 2 or n <= d:
        raise ValueError(f"no simple {d}-regular graph on {n} vertices")
    for _ in range(100000):
        stubs = [v for v in range(n) for _ in range(d)]
        rng.shuffle(stubs)
        edges = set()
        for a, b in zip(stubs[::2], stubs[1::2]):
            e = (min(a, b), max(a, b))
            if a == b or e in edges:
                break
            edges.add(e)
        else:
            return sorted(edges)
    raise RuntimeError("pairing model did not converge")


def build_qaoa(n, layers, gammas, betas, seed):
    rng = random.Random(seed)
    edges = random_regular_graph(n, 3, rng)
    gammas, betas = gammas[:layers], betas[:layers]
    if len(gammas) < layers or len(betas) < layers:
        raise ValueError("give at least --layers values of --gammas and --betas")
    core = QuantumCircuit(n)
    core.h(range(n))
    for g, b in zip(gammas, betas):
        for i, j in edges:                 # exp(-i g (1 - Z_i Z_j)/2), mod phase
            core.rzz(-g, i, j)
        for i in range(n):                 # exp(-i b X_i)
            core.rx(2.0 * b, i)
    return core, [], {'layers': layers, 'gammas': gammas, 'betas': betas,
                      'edges': [list(e) for e in edges]}


# ---------------------------------------------------------------- synthesis

def word_matrix(word):
    m = np.eye(2, dtype=complex)
    for g in word:
        qc = QuantumCircuit(1)
        getattr(qc, g)(0)
        m = Operator(qc).data @ m
    return m


def phase_free_dist(u, v):
    """min_phi ||u - e^{i phi} v||_2, the error measure of gridsynth."""
    tr = np.trace(u.conj().T @ v)
    ph = tr.conjugate() / abs(tr) if abs(tr) > 1e-300 else 1.0
    return float(np.linalg.norm(u - ph * v, 2))


class Synth:
    def __init__(self, method, eps, depth, basic_depth):
        self.method, self.eps, self.depth = method, eps, depth
        self.sk = (SolovayKitaevDecomposition(basis_gates=['h', 't', 'tdg'],
                                              depth=basic_depth)
                   if method == 'sk' else None)
        self.cache = {}
        self.stats = []        # one entry per rz of the circuit

    def word(self, theta):
        th = math.remainder(theta, 2.0 * math.pi)            # (-pi, pi]
        k = th / (math.pi / 4.0)
        if abs(k - round(k)) < 1e-9:                          # rz(k pi/4) ~ t^k
            w = DIAG_WORD[round(k) % 8]
            self.stats.append(('exact', len(w), sum(g in ('t', 'tdg') for g in w), 0.0))
            return list(w)
        key = round(th, 13)
        if key not in self.cache:
            if self.method == 'gridsynth':
                qc = gridsynth_rz(th, self.eps)
            else:
                qc = self.sk.run(RZGate(th).to_matrix(), self.depth)
            w = [ci.operation.name for ci in qc.data]
            bad = set(w) - EMIT
            if bad:
                raise RuntimeError(f"synthesis emitted {bad}")
            w = peephole([(g, (0,)) for g in w])
            w = [g for g, _ in w]
            err = phase_free_dist(RZGate(th).to_matrix(), word_matrix(w)) if w else \
                phase_free_dist(RZGate(th).to_matrix(), np.eye(2))
            self.cache[key] = (w, err)
        w, err = self.cache[key]
        self.stats.append(('synth', len(w), sum(g in ('t', 'tdg') for g in w), err))
        return list(w)


def peephole(gates):
    """Merge diagonal 1q runs mod T^8 and cancel h h, x x, y y per wire."""
    out, wire = [], {}                    # out: [name, qubits, exp|None, alive]
    for name, qs in gates:
        if len(qs) == 1:
            q = qs[0]
            st = wire.setdefault(q, [])
            top = out[st[-1]] if st else None
            if name in DIAG:
                if top is not None and top[0] == 'D':
                    top[2] = (top[2] + DIAG[name]) % 8
                    if top[2] == 0:
                        top[3] = False
                        st.pop()
                    continue
                out.append(['D', qs, DIAG[name], True])
                st.append(len(out) - 1)
                continue
            if name in INVOLUTION and top is not None and top[0] == name:
                top[3] = False
                st.pop()
                continue
            out.append([name, qs, None, True])
            st.append(len(out) - 1)
        else:
            out.append([name, qs, None, True])
            for q in qs:
                wire.setdefault(q, []).append(len(out) - 1)
    res = []
    for name, qs, e, alive in out:
        if not alive:
            continue
        if name == 'D':
            res += [(g, qs) for g in DIAG_WORD[e]]
        else:
            res.append((name, qs))
    return res


def compile_core(core, synth, seed):
    tc = transpile(core, basis_gates=BASIS, optimization_level=1, seed_transpiler=seed)
    gates = []
    for ci in tc.data:
        name = ci.operation.name
        qs = tuple(tc.find_bit(q).index for q in ci.qubits)
        if name == 'rz':
            gates += [(g, qs) for g in synth.word(float(ci.operation.params[0]))]
        elif name in EMIT:
            gates.append((name, qs))
        elif name == 'barrier':
            continue
        else:
            raise RuntimeError(f"transpiler left {name}")
    return tc, gates


def inverse(gates):
    return [(INVERSE.get(g, g), qs) for g, qs in reversed(gates)]


def to_qasm(n, gates):
    body = [f"{g} " + ",".join(f"q[{q}]" for q in qs) + ";" for g, qs in gates]
    return "\n".join(["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{n}];"] + body)


def tag_of(a):
    if a.synth == 'sk':
        return f"sk{a.depth}"
    k = -math.log10(a.eps)
    return f"e{round(k)}" if abs(k - round(k)) < 1e-9 else f"eps{a.eps:g}"


def generate(alg, n, a):
    t0 = time.time()
    synth = Synth(a.synth, a.eps, a.depth, a.sk_basic_depth)
    seed = a.seed + 1000 * n
    if alg in ('qft', 'qft_mirror'):
        core, prep, meta = build_qft(n, seed)
        name = f"synth_{alg}_n{n}_{tag_of(a)}"
    elif alg == 'tfim':
        core, prep, meta = build_tfim(n, a.steps, a.dt, a.J, a.hx)
        name = f"synth_tfim_n{n}_r{a.steps}_{tag_of(a)}"
    else:
        core, prep, meta = build_qaoa(n, a.layers, a.gammas, a.betas, seed)
        name = f"synth_qaoa_n{n}_p{a.layers}_{tag_of(a)}"
    tc, word = compile_core(core, synth, seed)
    if not a.no_peephole:
        word = peephole(word)
    gates = list(prep) + word + (inverse(word) if alg == 'qft_mirror' else [])
    txt = to_qasm(n, gates)
    t_gen = time.time() - t0

    syn = [s for s in synth.stats if s[0] == 'synth']
    nsyn = len(syn) or 1
    info = {
        'family': f"synth_{alg}", 'algorithm': alg, 'n': n, 'seed': seed,
        'params': meta,
        'synth': a.synth,
        'eps': a.eps if a.synth == 'gridsynth' else None,
        'sk_depth': a.depth if a.synth == 'sk' else None,
        'sk_basic_depth': a.sk_basic_depth if a.synth == 'sk' else None,
        'peephole': not a.no_peephole,
        'transpiled_ops': dict(tc.count_ops()),
        'rz_total': len(synth.stats) * (2 if alg == 'qft_mirror' else 1),
        'rz_exact_cliffT': sum(s[0] == 'exact' for s in synth.stats),
        'rz_synthesised': len(syn),
        'rz_synth_identity': sum(s[1] == 0 for s in syn),
        'distinct_angles_synthesised': len(synth.cache),
        'gates_per_rotation_mean': sum(s[1] for s in syn) / nsyn,
        't_per_rotation_mean': sum(s[2] for s in syn) / nsyn,
        'gates_per_rotation_max': max((s[1] for s in syn), default=0),
        'rotation_error_max': max((s[3] for s in syn), default=0.0),
        'gates': len(gates),
        't_count': sum(g in ('t', 'tdg') for g, _ in gates),
        'counts': dict(sorted(Counter(g for g, _ in gates).items())),
        'generation_time_s': round(t_gen, 3),
    }
    if n <= a.fid_max:
        if alg == 'qft_mirror':
            ex = Statevector.from_int(meta['x'], 2 ** n)
        else:
            pc = QuantumCircuit(n)
            for g, qs in prep:
                getattr(pc, g)(*qs)
            ex = Statevector(pc.compose(tc))
        sy = Statevector(QuantumCircuit.from_qasm_str(txt))
        info['fidelity'] = float(abs(np.vdot(ex.data, sy.data)) ** 2)
        info['fidelity_reference'] = ('basis state |x>' if alg == 'qft_mirror'
                                      else 'transpiled circuit with exact rz')
    return name, txt, info


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--alg", required=True, choices=['qft', 'qft_mirror', 'tfim', 'qaoa'])
    ap.add_argument("--qubits", type=int, nargs="+", required=True)
    ap.add_argument("--synth", choices=['gridsynth', 'sk'], default='gridsynth')
    ap.add_argument("--eps", type=float, default=1e-3, help="gridsynth error bound")
    ap.add_argument("--depth", type=int, default=2, help="Solovay-Kitaev recursion degree")
    ap.add_argument("--sk-basic-depth", type=int, default=16,
                    help="length of the basic h,t,tdg approximations (Qiskit default 12)")
    ap.add_argument("--steps", type=int, default=2, help="tfim: Trotter steps r")
    ap.add_argument("--dt", type=float, default=0.1, help="tfim: Trotter step")
    ap.add_argument("--J", type=float, default=1.0, help="tfim: ZZ coupling")
    ap.add_argument("--hx", type=float, default=1.0, help="tfim: transverse field")
    ap.add_argument("--layers", type=int, default=1, help="qaoa: p")
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.7, 0.45])
    ap.add_argument("--betas", type=float, nargs="+", default=[0.3, 0.55])
    ap.add_argument("--fid-max", type=int, default=12,
                    help="compute the fidelity for n up to this")
    ap.add_argument("--no-peephole", action="store_true")
    ap.add_argument("--seed", type=int, default=1234)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    print(f"{'file':36s}{'n':>4}{'gates':>9}{'T':>8}{'rot':>6}{'g/rot':>7}{'T/rot':>7}"
          f"{'fidelity':>12}{'sec':>7}")
    for n in a.qubits:
        name, txt, info = generate(a.alg, n, a)
        with open(os.path.join(a.out, name + ".qasm"), "w") as fh:
            fh.write(txt)
        with open(os.path.join(a.out, name + ".json"), "w") as fh:
            json.dump(info, fh, indent=1)
        fid = f"{info['fidelity']:.8f}" if 'fidelity' in info else '-'
        print(f"{name + '.qasm':36s}{n:>4}{info['gates']:>9}{info['t_count']:>8}"
              f"{info['rz_synthesised']:>6}{info['gates_per_rotation_mean']:>7.1f}"
              f"{info['t_per_rotation_mean']:>7.1f}{fid:>12}{info['generation_time_s']:>7.1f}")
