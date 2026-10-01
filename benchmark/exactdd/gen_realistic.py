#!/usr/bin/env python3
"""Clifford+T circuits that stand for something, rather than random words.

Every family in MQTBench that does real work -- QFT, QPE, VQE, QAOA, Shor,
portfolio optimisation, ground states -- uses arbitrary-angle rotations, and
so cannot be written exactly in Q[i,sqrt2]. Only GHZ and graph states survive
the filter there, and both have no T gates at all. These three do real work
AND are natively Clifford+T:

  hidden-shift   The standard benchmark for simulators of Clifford-dominated
                 circuits (Bravyi and Gosset). The oracle is a bent function
                 built from CZ and CCZ, and the T count is set by how many CCZ
                 it uses, so the T count is a dial rather than an accident.

  adder          Cuccaro ripple-carry addition. Toffoli-based, hence exactly
                 Clifford+T once decomposed, and it is the inner loop of the
                 modular exponentiation in Shor's algorithm -- arguably the
                 most-studied circuit that is Clifford+T by nature.

  ising          One Trotter step of a transverse-field Ising Hamiltonian at
                 angle pi/4, where exp(-i pi/4 ZZ) is CX, T, CX exactly. This
                 is Hamiltonian simulation, and the angle is the one that
                 keeps it in the ring.

All three scale in n, and the T count grows with n, which is what makes them
useful for checking a bound stated in the T count.
"""
import argparse, os, sys

from qiskit import QuantumCircuit, QuantumRegister, transpile
from qiskit.qasm2 import dumps

# what the LIMDD path accepts; ccx/ccz are decomposed into this by transpile
BASIS = ['h', 's', 'sdg', 't', 'tdg', 'x', 'y', 'z', 'cx', 'cy', 'cz', 'swap']


def hidden_shift(n, ccz_frac=0.5, seed=0):
    """Bravyi-Gosset hidden shift on n qubits (n even)."""
    import random
    rng = random.Random(seed)
    n -= n % 2
    shift = [rng.randint(0, 1) for _ in range(n)]
    qc = QuantumCircuit(n)

    def oracle():
        # f(x) = sum x_{2i} x_{2i+1}  (bent), plus CCZ terms for nonlinearity
        for i in range(0, n - 1, 2):
            qc.cz(i, i + 1)
        for i in range(0, n - 2):
            if rng.random() < ccz_frac:
                qc.ccz(i, i + 1, i + 2)

    for q in range(n):
        if shift[q]: qc.x(q)
    qc.h(range(n))
    oracle()
    qc.h(range(n))
    oracle()
    qc.h(range(n))
    for q in range(n):
        if shift[q]: qc.x(q)
    return qc


def _maj(qc, a, b, c):
    qc.cx(c, b); qc.cx(c, a); qc.ccx(a, b, c)


def _uma(qc, a, b, c):
    qc.ccx(a, b, c); qc.cx(c, a); qc.cx(a, b)


def adder(nbits):
    """Cuccaro ripple-carry adder: a += b, with one carry qubit."""
    a = QuantumRegister(nbits, 'a'); b = QuantumRegister(nbits, 'b')
    c = QuantumRegister(1, 'c')
    qc = QuantumCircuit(a, b, c)
    # a non-trivial input, so the circuit is not acting on |0..0>
    for i in range(0, nbits, 2): qc.x(a[i])
    for i in range(0, nbits, 3): qc.x(b[i])
    qc.h(a[0])                       # and a superposition, so it is quantum
    _maj(qc, c[0], b[0], a[0])
    for i in range(1, nbits):
        _maj(qc, a[i - 1], b[i], a[i])
    for i in range(nbits - 1, 0, -1):
        _uma(qc, a[i - 1], b[i], a[i])
    _uma(qc, c[0], b[0], a[0])
    return qc


def ising(n, steps=1, ring=True):
    """Trotter steps of a transverse-field Ising model at angle pi/4.

    exp(-i pi/4 Z_i Z_j) = CX_{ij} . T_j . CX_{ij} up to a global phase, and
    exp(-i pi/4 X) = H . T . H, so the whole step is Clifford+T exactly.
    """
    qc = QuantumCircuit(n)
    qc.h(range(n))                   # start from |+>^n, the usual quench
    for _ in range(steps):
        for i in range(n - 1):       # ZZ couplings along a chain
            qc.cx(i, i + 1); qc.t(i + 1); qc.cx(i, i + 1)
        if ring and n > 2:
            qc.cx(n - 1, 0); qc.t(0); qc.cx(n - 1, 0)
        for i in range(n):           # transverse field
            qc.h(i); qc.t(i); qc.h(i)
    return qc


BUILD = {'hidden-shift': hidden_shift, 'adder': adder, 'ising': ising}

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('family', choices=sorted(BUILD))
    ap.add_argument('n', type=int)
    ap.add_argument('--steps', type=int, default=1)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', default='.')
    a = ap.parse_args()

    if a.family == 'hidden-shift':  qc = hidden_shift(a.n, seed=a.seed)
    elif a.family == 'adder':       qc = adder(a.n)
    else:                           qc = ising(a.n, steps=a.steps)

    t = transpile(qc, basis_gates=BASIS, optimization_level=1, seed_transpiler=a.seed)
    ops = t.count_ops()
    tc = ops.get('t', 0) + ops.get('tdg', 0)
    os.makedirs(a.out, exist_ok=True)
    name = f"{a.family}_n{t.num_qubits}" + (f"_s{a.steps}" if a.family == 'ising' else "")
    path = os.path.join(a.out, name + ".qasm")
    open(path, 'w').write(dumps(t))
    print(f"{name}: {t.num_qubits} qubits, {sum(ops.values())} gates, "
          f"T={tc}, depth={t.depth()}")
