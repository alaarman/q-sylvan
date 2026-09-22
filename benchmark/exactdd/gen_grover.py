"""Grover in the LIMDD Clifford+T basis: n data qubits, n-2 v-chain ancillas."""
import sys
from qiskit import QuantumCircuit, QuantumRegister, transpile
from qiskit.qasm2 import dumps

BASIS = ['h','s','sdg','t','tdg','x','y','z','cx','cy','cz','swap']

def mcz(qc, ctrls, targ, anc):
    """Multi-controlled Z on targ, controls ctrls, using v-chain ancillas anc (len=len(ctrls)-1)."""
    k = len(ctrls)
    if k == 0: qc.z(targ); return
    if k == 1: qc.cz(ctrls[0], targ); return
    qc.ccx(ctrls[0], ctrls[1], anc[0])
    for i in range(2, k):
        qc.ccx(ctrls[i], anc[i-2], anc[i-1])
    qc.cz(anc[k-2], targ)
    for i in reversed(range(2, k)):
        qc.ccx(ctrls[i], anc[i-2], anc[i-1])
    qc.ccx(ctrls[0], ctrls[1], anc[0])

def grover(n, iterations, flag=None):
    if flag is None: flag = [1]*n
    q = QuantumRegister(n, 'q')
    a = QuantumRegister(max(n-2,1), 'a')
    qc = QuantumCircuit(q, a)
    qc.h(q)
    for _ in range(iterations):
        # oracle: phase flip on |flag>
        for i,b in enumerate(flag):
            if not b: qc.x(q[i])
        mcz(qc, list(q[:-1]), q[n-1], list(a))
        for i,b in enumerate(flag):
            if not b: qc.x(q[i])
        # diffusion
        qc.h(q); qc.x(q)
        mcz(qc, list(q[:-1]), q[n-1], list(a))
        qc.x(q); qc.h(q)
    return qc

if __name__ == '__main__':
    n = int(sys.argv[1]); it = int(sys.argv[2])
    qc = grover(n, it)
    t = transpile(qc, basis_gates=BASIS, optimization_level=0)
    s = dumps(t)
    open(f'grover_n{n}_it{it}.qasm','w').write(s)
    print(n, it, 'total qubits', t.num_qubits, 'gates', sum(t.count_ops().values()),
          'depth', t.depth(), 'bytes', len(s), dict(t.count_ops()))
