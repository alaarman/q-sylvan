"""W-state on K=2^m qubits, Clifford+T, m scratch qubits returned to |0>."""
import sys
from qiskit import QuantumCircuit, QuantumRegister, transpile
from qiskit.qasm2 import dumps
BASIS=['h','s','sdg','t','tdg','x','y','z','cx','cy','cz','swap']
def wstate(K):
    m=K.bit_length()-1; assert 1<<m==K
    c=QuantumRegister(m,'c'); t=QuantumRegister(K,'t')
    qc=QuantumCircuit(t,c)
    qc.x(t[0])
    for j in range(m):                      # spread onto one-hot, index in c
        qc.h(c[j]); blk=1<<j
        for i in range(blk): qc.cswap(c[j], t[i], t[i+blk])
    for j in range(m):                      # uncompute index from one-hot
        for k in range(K):
            if (k>>j)&1: qc.cx(t[k], c[j])
    return qc
if __name__=='__main__':
    K=int(sys.argv[1]); qc=wstate(K)
    tq=transpile(qc,basis_gates=BASIS,optimization_level=1)
    n=tq.num_qubits
    open(f'w-state_{K}_n{n}.qasm','w').write(dumps(tq))
    print('K',K,'n',n,'gates',sum(tq.count_ops().values()),'depth',tq.depth())
