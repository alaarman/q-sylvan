"""True W-state on K=2^m qubits, no ancillas, Clifford only: x + ch + cx."""
import sys
def wstate_qasm(K):
    m=K.bit_length()-1; assert 1<<m==K
    L=['OPENQASM 2.0;','include "qelib1.inc";',f'qreg q[{K}];','x q[0];']
    for j in range(m):
        blk=1<<j
        for i in range(blk):
            L.append(f'ch q[{i}],q[{i+blk}];')
            L.append(f'cx q[{i+blk}],q[{i}];')
    return '\n'.join(L)+'\n'
K=int(sys.argv[1]); open(f'wclif_{K}_n{K}.qasm','w').write(wstate_qasm(K)); print(K,'gates',2*(K-1))
