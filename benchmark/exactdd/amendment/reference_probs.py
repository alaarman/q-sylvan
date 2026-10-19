#!/usr/bin/env python3
"""Reference top-qubit probabilities by state-vector simulation, for the
circuits on which no exact run (EVDD or LIMDD, qisq2) yields one.

  reference_probs.py OUT.csv QASM...

For each circuit: the probability that qubit 0 (the first qubit of the first
qreg, the runner's "top qubit") measures 0, from a double-precision
simulation of the circuit on |0...0>. The state is kept sparse (index ->
amplitude, amplitudes below 1e-15 in magnitude dropped), which is exact up to
rounding and fast when most qubits stay in a basis state (ancillas), as in the
quantum walks. Qiskit (qiskit.qasm2) parses the file and gives the gate
matrices; measurements and barriers are skipped.

Writes circuit,qubits,gates,p0,support_max,seconds. The figures use p0 only
where no exact probability exists (make_figures.py --reference), so it never
overrides an exact run; on circuits where an exact run has one, the two can
be compared to validate this script (--check in make_figures.py prints it).
"""
import csv, os, sys, time
import numpy as np
from qiskit import qasm2

DROP = 1e-15

def apply(idx, amp, qubits, U):
    """Apply the k-qubit matrix U to qubits (qubits[0] is U's least
    significant bit, as in Qiskit) of the sparse state (idx, amp)."""
    k = len(qubits)
    mask = 0
    for q in qubits: mask |= 1 << q
    base = idx & ~np.uint64(mask)
    local = np.zeros(len(idx), dtype=np.int64)
    for j, q in enumerate(qubits):
        local |= ((idx >> np.uint64(q)) & np.uint64(1)).astype(np.int64) << j
    ub, inv = np.unique(base, return_inverse=True)
    dense = np.zeros((len(ub), 1 << k), dtype=np.complex128)
    dense[inv, local] = amp
    out = dense @ U.T
    nz = np.abs(out) > DROP
    rows, cols = np.nonzero(nz)
    new_idx = ub[rows].copy()
    for j, q in enumerate(qubits):
        new_idx |= ((cols >> j) & 1).astype(np.uint64) << np.uint64(q)
    return new_idx, out[rows, cols]

# gates the runner reads that qelib1.inc lacks
from qiskit.circuit.library import CCZGate, CSGate, CSdgGate
EXTRA = qasm2.LEGACY_CUSTOM_INSTRUCTIONS + (
    qasm2.CustomInstruction("ccz", 0, 3, CCZGate, builtin=True),
    qasm2.CustomInstruction("cs", 0, 2, CSGate, builtin=True),
    qasm2.CustomInstruction("csdg", 0, 2, CSdgGate, builtin=True))

def p0(path):
    qc = qasm2.load(path, custom_instructions=EXTRA)
    pos = {q: i for i, q in enumerate(qc.qubits)}
    idx = np.zeros(1, dtype=np.uint64); amp = np.ones(1, dtype=np.complex128)
    cache, gates, smax = {}, 0, 1
    for ins in qc.data:
        op = ins.operation
        if op.name in ("measure", "barrier", "reset"):
            if op.name == "reset": raise ValueError("reset")
            continue
        key = (op.name, tuple(float(x) for x in op.params))
        U = cache.get(key)
        if U is None:
            U = cache[key] = np.asarray(op.to_matrix(), dtype=np.complex128)
        idx, amp = apply(idx, amp, [pos[q] for q in ins.qubits], U)
        gates += 1; smax = max(smax, len(idx))
    w = np.abs(amp) ** 2
    norm = w.sum()
    return float(w[(idx & np.uint64(1)) == 0].sum() / norm), qc.num_qubits, gates, smax, float(norm)

def main():
    out, paths = sys.argv[1], sys.argv[2:]
    done = set()
    if os.path.exists(out):
        done = {r["circuit"] for r in csv.DictReader(open(out))}
    new = not os.path.exists(out)
    with open(out, "a", newline="") as fh:
        w = csv.writer(fh)
        if new: w.writerow(["circuit", "qubits", "gates", "p0", "support_max", "norm", "seconds"])
        for p in paths:
            c = os.path.basename(p)[:-5]
            if c in done: continue
            t = time.time()
            pr, n, g, s, nm = p0(p)
            w.writerow([c, n, g, f"{pr:.12g}", s, f"{nm:.15g}", f"{time.time() - t:.1f}"]); fh.flush()
            print(f"{c:28s} n={n:2d} gates={g:6d} p0={pr:.9f} support<={s} norm={nm:.12f} {time.time() - t:.1f}s", flush=True)

if __name__ == "__main__":
    main()
