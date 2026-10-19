#!/bin/bash
# The 91 circuits of the new families (mirror, IQP, Grover with many
# iterations, quantum walks, synthesised circuits), one directory per family,
# as the data's pilot/qasm/ has them (amendment/REPRODUCE.md).
#   ./gen_new_families.sh OUTDIR RANDOMDIR
#   OUTDIR     where <family>/<circuit>.qasm go (created)
#   RANDOMDIR  the paper's random circuits (ROOT/ldd/qasm), for the two mirrors
#              built on them
# Python from $PYQ (a venv with qiskit 2.5.2 and pygridsynth 2.0.0, for
# gen_grover.py and gen_synth.py), else python3; the other generators need no
# packages. Every .qasm equals the archived one byte for byte (checked 8 Oct
# 2026 with diff -r against the data's pilot/qasm); the .json sidecars of the
# synthesised circuits differ only in generation_time_s.
set -eu
[ $# -eq 2 ] || { echo "usage: $0 OUTDIR RANDOMDIR" >&2; exit 2; }
O=$1; R=$2
PYQ=${PYQ:-python3}
G=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$O"/{mirror,iqp_close,grover_iter,qwalk,qwalk_echo,synth_qft,synth_qft_mirror,synth_qaoa,synth_tfim}

# mirror: random Clifford+T U, then U^dagger; 'plus' starts from |+...+>, x0 flips qubit 0 in the middle
python3 "$G/gen_mirror.py" --out "$O/mirror" --middle none --pairs 8:150 12:300 16:500 20:700
python3 "$G/gen_mirror.py" --out "$O/mirror" --middle none --plus --pairs 8:150 12:300 14:400 16:500
python3 "$G/gen_mirror.py" --out "$O/mirror" --middle x0 --pairs 12:300 20:700
python3 "$G/gen_mirror.py" --out "$O/mirror" --middle x0 --plus --pairs 12:300 16:500
python3 "$G/gen_mirror.py" --out "$O/mirror" --middle none \
    --base "$R/clifford_T_circuit_30_700.qasm" "$R/clifford_T_circuit_40_750.qasm"

# IQP, closed with a final Hadamard layer: qubits:gates pairs
for p in 8:16 12:24 16:16 16:32 20:40 20:80 24:24 24:48 32:32 40:40; do
  python3 "$G/gen_iqp.py" --out "$O/iqp_close" --close --qubits "${p%:*}" --gates "${p#*:}"
done

# Grover with k iterations (gen_grover.py n k writes grover_n<n>_it<k>.qasm into the working directory)
for p in 8:4 8:12 8:13 10:25 10:26 12:50 12:51 16:32 16:128 16:201 20:64 20:256 20:512; do
  n=${p%:*}; k=${p#*:}
  (cd "$O/grover_iter" && "$PYQ" "$G/gen_grover.py" "$n" "$k" > /dev/null \
     && mv "grover_n${n}_it${k}.qasm" "grover_iter_n${n}_it${k}.qasm")
done

# quantum walks on a cycle of 2^m positions, t steps; echo walks t steps and back
python3 "$G/gen_qwalk.py" --out "$O/qwalk" --pairs 4:8 4:32 6:16 6:64 8:32 8:128 10:64 10:256 10:512 12:128 12:512 12:2048
python3 "$G/gen_qwalk.py" --out "$O/qwalk_echo" --echo --pairs 6:32 8:128 10:256 12:512

# synthesised: gridsynth with eps 1e-3 (e3) or 1e-6 (e6), or Solovay-Kitaev depth 2 (sk2)
for alg in qft qft_mirror; do
  "$PYQ" "$G/gen_synth.py" --out "$O/synth_$alg" --alg $alg --qubits 3 4 5 6 7 8 12 16
  "$PYQ" "$G/gen_synth.py" --out "$O/synth_$alg" --alg $alg --qubits 4 6 8 --eps 1e-6
  "$PYQ" "$G/gen_synth.py" --out "$O/synth_$alg" --alg $alg --qubits 4 8 --synth sk
done
"$PYQ" "$G/gen_synth.py" --out "$O/synth_qaoa" --alg qaoa --qubits 4 6 8 10 12 16
"$PYQ" "$G/gen_synth.py" --out "$O/synth_tfim" --alg tfim --qubits 4 6 8 10 12 16
