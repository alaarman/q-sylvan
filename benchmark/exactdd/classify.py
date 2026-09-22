#!/usr/bin/env python3
"""Paper's correct/wrong rule: float is WRONG if its top-qubit P(|0>) deviates
by more than 5% from the algebraic (exact) value.

  python3 classify.py <resultdir> [float_cfg] [exact_cfg]     # default A vs B
Prints: circuit, p_exact, p_float, rel_dev, verdict, runtimes, node counts.
"""
import glob, json, os, sys
d = sys.argv[1]
fc = sys.argv[2] if len(sys.argv) > 2 else 'A'
ec = sys.argv[3] if len(sys.argv) > 3 else 'B'
print(f"{'circuit':40s} {'p_'+ec:>10s} {'p_'+fc:>10s} {'reldev':>8s}  verdict   "
      f"{'t_'+ec:>9s} {'t_'+fc:>9s} {'fin_'+ec:>8s} {'fin_'+fc:>8s} {'pk_'+ec:>8s} {'pk_'+fc:>8s}")
for f in sorted(glob.glob(os.path.join(d, f'*_{ec}.json'))):
    name = os.path.basename(f)[:-len(f'_{ec}.json')]
    g = os.path.join(d, f'{name}_{fc}.json')
    if not os.path.exists(g):
        print(f'{name:40s} (no {fc} result)'); continue
    e = json.load(open(f))['statistics']; x = json.load(open(g))['statistics']
    pe, pf = e['first_qubit_measurement_prob'], x['first_qubit_measurement_prob']
    dev = abs(pf - pe) / abs(pe) if pe else (0.0 if pf == 0 else float('inf'))
    print(f'{name:40s} {pe:10.6f} {pf:10.6f} {dev:8.2%}  '
          f'{"WRONG  " if dev > 0.05 else "correct"}   '
          f"{e['simulation_time']:9.3f} {x['simulation_time']:9.3f} "
          f"{e['final_nodes']:8d} {x['final_nodes']:8d} {e['max_nodes']:8d} {x['max_nodes']:8d}")
