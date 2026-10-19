#!/usr/bin/env python3
"""Discrete-time coined quantum walk on a cycle, as OPENQASM 2.0, dependency-free.

The walk lives on N = 2^m positions on a cycle and one coin qubit.  One step is
a Hadamard on the coin followed by the conditional shift

    S |x>|1> = |x+1 mod N>|1>,      S |x>|0> = |x-1 mod N>|0>,

and the walk starts at position 0 with the coin in |0> (or, with --sym-coin,
in (|0> + i|1>)/sqrt2, prepared by h; s, the start that makes the Hadamard
walk's distribution symmetric).  After t steps the state is

    sum_x |x> (a_t(x)|0> + b_t(x)|1>),   a_t(x), b_t(x) in Z[1/sqrt2]
                                          (Z[i, 1/sqrt2] with --sym-coin)

and the Hadamard walk spreads ballistically, with strong interference between
the two coin components.  Why the family is interesting for exact edge weights:

  * Every gate is in Clifford+T, so qisq2 (exact arithmetic in Q(i, sqrt2)) is
    exact; the amplitudes have denominator sqrt2^t, so the numerators grow
    linearly in bits with t: the exact cost should grow with t, not only with
    the diagram's size.
  * Exact zeros by parity: after t steps only positions x = t mod 2 are
    populated, and every Toffoli of the shift (h; ccz; h) creates and then
    cancels a superposition of its target.  Float rounding can leave residues
    where the exact amplitude is zero, and can keep apart nodes that are
    exactly equal (the walk's ancillas return to |0> only up to rounding).
  * For m = 2, 3 (cycles of 4 and 8) the Hadamard walk is periodic, with
    period 8 and 24: the exact diagram returns to a single path, float drifts.
    --echo runs t steps and then their inverse, so the final state is exactly
    the start state (one path in the diagram) for every m; it measures how
    much rounding the round trip leaves.

Qubit layout (q[0] is the top variable of the diagram, and the runner's -m
measures q[0]):

    q[0]                 coin                         (--coin first, default)
    q[1] .. q[m]         position, most significant bit first
    q[m+1] ..            clean ancillas of the Toffoli ladder (m-2 of them)

With --coin last the position comes first (q[0] = MSB) and the coin sits
between the position and the ancillas.

The shift.  --shift ctrl (default) is the literal construction: increment of
the position controlled on coin = 1, then x on the coin, decrement controlled
on coin = 1, and x on the coin again.  The decrement is the increment's gate
list reversed (every gate in it is self-inverse).  --shift flip uses the
identity dec = X^m inc X^m: anti-controlled X from the coin onto every
position bit, an uncontrolled increment, and the same anti-controlled X again;
the same unitary with fewer Toffolis and one ancilla less.

The increment flips bit j of the position when the coin and the bits below j
are all 1, most significant bit first.  It is built from multi-controlled X
gates on a clean-ancilla Toffoli ladder:

  --ladder shared (default): one ladder anc_r = AND of the first r+2 controls
      is computed once, the top bit is flipped by a Toffoli, and each lower bit
      by a cx from its ancilla, which is uncomputed right after it is used:
      2m-3 Toffolis and m-1 cx per controlled increment, m-2 ancillas.
  --ladder per-mcx: every multi-controlled X computes and uncomputes its own
      ladder, O(m^2) Toffolis per increment.

Each Toffoli is written as h t; ccz a,b,t; h t (--toffoli ccz, default; the
runner applies ccz as one gate and counts it as T-count 7), or fully in
Clifford+T with the standard 7-T decomposition (--toffoli ct).  The emitted
gates are x h s cx ccz (ccz mode) or x h s t tdg cx (ct mode), all of which
the EVDD and LIMDD paths accept.

Files are <stem>_m<m>_t<t>.qasm, where the stem is "qwalk" followed by the
variant tags that differ from the defaults (e.g. qwalk_ct_sym_m6_t32).  The
circuit is deterministic: no seed is needed.  No creg, no measure, no trailing
newline, like the other generators in this directory.
"""
import argparse, os


# --------------------------------------------------------------------------
# abstract reversible ops: ('x', q), ('cx', c, t), ('ccx', a, b, t)
# --------------------------------------------------------------------------

def _inc_shared(chain, bits, anc):
    """Increment with the shared ladder.

    chain: the control sequence s_0, s_1, ...  (coin first if controlled, then
    position bits from the LSB).  bits: position bits, LSB first.  Bit j is
    flipped iff s_0..s_{k_j - 1} are all 1, where k_j = j + len(chain) -
    len(bits); targets are processed from the MSB down.
    """
    off = len(chain) - len(bits)          # 1 if controlled on the coin, else 0
    K = len(bits) - 1 + off               # controls of the top bit
    ops = []
    # compute anc_r = AND(s_0 .. s_{r+1}), r = 0 .. K-3
    for r in range(K - 2):
        a = chain[0] if r == 0 else anc[r - 1]
        ops.append(('ccx', a, chain[r + 1], anc[r]))
    for j in range(len(bits) - 1, -1, -1):
        k = j + off
        b = bits[j]
        if k == 0:
            ops.append(('x', b))
        elif k == 1:
            ops.append(('cx', chain[0], b))
        elif k == K:                      # the top bit: Toffoli from the ladder
            if k == 2:
                ops.append(('ccx', chain[0], chain[1], b))
            else:
                ops.append(('ccx', anc[k - 3], chain[k - 1], b))
        else:                             # 2 <= k < K: anc_{k-2} holds the AND
            ops.append(('cx', anc[k - 2], b))
            a = chain[0] if k == 2 else anc[k - 3]
            ops.append(('ccx', a, chain[k - 1], anc[k - 2]))   # uncompute
    return ops


def _mcx(ctrls, tgt, anc):
    """Multi-controlled X with its own clean-ancilla ladder."""
    k = len(ctrls)
    if k == 0:
        return [('x', tgt)]
    if k == 1:
        return [('cx', ctrls[0], tgt)]
    if k == 2:
        return [('ccx', ctrls[0], ctrls[1], tgt)]
    up = [('ccx', ctrls[0], ctrls[1], anc[0])]
    for r in range(1, k - 2):
        up.append(('ccx', anc[r - 1], ctrls[r + 1], anc[r]))
    return up + [('ccx', anc[k - 3], ctrls[k - 1], tgt)] + up[::-1]


def _inc_per_mcx(chain, bits, anc):
    off = len(chain) - len(bits)
    ops = []
    for j in range(len(bits) - 1, -1, -1):
        ops += _mcx(chain[:j + off], bits[j], anc)
    return ops


def step_shift(coin, bits, anc, shift, ladder):
    """The conditional shift as abstract ops.  bits: position, LSB first."""
    inc = _inc_shared if ladder == 'shared' else _inc_per_mcx
    if shift == 'ctrl':
        up = inc([coin] + bits, bits, anc)              # +1 on coin = 1
        return up + [('x', coin)] + up[::-1] + [('x', coin)]   # -1 on coin = 0
    # flip: (anti-cx coin -> bits) inc (anti-cx coin -> bits)
    flip = [('x', coin)] + [('cx', coin, b) for b in bits] + [('x', coin)]
    return flip + inc(bits, bits, anc) + flip


def n_ancillas(m, shift):
    K = m if shift == 'ctrl' else m - 1     # controls of the top bit
    return max(0, K - 2)


# --------------------------------------------------------------------------
# lowering to QASM
# --------------------------------------------------------------------------

def _toffoli_ct(a, b, c):
    """Exact Toffoli in Clifford+T (Nielsen & Chuang Fig. 4.9), T-count 7."""
    q = lambda i: f"q[{i}]"
    return [f"h {q(c)};", f"cx {q(b)},{q(c)};", f"tdg {q(c)};", f"cx {q(a)},{q(c)};",
            f"t {q(c)};", f"cx {q(b)},{q(c)};", f"tdg {q(c)};", f"cx {q(a)},{q(c)};",
            f"t {q(b)};", f"t {q(c)};", f"h {q(c)};", f"cx {q(a)},{q(b)};",
            f"t {q(a)};", f"tdg {q(b)};", f"cx {q(a)},{q(b)};"]


def lower(ops, toffoli):
    out = []
    for op in ops:
        if op[0] == 'x':
            out.append(f"x q[{op[1]}];")
        elif op[0] == 'cx':
            out.append(f"cx q[{op[1]}],q[{op[2]}];")
        elif toffoli == 'ccz':
            a, b, c = op[1:]
            out += [f"h q[{c}];", f"ccz q[{a}],q[{b}],q[{c}];", f"h q[{c}];"]
        else:
            out += _toffoli_ct(*op[1:])
    return out


def layout(m, coin_pos):
    """Returns (coin, position bits LSB first, first ancilla index)."""
    if coin_pos == 'first':
        coin = 0
        msb_first = list(range(1, m + 1))
    else:
        coin = m
        msb_first = list(range(0, m))
    return coin, msb_first[::-1], m + 1


def generate(m, t, toffoli='ccz', shift='ctrl', ladder='shared', sym=False,
             coin_pos='first', echo=False, start=0):
    """Returns (qasm text, number of qubits, number of Toffolis)."""
    assert m >= 1 and t >= 0
    coin, bits, a0 = layout(m, coin_pos)
    na = n_ancillas(m, shift)
    anc = list(range(a0, a0 + na))
    n = a0 + na
    sh = step_shift(coin, bits, anc, shift, ladder)
    ntof = sum(1 for op in sh if op[0] == 'ccx')
    body = [f"x q[{bits[j]}];" for j in range(m) if (start >> j) & 1]
    if sym:
        body += [f"h q[{coin}];", f"s q[{coin}];"]
    fwd = [f"h q[{coin}];"] + lower(sh, toffoli)
    body += fwd * t
    if echo:
        # inverse step: shift^-1 (reversed ops; each is self-inverse), then h
        # (the 7-T decomposition is an exact Toffoli, so it is self-inverse too)
        bwd = lower(sh[::-1], toffoli) + [f"h q[{coin}];"]
        body += bwd * t
    head = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{n}];"]
    return "\n".join(head + body), n, ntof * t * (2 if echo else 1)


def stem_for(a):
    s = "qwalk"
    if a.toffoli != 'ccz':
        s += "_ct"
    if a.shift != 'ctrl':
        s += "_flip"
    if a.ladder != 'shared':
        s += "_mcx"
    if a.sym_coin:
        s += "_sym"
    if a.coin != 'first':
        s += "_coinlast"
    if a.start:
        s += f"_x{a.start}"
    if a.echo:
        s += "_echo"
    return s


def _count(txt):
    g = t = c = 0
    for ln in txt.split("\n"):
        if not ln or ln.startswith(("OPENQASM", "include", "qreg")):
            continue
        g += 1
        if ln.startswith(("t ", "tdg ")):
            t += 1
        elif ln.startswith("ccz "):
            c += 1
    return g, t, c


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--m", type=int, nargs="+",
                    help="log2 of the cycle length (cross product with --t)")
    ap.add_argument("--t", type=int, nargs="+", help="number of walk steps")
    ap.add_argument("--pairs", nargs="+", metavar="M:T",
                    help="explicit (m, t) instances instead of --m x --t")
    ap.add_argument("--toffoli", choices=['ccz', 'ct'], default='ccz')
    ap.add_argument("--shift", choices=['ctrl', 'flip'], default='ctrl')
    ap.add_argument("--ladder", choices=['shared', 'per-mcx'], default='shared')
    ap.add_argument("--sym-coin", action="store_true",
                    help="start the coin in (|0>+i|1>)/sqrt2 (h; s)")
    ap.add_argument("--coin", choices=['first', 'last'], default='first',
                    help="coin above the position (q[0]) or below it")
    ap.add_argument("--start", type=int, default=0, help="start position")
    ap.add_argument("--echo", action="store_true",
                    help="t steps forward, then t steps backward (final state |0>)")
    ap.add_argument("--name", default=None,
                    help="filename stem (default: qwalk plus variant tags)")
    a = ap.parse_args()
    if a.pairs:
        inst = [tuple(int(v) for v in p.split(":")) for p in a.pairs]
    elif a.m and a.t:
        inst = [(m, t) for m in a.m for t in a.t]
    else:
        ap.error("give --pairs, or both --m and --t")
    stem = a.name or stem_for(a)
    os.makedirs(a.out, exist_ok=True)
    print(f"{'file':36s}{'qubits':>7}{'gates':>9}{'t+tdg':>8}{'ccz':>7}{'T(ccz=7)':>10}")
    for m, t in inst:
        txt, n, ntof = generate(m, t, a.toffoli, a.shift, a.ladder, a.sym_coin,
                                a.coin, a.echo, a.start % (1 << m))
        f = os.path.join(a.out, f"{stem}_m{m}_t{t}.qasm")
        with open(f, "w") as fh:
            fh.write(txt)
        g, tc, cc = _count(txt)
        assert tc + cc == 7 * ntof if a.toffoli == 'ct' else cc == ntof
        print(f"{os.path.basename(f):36s}{n:>7}{g:>9}{tc:>8}{cc:>7}{tc + 7 * cc:>10}")
