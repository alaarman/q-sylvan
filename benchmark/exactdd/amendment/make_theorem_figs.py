#!/usr/bin/env python3
"""Figures relating the theorems to the runs.

The obvious plot does not work: 2^t reaches 10^72 while the widths are in the
thousands, so on a log axis the bound is a ceiling at the top of the frame and
every trajectory is a flat line on the floor. Plotting log2(width) instead
turns the bound into the diagonal y = t, the trivial ceiling into the
horizontal y = n, and the gap between a trajectory and those two lines becomes
the readable quantity -- which is what the reviewers are asking to see.

  make_theorem_figs.py TRACEDIR OUTDIR

TRACEDIR holds the per-gate traces <qubits>_<gates>.csv of run_trace.sh
(columns t_count, width, wgt_bits among others); the qubit count is read off
the file name.
"""
import csv, glob, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

if len(sys.argv) != 3: sys.exit("usage: make_theorem_figs.py TRACEDIR OUTDIR")
IN, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)

runs = {}
for f in sorted(glob.glob(os.path.join(IN, "*.csv"))):
    name = os.path.basename(f)[:-4]
    try: nq = int(name.split("_")[0])
    except ValueError: continue
    t, w, b = [], [], []
    for r in csv.DictReader(open(f)):
        try:
            t.append(int(r["t_count"])); w.append(int(r["width"])); b.append(int(r["wgt_bits"]))
        except (TypeError, ValueError):
            break
    if t: runs[name] = (nq, t, w, b)
if not runs: print("no usable traces"); sys.exit(0)

cmap = plt.get_cmap("viridis")
nqs  = sorted({v[0] for v in runs.values()})
col  = {n: cmap(i / max(1, len(nqs) - 1)) for i, n in enumerate(nqs)}
lg2  = lambda v: (v.bit_length() - 1) if v > 0 else 0
tmax = max(max(v[1]) for v in runs.values())

# --- width: log2 on y, so the bound is the diagonal -------------------------
fig, ax = plt.subplots(figsize=(4.6, 3.5))
ax.plot([0, tmax], [0, tmax], "k--", lw=1.2, zorder=5)
for name, (nq, t, w, b) in sorted(runs.items()):
    ax.plot(t, [lg2(x) for x in w], lw=1.0, alpha=0.9, color=col[nq])
for nq in nqs:
    ax.axhline(nq, color=col[nq], ls=":", lw=0.8, alpha=0.7)
ax.set_xlabel("T gates applied"); ax.set_ylabel(r"$\log_2$ (LIMDD width)")
ax.set_xlim(0, tmax); ax.set_ylim(0, max(nqs) + 2); ax.tick_params(labelsize=8)
hs  = [plt.Line2D([], [], color="k", ls="--", lw=1.2, label=r"bound $\log_2 2^t=t$")]
hs += [plt.Line2D([], [], color=col[n], lw=1.2, label=f"{n} qubits") for n in nqs]
hs += [plt.Line2D([], [], color="0.4", ls=":", lw=0.8, label=r"ceiling $n$")]
ax.legend(handles=hs, fontsize=6.5, loc="upper left", ncol=2, framealpha=0.9)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "width_vs_tcount.pdf"), bbox_inches="tight")
plt.close(fig); print("  wrote width_vs_tcount.pdf")

# --- coefficient size: bits are already logarithmic -------------------------
fig, ax = plt.subplots(figsize=(4.6, 3.5))
for name, (nq, t, w, b) in sorted(runs.items()):
    ax.plot(t, b, lw=1.0, alpha=0.9, color=col[nq])
ax.set_xlabel("T gates applied"); ax.set_ylabel("largest coefficient (bits)")
ax.set_xlim(0, tmax); ax.set_ylim(bottom=0); ax.tick_params(labelsize=8)
ax.legend(handles=[plt.Line2D([], [], color=col[n], lw=1.2, label=f"{n} qubits") for n in nqs],
          fontsize=6.5, loc="lower right", ncol=2, framealpha=0.9)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "bits_vs_tcount.pdf"), bbox_inches="tight")
plt.close(fig); print("  wrote bits_vs_tcount.pdf")

print(f"\n  {'circuit':<10} {'n':>3} {'t':>5} {'width':>8} {'log2 w':>7} {'bits':>5}  slack (t - log2 w)")
for name, (nq, t, w, b) in sorted(runs.items()):
    mw = max(w)
    print(f"  {name:<10} {nq:>3} {t[-1]:>5} {mw:>8} {lg2(mw):>7} {max(b):>5}  {t[-1]-lg2(mw):>6}")
