#!/usr/bin/env python3
"""The clean sweep of the exact-DD paper: every run behind its Figures 6 and 7
and the T-count figure, scheduled so that the runtimes are usable.

  QSY_W1=... QSY_W2=... ./sweep.py DATA OUT [--dry] [--phases A,T,B,C,D,E]
        [--arms a,b,...] [--lazy-arms a,b] [--only REGEX] [--timeout S]

DATA  the data root (README.md, "Data layout"): the circuits in ldd/qasm,
      amend/qasm and pilot/qasm/<family>/; where it also holds runs (the
      stored clean sweep: amend/out_clean, pilot/out_clean, ...), those set
      the order (longest first) and the first memory estimates, and the
      footprints in DATA/sweep/clean_audit.csv refine them. Nothing in DATA
      is written.
OUT   where the runs go, in the same layout, so that OUT can be given to
      regen_amendment.sh as its ROOT: OUT/amend/out_clean, OUT/amend/trace_clean,
      OUT/pilot/out_clean, OUT/amend/out_clean_lazy, OUT/pilot/out_clean_lazy.
      OUT/ldd/qasm, OUT/amend/qasm, OUT/pilot/qasm and
      OUT/pilot/reference_probs.csv are made links to DATA's. The scheduler's
      own files are OUT/sweep/: audit.csv, mem.csv, and sweep.conf and
      sweep.pause (below).
QSY_W1, QSY_W2  the runners, Release builds with LIMDD_PAULI_WORDS=1 and 2
      (run_arm.sh). Every run is run_arm.sh ARM CIRCUIT OUTDIR TIMEOUT, or
      run_trace.sh for a trace.
--arms        the arms of phases A, B and C (default: the six of the paper)
--lazy-arms   the arms of phases D and E (default limdd_qisq2_lazy,limdd_lazy)
--phases      which phases to run (default all)
--only        only circuits whose name matches this regular expression (re.search)
--timeout     per run, in seconds (default 1800, the paper's)
--dry         print the plan and stop

The paper's runs: an Apple M1 Pro (8 performance and 2 efficiency cores,
32 GB), macOS, on power; at most six runs at a time (so that they stay on the
performance cores), each single-threaded (-w 1). The 6 October first attempt
showed why the rules below exist: with other load on the machine, runs that
reproduced earlier node counts took up to 12 times longer.

Every run: at most maxjobs simulators on the machine (every run_qasm_on_qmdd
counts, also ones this script did not start), and memory-aware:
  - a run starts only when the simulators' footprints (for a run of this sweep
    at least its estimate) plus the estimate of the new run stay within
    budget_gb; a run alone always starts;
  - when the simulators' footprints exceed hard_gb all the same (two samples
    in a row), a run of this sweep is stopped: the youngest whose footprint
    alone brings the total under hard_gb, else the largest. Its files are
    removed, and it is queued again with 1.5 times the footprint it had as its
    estimate;
  - nothing starts while the machine swaps out (more than SWAPOUT_HOLD bytes
    in the last minute).
The footprint is macOS's phys_footprint (compressed memory included), on
Linux resident plus swapped memory (/proc/PID/smaps_rollup); elsewhere it is
not measured, and only maxjobs limits the sweep. budget_gb, hard_gb, maxjobs,
the quality gate's limits and heavy_gb come from OUT/sweep/sweep.conf when it
exists (JSON), read again all the time, over the defaults below. The paper's
sweep used budget 9 to 18 GB and hard 11 to 15 GB over its course (what the
other applications open allowed), maxjobs 6, max_sys 0.06, min_cpu 0.85,
heavy_gb 5 to 6. A run whose estimate is heavy_gb or more is heavy: at most
one heavy run at a time, counted with what it uses now but at least half the
budget (not its whole estimate), so light runs fill the other half; when it
grows, the backstop stops light runs first. Nothing else should run on the
machine meanwhile.

The quality gate: a run of more than 20 s whose system time exceeds max_sys
of its CPU time (on macOS a page decompressed counts as system time of the
process that touches it), or whose CPU time is below min_cpu of its wall time,
is moved to <outdir>_rejected/ and run once more, with a larger memory
estimate. A second unclean run is kept and flagged in the audit. A run of an
earlier instance of this script, which this one cannot measure, is moved to
<outdir>_unmeasured/ once it ends and run again.

Phases, by priority. No run of a later phase starts while a run of phase A
waits to start; after that, the first waiting run in this order that fits starts.
  A  the 72 circuits of the paper's figures (49 random, 12 Grover, 11 W state)
     x the arms -> amend/out_clean/, longest first (by the earlier run)
  T  the ten per-gate traces of the T-count figure -> amend/trace_clean/
  B  the other 19 circuits of amend/qasm (hidden shift, adder, Ising) x the
     arms -> amend/out_clean/, longest first
  C  the new families (pilot/qasm, 91 circuits in ladders) x the arms ->
     pilot/out_clean/. A ladder (one variant growing in size) runs smallest
     first per arm; after a timeout or a full table, its larger instances are
     recorded as SKIPPED_PREDICTED for that arm. Longest ladders first.
  D  the lazy canonical form, limdd_qisq2_lazy and limdd_lazy, on the 79
     circuits of A and hidden shift -> amend/out_clean_lazy/, longest first
  E  the same two arms on the new families -> pilot/out_clean_lazy/, ladders
     as in C
Phase A done and the end are logged as "phase A done" and "end".

OUT/sweep/audit.csv gets one line per run:
  - start and end;
  - the wall time, and the CPU time of the run (wait4 on run_arm.sh: the
    simulator plus a few short helpers);
  - the simulator's peak footprint, its maximum RSS and its major faults;
  - the simulators' peak total footprint while it ran, and the machine's
    swapouts (bytes) meanwhile.
A CPU time well below the wall time means that something slowed the run.
OUT/sweep/mem.csv samples the machine every minute.

A restart skips every run whose .status exists, and waits for the runs of an
earlier instance that are still in flight. While OUT/sweep/sweep.pause
exists, no run starts (the runs in flight go on).

This is the driver of the paper's sweep (l2run/clean_sweep.py, kept with the
data as sweep/clean_sweep.py) with its paths made arguments, the arms and
phases options, and the memory probes made portable; the scheduling is the
same.
"""
import argparse, collections, ctypes, csv, glob, json, os, platform, re, signal, subprocess, sys, threading, time

HERE = os.path.dirname(os.path.abspath(__file__))
RUN, TRACE = HERE + "/run_arm.sh", HERE + "/run_trace.sh"
PAPER_ARMS = ["qisq2_low", "float_low", "float_l2", "limdd_qisq2", "limdd_float", "limdd_float_l2"]
PAPER_LAZY = ["limdd_qisq2_lazy", "limdd_lazy"]
EXACT = ("qisq2_low", "limdd_qisq2", "limdd_qisq2_lazy", "trace")   # exact weights: big memory
TRACES = ["20_700", "20_750", "30_700", "30_750", "40_700", "40_800", "50_700", "60_700", "70_700", "80_700"]
GB = 1 << 30
SWAPOUT_HOLD = 50000 * 16384                    # bytes swapped out in the last minute (50000 macOS pages)
DEFAULTS = dict(budget_gb=10, hard_gb=14, maxjobs=6, max_sys=0.15, min_cpu=0.85, heavy_gb=6)

ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("data"); ap.add_argument("out")
ap.add_argument("--dry", action="store_true")
ap.add_argument("--phases", default="A,T,B,C,D,E")
ap.add_argument("--arms", default=",".join(PAPER_ARMS))
ap.add_argument("--lazy-arms", default=",".join(PAPER_LAZY))
ap.add_argument("--only", default=None)
ap.add_argument("--timeout", default="1800")
ARGS = ap.parse_args()
DATA, OUTROOT = os.path.abspath(ARGS.data), os.path.abspath(ARGS.out)
ARMS = [a for a in ARGS.arms.split(",") if a]
LAZY = [a for a in ARGS.lazy_arms.split(",") if a]
PHASES = set(ARGS.phases.split(","))
ONLY = re.compile(ARGS.only) if ARGS.only else None
TO = ARGS.timeout
L = OUTROOT + "/sweep"
for v in ("QSY_W1", "QSY_W2"):
    if not os.access(os.environ.get(v, ""), os.X_OK) and not ARGS.dry:
        sys.exit(f"{v} must name a run_qasm_on_qmdd (run_arm.sh)")
ENV = dict(os.environ)
CONF = L + "/sweep.conf"                        # JSON, read again all the time, overrides these

def conf():
    c = dict(DEFAULTS)
    try: c.update(json.load(open(CONF)))
    except (OSError, ValueError): pass
    return c
OUT_A, OUT_T = OUTROOT + "/amend/out_clean", OUTROOT + "/amend/trace_clean"
OUT_C, OUT_D, OUT_E = OUTROOT + "/pilot/out_clean", OUTROOT + "/amend/out_clean_lazy", OUTROOT + "/pilot/out_clean_lazy"
DEAD = re.compile(r"\b(TIMEOUT|TAB_FULL|LIM_FULL|STAB_FULL|SKIPPED_PREDICTED)\b")
LOG = sys.stdout

def say(msg):
    print(f"{time.strftime('%H:%M:%S')} {msg}", file=LOG, flush=True)

# ---------------------------------------------------------------- earlier runs
def prev_info(base, dirs):
    """(status, wall) of the first earlier record, for the order and the estimate."""
    for d in dirs:
        st = os.path.join(d, base + ".status")
        if not os.path.exists(st): continue
        stxt = open(st).read()
        txt = stxt.split()
        status = txt[2] if len(txt) > 2 else "?"
        m = re.search(r"wall=(\d+)s", stxt)              # run_trace.sh puts it in the .status
        try: m = m or re.search(r"wall=(\d+)s", open(os.path.join(d, base + ".out")).read())
        except OSError: pass
        wall = int(m.group(1)) if m else 0
        if status == "TIMEOUT": wall = max(wall, 1800)
        return status, wall
    return None, 0

# earlier records in DATA (the stored clean sweep first, then older runs where present)
def existing(*ds): return [DATA + d for d in ds if os.path.isdir(DATA + d)]
PREV_AMEND = existing("/amend/out_clean", "/amend/out_rnd", "/amend/out", "/ldd/out")
PREV_LAZY = existing("/amend/out_clean_lazy")
PREV_PILOT = existing("/pilot/out_clean", "/pilot/out")
PREV_TRACE = existing("/amend/trace_clean", "/amend/trace")
OLD_AUDIT = DATA + "/sweep/clean_audit.csv"

def family(path):
    b = os.path.basename(path)
    if b.startswith("clifford_T_circuit"): return "random"
    if path.startswith(DATA + "/pilot/"): return os.path.basename(os.path.dirname(path))
    return re.match(r"[a-z]+", b).group(0)

def estimate(arm, fam, status, wall):
    if arm in ("float_low", "float_l2"): return 1 * GB
    if arm in ("limdd_float", "limdd_float_l2", "limdd_lazy"): return 2 * GB
    if arm == "trace": return 3 * GB
    if fam.startswith("synth"): return 3 * GB   # a ladder learns from its smaller instances (Ladder.current)
    hard = status in ("TIMEOUT", "TAB_FULL") or wall >= 600
    if arm in ("limdd_qisq2", "limdd_qisq2_lazy"): return 10 * GB if hard else 3 * GB
    return 6 * GB if hard else 3 * GB // 2

def learned_peaks():
    """(arm, circuit) -> the largest footprint an earlier attempt reached
    (clean_audit.csv: finished, rejected and stopped runs alike)."""
    out = {}
    for path in (OLD_AUDIT, AUDIT):
      try:
        for r in csv.DictReader(open(path)):
            try: p = float(r["peak_footprint_gb"]) * GB
            except (KeyError, ValueError): continue
            k = (r["arm"], r["circuit"])
            out[k] = max(out.get(k, 0), p)
      except OSError: pass
    return out
PEAKS = {}

# ---------------------------------------------------------------- jobs
class Job:
    def __init__(self, phase, arm, qasm, out, prev=(None, 0), ladder=None, trace=None):
        self.phase, self.arm, self.qasm, self.out, self.ladder, self.trace = phase, arm, qasm, out, ladder, trace
        self.base = trace if trace else os.path.basename(qasm)[:-5] + "_" + arm
        self.prev_status, self.prev_wall = prev
        self.est = estimate(arm, family(qasm), *prev)
        seen = PEAKS.get((arm, os.path.basename(qasm)[:-5]), 0)
        if seen: self.est = max(self.est, int(1.2 * seen))
        self.state = "wait"                       # wait, run, done
        self.root = self.binpid = None
        self.start = self.peak = self.peak_total = 0
        self.swap0 = 0
        self.killed = False
        self.tries = 0                           # rejected by the quality gate so far

    def status_file(self): return os.path.join(self.out, self.base + ".status")

    def cmd(self):
        if self.trace: return [TRACE, self.trace, DATA + "/ldd/qasm", self.out]
        return [RUN, self.arm, self.qasm, self.out, TO]

class Ladder:
    """One variant growing in size, one arm: the next instance runs once the
    previous one is done; after a dead one the rest is SKIPPED_PREDICTED."""
    def __init__(self, phase, key, arm, files, out, prevdirs, proxy=None):
        self.phase, self.key, self.arm, self.out = phase, key, arm, out
        self.files, self.i, self.dead, self.last_peak = files, 0, False, 0
        pa = proxy or arm
        self.prevs = [prev_info(os.path.basename(f)[:-5] + "_" + pa, prevdirs) for f in files]
        self.total = sum(w for _, w in self.prevs)
        self.job = None

    def current(self):
        """The job of the instance to run next, or None when the ladder is done."""
        while self.i < len(self.files):
            if self.job is None:
                self.job = Job(self.phase, self.arm, self.files[self.i], self.out, self.prevs[self.i], ladder=self)
                if self.last_peak: self.job.est = max(self.job.est, int(1.5 * self.last_peak))
            j = self.job
            if self.dead and not os.path.exists(j.status_file()):
                os.makedirs(self.out, exist_ok=True)
                base = os.path.basename(j.qasm)[:-5]
                with open(j.status_file(), "w") as fh: fh.write(f"{self.arm:<12} {base:<26} SKIPPED_PREDICTED\n")
            if os.path.exists(j.status_file()) and j.state != "run":
                if key(j) in BUSY: unmeasured(j); return j
                self.advance(j); continue
            return j
        return None

    def advance(self, j):
        try: line = open(j.status_file()).read()
        except OSError: line = ""
        if DEAD.search(line): self.dead = True
        if j.peak: self.last_peak = j.peak
        self.i += 1; self.job = None

def ladder_key(fam, name):
    nums = [int(x) for x in re.findall(r"\d+", name)]
    if fam == "mirror":
        if "clifford_T" in name: return "mirror-base", nums[-2:]
        return name.split("_")[0], nums
    if fam == "iqp_close":
        n, g = nums[-2:]; return f"iqp_close_g{g // n}n", nums[-2:]
    if fam == "grover_iter":
        n, k = nums[-2:]; return f"grover_iter_n{n}", [k]
    if fam.startswith("synth"):
        return f"{fam}_{name.rsplit('_', 1)[1]}", nums
    if fam == "qwalk":
        m, t = nums[-2:]; return f"qwalk_m{m}", [t]
    return fam, nums

def ladders(phase, arms, out, proxy=None):
    groups = {}
    for f in sorted(glob.glob(DATA + "/pilot/qasm/*/*.qasm")):
        if ONLY and not ONLY.search(os.path.basename(f)[:-5]): continue
        fam = os.path.basename(os.path.dirname(f)); key, size = ladder_key(fam, os.path.basename(f)[:-5])
        for arm in arms: groups.setdefault((key, arm), []).append((size, f))
    lad = [Ladder(phase, key, arm, [f for _, f in sorted(v)], out, PREV_PILOT, proxy and proxy[arm])
           for (key, arm), v in sorted(groups.items())]
    return sorted(lad, key=lambda x: -x.total)

def plan():
    paper = sorted(glob.glob(DATA + "/ldd/qasm/*.qasm")) + sorted(glob.glob(DATA + "/amend/qasm/*.qasm"))
    if ONLY: paper = [f for f in paper if ONLY.search(os.path.basename(f)[:-5])]
    figs = [f for f in paper if family(f) in ("random", "grover", "w", "wclif")]
    rest = [f for f in paper if f not in figs]
    def singles(phase, files, arms, out, prevdirs):
        js = [Job(phase, a, f, out, prev_info(os.path.basename(f)[:-5] + "_" + a, prevdirs)) for f in files for a in arms]
        return sorted(js, key=lambda j: -j.prev_wall)
    lazy79 = figs + [f for f in rest if family(f) == "hidden"]
    traces = [t for t in TRACES if not ONLY or ONLY.search(f"clifford_T_circuit_{t}")]
    q = []
    if "A" in PHASES: q += singles("A", figs, ARMS, OUT_A, PREV_AMEND)
    if "T" in PHASES:
        q += [Job("T", "trace", DATA + f"/ldd/qasm/clifford_T_circuit_{t}.qasm", OUT_T,
                  prev_info(t, PREV_TRACE), trace=t) for t in traces]
    if "B" in PHASES: q += singles("B", rest, ARMS, OUT_A, PREV_AMEND)
    if "C" in PHASES: q += ladders("C", ARMS, OUT_C)
    if "D" in PHASES: q += singles("D", lazy79, LAZY, OUT_D, PREV_LAZY)
    if "E" in PHASES:
        q += ladders("E", LAZY, OUT_E, proxy={"limdd_qisq2_lazy": "limdd_qisq2", "limdd_lazy": "limdd_float"})
    return q, len(figs), len(rest), len(lazy79)

# ---------------------------------------------------------------- the machine
SYSTEM = platform.system()
PAGE = os.sysconf("SC_PAGE_SIZE")

if SYSTEM == "Darwin":
    _libproc = ctypes.CDLL("/usr/lib/libproc.dylib")
    _libproc.proc_pid_rusage.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
    class _RI(ctypes.Structure):
        _fields_ = [("uuid", ctypes.c_uint8 * 16)] + [(n, ctypes.c_uint64) for n in (
            "user", "system", "pkg_idle", "interrupt", "pageins", "wired", "resident", "footprint",
            "start", "exit", "c_user", "c_system", "c_pkg_idle", "c_interrupt", "c_pageins",
            "c_elapsed", "dio_read", "dio_write")]

def footprint(pid):
    """Bytes: macOS phys_footprint (compressed memory included); Linux resident
    plus swapped; elsewhere 0 (not measured)."""
    if SYSTEM == "Darwin":
        ri = _RI()
        return ri.footprint if _libproc.proc_pid_rusage(pid, 2, ctypes.byref(ri)) == 0 else 0
    if SYSTEM == "Linux":
        kb = 0
        try:
            with open(f"/proc/{pid}/smaps_rollup") as fh:
                for line in fh:
                    if line.startswith(("Rss:", "Swap:")): kb += int(line.split()[1])
        except (OSError, ValueError): return 0
        return kb * 1024
    return 0

def swapouts():
    """Bytes swapped out since boot (0 where unknown)."""
    if SYSTEM == "Darwin":
        out = subprocess.run(["vm_stat"], capture_output=True, text=True).stdout
        m = re.search(r"page size of (\d+) bytes", out)
        page = int(m.group(1)) if m else PAGE
        for line in out.splitlines():
            if line.startswith("Swapouts:"): return int(line.split()[1].rstrip(".")) * page
        return 0
    if SYSTEM == "Linux":
        try:
            for line in open("/proc/vmstat"):
                if line.startswith("pswpout "): return int(line.split()[1]) * PAGE
        except OSError: pass
    return 0

MAXRSS_UNIT = 1 if SYSTEM == "Darwin" else 1024      # ru_maxrss: bytes on macOS, KiB on Linux

def proc_table():
    t = {}
    for line in subprocess.run(["ps", "-A", "-o", "pid=,ppid=,command="], capture_output=True, text=True).stdout.splitlines():
        p = line.split(None, 2)
        if len(p) == 3: t[int(p[0])] = (int(p[1]), p[2])
    return t

def descendants(t, root):
    kids = collections.defaultdict(list)
    for pid, (pp, _) in t.items(): kids[pp].append(pid)
    out, todo = [], [root]
    while todo:
        p = todo.pop(0); out.append(p); todo += kids.get(p, [])
    return out                                   # top-down

def is_sim(cmd): return cmd.split()[0].endswith("/run_qasm_on_qmdd")

# ---------------------------------------------------------------- the scheduler
lock = threading.Lock()
PAUSE = L + "/sweep.pause"                 # while it exists, nothing new starts
BUSY = set()                                     # keys of runs in flight from before
ALL = []                                         # every job object, for the estimates
running = []                                     # jobs of this sweep
foreign = {}                                     # pid -> footprint of other simulators
swaplog = collections.deque()                    # (time, swapouts)
last_start = 0.0
sampled = threading.Event()                      # the monitor has looked at the machine once
audit_lock = threading.Lock()
AUDIT = L + "/audit.csv"
AUDIT_COLS = ["phase", "arm", "circuit", "status", "start", "end", "wall_s", "sim_time_s", "cpu_s",
              "utime_s", "stime_s", "peak_footprint_gb", "maxrss_gb", "majflt", "peak_total_gb", "swapouts",
              "est_gb", "note"]

def set_aside(j, sub):
    """Move the files of j's run to <out>_<sub>/ (the run is then made again)."""
    d = j.out.rstrip("/") + "_" + sub
    os.makedirs(d, exist_ok=True)
    for ext in (".status", ".json", ".out", ".csv"):
        f = os.path.join(j.out, j.base + ext)
        if os.path.exists(f): os.replace(f, os.path.join(d, j.base + ext))

def audit(row):
    with audit_lock:
        new = not os.path.exists(AUDIT)
        with open(AUDIT, "a", newline="") as fh:
            w = csv.writer(fh)
            if new: w.writerow(AUDIT_COLS)
            w.writerow([row.get(c, "") for c in AUDIT_COLS])

def monitor():
    last_mem, over = 0, 0
    while True:
        t = proc_table()
        now = time.time()
        so = swapouts()
        with lock:
            swaplog.append((now, so))
            while swaplog and now - swaplog[0][0] > 60: swaplog.popleft()
            mine = set()
            for j in running:
                if j.root is None: continue
                d = descendants(t, j.root)
                mine.update(d)
                sims = [p for p in d if p in t and is_sim(t[p][1])]
                if sims: j.binpid = sims[0]
                if j.binpid: j.peak = max(j.peak, footprint(j.binpid))
            foreign.clear()
            for pid, (_, cmd) in t.items():
                if is_sim(cmd) and pid not in mine: foreign[pid] = footprint(pid)
            fp = {j: footprint(j.binpid) for j in running if j.binpid}
            total = sum(fp.values()) + sum(foreign.values())
            for j in running: j.peak_total = max(j.peak_total, total)
            victim = None
            HARD = int(conf()["hard_gb"] * GB)
            over = over + 1 if total > HARD else 0
            if over >= 2:                            # over HARD for two samples in a row
                mine_live = [j for j in running if j.root and not j.killed and fp.get(j, 0) > 0]
                if len(mine_live) > 1 or (mine_live and foreign):
                    # light runs before the heavy one, phase A last; among them the
                    # youngest whose footprint alone brings the total under HARD,
                    # else the largest
                    pool = ([j for j in mine_live if j.phase != "A" and not heavy(j)]
                            or [j for j in mine_live if j.phase != "A"] or mine_live)
                    enough = [j for j in pool if fp[j] >= total - HARD]
                    victim = max(enough, key=lambda j: j.start) if enough else max(pool, key=lambda j: fp[j])
                    victim.killed = True
                    over = 0
            if now - last_mem >= 60:
                last_mem = now
                with open(L + "/mem.csv", "a") as fh:
                    fh.write(f"{time.strftime('%H:%M:%S')},{len(running)},{len(foreign)},{total / GB:.2f},"
                             f"{committed() / GB:.2f},{so},{int(swapping())}\n")
        sampled.set()
        if victim:
            say(f"over {HARD / GB:.0f} GB ({total / GB:.1f}): stopping {victim.base} "
                f"({footprint(victim.binpid) / GB:.1f} GB), queued again")
            for pid in descendants(proc_table(), victim.root):
                try: os.kill(pid, signal.SIGTERM)
                except ProcessLookupError: pass
        time.sleep(5)

def heavy(j): return j.est >= conf()["heavy_gb"] * GB

def committed():
    """A light run counts with its estimate until it outgrows it; the heavy run
    (at most one) counts with what it uses now: when it grows, the backstop
    stops light runs first."""
    half = conf()["budget_gb"] * GB // 2
    return sum(max(footprint(j.binpid) if j.binpid else 0, min(j.est, half) if heavy(j) else j.est)
               for j in running) + sum(foreign.values())

def swapping():
    return len(swaplog) > 1 and swaplog[-1][1] - swaplog[0][1] > SWAPOUT_HOLD

def run_job(j):
    os.makedirs(j.out, exist_ok=True)
    p = subprocess.Popen(j.cmd(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=ENV,
                         start_new_session=True)
    with lock: j.root = p.pid
    _, wstat, ru = os.wait4(p.pid, 0)
    p.returncode = os.waitstatus_to_exitcode(wstat)
    end = time.time()
    so = swapouts()
    if j.killed:
        for ext in (".status", ".json", ".out", ".csv"):
            f = os.path.join(j.out, j.base + ext)
            if os.path.exists(f): os.remove(f)
        audit(dict(phase=j.phase, arm=j.arm, circuit=os.path.basename(j.qasm)[:-5], status="STOPPED",
                   start=time.strftime("%H:%M:%S", time.localtime(j.start)), end=time.strftime("%H:%M:%S"),
                   wall_s=f"{end - j.start:.0f}", peak_footprint_gb=f"{j.peak / GB:.2f}",
                   peak_total_gb=f"{j.peak_total / GB:.2f}", est_gb=f"{j.est / GB:.1f}", note="over HARD"))
        with lock:
            running.remove(j)
            j.est = max(int(1.5 * j.peak), 2 * j.est)
            j.state, j.killed, j.root, j.binpid, j.peak, j.peak_total = "wait", False, None, None, 0, 0
        return
    try: line = open(j.status_file()).read().strip()
    except OSError: line = f"{j.arm} {j.base} NOSTATUS rc={p.returncode}"
    f = line.split()
    simt = f[3] if len(f) > 3 and not j.trace else ""
    status = f[2] if len(f) > 2 else "?"            # a trace: rc=<exit code of the simulator>
    peak = j.peak
    # the quality gate: a single-threaded run that was not slowed down has CPU time
    # close to its wall time and little system time (decompressing a page counts
    # as system time of the process that touches it)
    c = conf()
    cpu, wall = ru.ru_utime + ru.ru_stime, end - j.start
    why = ""
    if wall > 20:
        if ru.ru_stime > c["max_sys"] * cpu: why = f"sys {ru.ru_stime:.0f} s of cpu {cpu:.0f} s"
        elif cpu < c["min_cpu"] * wall: why = f"cpu {cpu:.0f} s of wall {wall:.0f} s"
    retry = bool(why) and j.tries == 0
    audit(dict(phase=j.phase, arm=j.arm, circuit=os.path.basename(j.qasm)[:-5], status=status,
               start=time.strftime("%H:%M:%S", time.localtime(j.start)), end=time.strftime("%H:%M:%S"),
               wall_s=f"{wall:.1f}", sim_time_s=simt, cpu_s=f"{cpu:.1f}", utime_s=f"{ru.ru_utime:.1f}",
               stime_s=f"{ru.ru_stime:.1f}", peak_footprint_gb=f"{j.peak / GB:.2f}",
               maxrss_gb=f"{ru.ru_maxrss * MAXRSS_UNIT / GB:.2f}", majflt=ru.ru_majflt, peak_total_gb=f"{j.peak_total / GB:.2f}",
               swapouts=so - j.swap0, est_gb=f"{j.est / GB:.1f}",
               note=("rejected, run again: " if retry else "kept, second try: " if why else "") + why))
    if retry:
        set_aside(j, "rejected")
        say(f"[{j.phase}] {line}  REJECTED ({why}), queued again")
        with lock:
            running.remove(j)
            j.tries += 1
            j.est = max(j.est, int(1.5 * j.peak))
            j.state, j.root, j.binpid, j.peak, j.peak_total = "wait", None, None, 0, 0
        return
    say(f"[{j.phase}] {line}" + (f"  (unclean twice: {why})" if why else ""))
    with lock:
        running.remove(j)
        j.state = "done"
        if j.ladder: j.ladder.advance(j)
        if j.arm in EXACT and peak:                 # the other exact arms of this circuit
            for k in ALL:
                if k.state == "wait" and k.qasm == j.qasm and k.arm in EXACT:
                    k.est = max(k.est, int(1.2 * peak))

def key(j): return ("trace", j.trace, j.out) if j.trace else (j.arm, j.qasm, j.out)

def unmeasured(j):
    """The run of j by an earlier instance of this script has ended: it was not
    measured here (no CPU time), so its files are set aside and j runs again."""
    f = open(j.status_file()).read().split()
    set_aside(j, "unmeasured")
    BUSY.discard(key(j))
    audit(dict(phase=j.phase, arm=j.arm, circuit=os.path.basename(j.qasm)[:-5],
               status=f[2] if len(f) > 2 else "?", end=time.strftime("%H:%M:%S"),
               sim_time_s=f[3] if len(f) > 3 and not j.trace else "",
               note="run of an earlier instance, not measured: set aside, run again"))
    say(f"[{j.phase}] {' '.join(f)}  (earlier instance, not measured: set aside, queued again)")

def waiting(queue, phases=None, blocked=False):
    """The waiting jobs in priority order (a ladder offers its current instance).
    A job whose run is in flight from an earlier instance of this script waits
    for its .status, and is only listed with blocked=True."""
    for e in queue:
        if isinstance(e, Ladder):
            if e.job is not None and e.job.state == "run": continue
            j = e.current()
            if j is None or j.state != "wait": continue
        else:
            j = e
            if j.state != "wait": continue
            if os.path.exists(j.status_file()):
                if key(j) not in BUSY:
                    j.state = "done"; continue
                unmeasured(j)                        # and it waits to run again
        if phases and j.phase not in phases: continue
        if key(j) in BUSY and not blocked: continue
        yield j

def inflight():
    """(arm, qasm) of run_arm.sh runs already on the machine, e.g. of an earlier
    instance of this script, and the trace ids of run_trace.sh runs."""
    out = set()
    for _, cmd in proc_table().values():
        t = cmd.split()
        for i, x in enumerate(t):
            if x.endswith("/run_arm.sh") and len(t) > i + 3: out.add((t[i + 1], t[i + 2], t[i + 3]))
            if x.endswith("/run_trace.sh") and len(t) > i + 3: out.add(("trace", t[i + 1], t[i + 3]))
    return out

def link_inputs():
    """OUT gets links to DATA's circuits and reference probabilities, so that
    regen_amendment.sh can take OUT as its ROOT."""
    for d in ("ldd/qasm", "amend/qasm", "pilot/qasm", "pilot/reference_probs.csv"):
        src, dst = os.path.join(DATA, d), os.path.join(OUTROOT, d)
        if os.path.exists(src) and not os.path.lexists(dst) and os.path.realpath(src) != os.path.realpath(dst):
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            os.symlink(src, dst)

def main():
    global last_start
    PEAKS.update(learned_peaks())
    queue, nf, nr, nl = plan()
    jobs = lambda: [e for e in queue if isinstance(e, Job)]
    lads = [e for e in queue if isinstance(e, Ladder)]
    count = collections.Counter(j.phase for j in jobs())
    for e in lads: count[e.phase] += len(e.files)
    say(f"plan: {nf} figure circuits, {nr} others, {nl} for the lazy arms; runs per phase "
        + ", ".join(f"{k}={v}" for k, v in sorted(count.items())))
    slot = collections.Counter()
    for j in jobs(): slot[j.phase] += j.prev_wall
    for e in lads: slot[e.phase] += e.total
    say("earlier wall time per phase (h): " + ", ".join(f"{k}={v / 3600:.1f}" for k, v in sorted(slot.items())))
    if ARGS.dry:
        for j in jobs()[:12]: say(f"  {j.phase} {j.base:45s} prev {j.prev_status} {j.prev_wall}s est {j.est / GB:.1f} GB")
        for e in lads[:6]: say(f"  {e.phase} ladder {e.key} {e.arm}: {len(e.files)} instances, prev {e.total}s")
        return
    os.makedirs(L, exist_ok=True)
    link_inputs()
    if SYSTEM == "Darwin":                       # no sleep while the sweep runs
        subprocess.Popen(["caffeinate", "-i", "-s", "-w", str(os.getpid())])
    if os.path.exists(AUDIT) and open(AUDIT).readline().strip() != ",".join(AUDIT_COLS):
        os.replace(AUDIT, AUDIT + ".old")            # written with other columns
    say(f"settings {conf()} (from {CONF} where it exists)")
    BUSY.update(inflight())
    ALL.extend(jobs())
    threading.Thread(target=monitor, daemon=True).start()
    sampled.wait()                               # no start before the first sample
    phase_a_said = False
    last_busy = time.time()
    while True:
        if BUSY and time.time() - last_busy > 30:    # runs of an earlier instance that were stopped
            last_busy = time.time()
            gone = BUSY - inflight()
            with lock:
                for k in gone:
                    if not os.path.exists(k[2] + "/" + (k[1] if k[0] == "trace" else
                                          os.path.basename(k[1])[:-5] + "_" + k[0]) + ".status"):
                        BUSY.discard(k)
        with lock:
            for e in lads:                           # a ladder's job objects, for the estimates
                if e.job is not None and e.job not in ALL: ALL.append(e.job)
            pending_a = any(True for _ in waiting(queue, {"A"}, blocked=True))
            if not phase_a_said and not pending_a and not any(j.phase == "A" for j in running):
                say("phase A done"); phase_a_said = True
            nwait = sum(1 for _ in waiting(queue, blocked=True))
            if nwait == 0 and not running:
                break
            cf = conf()
            free = cf["maxjobs"] - len(running) - len(foreign)
            ready = (free > 0 and time.time() - last_start >= 3 and not swapping()
                     and not os.path.exists(PAUSE))
            pick = None
            if ready:
                c = committed()
                alone = not running and not foreign
                busy_heavy = any(heavy(k) for k in running)
                for j in waiting(queue, {"A"} if pending_a else None):
                    if heavy(j):                     # one heavy run at a time
                        if busy_heavy: continue
                        if alone or c + min(j.est, cf["budget_gb"] * GB // 2) <= cf["budget_gb"] * GB:
                            pick = j; break
                    elif alone or c + j.est <= cf["budget_gb"] * GB: pick = j; break
            if pick:
                pick.state, pick.start, pick.swap0 = "run", time.time(), swapouts()
                running.append(pick)
                last_start = pick.start
        if pick:
            threading.Thread(target=run_job, args=(pick,), daemon=True).start()
        time.sleep(2)
    say("end")

if __name__ == "__main__":
    main()
