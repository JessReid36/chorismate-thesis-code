#!/usr/bin/env python3
"""
sp_stage.py - the in-job checks of the stationary-point barriers (results 05_qmmm/24_sp_barriers/SP_PLAN.txt). Called
by the job scripts on hpc1 (job_*.pbs, written by sp_prepare.py) before and after each ORCA step. Python 3.6, standard
library; deterministic. Every check writes RD/STEP/STAGE_RESULT.txt (key=value lines; status=PASS or FAIL) and exits 1
on FAIL, so the job stops and the jobs that depend on it (PBS afterok) never start.

  prep RD STEP            before a step: its start.xyz. ts, r0: the generator's file, checked against starts.sha256.
                          ts_freq: the TS's last geometry (ts/ts.xyz); irc: ts_freq's start.xyz and freq.hess;
                          r_irc, p_irc: written by check_irc (checked against its record); r_irc_freq, p_irc_freq,
                          r0_freq: the minimisation's last geometry (opt.xyz). Each source step must have passed.
  check_opt RD STEP KIND  after an optimisation (KIND ts or min): ORCA terminated normally and printed "THE OPTIMIZATION
                          HAS CONVERGED"; the last geometry holds the full system with the start's elements; no atom
                          outside region.txt moved (2e-6 A); the substrate's breaking O-C and forming C-C distances
                          (r_irc and r0 must be reactant-like, p_irc product-like); final QM/MM energy; cycles; sizes.
                          Writes final_region.xyz (the region atoms, 5th column the 0-based index).
  check_freq RD STEP KIND after the partial-Hessian frequencies (24 substrate atoms displaced): counted imaginary modes
                          (|nu| >= 20 cm-1): exactly 1 for KIND ts, 0 for KIND min; smaller imaginary values listed.
  check_irc RD            after the IRC: ORCA terminated normally; both end files hold the full system; no atom other
                          than the substrate moved; one end reactant-like and the other product-like; writes
                          irc/r_irc_start.xyz and irc/p_irc_start.xyz (verbatim copies of those ends), which prep
                          copies to r_irc/start.xyz and p_irc/start.xyz.
Reactant-like: breaking O-C (6214-6215) <= 1.80 A and forming C-C (6207-6219) >= 2.40 A; product-like: forming C-C
<= 1.80 A and breaking O-C >= 2.40 A (0-based atom indices, as in ORCA).
"""
import hashlib, os, re, shutil, sys

QM0, NQM, NATOMS = 6207, 24, 55680
QMS = list(range(QM0, QM0 + NQM))
BRK, FRM = (6214, 6215), (6207, 6219)
R_BRK_MAX, R_FRM_MIN = 1.80, 2.40
IMAG_MIN = 20.0           # cm-1: imaginary modes at or above this magnitude are counted
MOVE_TOL = 2e-6           # A
BASE = {"ts": "ts", "ts_freq": "freq", "irc": "irc", "r_irc": "opt", "p_irc": "opt", "r0": "opt",
        "r_irc_freq": "freq", "p_irc_freq": "freq", "r0_freq": "freq"}

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
def result(rd, step, status, items, reason=""):
    lines = ["status=" + status] + (["reason=" + reason] if reason else []) + ["%s=%s" % kv for kv in items]
    p = os.path.join(rd, step, "STAGE_RESULT.txt")
    if os.path.exists(p):      # a step's record is written once; a second job on the same step changes nothing
        print("%s %s (not recorded: STAGE_RESULT.txt exists already): %s" % (step, status, reason)); sys.exit(1)
    with open(p, "w") as fh: fh.write("\n".join(lines) + "\n")
    print("%s %s: %s%s" % (step, status, reason if reason else "", "" if not items else
                           ("; " if reason else "") + ", ".join("%s=%s" % kv for kv in items)))
    sys.exit(0 if status == "PASS" else 1)
def fail(rd, step, reason, items=()): result(rd, step, "FAIL", list(items), reason)
def read_result(rd, step):
    p = os.path.join(rd, step, "STAGE_RESULT.txt")
    if not os.path.isfile(p): return {}
    return dict(l.rstrip("\n").split("=", 1) for l in open(p) if "=" in l)
def read_xyz(p):
    L = open(p).read().splitlines()
    n = int(L[0].split()[0]); R = [l.split()[:4] for l in L[2:2 + n]]
    if n != NATOMS or len(R) != NATOMS or any(len(r) < 4 for r in R): return None
    return R
def fl(r): return (float(r[1]), float(r[2]), float(r[3]))
def dist(a, b): return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5
def rc(R):
    X = lambda i: fl(R[i])
    return dist(X(BRK[0]), X(BRK[1])), dist(X(FRM[0]), X(FRM[1]))
def kind_of(b, f):
    if b <= R_BRK_MAX and f >= R_FRM_MIN: return "reactant"
    if f <= R_BRK_MAX and b >= R_FRM_MIN: return "product"
    return "neither"
def region(rd): return [int(x) for x in open(os.path.join(rd, "region.txt")).read().split()]
def moved(R, S, keep):
    """largest coordinate change of an atom not in keep, and whether every element matches"""
    m = 0.0; same = True
    for i in range(NATOMS):
        if R[i][0].lower() != S[i][0].lower(): same = False
        if i in keep: continue
        for k in (1, 2, 3):
            d = abs(float(R[i][k]) - float(S[i][k]))
            if d > m: m = d
    return m, same
def run_seconds(t):
    m = re.search(r"TOTAL RUN TIME: (\d+) days (\d+) hours (\d+) minutes (\d+) seconds", t)
    return "" if not m else str(((int(m.group(1)) * 24 + int(m.group(2))) * 60 + int(m.group(3))) * 60 + int(m.group(4)))
def energies(t): return re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", t)
def output(rd, step, name):
    p = os.path.join(rd, step, name)
    if not os.path.isfile(p): fail(rd, step, "no %s (ORCA did not start or did not write it)" % name)
    return open(p, errors="replace").read()

def prep(rd, step):
    d = os.path.join(rd, step); st = os.path.join(d, "start.xyz")
    def need(src_step):
        r = read_result(rd, src_step)
        if r.get("status") != "PASS": fail(rd, step, "%s has not passed its checks" % src_step)
        return r
    if step in ("ts", "r0"):
        want = dict((l.split()[1], l.split()[0]) for l in open(os.path.join(rd, "starts.sha256")) if l.strip())
        if not os.path.isfile(st) or sha(st) != want.get(step + "/start.xyz"):
            fail(rd, step, "start.xyz is missing or is not the generator's (starts.sha256)")
        print("%s: start.xyz is the generator's" % step); return
    if os.path.exists(st): fail(rd, step, "start.xyz exists already")
    if step == "ts_freq":
        r = need("ts"); shutil.copyfile(os.path.join(rd, "ts", "ts.xyz"), st)
        if sha(st) != r.get("xyz_sha256"): fail(rd, step, "ts/ts.xyz is not the geometry check_opt recorded")
    elif step == "irc":
        r = need("ts_freq"); shutil.copyfile(os.path.join(rd, "ts_freq", "start.xyz"), st)
        if sha(st) != r.get("start_sha256"): fail(rd, step, "ts_freq/start.xyz changed after its check")
        h = os.path.join(rd, "ts_freq", "freq.hess")
        if not os.path.isfile(h) or sha(h) != r.get("hess_sha256"): fail(rd, step, "ts_freq/freq.hess is missing or changed")
        shutil.copyfile(h, os.path.join(d, "freq.hess"))
    elif step in ("r_irc", "p_irc"):
        r = need("irc")
        p = os.path.join(rd, "irc", step + "_start.xyz")
        if not os.path.isfile(p) or sha(p) != r.get(step + "_start_sha256"): fail(rd, step, "the IRC end recorded by check_irc is missing or changed")
        shutil.copyfile(p, st)
    elif step in ("r_irc_freq", "p_irc_freq", "r0_freq"):
        src_step = step[:-5]; r = need(src_step)
        shutil.copyfile(os.path.join(rd, src_step, "opt.xyz"), st)
        if sha(st) != r.get("xyz_sha256"): fail(rd, step, "%s/opt.xyz is not the geometry check_opt recorded" % src_step)
    else:
        fail(rd, step, "unknown step")
    print("%s: start.xyz written (sha256 %s)" % (step, sha(st)))

def check_opt(rd, step, kind):
    d = os.path.join(rd, step); b = BASE[step]
    t = output(rd, step, b + ".out")
    it = [("cycles", str(t.count("GEOMETRY OPTIMIZATION CYCLE"))), ("run_seconds", run_seconds(t))]
    m = re.search(r"^Number of atoms\s+\.\.\.\.\s+(\d+)", t, re.M)
    it += [("region_atoms", str(len(region(rd)))), ("orca_optimiser_atoms", m.group(1) if m else "")]
    m = re.search(r"Geometry relaxation\s+\.\.\.\s+(\d+\.\d+) sec", t); it.append(("geometry_seconds", m.group(1) if m else ""))
    E = energies(t)
    it.append(("energy_count", str(len(E))))
    if "ORCA TERMINATED NORMALLY" not in t: fail(rd, step, "ORCA did not terminate normally", it)
    if "THE OPTIMIZATION HAS CONVERGED" not in t: fail(rd, step, "the optimisation did not converge", it)
    if not E: fail(rd, step, "no FINAL SINGLE POINT ENERGY (QM/MM)", it)
    it.append(("energy_Eh", E[-1]))
    px = os.path.join(d, b + ".xyz")
    if not os.path.isfile(px): fail(rd, step, "no %s.xyz" % b, it)
    R = read_xyz(px); S = read_xyz(os.path.join(d, "start.xyz"))
    if R is None or S is None: fail(rd, step, "%s.xyz or start.xyz is not the full system" % b, it)
    reg = region(rd); RS = set(reg)
    mv, same = moved(R, S, RS)
    it.append(("outside_region_max_move_A", "%.1e" % mv))
    if not same: fail(rd, step, "elements differ between start.xyz and %s.xyz" % b, it)
    if mv > MOVE_TOL: fail(rd, step, "an atom outside the region moved", it)
    br, fo = rc(R); k = kind_of(br, fo)
    it += [("break_A", "%.3f" % br), ("form_A", "%.3f" % fo), ("substrate", k)]
    want = {"r_irc": "reactant", "r0": "reactant", "p_irc": "product"}.get(step)
    if want and k != want: fail(rd, step, "the optimised substrate is %s, expected %s" % (k, want), it)
    with open(os.path.join(d, "final_region.xyz"), "w") as fh:
        fh.write("%d\n%s/%s: last geometry (%s.xyz), region atoms; 5th column: 0-based atom index\n" % (
            len(reg), os.path.basename(os.path.abspath(rd)), step, b)
            + "".join("%s %s %s %s %d\n" % (R[i][0], R[i][1], R[i][2], R[i][3], i) for i in reg))
    it += [("xyz_sha256", sha(px)), ("start_sha256", sha(os.path.join(d, "start.xyz")))]
    result(rd, step, "PASS", it)

def check_freq(rd, step, kind):
    d = os.path.join(rd, step)
    t = output(rd, step, "freq.out")
    it = [("run_seconds", run_seconds(t))]
    m = re.search(r"Number of displaced atoms\s+\.\.\.\s+(\d+)", t); it.append(("displaced_atoms", m.group(1) if m else ""))
    E = energies(t); it.append(("reference_energy_Eh", E[0] if E else ""))
    if "ORCA TERMINATED NORMALLY" not in t: fail(rd, step, "ORCA did not terminate normally", it)
    if not m or int(m.group(1)) != NQM: fail(rd, step, "not a partial Hessian over the 24 substrate atoms", it)
    k = t.rfind("VIBRATIONAL FREQUENCIES")
    if k < 0: fail(rd, step, "no VIBRATIONAL FREQUENCIES block", it)
    F = []
    for l in t[k:].splitlines()[1:]:
        mm = re.match(r"^\s*(\d+):\s+(-?\d+\.\d+) cm\*\*-1", l)
        if mm: F.append(float(mm.group(2)))
        elif F and not l.strip(): break
    it.append(("modes", str(len(F))))
    cnt = sorted(v for v in F if v < 0 and -v >= IMAG_MIN); small = sorted(v for v in F if v < 0 and -v < IMAG_MIN)
    it += [("imaginary_counted", " ".join("%.2f" % v for v in cnt) or "none"),
           ("imaginary_below_20", " ".join("%.2f" % v for v in small) or "none"),
           ("lowest_real_above_20_cm-1", "%.2f" % min([v for v in F if v >= IMAG_MIN] or [0.0]))]
    h = os.path.join(d, "freq.hess")
    if not os.path.isfile(h): fail(rd, step, "no freq.hess", it)
    it += [("start_sha256", sha(os.path.join(d, "start.xyz"))), ("hess_sha256", sha(h))]
    want = 1 if kind == "ts" else 0
    if len(cnt) != want: fail(rd, step, "%d counted imaginary modes, expected %d" % (len(cnt), want), it)
    result(rd, step, "PASS", it)

def check_irc(rd):
    step = "irc"; d = os.path.join(rd, step)
    t = output(rd, step, "irc.out")
    it = [("run_seconds", run_seconds(t))]
    fw, bw = t.find("FORWARD IRC"), t.find("BACKWARD IRC")
    def conv(a, z): return "yes" if a >= 0 and "THE IRC HAS CONVERGED" in t[a:z if z > a else len(t)] else "no"
    it += [("forward_converged", conv(fw, bw)), ("backward_converged", conv(bw, len(t)))]
    k = t.rfind("IRC PATH SUMMARY")
    rows = re.findall(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+\d+\.\d+\s+\d+\.\d+(.*)$", t[k:], re.M) if k >= 0 else []
    if rows:
        ts_rows = [r for r in rows if "TS" in r[3]]
        it += [("path_points", str(len(rows))), ("E_path_first_Eh", rows[0][1]), ("E_path_last_Eh", rows[-1][1]),
               ("E_ts_Eh", ts_rows[0][1] if ts_rows else "")]
    it.append(("hessian_atom_warning", "yes" if "corresponding active atoms in Hessian file" in t else "no"))
    if "ORCA TERMINATED NORMALLY" not in t: fail(rd, step, "ORCA did not terminate normally", it)
    S = read_xyz(os.path.join(d, "start.xyz")); QS = set(QMS); ends = {}
    for e in ("F", "B"):
        p = os.path.join(d, "irc_IRC_%s.xyz" % e)
        if not os.path.isfile(p): fail(rd, step, "no irc_IRC_%s.xyz" % e, it)
        R = read_xyz(p)
        if R is None: fail(rd, step, "irc_IRC_%s.xyz is not the full system" % e, it)
        mv, same = moved(R, S, QS)
        if not same or mv > MOVE_TOL: fail(rd, step, "irc_IRC_%s.xyz: elements differ or a non-substrate atom moved (%.1e A)" % (e, mv), it)
        br, fo = rc(R); ends[e] = (kind_of(br, fo), br, fo, p)
        it += [("end_%s" % e, "%s break %.3f form %.3f" % (ends[e][0], br, fo))]
    kinds = sorted(v[0] for v in ends.values())
    if kinds != ["product", "reactant"]: fail(rd, step, "the IRC ends are %s, not one reactant and one product" % " and ".join(kinds), it)
    for e, v in sorted(ends.items()):
        tgt = "r_irc" if v[0] == "reactant" else "p_irc"
        dst = os.path.join(d, tgt + "_start.xyz")
        shutil.copyfile(v[3], dst)
        it += [(tgt + "_from", "irc_IRC_%s.xyz" % e), (tgt + "_start_sha256", sha(dst))]
    result(rd, step, "PASS", it)

if __name__ == "__main__":
    a = sys.argv[1:]
    try:
        if len(a) == 3 and a[0] == "prep": prep(a[1], a[2])
        elif len(a) == 4 and a[0] == "check_opt" and a[3] in ("ts", "min"): check_opt(a[1], a[2], a[3])
        elif len(a) == 4 and a[0] == "check_freq" and a[3] in ("ts", "min"): check_freq(a[1], a[2], a[3])
        elif len(a) == 2 and a[0] == "check_irc": check_irc(a[1])
        else: sys.exit(__doc__)
    except Exception as e:      # an unreadable or truncated ORCA file: recorded as a failed check, never as "not run"
        fail(a[1], "irc" if a[0] == "check_irc" else a[2], "check could not read its files: %s: %s" % (type(e).__name__, e))
