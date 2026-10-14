#!/usr/bin/env python3
"""
spa_stage.py - the in-job checks of the stationary-point barriers under SP Amendment 1, the alternating scheme (results
05_qmmm/24_sp_barriers/SP_AMENDMENT1.txt). Called by the job scripts on hpc1 (job_*.pbs, written by spa_prepare.py)
before and after each ORCA step. Python 3.6, standard library; deterministic. Every check writes STAGE_RESULT.txt in
its folder (key=value lines; status=PASS or FAIL; never overwritten) and exits 1 on FAIL, so the job stops and the jobs
that depend on it (PBS afterok) never start. prep and alt_prep write no record: when they refuse (a folder or start
file exists already, a start file is missing or not the recorded one, the source step has not passed) the job stops
as a machinery stop (SP_AMENDMENT1.txt, A1.5).

Stationary points ts, r_irc, p_irc, r0 come from the alternating scheme: cycle K = environment step (folder
STEP/cKK_env: L-Opt of the region minus the 24 substrate atoms, substrate fixed) then core step (STEP/cKK_core: the 24
substrate atoms alone, environment fixed; OptTS with an exact substrate Hessian for ts, Opt for the minima).

  prep RD STEP              start.xyz of a step. ts, r0: the generator's file (starts.sha256). ts_freq: ts/final.xyz.
                            irc: ts_freq's start.xyz and freq.hess. r_irc, p_irc: the IRC end check_irc recorded.
                            r_irc_freq, p_irc_freq, r0_freq: STEP/final.xyz. The source step must have passed.
  alt_prep RD STEP K PART   folder STEP/cKK_PART (PART env or core): its input (STEP/PART.inp) and start.xyz - for env
                            in cycle 1 STEP/start.xyz, later the previous core step's core.xyz; for core this cycle's
                            env.xyz - each checked against the record of the step that wrote it
  alt_check RD STEP K PART  env: ORCA terminated normally, "The minimization has converged" with its printed largest
                            and RMS gradients at or below their thresholds, energy over the last 50 optimiser steps
                            fell by at most 2.5e-4 Eh (not tested below 51 steps), the substrate and every atom outside
                            the region unmoved; also records the gradient where the step began (Zhang's criterion,
                            reported only). core: ORCA terminated normally, "THE OPTIMIZATION HAS CONVERGED", every atom
                            other than the substrate unmoved, the substrate still TS-like (neither reactant- nor
                            product-like) for ts, reactant-like for r_irc and r0, product-like for p_irc
  alt_decide RD STEP K      after cycle K >= 2: converged when the environment step of cycle K changed the energy by at
                            most 0.10 kcal/mol, the core energy changed by at most 0.10 kcal/mol from cycle K-1, and no
                            substrate atom moved more than 0.02 A between the two cycles' core geometries. Then
                            STEP/final.xyz = cycle K's core.xyz and STEP/STAGE_RESULT.txt is PASS (exit 0). Not yet: exit
                            3. Not converged after 8 cycles: FAIL (exit 1)
  alt_fail RD STEP REASON   records a failed cycle as the step's FAIL
  check_freq RD STEP KIND   partial-Hessian frequencies (24 substrate atoms displaced): counted imaginary modes
                            (|nu| >= 20 cm-1): exactly 1 for KIND ts, 0 for KIND min
  check_irc RD              IRC: ORCA terminated normally; both end files hold the full system; nothing but the
                            substrate moved; one end reactant-like and the other product-like; writes
                            irc/r_irc_start.xyz and irc/p_irc_start.xyz
Reactant-like: breaking O-C (6214-6215) <= 1.80 A and forming C-C (6207-6219) >= 2.40 A; product-like: C-C <= 1.80 A
and O-C >= 2.40 A (0-based atom indices, as in ORCA). Energies: kcal/mol = Eh x 627.5094740631.
"""
import hashlib, os, re, shutil, sys

QM0, NQM, NATOMS = 6207, 24, 55680
QMS = list(range(QM0, QM0 + NQM))
BRK, FRM = (6214, 6215), (6207, 6219)
R_BRK_MAX, R_FRM_MIN = 1.80, 2.40
IMAG_MIN = 20.0              # cm-1
MOVE_TOL = 2e-6              # A
H = 627.5094740631
KMAX = 8                     # cycles (Zhang, Liu & Yang 2000: "no more than eight cycles")
ALT_DE = 0.10                # kcal/mol
ALT_DX = 0.02                # A
PLATEAU_EH, PLATEAU_N = 2.5e-4, 50
GMAX_KJ, GRMS_KJ = 0.4961, 0.1488   # kJ/mol/A: TolMaxG 1e-4 and TolRMSG 3e-5 Eh/bohr as ORCA's L-Opt prints them
ZHANG_RMS_KJ = 0.4184               # kJ/mol/A: Zhang, Liu & Yang 2000's MM criterion, RMS gradient 0.1 kcal/mol/A (reported only)
ALT_STEPS = ("ts", "r_irc", "p_irc", "r0")
WANT = {"ts": "neither", "r_irc": "reactant", "r0": "reactant", "p_irc": "product"}

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
def result(rd, step, status, items, reason="", code=None):
    lines = ["status=" + status] + (["reason=" + reason] if reason else []) + ["%s=%s" % kv for kv in items]
    p = os.path.join(rd, step, "STAGE_RESULT.txt")
    if os.path.exists(p):      # a record is written once; a second job on the same step changes nothing
        print("%s %s (not recorded: STAGE_RESULT.txt exists already): %s" % (step, status, reason)); sys.exit(1)
    with open(p, "w") as fh: fh.write("\n".join(lines) + "\n")
    print("%s %s: %s%s" % (step, status, reason, "" if not items else ("; " if reason else "") +
                           ", ".join("%s=%s" % kv for kv in items if not kv[0].endswith("sha256"))))
    sys.exit(code if code is not None else (0 if status == "PASS" else 1))
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
def rc(R): return dist(fl(R[BRK[0]]), fl(R[BRK[1]])), dist(fl(R[FRM[0]]), fl(R[FRM[1]]))
def kind_of(b, f):
    if b <= R_BRK_MAX and f >= R_FRM_MIN: return "reactant"
    if f <= R_BRK_MAX and b >= R_FRM_MIN: return "product"
    return "neither"
def region(rd): return [int(x) for x in open(os.path.join(rd, "region.txt")).read().split()]
def moved(R, S, keep):
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
def cyc(k, part): return "c%02d_%s" % (k, part)

def stop(step, why): print("SP_%s_FAIL machinery stop (prep): %s; no record written" % (step, why)); sys.exit(1)
def prep(rd, step):
    d = os.path.join(rd, step); st = os.path.join(d, "start.xyz")
    fail = lambda rd_, step_, why: stop(step_, why)          # prep's refusals are machinery stops, never step records
    def need(src_step):
        r = read_result(rd, src_step)
        if r.get("status") != "PASS": fail(rd, step, "%s has not passed its checks" % src_step)
        return r
    if step in ("ts", "r0"):
        want = dict((l.split()[1], l.split()[0]) for l in open(os.path.join(rd, "starts.sha256")) if l.strip())
        if not os.path.isfile(st) or sha(st) != want.get(step + "/start.xyz"):
            fail(rd, step, "start.xyz is missing or is not the generator's (starts.sha256)")
    else:
        if os.path.exists(st): fail(rd, step, "start.xyz exists already")
        if step == "ts_freq":
            r = need("ts"); shutil.copyfile(os.path.join(rd, "ts", "final.xyz"), st)
            if sha(st) != r.get("xyz_sha256"): fail(rd, step, "ts/final.xyz is not the geometry the alternation recorded")
        elif step == "irc":
            r = need("ts_freq"); shutil.copyfile(os.path.join(rd, "ts_freq", "start.xyz"), st)
            if sha(st) != r.get("start_sha256"): fail(rd, step, "ts_freq/start.xyz changed after its check")
            h = os.path.join(rd, "ts_freq", "freq.hess")
            if not os.path.isfile(h) or sha(h) != r.get("hess_sha256"): fail(rd, step, "ts_freq/freq.hess is missing or changed")
            shutil.copyfile(h, os.path.join(d, "freq.hess"))
        elif step in ("r_irc", "p_irc"):
            r = need("irc"); p = os.path.join(rd, "irc", step + "_start.xyz")
            if not os.path.isfile(p) or sha(p) != r.get(step + "_start_sha256"): fail(rd, step, "the IRC end recorded by check_irc is missing or changed")
            shutil.copyfile(p, st)
        elif step in ("r_irc_freq", "p_irc_freq", "r0_freq"):
            src_step = step[:-5]; r = need(src_step)
            shutil.copyfile(os.path.join(rd, src_step, "final.xyz"), st)
            if sha(st) != r.get("xyz_sha256"): fail(rd, step, "%s/final.xyz is not the geometry the alternation recorded" % src_step)
        else:
            fail(rd, step, "unknown step")
    with open(os.path.join(d, "start.sha256"), "w") as fh: fh.write(sha(st) + "\n")
    print("%s: start.xyz ready (sha256 %s)" % (step, sha(st)))

def alt_prep(rd, step, k, part):
    if step not in ALT_STEPS or part not in ("env", "core") or not 1 <= k <= KMAX: sys.exit("alt_prep: bad arguments")
    sub = os.path.join(step, cyc(k, part)); d = os.path.join(rd, sub)
    if os.path.exists(d): print("alt_prep: %s exists already" % sub); sys.exit(1)
    if part == "env" and k == 1:
        src = os.path.join(rd, step, "start.xyz"); want = open(os.path.join(rd, step, "start.sha256")).read().split()[0]
    else:
        prev = os.path.join(step, cyc(k - 1, "core") if part == "env" else cyc(k, "env"))
        r = read_result(rd, prev)
        if r.get("status") != "PASS": print("alt_prep: %s has not passed its checks" % prev); sys.exit(1)
        src = os.path.join(rd, prev, ("core" if part == "env" else "env") + ".xyz"); want = r.get("xyz_sha256")
    if not os.path.isfile(src) or sha(src) != want: print("alt_prep: %s is missing or not the recorded geometry" % src); sys.exit(1)
    os.makedirs(d)
    shutil.copyfile(src, os.path.join(d, "start.xyz"))
    shutil.copyfile(os.path.join(rd, step, part + ".inp"), os.path.join(d, part + ".inp"))
    os.symlink("../../../complex_solvated.ORCAFF.prms", os.path.join(d, "complex_solvated.ORCAFF.prms"))
    print("%s: start.xyz from %s (sha256 %s)" % (sub, os.path.relpath(src, rd), want))

def start_gradient(csv, table, src="table"):
    """the environment gradient where the step began (SP_AMENDMENT1.txt A1.3: Zhang et al.'s criterion, reported only).
    ORCA leaves line 1 of the CSV without gradients; line 2 is after the optimiser's first step, which in the pilot's
    twelve runs was 5e-6 to 3e-4 A. A step that converged at its first point has one line: its printed table is then the
    start. Never fails the check: anything unreadable is recorded as n/a."""
    try:
        G = [l.split(";") for l in open(csv) if l.strip() and not l.lstrip().startswith("#")]
        if len(G) >= 2 and len(G[1]) > 8 and G[1][7].strip() and G[1][8].strip():
            rms, mx = float(G[1][7]), float(G[1][8])
            note = "CSV line 2, %s A from the start" % (G[1][9].strip() if len(G[1]) > 9 and G[1][9].strip() else "?")
        elif len(G) == 1 and len(table) == 2:
            rms, mx = table["RMS Grad."][0], table["Max. Grad."][0]
            note = ("printed table" if src == "table" else "CSV last line") + " (converged at its first point)"
        else:
            return [("start_rms_grad_kJmolA", "n/a"), ("start_max_grad_kJmolA", "n/a"), ("start_grad_from", "n/a")]
        return [("start_rms_grad_kJmolA", "%.4f" % rms), ("start_max_grad_kJmolA", "%.4f" % mx), ("start_grad_from", note),
                ("start_meets_zhang_rms", "yes" if rms <= ZHANG_RMS_KJ else "no")]
    except (OSError, IndexError, ValueError):
        return [("start_rms_grad_kJmolA", "n/a"), ("start_max_grad_kJmolA", "n/a"), ("start_grad_from", "n/a")]

def alt_check(rd, step, k, part):
    sub = os.path.join(step, cyc(k, part)); d = os.path.join(rd, sub)
    t = output(rd, sub, part + ".out")
    reg = set(region(rd)); keep = (reg - set(QMS)) if part == "env" else set(QMS)
    it = [("run_seconds", run_seconds(t))]
    E = energies(t)
    if part == "env":
        csv = os.path.join(d, "env-minimize-ener.csv")
        EP = [float(l.split(";")[5]) for l in open(csv) if l.strip() and not l.lstrip().startswith("#")] if os.path.isfile(csv) else []
        it += [("steps", str(len(EP)))]
        it += start_gradient(csv, {})            # reported only; recorded before any check can fail the step
        if "ORCA TERMINATED NORMALLY" not in t: fail(rd, sub, "ORCA did not terminate normally", it)
        k = t.rfind("The minimization has converged")
        if k < 0: fail(rd, sub, "the minimisation did not converge", it)
        g = dict((m.group(1), (float(m.group(2)), float(m.group(3)))) for m in
                 re.finditer(r"^\s*(Max\. Grad\.|RMS Grad\.)\s*\|\s*(\d+\.\d+)\s*\|\s*(\d+\.\d+)", t[k:k + 2000], re.M))
        src = "table"
        if len(g) != 2:    # no printed table: the CSV's last row (RMS and Max gradient, kJ/mol/A) against TolRMSG/TolMaxG
            G = [l.split(";") for l in open(csv) if l.strip() and not l.lstrip().startswith("#")] if os.path.isfile(csv) else []
            if not G or len(G[-1]) < 9 or not G[-1][8].strip() or not G[-1][7].strip():
                fail(rd, sub, "no gradient table after \"The minimization has converged\" and no gradient in the CSV's last line", it)
            g = {"Max. Grad.": (float(G[-1][8]), GMAX_KJ), "RMS Grad.": (float(G[-1][7]), GRMS_KJ)}; src = "csv"
        it += [("max_grad_kJmolA", "%.4f/%.4f" % g["Max. Grad."]), ("rms_grad_kJmolA", "%.4f/%.4f" % g["RMS Grad."]),
               ("gradients_from", src)]
        if g["Max. Grad."][0] > g["Max. Grad."][1] or g["RMS Grad."][0] > g["RMS Grad."][1]:
            fail(rd, sub, "a printed gradient is above its threshold", it)
        if len(EP) < 1: fail(rd, sub, "no optimiser energies (env-minimize-ener.csv)", it)
        if len(EP) == 1:                         # converged at its first point: the start gradient is the final one
            it = [kv for kv in it if not kv[0].startswith("start_")] + start_gradient(csv, g, src)
        slide = EP[-1 - PLATEAU_N] - EP[-1] if len(EP) > PLATEAU_N else None
        it += [("E_first_Eh", "%.8f" % EP[0]), ("E_last_Eh", "%.8f" % EP[-1]), ("drop_kcal", "%.4f" % ((EP[0] - EP[-1]) * H)),
               ("last50_fall_Eh", "n/a (%d steps)" % len(EP) if slide is None else "%.2e" % slide)]
        if slide is not None and slide > PLATEAU_EH:
            fail(rd, sub, "energy still falling: %.2e Eh over the last %d steps" % (slide, PLATEAU_N), it)
    else:
        it += [("cycles", str(t.count("GEOMETRY OPTIMIZATION CYCLE")))]
        if "ORCA TERMINATED NORMALLY" not in t: fail(rd, sub, "ORCA did not terminate normally", it)
        if "THE OPTIMIZATION HAS CONVERGED" not in t: fail(rd, sub, "the optimisation did not converge", it)
    if not E: fail(rd, sub, "no FINAL SINGLE POINT ENERGY (QM/MM)", it)
    it.append(("energy_Eh", E[-1]))
    px = os.path.join(d, part + ".xyz")
    if not os.path.isfile(px): fail(rd, sub, "no %s.xyz" % part, it)
    R = read_xyz(px); S = read_xyz(os.path.join(d, "start.xyz"))
    if R is None or S is None: fail(rd, sub, "%s.xyz or start.xyz is not the full system" % part, it)
    mv, same = moved(R, S, keep)
    it.append(("fixed_atoms_max_move_A", "%.1e" % mv))
    if not same: fail(rd, sub, "elements differ between start.xyz and %s.xyz" % part, it)
    if mv > MOVE_TOL: fail(rd, sub, "an atom that should be fixed moved (%s step)" % part, it)
    br, fo = rc(R); kd = kind_of(br, fo)
    it += [("break_A", "%.3f" % br), ("form_A", "%.3f" % fo), ("substrate", kd)]
    if part == "core" and kd != WANT[step]: fail(rd, sub, "the substrate is %s, expected %s" % (kd, WANT[step]), it)
    it += [("xyz_sha256", sha(px)), ("start_sha256", sha(os.path.join(d, "start.xyz")))]
    result(rd, sub, "PASS", it)

def history(rd, step):
    rows = []
    for k in range(1, KMAX + 1):
        e, c = read_result(rd, os.path.join(step, cyc(k, "env"))), read_result(rd, os.path.join(step, cyc(k, "core")))
        if not e: break
        rows.append("%d:%s/%s,%s,%s" % (k, e.get("status", "-"), c.get("status", "-"), e.get("drop_kcal", "-"), c.get("energy_Eh", "-")))
    return " ".join(rows)

def alt_decide(rd, step, k):
    ce = read_result(rd, os.path.join(step, cyc(k, "env"))); cc = read_result(rd, os.path.join(step, cyc(k, "core")))
    if ce.get("status") != "PASS" or cc.get("status") != "PASS": fail(rd, step, "cycle %d did not pass" % k, [("history", history(rd, step))])
    msg = "cycle %d: environment step lowered the energy by %s kcal/mol; core energy %s Eh" % (k, ce["drop_kcal"], cc["energy_Eh"])
    if float(ce["drop_kcal"]) < 0: msg += " (the environment step RAISED the energy)"
    if k >= 2:
        pc = read_result(rd, os.path.join(step, cyc(k - 1, "core")))
        dE = abs(float(cc["energy_Eh"]) - float(pc["energy_Eh"])) * H
        A = read_xyz(os.path.join(rd, step, cyc(k - 1, "core"), "core.xyz")); B = read_xyz(os.path.join(rd, step, cyc(k, "core"), "core.xyz"))
        dx = max(dist(fl(A[i]), fl(B[i])) for i in QMS)
        msg += "; core energy change from cycle %d %.4f kcal/mol; largest substrate-atom shift %.4f A" % (k - 1, dE, dx)
        if abs(float(ce["drop_kcal"])) <= ALT_DE and dE <= ALT_DE and dx <= ALT_DX:
            src = os.path.join(rd, step, cyc(k, "core"), "core.xyz"); dst = os.path.join(rd, step, "final.xyz")
            shutil.copyfile(src, dst)
            if sha(dst) != cc["xyz_sha256"]: fail(rd, step, "final.xyz copy differs from cycle %d's core.xyz" % k)
            R = read_xyz(dst); reg = region(rd)
            with open(os.path.join(rd, step, "final_region.xyz"), "w") as fh:
                fh.write("%d\n%s/%s: converged alternation, cycle %d core geometry, region atoms; 5th column: 0-based atom index\n" % (
                    len(reg), os.path.basename(os.path.abspath(rd)), step, k)
                    + "".join("%s %s %s %s %d\n" % (R[i][0], R[i][1], R[i][2], R[i][3], i) for i in reg))
            print(msg)
            result(rd, step, "PASS", [("cycles", str(k)), ("energy_Eh", cc["energy_Eh"]), ("last_env_drop_kcal", ce["drop_kcal"]),
                                      ("last_core_change_kcal", "%.4f" % dE), ("last_core_shift_A", "%.4f" % dx),
                                      ("break_A", cc["break_A"]), ("form_A", cc["form_A"]), ("substrate", cc["substrate"]),
                                      ("history", history(rd, step)), ("xyz_sha256", cc["xyz_sha256"]),
                                      ("start_sha256", open(os.path.join(rd, step, "start.sha256")).read().split()[0])])
    print(msg)
    if k >= KMAX: fail(rd, step, "the alternation did not converge in %d cycles" % KMAX, [("history", history(rd, step))])
    print("%s: not converged after cycle %d; next cycle" % (step, k)); sys.exit(3)

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
    it.append(("hessian_atom_warning", "yes" if "corresponding active atoms in Hessian file" in t else "no"))
    k = t.rfind("IRC PATH SUMMARY")
    rows = re.findall(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+\d+\.\d+\s+\d+\.\d+(.*)$", t[k:], re.M) if k >= 0 else []
    if rows:
        ts_rows = [r for r in rows if "TS" in r[3]]
        it += [("path_points", str(len(rows))), ("E_path_first_Eh", rows[0][1]), ("E_path_last_Eh", rows[-1][1]),
               ("E_ts_Eh", ts_rows[0][1] if ts_rows else "")]
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
        elif len(a) == 5 and a[0] == "alt_prep": alt_prep(a[1], a[2], int(a[3]), a[4])
        elif len(a) == 5 and a[0] == "alt_check" and a[2] in ALT_STEPS and a[4] in ("env", "core"): alt_check(a[1], a[2], int(a[3]), a[4])
        elif len(a) == 4 and a[0] == "alt_decide" and a[2] in ALT_STEPS: alt_decide(a[1], a[2], int(a[3]))
        elif len(a) == 4 and a[0] == "alt_fail" and a[2] in ALT_STEPS: fail(a[1], a[2], a[3], [("history", history(a[1], a[2]))])
        elif len(a) == 4 and a[0] == "check_freq" and a[3] in ("ts", "min"): check_freq(a[1], a[2], a[3])
        elif len(a) == 2 and a[0] == "check_irc": check_irc(a[1])
        else: sys.exit(__doc__)
    except Exception as e:      # an unreadable or truncated file: recorded as a failed check, never as "not run"
        if a and a[0] in ("alt_prep", "prep"): print("%s failed (machinery stop, no record): %s: %s" % (a[0], type(e).__name__, e)); sys.exit(1)
        step = "irc" if a[0] == "check_irc" else (os.path.join(a[2], cyc(int(a[3]), a[4])) if a[0] == "alt_check" else a[2])
        fail(a[1], step, "check could not read its files: %s: %s" % (type(e).__name__, e))
