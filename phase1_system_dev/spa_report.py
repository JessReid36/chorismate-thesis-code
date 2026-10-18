#!/usr/bin/env python3
"""
spa_report.py - read-only report of the stationary-point barriers under SP Amendment 1, the alternating scheme (results
05_qmmm/24_sp_barriers/SP_AMENDMENT1.txt and SP_PLAN.txt), reading up to SP Amendment 2's 20 cycles (A2.6). Runs on hpc1 or on a fetched copy. Python 3.6, standard
library; deterministic; writes nothing.

  spa_report.py RD [RD ...]       per run folder (f<FRAME>_r<R>[_ci]_a1): every job, every step, every alternation
                                  cycle, and the quantities of SP_PLAN.txt
  spa_report.py smoke SD FRAME    the region rule (SP_AMENDMENT1.txt) on SD/f<FRAME>_r8_a1 and SD/f<FRAME>_r6_a1; it is
                                  Amendment 1's rule, which smoke test 2 decided, and its text keeps Amendment 1's
                                  8 cycles (SP_AMENDMENT2.txt A2.6); not for the package

Energies are re-read from the ORCA outputs (.out or .out.xz) and must equal the in-job records (STAGE_RESULT.txt);
any difference is printed as PROBLEM. kcal/mol = Eh x 627.5094740631.
"""
import os, re, subprocess, sys

H = 627.5094740631
QM0, NQM = 6207, 24
QMS = list(range(QM0, QM0 + NQM))
KMAX = 20                    # SP_AMENDMENT2.txt A2.6 (Amendment 1: 8); the smoke verdict below is Amendment 1's
STEPS = ["ts", "ts_freq", "irc", "r_irc", "r_irc_freq", "p_irc", "p_irc_freq", "r0", "r0_freq"]
BASE = {"ts_freq": "freq", "irc": "irc", "r_irc_freq": "freq", "p_irc_freq": "freq", "r0_freq": "freq"}
JOBS = [("job_ts", ["ts"]), ("job_ts_freq_irc", ["ts_freq", "irc"]), ("job_r_irc", ["r_irc"]),
        ("job_r_irc_freq", ["r_irc_freq"]), ("job_p_irc", ["p_irc"]), ("job_p_irc_freq", ["p_irc_freq"]),
        ("job_r0", ["r0"]), ("job_r0_freq", ["r0_freq"])]
ALT = ("ts", "r_irc", "p_irc", "r0")
PROBLEMS = []

def text(p):
    if os.path.isfile(p): return open(p, errors="replace").read()
    if os.path.isfile(p + ".xz"):
        try:
            import lzma
            return lzma.open(p + ".xz", "rt", errors="replace").read()
        except ImportError:
            return subprocess.run(["xz", "-dc", p + ".xz"], stdout=subprocess.PIPE).stdout.decode(errors="replace")
    return None
def result(rd, step):
    p = os.path.join(rd, step, "STAGE_RESULT.txt")
    if not os.path.isfile(p): return {}
    return dict(l.rstrip("\n").split("=", 1) for l in open(p) if "=" in l)
def cyc(k, part): return "c%02d_%s" % (k, part)
def job_state(rd, job):
    p = os.path.join(rd, job + ".pbs.out")
    if not os.path.isfile(p): return "not started", {}, "", []
    t = open(p, errors="replace").read()
    use = dict(re.findall(r"^pbs_resources_used\.(\w+) = (\S+)", t, re.M))
    kb = " ".join(re.findall(r"folder_kB_after_cleanup=(\d+)", t))
    fails = [l.strip() for l in t.splitlines() if re.match(r"^SP_\S+_FAIL", l.strip())]
    if "PBS: job killed" in t: return "killed by PBS (" + t[t.find("PBS: job killed"):].splitlines()[0] + ")", use, kb, fails
    if "end_job=" in t: return "ended", use, kb, fails
    return "running, or ended without its last line", use, kb, fails
def size_before(rd, sub):
    p = os.path.join(rd, sub, "FILES_BEFORE_CLEANUP.txt")
    if not os.path.isfile(p): return None, None
    tot = 0; big = (0, "")
    for l in open(p):
        f = l.split(None, 1)
        if len(f) == 2 and f[0].isdigit():
            n = int(f[0]); tot += n
            if n > big[0]: big = (n, f[1].strip())
    return tot, big
def last_energy(rd, sub, base):
    t = text(os.path.join(rd, sub, base + ".out"))
    if t is None: return None
    E = re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", t)
    return E[-1] if E else None
def core_xyz(rd, step, k):
    p = os.path.join(rd, step, cyc(k, "core"), "core.xyz")
    try:
        L = open(p).read().splitlines()
        return [tuple(float(v) for v in L[2 + i].split()[1:4]) for i in QMS]
    except (OSError, IndexError, ValueError):
        return None
def csv_line2_rms(rd, sub):
    p = os.path.join(rd, sub, "env-minimize-ener.csv")
    try:
        G = [l.split(";") for l in open(p) if l.strip() and not l.lstrip().startswith("#")]
        return float(G[1][7]) if len(G) >= 2 and len(G[1]) > 8 and G[1][7].strip() and G[1][8].strip() else None
    except (OSError, IndexError, ValueError): return None
def csv_energies(rd, sub):
    p = os.path.join(rd, sub, "env-minimize-ener.csv")
    try: return [float(l.split(";")[5]) for l in open(p) if l.strip() and not l.lstrip().startswith("#")]
    except (OSError, IndexError, ValueError): return None
def dist(a, b): return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5
def hours(s):
    try: return "%.1f h" % (int(s) / 3600.0)
    except ValueError: return "-"

def cycles(rd, name, step):
    """the alternation table of one step; checks every recorded energy against its output"""
    rows = []; prev = None; prevx = None
    ss = os.path.join(rd, step, "start.sha256")
    chain = open(ss).read().split()[0] if os.path.isfile(ss) else None     # the sha256 the next part must start from
    for k in range(1, KMAX + 1):
        e, c = result(rd, os.path.join(step, cyc(k, "env"))), result(rd, os.path.join(step, cyc(k, "core")))
        if not e and not os.path.isdir(os.path.join(rd, step, cyc(k, "env"))): break
        for part, r in (("env", e), ("core", c)):
            sub = os.path.join(step, cyc(k, part))
            if r.get("energy_Eh"):
                eo = last_energy(rd, sub, part)
                if eo != r["energy_Eh"]:
                    PROBLEMS.append("%s %s: energy in STAGE_RESULT %s, in the output %s" % (name, sub, r["energy_Eh"], eo))
            if r.get("start_sha256"):
                if r["start_sha256"] != chain:
                    PROBLEMS.append("%s %s: started from sha256 %s, the previous part recorded %s" % (name, sub, r["start_sha256"], chain))
                chain = r.get("xyz_sha256")
            if part == "env" and r.get("drop_kcal"):
                EP = csv_energies(rd, sub)
                if not EP:
                    PROBLEMS.append("%s %s: env-minimize-ener.csv missing or unreadable" % (name, sub))
                else:
                    sl = EP[-51] - EP[-1] if len(EP) > 50 else None
                    if "%.4f" % ((EP[0] - EP[-1]) * H) != r["drop_kcal"] or str(len(EP)) != r.get("steps") or \
                            ("n/a (%d steps)" % len(EP) if sl is None else "%.2e" % sl) != r.get("last50_fall_Eh"):
                        PROBLEMS.append("%s %s: steps, drop or last-50 fall in STAGE_RESULT differ from the CSV" % (name, sub))
                g2 = csv_line2_rms(rd, sub)
                if g2 is not None and r.get("start_rms_grad_kJmolA") not in ("%.4f" % g2, None, "n/a"):
                    PROBLEMS.append("%s %s: start_rms_grad_kJmolA %s, CSV line 2 %.4f" % (name, sub, r.get("start_rms_grad_kJmolA"), g2))
        dE = dx = j = "-"
        if prev is not None and e.get("E_first_Eh"): j = "%+.4f" % ((float(e["E_first_Eh"]) - prev) * H)
        if c.get("energy_Eh"):
            if prev is not None: dE = "%+.4f" % ((float(c["energy_Eh"]) - prev) * H)
            prev = float(c["energy_Eh"])
            x = core_xyz(rd, step, k)
            if x is not None and prevx is not None: dx = "%.4f" % max(dist(a, b) for a, b in zip(x, prevx))
            prevx = x
        rows.append("    %2d  env %-4s %5s steps start %8s grad %s/%s drop %9s last50 %8s %7s | core %-4s %4s cycles E %s dE %10s"
                    " dx %7s %s/%s %s %7s%s" % (
                        k, e.get("status", "-"), e.get("steps", "-"), j, e.get("start_rms_grad_kJmolA", "-"),
                        e.get("start_max_grad_kJmolA", "-"), e.get("drop_kcal", "-"), e.get("last50_fall_Eh", "-"),
                        hours(e.get("run_seconds", "")), c.get("status", "-"), c.get("cycles", "-"), c.get("energy_Eh", "-"),
                        dE, dx, c.get("break_A", "-"), c.get("form_A", "-"), c.get("substrate", "-"),
                        hours(c.get("run_seconds", "")),
                        "".join("  [%s: %s]" % (p, r["reason"]) for p, r in (("env", e), ("core", c)) if r.get("reason"))))
        sz = []
        for part in ("env", "core"):
            tot, big = size_before(rd, os.path.join(step, cyc(k, part)))
            if tot is not None: sz.append("%s %.1f MB (largest %s %.1f MB)" % (part, tot / 1e6, big[1], big[0] / 1e6))
        if sz: rows.append("        files before clean-up: " + "; ".join(sz))
    return rows

def report(rd):
    rd = os.path.abspath(rd); name = os.path.basename(rd)
    print("=" * 120)
    nreg = len(open(os.path.join(rd, "region.txt")).read().split()) if os.path.isfile(os.path.join(rd, "region.txt")) else 0
    print("%s: region %d atoms (environment %d, core %d)" % (name, nreg, max(nreg - NQM, 0), NQM))
    print("jobs:")
    for job, steps in JOBS:
        st, use, kb, fails = job_state(rd, job)
        print("  %-16s %-40s mem %-10s walltime %-9s folder after clean-up %s kB" % (
            job, st, use.get("mem", "-"), use.get("walltime", "-"), kb or "-"))
        for l in fails: print("                   " + l)
    print("steps:")
    R = {}
    ended = dict((s, job_state(rd, job)[0] != "running, or ended without its last line") for job, steps in JOBS for s in steps)
    for s in STEPS:
        r = result(rd, s); R[s] = r
        if not r and not (s in ALT and os.path.isdir(os.path.join(rd, s, cyc(1, "env")))):
            print("  %-11s not run" % s); continue
        keys = [k for k in r if k not in ("status", "reason", "history") and not k.endswith("sha256")]
        st = r.get("status") or ("stopped (its job has ended without a record: see the job's lines above)" if ended[s] else "running")
        print("  %-11s %-5s %s%s" % (s, st, ("(" + r["reason"] + ") ") if r.get("reason") else "",
                                    ", ".join("%s=%s" % (k, r[k]) for k in keys)))
        if s in ALT:
            print("    cycle: environment step (L-Opt, substrate fixed): optimiser steps, start = its first energy minus the previous"
                  " core energy (kcal/mol), grad = RMS/largest environment gradient where it began (kJ/mol/A; Zhang et al.'s"
                  " criterion is RMS <= 0.4184, reported only), drop (kcal/mol), last50 = fall over the last 50 steps (Eh),"
                  " wall time | core step"
                  " (substrate alone): optimiser cycles, E (Eh), dE from the previous core (kcal/mol), largest substrate-atom"
                  " shift dx (A), break/form (A), wall time")
            for row in cycles(rd, name, s): print(row)
        else:
            tot, big = size_before(rd, s)
            if tot is not None: print("              files before clean-up %.1f MB (largest %s, %.1f MB)" % (tot / 1e6, big[1], big[0] / 1e6))
    # the step record against its final cycle
    for s in ALT:
        r = R[s]
        if r.get("energy_Eh"):
            k = int(r.get("cycles", "0"))
            e = last_energy(rd, os.path.join(s, cyc(k, "core")), "core")
            if e != r["energy_Eh"]: PROBLEMS.append("%s %s: energy in STAGE_RESULT %s, in %s/core.out %s" % (name, s, r["energy_Eh"], cyc(k, "core"), e))
            a, b = core_xyz(rd, s, k - 1), core_xyz(rd, s, k)
            if a is None or b is None or "%.4f" % max(dist(x, y) for x, y in zip(a, b)) != r.get("last_core_shift_A"):
                PROBLEMS.append("%s %s: last_core_shift_A %s not reproduced from the core geometries" % (name, s, r.get("last_core_shift_A")))
    ok = lambda s: R[s].get("status") == "PASS"
    E = dict((s, float(R[s]["energy_Eh"])) for s in ALT if ok(s))
    print("quantities (kcal/mol):")
    complete = ok("ts") and ok("ts_freq") and ok("irc") and ok("r_irc") and ok("r_irc_freq")
    if complete:
        print("  B_sp = E(TS) - E(R_IRC)               %8.2f   complete (TS: one imaginary mode; IRC: reactant and"
              " product ends; R_IRC: a minimum)" % ((E["ts"] - E["r_irc"]) * H))
    else:
        miss = [s for s in ("ts", "ts_freq", "irc", "r_irc", "r_irc_freq") if not ok(s)]
        print("  B_sp: none - the frame is not complete (not passed: %s)" % ", ".join(miss))
    if ok("ts_freq"): print("  TS imaginary mode (cm-1)               %s" % R["ts_freq"].get("imaginary_counted"))
    if "r0" in E and "r_irc" in E: print("  dR0 = E(R0) - E(R_IRC)                 %8.2f" % ((E["r0"] - E["r_irc"]) * H))
    if "ts" in E and "r0" in E: print("  E(TS) - E(R0)                          %8.2f   (reported only)" % ((E["ts"] - E["r0"]) * H))
    if "p_irc" in E and "r_irc" in E: print("  dE_rxn = E(P_IRC) - E(R_IRC)           %8.2f" % ((E["p_irc"] - E["r_irc"]) * H))
    for s in ("r_irc_freq", "p_irc_freq", "r0_freq"):
        if R[s]: print("  %-11s counted imaginary modes: %s" % (s, R[s].get("imaginary_counted", "-")))
    # the SP_PLAN run of the same frame and radius (exploratory since Amendment 1): same start files, reported only
    if name.endswith("_a1"):
        od = rd[:-3]
        for s in ("ts", "r0"):
            o = result(od, s)
            if s in E and o.get("status") == "PASS" and o.get("energy_Eh"):
                print("  %-3s this run minus the SP_PLAN run %s (same start file; reported only): %8.2f" % (
                    s.upper(), os.path.basename(od), (E[s] - float(o["energy_Eh"])) * H))
    return R

DISK = ("no space left on device", "disk quota exceeded")
def disk_failure(rd, step, job):
    """True when a failed step's ORCA output, or its job's output, says the disk or quota was full (SP_AMENDMENT1.txt A1.5)"""
    folders = [os.path.join(rd, step)] + sorted(os.path.join(rd, step, d) for d in (os.listdir(os.path.join(rd, step))
                                                                                  if os.path.isdir(os.path.join(rd, step)) else [])
                                                 if re.match(r"^c\d\d_(env|core)$", d))
    texts = [text(os.path.join(f, b + ".out")) for f in folders for b in ("env", "core", "freq", "irc")]
    texts.append(text(os.path.join(rd, job + ".pbs.out")))
    return any(x and any(s in x.lower() for s in DISK) for x in texts)
def smoke(sd, frame):
    verdict = {}
    for rad in (8, 6):
        rd = os.path.join(sd, "f%s_r%d_a1" % (frame, rad))
        R = report(rd); why = []; pending = []; stopped = []
        for job, steps in JOBS:
            st, use, kb, fails = job_state(rd, job)
            if st.startswith("killed"): why.append(job + " " + st)
            elif st != "ended": pending.append(job + ": " + st)
            elif fails and not any(R[s].get("status") == "FAIL" for s in steps):
                stopped.append(job + ": " + "; ".join(fails))
        for s in STEPS:
            r = R[s]
            if not r: continue
            if r.get("status") == "PASS": continue
            rs = r.get("reason", "")
            if s in ("r_irc_freq", "p_irc_freq", "r0_freq") and "counted imaginary modes" in rs: continue   # a result
            job = [j for j, steps in JOBS if s in steps][0]
            if disk_failure(rd, s, job): stopped.append("%s: %s, and its output reports a full disk or quota" % (s, rs)); continue
            why.append("%s: %s" % (s, rs))
        again = os.path.isdir(rd + "_stop1")
        if why: verdict[rad] = "fails (" + "; ".join(why) + ")"
        elif stopped and again:
            verdict[rad] = "fails (a second machinery stop, after %s_stop1: %s)" % (os.path.basename(rd), "; ".join(stopped))
        elif stopped: verdict[rad] = ("not decided: machinery stop (SP_AMENDMENT1.txt A1.5: rename the folder to %s_stop1,"
                                      " generate it again and resubmit once): " % os.path.basename(rd)) + "; ".join(stopped)
        elif pending or any(not R[s] for s in STEPS): verdict[rad] = "not decided yet (" + "; ".join(pending or ["steps not run"]) + ")"
        else: verdict[rad] = "meets the rule"
    print("=" * 120)
    print("REGION RULE (SP_AMENDMENT1.txt), frame %s: every job ended inside its walltime with ORCA terminating normally;"
          " the TS, R_IRC, P_IRC and R0 alternations converged within 8 cycles with every cycle passing its checks (R_IRC"
          " and R0 reactant-like, P_IRC product-like, the TS neither); the TS has one counted imaginary mode; the IRC ends"
          " are one reactant and one product" % frame)
    for rad in (8, 6): print("  %d A: %s" % (rad, verdict[rad]))
    if verdict[8] == "meets the rule": print("  -> the 30 frames use 8 A")
    elif verdict[8].startswith("fails") and verdict[6] == "meets the rule": print("  -> the 30 frames use 6 A")
    elif verdict[8].startswith("fails") and verdict[6].startswith("fails"): print("  -> neither radius meets the rule: stop and report")
    else: print("  -> not decided yet")

if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) == 3 and a[0] == "smoke": smoke(a[1], a[2])
    elif a and a[0] != "smoke":
        for rd in a: report(rd)
    else: sys.exit(__doc__)
    print("=" * 120)
    print("PROBLEM: " + "\nPROBLEM: ".join(PROBLEMS) if PROBLEMS else "records against outputs: no problem")
