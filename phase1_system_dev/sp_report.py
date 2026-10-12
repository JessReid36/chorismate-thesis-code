#!/usr/bin/env python3
"""
sp_report.py - read-only report of the stationary-point barriers (results 05_qmmm/24_sp_barriers/SP_PLAN.txt). Runs on
hpc1 or on a fetched copy. Python 3.6, standard library; deterministic; writes nothing.

  sp_report.py RD [RD ...]       per run folder (f<FRAME>_r<R>[_ci]): every job and step, the quantities of SP_PLAN.txt
  sp_report.py smoke SD FRAME    the region rule of SP_PLAN.txt on SD/f<FRAME>_r8 and SD/f<FRAME>_r6

Energies are re-read from the ORCA outputs (.out or .out.xz) and must equal the in-job record (STAGE_RESULT.txt);
any difference is printed as PROBLEM. kcal/mol = Eh x 627.5094740631.
"""
import os, re, subprocess, sys

H = 627.5094740631
STEPS = ["ts", "ts_freq", "irc", "r_irc", "r_irc_freq", "p_irc", "p_irc_freq", "r0", "r0_freq"]
BASE = {"ts": "ts", "ts_freq": "freq", "irc": "irc", "r_irc": "opt", "p_irc": "opt", "r0": "opt",
        "r_irc_freq": "freq", "p_irc_freq": "freq", "r0_freq": "freq"}
JOBS = [("job_ts", ["ts"]), ("job_ts_freq_irc", ["ts_freq", "irc"]), ("job_r_irc", ["r_irc"]),
        ("job_r_irc_freq", ["r_irc_freq"]), ("job_p_irc", ["p_irc"]), ("job_p_irc_freq", ["p_irc_freq"]),
        ("job_r0", ["r0"]), ("job_r0_freq", ["r0_freq"])]
OPTS = ("ts", "r_irc", "p_irc", "r0")
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
def job_state(rd, job):
    p = os.path.join(rd, job + ".pbs.out")
    if not os.path.isfile(p): return "not started", {}, ""
    t = open(p, errors="replace").read()
    use = dict(re.findall(r"^pbs_resources_used\.(\w+) = (\S+)", t, re.M))
    kb = " ".join(re.findall(r"folder_kB_after_cleanup=(\d+)", t))
    if "PBS: job killed" in t: return "killed by PBS (" + t[t.find("PBS: job killed"):].splitlines()[0] + ")", use, kb
    if "end_job=" in t: return "ended", use, kb
    return "running, or ended without its last line", use, kb
def size_before(rd, step):
    p = os.path.join(rd, step, "FILES_BEFORE_CLEANUP.txt")
    if not os.path.isfile(p): return None, None
    tot = 0; big = (0, "")
    for l in open(p):
        f = l.split(None, 1)
        if len(f) == 2 and f[0].isdigit():
            n = int(f[0]); tot += n
            if n > big[0]: big = (n, f[1].strip())
    return tot, big
def energy_from_out(rd, step, first=False):
    t = text(os.path.join(rd, step, BASE[step] + ".out"))
    if t is None: return None
    E = re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", t)
    return (E[0] if first else E[-1]) if E else None

def report(rd):
    rd = os.path.abspath(rd); name = os.path.basename(rd)
    print("=" * 100)
    nreg = len(open(os.path.join(rd, "region.txt")).read().split()) if os.path.isfile(os.path.join(rd, "region.txt")) else 0
    print("%s: region %d atoms" % (name, nreg))
    print("jobs:")
    for job, steps in JOBS:
        st, use, kb = job_state(rd, job)
        print("  %-16s %-40s mem %-10s walltime %-9s folder after clean-up %s kB" % (
            job, st, use.get("mem", "-"), use.get("walltime", "-"), kb or "-"))
    print("steps:")
    R = {}
    for s in STEPS:
        r = result(rd, s); R[s] = r
        if not r: print("  %-11s not run" % s); continue
        keys = [k for k in r if k not in ("status", "reason") and not k.endswith("sha256")]
        print("  %-11s %-5s %s%s" % (s, r.get("status"), ("(" + r["reason"] + ") ") if r.get("reason") else "",
                                    ", ".join("%s=%s" % (k, r[k]) for k in keys)))
        tot, big = size_before(rd, s)
        if tot is not None: print("              files before clean-up %.1f MB (largest %s, %.1f MB)" % (tot / 1e6, big[1], big[0] / 1e6))
    # in-job record against the outputs
    for s in OPTS:
        r = R[s]
        if r.get("energy_Eh"):
            e = energy_from_out(rd, s)
            if e != r["energy_Eh"]: PROBLEMS.append("%s %s: energy in STAGE_RESULT %s, in the output %s" % (name, s, r["energy_Eh"], e))
    ok = lambda s: R[s].get("status") == "PASS"
    E = dict((s, float(R[s]["energy_Eh"])) for s in OPTS if ok(s))
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
    return R

def smoke(sd, frame):
    verdict = {}
    for rad in (8, 6):
        rd = os.path.join(sd, "f%s_r%d" % (frame, rad))
        R = report(rd); why = []; pending = []
        for job, steps in JOBS:
            st, use, kb = job_state(rd, job)
            if st.startswith("killed"): why.append(job + " " + st)
            elif st != "ended": pending.append(job + ": " + st)
        for s in STEPS:
            r = R[s]
            if not r: continue
            if r.get("status") == "PASS": continue
            rs = r.get("reason", "")
            if s in ("r_irc_freq", "p_irc_freq", "r0_freq") and "counted imaginary modes" in rs: continue   # a result
            why.append("%s: %s" % (s, rs))
        if why: verdict[rad] = "fails (" + "; ".join(why) + ")"
        elif pending or any(not R[s] for s in STEPS): verdict[rad] = "not decided yet (" + "; ".join(pending or ["steps not run"]) + ")"
        else: verdict[rad] = "meets the rule"
    print("=" * 100)
    print("REGION RULE (SP_PLAN.txt), frame %s: every job ended inside its walltime with ORCA terminating normally; the TS,"
          " R_IRC, P_IRC and R0 optimisations passed their checks (converged; R_IRC and R0 reactant-like, P_IRC"
          " product-like); the TS has one counted imaginary mode; the IRC ends are one reactant and one product" % frame)
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
    print("=" * 100)
    print("PROBLEM: " + "\nPROBLEM: ".join(PROBLEMS) if PROBLEMS else "record against outputs: no problem")
