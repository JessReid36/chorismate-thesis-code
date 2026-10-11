#!/usr/bin/env python3
"""
c4a2_pilot_summary.py - one line per pilot frame of the C4 Amendment 2 pilot (exploratory; PILOT_A2.txt and
PILOT_A2_ADDENDUM1.txt) and the addendum's decision tally. Read-only; Python 3.6, standard library. Works while runs are
going. The full per-frame report is c4a2_pilot_report.py FRAME.
Reads, in the C4 folder D: c4_bands.tsv, active_atoms.txt, S2_F/, pilot_a2/T/T1_F/, pilot_a2/R/R2_F/ for every frame F
with a pilot_a2/T/T1_F folder.
USAGE  python3 c4a2_pilot_summary.py [D]
"""
import glob, lzma, os, re, sys
D = sys.argv[1] if len(sys.argv) > 1 else "/home/18660916/system_development/05_qmmm/23_c4_reactant"
H = 627.5094740631
QM0, NQM, NATOMS = 6207, 24, 55680
PLATEAU_EH, PLATEAU_N = 2.5e-4, 50
TOL_START, TOL_DR, NSTOP = 0.01, 1.0, 3
def text(p):
    if os.path.exists(p): return open(p, errors="replace").read()
    if os.path.exists(p + ".xz"): return lzma.open(p + ".xz", "rt", errors="replace").read()
    return None
def state(d, nact):
    """status code, energies and machinery problems of one optimisation folder"""
    r = {"code": "-", "E": [], "prob": []}
    if not os.path.isfile(os.path.join(d, "reopt.inp")): return r
    t = text(os.path.join(d, "reopt.out"))
    if t is None: r["code"] = "queued"; return r
    m = re.search(r"The system contains \d+ active degrees of freedom \((\d+) active atoms", t)
    if m and int(m.group(1)) != nact: r["prob"].append("%s active atoms" % m.group(1))
    if "ABORTING THE RUN" in t or "error termination" in t: r["prob"].append("ORCA error")
    cp = os.path.join(d, "reopt-minimize-ener.csv")
    if os.path.exists(cp):
        r["E"] = [float(l.split(";")[5]) for l in open(cp) if l.strip() and not l.lstrip().startswith("#")]
    po = text(os.path.join(d, "reopt.pbs.out")) or ""
    grad = "ORCA TERMINATED NORMALLY" in t and "The minimization has converged." in t
    E = r["E"]; slide = (E[-1 - PLATEAU_N] - E[-1]) if len(E) > PLATEAU_N else None
    plateau = slide is None or slide <= PLATEAU_EH
    if "end=" not in po: r["code"] = "run %d" % len(E)
    elif grad and plateau and "extract=ok" in po: r["code"] = "done %d" % len(E); r["done"] = True
    else:
        r["code"] = "ENDED %d" % len(E)
        if "extract=FAILED" in po: r["prob"].append("c4_extract failed (an atom outside the active atoms moved?)")
        elif not grad: r["prob"].append("gradients not met")
        elif not plateau: r["prob"].append("no energy plateau (continue)")
    return r
BL = [l.rstrip("\n").split("\t") for l in open(os.path.join(D, "c4_bands.tsv")) if l.strip() and not l.startswith("#")]
BAND = {r[0]: dict(zip(BL[0], r)) for r in BL[1:]}
ACT = [int(x) for x in open(os.path.join(D, "active_atoms.txt")).read().split()]
NENV = len([i for i in ACT if not QM0 <= i < QM0 + NQM])
frames = sorted(os.path.basename(p)[3:] for p in glob.glob(os.path.join(D, "pilot_a2", "T", "T1_*")))
print("C4 Amendment 2 pilot - summary (c4a2_pilot_summary.py); exploratory, not part of C4's result. kcal/mol;")
print("fall_T = E_ci - E_min(T1); dR = E_min(R2) - E_min(S2), test |dR| <= %.1f; B = E_ci - E_klow (band); B3 = E_min(T1) - E_min(R2)" % TOL_DR)
print("%-6s %-10s %-11s %7s  %-11s %7s  %-8s %6s %6s %6s" % ("frame", "S2", "T1", "fall_T", "R2", "dR", "test", "B", "B3", "B3-B"))
met = notmet = pend = 0; probs = []
for f in frames:
    b = BAND[f]; eci, ekl = float(b["E_ci"]), float(b["E_klow"]); B = (eci - ekl) * H
    s2 = state(os.path.join(D, "S2_" + f), len(ACT))
    t1 = state(os.path.join(D, "pilot_a2", "T", "T1_" + f), NENV)
    r2 = state(os.path.join(D, "pilot_a2", "R", "R2_" + f), len(ACT))
    for nm, r in (("S2", s2), ("T1", t1), ("R2", r2)):
        for p in r["prob"]: probs.append("%s %s: %s" % (f, nm, p))
    if t1["E"] and abs((t1["E"][0] - eci) * H) > TOL_START:
        probs.append("%s T1: E_start - E_ci = %+.4f (limit %.2f)" % (f, (t1["E"][0] - eci) * H, TOL_START))
    if t1.get("done"):
        S = text(os.path.join(D, "pilot_a2", "T", "T1_" + f, "start.xyz")); X = text(os.path.join(D, "pilot_a2", "T", "T1_" + f, "reopt.xyz"))
        if S and X:
            SL, XL = S.splitlines(), X.splitlines()
            mv = max(abs(float(a) - float(c)) for i in range(QM0, QM0 + NQM) for a, c in zip(SL[2 + i].split()[1:4], XL[2 + i].split()[1:4]))
            if mv > 2e-6: probs.append("%s T1: substrate moved %.1e A" % (f, mv))
    fall = "%7.2f" % ((eci - t1["E"][-1]) * H) if t1["E"] else "%7s" % "-"
    dr = test = b3 = b3b = None
    if r2["E"] and s2["E"]:
        dr = (r2["E"][-1] - s2["E"][-1]) * H
    if r2.get("done") and s2.get("done"):
        test = "met" if abs(dr) <= TOL_DR else "NOT met"
        if test == "met": met += 1
        else: notmet += 1
    else: pend += 1
    if t1.get("done") and r2.get("done"):
        b3 = (t1["E"][-1] - r2["E"][-1]) * H; b3b = b3 - B
    print("%-6s %-10s %-11s %s  %-11s %7s  %-8s %6.2f %6s %6s" % (
        f, s2["code"], t1["code"], fall, r2["code"], "-" if dr is None else "%+.2f" % dr, test or "pending",
        B, "-" if b3 is None else "%.2f" % b3, "-" if b3b is None else "%+.2f" % b3b))
print("\nfirst-cycle test over %d frames: met %d, not met %d, pending %d" % (len(frames), met, notmet, pend))
print("addendum 1 rule: revisit the design if %d or more of the six do not meet it; any machinery problem -> fix first" % NSTOP)
if notmet >= NSTOP: print("-> %d not met: the rule's threshold is reached" % notmet)
print("machinery problems: " + ("none" if not probs else "\n  " + "\n  ".join(probs)))
