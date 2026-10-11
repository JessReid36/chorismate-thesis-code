#!/usr/bin/env python3
"""
c4a2_pilot_report.py - report of the C4 Amendment 2 pilot (exploratory; pilot_a2/PILOT_A2.txt). Python 3.6, standard
library. Read-only: prints to the screen, writes nothing. Works while runs are going (progress) and after (results).
Reads, in the C4 folder D: c4_bands.tsv, active_atoms.txt, frames/frame_F/md_substrate.xyz, S2_F/ (the step-2 run),
pilot_a2/T/T1_F/, pilot_a2/R/R2_F/ and pilot_a2/{T,R}/active_atoms.txt.
USAGE  python3 c4a2_pilot_report.py FRAME [D]
"""
import lzma, os, re, sys
F = sys.argv[1]
D = sys.argv[2] if len(sys.argv) > 2 else "/home/18660916/system_development/05_qmmm/23_c4_reactant"
H = 627.5094740631
QM0, NQM, NATOMS = 6207, 24, 55680
EH_BOHR = 2625.49964 / 0.529177210903          # kJ/mol/A per Eh/bohr
PLATEAU_EH, PLATEAU_N = 2.5e-4, 50
TOL_START, TOL_DR = 0.01, 1.0                  # kcal/mol: T1's start = the climbing image's energy; draft Amendment 2 test
def text(p):
    if os.path.exists(p): return open(p, errors="replace").read()
    if os.path.exists(p + ".xz"): return lzma.open(p + ".xz", "rt", errors="replace").read()
    return None
def kabsch_rmsd(A, B):
    n = len(A)
    ca = [sum(p[k] for p in A) / n for k in range(3)]; cb = [sum(p[k] for p in B) / n for k in range(3)]
    a = [[p[k] - ca[k] for k in range(3)] for p in A]; b = [[p[k] - cb[k] for k in range(3)] for p in B]
    S = [[sum(a[i][r] * b[i][c] for i in range(n)) for c in range(3)] for r in range(3)]
    (xx, xy, xz), (yx, yy, yz), (zx, zy, zz) = S
    N = [[xx + yy + zz, yz - zy, zx - xz, xy - yx], [yz - zy, xx - yy - zz, xy + yx, zx + xz],
         [zx - xz, xy + yx, -xx + yy - zz, yz + zy], [xy - yx, zx + xz, yz + zy, -xx - yy + zz]]
    lam = max(abs(v) for row in N for v in row) * 4 + 1.0
    v = [1.0, 0.0, 0.0, 0.0]
    for _ in range(500):
        w = [sum((N[i][j] + (lam if i == j else 0.0)) * v[j] for j in range(4)) for i in range(4)]
        s = sum(x * x for x in w) ** 0.5; v = [x / s for x in w]
    emax = sum(v[i] * sum(N[i][j] * v[j] for j in range(4)) for i in range(4))
    ga = sum(x * x for p in a for x in p); gb = sum(x * x for p in b for x in p)
    return max(0.0, (ga + gb - 2.0 * emax) / n) ** 0.5
def echo(t):
    i = t.find("INPUT FILE"); j = t.find("****END OF INPUT****", i)
    if i < 0 or j < 0: return None
    e = [m.group(1).rstrip() for m in re.finditer(r"^\|\s*\d+> ?(.*)$", t[i:j], re.M)]
    while e and not e[-1]: e.pop()
    return e
def xyz(p):
    t = text(p)
    if t is None: return None
    L = t.splitlines(); n = int(L[0].split()[0])
    return [tuple(float(v) for v in l.split()[1:4]) for l in L[2:2 + n]]
def run_state(d, nact):
    """one optimisation folder: status, energies, checks"""
    r = {"dir": d, "notes": []}
    if not os.path.isfile(os.path.join(d, "reopt.inp")): r["status"] = "no inputs"; return r
    t = text(os.path.join(d, "reopt.out"))
    if t is None: r["status"] = "not started"; return r
    raw = open(os.path.join(d, "reopt.inp")).read()
    inp = [l.rstrip() for l in raw.splitlines()]
    while inp and not inp[-1]: inp.pop()
    e = echo(t)
    if e is not None and e != inp: r["notes"].append("PROBLEM: echoed input differs from reopt.inp")
    m = re.search(r"The system contains (\d+) active degrees of freedom \((\d+) active atoms", t)
    r["nact"] = int(m.group(2)) if m else None
    if m and int(m.group(2)) != nact: r["notes"].append("PROBLEM: ORCA reports %s active atoms, expected %d" % (m.group(2), nact))
    want = {k: float(re.search(k + r"\s+(\S+)", raw).group(1)) * EH_BOHR for k in ("TolMaxG", "TolRMSG")}
    tr = re.findall(r"RMS\(Grad\) Convergence:\s+(\d+\.\d+)", t); tm = re.findall(r"Max\(Grad\) Convergence:\s+(\d+\.\d+)", t)
    if tr and tm and (abs(float(tr[0]) - want["TolRMSG"]) > max(6e-4, 2e-3 * want["TolRMSG"])
                      or abs(float(tm[0]) - want["TolMaxG"]) > max(6e-4, 2e-3 * want["TolMaxG"])):
        r["notes"].append("PROBLEM: thresholds %s/%s kJ/mol/A" % (tm[0], tr[0]))
    r["thr"] = (tm[0], tr[0]) if tr and tm else None
    for pat in ("ABORTING THE RUN", "error termination"):
        if pat in t: r["notes"].append("PROBLEM: '%s' in reopt.out" % pat); break
    cp = os.path.join(d, "reopt-minimize-ener.csv")
    rows = [l.split(";") for l in open(cp) if l.strip() and not l.lstrip().startswith("#")] if os.path.exists(cp) else []
    E = [float(x[5]) for x in rows]
    r["E"] = E
    r["grad"] = "ORCA TERMINATED NORMALLY" in t and "The minimization has converged." in t
    r["slide"] = (E[-1 - PLATEAU_N] - E[-1]) if len(E) > PLATEAU_N else None
    r["plateau"] = r["slide"] is None or r["slide"] <= PLATEAU_EH
    po = text(os.path.join(d, "reopt.pbs.out")) or ""
    r["extract"] = "extract=ok" in po
    ended = "end=" in po
    if not E: r["status"] = "started, no optimiser step yet"
    elif not ended: r["status"] = "running: %d steps so far" % len(E)
    elif r["grad"] and r["plateau"] and r["extract"]: r["status"] = "finished: %d steps, gradients met, energy plateau, c4_extract ok" % len(E)
    else: r["status"] = "ended without finishing (gradients %s, plateau %s, extract %s): %d steps" % (r["grad"], r["plateau"], r["extract"], len(E))
    r["done"] = ended and r["grad"] and r["plateau"] and r["extract"]
    sp = re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", t)
    r["Esp"] = float(sp[-1]) if sp else None
    if r["done"]:
        nt = re.findall(r"N\(Total\)\s*:\s*(\d+\.\d+)\s+electrons", t)
        if not nt or round(float(nt[-1])) != 118: r["notes"].append("PROBLEM: electrons %s" % (nt[-1] if nt else "not found"))
    return r
BL = [l.rstrip("\n").split("\t") for l in open(os.path.join(D, "c4_bands.tsv")) if l.strip() and not l.startswith("#")]
B = {r[0]: dict(zip(BL[0], r)) for r in BL[1:]}[F]
E_ci, E_kl = float(B["E_ci"]), float(B["E_klow"])
ACT = [int(x) for x in open(os.path.join(D, "active_atoms.txt")).read().split()]
ENV = [i for i in ACT if not QM0 <= i < QM0 + NQM]
PD = os.path.join(D, "pilot_a2")
s2 = run_state(os.path.join(D, "S2_" + F), len(ACT))
t1 = run_state(os.path.join(PD, "T", "T1_" + F), len(ENV))
r2 = run_state(os.path.join(PD, "R", "R2_" + F), len(ACT))
out = []; say = out.append
say("C4 Amendment 2 pilot, frame %s (c4a2_pilot_report.py) - exploratory, not part of C4's result (pilot_a2/PILOT_A2.txt)" % F)
say("energies: the optimiser's (reopt-minimize-ener.csv), kcal/mol unless stated; band energies from c4_bands.tsv (5 decimals)")
say("band %s: climbing image %s, k_low %s; B = E_ci - E_klow = %.2f" % (B["band"], B["ci"], B["k_low"], (E_ci - E_kl) * H))
say("\nS2_%s (step 2, the reactant re-optimised from k_low): %s" % (F, s2["status"]))
if s2.get("E"):
    say("  E_min %.8f Eh; dE_reopt = E_min - E_klow = %.2f" % (s2["E"][-1], (s2["E"][-1] - E_kl) * H))
for n in s2["notes"]: say("  " + n)
say("\nT1_%s (environment relaxed around the climbing-image substrate; substrate fixed): %s" % (F, t1["status"]))
if t1.get("nact") is not None: say("  ORCA: %d active atoms (expected %d)%s" % (t1["nact"], len(ENV), "" if t1["nact"] == len(ENV) else "  <- PROBLEM"))
if t1.get("thr"): say("  ORCA thresholds (kJ/mol/A): max %s, RMS %s" % t1["thr"])
if t1.get("E"):
    d0 = (t1["E"][0] - E_ci) * H
    say("  E_start - E_ci = %+.4f (start = the converged climbing image; must be within %.2f) -> %s" % (d0, TOL_START, "ok" if abs(d0) <= TOL_START else "PROBLEM"))
    say("  fall_T = E_ci - E_min(T1) = %.2f%s" % ((E_ci - t1["E"][-1]) * H, "" if t1["done"] else " (so far)"))
    if t1["slide"] is not None: say("  energy fall over the last %d steps: %.2e Eh" % (PLATEAU_N, t1["slide"]))
if t1.get("done"):
    S = xyz(os.path.join(t1["dir"], "start.xyz")); Fx = xyz(os.path.join(t1["dir"], "reopt.xyz"))
    if S and Fx and len(S) == NATOMS and len(Fx) == NATOMS:
        sq = [(sum((Fx[i][k] - S[i][k]) ** 2 for k in range(3)), i) for i in ENV]; mx = max(sq)
        sub = max(max(abs(Fx[i][k] - S[i][k]) for k in range(3)) for i in range(QM0, QM0 + NQM))
        say("  environment displacement (start -> end, %d atoms): RMS %.3f A, largest %.2f A (atom %d%s)" % (
            len(ENV), (sum(x[0] for x in sq) / len(sq)) ** 0.5, mx[0] ** 0.5, mx[1], ", a water" if mx[1] >= 6342 else ""))
        say("  substrate displacement: largest coordinate change %.1e A (must be 0)%s" % (sub, "" if sub <= 2e-6 else "  <- PROBLEM"))
    if t1["Esp"] is not None: say("  final single point - last optimiser step = %+.4f" % ((t1["Esp"] - t1["E"][-1]) * H))
for n in t1["notes"]: say("  " + n)
say("\nR2_%s (back to the reactant from T1's environment; all 1959 atoms active): %s" % (F, r2["status"]))
if r2.get("nact") is not None: say("  ORCA: %d active atoms (expected %d)%s" % (r2["nact"], len(ACT), "" if r2["nact"] == len(ACT) else "  <- PROBLEM"))
if r2.get("E") and s2.get("E"):
    dr = (r2["E"][-1] - s2["E"][-1]) * H
    say("  dR = E_min(R2) - E_min(S2) = %+.2f%s" % (dr, "" if r2["done"] else " (so far)"))
    if r2["done"]: say("  draft Amendment 2 test |dR| <= %.1f: %s" % (TOL_DR, "met" if abs(dr) <= TOL_DR else "not met (a second cycle would follow)"))
if r2.get("done"):
    Fx = xyz(os.path.join(r2["dir"], "reopt.xyz"))
    md = xyz(os.path.join(D, "frames", "frame_" + F, "md_substrate.xyz"))
    if Fx and md and len(Fx) == NATOMS:
        say("  substrate RMSD to its MD structure: %.3f A (A2 limit 0.30)" % kabsch_rmsd([Fx[i] for i in range(QM0, QM0 + NQM)], md))
    if r2["Esp"] is not None: say("  final single point - last optimiser step = %+.4f" % ((r2["Esp"] - r2["E"][-1]) * H))
for n in r2["notes"]: say("  " + n)
say("\nsummary")
say("  B (band, k_low to climbing image)            %.2f" % ((E_ci - E_kl) * H))
if t1.get("done") and r2.get("done"):
    say("  B3 = E_min(T1) - E_min(R2)                   %.2f" % ((t1["E"][-1] - r2["E"][-1]) * H))
if t1.get("done") and s2.get("E"):
    say("  E_min(T1) - E_min(S2) (cross-lineage, check)  %.2f" % ((t1["E"][-1] - s2["E"][-1]) * H))
if not (t1.get("done") and r2.get("done")): say("  (B3 when both T1 and R2 have finished)")
print("\n".join(out))
