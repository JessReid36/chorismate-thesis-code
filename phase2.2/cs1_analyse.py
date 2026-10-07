#!/usr/bin/env python3
"""
cs1_analyse.py - Stage 1 (charged-sphere implementation checks): applies STAGE1_CRITERIA.txt to the seven
single points in phase2.2/cs_stage1/. Python 3.6, standard library only (hpc1 login node).
Writes STAGE1_REPORT.txt into the folder. USAGE  python3 cs1_analyse.py [cs_stage1 folder]
"""
import math, os, re, sys
D = sys.argv[1] if len(sys.argv) > 1 else "/home/18660916/system_development/phase2.2/cs_stage1"
H2K = 627.508107          # ORCA's energy constant, measured in WP1
H2KCAL = 627.5094740631   # project constant for reporting
RUNS = ["s1a_q1_N", "s1b_q1_C", "s1c_q1_Ar", "s1d_q0_N", "s1e_q05_N", "s1f_two_R", "s1g_two_TS"]
WP4_CYCLE1 = -836.340705524819   # WP4 armb_anch, first single point, same geometry and force field (8 processes)

def parse(run):
    p = os.path.join(D, run, "job.out")
    if not os.path.exists(p):
        return None
    t = open(p, errors="replace").read()
    r = {"normal": "ORCA TERMINATED NORMALLY" in t}
    m = re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", t)
    r["E_str"] = m[-1] if m else None
    r["E"] = float(m[-1]) if m else None
    r["LJ_str"] = (re.findall(r"LJ interaction\s+(-?\d+\.\d+)", t) or [None])[-1]
    r["CO_str"] = (re.findall(r"Coulomb interaction\s+(-?\d+\.\d+)", t) or [None])[-1]
    r["charge"] = (re.findall(r"Total Charge\s+Charge\s+\.\.\.\.\s+(-?\d+)", t) or [None])[-1]
    r["mult"] = (re.findall(r"Multiplicity\s+Mult\s+\.\.\.\.\s+(\d+)", t) or [None])[-1]
    g = []
    i = t.rfind("CARTESIAN GRADIENT (QM/MM)")
    if i >= 0:
        for l in t[i:].splitlines()[1:]:
            mm = re.match(r"^\s*(\d+)\s+(\w+)\s*:\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", l)
            if mm:
                g.append(l.strip())
            elif g:
                break
    r["grad_qm"] = g[:24]
    return r

def geometry(run):
    L = open(os.path.join(D, run, "combined.xyz")).read().splitlines()
    n = int(L[0])
    return [[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]]

def ff(run):
    L = open(os.path.join(D, run, "site.ORCAFF.prms")).read().splitlines()
    i = L.index("$atoms")
    n = int(L[i + 1].split()[0])
    return [(float(l.split()[2]), abs(float(l.split()[3])), float(l.split()[4])) for l in L[i + 2:i + 2 + n]]

def lj_qm_mm(run, ron=10.0, roff=12.0):
    X, P = geometry(run), ff(run)
    E = 0.0
    for j in range(24, len(X)):
        for i in range(24):
            e = math.sqrt(P[i][1] * P[j][1])
            if e == 0.0:
                continue
            Rij = (P[i][2] + P[j][2]) / 2.0
            r = math.sqrt(sum((X[i][k] - X[j][k]) ** 2 for k in range(3)))
            A, B = e * Rij ** 12, 2.0 * e * Rij ** 6
            E += A / r ** 12 - B / r ** 6
            if r < ron:
                E += -A / (ron * roff) ** 6 + B / (ron * roff) ** 3
    return E / H2K

R = {k: parse(k) for k in RUNS}
rep = []
say = rep.append
say("Stage 1 report (cs1_analyse.py) - criteria from STAGE1_CRITERIA.txt")
for k in RUNS:
    r = R[k]
    if r is None:
        say("  %-11s no output" % k); continue
    say("  %-11s %s  E(QM/MM) %s  LJ %s  Coulomb(MM-MM) %s  QM charge %s mult %s" % (
        k, "normal" if r["normal"] else "NOT NORMAL", r["E_str"], r["LJ_str"], r["CO_str"], r["charge"], r["mult"]))
say("")
# S1.1 element label
ok_runs = [k for k in ("s1a_q1_N", "s1b_q1_C", "s1c_q1_Ar") if R[k] and R[k]["normal"] and R[k]["E_str"]]
same_E = len({R[k]["E_str"] for k in ok_runs}) == 1
same_G = len({tuple(R[k]["grad_qm"]) for k in ok_runs}) == 1 and all(len(R[k]["grad_qm"]) == 24 for k in ok_runs)
same_LJ = len({R[k]["LJ_str"] for k in ok_runs}) == 1
say("S1.1 element label (runs compared: %s): energy %s, QM-atom gradient %s, LJ %s -> %s" % (
    ", ".join(ok_runs), "identical" if same_E else "DIFFERS", "identical" if same_G else "DIFFERS",
    "identical" if same_LJ else "DIFFERS", "PASS" if (len(ok_runs) >= 2 and same_E and same_G and same_LJ) else "FAIL"))
if "s1c_q1_Ar" not in ok_runs:
    say("     note: the Ar-labelled run did not complete normally; see its job.out")
# S1.2 fractional charge
r0, r5, r1 = R["s1d_q0_N"], R["s1e_q05_N"], R["s1a_q1_N"]
if all(x and x["normal"] for x in (r0, r5, r1)):
    d1 = (r1["E"] - r0["E"]) * H2KCAL
    d5 = (r5["E"] - r0["E"]) * H2KCAL
    b = 2.0 * (d1 - 2.0 * d5)
    a = d1 - b
    acc = r5["charge"] == "-2" and r5["mult"] == "1"
    say("S1.2 fractional charge: q=+0.5 run %s (QM charge %s, multiplicity %s); E(q)-E(0) = a q + b q^2 with "
        "a = %.3f, b = %.3f kcal/mol -> %s" % ("accepted" if acc else "NOT as expected", r5["charge"], r5["mult"],
                                               a, b, "PASS" if (acc and b < 0) else "FAIL"))
else:
    say("S1.2 fractional charge: a run did not complete normally -> FAIL")
# S1.3 frozen sphere-sphere terms
f, g = R["s1f_two_R"], R["s1g_two_TS"]
if f and g and f["normal"] and g["normal"]:
    mmlj_f = float(f["LJ_str"]) - lj_qm_mm("s1f_two_R")
    mmlj_g = float(g["LJ_str"]) - lj_qm_mm("s1g_two_TS")
    co_same = f["CO_str"] == g["CO_str"]
    lj_same = abs(mmlj_f - mmlj_g) <= 2e-6
    say("S1.3 frozen sphere-sphere terms: MM-MM Coulomb R %s / TS %s (%s); MM-MM LJ (ORCA LJ minus QM-MM LJ) "
        "R %.8f / TS %.8f Eh (difference %.1e) -> %s" % (f["CO_str"], g["CO_str"], "identical" if co_same else "DIFFERS",
                                                         mmlj_f, mmlj_g, abs(mmlj_f - mmlj_g),
                                                         "PASS" if (co_same and lj_same and float(f["CO_str"]) != 0.0) else "FAIL"))
    if float(f["CO_str"]) == 0.0:
        say("     note: the MM-MM Coulomb term is zero, so the test is not informative")
else:
    say("S1.3 frozen sphere-sphere terms: a run did not complete normally -> FAIL")
# S1.4 QM-MM Lennard-Jones reproduced by the rigid-model formula
dev = []
for k in ("s1a_q1_N", "s1b_q1_C", "s1c_q1_Ar", "s1d_q0_N", "s1e_q05_N"):
    if R[k] and R[k]["normal"]:
        dev.append((k, float(R[k]["LJ_str"]) - lj_qm_mm(k)))
mx = max(abs(d) for k, d in dev) if dev else float("nan")
say("S1.4 QM-MM Lennard-Jones from the force-field file (Lorentz-Berthelot on R_min, force-switch shift) vs ORCA: "
    "max |difference| %.1e Eh over %d runs -> %s" % (mx, len(dev), "PASS" if dev and mx <= 1e-6 else "FAIL"))
# informative
if r1 and r1["E"] is not None:
    say("info: s1a vs WP4 arm b's first single point (same geometry and force field; 4 vs 8 processes): %.1e Eh"
        % (r1["E"] - WP4_CYCLE1))
open(os.path.join(D, "STAGE1_REPORT.txt"), "w").write("\n".join(rep) + "\n")
print("\n".join(rep))
