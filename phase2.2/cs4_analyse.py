#!/usr/bin/env python3
"""
cs4_analyse.py - Stage 4 (spill-out check): applies STAGE4_CRITERIA.txt to phase2.2/cs_stage4/.
Python 3.6, standard library only (hpc1 login node). Writes STAGE4_REPORT.txt.
USAGE  python3 cs4_analyse.py [cs_stage4 folder]
"""
import math, os, re, sys
D = sys.argv[1] if len(sys.argv) > 1 else "/home/18660916/system_development/phase2.2/cs_stage4"
H2KCAL = 627.5094740631
NAMES = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
OXY = ["O1", "O2", "O3", "O4", "O5", "O6"]
BONDS = [(0, 1), (0, 2), (0, 3), (3, 4), (4, 5), (4, 6), (3, 7), (7, 8), (8, 9), (8, 10), (10, 11), (10, 12),
         (12, 13), (13, 14), (13, 15), (15, 16), (15, 17), (17, 18), (17, 19), (19, 20), (17, 8), (12, 21),
         (21, 22), (21, 23)]
SVP_CONTACT = {"o1_p1_svpd": ("O2", 2.674), "o2_q1_svpd": ("O2", 2.711)}   # Stage 3 p1 and WP4 armb_anch, def2-SVP
def dist(a, b): return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))
def read_xyz(p):
    L = open(p).read().splitlines(); n = int(L[0].split()[0])
    return [[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]]
def text(p): return open(p, errors="replace").read() if os.path.exists(p) else ""
def energy(t):
    m = re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", t)
    return float(m[-1]) if m and "ORCA TERMINATED NORMALLY" in t else None
def charges(t, header, pat):
    i = t.rfind(header)
    if i < 0: return None
    out = {}
    for l in t[i:].splitlines()[1:200]:
        m = re.match(pat, l)
        if m: out[int(m.group(1))] = float(m.group(3))
        elif out and l.strip() and not l.strip().startswith("-"): 
            if len(out) >= 24: break
    return out if len(out) >= 24 else None
HIRSH = r"^\s*(\d+)\s+([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$"
MULL = r"^\s*(\d+)\s+([A-Z][a-z]?)\s*:\s*(-?\d+\.\d+)"
rep = []; say = rep.append
say("Stage 4 report (cs4_analyse.py) - criteria from STAGE4_CRITERIA.txt")
ok1 = True
for run in ("o1_p1_svpd", "o2_q1_svpd"):
    F = os.path.join(D, run); t = text(os.path.join(F, "job.out"))
    conv = "THE OPTIMIZATION HAS CONVERGED" in t
    if not os.path.exists(os.path.join(F, "job.xyz")):
        say("== %s: no final geometry -> S4.1 not assessable" % run); ok1 = False; continue
    X0, X1 = read_xyz(os.path.join(F, "combined.xyz")), read_xyz(os.path.join(F, "job.xyz"))
    s = X1[24]
    do = sorted((dist(s, X1[NAMES.index(o)]), o) for o in OXY)
    da = min(dist(s, X1[j]) for j in range(24))
    db = max(abs(dist(X1[a], X1[b]) - dist(X0[a], X0[b])) for a, b in BONDS)
    ref_o, ref_d = SVP_CONTACT[run]
    shift = do[0][0] - ref_d
    good = conv and do[0][1] == ref_o and abs(shift) < 0.05 and da >= 2.0
    ok1 = ok1 and good
    say("== %s (%d cycles, %s): nearest O %s %.3f A with def2-SVPD vs %s %.3f A with def2-SVP: shift %+.3f A; "
        "closest atom %.3f A; largest bond change from the SVP minimum %.3f A -> %s" % (
            run, t.count("GEOMETRY OPTIMIZATION CYCLE"), "converged" if conv else "NOT converged", do[0][1], do[0][0],
            ref_o, ref_d, shift, da, db, "within 0.05 A" if good else "NOT within 0.05 A (or other condition failed)"))
say("S4.1 contact distance moves less than 0.05 A with diffuse functions, same nearest oxygen, both runs -> %s" % ("PASS" if ok1 else "FAIL"))
say("")
res = {}
for tag in ("svp", "svpd"):
    E, Hq, Mq = {}, {}, {}
    for qt, q in (("q0", 0.0), ("qp05", 0.5), ("qm05", -0.5), ("qp1", 1.0), ("qm1", -1.0)):
        t = text(os.path.join(D, "sp_%s_%s" % (tag, qt), "job.out"))
        E[q] = energy(t)
        Hq[q] = charges(t, "HIRSHFELD ANALYSIS", HIRSH)
        Mq[q] = charges(t, "MULLIKEN ATOMIC CHARGES", MULL)
    if any(v is None for v in E.values()):
        say("%s: a single point is missing or did not terminate normally" % tag); res[tag] = None; continue
    D1 = (E[1.0] - E[-1.0]) * H2KCAL; D5 = (E[0.5] - E[-0.5]) * H2KCAL
    S1 = (E[1.0] + E[-1.0] - 2 * E[0.0]) * H2KCAL; S5 = (E[0.5] + E[-0.5] - 2 * E[0.0]) * H2KCAL
    a = (8 * D5 - D1) / 6; c = (D1 - 2 * a) / 2
    b = (16 * S5 - S1) / 6; d4 = (S1 - 2 * b) / 2
    res[tag] = (a, b, c, d4, Hq, Mq, E)
    say("%-5s a = %10.4f  b = %8.4f  (q^3 %+.4f, q^4 %+.4f) kcal/mol" % ("def2-" + tag.upper().replace("SVPD", "SVPD"), a, b, c, d4))
if res.get("svp") and res.get("svpd"):
    a0, b0 = res["svp"][0], res["svp"][1]; a1, b1 = res["svpd"][0], res["svpd"][1]
    rel = abs(b1 - b0) / abs(b0)
    say("S4.2 polarisation term at contact: b(SVPD) %.4f vs b(SVP) %.4f kcal/mol, change %.1f%% (criterion < 10%%) -> %s" % (
        b1, b0, 100 * rel, "PASS" if rel < 0.10 else "FAIL"))
    say("informative: first-order a changes %+.4f kcal/mol (%.2f%%) with the basis (it is the potential of a different density)" % (
        a1 - a0, 100 * (a1 - a0) / abs(a0)))
    for lab, k in (("Hirshfeld", 4), ("Mulliken", 5)):
        q0s, q1s = res["svp"][k][0.0], res["svp"][k][1.0]
        q0d, q1d = res["svpd"][k][0.0], res["svpd"][k][1.0]
        if q0s and q1s and q0d and q1d:
            o2 = NAMES.index("O2")
            say("informative %s charge on O2 (nearest the sphere): q=0 -> q=+1 change %+.4f e (SVP) vs %+.4f e (SVPD)" % (
                lab, q1s[o2] - q0s[o2], q1d[o2] - q0d[o2]))
        else:
            say("informative %s charges: not found in every output" % lab)
    p1E = re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", text(os.path.join(D, "..", "cs_stage3", "p1_41786_qp05", "job.out")))
    if p1E:
        say("check: def2-SVP q=+0.5 single point vs Stage 3 p1's final energy: %.1e Eh" % (res["svp"][6][0.5] - float(p1E[-1])))
else:
    say("S4.2 -> FAIL (single points missing)")
open(os.path.join(D, "STAGE4_REPORT.txt"), "w").write("\n".join(rep) + "\n")
print("\n".join(rep))
