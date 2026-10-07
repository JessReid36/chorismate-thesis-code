#!/usr/bin/env python3
"""
cs3_analyse.py - Stage 3: applies STAGE3_CRITERIA.txt to the four QM/MM optimisations in phase2.2/cs_stage3/.
Reads per run: combined.xyz (start), job.xyz (final), job.out (convergence), job_trj.xyz (trajectory, for the
extremes along the path and the stability rule if a run stopped early). Writes STAGE3_REPORT.txt.
Python 3.6, standard library only (hpc1 login node). USAGE  python3 cs3_analyse.py [cs_stage3 folder]
"""
import math, os, sys
D = sys.argv[1] if len(sys.argv) > 1 else "/home/18660916/system_development/phase2.2/cs_stage3"
NAMES = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
BONDS = [(0, 1), (0, 2), (0, 3), (3, 4), (4, 5), (4, 6), (3, 7), (7, 8), (8, 9), (8, 10), (10, 11), (10, 12),
         (12, 13), (13, 14), (13, 15), (15, 16), (15, 17), (17, 18), (17, 19), (19, 20), (17, 8), (12, 21),
         (21, 22), (21, 23)]
OXY = {"O1": "carbA", "O2": "carbA", "O5": "carbB", "O6": "carbB", "O3": "ether", "O4": "hydroxyl"}
RIGID_QP05 = {"carbA": 2.7441, "carbB": 2.9696}   # Stage 2, frame 41786 R, N sphere, q +0.5, own basin (cs2_contacts.tsv)
WP4_QP1 = 2.711                                    # WP4 armb_anch, same start and set-up at q +1: nearest O2
O4, H8 = NAMES.index("O4"), NAMES.index("H8")
def dist(a, b): return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))
def read_xyz(p):
    L = open(p).read().splitlines(); n = int(L[0].split()[0])
    return [[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]]
def read_traj(p):
    if not os.path.exists(p): return []
    L, out, i = open(p).read().splitlines(), [], 0
    while i < len(L):
        if not L[i].strip(): i += 1; continue
        n = int(L[i].split()[0]); out.append([[float(v) for v in l.split()[1:4]] for l in L[i + 2:i + 2 + n]]); i += n + 2
    return out
def status(p):
    t = open(p, errors="replace").read()
    return "THE OPTIMIZATION HAS CONVERGED" in t, "ORCA TERMINATED NORMALLY" in t, t.count("GEOMETRY OPTIMIZATION CYCLE")
rep = []; say = rep.append
say("Stage 3 report (cs3_analyse.py) - criteria from STAGE3_CRITERIA.txt")
for run in ("p1_41786_qp05", "n1_41786_qm1", "n2_41786_qm05", "n3_55446_qm1"):
    F = os.path.join(D, run)
    if not os.path.exists(os.path.join(F, "job.xyz")):
        say("\n== %s: no final geometry -> not assessable" % run); continue
    X0, X1 = read_xyz(os.path.join(F, "combined.xyz")), read_xyz(os.path.join(F, "job.xyz"))
    conv, term, ncyc = status(os.path.join(F, "job.out"))
    T = read_traj(os.path.join(F, "job_trj.xyz"))
    s = X1[24]
    say("\n== %s  (%d cycles, %s)" % (run, ncyc, "converged" if conv else "NOT converged"))
    held = max(dist(X0[k], X1[k]) for k in (10, 12, 13, 15, 24))
    say("   held atoms (C5 C6 C7 C8, sphere) moved at most %.4f A" % held)
    other = [b for b in BONDS if b != (19, 20)] if run.startswith("n") else BONDS
    wb = max(other, key=lambda b: abs(dist(X1[b[0]], X1[b[1]]) - dist(X0[b[0]], X0[b[1]])))
    dbond = abs(dist(X1[wb[0]], X1[wb[1]]) - dist(X0[wb[0]], X0[wb[1]]))
    integ = dbond <= 0.15
    say("   integrity: largest change %s-%s %.3f A -> %s" % (NAMES[wb[0]], NAMES[wb[1]], dbond, "PASS" if integ else "FAIL"))
    stable = True
    if not conv:
        key = (lambda X: min(dist(X[24], X[j]) for j in range(24))) if run.startswith("p") else (lambda X: dist(X[24], X[H8]))
        tail = [key(X) for X in T[-max(1, len(T) // 3):]] if len(T) >= 3 else []
        stable = bool(tail) and max(tail) - min(tail) <= 0.05
        say("   stopped early: contact moved %s over the last third -> stability rule %s" % (
            "%.3f A" % (max(tail) - min(tail)) if tail else "n/a", "MET" if stable else "NOT MET"))
    if run.startswith("p"):
        do = sorted((dist(s, X1[NAMES.index(o)]), o) for o in OXY)
        da = min((dist(s, X1[j]), NAMES[j]) for j in range(24))
        g = OXY[do[0][1]]
        ref = RIGID_QP05.get(g)
        say("   sphere to oxygens: " + ", ".join("%s %.3f" % (o, x) for x, o in do) + "; closest atom %s %.3f A" % (da[1], da[0]))
        ok_c = ref is not None and abs(do[0][0] - ref) <= 0.3
        say("   S3.1 contact: nearest O %s (%s) %.3f A vs rigid own-basin %s -> %s" % (
            do[0][1], g, do[0][0], "%.3f A (difference %+.3f)" % (ref, do[0][0] - ref) if ref else "none (rigid model holds no sphere there)",
            "within 0.3 A" if ok_c else "OUTSIDE"))
        nc = da[0] >= 2.0
        say("   S3.1 no collapse (closest atom >= 2.0 A): %s" % ("yes" if nc else "NO"))
        say("   S3.1 -> %s" % ("PASS" if (ok_c and nc and integ and (conv or stable)) else "FAIL"))
        say("   informative: q +0.5 contact %.3f A vs q +1 (WP4 armb_anch, same start and set-up) %.3f A: %+.3f A" % (
            do[0][0], WP4_QP1, do[0][0] - WP4_QP1))
    else:
        oh, hs = dist(X1[O4], X1[H8]), dist(s, X1[H8])
        oh0, hs0 = dist(X0[O4], X0[H8]), dist(s, X0[H8])
        tmin = min(dist(X[24], X[H8]) for X in T) if T else float("nan")
        tmax = max(dist(X[O4], X[H8]) for X in T) if T else float("nan")
        say("   O4-H8 %.3f -> %.3f A (largest along the path %.3f); H8..sphere %.3f -> %.3f A (smallest along the path %.3f)" % (
            oh0, oh, tmax, hs0, hs, tmin))
        say("   H8 to O2 %.3f -> %.3f A; sphere to its closest heavy atom %.3f A" % (
            dist(X0[H8], X0[6]), dist(X1[H8], X1[6]), min(dist(s, X1[j]) for j in range(24) if NAMES[j][0] != "H")))
        ok = oh <= 1.20 and hs >= 1.50
        say("   S3.2 no H8 artefact (O4-H8 <= 1.20 A and H8..sphere >= 1.50 A, final geometry) -> %s" % (
            "PASS" if (ok and integ and (conv or stable)) else "FAIL"))
open(os.path.join(D, "STAGE3_REPORT.txt"), "w").write("\n".join(rep) + "\n")
print("\n".join(rep))
