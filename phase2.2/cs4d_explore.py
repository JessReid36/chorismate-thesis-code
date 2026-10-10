#!/usr/bin/env python3
"""
cs4d_explore.py - EXPLORATORY analysis of Stage 4d (not pre-registered; written 10 October 2026 after
STAGE4D_REPORT.txt was seen). It changes no verdict. It looks for the sources of PC1's and PC2's failures:
  1. the enzyme's barrier change at frame 41786 split as Claeyssens et al. 2005 split theirs (electrostatic incl.
     polarisation; QM/MM van der Waals; the rest), and the frame's place in the Phase 1 ensemble;
  2. the steric (LJ) term of the frozen environment, by residue;
  3. the residues' hydrogen bonds to the substrate: reactant in env_R, TS in env_R (frozen), TS in env_TS (own);
  4. the first-order field (A_v2: V_TS - V_R per unit charge) at the design-grid sites nearest each residue's charge
     centre, in frame 41786 and across all 30 A_v2 frames (is a residue's sign specific to this frame?);
  5. the substrate's reactant and TS against the structures Szefczyk et al. 2004 analysed.
Reads the run outputs (job.out or job.out.xz), pc/*.pc(.xz), lj.tsv, residues.tsv and each run's substrate.xyz.
Python 3.6, standard library.
USAGE  python3 cs4d_explore.py <cs_stage4d folder> <phase2.2 folder> <ensemble_barriers.tsv>
       -> writes <cs_stage4d folder>/exploratory/STAGE4D_EXPLORATORY_REPORT.txt
"""
import lzma, math, os, re, statistics, sys
D, P2, ENS = sys.argv[1:4]
H = 627.5094740631
N = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
NAMES = {"ARG": "N H CA HA CB HB2 HB3 CG HG2 HG3 CD HD2 HD3 NE HE CZ NH1 HH11 HH12 NH2 HH21 HH22 C O".split(),
         "LYS": "N H CA HA CB HB2 HB3 CG HG2 HG3 CD HD2 HD3 CE HE2 HE3 NZ HZ1 HZ2 HZ3 C O".split(),
         "GLU": "N H CA HA CB HB2 HB3 CG HG2 HG3 CD OE1 OE2 C O".split()}       # Amber atom order
CENTRE = {"ARG": "CZ", "LYS": "NZ", "GLU": "CD"}
def opn(p):
    return open(p) if os.path.exists(p) else lzma.open(p + ".xz", "rt")
def energy(run):
    t = opn(os.path.join(D, run, "job.out")).read()
    if "ORCA TERMINATED NORMALLY" not in t: raise SystemExit("not terminated normally: " + run)
    return float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)[-1])
def sub(g):
    L = open(os.path.join(D, "bare_" + g, "substrate.xyz")).read().splitlines()
    return [tuple(float(v) for v in l.split()[1:4]) for l in L[2:26]]
def pc(name):
    L = opn(os.path.join(D, "pc", name + ".pc")).read().splitlines()
    return [tuple(float(v) for v in l.split()) for l in L[1:1 + int(L[0])]]
def dist(a, b): return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))
X = {"R": sub("R"), "TS": sub("TS")}
E = {r: energy(r) for r in ("bare_R", "bare_TS", "ENZ_R", "ENZ_TS", "ENZTS_R", "ENZTS_TS")}
LJ = {}
for l in open(os.path.join(D, "lj.tsv")):
    if not l.startswith("#") and l.strip():
        f = l.split("\t"); LJ[f[0]] = (float(f[2]), float(f[3]))
RES = [l.rstrip("\n").split("\t") for l in open(os.path.join(D, "residues.tsv")) if not l.startswith("#") and l.strip()]
out = []; say = out.append
say("Stage 4d - EXPLORATORY report (cs4d_explore.py; not pre-registered; verdicts unchanged)")
# ---- 1
qm_own = ((E["ENZTS_TS"] - E["bare_TS"]) - (E["ENZ_R"] - E["bare_R"])) * H
lj_own = LJ["envTS_all"][1] - LJ["envR_all"][0]
rows = [l.rstrip("\n").split("\t") for l in open(ENS) if not l.startswith("#")]
i = rows[0].index("stab_TS"); st = sorted((float(r[i]), r[0]) for r in rows[1:] if r[i].strip())
s41 = [v for v, f in st if f == "41786"][0]
say("\n1. Frame 41786's TS stabilisation by the enzyme, split as Claeyssens et al. 2005 split theirs (kcal/mol):")
say("   electrostatic incl. polarisation (QM part, R in env_R, TS in env_TS)   %+7.3f   (Claeyssens: 4.7 on average, 16 paths)" % qm_own)
say("   substrate-environment van der Waals (LJ, same geometries)              %+7.3f" % lj_own)
say("   the rest (Phase 1 stab_TS minus the two; mainly the protein's own MM energy change) %+7.3f" % (s41 - qm_own - lj_own))
say("   stab_TS (Phase 1)                                                       %+7.3f   (Claeyssens: 4.2 on average)" % s41)
say("   frame 41786 in the Phase 1 ensemble: %d of %d frames stabilise the TS less or equally (mean %+.2f, sd %.2f)" % (
    sum(1 for v, f in st if v >= s41), len(st), statistics.mean(v for v, f in st), statistics.stdev(v for v, f in st)))
# ---- 2
say("\n2. The frozen environment's steric term (LJ(TS) - LJ(R), substrate against the reactant's environment, kcal/mol):")
say("   whole environment %+7.3f; of which, by residue alone: %s" % (LJ["envR_all"][1] - LJ["envR_all"][0], ", ".join(
    "%s %+.3f" % (r[0], LJ["res_" + r[0]][1] - LJ["res_" + r[0]][0]) for r in RES)))
say("   environment of the TS (env_TS): %+7.3f" % (LJ["envTS_all"][1] - LJ["envTS_all"][0]))
# ---- 3
ENZ, ENZTS = pc("ENZ"), pc("ENZTS")
say("\n3. Hydrogen bonds between residues and substrate (A; H...O or O...H below 2.6 A in any of the three):")
say("   %-5s %-11s %8s %14s %12s" % ("", "contact", "R, env_R", "TS, env_R", "TS, env_TS"))
CEN = {}
for r in RES:
    lab, rn = r[0], r[2]; A = pc("ALONE_" + lab); nm = NAMES[rn]
    if len(nm) != len(A): raise SystemExit("atom count of %s" % lab)
    idx = []
    for a in A:
        j = min(range(len(ENZ)), key=lambda k: max(abs(ENZ[k][m] - a[m]) for m in (1, 2, 3)))
        if max(abs(ENZ[j][m] - a[m]) for m in (1, 2, 3)) > 1e-6: raise SystemExit("residue %s not found in ENZ.pc" % lab)
        idx.append(j)
    CEN[lab] = A[nm.index(CENTRE[rn])][1:]
    shown = 0
    for k, an in enumerate(nm):
        if an[0] not in "HO": continue
        pr, pt = A[k][1:], ENZTS[idx[k]][1:]
        for j, sn in enumerate(N):
            if (an[0] == "H" and sn[0] == "O") or (an[0] == "O" and sn[0] == "H"):
                d = (dist(pr, X["R"][j]), dist(pr, X["TS"][j]), dist(pt, X["TS"][j]))
                if min(d) < 2.6:
                    say("   %-5s %-11s %8.2f %14.2f %12.2f" % (lab, an + "-" + sn, d[0], d[1], d[2])); shown += 1
    if not shown: say("   %-5s none" % lab)
# ---- 4
T = [[float(v) for v in l.split()] for l in open(os.path.join(P2, "aligned/frame_41786/transform.txt")) if not l.startswith("#")]
Rm, t = T[:3], T[3]
G = [[float(v) for v in l.split("\t")[1:4]] for l in open(os.path.join(P2, "grid_v2.tsv")) if not l.startswith("#") and not l.startswith("idx")]
fr = [l for l in open(os.path.join(P2, "A_v2.tsv")) if l.startswith("# frames")][0].split()[3:]
A2 = [[float(x) for x in l.split("\t")] for l in open(os.path.join(P2, "A_v2.tsv")) if not l.startswith("#")]
c41 = fr.index("41786")
say("\n4. First-order field a = V_TS - V_R (kcal/mol per unit charge; a < 0: a cation there lowers the barrier) at the three")
say("   grid_v2 sites nearest each residue's charge centre; frame 41786, and over all %d A_v2 frames:" % len(fr))
for r in RES:
    lab = r[0]; c = CEN[lab]
    ca = [sum(c[j] * Rm[j][i] for j in range(3)) + t[i] for i in range(3)]        # aligned = original @ R + t
    near = sorted(range(len(G)), key=lambda s: dist(G[s], ca))[:3]
    for s in near:
        col = A2[s]
        say("   %-5s (%s %+d) site %3d, %.1f A: frame 41786 %+6.2f; all frames median %+6.2f, range %+6.2f to %+6.2f, negative in %d/%d"
            % (lab, r[2], round(float(r[4])), s, dist(G[s], ca), col[c41], statistics.median(col), min(col), max(col),
               sum(1 for v in col if v < 0), len(col)))
# ---- 5
def dd(g, a, b): return dist(X[g][N.index(a)], X[g][N.index(b)])
say("\n5. Substrate geometry: breaking C4-O3, forming C1-C6 (A). Szefczyk et al. 2004 (their numbering C-O3, C3-C9):")
say("   substrate 1.45, 3.44; TS 2.17, 2.58. Frame 41786: R %.3f, %.3f; TS %.3f, %.3f" % (
    dd("R", "C4", "O3"), dd("R", "C1", "C6"), dd("TS", "C4", "O3"), dd("TS", "C1", "C6")))
os.makedirs(os.path.join(D, "exploratory"), exist_ok=True)
open(os.path.join(D, "exploratory", "STAGE4D_EXPLORATORY_REPORT.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
