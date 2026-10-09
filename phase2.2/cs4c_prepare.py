#!/usr/bin/env python3
"""
cs4c_prepare.py - Stage 4c of the charged-sphere development (pre-registered; STAGE4C_CRITERIA.txt): does the basis
set change the barrier change a sphere arrangement produces, or the ranking of arrangements?
Built from committed files only:
  aligned/frame_41786/{reactant_qm,transition_state_qm}.xyz  frame 41786's R and TS (aligned frame; identical
                                                             internal geometry to Stage 1's)
  grid_v2.tsv, A_v2.tsv                                      the design grid and A matrix (column for frame 41786)
  cs_stage4/sp_svp_q0/{job.inp,site.ORCAFF.prms}              the production QM/MM input and force-field file
Eight arrangements, chosen by a stated rule from frame 41786's a (kcal/mol per unit charge; a < 0 means a positive
charge there lowers the barrier), sites at least 3.25 A apart (the approach document's sphere spacing):
  E1 most negative a, q = +1        E2 the same site, q = +0.5     E3 most positive a, q = -1
  E4 the site nearest a carboxylate oxygen (O1, O2, O5, O6), q = +1
  E5 three sites taken greedily by |a|, q = -sign(a) x 1.0        E6 five sites, the same rule, x 0.5
  E7 E5 with every sign reversed      E8 three sites taken greedily from the most negative a, q = +1 each
Levels (every arrangement at R and TS):
  L1 production: ORCA QM/MM, B3LYP-D3BJ/def2-SVP, spheres as MM atoms with Amber N LJ (as Stages 1-4)
  L0, L2, L3: plain ORCA, spheres as coreless ECP centres Ne> carrying q and the Ne-type pseudopotential (Marefat
      Khah et al. 2020; given inline, which Stage 4b showed equal to the %basis form), def2-SVP / def2-SVPD /
      def2-TZVPD; bare substrate (no spheres) at R and TS in each basis as the reference.
USAGE  python3 cs4c_prepare.py <results repo> <output folder>
"""
import os, sys
import numpy as np
REPO, OUT = sys.argv[1], sys.argv[2]
P2 = os.path.join(REPO, "phase2.2")
SD = "/home/18660916/system_development/phase2.2/cs_stage4c"
def xyz(p):
    L = open(p).read().splitlines(); n = int(L[0])
    return [l.split()[0] for l in L[2:2 + n]], np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])
LAB, GEOM = None, {}
for g, f in (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz")):
    lab, X = xyz(os.path.join(P2, "aligned/frame_41786", f)); assert len(lab) == 24
    LAB = LAB or lab; assert lab == LAB; GEOM[g] = X
fr = [l for l in open(os.path.join(P2, "A_v2.tsv")) if l.startswith("# frames")][0].split()[3:]
A = np.array([[float(x) for x in l.split("\t")] for l in open(os.path.join(P2, "A_v2.tsv")) if not l.startswith("#")])
a = A[:, fr.index("41786")]
G = np.array([[float(v) for v in l.split("\t")[1:4]] for l in open(os.path.join(P2, "grid_v2.tsv"))
              if not l.startswith("#") and not l.startswith("idx")])
assert G.shape == (233, 3) and A.shape[0] == 233
def greedy(order, n):
    pick = []
    for i in order:
        if all(np.linalg.norm(G[i] - G[j]) >= 3.25 for j in pick): pick.append(int(i))
        if len(pick) == n: return pick
    raise SystemExit("greedy selection could not place %d sites" % n)
i1, i3 = int(a.argmin()), int(a.argmax())
N = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
carb = [N.index(x) for x in ("O1", "O2", "O5", "O6")]
dR = np.linalg.norm(G[:, None] - GEOM["R"][None], axis=2)
i4 = int(dR[:, carb].min(axis=1).argmin())
o = np.argsort(-abs(a), kind="stable")
s5, s6, s8 = greedy(o, 3), greedy(o, 5), greedy(np.argsort(a, kind="stable"), 3)
ENV = {"E1": [(i1, 1.0)], "E2": [(i1, 0.5)], "E3": [(i3, -1.0)], "E4": [(i4, 1.0)],
       "E5": [(i, -float(np.sign(a[i]))) for i in s5], "E6": [(i, -0.5 * float(np.sign(a[i]))) for i in s6],
       "E7": [(i, float(np.sign(a[i]))) for i in s5], "E8": [(i, 1.0) for i in s8]}
for e, S in ENV.items():
    for i, _ in S: assert G[i].shape == (3,)
    for x in range(len(S)):
        for y in range(x + 1, len(S)): assert np.linalg.norm(G[S[x][0]] - G[S[y][0]]) >= 3.25
SPIN = open(os.path.join(P2, "cs_stage4/sp_svp_q0/job.inp")).read()
assert SPIN.startswith("! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF QMMM")
PR = open(os.path.join(P2, "cs_stage4/sp_svp_q0/site.ORCAFF.prms")).read().splitlines()
ia = PR.index("$atoms"); ib = PR.index("$bonds")
SUBROWS = PR[ia + 2: ia + 2 + 24]; TAIL = PR[ib:]
HEAD = """! B3LYP D3BJ %s def2/J RIJCOSX TightSCF
%%maxcore 3000
%%pal nprocs 8 end
%%scf MaxIter 300 end
"""
ECP = """  N_core 10
  lmax f
  s 1
    1  2.0475   54.5100  2
  p 1
    1  0.44815   1.46500  2
  d 1
    1  0.49205  -0.83800  2
  f 1
    1  1.0000    0.0000  2
end
"""
os.makedirs(OUT)
batches = {"svp": [], "svpd": [], "tzvpd_a": [], "tzvpd_b": []}
def sub(g):
    return "".join("  %-2s %15.8f%15.8f%15.8f\n" % (LAB[k], GEOM[g][k, 0], GEOM[g][k, 1], GEOM[g][k, 2]) for k in range(24))
def plain(name, basis, g, spheres, batch):
    d = os.path.join(OUT, name); os.makedirs(d)
    s = HEAD % basis + "* xyz -2 1\n" + sub(g)
    for i, q in spheres:
        s += "  Ne>  %.6f %15.8f%15.8f%15.8f  NewECP\n" % (q, G[i, 0], G[i, 1], G[i, 2]) + ECP
    open(os.path.join(d, "job.inp"), "w").write(s + "*\n"); batches[batch].append(name)
def qmmm(name, g, spheres, batch):
    d = os.path.join(OUT, name); os.makedirs(d)
    with open(os.path.join(d, "combined.xyz"), "w") as f:
        f.write("%d\nframe 41786 %s (aligned), Stage 4c %s\n" % (24 + len(spheres), g, name))
        f.write("".join("%-2s %15.8f%15.8f%15.8f\n" % (LAB[k], GEOM[g][k, 0], GEOM[g][k, 1], GEOM[g][k, 2]) for k in range(24)))
        f.write("".join("%-2s %15.8f%15.8f%15.8f\n" % ("N", G[i, 0], G[i, 1], G[i, 2]) for i, _ in spheres))
    rows = SUBROWS + ["%6d   %-2s%13.6f%13.6f%13.6f%13.6f%13.6f" % (25 + k, "N", q, -0.17, 3.648, -0.085, 3.648)
                      for k, (_, q) in enumerate(spheres)]
    open(os.path.join(d, "site.ORCAFF.prms"), "w").write("\n".join(PR[:ia + 1] + ["%d 1 4" % (24 + len(spheres))] + rows + TAIL) + "\n")
    open(os.path.join(d, "job.inp"), "w").write(SPIN); batches[batch].append(name)
for g in ("R", "TS"):
    for basis, bt, batch in (("def2-SVP", "svp", "svp"), ("def2-SVPD", "svpd", "svpd"), ("def2-TZVPD", "tzvpd", "tzvpd_a")):
        plain("bare_%s_%s" % (bt, g), basis, g, [], batch)
for k, (e, S) in enumerate(ENV.items()):
    for g in ("R", "TS"):
        qmmm("L1_%s_%s" % (e, g), g, S, "svp")
        plain("L0_%s_%s" % (e, g), "def2-SVP", g, S, "svp")
        plain("L2_%s_%s" % (e, g), "def2-SVPD", g, S, "svpd")
        plain("L3_%s_%s" % (e, g), "def2-TZVPD", g, S, "tzvpd_a" if k < 4 else "tzvpd_b")
for b, runs in batches.items():
    with open(os.path.join(OUT, "batch_%s.pbs" % b), "w") as f:
        f.write("""#!/bin/bash
#PBS -N cm34c_%s
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -l walltime=168:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %s/batch_%s.pbs.out
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
echo "host=$(hostname) start=$(date)"
""" % (b, SD, b))
        f.write("for d in %s; do\n  cd %s/$d\n  /home/apps2/ORCA/6.0.1/orca job.inp > job.out 2>&1 </dev/null\n" % (" ".join(runs), SD))
        f.write('  echo "$d orca_exit=$? $(grep -c \'ORCA TERMINATED NORMALLY\' job.out) $(date)"\n  rm -f job.gbw job*.tmp 2>/dev/null\ndone\necho "end=$(date)"\n')
with open(os.path.join(OUT, "arrangements.tsv"), "w") as f:
    f.write("# arrangement\tsite (grid_v2 row)\tq\ta_41786 (kcal/mol per e)\tnearest atom at R (A)\tx\ty\tz\n")
    for e, S in ENV.items():
        for i, q in S:
            j = int(dR[i].argmin())
            f.write("%s\t%d\t%+.2f\t%.6f\t%s %.3f\t%.8f\t%.8f\t%.8f\n" % (e, i, q, a[i], N[j], dR[i, j], G[i, 0], G[i, 1], G[i, 2]))
    for e, S in ENV.items():
        f.write("# %s first-order prediction sum(a q) = %+.4f kcal/mol\n" % (e, sum(a[i] * q for i, q in S)))
print({b: len(r) for b, r in batches.items()}, "runs;", sum(len(r) for r in batches.values()), "total")
