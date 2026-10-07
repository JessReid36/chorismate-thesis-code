#!/usr/bin/env python3
"""
cs3_prepare.py - Stage 3 of the charged-sphere development: QM/MM confirmation at selected charges.
Four QM/MM optimisations, built from committed files only, with exactly WP4 arm b's anchored set-up
(phase2.2/wp4_s13_rerun/armb_anch/job.inp: B3LYP-D3BJ/def2-SVP RIJCOSX TightSCF QMMM Opt; C5 C6 C7 C8 and the
sphere held; every other substrate atom free) and Amber N spheres (R_min 3.648 A full, eps 0.170 kcal/mol):
  p1_41786_qp05  frame 41786 R, q = +0.5, sphere at WP4 arm b's start (40.832, 18.750, 42.551)
                 (q = +1 at this start and set-up is WP4 armb_anch itself: nearest O2 2.711 A)
  n1_41786_qm1   frame 41786 R, q = -1,   sphere at the closest approach to H8 that every other atom's
  n2_41786_qm05  frame 41786 R, q = -0.5, Lennard-Jones wall allows (distance >= R_ij = (R_sphere + R_j)/2 to
  n3_55446_qm1   frame 55446 R, q = -1,   every atom with LJ; H8 has none), found exactly by excluded intervals
                                          along rays from H8 (20,000 directions)
Substrate geometries: 05_qmmm/20_invacuo/<frame>_R.xyz (the QM/MM-optimised reactants; 41786's equals Stage 1's).
Substrate force-field rows: phase2.2/cs_stage1/s1a_q1_N/site.ORCAFF.prms (identical to WP4 arm b's).
USAGE  python3 cs3_prepare.py <results repo> <output folder>
"""
import os, sys
import numpy as np
REPO, OUT = sys.argv[1], sys.argv[2]
NAMES = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
P = open(os.path.join(REPO, "phase2.2/cs_stage1/s1a_q1_N/site.ORCAFF.prms")).read().splitlines()
i = P.index("$atoms")
ROWS = P[i + 2:i + 26]
RC = np.array([float(l.split()[4]) for l in ROWS])
RS = 3.648
INP = open(os.path.join(REPO, "phase2.2/wp4_s13_rerun/armb_anch/job.inp")).read()
assert "ActiveAtoms {0 1 2 3 4 5 6 7 8 9 11 14 16 17 18 19 20 21 22 23}" in INP and "QMMM Opt" in INP
PBS = """#!/bin/bash
#PBS -N cm33_{name}
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -l walltime=168:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o /home/18660916/system_development/phase2.2/cs_stage3/{name}/job.pbs.out
cd /home/18660916/system_development/phase2.2/cs_stage3/{name}
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${{LD_LIBRARY_PATH:-}}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
echo "host=$(hostname) start=$(date)"
/home/apps2/ORCA/6.0.1/orca job.inp > job.out 2>&1 </dev/null
echo "orca_exit=$? end=$(date)"
grep -c 'GEOMETRY OPTIMIZATION CYCLE' job.out; grep -E 'THE OPTIMIZATION HAS CONVERGED|ORCA TERMINATED NORMALLY' job.out | tail -2
rm -f job.gbw job*.tmp 2>/dev/null
"""
def substrate(frame):
    L = open(os.path.join(REPO, "05_qmmm/20_invacuo/%s_R.xyz" % frame)).read().splitlines()
    S = [l.split() for l in L[2:26]]
    assert "".join(a[0][0] for a in S) == "CHHCCOOOCHCHCCHCHCHOHCOO"
    return S, np.array([[float(v) for v in a[1:4]] for a in S])

def h8_closest(X):
    h8 = X[NAMES.index("H8")]
    m = RC > 0
    Xm, R = X[m], 0.5 * (RS + RC[m])
    n = 20000
    k = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * k / n); th = np.pi * (1 + 5 ** 0.5) * k
    U = np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)
    best_t, best_u = 1e9, None
    for u in U:
        w = h8 - Xm; b = w @ u; c = (w * w).sum(1) - R * R; disc = b * b - c
        iv = [(-b[j] - np.sqrt(disc[j]), -b[j] + np.sqrt(disc[j])) for j in range(len(Xm)) if disc[j] > 0]
        t, changed = 0.0, True
        while changed:
            changed = False
            for a, z in iv:
                if a < t < z:
                    t, changed = z, True
        if t < best_t:
            best_t, best_u = t, u
    p = h8 + best_t * best_u
    d = np.linalg.norm(X - p, axis=1)
    assert all(d[j] >= 0.5 * (RS + RC[j]) - 1e-9 for j in range(24) if RC[j] > 0)
    return p, best_t

def write(name, S, p, q, title):
    d = os.path.join(OUT, name); os.makedirs(d)
    with open(os.path.join(d, "combined.xyz"), "w") as f:
        f.write("25\n%s\n" % title)
        for a in S: f.write("%-2s %15.8f%15.8f%15.8f\n" % (a[0], float(a[1]), float(a[2]), float(a[3])))
        f.write("%-2s %15.8f%15.8f%15.8f\n" % ("N", p[0], p[1], p[2]))
    with open(os.path.join(d, "site.ORCAFF.prms"), "w") as f:
        f.write("$fftype\nAMBER\n$atoms\n25 1 4\n")
        for l in ROWS: f.write(l + "\n")
        f.write("%6d   %-2s%13.6f%13.6f%13.6f%13.6f%13.6f\n" % (25, "N", q, -0.17, 3.648, -0.085, 3.648))
        f.write("$bonds\n0 2 2\n$angles\n0 3 2\n$dihedrals\n0 4 3\n")
    open(os.path.join(d, "job.inp"), "w").write(INP)
    open(os.path.join(d, "job.pbs"), "w").write(PBS.format(name=name))

S41, X41 = substrate("41786")
S55, X55 = substrate("55446")
p41, t41 = h8_closest(X41)
p55, t55 = h8_closest(X55)
write("p1_41786_qp05", S41, np.array([40.832, 18.750, 42.551]), 0.5, "frame 41786 R + sphere q=+0.5 at WP4 arm b start")
write("n1_41786_qm1", S41, p41, -1.0, "frame 41786 R + sphere q=-1 at closest LJ-allowed approach to H8 (%.3f A)" % t41)
write("n2_41786_qm05", S41, p41, -0.5, "frame 41786 R + sphere q=-0.5 at closest LJ-allowed approach to H8 (%.3f A)" % t41)
write("n3_55446_qm1", S55, p55, -1.0, "frame 55446 R + sphere q=-1 at closest LJ-allowed approach to H8 (%.3f A)" % t55)
with open(os.path.join(OUT, "placements.tsv"), "w") as f:
    f.write("run\tframe\tq\tx_A\ty_A\tz_A\tto_H8_A\tto_O4_A\n")
    for name, fr, q, p, X in (("p1_41786_qp05", "41786", 0.5, np.array([40.832, 18.750, 42.551]), X41),
                              ("n1_41786_qm1", "41786", -1.0, p41, X41), ("n2_41786_qm05", "41786", -0.5, p41, X41),
                              ("n3_55446_qm1", "55446", -1.0, p55, X55)):
        f.write("%s\t%s\t%+.1f\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\n" % (name, fr, q, p[0], p[1], p[2],
                np.linalg.norm(p - X[NAMES.index("H8")]), np.linalg.norm(p - X[NAMES.index("O4")])))
print("4 runs and placements.tsv written to", OUT)
