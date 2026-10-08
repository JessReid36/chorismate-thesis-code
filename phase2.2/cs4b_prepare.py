#!/usr/bin/env python3
"""
cs4b_prepare.py - Stage 4b of the charged-sphere development (pre-registered; criteria in STAGE4B_CRITERIA.txt).
Built from committed files only: Stage 3 p1's converged geometry (sphere 2.674 A from O2), its force-field file,
Stage 4's single-point input (B3LYP-D3BJ, def2/J, RIJCOSX, TightSCF, Hirshfeld printed).

Arm A - safe-region map for the bare sphere (ORCA QM/MM, sphere = MM atom with q + Amber N LJ, as before).
  The sphere is moved out along the line O2 -> sphere to 2.90, 3.20 and 3.50 A from O2 (nearest atom H8 2.77, 3.02,
  3.28 A); at each, q = -1, -0.5, 0, +0.5, +0.75, +1 in def2-SVP and def2-SVPD. At 2.674 A only q = +0.75 is new
  (Stage 4's committed single points supply the rest).
Arm B - pseudopotential sphere (plain ORCA, no QM/MM: the sphere is a coreless ECP centre "Ne>" with charge q and the
  Ne-type all-electron pseudopotential of Marefat Khah et al. 2020, SI Table S2, at p1's sphere position).
  q = -1, -0.5, 0, +0.5, +0.75, +1 in def2-SVP and def2-SVPD, the ECP assigned in the %basis block (NewECP Ne);
  controls: the same centre with no ECP at q = +1 and -1 (def2-SVP), which must reproduce Stage 4's QM/MM energies
  minus their printed LJ term; and one run with the ECP given inline on the centre instead (q = +1, def2-SVP), which
  must equal the %basis run.
ECP translation (TURBOMOLE columns coefficient / power / exponent -> ORCA exponent / coefficient / power; power 2 = r^0
in both): s 2.0475 / 54.5100, p 0.44815 / 1.46500, d 0.49205 / -0.8380, local f 1.0 / 0.0; N_core 10.
USAGE  python3 cs4b_prepare.py <results repo> <output folder>
"""
import os, sys
import numpy as np
REPO, OUT = sys.argv[1], sys.argv[2]
SD = "/home/18660916/system_development/phase2.2/cs_stage4b"
P1 = os.path.join(REPO, "phase2.2/cs_stage3/p1_41786_qp05")
S4 = os.path.join(REPO, "phase2.2/cs_stage4")
L = open(os.path.join(P1, "job.xyz")).read().splitlines()
assert int(L[0]) == 25
ATOMS = [l.split() for l in L[2:27]]
X = np.array([[float(v) for v in a[1:4]] for a in ATOMS])
O2 = X[6]; S0 = X[24]
U = (S0 - O2) / np.linalg.norm(S0 - O2)
PR = open(os.path.join(P1, "site.ORCAFF.prms")).read().splitlines()
SPIN = open(os.path.join(S4, "sp_svp_q0/job.inp")).read()
assert SPIN.startswith("! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF QMMM")
ENV = """export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1"""
HDR = """#!/bin/bash
#PBS -N cm34b_%s
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -l walltime=168:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %s/%s.pbs.out
"""
QS = ((-1.0, "qm1"), (-0.5, "qm05"), (0.0, "q0"), (0.5, "qp05"), (0.75, "qp075"), (1.0, "qp1"))
BASES = (("def2-SVP", "svp"), ("def2-SVPD", "svpd"))
ECP_BLOCK = """  N_core 10
  lmax f
  s 1
    1  2.0475   54.5100  2
  p 1
    1  0.44815   1.46500  2
  d 1
    1  0.49205  -0.83800  2
  f 1
    1  1.0000    0.0000  2
end"""
os.makedirs(OUT)
batches = {"A_svp": [], "A_svpd": [], "B": []}

def prms(q):
    P = list(PR); i = P.index("$atoms")
    P[i + 2 + 24] = "%6d   %-2s%13.6f%13.6f%13.6f%13.6f%13.6f" % (25, "N", q, -0.17, 3.648, -0.085, 3.648)
    return "\n".join(P) + "\n"

# ---------------- Arm A
for basis, bt in BASES:
    for d, dt in ((2.674, "d2674"), (2.90, "d290"), (3.20, "d320"), (3.50, "d350")):
        p = S0 if dt == "d2674" else O2 + d * U
        for q, qt in QS:
            if dt == "d2674" and qt != "qp075":
                continue
            name = "A_%s_%s_%s" % (bt, dt, qt)
            dd = os.path.join(OUT, name); os.makedirs(dd)
            with open(os.path.join(dd, "combined.xyz"), "w") as f:
                f.write("25\np1 geometry, sphere %.3f A from O2, q=%+.2f, %s\n" % (d, q, basis))
                for a in ATOMS[:24]: f.write("%-2s %15.8f%15.8f%15.8f\n" % (a[0], float(a[1]), float(a[2]), float(a[3])))
                f.write("%-2s %15.8f%15.8f%15.8f\n" % ("N", p[0], p[1], p[2]))
            open(os.path.join(dd, "site.ORCAFF.prms"), "w").write(prms(q))
            open(os.path.join(dd, "job.inp"), "w").write(SPIN.replace("def2-SVP def2/J", basis + " def2/J", 1))
            batches["A_" + bt].append(name)

# ---------------- Arm B (plain ORCA, sphere = coreless ECP centre)
HEAD = """! B3LYP D3BJ %s def2/J RIJCOSX TightSCF PrintBasis
%%maxcore 3000
%%pal nprocs 8 end
%%scf MaxIter 300 end
%%output
  Print[P_Hirshfeld] 1
end
"""
def sub_lines():
    return "".join("  %-2s %15.8f%15.8f%15.8f\n" % (a[0], float(a[1]), float(a[2]), float(a[3])) for a in ATOMS[:24])
def write_B(name, basis, q, ecp):   # ecp: "basis" | "inline" | None
    dd = os.path.join(OUT, name); os.makedirs(dd)
    s = HEAD % basis
    if ecp == "basis":
        s += "%basis\n  NewECP Ne\n" + "\n".join("  " + l for l in ECP_BLOCK.splitlines()[:-1]) + "\n  end\nend\n"
    s += "* xyz -2 1\n" + sub_lines()
    if ecp == "inline":
        s += "  Ne>  %.6f %15.8f%15.8f%15.8f  NewECP\n" % (q, S0[0], S0[1], S0[2]) + ECP_BLOCK + "\n"
    else:
        s += "  Ne>  %.6f %15.8f%15.8f%15.8f\n" % (q, S0[0], S0[1], S0[2])
    s += "*\n"
    open(os.path.join(dd, "job.inp"), "w").write(s)
    batches["B"].append(name)
for basis, bt in BASES:
    for q, qt in QS:
        write_B("B_ecp_%s_%s" % (bt, qt), basis, q, "basis")
write_B("B_noecp_svp_qp1", "def2-SVP", 1.0, None)
write_B("B_noecp_svp_qm1", "def2-SVP", -1.0, None)
write_B("B_inline_svp_qp1", "def2-SVP", 1.0, "inline")

for b, runs in batches.items():
    with open(os.path.join(OUT, "batch_%s.pbs" % b), "w") as f:
        f.write(HDR % (b, SD, "batch_" + b) + ENV + "\n" + 'echo "host=$(hostname) start=$(date)"\n')
        f.write("for d in %s; do\n  cd %s/$d\n  /home/apps2/ORCA/6.0.1/orca job.inp > job.out 2>&1 </dev/null\n" % (" ".join(runs), SD))
        f.write('  echo "$d orca_exit=$? $(grep -c \'ORCA TERMINATED NORMALLY\' job.out) $(date)"\n  rm -f job.gbw job*.tmp 2>/dev/null\ndone\necho "end=$(date)"\n')
with open(os.path.join(OUT, "sphere_positions.tsv"), "w") as f:
    f.write("label\td_O2_A\tx\ty\tz\tnearest_atom_A\n")
    for d, dt in ((2.674, "d2674"), (2.90, "d290"), (3.20, "d320"), (3.50, "d350")):
        p = S0 if dt == "d2674" else O2 + d * U
        f.write("%s\t%.4f\t%.8f\t%.8f\t%.8f\t%.4f\n" % (dt, np.linalg.norm(p - O2), p[0], p[1], p[2], np.linalg.norm(X[:24] - p, axis=1).min()))
print({k: len(v) for k, v in batches.items()}, "runs written to", OUT)
