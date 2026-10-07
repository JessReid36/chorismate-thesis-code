#!/usr/bin/env python3
"""
cs4_prepare.py - Stage 4 of the charged-sphere development: spill-out check, def2-SVPD against def2-SVP, for
positive spheres at contact distance (negative spheres push density away and are not tested). Built from
committed files only.
  o1_p1_svpd   restart Stage 3's converged p1 (frame 41786, q +0.5, nearest O2 2.674 A with def2-SVP)
               with def2-SVPD; same anchored set-up (C5-C8 and the sphere held). Starting from the SVP minimum
               measures how that minimum moves, without the risk of a different binding mode.
  o2_q1_svpd   the same for WP4 armb_anch (q +1, nearest O2 2.711 A with def2-SVP)
  sp_<basis>_q<..>  ten single points on p1's converged geometry (sphere 2.674 A from O2), q = 0, +-0.5, +-1,
               with def2-SVP and def2-SVPD, Hirshfeld charges printed: first-order term a and polarisation term
               b at true contact distance in both basis sets
Method otherwise as Stage 3 / WP4: B3LYP-D3BJ, def2/J, RIJCOSX, TightSCF, QM/MM, Amber N sphere.
USAGE  python3 cs4_prepare.py <results repo> <output folder>
"""
import os, sys
REPO, OUT = sys.argv[1], sys.argv[2]
SD = "/home/18660916/system_development/phase2.2/cs_stage4"
P1 = os.path.join(REPO, "phase2.2/cs_stage3/p1_41786_qp05")
WB = os.path.join(REPO, "phase2.2/wp4_s13_rerun/armb_anch")
INP_OPT = open(os.path.join(P1, "job.inp")).read()
assert INP_OPT.startswith("! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF QMMM Opt") and open(os.path.join(WB, "job.inp")).read() == INP_OPT
ENV = """export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1"""
HDR = """#!/bin/bash
#PBS -N cm34_%s
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -l walltime=168:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %s
"""
def final_xyz(folder):
    L = open(os.path.join(folder, "job.xyz")).read().splitlines()
    assert int(L[0]) == 25
    return [l.split() for l in L[2:27]]
def prms_with(q, folder):
    L = open(os.path.join(folder, "site.ORCAFF.prms")).read().splitlines()
    i = L.index("$atoms")
    L[i + 2 + 24] = "%6d   %-2s%13.6f%13.6f%13.6f%13.6f%13.6f" % (25, "N", q, -0.17, 3.648, -0.085, 3.648)
    return "\n".join(L) + "\n"
def write_xyz(path, atoms, title):
    with open(path, "w") as f:
        f.write("25\n%s\n" % title)
        for a in atoms: f.write("%-2s %15.8f%15.8f%15.8f\n" % (a[0], float(a[1]), float(a[2]), float(a[3])))
os.makedirs(OUT)
for name, src, q in (("o1_p1_svpd", P1, 0.5), ("o2_q1_svpd", WB, 1.0)):
    d = os.path.join(OUT, name); os.makedirs(d)
    write_xyz(os.path.join(d, "combined.xyz"), final_xyz(src), "restart of the def2-SVP minimum (%s) with def2-SVPD" % os.path.basename(src))
    open(os.path.join(d, "site.ORCAFF.prms"), "w").write(prms_with(q, src))
    open(os.path.join(d, "job.inp"), "w").write(INP_OPT.replace("def2-SVP def2/J", "def2-SVPD def2/J", 1))
    open(os.path.join(d, "job.pbs"), "w").write(HDR % (name, "%s/%s/job.pbs.out" % (SD, name)) + "cd %s/%s\n%s\n" % (SD, name, ENV)
        + 'echo "host=$(hostname) start=$(date)"\n/home/apps2/ORCA/6.0.1/orca job.inp > job.out 2>&1 </dev/null\n'
        + 'echo "orca_exit=$? end=$(date)"\nrm -f job.gbw job*.tmp 2>/dev/null\n')
geo = final_xyz(P1)
SP = """! B3LYP D3BJ %s def2/J RIJCOSX TightSCF QMMM
%%maxcore 3000
%%pal nprocs 8 end
%%scf MaxIter 300 end
%%output
  Print[P_Hirshfeld] 1
end
%%qmmm
  QMAtoms {0:23} end
  ORCAFFFilename "site.ORCAFF.prms"
end
* xyzfile -2 1 combined.xyz
"""
runs = []
for basis, tag in (("def2-SVP", "svp"), ("def2-SVPD", "svpd")):
    for q, qt in ((0.0, "q0"), (0.5, "qp05"), (-0.5, "qm05"), (1.0, "qp1"), (-1.0, "qm1")):
        name = "sp_%s_%s" % (tag, qt); runs.append(name)
        d = os.path.join(OUT, name); os.makedirs(d)
        write_xyz(os.path.join(d, "combined.xyz"), geo, "p1 converged geometry (sphere 2.674 A from O2), q=%+.1f, %s" % (q, basis))
        open(os.path.join(d, "site.ORCAFF.prms"), "w").write(prms_with(q, P1))
        open(os.path.join(d, "job.inp"), "w").write(SP % basis)
with open(os.path.join(OUT, "sp_batch.pbs"), "w") as f:
    f.write(HDR % ("sp_batch", "%s/sp_batch.pbs.out" % SD) + ENV + "\n" + 'echo "host=$(hostname) start=$(date)"\n')
    f.write("for d in %s; do\n  cd %s/$d\n  /home/apps2/ORCA/6.0.1/orca job.inp > job.out 2>&1 </dev/null\n" % (" ".join(runs), SD))
    f.write('  echo "$d orca_exit=$? $(grep -c \'ORCA TERMINATED NORMALLY\' job.out) $(date)"\n  rm -f job.gbw job*.tmp 2>/dev/null\ndone\necho "end=$(date)"\n')
print("2 optimisations, 10 single points and sp_batch.pbs written to", OUT)
