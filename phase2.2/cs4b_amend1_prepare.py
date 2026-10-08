#!/usr/bin/env python3
"""
cs4b_amend1_prepare.py - Stage 4b, Amendment 1 (STAGE4B_AMENDMENT1.txt). Writes the eight replacement runs for Arm B
into phase2.2/cs_stage4b/, beside the pre-registered inputs (which are not touched), from the same committed files
cs4b_prepare.py used (Stage 3 p1's converged geometry):
  B2_inline_{svp,svpd}_q0        the Ne-type pseudopotential given inline on the Ne> centre, q = 0
  B2_inline_{svp,svpd}_q{p,m}00001  the same at q = +-0.0001 (fallback if ORCA refuses a zero-charge centre)
  B2_pc_svp_q{p,m}1              no pseudopotential: the sphere as an external point-charge file (%pointcharges),
                                 q = +-1, the control for check 1
Method line, substrate lines and ECP block are formatted exactly as in cs4b_prepare.py's Arm B inputs.
USAGE  python3 cs4b_amend1_prepare.py <results repo> <cs_stage4b folder>
"""
import os, sys
REPO, OUT = sys.argv[1], sys.argv[2]
SD = "/home/18660916/system_development/phase2.2/cs_stage4b"
P1 = os.path.join(REPO, "phase2.2/cs_stage3/p1_41786_qp05")
L = open(os.path.join(P1, "job.xyz")).read().splitlines()
assert int(L[0]) == 25
ATOMS = [l.split() for l in L[2:27]]
S0 = [float(v) for v in ATOMS[24][1:4]]
assert os.path.isdir(OUT) and os.path.exists(os.path.join(OUT, "STAGE4B_CRITERIA.txt")), "cs_stage4b must already exist"
HEAD = """! B3LYP D3BJ %s def2/J RIJCOSX TightSCF PrintBasis
%%maxcore 3000
%%pal nprocs 8 end
%%scf MaxIter 300 end
%%output
  Print[P_Hirshfeld] 1
end
"""
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
SUB = "".join("  %-2s %15.8f%15.8f%15.8f\n" % (a[0], float(a[1]), float(a[2]), float(a[3])) for a in ATOMS[:24])
runs = []
def new(name):
    d = os.path.join(OUT, name); os.makedirs(d); runs.append(name); return d
for basis, bt in (("def2-SVP", "svp"), ("def2-SVPD", "svpd")):
    for q, qt in ((0.0, "q0"), (0.0001, "qp00001"), (-0.0001, "qm00001")):
        d = new("B2_inline_%s_%s" % (bt, qt))
        s = HEAD % basis + "* xyz -2 1\n" + SUB
        s += "  Ne>  %.6f %15.8f%15.8f%15.8f  NewECP\n" % (q, S0[0], S0[1], S0[2]) + ECP_BLOCK + "\n*\n"
        open(os.path.join(d, "job.inp"), "w").write(s)
for q, qt in ((1.0, "qp1"), (-1.0, "qm1")):
    d = new("B2_pc_svp_%s" % qt)
    s = HEAD % "def2-SVP" + '%pointcharges "sphere.pc"\n' + "* xyz -2 1\n" + SUB + "*\n"
    open(os.path.join(d, "job.inp"), "w").write(s)
    open(os.path.join(d, "sphere.pc"), "w").write("1\n%10.6f %15.8f%15.8f%15.8f\n" % (q, S0[0], S0[1], S0[2]))
with open(os.path.join(OUT, "batch_B2.pbs"), "w") as f:
    f.write("""#!/bin/bash
#PBS -N cm34b_B2
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -l walltime=168:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %s/batch_B2.pbs.out
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
echo "host=$(hostname) start=$(date)"
""" % SD)
    f.write("for d in %s; do\n  cd %s/$d\n  /home/apps2/ORCA/6.0.1/orca job.inp > job.out 2>&1 </dev/null\n" % (" ".join(runs), SD))
    f.write('  echo "$d orca_exit=$? $(grep -c \'ORCA TERMINATED NORMALLY\' job.out) $(date)"\n  rm -f job.gbw job*.tmp 2>/dev/null\ndone\necho "end=$(date)"\n')
print(len(runs), "replacement runs written to", OUT)
