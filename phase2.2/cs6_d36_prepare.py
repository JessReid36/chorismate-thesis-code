#!/usr/bin/env python3
"""cs6_d36_prepare.py - DIAGNOSTIC 36: a correlated reference for the in vacuo barrier.

Six single points on frame 42436's committed bare substrate. Arm C is D28's own input
copied byte for byte and is the control. Arms H and B are written here, because the
record contains no Hartree-Fock or coupled-cluster precedent to copy.

Usage:  cs6_d36_prepare.py <results-repo> <out-dir> <hpc-dir>
"""
import hashlib, os, shutil, sys

S6 = "phase2.2/cs_stage6"
GEOM = {"R": "d33a/sp_42436_R/substrate.xyz", "TS": "d33a/sp_42436_TS/substrate.xyz"}
CTL = "d28/CAM-B3LYP/def2-TZVPD/ctl_%s/job.inp"

ARM_H = """! HF def2-TZVPD def2/J RIJCOSX TightSCF
%maxcore 3000
%pal nprocs 8 end
%scf MaxIter 300 end
* xyzfile -2 1 substrate.xyz
"""

ARM_B = """! DLPNO-CCSD(T) def2-TZVPD AutoAux TightSCF TightPNO
%maxcore 5000
%pal nprocs 8 end
%scf MaxIter 300 end
* xyzfile -2 1 substrate.xyz
"""

# arm, leg -> (how the input is made, a label for the job script)
PLAN = [("C", "R"), ("C", "TS"), ("H", "R"), ("H", "TS"), ("B", "R"), ("B", "TS")]
WALL = "24:00:00"
MEM = "48gb"


def stop(m):
    sys.stderr.write("REFUSED  %s\n" % m)
    raise SystemExit(2)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def main():
    if len(sys.argv) != 4:
        stop("usage: cs6_d36_prepare.py <results-repo> <out-dir> <hpc-dir>")
    res, out, hpc = sys.argv[1], sys.argv[2], sys.argv[3]
    s6 = os.path.join(res, S6)

    for leg, g in GEOM.items():
        if not os.path.isfile(os.path.join(s6, g)):
            stop("%s is missing; it is the geometry D36 runs on." % g)
    for leg in ("R", "TS"):
        p = os.path.join(s6, CTL % leg)
        if not os.path.isfile(p):
            stop("%s is missing; arm C is a byte-for-byte repeat of it." % (CTL % leg))
        # the control's own geometry must already be the one D36 uses, or the repeat
        # is not a repeat
        gsrc = os.path.join(os.path.dirname(p), "substrate.xyz")
        if not os.path.isfile(gsrc):
            stop("%s has no substrate.xyz beside it." % (CTL % leg))
        if sha(gsrc) != sha(os.path.join(s6, GEOM[leg])):
            stop("D28's control geometry at %s is not frame 42436's committed substrate.\n"
                 "         Arm C would not be a repeat of the committed calculation." % leg)
        first = open(p).read().splitlines()[0]
        if "def2-TZVPD" not in first or "CAM-B3LYP" not in first:
            stop("%s does not carry the method line arm C expects: %s" % (CTL % leg, first))

    os.makedirs(out, exist_ok=True)
    made, sums = [], []
    for arm, leg in PLAN:
        d = os.path.join(out, "%s_%s" % (arm, leg))
        os.makedirs(d, exist_ok=True)
        dst = os.path.join(d, "job.inp")
        if arm == "C":
            src = os.path.join(s6, CTL % leg)
            shutil.copyfile(src, dst)
            if sha(src) != sha(dst):
                stop("%s did not copy byte for byte into %s" % (src, d))
        else:
            open(dst, "w").write(ARM_H if arm == "H" else ARM_B)
        g = os.path.join(s6, GEOM[leg])
        shutil.copyfile(g, os.path.join(d, "substrate.xyz"))
        if sha(g) != sha(os.path.join(d, "substrate.xyz")):
            stop("%s did not copy byte for byte into %s" % (g, d))
        t = open(dst).read()
        if "* xyzfile -2 1 substrate.xyz" not in t:
            stop("%s/job.inp does not read substrate.xyz at charge -2, multiplicity 1." % d)
        made.append("%s_%s" % (arm, leg))
        for fn in sorted(os.listdir(d)):
            sums.append((sha(os.path.join(d, fn)), "%s_%s/%s" % (arm, leg, fn)))
        print("   %-5s %s" % ("%s_%s" % (arm, leg), t.splitlines()[0]))

    if len(made) != 6:
        stop("built %d run directories, expected 6." % len(made))

    pbs = os.path.join(out, "run.pbs")
    with open(pbs, "w") as f:
        f.write("#!/bin/bash\n#PBS -N cm_d36\n")
        f.write("#PBS -l select=1:ncpus=8:mem=%s\n#PBS -l walltime=%s\n" % (MEM, WALL))
        f.write("#PBS -m ae\n#PBS -M 18660916@sun.ac.za\n#PBS -j oe\n")
        f.write("#PBS -o %s/run.pbs.out\n" % hpc)
        f.write("# Stage 6 DIAGNOSTIC 36: a correlated reference for the in vacuo barrier.\n")
        f.write("# Arm C is D28's own CAM-B3LYP/def2-TZVPD input byte for byte and is the CONTROL:\n")
        f.write("#   it must return -837.20938385 and -837.16842778 Eh (G36.3). It runs first and\n")
        f.write("#   nothing else is readable without it.\n")
        f.write("# Arm H is the Hartree-Fock reference the coupled-cluster arm is built on.\n")
        f.write("# Arm B is DLPNO-CCSD(T)/def2-TZVPD with AutoAux and TightPNO, and runs LAST, so\n")
        f.write("#   a walltime or memory failure there costs nothing already measured.\n")
        f.write("# No -q: the allocation is automatic and walltime %s routes to day. Restartable.\n" % WALL)
        f.write('export PATH="/apps/openmpi/4.1.1/bin:$PATH"\n')
        f.write('export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:'
                '/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"\n')
        f.write("export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1\n")
        f.write("ORCA=/home/apps2/ORCA/6.0.1\n")
        f.write('W=%s\ncd "$W" || exit 1\n' % hpc)
        f.write('if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running in $W"; exit 1; fi\n')
        f.write("trap 'rmdir \"$W/.running\" 2>/dev/null' EXIT\n")
        f.write('echo "host=$(hostname) start=$(date)"\n')
        f.write('echo "free memory at start:"; free -g | head -2\n')
        f.write("run() {\n")
        f.write('  d="$W/$1"\n')
        f.write('  if [ -s "$d/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$d/job.out"; then echo "$1 done before"\n')
        f.write('  else ( cd "$d" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); '
                'echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi\n')
        f.write("}\n")
        for d in made:
            f.write("run %s\n" % d)
        f.write('echo "end=$(date)"\n')
        f.write('echo "finished: $(grep -l \'ORCA TERMINATED NORMALLY\' */job.out 2>/dev/null | wc -l) of 6"\n')
        f.write('echo "ccsd(t) energies found: '
                '$(grep -l \'E(CCSD(T))\' B_*/job.out 2>/dev/null | wc -l) of 2"\n')
    os.chmod(pbs, 0o755)
    sums.append((sha(pbs), "run.pbs"))

    with open(os.path.join(out, "INPUT_CHECKSUMS.sha256"), "w") as f:
        for h, n in sums:
            f.write("%s  %s\n" % (h, n))
    print("   6 run directories, %d files checksummed, %s, walltime %s" % (len(sums), MEM, WALL))


if __name__ == "__main__":
    main()
