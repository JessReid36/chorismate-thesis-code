#!/usr/bin/env python3
"""cs6_d36_amend1_prepare.py - DIAGNOSTIC 36 AMENDMENT 1: arm B again, with resources
that do not starve it onto disk.

Copies the committed d36/B_{R,TS} inputs and changes exactly one line, %maxcore. Writes
one PBS script per leg so that B_TS starts whatever B_R does. Nothing about the method,
the auxiliary basis, the PNO thresholds or the geometry is touched.

Usage:  cs6_d36_amend1_prepare.py <results-repo> <out-dir> <hpc-dir>
"""
import hashlib, os, shutil, sys

S6 = "phase2.2/cs_stage6"
# two resource lines change, and only these two: the memory per core and the number
# of cores ORCA itself uses. PBS granting 16 cores is no use if %pal still says 8.
SUBS = (("%maxcore 5000", "%maxcore 12000"),
        ("%pal nprocs 8 end", "%pal nprocs 16 end"))
METHOD = "! DLPNO-CCSD(T) def2-TZVPD AutoAux TightSCF TightPNO"
NCPUS, MEM, WALL = 16, "240gb", "168:00:00"
LEGS = ("R", "TS")


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
        stop("usage: cs6_d36_amend1_prepare.py <results-repo> <out-dir> <hpc-dir>")
    res, out, hpc = sys.argv[1], sys.argv[2], sys.argv[3]
    s6 = os.path.join(res, S6)

    for leg in LEGS:
        src = os.path.join(s6, "d36", "B_%s" % leg, "job.inp")
        if not os.path.isfile(src):
            stop("%s is missing; amendment 1 copies attempt 1's own input." % src)
        t = open(src).read()
        if t.splitlines()[0] != METHOD:
            stop("%s: the method line is %r, not the pre-registered one.\n"
                 "         Amendment 1 changes resources only and will not proceed on a\n"
                 "         different calculation." % (src, t.splitlines()[0]))
        for old, _ in SUBS:
            if t.count(old + "\n") != 1:
                stop("%s: expected exactly one %r line, found %d."
                     % (src, old, t.count(old + "\n")))
        g = os.path.join(s6, "d33a", "sp_42436_%s" % leg, "substrate.xyz")
        if not os.path.isfile(g):
            stop("%s is missing; it is the geometry arm B runs on." % g)
        gsrc = os.path.join(s6, "d36", "B_%s" % leg, "substrate.xyz")
        if os.path.isfile(gsrc) and sha(gsrc) != sha(g):
            stop("attempt 1's %s geometry is not frame 42436's committed substrate." % leg)

    os.makedirs(out, exist_ok=True)
    sums = []
    for leg in LEGS:
        src = os.path.join(s6, "d36", "B_%s" % leg, "job.inp")
        d = os.path.join(out, "B_%s" % leg)
        os.makedirs(d, exist_ok=True)
        t = open(src).read()
        for old, new in SUBS:
            t = t.replace(old + "\n", new + "\n")
        open(os.path.join(d, "job.inp"), "w").write(t)
        g = os.path.join(s6, "d33a", "sp_42436_%s" % leg, "substrate.xyz")
        shutil.copyfile(g, os.path.join(d, "substrate.xyz"))
        if sha(g) != sha(os.path.join(d, "substrate.xyz")):
            stop("%s did not copy byte for byte into %s" % (g, d))
        a, b = open(src).read().splitlines(), t.splitlines()
        diff = [(x, y) for x, y in zip(a, b) if x != y]
        if len(a) != len(b) or sorted(diff) != sorted(SUBS):
            stop("B_%s: the new input differs from the committed one in %d lines, "
                 "expected exactly the two resource lines." % (leg, len(diff)))
        print("   B_%-3s %s" % (leg, t.splitlines()[0]))
        for old, new in SUBS:
            print("         %-20s -> %-21s" % (old, new))
        for fn in sorted(os.listdir(d)):
            sums.append((sha(os.path.join(d, fn)), "B_%s/%s" % (leg, fn)))

    for leg in LEGS:
        p = os.path.join(out, "run_B_%s.pbs" % leg)
        with open(p, "w") as f:
            f.write("#!/bin/bash\n#PBS -N cm_d36a1_%s\n" % leg)
            f.write("#PBS -l select=1:ncpus=%d:mem=%s\n#PBS -l walltime=%s\n" % (NCPUS, MEM, WALL))
            f.write("#PBS -m ae\n#PBS -M 18660916@sun.ac.za\n#PBS -j oe\n")
            f.write("#PBS -o %s/run_B_%s.pbs.out\n" % (hpc, leg))
            f.write("# DIAGNOSTIC 36 AMENDMENT 1: arm B, %s leg, one job of its own.\n" % leg)
            f.write("# Attempt 1 asked %%maxcore 5000 on 8 cores and ORCA put the PNO integrals\n")
            f.write("#   ON DISK, writing 36 GB onto NFS and streaming it through every CCSD\n")
            f.write("#   iteration. %%maxcore %d on %d cores is %d GB inside a %s request.\n"
                    % (12000, NCPUS, 12000 * NCPUS // 1000, MEM))
            f.write("# The method line, AutoAux, TightPNO, the basis and the geometry are\n")
            f.write("#   unchanged. Only resources change, and resources cannot move the answer.\n")
            f.write("# walltime %s routes to week automatically. No -q.\n" % WALL)
            f.write('export PATH="/apps/openmpi/4.1.1/bin:$PATH"\n')
            f.write('export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:'
                    '/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"\n')
            f.write("export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1\n")
            f.write("ORCA=/home/apps2/ORCA/6.0.1\n")
            f.write('W=%s/B_%s\ncd "$W" || exit 1\n' % (hpc, leg))
            f.write('if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running in $W"; exit 1; fi\n')
            f.write("trap 'rmdir \"$W/.running\" 2>/dev/null' EXIT\n")
            f.write('echo "host=$(hostname) start=$(date)"\n')
            f.write('echo "cores=$(nproc) ; memory:"; free -g | head -2\n')
            f.write('if [ -s job.out ] && grep -q "ORCA TERMINATED NORMALLY" job.out; then\n')
            f.write('  echo "B_%s done before"; exit 0; fi\n' % leg)
            f.write('"$ORCA/orca" job.inp > job.out 2>&1 </dev/null\n')
            f.write('echo "orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' job.out) $(date)"\n')
            f.write('echo "PNO integral storage: $(grep -m1 \'PNO-Integral Storage\' job.out)"\n')
            f.write('echo "ccsdt energy: $(grep -m1 \'E(CCSD(T))\' job.out)"\n')
            f.write('grep -m1 "TOTAL RUN TIME" job.out\n')
        os.chmod(p, 0o755)
        sums.append((sha(p), "run_B_%s.pbs" % leg))

    with open(os.path.join(out, "INPUT_CHECKSUMS.sha256"), "w") as f:
        for h, n in sums:
            f.write("%s  %s\n" % (h, n))
    print("   2 run directories, 2 job scripts, %d files checksummed" % len(sums))
    print("   %d cores, %s, walltime %s per leg" % (NCPUS, MEM, WALL))


if __name__ == "__main__":
    main()
