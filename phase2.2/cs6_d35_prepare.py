#!/usr/bin/env python3
"""cs6_d35_prepare.py - DIAGNOSTIC 35: Stage 4e's five levels on frame 42436.

Copies each of Stage 4e's bare inputs byte for byte and places frame 42436's committed
substrate beside it. Nothing inside an input is edited, so the only variable is the
geometry. The script refuses if a source input is missing, if the two legs of a level
are not the same file, or if a copy does not match its source by sha256.

Usage:  cs6_d35_prepare.py <results-repo> <out-dir> <hpc-dir>
"""
import hashlib, os, shutil, sys

S6 = "phase2.2/cs_stage6"
S4E = "phase2.2/cs_stage4e"
LEVELS = ("L1", "L2", "L3", "L4", "L5")
LEGS = ("R", "TS")
GEOM = {"R": "d33a/sp_42436_R/substrate.xyz", "TS": "d33a/sp_42436_TS/substrate.xyz"}
WALL = "24:00:00"


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
        stop("usage: cs6_d35_prepare.py <results-repo> <out-dir> <hpc-dir>")
    res, out, hpc = sys.argv[1], sys.argv[2], sys.argv[3]
    s6, s4e = os.path.join(res, S6), os.path.join(res, S4E)

    for leg, g in GEOM.items():
        if not os.path.isfile(os.path.join(s6, g)):
            stop("%s is missing; it is the geometry D35 runs on." % g)

    srcs = {}
    for lv in LEVELS:
        paths = [os.path.join(s4e, "%s_bare_%s" % (lv, leg), "job.inp") for leg in LEGS]
        for p in paths:
            if not os.path.isfile(p):
                stop("%s is missing; D35 copies Stage 4e's own inputs." % p)
        h = {sha(p) for p in paths}
        if len(h) != 1:
            stop("%s's reactant and transition-state inputs differ.\n"
                 "         D35 assumes the only difference between the two legs is the\n"
                 "         geometry file. Refusing to proceed on an assumption that is false." % lv)
        srcs[lv] = paths[0]
        first = open(paths[0]).read().splitlines()[0]
        if not first.startswith("! "):
            stop("%s: the first line is not an ORCA keyword line." % paths[0])
        print("   %-3s %s" % (lv, first))

    os.makedirs(out, exist_ok=True)
    made, sums = [], []
    for lv in LEVELS:
        for leg in LEGS:
            d = os.path.join(out, "%s_%s" % (lv, leg))
            os.makedirs(d, exist_ok=True)
            for src, dst in ((srcs[lv], os.path.join(d, "job.inp")),
                             (os.path.join(s6, GEOM[leg]), os.path.join(d, "substrate.xyz"))):
                shutil.copyfile(src, dst)
                if sha(src) != sha(dst):
                    stop("%s did not copy byte for byte into %s" % (src, d))
            made.append(os.path.basename(d))
            for fn in sorted(os.listdir(d)):
                sums.append((sha(os.path.join(d, fn)), "%s/%s" % (os.path.basename(d), fn)))

    if len(made) != len(LEVELS) * len(LEGS):
        stop("built %d run directories, expected %d." % (len(made), len(LEVELS) * len(LEGS)))

    pbs = os.path.join(out, "run.pbs")
    with open(pbs, "w") as f:
        f.write("#!/bin/bash\n#PBS -N cm_d35\n")
        f.write("#PBS -l select=1:ncpus=8:mem=24gb\n#PBS -l walltime=%s\n" % WALL)
        f.write("#PBS -m ae\n#PBS -M 18660916@sun.ac.za\n#PBS -j oe\n")
        f.write("#PBS -o %s/run.pbs.out\n" % hpc)
        f.write("# Stage 6 DIAGNOSTIC 35: Stage 4e's five levels on frame 42436's stationary points.\n")
        f.write("# Each job.inp is Stage 4e's own file, byte for byte; only substrate.xyz differs.\n")
        f.write("# L1 runs first and is the control: its level already exists at this frame and it\n")
        f.write("#   must reproduce 18.461014 kcal/mol (G35.3). Nothing else is readable without it.\n")
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
        f.write("run() {\n")
        f.write('  d="$W/$1"\n')
        f.write('  if [ -s "$d/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$d/job.out"; then echo "$1 done before"\n')
        f.write('  else ( cd "$d" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); '
                'echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi\n')
        f.write("}\n")
        for d in made:
            f.write("run %s\n" % d)
        f.write('echo "end=$(date)"\n')
        f.write('echo "finished: $(grep -l \'ORCA TERMINATED NORMALLY\' */job.out 2>/dev/null | wc -l) of %d"\n' % len(made))
    os.chmod(pbs, 0o755)
    sums.append((sha(pbs), "run.pbs"))

    with open(os.path.join(out, "INPUT_CHECKSUMS.sha256"), "w") as f:
        for h, n in sums:
            f.write("%s  %s\n" % (h, n))
    print("   %d run directories, %d files checksummed, walltime %s" % (len(made), len(sums), WALL))


if __name__ == "__main__":
    main()
