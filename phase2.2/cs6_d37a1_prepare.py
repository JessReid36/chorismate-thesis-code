#!/usr/bin/env python3
"""cs6_d37a1_prepare.py - DIAGNOSTIC 37, AMENDMENT 1: revDSD-PBEP86-D4, which D37 could not run.

ORCA 6.0.1 rejected D37's keyword RI-revDSD-PBEP86-D4/2021. This writes 10 single points on
the same committed geometries: a control pair (D28's CAM-B3LYP input at def2-SVPD, byte for
byte) and two arms of four runs each.

   arm RI   the functional's own keyword, with RI asked for in an %mp2 block
   arm K    the functional's own keyword and nothing else, as the manual's example writes it

Every input is the control input of the same basis with its first line replaced; arm RI
adds one line, "%mp2 RI true end", and nothing else. The script refuses to finish unless
that is true of every file it wrote.

Usage:  cs6_d37a1_prepare.py <results-repo> <out-dir> <hpc-dir>
"""
import hashlib, os, shutil, sys

S6 = "phase2.2/cs_stage6"
GEOM = {"R": "d33a/sp_42436_R/substrate.xyz", "TS": "d33a/sp_42436_TS/substrate.xyz"}
CTL = "d28/CAM-B3LYP/%s/ctl_%s/job.inp"
BASES = [("SV", "def2-SVPD"), ("TZ", "def2-TZVPD")]
CTL_LINE = "! CAM-B3LYP D3BJ %s def2/J RIJCOSX TightSCF"
LINE = "! revDSD-PBEP86-D4/2021 %s def2/J %s/C RIJCOSX TightSCF"
FAILED = "! RI-revDSD-PBEP86-D4/2021 %s def2/J %s/C RIJCOSX TightSCF"     # what D37 ran
MP2 = "%mp2 RI true end"
GEOLINE = "* xyzfile -2 1 substrate.xyz"
NRUN = 10
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
        stop("usage: cs6_d37a1_prepare.py <results-repo> <out-dir> <hpc-dir>")
    res, out, hpc = sys.argv[1], sys.argv[2], sys.argv[3]
    s6 = os.path.join(res, S6)
    if os.path.exists(out) and os.listdir(out):
        stop("%s exists and is not empty; this script writes into a new directory only." % out)

    ctl_text = {}
    for tag, basis in BASES:
        both = []
        for leg in ("R", "TS"):
            rel = CTL % (basis, leg)
            p = os.path.join(s6, rel)
            if not os.path.isfile(p):
                stop("%s is missing; every input here is built from it." % rel)
            g = os.path.join(os.path.dirname(p), "substrate.xyz")
            if not os.path.isfile(g) or sha(g) != sha(os.path.join(s6, GEOM[leg])):
                stop("D28's control geometry for %s at %s is not frame 42436's committed substrate." % (basis, leg))
            t = open(p).read()
            if t.splitlines()[0] != CTL_LINE % basis:
                stop("%s does not carry the method line the control expects." % rel)
            if t.count(GEOLINE) != 1:
                stop("%s does not read substrate.xyz at charge -2, multiplicity 1, exactly once." % rel)
            both.append(t)
        if both[0] != both[1]:
            stop("D28's two control inputs at %s differ from each other." % basis)
        ctl_text[tag] = both[0]
        # what D37 ran and ORCA rejected: this amendment's line is that line without the prefix
        for leg in ("R", "TS"):
            f = os.path.join(s6, "d37", "revDSD_%s_%s" % (tag, leg), "job.inp")
            if not os.path.isfile(f):
                stop("d37/revDSD_%s_%s/job.inp is missing; the amendment corrects that input." % (tag, leg))
            l0 = open(f).read().splitlines()[0]
            if l0 != FAILED % (basis, basis):
                stop("d37/revDSD_%s_%s does not carry the line ORCA rejected." % (tag, leg))
            if l0.replace("! RI-", "! ", 1) != LINE % (basis, basis):
                stop("the amended line is not D37's line with the RI- prefix removed.")

    os.makedirs(out, exist_ok=True)
    made, sums = [], []

    def place(name, leg, text, copy_from=None):
        d = os.path.join(out, name)
        os.makedirs(d)
        dst = os.path.join(d, "job.inp")
        if copy_from:
            shutil.copyfile(copy_from, dst)
            if sha(copy_from) != sha(dst):
                stop("%s did not copy byte for byte into %s" % (copy_from, name))
        else:
            with open(dst, "w") as f:
                f.write(text)
        g = os.path.join(s6, GEOM[leg])
        shutil.copyfile(g, os.path.join(d, "substrate.xyz"))
        if sha(g) != sha(os.path.join(d, "substrate.xyz")):
            stop("%s did not copy byte for byte into %s" % (g, name))
        made.append(name)
        for fn in sorted(os.listdir(d)):
            sums.append((sha(os.path.join(d, fn)), "%s/%s" % (name, fn)))
        print("   %-16s %s%s" % (name, open(dst).read().splitlines()[0], "   + " + MP2 if MP2 in open(dst).read() else ""))

    # the control, first
    for leg in ("R", "TS"):
        place("CAM_SV_%s" % leg, leg, None, os.path.join(s6, CTL % ("def2-SVPD", leg)))
    # both arms at the small basis, then both at the large one; arm RI before arm K
    for tag, basis in BASES:
        body = ctl_text[tag].split("\n", 1)[1]
        body_ri = body.replace(GEOLINE, MP2 + "\n" + GEOLINE, 1)
        for arm, b in (("RI", body_ri), ("K", body)):
            for leg in ("R", "TS"):
                place("revDSD_%s_%s_%s" % (arm, tag, leg), leg, (LINE % (basis, basis)) + "\n" + b)

    # ---- the comparability check, on what was actually written
    if len(made) != NRUN or len(set(made)) != NRUN:
        stop("built %d run directories, expected %d." % (len(made), NRUN))
    for tag, basis in BASES:
        c = ctl_text[tag].split("\n")
        for leg in ("R", "TS"):
            k = open(os.path.join(out, "revDSD_K_%s_%s" % (tag, leg), "job.inp")).read().split("\n")
            r = open(os.path.join(out, "revDSD_RI_%s_%s" % (tag, leg), "job.inp")).read().split("\n")
            if k[0] != LINE % (basis, basis) or r[0] != k[0]:
                stop("revDSD_*_%s_%s does not carry the amended method line." % (tag, leg))
            if k[1:] != c[1:]:
                stop("revDSD_K_%s_%s differs from the control below its first line." % (tag, leg))
            extra = [l for l in r[1:] if l not in c[1:]]
            if extra != [MP2] or [l for l in r[1:] if l != MP2] != c[1:]:
                stop("revDSD_RI_%s_%s differs from the control by more than the one %%mp2 line." % (tag, leg))
            if r.index(MP2) != r.index(GEOLINE) - 1:
                stop("revDSD_RI_%s_%s: the %%mp2 line is not directly above the geometry line." % (tag, leg))

    pbs = os.path.join(out, "run.pbs")
    with open(pbs, "w") as f:
        f.write("#!/bin/bash\n#PBS -N cm_d37a1\n")
        f.write("#PBS -l select=1:ncpus=8:mem=%s\n#PBS -l walltime=%s\n" % (MEM, WALL))
        f.write("#PBS -m ae\n#PBS -M 18660916@sun.ac.za\n#PBS -j oe\n")
        f.write("#PBS -o %s/run.pbs.out\n" % hpc)
        f.write("# Stage 6 DIAGNOSTIC 37, AMENDMENT 1: revDSD-PBEP86-D4, which D37 could not run.\n")
        f.write("# The CAM_* runs are D28's own CAM-B3LYP input byte for byte and are the CONTROL.\n")
        f.write("# Arm RI asks for RI in an %mp2 block; arm K carries the keyword alone.\n")
        f.write("# A run that fails does not stop the ones after it.\n")
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
        f.write('echo "finished: $(grep -l \'ORCA TERMINATED NORMALLY\' */job.out 2>/dev/null | wc -l) of %d"\n' % NRUN)
    os.chmod(pbs, 0o755)
    sums.append((sha(pbs), "run.pbs"))

    with open(os.path.join(out, "INPUT_CHECKSUMS.sha256"), "w") as f:
        for h, n in sums:
            f.write("%s  %s\n" % (h, n))
    print("   %d run directories, %d files checksummed, %s, walltime %s" % (NRUN, len(sums), MEM, WALL))


if __name__ == "__main__":
    main()
