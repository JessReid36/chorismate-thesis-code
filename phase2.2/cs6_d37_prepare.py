#!/usr/bin/env python3
"""cs6_d37_prepare.py - DIAGNOSTIC 37: six candidate functionals on the record's own two tests.

28 single points on frame 42436's committed bare substrate. The control arm is D28's own
CAM-B3LYP inputs at def2-SVPD and def2-TZVPD, copied byte for byte. Every candidate input is
the control input of the same basis with its first line replaced and nothing else touched,
and this script refuses to finish unless that is true of every file it wrote.

Usage:  cs6_d37_prepare.py <results-repo> <out-dir> <hpc-dir>
"""
import hashlib, os, shutil, sys

S6 = "phase2.2/cs_stage6"
GEOM = {"R": "d33a/sp_42436_R/substrate.xyz", "TS": "d33a/sp_42436_TS/substrate.xyz"}
CTL = "d28/CAM-B3LYP/%s/ctl_%s/job.inp"
BASES = [("SV", "def2-SVPD"), ("TZ", "def2-TZVPD")]
CTL_LINE = "! CAM-B3LYP D3BJ %s def2/J RIJCOSX TightSCF"

# tag -> the first line, with %s for the orbital basis. Everything after the functional
# and its dispersion correction is the control's own, in the control's own order; the two
# double hybrids add the correlation-fitting set of the same basis and nothing else.
CAND = [
    ("wB97MV", "! wB97M-V %s def2/J RIJCOSX TightSCF"),
    ("wB97XV", "! wB97X-V %s def2/J RIJCOSX TightSCF"),
    ("PW6B95", "! PW6B95 D3BJ %s def2/J RIJCOSX TightSCF"),
    ("M062X",  "! M062X D3ZERO %s def2/J RIJCOSX TightSCF"),
    ("PWPB95", "! RI-PWPB95 D3BJ %s def2/J %s/C RIJCOSX TightSCF"),
    ("revDSD", "! RI-revDSD-PBEP86-D4/2021 %s def2/J %s/C RIJCOSX TightSCF"),
]
NRUN = 4 + len(CAND) * 4
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


def line1(tmpl, basis):
    return tmpl % ((basis,) * tmpl.count("%s"))


def main():
    if len(sys.argv) != 4:
        stop("usage: cs6_d37_prepare.py <results-repo> <out-dir> <hpc-dir>")
    res, out, hpc = sys.argv[1], sys.argv[2], sys.argv[3]
    s6 = os.path.join(res, S6)
    if os.path.exists(out) and os.listdir(out):
        stop("%s exists and is not empty; this script writes into a new directory only." % out)

    # ---- the sources, before anything is written
    for leg, g in GEOM.items():
        if not os.path.isfile(os.path.join(s6, g)):
            stop("%s is missing; it is the geometry D37 runs on." % g)
    ctl_text = {}
    for tag, basis in BASES:
        both = []
        for leg in ("R", "TS"):
            rel = CTL % (basis, leg)
            p = os.path.join(s6, rel)
            if not os.path.isfile(p):
                stop("%s is missing; the control arm is a byte-for-byte repeat of it." % rel)
            gsrc = os.path.join(os.path.dirname(p), "substrate.xyz")
            if not os.path.isfile(gsrc):
                stop("%s has no substrate.xyz beside it." % rel)
            if sha(gsrc) != sha(os.path.join(s6, GEOM[leg])):
                stop("D28's control geometry for %s at %s is not frame 42436's committed substrate." % (basis, leg))
            t = open(p).read()
            if t.splitlines()[0] != CTL_LINE % basis:
                stop("%s does not carry the method line the control expects: %s" % (rel, t.splitlines()[0]))
            if "* xyzfile -2 1 substrate.xyz" not in t:
                stop("%s does not read substrate.xyz at charge -2, multiplicity 1." % rel)
            both.append(t)
        if both[0] != both[1]:
            stop("D28's two control inputs at %s differ from each other; one template cannot serve both legs." % basis)
        ctl_text[tag] = both[0]

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
        print("   %-13s %s" % (name, open(dst).read().splitlines()[0]))

    # ---- the control, first: both bases, reactant before transition state
    for tag, basis in BASES:
        for leg in ("R", "TS"):
            place("CAM_%s_%s" % (tag, leg), leg, None, os.path.join(s6, CTL % (basis, leg)))
    # ---- every candidate at the small basis, then every candidate at the large one
    for tag, basis in BASES:
        body = ctl_text[tag].split("\n", 1)[1]
        for name, tmpl in CAND:
            for leg in ("R", "TS"):
                place("%s_%s_%s" % (name, tag, leg), leg, line1(tmpl, basis) + "\n" + body)

    # ---- the comparability check, on what was actually written
    if len(made) != NRUN or len(set(made)) != NRUN:
        stop("built %d run directories, expected %d." % (len(made), NRUN))
    first = set()
    for tag, basis in BASES:
        c = open(os.path.join(out, "CAM_%s_R" % tag, "job.inp")).read().split("\n")
        for name, tmpl in CAND:
            for leg in ("R", "TS"):
                t = open(os.path.join(out, "%s_%s_%s" % (name, tag, leg), "job.inp")).read().split("\n")
                if len(t) != len(c) or t[1:] != c[1:]:
                    stop("%s_%s_%s differs from the control below its first line." % (name, tag, leg))
                if t[0] == c[0]:
                    stop("%s_%s_%s carries the control's own method line." % (name, tag, leg))
                want = line1(tmpl, basis)
                if t[0] != want:
                    stop("%s_%s_%s: first line is %r, expected %r" % (name, tag, leg, t[0], want))
                # what the line shares with the control's: basis, Coulomb fitting, COSX, convergence
                for tok in (basis, "def2/J", "RIJCOSX", "TightSCF"):
                    if tok not in t[0].split():
                        stop("%s_%s_%s: first line lacks %s" % (name, tag, leg, tok))
                first.add((name, tag, t[0]))
    if len(first) != len(CAND) * len(BASES):
        stop("expected %d distinct candidate method lines, found %d." % (len(CAND) * len(BASES), len(first)))

    pbs = os.path.join(out, "run.pbs")
    with open(pbs, "w") as f:
        f.write("#!/bin/bash\n#PBS -N cm_d37\n")
        f.write("#PBS -l select=1:ncpus=8:mem=%s\n#PBS -l walltime=%s\n" % (MEM, WALL))
        f.write("#PBS -m ae\n#PBS -M 18660916@sun.ac.za\n#PBS -j oe\n")
        f.write("#PBS -o %s/run.pbs.out\n" % hpc)
        f.write("# Stage 6 DIAGNOSTIC 37: six candidate functionals on the record's own two tests.\n")
        f.write("# The CAM_* runs are D28's own CAM-B3LYP inputs byte for byte and are the CONTROL:\n")
        f.write("#   they must return D28's committed energies (G37.4). They run first and nothing\n")
        f.write("#   else is readable without them.\n")
        f.write("# Then every candidate at def2-SVPD, then every candidate at def2-TZVPD, so each\n")
        f.write("#   keyword this record has never run is first tried in a run of a few minutes.\n")
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
