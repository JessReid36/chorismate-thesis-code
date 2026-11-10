#!/usr/bin/env python3
"""cs6_d34_prepare.py - DIAGNOSTIC 34 stage B: lift the committed frozen single points
to CAM-B3LYP/def2-SVPD, in two arms.

Arm 1 changes exactly two lines of each committed input, the method line and the
iteration cap, and copies everything else byte for byte.  Arm 2 is arm 1 with the
Amendment 5 stability setting added, and differs from arm 1 in the %scf block alone.

The run set is taken from the committed outputs tarball, not from a list written here,
so it cannot drift from what the frozen half actually measured.

Usage:  cs6_d34_prepare.py <results-repo> <out-dir> <hpc-dir>
"""
import hashlib, os, re, shutil, sys, tarfile

S6 = "phase2.2/cs_stage6"

OLD_QMMM = "! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF QMMM EnGrad"
NEW_QMMM = "! CAM-B3LYP D3BJ def2-SVPD def2/J RIJCOSX TightSCF QMMM EnGrad"
OLD_BARE = "! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF EnGrad"
NEW_BARE = "! CAM-B3LYP D3BJ def2-SVPD def2/J RIJCOSX TightSCF EnGrad"
METHOD = ((OLD_QMMM, NEW_QMMM), (OLD_BARE, NEW_BARE))

OLD_SCF = "%scf MaxIter 250 end"
A1_SCF = "%scf MaxIter 300 end"
A2_SCF = "%scf\n  MaxIter 300\n  STABPerform true\nend"

WALL = "24:00:00"


def stop(msg):
    sys.stderr.write("REFUSED  %s\n" % msg)
    raise SystemExit(2)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def sphere_charges(prms):
    """The MM charges of the spheres: everything after the 24 QM atoms."""
    L = open(prms).read().splitlines()
    try:
        i = L.index("$atoms")
    except ValueError:
        stop("%s: no $atoms block." % prms)
    n = int(L[i + 1].split()[0])
    q = [float(L[i + 2 + k].split()[2]) for k in range(n)]
    if n < 24:
        stop("%s: %d atoms, fewer than the 24 QM atoms." % (prms, n))
    return q[24:]


def runs_from_outputs(tarpath):
    """(arrangement, leg) for every committed frozen output, read off the tarball."""
    got = set()
    with tarfile.open(tarpath) as tf:
        for n in tf.getnames():
            m = re.match(r"\./([^/]+)/(R|TS)_sp/job\.out$", n)
            if m:
                got.add((m.group(1), m.group(2)))
    if not got:
        stop("%s: no ./<arrangement>/<leg>_sp/job.out entries." % tarpath)
    return sorted(got)


def lift(src_inp, arm):
    t = open(src_inp).read()
    hits = [(o, n) for o, n in METHOD if t.count(o + "\n") == 1]
    if len(hits) != 1:
        stop("%s: expected exactly one of the two committed method lines, matched %d.\n"
             "         Refusing to guess which line to change." % (src_inp, len(hits)))
    old, new = hits[0]
    if t.count(OLD_SCF + "\n") != 1:
        stop("%s: expected exactly one %r line, found %d."
             % (src_inp, OLD_SCF, t.count(OLD_SCF + "\n")))
    t = t.replace(old + "\n", new + "\n")
    t = t.replace(OLD_SCF + "\n", (A1_SCF if arm == 1 else A2_SCF) + "\n")
    # nothing else may change: the only differences from the source are those two
    return t, old, new


def main():
    if len(sys.argv) != 4:
        stop("usage: cs6_d34_prepare.py <results-repo> <out-dir> <hpc-dir>")
    res, out, hpc = sys.argv[1], sys.argv[2], sys.argv[3]
    s6 = os.path.join(res, S6)
    tarpath = os.path.join(s6, "frozen", "sp_outputs.tar.xz")
    if not os.path.isfile(tarpath):
        stop("%s is missing; the run set is read from it." % tarpath)

    runs = runs_from_outputs(tarpath)
    arrs = sorted({a for a, _ in runs})
    print("   %d committed frozen outputs over %d arrangements" % (len(runs), len(arrs)))

    # the charge budget of each arrangement, from its own committed force field
    budget, nsph = {}, {}
    for a in arrs:
        prms = os.path.join(s6, "relax", a, "R_sp", "site.ORCAFF.prms")
        if os.path.isfile(prms):
            q = sphere_charges(prms)
            budget[a], nsph[a] = sum(abs(x) for x in q), len(q)
        else:
            budget[a], nsph[a] = 0.0, 0          # the charge-free control

    ctl = [a for a in arrs if nsph[a] == 0]
    if len(ctl) != 1:
        stop("expected exactly one charge-free arrangement as the control, found %d: %s"
             % (len(ctl), ", ".join(ctl) or "none"))
    ctl = ctl[0]
    order = [ctl] + sorted((a for a in arrs if a != ctl), key=lambda a: (budget[a], a))
    print("   control %s; then ascending Sigma|q|: %s"
          % (ctl, ", ".join("%s %.4f" % (a, budget[a]) for a in order[1:])))

    os.makedirs(out, exist_ok=True)
    made, sums = [], []
    for arm in (1, 2):
        for a in order:
            for leg in ("R", "TS"):
                if (a, leg) not in runs:
                    continue
                src = os.path.join(s6, "relax", a, "%s_sp" % leg)
                dst = os.path.join(out, "a%d_%s_%s" % (arm, a, leg))
                os.makedirs(dst, exist_ok=True)
                txt, old, new = lift(os.path.join(src, "job.inp"), arm)
                with open(os.path.join(dst, "job.inp"), "w") as f:
                    f.write(txt)
                for fn in ("combined.xyz", "site.ORCAFF.prms"):
                    p = os.path.join(src, fn)
                    if os.path.isfile(p):
                        shutil.copyfile(p, os.path.join(dst, fn))
                        if sha(p) != sha(os.path.join(dst, fn)):
                            stop("%s did not copy byte for byte into %s" % (p, dst))
                    elif fn == "combined.xyz":
                        stop("%s is missing." % p)
                    elif nsph[a] != 0:
                        stop("%s has %d spheres but no site.ORCAFF.prms." % (a, nsph[a]))
                made.append(os.path.basename(dst))
                for fn in sorted(os.listdir(dst)):
                    sums.append((sha(os.path.join(dst, fn)), "%s/%s" % (os.path.basename(dst), fn)))

    if len(made) != 2 * len(runs):
        stop("built %d run directories, expected %d." % (len(made), 2 * len(runs)))

    pbs = os.path.join(out, "run.pbs")
    with open(pbs, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("#PBS -N cm_d34\n")
        f.write("#PBS -l select=1:ncpus=8:mem=24gb\n")
        f.write("#PBS -l walltime=%s\n" % WALL)
        f.write("#PBS -m ae\n#PBS -M 18660916@sun.ac.za\n#PBS -j oe\n")
        f.write("#PBS -o %s/run.pbs.out\n" % hpc)
        f.write("# Stage 6 DIAGNOSTIC 34 stage B: %d runs at CAM-B3LYP/def2-SVPD.\n" % len(made))
        f.write("#   arm 1 (a1_*) is the lift and carries no stability analysis, matching the\n")
        f.write("#   committed runs, so arm 1 against the committed ddE changes one variable.\n")
        f.write("#   arm 2 (a2_*) adds STABPerform, which no committed QM/MM run has used; it is\n")
        f.write("#   ordered last so a rejection costs nothing already measured.\n")
        f.write("# The control runs first in each arm: it is the gate on everything else.\n")
        f.write("# Within each arm, ascending Sigma|q|, so the designs inside D30's confirmed\n")
        f.write("#   basis window land before the three outside it.\n")
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

    print("   %d run directories, %d input files checksummed, run.pbs walltime %s"
          % (len(made), len(sums), WALL))
    print("   arm 1: %d   arm 2: %d" % (sum(1 for d in made if d.startswith("a1_")),
                                        sum(1 for d in made if d.startswith("a2_"))))


if __name__ == "__main__":
    main()
