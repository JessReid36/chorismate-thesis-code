#!/usr/bin/env python3
"""cs6_d40_prepare.py - DIAGNOSTIC 40: detachment energies at the density-functional level.

21 single points on frame 42436. Twenty are monoanions (charge -1, multiplicity 2); each is a
COMMITTED dianion input with its charge-and-multiplicity line changed and nothing else touched,
beside byte-for-byte copies of that run's own geometry, charge and force-field files. One is a
control: D28's CAM-B3LYP dianion input, byte for byte. This script refuses to finish unless all
of that is true of every file it wrote.

Part A  the bare substrate at the reactant, def2-SVPD, under eight functionals
Part B  CAM-B3LYP-D3BJ/def2-SVPD with the stability analysis, reactant and transition state:
        the bare substrate, D34's three RELIABLE arrangements, D33 stage B's two designs

Usage:  cs6_d40_prepare.py <results-repo> <out-dir> <hpc-dir>
"""
import hashlib, os, shutil, sys

S6 = "phase2.2/cs_stage6"
GEOM = {"R": "d33a/sp_42436_R/substrate.xyz", "TS": "d33a/sp_42436_TS/substrate.xyz"}
DI, MONO = "-2 1", "-1 2"

# name, source run directory (relative to cs_stage6), state, files copied beside job.inp
CONTROL = ("C_CAM_R", "d28/CAM-B3LYP/def2-SVPD/ctl_R", "R", ["substrate.xyz"])
PART_A = [
    ("A_B3LYP",   "d27/def2-SVPD/ctl_R",            "R", ["substrate.xyz"]),
    ("A_CAM",     "d28/CAM-B3LYP/def2-SVPD/ctl_R",  "R", ["substrate.xyz"]),
    ("A_wB97XD3", "d28/wB97X-D3/def2-SVPD/ctl_R",   "R", ["substrate.xyz"]),
    ("A_wB97MV",  "d37/wB97MV_SV_R",                "R", ["substrate.xyz"]),
    ("A_wB97XV",  "d37/wB97XV_SV_R",                "R", ["substrate.xyz"]),
    ("A_PW6B95",  "d37/PW6B95_SV_R",                "R", ["substrate.xyz"]),
    ("A_M062X",   "d37/M062X_SV_R",                 "R", ["substrate.xyz"]),
    ("A_PWPB95",  "d37/PWPB95_SV_R",                "R", ["substrate.xyz"]),
]
PART_B = [("B_bare_%s" % st, "d31/ctl_%s" % st, st, ["substrate.xyz"]) for st in ("R", "TS")]
for arr in ("ADD-s540_qm06500", "ADD-s678_qm07625", "ADD-s678_qm10000"):
    PART_B += [("B_%s_%s" % (arr, st), "d34/a2_%s_%s" % (arr, st), st, ["combined.xyz", "site.ORCAFF.prms"]) for st in ("R", "TS")]
for des in ("bare_n05_q0.50", "bare_n11_q0.25"):
    PART_B += [("B_%s_%s" % (des, st), "d33b/%s_%s" % (des, st), st, ["substrate.xyz", "pc.pc"]) for st in ("R", "TS")]
# the order of running: the control; the two bare runs that reproduce or anchor; the production-level
# runs; then the six functionals whose open-shell keywords this record has not run before
ORDER = ([CONTROL[0], "A_B3LYP", "A_CAM"] + [x[0] for x in PART_B] +
         [x[0] for x in PART_A if x[0] not in ("A_B3LYP", "A_CAM")])
NRUN = 1 + len(PART_A) + len(PART_B)
C5_MONO = "qscan/VDE/C5_mono_R/job.inp"        # Probe 1's committed monoanion input at B3LYP-D3BJ/def2-SVPD
WALL, MEM = "24:00:00", "48gb"


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
        stop("usage: cs6_d40_prepare.py <results-repo> <out-dir> <hpc-dir>")
    res, out, hpc = sys.argv[1], sys.argv[2], sys.argv[3]
    s6 = os.path.join(res, S6)
    if os.path.exists(out) and os.listdir(out):
        stop("%s exists and is not empty; this script writes into a new directory only." % out)
    if len(ORDER) != NRUN or len(set(ORDER)) != NRUN:
        stop("the run order names %d runs, expected %d." % (len(set(ORDER)), NRUN))

    # ---- every source, before anything is written
    for name, src, st, extra in [CONTROL] + PART_A + PART_B:
        d = os.path.join(s6, src)
        for fn in ["job.inp", "job.out"] + extra:
            if not os.path.isfile(os.path.join(d, fn)):
                stop("%s/%s is missing; %s is built from it." % (src, fn, name))
        t = open(os.path.join(d, "job.inp")).read()
        if sum(1 for l in t.split("\n") if l.startswith("* xyzfile ")) != 1 or ("* xyzfile %s " % DI) not in t:
            stop("%s/job.inp does not have exactly one geometry line at charge -2, multiplicity 1." % src)
        if "****ORCA TERMINATED NORMALLY****" not in open(os.path.join(d, "job.out"), errors="replace").read():
            stop("%s/job.out, the committed dianion this run is read against, did not terminate normally." % src)
        if "substrate.xyz" in extra and sha(os.path.join(d, "substrate.xyz")) != sha(os.path.join(s6, GEOM[st])):
            stop("%s/substrate.xyz is not frame 42436's committed %s geometry." % (src, st))

    os.makedirs(out, exist_ok=True)
    sums, first = [], {}

    def place(name, src, extra, mono):
        d = os.path.join(out, name)
        os.makedirs(d)
        sd = os.path.join(s6, src)
        for fn in extra:
            shutil.copyfile(os.path.join(sd, fn), os.path.join(d, fn))
            if sha(os.path.join(sd, fn)) != sha(os.path.join(d, fn)):
                stop("%s/%s did not copy byte for byte into %s" % (src, fn, name))
        dst = os.path.join(d, "job.inp")
        if not mono:
            shutil.copyfile(os.path.join(sd, "job.inp"), dst)
            if sha(os.path.join(sd, "job.inp")) != sha(dst):
                stop("%s/job.inp did not copy byte for byte into %s" % (src, name))
        else:
            lines = open(os.path.join(sd, "job.inp")).read().split("\n")
            new = [l.replace("* xyzfile %s " % DI, "* xyzfile %s " % MONO) if l.startswith("* xyzfile ") else l for l in lines]
            with open(dst, "w") as f:
                f.write("\n".join(new))
        for fn in sorted(os.listdir(d)):
            sums.append((sha(os.path.join(d, fn)), "%s/%s" % (name, fn)))
        first[name] = open(dst).read().splitlines()[0]

    spec = {x[0]: x for x in [CONTROL] + PART_A + PART_B}
    for name in ORDER:
        _n, src, st, extra = spec[name]
        place(name, src, extra, mono=(name != CONTROL[0]))

    # ---- the comparability check, on what was actually written
    for name in ORDER:
        _n, src, st, extra = spec[name]
        a = open(os.path.join(s6, src, "job.inp")).read().split("\n")
        b = open(os.path.join(out, name, "job.inp")).read().split("\n")
        diff = [(x, y) for x, y in zip(a, b) if x != y]
        if len(a) != len(b):
            stop("%s: the input has a different number of lines from its committed source." % name)
        if name == CONTROL[0]:
            if diff:
                stop("%s is not D28's input byte for byte." % name)
            continue
        if len(diff) != 1:
            stop("%s differs from %s/job.inp in %d lines; exactly one is allowed." % (name, src, len(diff)))
        x, y = diff[0]
        if not (x.startswith("* xyzfile %s " % DI) and y == x.replace(DI, MONO, 1)):
            stop("%s: the one changed line is not the charge-and-multiplicity line: %r -> %r" % (name, x, y))
        if sorted(os.listdir(os.path.join(out, name))) != sorted(["job.inp"] + extra):
            stop("%s holds files other than its input and the %d files copied beside it." % (name, len(extra)))
    # the B3LYP monoanion must be Probe 1's committed monoanion input, to the byte
    c5 = os.path.join(s6, C5_MONO)
    if not os.path.isfile(c5) or sha(c5) != sha(os.path.join(out, "A_B3LYP", "job.inp")):
        stop("A_B3LYP/job.inp is not byte-identical to the committed %s." % C5_MONO)
    for name, src, st, extra in PART_B:
        t = open(os.path.join(out, name, "job.inp")).read()
        if "STABPerform true" not in t or first[name] not in ("! CAM-B3LYP D3BJ def2-SVPD def2/J RIJCOSX TightSCF",
                                                              "! CAM-B3LYP D3BJ def2-SVPD def2/J RIJCOSX TightSCF QMMM EnGrad"):
            stop("%s is not a CAM-B3LYP-D3BJ/def2-SVPD input with the stability analysis." % name)
    if len(set(first[x[0]] for x in PART_A)) != len(PART_A):
        stop("two runs of part A carry the same method line.")

    pbs = os.path.join(out, "run.pbs")
    with open(pbs, "w") as f:
        f.write("#!/bin/bash\n#PBS -N cm_d40\n")
        f.write("#PBS -l select=1:ncpus=8:mem=%s\n#PBS -l walltime=%s\n" % (MEM, WALL))
        f.write("#PBS -m ae\n#PBS -M 18660916@sun.ac.za\n#PBS -j oe\n")
        f.write("#PBS -o %s/run.pbs.out\n" % hpc)
        f.write("# Stage 6 DIAGNOSTIC 40: detachment energies at the density-functional level.\n")
        f.write("# C_CAM_R is D28's own CAM-B3LYP dianion input byte for byte and is the CONTROL (G40.4).\n")
        f.write("# Every other run is a monoanion, charge -1 multiplicity 2: a committed dianion input\n")
        f.write("#   with that one line changed. A_B3LYP is Probe 1's committed monoanion input and must\n")
        f.write("#   return its committed detachment energy (G40.5).\n")
        f.write("# Order: the control, the two anchors, the production-level runs, then the six\n")
        f.write("#   functionals whose open-shell runs this record has not made before.\n")
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
        for d in ORDER:
            f.write("run %s\n" % d)
        f.write('echo "end=$(date)"\n')
        f.write('echo "finished: $(grep -l \'ORCA TERMINATED NORMALLY\' */job.out 2>/dev/null | wc -l) of %d"\n' % NRUN)
    os.chmod(pbs, 0o755)
    sums.append((sha(pbs), "run.pbs"))

    with open(os.path.join(out, "INPUT_CHECKSUMS.sha256"), "w") as f:
        for h, n in sums:
            f.write("%s  %s\n" % (h, n))
    for name in ORDER:
        print("   %-26s %s" % (name, first[name]))
    print("   %d run directories, %d files checksummed, %s, walltime %s" % (NRUN, len(sums), MEM, WALL))


if __name__ == "__main__":
    main()
