#!/usr/bin/env python3
"""
cs6_diffuse.py - the diffuse-basis diagnostic of STAGE6_DIAGNOSTIC1_20261022.txt.

  python3 cs6_diffuse.py <relax_dir> <OUT> <hpc_dir>

A DIAGNOSTIC, NOT A METHOD CHANGE. STAGE6_CRITERIA.txt fixes B3LYP-D3BJ/def2-SVP for everything Stage 6 reports and
that stays fixed; nothing these runs produce replaces a reported number. They exist because
STAGE6_FINDINGS_20261022.md Finding 8 showed the bare dianion has a HOMO of +2.28 eV at def2-SVP - the highest
electrons are not bound, and the calculation only looks well-behaved because the basis has no diffuse functions to
let one leave. An applied field acts on loosely-held density exactly where that artefact bites hardest, so part of
what Finding 1 reads as polarisation may be the basis holding unbound density differently in different fields.

The runs are the committed frozen single points of four arrangements at both geometries, byte-identical except for
the basis set: def2-SVP -> ma-def2-SVP, minimally augmented, which adds diffuse s and p functions on the heavy
atoms and is the standard cheap choice for anions. def2/J stays as the auxiliary, which is what ORCA documents for
the ma-def2 sets. Eight single points of a few minutes each, against the fourteen six-hour relaxations they are
meant to inform.

The four: S6-C (the bare dianion, the worst case for an unbound HOMO), S6-P (the only design that binds it), S6-2
(the smallest sphere pair) and S6-20 (the largest design and the one with the most added gradient). They span the
range of sum q*V_R from -120 to +44 kcal/mol.

Checks, each a stop:
  (a) every source folder exists and its three input files are the committed ones, by relax/INPUT_CHECKSUMS.sha256;
  (b) each generated job.inp differs from its source in EXACTLY the basis-set token and nothing else;
  (c) combined.xyz and site.ORCAFF.prms are copied byte for byte, so the geometry and the charges are identical and
      the basis is the only thing that changed;
  (d) the output folder does not exist, and nothing is written anywhere else.
"""
import hashlib, os, shutil, sys
from pathlib import Path

RELAX, OUT, HPC = sys.argv[1], sys.argv[2], sys.argv[3]
DESIGNS = ("S6-C", "S6-P", "S6-2", "S6-20")
STATES = ("R", "TS")
OLD, NEW = "def2-SVP", "ma-def2-SVP"
WALLTIME = "02:00:00"


def stop(m):
    sys.exit("STOP: " + m)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


relax, out = Path(RELAX), Path(OUT)
if out.exists():
    stop("%s exists" % out)
# (a) the sources are the committed inputs
rec = {}
for line in (relax / "INPUT_CHECKSUMS.sha256").read_text().splitlines():
    h, rel = line.split("  ", 1)
    rec[rel] = h
made = []
for name in DESIGNS:
    for st in STATES:
        src = relax / name / ("%s_sp" % st)
        if not src.is_dir():
            stop("%s is missing" % src)
        for fn in ("job.inp", "combined.xyz", "site.ORCAFF.prms"):
            f = src / fn
            if not f.is_file():
                if fn == "site.ORCAFF.prms" and name == "S6-C":
                    continue          # the control is a plain QM job and has no force field
                stop("%s is missing" % f)
            rel = str(f.relative_to(relax))
            if rel not in rec:
                stop("%s is not in relax/INPUT_CHECKSUMS.sha256" % rel)
            if sha(f) != rec[rel]:
                stop("%s has changed since it was generated" % rel)
        d = out / ("DIFF-%s" % name) / ("%s_sp" % st)
        d.mkdir(parents=True)
        inp = (src / "job.inp").read_text()
        if inp.count(OLD) != 1:
            stop("%s/job.inp contains %r %d times, expected once" % (src, OLD, inp.count(OLD)))
        new = inp.replace(OLD, NEW)
        # (b) exactly the basis token differs
        a, b = inp.splitlines(), new.splitlines()
        diff = [i for i in range(len(a)) if a[i] != b[i]]
        if len(a) != len(b) or len(diff) != 1 or a[diff[0]].replace(OLD, NEW) != b[diff[0]]:
            stop("%s/job.inp would differ in more than the basis set" % src)
        (d / "job.inp").write_text(new)
        # (c) everything else byte for byte
        for fn in ("combined.xyz", "site.ORCAFF.prms"):
            if (src / fn).is_file():
                shutil.copyfile(src / fn, d / fn)
                if sha(d / fn) != sha(src / fn):
                    stop("%s was not copied byte for byte" % fn)
        remote = "%s/DIFF-%s/%s_sp" % (HPC.rstrip("/"), name, st)
        (d / "job.pbs").write_text(
            "#!/bin/bash\n#PBS -N cm_s6d_%s_%s\n#PBS -l select=1:ncpus=8:mem=24gb\n#PBS -l walltime=%s\n"
            "#PBS -m ae\n#PBS -M 18660916@sun.ac.za\n#PBS -j oe\n#PBS -o %s/job.pbs.out\ncd %s\n"
            'export PATH="/apps/openmpi/4.1.1/bin:$PATH"\n'
            'export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:'
            '/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"\n'
            "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1\n"
            'echo "host=$(hostname) start=$(date)"\n'
            "/home/apps2/ORCA/6.0.1/orca job.inp > job.out 2>&1 </dev/null\n"
            'echo "orca_exit=$? end=$(date)"\n'
            "rm -f job.gbw job*.tmp 2>/dev/null || true\n"
            % (name.replace("S6-", ""), st, WALLTIME, remote, remote))
        made.append(str(d.relative_to(out)))
        print("   %-26s basis %s -> %s, geometry and charges byte-identical" % (d.relative_to(out), OLD, NEW))
lines = []
for root, _s, files in os.walk(out):
    for fn in sorted(files):
        p = Path(root) / fn
        lines.append("%s  %s" % (sha(p), p.relative_to(out)))
(out / "INPUT_CHECKSUMS.sha256").write_text("\n".join(sorted(lines, key=lambda l: l.split("  ", 1)[1])) + "\n")
print("   %d folders, %d files, with a manifest" % (len(made), len(lines)))
