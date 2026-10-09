#!/usr/bin/env python3
"""
cs4c_amend1_prepare.py - Stage 4c, Amendment 1 (STAGE4C_AMENDMENT1.txt). Writes the zero-charge runs that remove the
pseudopotential's own (charge-independent) reactant/TS energy from the pseudopotential levels' ddE, beside the
committed inputs and outputs in phase2.2/cs_stage4c/ (nothing existing is touched).
Each zero-charge input is the committed input of the arrangement named below with every Ne> charge set to 0.000000
and nothing else changed (method line, substrate lines, sphere positions and ECP blocks identical):
  site set   from   sites                      used for
  S1         E1     117                        E1, E2
  S3         E3     232                        E3
  S4         E4     189                        E4
  S5         E5     117, 232, 60               E5, E7
  S6         E6     117, 232, 60, 88, 228      E6
  S8         E8     117, 60, 16                E8
Levels L0 (def2-SVP), L2 (def2-SVPD), L3 (def2-TZVPD); R and TS: 36 runs, Z_<level>_<set>_<R|TS>.
USAGE  python3 cs4c_amend1_prepare.py <cs_stage4c folder>
"""
import os, re, sys
D = sys.argv[1]
SD = "/home/18660916/system_development/phase2.2/cs_stage4c"
SETS = (("S1", "E1"), ("S3", "E3"), ("S4", "E4"), ("S5", "E5"), ("S6", "E6"), ("S8", "E8"))
SHARED = (("E2", "E1"), ("E7", "E5"))
assert os.path.exists(os.path.join(D, "STAGE4C_CRITERIA.txt")), "not the cs_stage4c folder"
NSPH = {}
for l in open(os.path.join(D, "arrangements.tsv")):
    if not l.startswith("#") and l.strip():
        NSPH[l.split("\t")[0]] = NSPH.get(l.split("\t")[0], 0) + 1
NE = re.compile(r"^(  Ne>  )(-?\d+\.\d+)( .*)$")
def ne_lines(path):
    return [l for l in open(path).read().splitlines() if NE.match(l)]
for L in ("L0", "L2", "L3"):          # sets shared between arrangements must be the same positions
    for g in ("R", "TS"):
        for e, f in SHARED:
            a = [NE.match(l).group(3) for l in ne_lines(os.path.join(D, "%s_%s_%s" % (L, e, g), "job.inp"))]
            b = [NE.match(l).group(3) for l in ne_lines(os.path.join(D, "%s_%s_%s" % (L, f, g), "job.inp"))]
            assert a == b, "%s and %s do not share sphere positions" % (e, f)
batches = {"L0": [], "L2": [], "L3": []}
for L in ("L0", "L2", "L3"):
    for s, e in SETS:
        for g in ("R", "TS"):
            src = open(os.path.join(D, "%s_%s_%s" % (L, e, g), "job.inp")).read().splitlines()
            new, n = [], 0
            for l in src:
                m = NE.match(l)
                if m:
                    l = m.group(1) + "0.000000" + m.group(3); n += 1
                new.append(l)
            assert n == NSPH[e], "%s_%s_%s: %d Ne> lines, expected %d" % (L, e, g, n, NSPH[e])
            assert sum(1 for x, y in zip(src, new) if x != y) == n and len(src) == len(new)
            name = "Z_%s_%s_%s" % (L, s, g)
            os.makedirs(os.path.join(D, name))
            open(os.path.join(D, name, "job.inp"), "w").write("\n".join(new) + "\n")
            batches[L].append(name)
for L, b in (("L0", "svp"), ("L2", "svpd"), ("L3", "tzvpd")):
    with open(os.path.join(D, "batch_Z_%s.pbs" % b), "w") as f:
        f.write("""#!/bin/bash
#PBS -N cm34c_Z%s
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -l walltime=168:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %s/batch_Z_%s.pbs.out
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
echo "host=$(hostname) start=$(date)"
""" % (b, SD, b))
        f.write("for d in %s; do\n  cd %s/$d\n  /home/apps2/ORCA/6.0.1/orca job.inp > job.out 2>&1 </dev/null\n" % (" ".join(batches[L]), SD))
        f.write('  echo "$d orca_exit=$? $(grep -c \'ORCA TERMINATED NORMALLY\' job.out) $(date)"\n  rm -f job.gbw job*.tmp 2>/dev/null\ndone\necho "end=$(date)"\n')
print({L: len(r) for L, r in batches.items()}, "zero-charge runs written")
