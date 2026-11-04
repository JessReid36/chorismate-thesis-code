#!/usr/bin/env python3
"""DIAGNOSTIC 27 -- build the 24 single points that settle the basis question.

Three bases x three charge vectors x two states, plus a control per basis. Generates
inputs only; submits nothing. Every geometry and every charge file is COPIED from the
committed DIAGNOSTIC 24 inputs, never regenerated, so the vectors are identical by
construction and not by coincidence.
"""
import argparse, hashlib, shutil, sys
from pathlib import Path

S6 = "phase2.2/cs_stage6/"
BASES = ["def2-SVP", "def2-SVPD", "def2-TZVPD"]
VECTORS = ["1.3000e_n10", "2.8800e_n10", "5.6433e_n20"]
STATES = {"R": "reactant_qm.xyz", "TS": "transition_state_qm.xyz"}
REMOTE = "/home/18660916/system_development/phase2.2/cs_stage6/d27"

INP = """! B3LYP D3BJ {basis} def2/J RIJCOSX TightSCF
%maxcore 3000
%pal nprocs 8 end
%scf MaxIter 300 end
{pc}* xyzfile -2 1 substrate.xyz
"""

PBS = """#!/bin/bash
#PBS -N cm_d27
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -l walltime=24:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o {remote}/run.pbs.out
# Stage 6 DIAGNOSTIC 27: {n} single points over three bases. Restartable.
# No -q is given: the allocation is automatic, and walltime 24:00:00 routes through workq to day.
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${{LD_LIBRARY_PATH:-}}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
ORCA=/home/apps2/ORCA/6.0.1
W={remote}
cd "$W" || exit 1
if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running in $W"; exit 1; fi
trap 'rmdir "$W/.running" 2>/dev/null' EXIT
echo "host=$(hostname) start=$(date)"
run() {{
  d="$W/$1"
  if [ -s "$d/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$d/job.out"; then echo "$1 done before"
  else ( cd "$d" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); echo "$1 orca_exit=$? normal=$(grep -c 'ORCA TERMINATED NORMALLY' "$d/job.out") $(date)"; fi
}}
{runs}
echo "end=$(date)"
"""


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    src = Path(a.results) / S6 / "d24"
    out = Path(a.out)
    if out.exists(): print("STOP  %s already exists." % out); return 2
    for f in list(STATES.values()) + ["pc_%s.pc" % v for v in VECTORS]:
        if not (src / f).exists(): print("STOP  committed input %s is missing." % f); return 2

    names = []
    # cheapest basis first, so an early look already answers something
    for basis in BASES:
        for st in ("R", "TS"):
            d = out / basis / ("ctl_" + st); d.mkdir(parents=True)
            shutil.copy2(src / STATES[st], d / "substrate.xyz")
            (d / "job.inp").write_text(INP.format(basis=basis, pc=""))
            names.append("%s/ctl_%s" % (basis, st))
        for v in VECTORS:
            for st in ("R", "TS"):
                d = out / basis / ("%s_%s" % (v, st)); d.mkdir(parents=True)
                shutil.copy2(src / STATES[st], d / "substrate.xyz")
                shutil.copy2(src / ("pc_%s.pc" % v), d / "pc.pc")
                (d / "job.inp").write_text(INP.format(basis=basis, pc='%pointcharges "pc.pc"\n'))
                names.append("%s/%s_%s" % (basis, v, st))

    (out / "run.pbs").write_text(PBS.format(remote=REMOTE, n=len(names),
                                            runs="\n".join("run " + n for n in names)))
    lines = []
    for p in sorted(out.rglob("*")):
        if p.is_file(): lines.append("%s  %s" % (sha(p), p.relative_to(out)))
    (out / "INPUT_CHECKSUMS.sha256").write_text("\n".join(lines) + "\n")

    print("built %d single points in %s" % (len(names), out))
    for basis in BASES:
        print("   %-11s %d runs" % (basis, sum(1 for n in names if n.split("/")[0] == basis)))
    print("   run.pbs, INPUT_CHECKSUMS.sha256")
    print()
    print("geometries and charge files copied from the committed d24 inputs:")
    for f in list(STATES.values()) + ["pc_%s.pc" % v for v in VECTORS]:
        print("   %s  %s" % (sha(src / f)[:16], f))
    return 0


if __name__ == "__main__":
    sys.exit(main())
