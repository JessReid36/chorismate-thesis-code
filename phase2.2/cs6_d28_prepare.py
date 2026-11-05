#!/usr/bin/env python3
"""DIAGNOSTIC 28: the same conditions as DIAGNOSTIC 27, with two range-separated functionals.

Every input is generated from the SAME committed geometries and charge files as d27, with the SAME
%maxcore, %pal, %scf MaxIter, basis, def2/J, RIJCOSX, TightSCF, charge and multiplicity. The ONLY
difference from the corresponding d27 input is the functional and its dispersion token. The checker
verifies that against d27's own committed inputs rather than taking it on trust.
"""
import argparse, hashlib, shutil, sys
from pathlib import Path

S6 = "phase2.2/cs_stage6"
# keyword line fragment per functional. wB97X-D3 carries its own D3, so NO separate D3BJ:
# adding one would double-count. CAM-B3LYP has no built-in dispersion, so D3BJ is added.
FUNCS = {"wB97X-D3": "wB97X-D3", "CAM-B3LYP": "CAM-B3LYP D3BJ"}
VECTORS = ["1.3000e_n10", "2.8800e_n10", "5.6433e_n20"]
STATES = {"R": "reactant_qm.xyz", "TS": "transition_state_qm.xyz"}
REMOTE = "/home/18660916/system_development/" + S6 + "/d28"

INP = """! {kw} {basis} def2/J RIJCOSX TightSCF
%maxcore 3000
%pal nprocs 8 end
%scf MaxIter 300 end
{pc}* xyzfile -2 1 substrate.xyz
"""
# a few-second job whose only purpose is to prove ORCA accepts the keyword before 20 real jobs run
PROBE = """! {kw} def2-SVP def2/J RIJCOSX TightSCF
%maxcore 1000
%pal nprocs 2 end
* xyz 0 1
O   0.000000   0.000000   0.117300
H   0.000000   0.757200  -0.469200
H   0.000000  -0.757200  -0.469200
*
"""

PBS = """#!/bin/bash
#PBS -N cm_d28
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -l walltime=24:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o {remote}/run.pbs.out
# Stage 6 DIAGNOSTIC 28: {n} single points, two range-separated functionals. Restartable.
# No -q is given: the allocation is automatic, and walltime 24:00:00 routes through workq to day.
# Order is deliberate: keyword probes first, then the controls that carry the primary reading,
# then the charged arms. An early look answers the main question.
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


def write(d, kw, basis, geom, pc, src):
    d.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src / geom, d / "substrate.xyz")
    if pc:
        shutil.copy2(src / pc, d / "pc.pc")
    (d / "job.inp").write_text(INP.format(kw=kw, basis=basis,
                                          pc='%pointcharges "pc.pc"\n' if pc else ""))


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
    for fn, kw in FUNCS.items():                                   # the keyword probes
        d = out / "probe" / fn; d.mkdir(parents=True)
        (d / "job.inp").write_text(PROBE.format(kw=kw))
        names.append("probe/%s" % fn)
    for basis in ("def2-SVPD", "def2-TZVPD"):                      # the controls, both bases
        for fn, kw in FUNCS.items():
            for st, geom in STATES.items():
                write(out / fn / basis / ("ctl_%s" % st), kw, basis, geom, None, src)
                names.append("%s/%s/ctl_%s" % (fn, basis, st))
    for fn, kw in FUNCS.items():                                   # charged arms, def2-SVPD only
        for v in VECTORS:
            for st, geom in STATES.items():
                write(out / fn / "def2-SVPD" / ("%s_%s" % (v, st)), kw, "def2-SVPD", geom,
                      "pc_%s.pc" % v, src)
                names.append("%s/def2-SVPD/%s_%s" % (fn, v, st))

    (out / "run.pbs").write_text(PBS.format(remote=REMOTE, n=len(names),
                                            runs="\n".join("run " + n for n in names)))
    lines = []
    for p in sorted(out.rglob("*")):
        if p.is_file(): lines.append("%s  %s" % (sha(p), p.relative_to(out)))
    (out / "INPUT_CHECKSUMS.sha256").write_text("\n".join(lines) + "\n")
    print("built %d jobs in %s" % (len(names), out))
    print("   %-12s %d keyword probes" % ("probe", sum(1 for n in names if n.startswith("probe"))))
    for fn in FUNCS:
        print("   %-12s %d runs" % (fn, sum(1 for n in names if n.startswith(fn + "/"))))
    print("   run.pbs, INPUT_CHECKSUMS.sha256")
    print()
    print("geometries and charge files copied from the committed d24 inputs, the same files d27 used:")
    for f in list(STATES.values()) + ["pc_%s.pc" % v for v in VECTORS]:
        print("   %s  %s" % (sha(src / f)[:16], f))
    return 0


if __name__ == "__main__":
    sys.exit(main())
