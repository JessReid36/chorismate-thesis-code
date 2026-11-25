#!/usr/bin/env python3
"""cs6_d60_prepare.py - form the inputs of STAGE 6, DIAGNOSTIC 60 from committed files.

AMENDMENT 13, A13.4 items 7 and 8, as AMENDMENT 14 has them. Four single points of the substrate alone (the
dianion, frame 42436), each followed by orca_vpot at the 2626 sites of d33a/points.xyz:
  WXV_R, WXV_TS   wB97X-V/def2-SVPD: the inputs and geometries of D37's committed d37/wB97XV_SV_<state>, byte
                  for byte; the second functional's map is V_TS - V_R at the sites
  CAM_R, CAM_TS   CAM-B3LYP-D3BJ/def2-SVPD: the inputs and geometries of D33 stage A's committed
                  d33a/sp_42436_<state>, byte for byte; the control, which must give again D33's potentials
This program writes the inputs; it computes nothing and reads no output.

usage:
   python3 cs6_d60_prepare.py --res <results repository> --out <folder to write d60 into>
The reader cs6_d60.py calls build() and stops unless the committed d60/ is what build() gives.
"""
import argparse, hashlib, os, sys

S6REL = "phase2.2/cs_stage6"
DIR = "d60"
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/" + DIR
STATES = ("R", "TS")
SRC = {"WXV": "d37/wB97XV_SV_%s", "CAM": "d33a/sp_42436_%s"}
POINTS = "d33a/points.xyz"
CAP = 1800                     # s, D56's cap at def2-SVPD
JOB = "V1"


def rd(p):
    with open(p) as f:
        return f.read()


def runs():
    """[(run, cap)] in the order the job makes them"""
    return [("%s_%s" % (lv, st), CAP) for lv in ("WXV", "CAM") for st in STATES]


def source(run):
    lv, st = run.split("_")
    return SRC[lv] % st


def build(S6):
    """every file of d60/, as {relative path: text}"""
    F = {}
    for r, c in runs():
        s = source(r)
        F[r + "/job.inp"] = rd(os.path.join(S6, s, "job.inp"))
        F[r + "/substrate.xyz"] = rd(os.path.join(S6, s, "substrate.xyz"))
    F["points.xyz"] = rd(os.path.join(S6, POINTS))
    F["RUNS.tsv"] = "run\tjob\tcap_s\tsource\n" + "".join("%s\t%s\t%d\t%s\n" % (r, JOB, c, source(r)) for r, c in runs())
    F["run_%s.pbs" % JOB] = pbs(JOB, runs())
    F["INPUT_CHECKSUMS.sha256"] = "".join("%s  %s\n" % (hashlib.sha256(F[p].encode()).hexdigest(), p) for p in sorted(F))
    return F


def pbs(j, mine):
    t = """#!/bin/bash
#PBS -N cm_d60_%s
#PBS -l select=1:ncpus=8:mem=48gb
#PBS -l walltime=24:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %s/run_%s.pbs.out
# Stage 6 DIAGNOSTIC 60: the map at a second functional, wB97X-V/def2-SVPD, and its control (AMENDMENTS 13 and 14).
# Job %s of 1. It makes %d runs in the order below (name, cap in s), each followed by orca_vpot at the 2626 sites; a run
# that fails does not stop the ones after it.
# No -q: the allocation is automatic and walltime 24:00:00 routes to day. The job script is D58's committed run_E1.pbs
# with its names, runs and caps changed. No run is made twice:
#   submitted again, it makes only the runs whose folder holds no output yet (D44 NOTE 3).
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
ORCA=/home/apps2/ORCA/6.0.1
W=%s
cd "$W" || exit 1
if ! mkdir .running_%s 2>/dev/null; then echo "FAIL: another job holds .running_%s in $W"; exit 1; fi
trap 'rmdir "$W/.running_%s" 2>/dev/null' EXIT
echo "host=$(hostname) start=$(date)"
echo "free memory at start:"; free -g | head -2
run() {
  d="$W/$1"
  if [ -e "$d/job.out" ]; then echo "$1 has an output already: not run again"
  else ( cd "$d" && timeout -k 10 "$2" "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); echo "$1 orca_exit=$? normal=$(grep -c 'ORCA TERMINATED NORMALLY' "$d/job.out") $(date)"; fi
}
vpot() {
  d="$W/$1"; pts="$W/points.xyz"
  if [ -e "$d/vpot.log" ]; then echo "$1 vpot has a log already: not run again"; return; fi
  if ! grep -q "ORCA TERMINATED NORMALLY" "$d/job.out" 2>/dev/null; then echo "$1 vpot skipped, run not normal"; return; fi
  ( cd "$d" && "$ORCA/orca_vpot" job.gbw job.scfp "$pts" vpot.out job > vpot.log 2>&1 )
  echo "$1 vpot_exit=$? lines=$(wc -l < "$d/vpot.out" 2>/dev/null || echo 0) expected=2627"
}
""" % (j.lower(), HPC, j, j, len(mine), HPC, j, j, j)
    for r, c in mine:
        t += "run %s %d\n" % (r, c)
        t += "vpot %s\n" % r
    t += 'echo "end=$(date)"\n'
    t += "n=0; for r in %s; do grep -q 'ORCA TERMINATED NORMALLY' \"$W/$r/job.out\" 2>/dev/null && n=$((n+1)); done\n" % " ".join(r for r, _ in mine)
    t += 'echo "finished: $n of %d"\n' % len(mine)
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    S6 = os.path.join(a.res, S6REL)
    F = build(S6)
    dst = os.path.join(a.out, DIR)
    if os.path.exists(dst):
        sys.exit("REFUSED: %s exists already." % dst)
    for p in sorted(F):
        os.makedirs(os.path.dirname(os.path.join(dst, p)) or dst, exist_ok=True)
        with open(os.path.join(dst, p), "w") as f:
            f.write(F[p])
    print("wrote %d files into %s: %d runs in 1 job" % (len(F), dst, len(runs())))


if __name__ == "__main__":
    main()
