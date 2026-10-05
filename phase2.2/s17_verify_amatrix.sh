#!/usr/bin/env bash
#PBS -N cm_s17chk
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -l walltime=04:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o /home/18660916/system_development/phase2.2/s17_check/s17_check.pbs.out
#
# s17_verify_amatrix.sh - does the linear-response a-matrix agree with direct DFT?
#
# WHY THIS IS NEEDED
# The first certified design puts a +1 at site 117 and the a-matrix says that is worth
# -9.15 kcal/mol averaged over thirty frames. Three independent numbers say that is
# large:
#     probe scan, +1 at 2.8 A on frame 20000      -2.90 kcal/mol
#     probe scan, +1 at 3.2 A, mean of 3 frames   -1.53 kcal/mol
#     Burschowsky's real arginine, experiment     -5.30 kcal/mol
# Site 117 is 2.80 A from the ether oxygen, so distance alone does not explain a factor
# of three against the probe at the same distance. Direction might: the probe ran along
# one specific ray while the grid site is wherever the surface put it. But a factor of
# three has to be demonstrated, not assumed.
#
# WHAT THIS MEASURES
# For each test frame, the differential stabilisation computed the way the design
# claims it, directly:
#     ddE = [E_TS(with charge) - E_TS(bare)] - [E_R(with charge) - E_R(bare)]
# at exactly the position the optimiser chose, with the charge mapped back into that
# frame's own coordinates by the inverse alignment transform. The bare energies already
# exist in 20_invacuo and are reused, so only the two field calculations per frame are
# new.
#
# HOW TO READ THE RESULT
#   agreement within a few per cent -> the linear model is sound and -9.15 is a real
#     geometric result: the optimiser found a better direction than the probe ray.
#   DFT much smaller than predicted -> the a-matrix is wrong, and the likely causes are
#     the reference density, a sign or unit convention in the potential, or the grid
#     point not being where the transform says it is.
# Either way the answer is decided by one number per frame, not by argument.
#
# NOTE ON min_approach. Every site in the first design sits at exactly 2.50 A from the
# nearest atom, which is the min_approach floor. The optimiser pushes every charge as
# close as the grid allows, so that parameter is setting the answer and it was chosen
# as a placeholder. That is a separate issue from this check and is recorded in the plan.

set -uo pipefail
ORCA=/home/apps2/ORCA/6.0.1
ROOT=/home/18660916/system_development
VAC="$ROOT/05_qmmm/20_invacuo"
WORK="$ROOT/phase2.2/s17_check"
SITE=117
FRAMES="20000 21634 23268 24883 26495"

export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="$ORCA/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

mkdir -p "$WORK"; cd "$WORK" || exit 1
if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running"; exit 1; fi
SCRATCH="${TMPDIR:-/tmp}/orca_$$"; mkdir -p "$SCRATCH"
trap 'rm -rf "$SCRATCH"; rmdir "$WORK/.running" 2>/dev/null' EXIT
echo "host=$(hostname) start=$(date)"
echo "site $SITE, frames: $FRAMES"
echo

# write the point-charge files, mapped per frame
python3 - <<PYEOF
import numpy as np, os
ROOT = "$ROOT"
site = $SITE
frames = "$FRAMES".split()
g = [l.split("\t") for l in open(f"{ROOT}/phase2.2/grid_v2.tsv")
     if not l.startswith(("#", "idx"))]
xyz = np.array([[float(c) for c in r[1:4]] for r in g])
for fr in frames:
    M = np.loadtxt(f"{ROOT}/phase2.2/aligned/frame_{fr}/transform.txt")
    R, t = M[:3], M[3]
    p = (xyz[site] - t) @ R.T
    open(f"site{site}_{fr}.pc", "w").write(f"1\n1.0000 {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}\n")
print(f"  wrote {len(frames)} point-charge files for site {site}")
PYEOF
echo

ok=0; bad=0
for f in $FRAMES; do
  for tag in R TS; do
    base="${f}_${tag}"
    out="field_${SITE}_${base}.out"
    [ -s "$out" ] && grep -q "TERMINATED NORMALLY" "$out" && { ok=$((ok+1)); continue; }
    [ -s "$VAC/${base}.xyz" ] || { echo "  MISSING $VAC/${base}.xyz"; bad=$((bad+1)); continue; }

    d="$SCRATCH/$base"; rm -rf "$d"; mkdir -p "$d"
    cp "$VAC/${base}.xyz" "site${SITE}_${f}.pc" "$d/"
    cat > "$d/in.inp" <<EOF
! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF
%maxcore 3000
%pal nprocs 8 end
%scf MaxIter 300 end
%pointcharges "site${SITE}_${f}.pc"
* xyzfile -2 1 ${base}.xyz
EOF
    ( cd "$d" && "$ORCA/orca" in.inp > in.out 2>&1 </dev/null )
    cp "$d/in.out" "$out" 2>/dev/null
    rm -rf "$d"
    if grep -q "TERMINATED NORMALLY" "$out" 2>/dev/null; then ok=$((ok+1))
    else bad=$((bad+1)); echo "  FAILED $base"; fi
  done
done
echo "ran=$ok failed=$bad"
echo

python3 - <<PYEOF
import re, numpy as np
ROOT = "$ROOT"
site = $SITE
frames = "$FRAMES".split()
H = 627.5094740631

def energy(p):
    try:
        t = open(p, errors="replace").read()
    except OSError:
        return None
    if "TERMINATED NORMALLY" not in t:
        return None
    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    return float(m[-1]) if m else None

A = np.array([[float(x) for x in l.split()]
              for l in open(f"{ROOT}/phase2.2/A_v2.tsv") if not l.startswith("#")])
allfr = [x for x in open(f"{ROOT}/phase2.2/s15_dv2/frames.txt").read().split() if x]

print(f"{'frame':>8}{'linear model':>15}{'direct DFT':>13}{'difference':>13}{'ratio':>9}")
rows = []
for fr in frames:
    eR  = energy(f"{ROOT}/05_qmmm/20_invacuo/sp_{fr}_R.out")
    eT  = energy(f"{ROOT}/05_qmmm/20_invacuo/sp_{fr}_TS.out")
    fR  = energy(f"field_{site}_{fr}_R.out")
    fT  = energy(f"field_{site}_{fr}_TS.out")
    if None in (eR, eT, fR, fT):
        print(f"{fr:>8}   incomplete")
        continue
    dft = ((fT - eT) - (fR - eR)) * H
    lin = A[site, allfr.index(fr)]
    rows.append((fr, lin, dft))
    print(f"{fr:>8}{lin:>15.3f}{dft:>13.3f}{dft-lin:>13.3f}"
          f"{(dft/lin if lin else float('nan')):>9.3f}")

if rows:
    L = np.array([r[1] for r in rows]); D = np.array([r[2] for r in rows])
    rel = np.abs((D - L) / L) * 100
    print()
    print(f"  mean linear {L.mean():+.3f}, mean DFT {D.mean():+.3f} kcal/mol")
    print(f"  relative difference: min {rel.min():.1f}%, mean {rel.mean():.1f}%, "
          f"max {rel.max():.1f}%")
    print()
    if rel.mean() < 10:
        print("  AGREES. The linear-response model reproduces direct DFT at this site,")
        print("  so the magnitude is a real geometric result and not an artefact: the")
        print("  optimiser found a more effective direction than the probe ray sampled.")
    elif rel.mean() < 30:
        print("  PARTIAL AGREEMENT. The linear model is in the right range but the")
        print("  discrepancy exceeds what second-order polarisation would explain")
        print("  (measured at 20 per cent at q = 1). Worth tracing before quoting.")
    else:
        print("  DISAGREES. The a-matrix does not reproduce direct DFT and must not be")
        print("  used for a reported design until the cause is found. Check, in order:")
        print("    - the reference density: in vacuo, as the design environment requires")
        print("    - the sign and unit convention in the potential")
        print("    - whether the mapped point is where the transform says it is")
PYEOF
echo
echo "end=$(date)"
