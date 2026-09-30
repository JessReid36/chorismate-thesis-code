#!/usr/bin/env bash
#PBS -N cm_s13lj
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -l walltime=04:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o /home/18660916/system_development/phase2.2/s13_lj_test/s13_lj_test.pbs.out
#
# s13_lj_charge_site_test.sh - can a designed charge site carry a Pauli wall?
#
# THE PROBLEM, from step_7_relaxed_validation_and_point_charge_collapse.md
# Relaxed optimisations with designed charges as bare ORCA point charges collapse: a
# substrate oxygen migrates onto a designed +1 and ends INSIDE its own van der Waals
# radius, closest approach 0.879 A, 0.641 A beyond the vdW contact. Four of eight designs
# collapsed. The cause recorded there is the missing Pauli wall - a bare point charge is
# an unbounded 1/r attractor. Nothing in the standard output flags it; the converged run
# reported an energy 249 kcal/mol below what the frozen calculation implies.
#
# WHY THIS ROUTE WAS PREVIOUSLY REJECTED, AND WHY IT IS WORTH RE-TESTING
# step_7b rejected QM/MM with Lennard-Jones parameters because ORCA aborts with
#   "CPCM or SMD or ALPB or ddCOSMO or CPCMX requested together with QM/MM method.
#    This is not implemented."
# That is a real, quoted test result. But it is CONDITIONAL ON CPCM BEING REQUIRED, and
# the Phase 2.2 environment decision removes the continuum: [D2018] never runs charges
# and COSMO together, and [B2021] uses explicit molecules with no continuum at all.
# With CPCM gone the abort should not fire.
#
# step_7b also records that orca_mm -makeff "cannot generate a force field for
# free-floating charge sites: it runs an internal xtb geometry optimisation which fails
# on disconnected atoms". That is true of GENERATING a force field. The ORCA 6 manual
# (Molecular Mechanics) states that the ORCAFF.prms file's nonbonded parameters "can be
# conveniently modified directly within the ORCA Force Field File. This file can then be
# modified, the required values can be added, and the resulting file can be defined as
# input for the QMMM calculation." So the site is ADDED BY HAND to an existing file
# rather than generated.
#
# WHAT THIS TESTS, in order, stopping at the first failure
#   A. Does !QMMM run at all without CPCM on this system?
#   B. Is a hand-edited ORCAFF.prms with an added free-floating site accepted?
#   C. Under relaxed optimisation, where does a substrate oxygen stop against a +1 that
#      carries LJ parameters, compared with the 0.879 A of the bare point charge?
#
# THE CRITERION, fixed before running. step_7 establishes that inward motion is PARTLY
# PHYSICAL: a real guanidinium should pull a carboxylate into a salt bridge, and the
# structural consensus places Arg-carboxylate contacts at 2.6-3.0 A. So the target is not
# to prevent approach but to make it stop in the right place.
#     closest approach < 2.0 A   -> the wall is too weak, no better than bare
#     2.0 to 3.2 A               -> PHYSICAL, this is the salt-bridge range
#     > 4.0 A                    -> the wall is too hard, the charge cannot act
#
# NOT TESTED HERE: whether the LJ parameters chosen are the right ones. The sodium values
# are lifted verbatim from the committed complex_solvated.ORCAFF.prms that step13a
# generated, so they are sourced rather than invented - but a designed cation standing in
# for a guanidinium is not sodium-sized, and r_min 2.738 A is a choice this test does not
# justify. If the mechanism works, the parameter choice becomes the next question, and
# Burschowsky's 3.2 A citrulline contact is the obvious calibration target.

set -uo pipefail
ORCA=/home/apps2/ORCA/6.0.1
ROOT=/home/18660916/system_development
WORK="$ROOT/phase2.2/s13_lj_test"

export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="$ORCA/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

mkdir -p "$WORK"; cd "$WORK" || exit 1
if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running"; exit 1; fi
SCRATCH="${TMPDIR:-/tmp}/orca_$$"; mkdir -p "$SCRATCH"
trap 'rm -rf "$SCRATCH"; rmdir "$WORK/.running" 2>/dev/null' EXIT
echo "host=$(hostname)  start=$(date)  scratch=$SCRATCH"
echo

# ---------------------------------------------------------------- inputs
FRAME=24883
SRC="$ROOT/05_qmmm/19_ensemble_barriers/frame_$FRAME/reactant_qm.xyz"
[ -s "$SRC" ] || { echo "FAIL: $SRC missing"; exit 1; }
cp "$SRC" sub.xyz
echo "substrate: frame $FRAME reactant, $(head -1 sub.xyz) atoms"

# Place a +1 at 3.2 A from the ether oxygen O3 (index 7, 0-based), outward along the
# centroid->O3 ray. Same construction as s11_probe_distance_scan.py.
python3 - <<'PY'
import math
L=open('sub.xyz').read().splitlines(); n=int(L[0].split()[0])
at=[(f[0],float(f[1]),float(f[2]),float(f[3])) for f in (l.split() for l in L[2:2+n])]
cx=sum(a[1] for a in at)/n; cy=sum(a[2] for a in at)/n; cz=sum(a[3] for a in at)/n
ox,oy,oz=at[7][1:4]
vx,vy,vz=ox-cx,oy-cy,oz-cz
d=math.sqrt(vx*vx+vy*vy+vz*vz)
px,py,pz=ox+3.2*vx/d, oy+3.2*vy/d, oz+3.2*vz/d
open('probe.xyz','w').write(f"{px:.8f} {py:.8f} {pz:.8f}\n")
open('bare.pc','w').write(f"1\n1.0000 {px:.8f} {py:.8f} {pz:.8f}\n")
# combined structure: substrate then the site, as one xyz for the QM/MM run
with open('combined.xyz','w') as fh:
    fh.write(f"{n+1}\nsubstrate + one charge site\n")
    for e,x,y,z in at: fh.write(f"{e:<3}{x:15.8f}{y:15.8f}{z:15.8f}\n")
    fh.write(f"Na {px:15.8f}{py:15.8f}{pz:15.8f}\n")   # Na as the MM site carrier
print(f"  probe at {px:.3f} {py:.3f} {pz:.3f}, 3.200 A from O3")
PY
echo

# ------------------------------------------------- A: does QMMM run without CPCM?
echo "=== TEST A: !QMMM without CPCM"
cat > testA.inp <<'EOF'
! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF QMMM
%maxcore 3000
%pal nprocs 8 end
%qmmm
  QMAtoms {0:23} end
  ORCAFFFilename "site.ORCAFF.prms"
end
* xyzfile -2 1 combined.xyz
EOF

# Hand-write a minimal ORCAFF.prms: 25 atoms, the 24 QM ones plus the charge site.
# Columns after the index and element are: charge, epsilon, r_min, epsilon_14, r_min_14,
# read from the existing complex_solvated.ORCAFF.prms written by step13a.
# The site carries charge +1 and the LJ parameters of a sodium ion, which is the closest
# available analogue for a small monovalent cation and is a STARTING POINT, not a
# justified choice.
python3 - <<'PY'
L=open('sub.xyz').read().splitlines(); n=int(L[0].split()[0])
els=[l.split()[0] for l in L[2:2+n]]
with open('site.ORCAFF.prms','w') as fh:
    fh.write("$fftype\nAMBER\n$atoms\n")
    fh.write(f"{n+1} 1 4\n")
    for i,e in enumerate(els, start=1):
        fh.write(f"{i:>6}   {e:<3}  0.000000    -0.000000     0.000000    -0.000000     0.000000\n")
    # Sodium LJ parameters taken verbatim from the committed
    # 05_qmmm/13_bridge/complex_solvated.ORCAFF.prms, which step13a generated from the
    # AMBER topology. Line 6260 there reads:
    #   6256   Na     1.000000    -0.087439     2.738000    -0.043720     2.738000
    # These are the ff14SB/Joung-Cheatham sodium parameters as ORCA received them, not
    # invented values. They are still a PLACEHOLDER in the sense that a designed cation
    # is not necessarily sodium-sized; the point of this test is the MECHANISM.
    fh.write(f"{n+1:>6}   Na   1.000000    -0.087439     2.738000    -0.043720     2.738000\n")
    fh.write("$end\n")
print(f"  wrote site.ORCAFF.prms with {n+1} atoms, last one the +1 site")
PY

"$ORCA/orca" testA.inp > testA.out 2>&1 </dev/null
if grep -q "TERMINATED NORMALLY" testA.out; then
  echo "  PASS: QMMM ran without CPCM"
elif grep -qi "not implemented" testA.out; then
  echo "  FAIL: ORCA still refuses. The message:"
  grep -i -A2 "not implemented" testA.out | head -5 | sed 's/^/    /'
  echo "  -> the LJ route is unavailable regardless of CPCM; surrogates stand."
  echo "end=$(date)"; exit 0
else
  echo "  FAIL for another reason. Last 20 lines:"
  tail -20 testA.out | sed 's/^/    /'
  echo "end=$(date)"; exit 0
fi
echo

# --------------------------------- C: relaxed optimisation, where does the oxygen stop?
echo "=== TEST C: relaxed optimisation with the site carrying LJ"
sed 's/TightSCF QMMM/TightSCF QMMM Opt/' testA.inp > testC.inp
cat >> testC.inp <<'EOF'
%geom
  Constraints
    { C 24 C }      # freeze the charge site; the substrate relaxes around it
  end
end
EOF
"$ORCA/orca" testC.inp > testC.out 2>&1 </dev/null

echo "--- and the bare point charge, for comparison, same geometry and level"
cat > testBare.inp <<'EOF'
! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF Opt
%maxcore 3000
%pal nprocs 8 end
%pointcharges "bare.pc"
* xyzfile -2 1 sub.xyz
EOF
"$ORCA/orca" testBare.inp > testBare.out 2>&1 </dev/null
echo

# ------------------------------------------------------------------- read the result
python3 - <<'PY'
import re, math
px,py,pz = [float(x) for x in open('probe.xyz').read().split()]

def final_geom(path, nsub=24):
    t=open(path, errors='replace').read()
    if "TERMINATED NORMALLY" not in t: return None, "did not terminate"
    blocks=re.findall(r"CARTESIAN COORDINATES \(ANGSTROEM\)\n-+\n(.*?)\n\n", t, re.S)
    if not blocks: return None, "no coordinate block"
    at=[]
    for line in blocks[-1].splitlines():
        f=line.split()
        if len(f)>=4:
            try: at.append((f[0], float(f[1]), float(f[2]), float(f[3])))
            except ValueError: pass
    return at[:nsub], None

print(f"{'run':<22}{'closest O to site':>20}{'verdict':>34}")
for label, path in (("bare point charge","testBare.out"), ("site with LJ","testC.out")):
    at, err = final_geom(path)
    if at is None:
        print(f"{label:<22}{'-':>20}   {err}")
        continue
    ds=[(math.sqrt((x-px)**2+(y-py)**2+(z-pz)**2), e)
        for e,x,y,z in at if e.upper().startswith('O')]
    d=min(ds)[0]
    if   d < 2.0: v="TOO WEAK - no better than bare"
    elif d <= 3.2: v="PHYSICAL - salt-bridge range"
    elif d <= 4.0: v="slightly long"
    else: v="TOO HARD - the charge cannot act"
    print(f"{label:<22}{d:>17.3f} A   {v:>31}")

print()
print("step_7 recorded 0.879 A for a bare +1 with a designed charge set, 0.641 A inside")
print("the oxygen vdW contact. The criterion here is 2.0-3.2 A, the range the structural")
print("consensus gives for Arg-carboxylate salt bridges, because step_7's own reframing")
print("is that inward motion is correct down to that distance and unphysical only below.")
print()
print("If the LJ site stops in range, the mechanism works and the next question is which")
print("LJ parameters are justified - the sodium values used here are a placeholder.")
print("If it does not, molecular surrogates stand as step_7b concluded.")
PY
echo
echo "end=$(date)"
