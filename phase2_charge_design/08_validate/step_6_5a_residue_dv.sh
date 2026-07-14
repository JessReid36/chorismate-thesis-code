#!/usr/bin/env bash
# Phase 2 step 6.5a - evaluate Dv = V_TS - V_R at the enzyme residue charge-centers (for the
# enzyme-constrained comparison, 6.5). Same validated orca_vpot machinery as 3.1, pointed at
# residue_centers.xyz instead of the grid. Produces residue_dv.npz for the constrained MIQP.
set -euo pipefail
source "$HOME/system_development/phase2_charge_design/00_admin/phase2_paths.sh"
dvpot="$PHASE2_ROOT/03_dvpot"
stage="$PHASE2_ROOT/05_select"
cd "$stage"

ORCA=/home/apps2/ORCA/6.0.1
export PATH="/apps/openmpi/4.1.1/bin:${PATH:-}"
export LD_LIBRARY_PATH="$ORCA/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"

[[ -s residue_centers.xyz ]] || { echo "FAIL: residue_centers.xyz missing"; exit 1; }
for f in sp_reactant.gbw sp_reactant.densities sp_ts.gbw sp_ts.densities; do
  [[ -s "$dvpot/$f" ]] || { echo "FAIL: missing $dvpot/$f"; exit 1; }
done

# residue_centers.xyz: line1=N, then "X x y z" (Angstrom). Convert to orca_vpot points (Bohr).
python3 - residue_centers.xyz residue_points_bohr.xyz <<'PY'
import sys
from pathlib import Path
ANG2BOHR=1.8897259886
src,out=map(Path,sys.argv[1:3])
lines=src.read_text().splitlines()
n=int(lines[0].split()[0])
rows=[l.split()[1:4] for l in lines[1:1+n]]
with out.open("w") as fh:
    fh.write(f"{n}\n")
    for x,y,z in rows:
        fh.write(f"{float(x)*ANG2BOHR:.10f} {float(y)*ANG2BOHR:.10f} {float(z)*ANG2BOHR:.10f}\n")
print(f"converted {n} residue centers to Bohr")
PY

echo "=== orca_vpot at residue centers (reactant) ==="
"$ORCA/orca_vpot" "$dvpot/sp_reactant.gbw" sp_reactant.scfp residue_points_bohr.xyz vR_res.out "$dvpot/sp_reactant"
echo "=== orca_vpot at residue centers (TS) ==="
"$ORCA/orca_vpot" "$dvpot/sp_ts.gbw" sp_ts.scfp residue_points_bohr.xyz vTS_res.out "$dvpot/sp_ts"

echo "=== assemble residue Dv and pack with names/signs ==="
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 - <<'PY'
import numpy as np
from pathlib import Path
def read_vpot(f):
    vals=[]
    for ln in Path(f).read_text().splitlines():
        p=ln.split()
        if len(p)>=4:
            try: vals.append(float(p[-1]))
            except ValueError: pass
    return np.array(vals)
vR=read_vpot("vR_res.out"); vTS=read_vpot("vTS_res.out")
meta=[l.split("\t") for l in Path("residue_meta.tsv").read_text().splitlines()[1:]]
names=np.array([m[1] for m in meta]); signs=np.array([int(m[2]) for m in meta])
centers=np.array([[float(m[3]),float(m[4]),float(m[5])] for m in meta])
assert len(vR)==len(vTS)==len(names), f"len mismatch {len(vR)}/{len(vTS)}/{len(names)}"
dv=vTS-vR
np.savez("residue_dv.npz", names=names, signs=signs, centers=centers, dv=dv)
print(f"{'residue':>8} {'sign':>4} {'Dv(Eh)':>11}  (Dv<0 => a cation there stabilizes the TS)")
for n,s,d in zip(names,signs,dv):
    print(f"{n:>8} {s:+4d} {d:+11.6f}")
print("\nsaved residue_dv.npz")
PY
echo "done 6.5a"
