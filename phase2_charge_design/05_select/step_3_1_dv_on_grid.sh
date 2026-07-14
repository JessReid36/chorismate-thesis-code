#!/usr/bin/env bash
# Phase 2 step 3.1 - evaluate the difference potential Dv = V_TS - V_R on the 326-point candidate
# grid. Reuses the validated 1.2 orca_vpot machinery (gbw + .densities container, points in BOHR,
# container basename as arg 5). Produces the objective vector dv_grid.tsv for the Stage 3 optimizer.
set -euo pipefail
source "$HOME/system_development/phase2_charge_design/00_admin/phase2_paths.sh"

dvpot="$PHASE2_ROOT/03_dvpot"     # densities live here (from 1.1b)
grid="$PHASE2_ROOT/04_grid"       # grid_final.xyz pushed here from the workstation
stage="$PHASE2_ROOT/05_select"
mkdir -p "$stage"; cd "$stage"

ORCA=/home/apps2/ORCA/6.0.1
export PATH="/apps/openmpi/4.1.1/bin:${PATH:-}"
export LD_LIBRARY_PATH="$ORCA/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"

gridxyz="$grid/grid_final.xyz"
[[ -s "$gridxyz" ]] || { echo "FAIL: grid_final.xyz not found at $gridxyz (scp it from the workstation)"; exit 1; }
for f in sp_reactant.gbw sp_reactant.densities sp_ts.gbw sp_ts.densities; do
  [[ -s "$dvpot/$f" ]] || { echo "FAIL: missing $dvpot/$f"; exit 1; }
done

# grid_final.xyz format: line1=N, line2=comment, then "X x y z shell=..". Convert coords (Angstrom)
# to the orca_vpot points file (Bohr): line1=N, then "x y z". Preserve row order (Dv matches by row).
python3 - "$gridxyz" grid_points_bohr.xyz grid_pts_ang.tsv <<'PY'
import sys
from pathlib import Path
ANG2BOHR = 1.8897259886
src, outbohr, outang = map(Path, sys.argv[1:4])
lines = src.read_text().splitlines()
n = int(lines[0].split()[0])
rows = []
for ln in lines[2:2+n]:
    p = ln.split()
    x, y, z = float(p[1]), float(p[2]), float(p[3])
    shell = ln.split("shell=")[-1].strip() if "shell=" in ln else "NA"
    rows.append((x, y, z, shell))
with outbohr.open("w") as fh:
    fh.write(f"{len(rows)}\n")
    for x,y,z,_ in rows:
        fh.write(f"{x*ANG2BOHR:.10f} {y*ANG2BOHR:.10f} {z*ANG2BOHR:.10f}\n")
with outang.open("w") as fh:
    fh.write("idx\tx_ang\ty_ang\tz_ang\tshell\n")
    for i,(x,y,z,s) in enumerate(rows):
        fh.write(f"{i}\t{x:.4f}\t{y:.4f}\t{z:.4f}\t{s}\n")
print(f"converted {len(rows)} grid points to Bohr")
PY

echo "=== orca_vpot on the grid (reactant) ==="
"$ORCA/orca_vpot" "$dvpot/sp_reactant.gbw" sp_reactant.scfp grid_points_bohr.xyz vR_grid.out "$dvpot/sp_reactant"
echo "=== orca_vpot on the grid (TS) ==="
"$ORCA/orca_vpot" "$dvpot/sp_ts.gbw" sp_ts.scfp grid_points_bohr.xyz vTS_grid.out "$dvpot/sp_ts"

echo "=== assemble Dv = V_TS - V_R ==="
python3 - <<'PY'
from pathlib import Path
def read_vpot(f):
    vals=[]
    for ln in Path(f).read_text().splitlines():
        p=ln.split()
        if len(p)>=4:
            try: vals.append(float(p[-1]))
            except ValueError: pass
    return vals
vR, vTS = read_vpot("vR_grid.out"), read_vpot("vTS_grid.out")
meta = Path("grid_pts_ang.tsv").read_text().splitlines()[1:]
if not (len(vR)==len(vTS)==len(meta)):
    print(f"FAIL: length mismatch vR={len(vR)} vTS={len(vTS)} grid={len(meta)}"); raise SystemExit(1)
rows=[]
nneg=0
with open("dv_grid.tsv","w") as fh:
    fh.write("idx\tx_ang\ty_ang\tz_ang\tshell\tV_R\tV_TS\tDv\n")
    for m,a,b in zip(meta,vR,vTS):
        idx,x,y,z,shell = m.split("\t")
        dv=b-a; nneg += (dv<0)
        fh.write(f"{idx}\t{x}\t{y}\t{z}\t{shell}\t{a:.6f}\t{b:.6f}\t{dv:.6f}\n")
        rows.append(dv)
import statistics as st
print(f"assembled Dv for {len(rows)} grid points")
print(f"  Dv range [{min(rows):+.6f}, {max(rows):+.6f}] Eh")
print(f"  Dv mean {st.mean(rows):+.6f}  std {st.pstdev(rows):.6f}")
print(f"  stabilizing points (Dv<0, a +1 cation lowers the barrier): {nneg}/{len(rows)}")
print("saved dv_grid.tsv (objective vector for Stage 3)")
PY
echo "done 3.1"
