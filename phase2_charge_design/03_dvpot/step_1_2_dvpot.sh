#!/usr/bin/env bash
# Phase 2 step 1.2 - difference-potential engine + sign check.
# Dv = V_TS - V_R at TEST points via orca_vpot; verifies a +1 charge near the breaking
# C4-O3 ether oxygen gives a stabilizing (negative) q*Dv (the Arg90 role).
# This ORCA build stores the SCF density as member <base>.scfp INSIDE the <base>.densities
# container, so orca_vpot is called as: GBW  PName(member)  XYZ(BOHR)  POTout  DensContainerBase.
set -euo pipefail
source "$HOME/system_development/phase2_charge_design/00_admin/phase2_paths.sh"
stage="$PHASE2_ROOT/03_dvpot"
cd "$stage"

ORCA=/home/apps2/ORCA/6.0.1
export PATH="/apps/openmpi/4.1.1/bin:${PATH:-}"
export LD_LIBRARY_PATH="$ORCA/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"

# .scfp is a container member (not a standalone file); check the container + gbw instead.
for f in sp_reactant.gbw sp_reactant.densities sp_ts.gbw sp_ts.densities; do
  [[ -s "$stage/$f" ]] || { echo "FAIL: missing $f (rerun 1.1b)"; exit 1; }
done

python3 - "$SUBSTRATE_DIR/$SUBSTRATE_REACTANT" "$PHASE2_ROOT/00_admin/reacting_atoms.tsv" test_points.xyz <<'PY'
import sys, math
from pathlib import Path
geo, mapping, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
ANG2BOHR = 1.8897259886
lines = geo.read_text().splitlines()
n = int(lines[0]); at = [(p.split()[0], *map(float, p.split()[1:4])) for p in lines[2:2+n]]
I = {"O3":7, "C4":8, "C1":0, "C6":12}
def vec(i): return [at[i][1], at[i][2], at[i][3]]
def sub(a,b): return [a[k]-b[k] for k in range(3)]
def add(a,b): return [a[k]+b[k] for k in range(3)]
def norm(a): return math.sqrt(sum(x*x for x in a))
def scale(a,s): return [x*s for x in a]
def cross(a,b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
O3, C4, C6, C1 = vec(I["O3"]), vec(I["C4"]), vec(I["C6"]), vec(I["C1"])
d = sub(O3, C4); u = scale(d, 1.0/norm(d)); pA = add(O3, scale(u, 2.0))
d2 = sub(C6, C1); u2 = scale(d2, 1.0/norm(d2)); pB = add(C6, scale(u2, 2.0))
mid = scale(add(add(O3,C4), add(C6,C1)), 0.25)
nrm = cross(sub(C4,C6), sub(C1,O3)); nrm = scale(nrm, 3.0/(norm(nrm) or 1.0)); pC = add(mid, nrm)
pts = {"A_etherO_C-O_dir": pA, "B_formingC6_dir": pB, "C_out_of_plane_mid": pC}
with out.open("w") as fh:
    fh.write(f"{len(pts)}\n")
    for name,p in pts.items():
        pb = [c*ANG2BOHR for c in p]
        fh.write(f"{pb[0]:.10f} {pb[1]:.10f} {pb[2]:.10f}\n")
Path("test_points_labels.tsv").write_text(
    "label\tx_ang\ty_ang\tz_ang\n" +
    "".join(f"{name}\t{p[0]:.4f}\t{p[1]:.4f}\t{p[2]:.4f}\n" for name,p in pts.items()))
print(f"wrote {len(pts)} test points (Bohr) to {out}")
PY

echo "=== running orca_vpot (reactant) ==="
"$ORCA/orca_vpot" sp_reactant.gbw sp_reactant.scfp test_points.xyz vR.out sp_reactant
echo "=== running orca_vpot (TS) ==="
"$ORCA/orca_vpot" sp_ts.gbw sp_ts.scfp test_points.xyz vTS.out sp_ts

echo "=== raw vpot output (first lines, for format check) ==="
echo "--- vR.out ---"; head -8 vR.out
echo "--- vTS.out ---"; head -8 vTS.out

echo "=== Dv = V_TS - V_R + sign check ==="
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
vR, vTS = read_vpot("vR.out"), read_vpot("vTS.out")
labels=[l.split("\t")[0] for l in Path("test_points_labels.tsv").read_text().splitlines()[1:]]
if len(vR)!=len(labels) or len(vTS)!=len(labels):
    print(f"WARN: parsed {len(vR)}/{len(vTS)} potentials but {len(labels)} points - check raw format above")
print(f"{'label':22s} {'V_R':>12s} {'V_TS':>12s} {'Dv=V_TS-V_R':>14s}  q*Dv(+1) sign")
etherO_dv=None
for lab,a,b in zip(labels,vR,vTS):
    dv=b-a
    if lab.startswith("A_etherO"): etherO_dv=dv
    sign = "STABILIZING (-)" if dv<0 else "destabilizing (+)"
    print(f"{lab:22s} {a:12.6f} {b:12.6f} {dv:14.6f}   {sign}")
print()
if etherO_dv is None:
    print("SIGN CHECK INCONCLUSIVE: ether-O point missing")
elif etherO_dv < 0:
    print(f"SIGN CHECK PASS: +1 charge at ether O gives q*Dv = {etherO_dv:.6f} Eh < 0 (stabilizing).")
    print("Convention correct: cation near the developing negative charge lowers the barrier (Arg90 role).")
else:
    print(f"SIGN CHECK FAIL: +1 charge at ether O gives q*Dv = {etherO_dv:.6f} Eh > 0.")
    print("Investigate sign/units before proceeding.")
PY
echo "done 1.2"
