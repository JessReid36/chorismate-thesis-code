#!/bin/bash
# Compute V_R and V_TS at the outer-shell points via orca_vpot, using the COMMITTED R and TS wavefunctions
# from 02_singlepoints. Convention (TIER2_KNOWN_ISSUES E3):
#   orca_vpot <arg1=.gbw> <arg2=.scfp member name, relative> <arg3=points(Bohr)> <arg4=out> <arg5=.densities base full path>
# Edit the paths below to your 02_singlepoints outputs, then run on the HPC (module-loaded ORCA).
set -e
SP=~/system_development/phase2b_charge_design/02_singlepoints   # <-- adjust if different
ORCA_VPOT=/home/apps2/ORCA/6.0.1/orca_vpot
# 1. convert outer_shell_points.tsv (Angstrom) -> points in Bohr for orca_vpot
python3 - <<'PY'
ANG2BOHR=1.8897259886
rows=open("../grid_ext/outer_shell_points.tsv").read().splitlines()[1:]
with open("points_bohr.xyz","w") as fh:
    fh.write("%d\n"%len(rows))
    for ln in rows:
        p=ln.split("\t"); fh.write("%.8f %.8f %.8f\n"%(float(p[1])*ANG2BOHR,float(p[2])*ANG2BOHR,float(p[3])*ANG2BOHR))
print("wrote points_bohr.xyz (%d pts)"%len(rows))
PY
# 2. run orca_vpot for R and TS (adjust member/base names to your 02_singlepoints files)
$ORCA_VPOT $SP/sp_reactant.gbw sp_reactant.scfp points_bohr.xyz vpot_R.out $SP/sp_reactant.densities
$ORCA_VPOT $SP/sp_ts.gbw       sp_ts.scfp       points_bohr.xyz vpot_TS.out $SP/sp_ts.densities
# 3. assemble dv_grid_ext.tsv (V_R, V_TS, Dv=V_TS-V_R) merged with existing 331-point grid
python3 - <<'PY'
def readvpot(f):
    v=[]
    for ln in open(f):
        t=ln.split()
        if len(t)>=4:
            try: v.append(float(t[3]))
            except: pass
    return v
VR=readvpot("vpot_R.out"); VT=readvpot("vpot_TS.out")
rows=open("../grid_ext/outer_shell_points.tsv").read().splitlines()[1:]
assert len(VR)==len(VT)==len(rows), "count mismatch: %d/%d/%d"%(len(VR),len(VT),len(rows))
out=["idx\tx_ang\ty_ang\tz_ang\tshell\tV_R\tV_TS\tDv"]
for ln,vr,vt in zip(rows,VR,VT):
    p=ln.split("\t"); out.append("%s\t%s\t%s\t%s\t%s\t%.8f\t%.8f\t%.8f"%(p[0],p[1],p[2],p[3],p[4],vr,vt,vt-vr))
open("dv_grid_ext.tsv","w").write("\n".join(out)+"\n")
print("wrote dv_grid_ext.tsv (%d outer-shell points with Dv)"%len(rows))
PY
echo "Now MERGE with the existing grid: cat ../../../05_optimise/dv_grid.tsv (keep header) + dv_grid_ext.tsv (skip header) -> dv_grid_full.tsv"
