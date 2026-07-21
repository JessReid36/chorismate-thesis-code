#!/usr/bin/env python3
# Export certified Tier-1 designs (K list, chosen objective) to ORCA external point-charge (.pc)
# files for Tier-2. .pc format: line1=count, then "q x y z" (charge a.u., position Angstrom, in the
# SAME box frame as the substrate geometry 01_geometry/reactant.xyz). Plus a manifest.
import argparse, math
from tier1_sweep import load_dv_grid, ether_o, maxlower, distributed, derive_band
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("dv_grid"); ap.add_argument("reactant_xyz")
    ap.add_argument("--bare",type=float,default=17.47); ap.add_argument("--menu",default="-1,1")
    ap.add_argument("--objective",choices=["maxlower","distributed"],default="maxlower")
    ap.add_argument("--K",default="1,2,3,4"); ap.add_argument("--outdir",default=".")
    a=ap.parse_args()
    grid=load_dv_grid(a.dv_grid); o3=ether_o(a.reactant_xyz)
    menu=[float(x) for x in a.menu.split(",")]; band=derive_band(grid,a.bare)
    man=open(a.outdir+"/designs_manifest.tsv","w")
    man.write("design\tobjective\tK\tddE_kcal\tbarrier_kcal\tgap\tn\tcharges(q@dist_to_etherO_A)\n")
    for K in [int(x) for x in a.K.split(",")]:
        if a.objective=="maxlower": ddE,gap,pl=maxlower(grid,K,menu)
        else: ddE,gap,pl=distributed(grid,K,menu,band)
        name="design_%s_K%d"%(a.objective,K)
        with open(a.outdir+"/"+name+".pc","w") as f:
            f.write("%d\n"%len(pl))
            for p in pl:
                x,y,z=p["pos"]; f.write("%+.4f  %.6f %.6f %.6f\n"%(p["q"],x,y,z))
        desc=";".join("%+g@%.2f"%(p["q"],math.sqrt(sum((p["pos"][k]-o3[k])**2 for k in range(3)))) for p in sorted(pl,key=lambda z:z["dv"]))
        man.write("%s\t%s\t%d\t%+.2f\t%+.2f\t%.4f\t%d\t%s\n"%(name,a.objective,K,ddE,a.bare+ddE,gap,len(pl),desc))
        print("%s: %d charges  ddE=%+.2f  barrier=%+.2f  gap=%.4f -> %s.pc"%(name,len(pl),ddE,a.bare+ddE,gap,name))
    man.close(); print("wrote designs_manifest.tsv")
if __name__=="__main__": main()
