#!/usr/bin/env python3
# K=10 discrete +-1 (net-free) re-optimised on a SINGLE shell. Certified MILP (HiGHS gap 0).
# Usage: optimise_shell.py <dv_grid_full.tsv> <shell e.g. 5.0> --K 10 [--min_dist 4.7]
# Writes <shell>_coords.xyz + <shell>_charges.txt (24 QM substrate + guan/formate surrogates).
import sys, argparse, numpy as np
import highspy
HA2KCAL=627.509
def load_grid(path, shell):
    g=[]
    for ln in open(path).read().splitlines()[1:]:
        p=ln.split("\t")
        if p[4]==shell or ("%.1f"%float(p[4]))==("%.1f"%float(shell)):
            g.append((int(p[0]),float(p[1]),float(p[2]),float(p[3]),float(p[7])))  # idx,x,y,z,Dv
    return g
def opt(grid,K,min_dist):
    n=len(grid); h=highspy.Highs(); h.setOptionValue("output_flag",False); h.setOptionValue("mip_rel_gap",0.0)
    y=[[h.addBinary() for _ in range(2)] for _ in range(n)]; V=[-1.0,1.0]
    for i in range(n): h.addConstr(sum(y[i][k] for k in range(2))<=1)
    h.addConstr(sum(y[i][k] for i in range(n) for k in range(2))==K)
    for i in range(n):
        xi,yi,zi=grid[i][1:4]
        for j in range(i+1,n):
            xj,yj,zj=grid[j][1:4]
            if (xi-xj)**2+(yi-yj)**2+(zi-zj)**2<min_dist*min_dist:
                h.addConstr(sum(y[i][k] for k in range(2))+sum(y[j][k] for k in range(2))<=1)
    h.minimize(sum(V[k]*grid[i][4]*HA2KCAL*y[i][k] for i in range(n) for k in range(2))); h.run()
    ch=[(grid[i][0],V[k],grid[i][1:4]) for i in range(n) for k in range(2) if h.variableValue(y[i][k])>0.5]
    b=sum(V[k]*grid[i][4]*HA2KCAL*h.variableValue(y[i][k]) for i in range(n) for k in range(2))
    return ch,b,h.getInfo().mip_gap
GUAN=np.array([[57.4994,32.313,60.9995],[57.649977,32.69822,59.735446],[57.33946,33.360722,60.4317],[57.846146,31.743182,59.999113],[57.090501,33.185404,61.916349],[57.28667,32.230366,62.180016],[57.204849,33.47794,60.956428],[57.757721,31.055376,61.346705],[57.87207,31.347911,60.386784],[57.447204,31.717878,62.042959]]); GT=GUAN-GUAN[0]; GEL=["C","N","H","H","N","H","H","N","H","H"]; GQ=[0.64,-0.80,0.46,0.46,-0.80,0.46,0.46,-0.80,0.46,0.46]
FORM=np.array([[63.0044,29.1163,52.2256],[64.188225,29.178076,51.829063],[62.293116,28.088917,52.258159],[62.546421,30.05215,52.578361]]); FT=FORM-FORM[0]; FEL=["C","O","O","H"]; FQ=[0.45,-0.80,-0.80,0.15]
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("grid"); ap.add_argument("shell"); ap.add_argument("reactant")
    ap.add_argument("--K",type=int,default=10); ap.add_argument("--min_dist",type=float,default=4.7); ap.add_argument("--bare",type=float,default=17.47)
    a=ap.parse_args()
    grid=load_grid(a.grid,a.shell)
    if len(grid)<a.K: print("shell %s: only %d sites, K=%d infeasible"%(a.shell,len(grid),a.K)); return
    ch,b,gap=opt(grid,a.K,a.min_dist)
    tag="shell%s"%a.shell.replace(".","p")
    L=open(a.reactant).read().splitlines(); qm=L[2:26]; lines=list(qm); charges=[0.0]*24
    for idx,q,pos in ch:
        c=np.array(pos)
        if q>0:
            for el,off in zip(GEL,GT): lines.append("%s  %12.6f %12.6f %12.6f"%(el,*(c+off)))
            charges+=GQ
        else:
            for el,off in zip(FEL,FT): lines.append("%s  %12.6f %12.6f %12.6f"%(el,*(c+off)))
            charges+=FQ
    open("%s_coords.xyz"%tag,"w").write("%d\n%s K=%d +-1 net-free proxy%+.2f\n%s\n"%(len(lines),tag,a.K,a.bare+b,"\n".join(lines)))
    open("%s_charges.txt"%tag,"w").write(",".join("%.3f"%c for c in charges))
    npos=sum(1 for _,q,_ in ch if q>0)
    print("%s: K=%d (%d+/%d-) proxy barrier %+.2f gap %.3f net %+d -> %s_coords.xyz"%(tag,a.K,npos,a.K-npos,a.bare+b,gap,2*npos-a.K,tag))
if __name__=="__main__": main()
