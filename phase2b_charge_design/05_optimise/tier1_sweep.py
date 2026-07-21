#!/usr/bin/env python3
# Tier-1 certified screen (HiGHS backend). Two objectives, K-sweep, net-free, configurable charge
# menu (any magnitudes, both signs), any shell mixed freely. HiGHS closes the branch-and-bound gap
# far faster than CBC (certified 0.000 where CBC stalled). Provable-gap = the novelty vs GOCAT's GA.
import sys, math, argparse, highspy
HA2KCAL=627.509
def load_dv_grid(path):
    grid=[]
    with open(path) as fh:
        next(fh)
        for ln in fh:
            p=ln.rstrip("\n").split("\t")
            grid.append((int(p[0]),float(p[1]),float(p[2]),float(p[3]),float(p[7]),p[4]))
    return grid
def ether_o(xyz, idx=7):
    L=open(xyz).read().splitlines(); return [float(t) for t in L[2+idx].split()[1:4]]
def _skel(grid,K,menu,min_dist,tmax):
    n=len(grid); V=list(menu); nv=len(V)
    h=highspy.Highs(); h.setOptionValue("output_flag",False)
    h.setOptionValue("time_limit",float(tmax)); h.setOptionValue("mip_rel_gap",0.0)
    y=[[h.addBinary() for _ in range(nv)] for _ in range(n)]
    for i in range(n): h.addConstr(sum(y[i][k] for k in range(nv))<=1)
    h.addConstr(sum(y[i][k] for i in range(n) for k in range(nv))==K)
    for i in range(n):
        xi,yi,zi=grid[i][1:4]
        for j in range(i+1,n):
            xj,yj,zj=grid[j][1:4]
            if (xi-xj)**2+(yi-yj)**2+(zi-zj)**2<min_dist*min_dist:
                h.addConstr(sum(y[i][k] for k in range(nv))+sum(y[j][k] for k in range(nv))<=1)
    return h,y,V,nv,n
def _place(grid,h,y,V,nv,n):
    out=[]
    for i in range(n):
        for k in range(nv):
            if h.variableValue(y[i][k])>0.5:
                out.append({"q":V[k],"shell":grid[i][5],"dv":grid[i][4],"pos":grid[i][1:4]})
    return out
def maxlower(grid,K,menu,min_dist=2.5,tmax=120):
    h,y,V,nv,n=_skel(grid,K,menu,min_dist,tmax)
    h.minimize(sum(V[k]*grid[i][4]*y[i][k] for i in range(n) for k in range(nv)))
    return h.getObjectiveValue()*HA2KCAL, h.getInfo().mip_gap, _place(grid,h,y,V,nv,n)
def distributed(grid,K,menu,band,min_dist=2.5,tmax=120):
    h,y,V,nv,n=_skel(grid,K,menu,min_dist,tmax)
    Mx=h.addVariable(lb=0)
    contribs=[sum(V[k]*grid[i][4]*HA2KCAL*y[i][k] for k in range(nv)) for i in range(n)]
    for c in contribs: h.addConstr(Mx-c>=0); h.addConstr(Mx+c>=0)
    tot=sum(contribs); h.addConstr(tot>=band[0]); h.addConstr(tot<=band[1])
    h.minimize(Mx)
    return h.getObjectiveValue(), h.getInfo().mip_gap, _place(grid,h,y,V,nv,n)
def derive_band(grid,bare,floor=5.0):
    single_best=min(g[4] for g in grid)*HA2KCAL; band_hi=single_best; band_lo=-(bare-floor)
    if band_lo>=band_hi: band_lo=band_hi-2.0
    return (band_lo,band_hi)
def show(grid,o3,pl):
    for p in sorted(pl,key=lambda z:z["dv"]):
        d=math.sqrt(sum((p["pos"][k]-o3[k])**2 for k in range(3)))
        tag="<- +1 @ ether-O valley (Arg90)" if (p["q"]>0 and d<4.5 and p["dv"]<0) else ""
        print("      q=%+g shell %s dv=%+.5f  %.2f A from ether-O %s"%(p["q"],p["shell"],p["dv"],d,tag))
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("dv_grid"); ap.add_argument("reactant_xyz")
    ap.add_argument("--bare",type=float,required=True)
    ap.add_argument("--menu",default="-1,1"); ap.add_argument("--Kmax",type=int,default=10)
    ap.add_argument("--place-K",default="1,2,3")
    a=ap.parse_args()
    grid=load_dv_grid(a.dv_grid); o3=ether_o(a.reactant_xyz)
    menu=[float(x) for x in a.menu.split(",")]; Ks=list(range(1,a.Kmax+1))
    band=derive_band(grid,a.bare)
    print("grid=%d pts, ether-O=%s, menu=%s, bare barrier=+%.2f (solver: HiGHS)"%(len(grid),[round(c,2) for c in o3],menu,a.bare))
    print("distributed band (total ddE): [%.2f, %.2f] kcal/mol\n"%band)
    print("K | MAX-LOWERING            | DISTRIBUTED (min max|contrib|)")
    print("  | ddE      barrier   gap  | max|contrib|  total    barrier   gap")
    designs={}
    for K in Ks:
        dl,gl,pl=maxlower(grid,K,menu)
        dd,gd,pd=distributed(grid,K,menu,band)
        totd=sum(p["q"]*p["dv"]*HA2KCAL for p in pd)
        print("%2d| %+7.2f  %+7.2f  %.3f | %10.2f  %+7.2f  %+7.2f  %.3f"%(K,dl,a.bare+dl,gl,dd,totd,a.bare+totd,gd))
        designs[K]=(pl,pd)
    for K in [int(x) for x in a.place_K.split(",") if int(x) in designs]:
        print("\n--- K=%d max-lowering ---"%K); show(grid,o3,designs[K][0])
        print("--- K=%d distributed ---"%K); show(grid,o3,designs[K][1])
