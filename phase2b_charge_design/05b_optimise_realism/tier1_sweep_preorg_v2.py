#!/usr/bin/env python3
# 05b Step 2 (PLAN.md): Tier-1 certified screen + PREORGANISATION PENALTY (penalty-sweep).
#
# Extends 05_optimise/tier1_sweep.py with ONE new idea: price in reactant preorganisation so the
# certified optimum is pushed off the C1-C6-distorting sites the frozen Dv proxy was blind to.
#
# NET-FREE (unchanged): at each K the MILP picks whatever +-1 arrangement is optimal - NO net-charge
# constraint. Realism ingredients already present/baked in: bounded +-1 (menu), excluded-volume
# (min_dist no-overlap), CPCM(eps=4) screening (already in dv_grid V_R/V_TS/Dv).
#
# NEW: distorting drive of a design  D = sum_i (sum_k V[k]*y[i][k]) * a_pot_i   [LINEAR in placements]
#   a_pot_i = V_i(C1) - V_i(C6) from Step 0 (preorg/distorting_coeffs.tsv, per unit +1 at site i).
#   |D| linearised with one aux var Dabs>=0 and two rows (Dabs>=D, Dabs>=-D) - same epigraph trick
#   the distributed objective already uses, so the problem stays a certifiable MILP (HiGHS gap 0.000).
# Objective (max-lowering + penalty):  minimize  sum_i sum_k V[k]*Dv_i*y[i][k]  +  mu * Dabs
#   mu in kcal/mol per (charge*Angstrom^-1) unit of distortion; swept. mu=0 reproduces plain max-lowering.
#
# Certificate note: at each FIXED mu the objective is linear -> HiGHS returns gap=0.000. Sweeping mu
# traces the certified trade-off between barrier-lowering (Dv) and preorganisation (low |D|). Each point
# on that curve is individually certified.

import sys, math, argparse, highspy
HA2KCAL=627.509

def load_dv_grid(path):
    grid=[]
    with open(path) as fh:
        next(fh)
        for ln in fh:
            p=ln.rstrip("\n").split("\t")
            # idx,x,y,z,shell,V_R,V_TS,Dv  -> keep idx,x,y,z,Dv,shell
            grid.append((int(p[0]),float(p[1]),float(p[2]),float(p[3]),float(p[7]),p[4]))
    return grid

def load_apot(path, col="a_pot"):
    # preorg/distorting_coeffs.tsv : idx, shell, a_pot, a_pot_norm, a_force_xcheck
    a={}
    with open(path) as fh:
        next(fh)
        for ln in fh:
            p=ln.rstrip("\n").split("\t")
            a[int(p[0])]=float(p[2]) if col=="a_pot" else float(p[4])   # col2=a_pot, col4=a_force
    return a

def ether_o(xyz, idx=7):
    L=open(xyz).read().splitlines(); return [float(t) for t in L[2+idx].split()[1:4]]

def solve(grid, apot, K, menu, mu, min_dist=2.5, tmax=120):
    n=len(grid); V=list(menu); nv=len(V)
    h=highspy.Highs(); h.setOptionValue("output_flag",False)
    h.setOptionValue("time_limit",float(tmax)); h.setOptionValue("mip_rel_gap",0.0)
    y=[[h.addBinary() for _ in range(nv)] for _ in range(n)]
    for i in range(n): h.addConstr(sum(y[i][k] for k in range(nv))<=1)
    h.addConstr(sum(y[i][k] for i in range(n) for k in range(nv))==K)
    # excluded-volume no-overlap (unchanged)
    for i in range(n):
        xi,yi,zi=grid[i][1:4]
        for j in range(i+1,n):
            xj,yj,zj=grid[j][1:4]
            if (xi-xj)**2+(yi-yj)**2+(zi-zj)**2<min_dist*min_dist:
                h.addConstr(sum(y[i][k] for k in range(nv))+sum(y[j][k] for k in range(nv))<=1)
    # distorting drive D = sum_i (sum_k V[k]*y[i][k]) * a_pot_i   (linear); |D| via aux
    Dexpr=sum(V[k]*apot[grid[i][0]]*y[i][k] for i in range(n) for k in range(nv))
    Dabs=h.addVariable(lb=0)
    h.addConstr(Dabs-Dexpr>=0); h.addConstr(Dabs+Dexpr>=0)
    # objective: max-lowering (Dv, in Ha) scaled to kcal + mu*|D|
    barrier=sum(V[k]*grid[i][4]*HA2KCAL*y[i][k] for i in range(n) for k in range(nv))
    h.minimize(barrier + mu*Dabs)
    # extract
    pl=[]
    for i in range(n):
        for k in range(nv):
            if h.variableValue(y[i][k])>0.5:
                pl.append({"q":V[k],"shell":grid[i][5],"dv":grid[i][4],"pos":grid[i][1:4],
                           "apot":apot[grid[i][0]],"idx":grid[i][0]})
    ddE=sum(p["q"]*p["dv"]*HA2KCAL for p in pl)
    Dval=sum(p["q"]*p["apot"] for p in pl)
    return ddE, Dval, h.getInfo().mip_gap, pl

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("dv_grid"); ap.add_argument("reactant_xyz"); ap.add_argument("apot_tsv")
    ap.add_argument("--bare",type=float,required=True)
    ap.add_argument("--menu",default="-1,1")
    ap.add_argument("--K",type=int,default=2)
    ap.add_argument("--mu",default="0,20,50,100,200,500")  # kcal/mol per unit |D|
    ap.add_argument("--coeff",default="a_pot",choices=["a_pot","a_force"])
    a=ap.parse_args()
    grid=load_dv_grid(a.dv_grid); apot=load_apot(a.apot_tsv, a.coeff); o3=ether_o(a.reactant_xyz)
    menu=[float(x) for x in a.menu.split(",")]; mus=[float(x) for x in a.mu.split(",")]
    print("grid=%d pts  menu=%s  K=%d  bare=+%.2f  (net-free, HiGHS certified)"%(len(grid),menu,a.K,a.bare))
    print("preorg penalty sweep: distorting drive D = sum q*a_COEFF (a_pot or a_force per --coeff)")
    print("  (D>0 = design drives C1-C6 apart = distorting; goal: mu pushes certified optimum toward D<=0)\n")
    print(" mu    | ddE(Dv)  barrier |    D(distort)  gap  | placed sites (idx:q@shell, apot)")
    print("-------+------------------+---------------------+-------------------------------------")
    for mu in mus:
        ddE,Dval,gap,pl=solve(grid,apot,a.K,menu,mu)
        sites=" ".join("%d:%+g@%s(a=%+.3f)"%(p["idx"],p["q"],p["shell"],p["apot"]) for p in sorted(pl,key=lambda z:z["idx"]))
        print("%6.0f | %+7.2f  %+7.2f | %+12.4e %.3f | %s"%(mu,ddE,a.bare+ddE,Dval,gap,sites))
