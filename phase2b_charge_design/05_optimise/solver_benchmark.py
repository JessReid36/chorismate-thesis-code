#!/usr/bin/env python3
# Solver benchmark: CBC (python-mip, the default) vs HiGHS (highspy) on the Tier-1 DISTRIBUTED
# min-max MILP at production scale (~330 candidate points, K=7). Both find the SAME optimum; the
# question is whether the branch-and-bound gap CLOSES (certified) within a time budget. Reproducible.
import math, random, time
random.seed(5); HA2KCAL=627.509
valley=(2.6,0,0); peak=(-2.2,1.5,0); grid=[]; gid=0
for shell,r in [("2.0",3.2),("3.0",4.2),("4.0",5.4)]:
    for _ in range(110):
        th,ph=random.uniform(0,math.pi),random.uniform(0,2*math.pi)
        x=r*math.sin(th)*math.cos(ph); y=r*math.sin(th)*math.sin(ph); z=r*math.cos(th)
        dv=-0.008*math.exp(-math.dist((x,y,z),valley)/2.0)+0.009*math.exp(-math.dist((x,y,z),peak)/2.0)+random.gauss(0,0.0004)
        grid.append((gid,x,y,z,dv,shell)); gid+=1
n=len(grid); menu=[-1.0,1.0]; nv=2; K=7; MIN_DIST=2.5; band=(-12.47,-5.27); TMAX=120
close=[(i,j) for i in range(n) for j in range(i+1,n)
       if (grid[i][1]-grid[j][1])**2+(grid[i][2]-grid[j][2])**2+(grid[i][3]-grid[j][3])**2 < MIN_DIST*MIN_DIST]
print("distributed min-max MILP: grid=%d pts, close-pairs=%d, K=%d, time budget=%ds/solver\n"%(n,len(close),K,TMAX))
from mip import Model, xsum, BINARY, CONTINUOUS, minimize
def cbc():
    m=Model(); m.verbose=0
    y=[[m.add_var(var_type=BINARY) for _ in range(nv)] for _ in range(n)]; Mx=m.add_var(var_type=CONTINUOUS,lb=0)
    for i in range(n): m += xsum(y[i][k] for k in range(nv))<=1
    m += xsum(y[i][k] for i in range(n) for k in range(nv))==K
    for i,j in close: m += xsum(y[i][k] for k in range(nv))+xsum(y[j][k] for k in range(nv))<=1
    contribs=[xsum(menu[k]*grid[i][4]*HA2KCAL*y[i][k] for k in range(nv)) for i in range(n)]
    for c in contribs: m += Mx>=c; m += Mx>=-c
    tot=xsum(contribs); m += tot>=band[0]; m += tot<=band[1]
    m.objective=minimize(Mx); t=time.time(); m.optimize(max_seconds=TMAX); return m.objective_value,m.gap,time.time()-t
import highspy
def highs():
    h=highspy.Highs(); h.setOptionValue("output_flag",False); h.setOptionValue("time_limit",float(TMAX)); h.setOptionValue("mip_rel_gap",0.0)
    y=[[h.addBinary() for _ in range(nv)] for _ in range(n)]; Mx=h.addVariable(lb=0)
    for i in range(n): h.addConstr(sum(y[i][k] for k in range(nv))<=1)
    h.addConstr(sum(y[i][k] for i in range(n) for k in range(nv))==K)
    for i,j in close: h.addConstr(sum(y[i][k] for k in range(nv))+sum(y[j][k] for k in range(nv))<=1)
    contribs=[sum(menu[k]*grid[i][4]*HA2KCAL*y[i][k] for k in range(nv)) for i in range(n)]
    for c in contribs: h.addConstr(Mx-c>=0); h.addConstr(Mx+c>=0)
    tot=sum(contribs); h.addConstr(tot>=band[0]); h.addConstr(tot<=band[1])
    t=time.time(); h.minimize(Mx); return h.getObjectiveValue(),h.getInfo().mip_gap,time.time()-t
o1,g1,t1=cbc(); print("CBC   (python-mip): obj=%.4f  gap=%.4f  time=%.1fs  -> %s"%(o1,g1,t1,"CERTIFIED" if g1<1e-4 else "NOT certified (gap open)"))
o2,g2,t2=highs(); print("HiGHS (highspy)  : obj=%.4f  gap=%.4f  time=%.1fs  -> %s"%(o2,g2,t2,"CERTIFIED" if g2<1e-4 else "NOT certified"))
print("\nsame optimum: %s | HiGHS speedup to certification: >%.0fx"%(abs(o1-o2)<1e-3, t1/max(t2,0.01)))
