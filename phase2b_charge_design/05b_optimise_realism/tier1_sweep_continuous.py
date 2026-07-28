#!/usr/bin/env python3
# tier1_sweep_continuous.py -- CONTINUOUS fractional charges q_i in [-1,+1] (GOCAT charge model).
# (New script; discrete tier1_sweep_preorg_surrsized.py retained unchanged.)
#
# WHY: GOCAT (Dittner 2019, line 1342) uses partial charges q_i in [-1,+1] (CONTINUOUS), min-distance, and
# (they use) net-neutrality; working designs = 81 such charges on a sphere. Our discrete +-1 model either
# under-constrains (few charges) or OVER-POLARISES (K10 fragmented) - because full +-1 charges near the -2
# dianion are a harsh field. This variant relaxes to CONTINUOUS magnitudes (the actual GOCAT model), which
# can shape a GENTLE distributed field. Tests Path C: can a GOCAT-scale fractional design hold the
# near-attack geometry WITHOUT the cradle?
#
# NET-FREE (user): NO net-neutrality constraint - charges free to sum to anything (unlike GOCAT's net-0).
# On OUR dv_grid (not a sphere). Certificate PRESERVED: Dv objective is LINEAR in q_i -> this is an LP/MILP
# (continuous q + binary use-indicators for the K-count and min-distance), HiGHS solves to gap 0.000.
#
# Formulation:
#   per site i:  q_i in [-1,1] (continuous) ; u_i in {0,1} (site used)
#   link:        -u_i <= q_i <= u_i         (q_i = 0 unless the site is used)
#   count:       sum_i u_i = K              (exactly K active charges)
#   min-dist:    u_i + u_j <= 1 for sites closer than min_dist
#   objective:   minimize  sum_i q_i * Dv_i * HA2KCAL   (max-lowering; Dv linear in q)  [+ optional mu*|D|]
# Net-free: no constraint on sum_i q_i.

import sys, argparse
import highspy
HA2KCAL=627.509

def load_grid(path):
    g=[]
    for ln in open(path).read().splitlines()[1:]:
        p=ln.split("\t"); g.append((int(p[0]),float(p[1]),float(p[2]),float(p[3]),float(p[5]),p[4]))
        # idx,x,y,z,Dv(col5),shell(col4)  -- match dv_grid.tsv columns: idx x y z shell V_R V_TS Dv
    return g

def load_grid_dv(path):
    # dv_grid.tsv header: idx x_ang y_ang z_ang shell V_R V_TS Dv
    g=[]
    for ln in open(path).read().splitlines()[1:]:
        p=ln.split("\t")
        g.append((int(p[0]),float(p[1]),float(p[2]),float(p[3]),p[4],float(p[7])))  # idx,x,y,z,shell,Dv
    return g

def solve(grid, K, min_dist=4.7, tmax=300):
    n=len(grid)
    h=highspy.Highs(); h.setOptionValue("output_flag",False)
    h.setOptionValue("time_limit",float(tmax)); h.setOptionValue("mip_rel_gap",0.0)
    q=[h.addVariable(lb=-1.0, ub=1.0) for _ in range(n)]   # continuous fractional charge
    u=[h.addBinary() for _ in range(n)]                    # site-used indicator
    for i in range(n):
        h.addConstr(q[i] - u[i] <= 0)     # q_i <= u_i
        h.addConstr(q[i] + u[i] >= 0)     # q_i >= -u_i   => q_i=0 unless u_i=1
    h.addConstr(sum(u[i] for i in range(n)) == K)
    for i in range(n):
        xi,yi,zi=grid[i][1:4]
        for j in range(i+1,n):
            xj,yj,zj=grid[j][1:4]
            if (xi-xj)**2+(yi-yj)**2+(zi-zj)**2 < min_dist*min_dist:
                h.addConstr(u[i]+u[j] <= 1)
    # objective: max-lowering, Dv (Hartree) * q, to kcal. minimize (most negative = most lowering)
    barrier=sum(q[i]*grid[i][5]*HA2KCAL for i in range(n))
    h.minimize(barrier)
    h.run()
    gap=h.getInfo().mip_gap
    placed=[]
    for i in range(n):
        qi=h.variableValue(q[i])
        if abs(qi)>1e-4:
            placed.append({"idx":grid[i][0],"q":qi,"shell":grid[i][4],"dv":grid[i][5],"pos":grid[i][1:4]})
    ddE=sum(p["q"]*p["dv"]*HA2KCAL for p in placed)
    return ddE,gap,placed

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("grid_tsv"); ap.add_argument("--bare",type=float,required=True)
    ap.add_argument("--K",type=int,required=True)
    ap.add_argument("--min_dist",type=float,default=4.7)
    a=ap.parse_args()
    grid=load_grid_dv(a.grid_tsv)
    print("grid=%d pts  K=%d  bare=+%.2f  min_dist=%.1f A  CONTINUOUS q in[-1,1] NET-FREE (HiGHS certified)"%(
        len(grid),a.K,a.bare,a.min_dist))
    ddE,gap,pl=solve(grid,a.K,a.min_dist)
    net=sum(p["q"] for p in pl)
    print("  proxy barrier = %+.2f  (ddE %+.2f)  gap %.3f  net charge %+.3f  (%d active sites)"%(
        a.bare+ddE,ddE,gap,net,len(pl)))
    for p in sorted(pl,key=lambda z:z["idx"]):
        print("    site %3d  q=%+.3f  shell %s  Dv=%+.5f"%(p["idx"],p["q"],p["shell"],p["dv"]))

if __name__=="__main__":
    main()
