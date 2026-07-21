#!/usr/bin/env python3
# The three certified charge-placement formulations, each a linear MILP solved to a
# provable branch-and-bound gap (the novelty vs GOCAT's genetic algorithm). q_i in {-1,0,+1}
# via two binaries p_i (+1) and q_i (-1); exactly K charges; min-distance MIN_DIST between
# placed charges. All operate on Dv values only -> geometry-independent.
import math
from mip import Model, xsum, BINARY, CONTINUOUS, minimize
HA2KCAL=627.509; MIN_DIST=2.5

def _base(pts, K):
    m=Model(); m.verbose=0
    p=[m.add_var(var_type=BINARY) for _ in pts]; q=[m.add_var(var_type=BINARY) for _ in pts]
    for i in range(len(pts)): m += p[i]+q[i] <= 1
    m += xsum(p[i]+q[i] for i in range(len(pts))) == K
    for i in range(len(pts)):
        for j in range(i+1, len(pts)):
            if math.dist(pts[i][1:4], pts[j][1:4]) < MIN_DIST:
                m += p[i]+q[i]+p[j]+q[j] <= 1
    return m, p, q

def certified_milp(grid, K):
    # maximise barrier lowering: minimise ddE = sum (p_i - q_i)*Dv_i. grid=(id,x,y,z,dv)
    m,p,q=_base(grid,K)
    m.objective=minimize(xsum((p[i]-q[i])*grid[i][4] for i in range(len(grid))))
    m.optimize(max_seconds=120)
    return m.objective_value*HA2KCAL, m.gap

def distributed_minmax(grid, K, band):
    # minimise the largest single-charge contribution, subject to total effect in a band
    # [band_lo, band_hi] (kcal/mol, both negative) that keeps the barrier positive by construction.
    m,p,q=_base(grid,K); Mx=m.add_var(var_type=CONTINUOUS, lb=0)
    contrib=[(p[i]-q[i])*grid[i][4]*HA2KCAL for i in range(len(grid))]
    for c in contrib:
        m += Mx >= c; m += Mx >= -c
    tot=xsum(contrib); m += tot >= band[0]; m += tot <= band[1]
    m.objective=minimize(Mx); m.optimize(max_seconds=120)
    return m.objective_value, m.gap

def wholepath_minmax(pts, b, dv, K):
    # minimise the MAX barrier across NEB images: min t s.t. t >= b_m + sum(p_i-q_i)*Dv_m(i).
    m,p,q=_base(pts,K); t=m.add_var(var_type=CONTINUOUS, lb=-1e6)
    for mm in range(len(b)):
        m += t >= b[mm] + xsum((p[i]-q[i])*dv[mm][i]*HA2KCAL for i in range(len(pts)))
    m.objective=minimize(t); m.optimize(max_seconds=120)
    return m.objective_value, m.gap

def true_path_barrier_from_singlepoint(pts, b, dv, K):
    # single-point optimises only the TS image; report its TS value AND the true max-over-images.
    ts=b.index(max(b)); m,p,q=_base(pts,K)
    m.objective=minimize(b[ts]+xsum((p[i]-q[i])*dv[ts][i]*HA2KCAL for i in range(len(pts))))
    m.optimize(max_seconds=120)
    sel=[(1 if p[i].x>0.5 else 0)-(1 if q[i].x>0.5 else 0) for i in range(len(pts))]
    tsval=b[ts]+sum(sel[i]*dv[ts][i]*HA2KCAL for i in range(len(pts)))
    truemax=max(b[mm]+sum(sel[i]*dv[mm][i]*HA2KCAL for i in range(len(pts))) for mm in range(len(b)))
    return tsval, truemax
