#!/usr/bin/env python3
# Phase 2 step 3.4 - DISCRETE MIQP design (option c): the buildable K-charge catalyst.
# Charges in {-1,0,+1}; at most K active; optionally net-neutral; PHYSICAL minimum-separation
# between any two active charges (default 2.5 A - the floor below which real charged moieties,
# ~2.5-3 A across, would overlap; a genuine salt bridge sits ~2.8-4.0 A). Objective LINEAR =>
# mixed-integer LINEAR program; CBC solves to PROVABLE optimality. Net-neutral & net-free
# variants; compared to the continuous L1-budgeted ceiling (3.2) at matching charge Q=K.
import sys, numpy as np
from mip import Model, xsum, minimize, BINARY, OptimizationStatus
from scipy.spatial import cKDTree
HART2KCAL = 627.509

def load_dv(path):
    idx, xyz, shell, dv = [], [], [], []
    for ln in open(path).read().splitlines()[1:]:
        p = ln.split("\t")
        idx.append(int(p[0])); xyz.append([float(p[1]),float(p[2]),float(p[3])])
        shell.append(float(p[4])); dv.append(float(p[7]))
    return np.array(idx), np.array(xyz), np.array(shell), np.array(dv)

def close_pairs(xyz, sep):
    return cKDTree(xyz).query_pairs(sep)

def solve_miqp(dv, xyz, K, net_neutral, sep, max_seconds=180):
    n = len(dv); md = Model(sense=minimize); md.verbose = 0
    p = [md.add_var(var_type=BINARY) for _ in range(n)]
    m = [md.add_var(var_type=BINARY) for _ in range(n)]
    a = [md.add_var(var_type=BINARY) for _ in range(n)]
    for i in range(n):
        md += p[i] + m[i] <= 1
        md += a[i] == p[i] + m[i]
    md += xsum(a) <= K
    if net_neutral:
        md += xsum(p[i] - m[i] for i in range(n)) == 0
    for (i,j) in close_pairs(xyz, sep):
        md += a[i] + a[j] <= 1
    md.objective = xsum(dv[i]*(p[i] - m[i]) for i in range(n))
    st = md.optimize(max_seconds=max_seconds)
    q = np.array([(1 if p[i].x>0.5 else 0) - (1 if m[i].x>0.5 else 0) for i in range(n)])
    gap = md.gap if md.gap is not None else float('nan')
    return q, md.objective_value, st, gap

def solve_continuous_at_Q(dv, Q, qmax=1.0):
    import cvxpy as cp
    n=len(dv); q=cp.Variable(n)
    cp.Problem(cp.Minimize(dv@q),[cp.norm1(q)<=Q, cp.sum(q)==0, q<=qmax, q>=-qmax]).solve()
    return float(dv@q.value)

def min_active_sep(xyz, q):
    sel = np.where(q!=0)[0]
    if len(sel)<2: return float('inf')
    d,_ = cKDTree(xyz[sel]).query(xyz[sel], k=2)
    return float(d[:,1].min())

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "dv_grid.tsv"
    sep  = float(sys.argv[2]) if len(sys.argv) > 2 else 2.5
    Ks   = [int(x) for x in (sys.argv[3:] or [2,4,6,8])]
    idx, xyz, shell, dv = load_dv(path)
    npairs = len(close_pairs(xyz, sep))
    print(f"loaded {len(dv)} grid points | Dv range [{dv.min():+.6f},{dv.max():+.6f}] Eh, {(dv<0).sum()} stabilizing")
    print(f"separation constraint: {sep} A  ({npairs} close pairs forbidden from co-selection)\n")
    print(f"{'K':>3} {'mode':>12} {'dDE(kcal)':>10} {'status':>9} {'gap':>6} {'(+/-)':>7} "
          f"{'minSep':>7} {'ceiling':>9} {'%ceil':>7}")
    best = {}
    for K in Ks:
        cont = solve_continuous_at_Q(dv, float(K)) * HART2KCAL
        for nn in (True, False):
            if nn and K % 2 != 0: continue
            q, obj, st, gap = solve_miqp(dv, xyz, K, nn, sep)
            kcal=obj*HART2KCAL; npos=int((q>0).sum()); nneg=int((q<0).sum())
            msep=min_active_sep(xyz,q); pct=100*kcal/cont if cont else float('nan')
            sts="OPT" if st==OptimizationStatus.OPTIMAL else str(st).split('.')[-1][:6]
            mode="neutral" if nn else "free"
            print(f"{K:3d} {mode:>12} {kcal:+10.3f} {sts:>9} {gap:6.3f} {npos:3d}/{nneg:<3d} "
                  f"{msep:7.2f} {cont:9.3f} {pct:6.1f}%")
            best[(K,nn)]=(q,obj)
    out={}
    for (K,nn),(q,obj) in best.items():
        tag=f"K{K}_{'neutral' if nn else 'free'}"; out[tag+"_q"]=q; out[tag+"_ddE"]=obj
    np.savez_compressed("sol_discrete.npz", idx=idx, xyz=xyz, shell=shell, dv=dv, sep=sep, **out)
    print(f"\nsaved sol_discrete.npz (sep={sep} A)")
    for tag,(K,nn) in [("K=2 net-free",(2,False)),("K=2 net-neutral",(2,True)),
                       ("K=6 net-neutral",(6,True))]:
        if (K,nn) in best:
            q=best[(K,nn)][0]; sel=np.where(q!=0)[0]
            print(f"\n{tag}:  (min active sep {min_active_sep(xyz,q):.2f} A)")
            for s in sel:
                print(f"  idx {idx[s]}: q={q[s]:+d}, Dv={dv[s]:+.6f}, shell {shell[s]:.1f}A, "
                      f"xyz=({xyz[s][0]:.2f},{xyz[s][1]:.2f},{xyz[s][2]:.2f})")
    print(f"\ngrid min Dv = {dv.min():+.6f} at idx {idx[np.argmin(dv)]}; "
          f"max Dv = {dv.max():+.6f} at idx {idx[np.argmax(dv)]}")
