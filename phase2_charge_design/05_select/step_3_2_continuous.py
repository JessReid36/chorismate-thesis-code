#!/usr/bin/env python3
# Phase 2 step 3.2 - CHARGE-BUDGETED CONTINUOUS optimizer (option a, done right).
# dDE = sum_i q_i*Dv_i is unbounded without a limit on deployed charge; the well-posed problem
# constrains total charge via an L1 budget sum|q_i| <= Q (net-neutral, bounded per point), which
# also induces sparsity (LASSO): the solver CHOOSES where a limited budget does the most good.
# Sweeping Q gives the barrier-lowering vs charge-deployed curve. Convex; global optimum per Q.
import sys, numpy as np, cvxpy as cp

def load_dv(path):
    idx, xyz, shell, dv = [], [], [], []
    for ln in open(path).read().splitlines()[1:]:
        p = ln.split("\t")
        idx.append(int(p[0])); xyz.append([float(p[1]),float(p[2]),float(p[3])])
        shell.append(float(p[4])); dv.append(float(p[7]))
    return np.array(idx), np.array(xyz), np.array(shell), np.array(dv)

def solve_budget(dv, Q, qmax=1.0):
    n = len(dv); q = cp.Variable(n)
    obj = cp.Minimize(dv @ q)
    cons = [cp.norm1(q) <= Q, cp.sum(q) == 0, q <= qmax, q >= -qmax]
    prob = cp.Problem(obj, cons); prob.solve()
    return q.value, prob.status

def describe(q, dv, Q):
    ddE = float(dv @ q)
    nnz = int((np.abs(q) > 1e-4).sum())
    pos = int((q > 1e-4).sum()); neg = int((q < -1e-4).sum())
    return dict(Q=Q, ddE_Eh=ddE, ddE_kcal=ddE*627.509, nnz=nnz, npos=pos, nneg=neg,
                l1_used=float(np.abs(q).sum()), maxq=float(np.abs(q).max()), sumq=float(q.sum()))

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "dv_grid.tsv"
    qmax = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    idx, xyz, shell, dv = load_dv(path)
    print(f"loaded {len(dv)} grid points | Dv range [{dv.min():+.6f},{dv.max():+.6f}] Eh, "
          f"{(dv<0).sum()} stabilizing\n")
    budgets = [1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0]
    print(f"{'Q(e)':>5} {'dDE(Eh)':>11} {'dDE(kcal/mol)':>14} {'nonzero':>8} {'(+/-)':>9} {'L1used':>7} {'max|q|':>7}")
    sols = {}
    for Q in budgets:
        q, status = solve_budget(dv, Q, qmax)
        if q is None:
            print(f"{Q:5.1f}  solve failed ({status})"); continue
        d = describe(q, dv, Q)
        print(f"{Q:5.1f} {d['ddE_Eh']:+11.6f} {d['ddE_kcal']:+14.3f} {d['nnz']:8d} "
              f"{d['npos']:4d}/{d['nneg']:<4d} {d['l1_used']:7.2f} {d['maxq']:7.3f}")
        sols[Q] = q
    np.savez_compressed("sol_continuous_budget.npz", idx=idx, xyz=xyz, shell=shell, dv=dv,
                        budgets=np.array(list(sols.keys())),
                        solutions=np.array([sols[Q] for Q in sols]))
    print("\nsaved sol_continuous_budget.npz (barrier-vs-budget sweep)")
    kcs = [describe(sols[Q],dv,Q)['ddE_kcal'] for Q in sols]
    print(f"monotonic (barrier lowers as budget grows): {all(kcs[i] >= kcs[i+1]-1e-6 for i in range(len(kcs)-1))}")
    if 2.0 in sols:
        q2 = sols[2.0]; toppos = np.argsort(q2)[::-1][:1]
        print(f"at Q=2e: {int((np.abs(q2)>1e-4).sum())} charges active; most-stabilizing cation at "
              f"grid idx {idx[toppos[0]]} (Dv={dv[toppos[0]]:+.6f}); min Dv over grid = {dv.min():+.6f} at idx {idx[np.argmin(dv)]}")
