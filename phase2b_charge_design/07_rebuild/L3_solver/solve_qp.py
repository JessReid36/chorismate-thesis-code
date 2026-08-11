#!/usr/bin/env python3
"""Tier-1 convex QP for VARIED-MAGNITUDE external point charges (certified, via CVXPY)."""
import os, sys
import numpy as np

HERE   = os.path.dirname(os.path.abspath(__file__))
GRID   = os.path.join(HERE, "..", "L2_grid", "grid_final.tsv")
INP    = os.path.join(HERE, "..", "inputs")
REACT  = "reactant.xyz"
TS     = "ts.xyz"
QMAX   = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
LAM    = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
MU     = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
SOLVER = sys.argv[4] if len(sys.argv) > 4 else "CLARABEL"
PRESELECT = 2000
CLEAR_OF_CRADLE = 2.0
REACTING_BONDS = [(7, 8), (0, 12)]


def load_grid():
    rows = open(GRID).read().splitlines()[1:]
    xyz, dv, md = [], [], []
    for ln in rows:
        p = ln.split("\t")
        xyz.append((float(p[1]), float(p[2]), float(p[3])))
        md.append(float(p[4])); dv.append(float(p[8]))
    return np.array(xyz), np.array(dv), np.array(md)


def load_sub(fname):
    L = open(os.path.join(INP, fname)).read().splitlines()[2:26]
    return np.array([[float(v) for v in l.split()[1:4]] for l in L])


def field_gram(sites):
    R = load_sub(REACT)
    B = np.zeros((len(REACTING_BONDS), len(sites)))
    for b, (i, j) in enumerate(REACTING_BONDS):
        m = 0.5 * (R[i] + R[j])
        u = R[j] - R[i]; u = u / np.linalg.norm(u)
        d = m - sites
        r = np.linalg.norm(d, axis=1)
        B[b] = (d @ u) / (r ** 3 + 1e-12)
    return B.T @ B


def load_cradle_positions():
    p = os.path.join(HERE, "design_cradle2_coords.xyz")
    if not os.path.exists(p):
        return None
    L = open(p).read().splitlines(); n = int(L[0].split()[0])
    return np.array([[float(v) for v in l.split()[1:4]] for l in L[2 + 24:2 + n]])


def main():
    import cvxpy as cp
    xyz, dv, md = load_grid()
    cradle = load_cradle_positions()
    if cradle is not None:
        clear = np.array([np.linalg.norm(cradle - p, axis=1).min() >= CLEAR_OF_CRADLE for p in xyz])
        xyz, dv, md = xyz[clear], dv[clear], md[clear]
    order = np.argsort(-np.abs(dv))[:PRESELECT]
    xs, ds, ms = xyz[order], dv[order], md[order]
    N = len(ds)
    Wf = field_gram(xs)
    Wf = 0.5 * (Wf + Wf.T)
    top = float(np.linalg.eigvalsh(Wf).max()) or 1.0
    GAMMA = 5.0
    W = np.eye(N) + GAMMA * (Wf / top)  # identity ridge regularises ALL charges; field term adds bond penalty
    q = cp.Variable(N)
    obj = ds @ q + LAM * cp.quad_form(q, cp.psd_wrap(W))
    if MU > 0:
        obj = obj + MU * cp.norm1(q)
    constraints = [cp.sum(q) == 0, q >= -QMAX, q <= QMAX]
    prob = cp.Problem(cp.Minimize(obj), constraints)
    prob.solve(solver=getattr(cp, SOLVER))
    print("solver=%s status=%s optval=%.6e" % (SOLVER, prob.status, prob.value))
    if q.value is None:
        sys.exit("no solution")
    qv = np.asarray(q.value).ravel()
    act = np.where(np.abs(qv) >= 1e-3)[0]
    qdV = float(ds[act] @ qv[act])
    print("active charges=%d | net_q=%+.2e | Sum q.dV=%.6e | |q| range [%.3f, %.3f]"
          % (len(act), qv.sum(), qdV, np.abs(qv[act]).min() if len(act) else 0, np.abs(qv[act]).max() if len(act) else 0))
    uq = np.round(np.sort(np.abs(qv[act]))[::-1], 3)
    print("charge magnitudes (varied, not all %.2f): %s" % (QMAX, uq[:20]))
    sub_lines = open(os.path.join(INP, REACT)).read().splitlines()[2:26]
    tag = "qp_q%s_l%s%s" % (("%.2f" % QMAX).replace(".", "p"),
                            ("%.4f" % LAM).replace(".", "p"),
                            ("_m%s" % ("%.4f" % MU).replace(".", "p")) if MU > 0 else "")
    ncat = len(act)
    with open(os.path.join(HERE, "design_%s_coords.xyz" % tag), "w") as f:
        ntot = 24 + (len(cradle) if cradle is not None else 0) + ncat
        f.write("%d\nTier-1 QP varied charges (lam=%.2f mu=%.2f, %s)\n" % (ntot, LAM, MU, SOLVER))
        for l in sub_lines:
            f.write(l.rstrip() + "\n")
        if cradle is not None:
            for p in cradle:
                f.write("C   %14.8f %14.8f %14.8f\n" % (p[0], p[1], p[2]))
        for i in act:
            f.write("C   %14.8f %14.8f %14.8f\n" % (xs[i, 0], xs[i, 1], xs[i, 2]))
    with open(os.path.join(HERE, "design_%s_charges.txt" % tag), "w") as f:
        vals = ["0.0000"] * 24
        if cradle is not None:
            vals += ["1.0000"] * len(cradle)
        vals += ["%.4f" % qv[i] for i in act]
        f.write(",".join(vals))
    with open(os.path.join(HERE, "design_%s_report.tsv" % tag), "w") as f:
        f.write("global_idx\tmin_dist\tdV\tq\n")
        for i in act:
            f.write("%d\t%.3f\t%+.6f\t%+.4f\n" % (order[i], ms[i], ds[i], qv[i]))
    print("wrote design_%s_* (%d catalytic charges, varied magnitudes)" % (tag, ncat))


if __name__ == "__main__":
    main()
