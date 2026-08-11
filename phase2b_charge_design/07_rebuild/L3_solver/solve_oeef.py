#!/usr/bin/env python3
"""OEEF optimiser: certified charges making a ~UNIFORM catalytic -z field over the reacting region."""
import os, sys
import numpy as np

HERE  = os.path.dirname(os.path.abspath(__file__))
GRID  = os.path.join(HERE, "..", "L2_grid", "grid_final.tsv")
INP   = os.path.join(HERE, "..", "inputs")
REACT = "reactant.xyz"
QMAX    = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
FMAX    = float(sys.argv[2]) if len(sys.argv) > 2 else 0.015
FTARGET = float(sys.argv[3]) if len(sys.argv) > 3 else 0.03
SOLVER  = sys.argv[4] if len(sys.argv) > 4 else "CLARABEL"
PRESELECT = 3000
CLEAR_OF_CRADLE = 2.0
BOHR = 0.529177
D_HAT = np.array([0.0, 0.0, -1.0])
REACTING = [7, 8, 0, 12]
REACTING_BONDS = [(7, 8), (0, 12)]


def load_grid():
    rows = open(GRID).read().splitlines()[1:]
    xyz, md = [], []
    for ln in rows:
        p = ln.split("\t")
        xyz.append((float(p[1]), float(p[2]), float(p[3]))); md.append(float(p[4]))
    return np.array(xyz), np.array(md)


def load_sub():
    L = open(os.path.join(INP, REACT)).read().splitlines()[2:26]
    return np.array([[float(v) for v in l.split()[1:4]] for l in L])


def field_matrix_at(sites, point):
    d = (point - sites) / BOHR
    r = np.linalg.norm(d, axis=1)
    return (d / (r[:, None] ** 3)).T


def field_matrix_bonds(sites):
    R = load_sub()
    B = np.zeros((len(REACTING_BONDS), len(sites)))
    for b, (i, j) in enumerate(REACTING_BONDS):
        m = 0.5 * (R[i] + R[j]); u = R[j] - R[i]; u = u / np.linalg.norm(u)
        d = (m - sites) / BOHR; r = np.linalg.norm(d, axis=1)
        B[b] = (d @ u) / (r ** 3)
    return B


def main():
    import cvxpy as cp
    xyz, md = load_grid()
    R = load_sub()
    cradle = None
    cp_file = os.path.join(HERE, "design_cradle2_coords.xyz")
    if os.path.exists(cp_file):
        L = open(cp_file).read().splitlines(); n = int(L[0].split()[0])
        cradle = np.array([[float(v) for v in l.split()[1:4]] for l in L[2 + 24:2 + n]])
        clear = np.array([np.linalg.norm(cradle - p, axis=1).min() >= CLEAR_OF_CRADLE for p in xyz])
        xyz, md = xyz[clear], md[clear]
    centroid = R[REACTING].mean(axis=0)
    dcen = np.linalg.norm(xyz - centroid, axis=1)
    order = np.argsort(dcen)[:PRESELECT]
    xs, ms = xyz[order], md[order]
    N = len(xs)

    Catoms = [field_matrix_at(xs, R[a]) for a in REACTING]
    czbar = sum(D_HAT @ Ca for Ca in Catoms) / len(Catoms)
    B = field_matrix_bonds(xs)

    q = cp.Variable(N)
    constraints = [cp.sum(q) == 0, q >= -QMAX, q <= QMAX, cp.norm(B @ q, 2) <= FMAX]
    for Ca in Catoms:
        constraints.append(cp.norm(Ca @ q, 2) <= FTARGET)
    prob = cp.Problem(cp.Maximize(czbar @ q), constraints)
    prob.solve(solver=getattr(cp, SOLVER))
    achieved = float(czbar @ q.value) if q.value is not None else 0.0
    print("solver=%s status=%s  avg -z field over reacting atoms = %.4f a.u. (target ~0.02)" %
          (SOLVER, prob.status, achieved))
    if q.value is None:
        sys.exit("no solution")
    qv = np.asarray(q.value).ravel()
    act = np.where(np.abs(qv) >= 1e-3)[0]
    print("-z field at each reacting atom: " +
          " ".join("%.4f" % float(D_HAT @ (Ca @ qv)) for Ca in Catoms))
    print("reacting-bond field norm ||Bq|| = %.4f (cap %.3f)" % (np.linalg.norm(B @ qv), FMAX))
    print("active charges=%d | net_q=%+.2e | |q| range [%.3f, %.3f]" %
          (len(act), qv.sum(), np.abs(qv[act]).min() if len(act) else 0,
           np.abs(qv[act]).max() if len(act) else 0))

    sub_lines = open(os.path.join(INP, REACT)).read().splitlines()[2:26]
    tag = "oeef_q%s_f%s" % (("%.2f" % QMAX).replace(".", "p"), ("%.3f" % FMAX).replace(".", "p"))
    ncr = len(cradle) if cradle is not None else 0
    with open(os.path.join(HERE, "design_%s_coords.xyz" % tag), "w") as f:
        f.write("%d\nOEEF uniform -z field + cradle(%d)\n" % (24 + ncr + len(act), ncr))
        for l in sub_lines:
            f.write(l.rstrip() + "\n")
        if cradle is not None:
            for p in cradle:
                f.write("C   %14.8f %14.8f %14.8f\n" % (p[0], p[1], p[2]))
        for i in act:
            f.write("C   %14.8f %14.8f %14.8f\n" % (xs[i, 0], xs[i, 1], xs[i, 2]))
    with open(os.path.join(HERE, "design_%s_charges.txt" % tag), "w") as f:
        vals = ["0.0000"] * 24 + (["1.0000"] * ncr) + ["%.4f" % qv[i] for i in act]
        f.write(",".join(vals))
    with open(os.path.join(HERE, "design_%s_report.tsv" % tag), "w") as f:
        f.write("global_idx\tmin_dist\tq\n")
        for i in act:
            f.write("%d\t%.3f\t%+.4f\n" % (order[i], ms[i], qv[i]))
    print("wrote design_%s_* (%d catalytic charges, uniform -z field)" % (tag, len(act)))


if __name__ == "__main__":
    main()
