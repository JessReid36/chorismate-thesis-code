#!/usr/bin/env python3
"""Layer 3 -- the certified OPTIMISER (our contribution). Clean-room rebuild.

Certified convex L1-sparse charge design on the Layer-2 standoff grid. Objective:
    minimise  sum_i q_i*dV_i + lam1*sum_i|q_i|   s.t.  |q_i|<=QMAX, sum q_i = 0 (if CHARGED).
Integrity is NOT pre-constrained; it is enforced GOCAT's way by SCREENING the relaxed structure
through the L1 oracle. No force-constraints / neutralisers (those compensated for the frozen-charge-
in-the-wall artifact, which the L2 standoff grid removes). Certified LP via native HiGHS.
"""
import os
import numpy as np
import highspy

HERE   = os.path.dirname(os.path.abspath(__file__))
GRID   = os.path.join(HERE, "..", "L2_grid", "grid_final.tsv")
INP    = os.path.join(HERE, "..", "inputs")
REACT  = "reactant.xyz"
QMAX   = 1.0
CHARGED = True
PRUNE  = 1.0e-3
NLAM   = 6
LAM_LO = 0.05
LAM_HI = 0.95


def load_grid():
    rows = open(GRID).read().splitlines()[1:]
    xyz, dv, md = [], [], []
    for ln in rows:
        p = ln.split("\t")
        xyz.append((float(p[1]), float(p[2]), float(p[3])))
        md.append(float(p[4])); dv.append(float(p[8]))
    return np.array(xyz), np.array(dv), np.array(md)


def solve_l1(dV, lam1, qmax, charged):
    N = dV.size; M = 2 * N
    cost = np.concatenate([dV + lam1, -dV + lam1]).astype(np.float64)
    h = highspy.Highs(); h.setOptionValue("output_flag", False)
    h.addVars(M, np.zeros(M), np.full(M, qmax))
    h.changeColsCost(M, np.arange(M, dtype=np.int32), cost)
    if charged:
        idx = np.arange(M, dtype=np.int32)
        h.addRow(0.0, 0.0, M, idx, np.concatenate([np.ones(N), -np.ones(N)]).astype(np.float64))
    h.run()
    st = h.getModelStatus()
    if st != highspy.HighsModelStatus.kOptimal:
        raise RuntimeError("HiGHS: %s" % h.modelStatusToString(st))
    x = np.array(h.getSolution().col_value)
    return x[:N] - x[N:], float(h.getInfo().objective_function_value)


def write_design(tag, sub_lines, xyz, dv, md, q):
    act = np.where(np.abs(q) >= PRUNE)[0]
    n = 24 + len(act)
    with open(os.path.join(HERE, "design_%s_coords.xyz" % tag), "w") as f:
        f.write("%d\ndesign %s: 24 substrate + %d catalytic\n" % (n, tag, len(act)))
        for l in sub_lines:
            f.write(l.rstrip() + "\n")
        for i in act:
            f.write("C   %14.8f %14.8f %14.8f\n" % (xyz[i, 0], xyz[i, 1], xyz[i, 2]))
    with open(os.path.join(HERE, "design_%s_charges.txt" % tag), "w") as f:
        f.write(",".join(["0.0000"] * 24 + ["%.4f" % q[i] for i in act]))
    with open(os.path.join(HERE, "design_%s_report.tsv" % tag), "w") as f:
        f.write("grid_idx\tmin_dist\tdV\tq\n")
        for i in act:
            f.write("%d\t%.3f\t%+.6f\t%+.4f\n" % (i, md[i], dv[i], q[i]))
    return act


def dist_hist(md, act):
    if len(act) == 0:
        return ""
    v, c = np.unique(np.floor(md[act]).astype(int), return_counts=True)
    return " ".join("%dA:%d" % (a, b) for a, b in zip(v, c))


def main():
    xyz, dv, md = load_grid()
    sub_lines = open(os.path.join(INP, REACT)).read().splitlines()[2:26]
    dVmax = abs(dv).max()
    lambdas = np.geomspace(LAM_LO, LAM_HI, NLAM) * dVmax

    print("grid: %d sites | min_dist %.2f-%.2f A (LJ standoff) | |dV|max %.5f | QMAX=%.1f CHARGED=%s"
          % (len(dv), md.min(), md.max(), dVmax, QMAX, CHARGED))
    print("certified L1-sparse LP (native HiGHS). integrity enforced downstream by the L1 oracle.\n")
    print("tag  lam1/|dV|max  obj           net_q     #active  min_dist histogram")

    summary = open(os.path.join(HERE, "explore_summary.tsv"), "w")
    summary.write("tag\tlam1_frac\tlam1\tobj\tnet_q\tn_active\tdist_hist\n")
    for k, lam1 in enumerate(lambdas):
        tag = "s%d" % k
        q, obj = solve_l1(dv, lam1, QMAX, CHARGED)
        act = write_design(tag, sub_lines, xyz, dv, md, q)
        hist = dist_hist(md, act)
        print("%-4s %.3f         %+.5e  %+.1e  %6d   %s" % (tag, lam1 / dVmax, obj, q.sum(), len(act), hist))
        summary.write("%s\t%.4f\t%.6e\t%.6e\t%.2e\t%d\t%s\n"
                      % (tag, lam1 / dVmax, lam1, obj, q.sum(), len(act), hist))
    summary.close()
    print("\nwrote design_s{0..%d}_{coords.xyz,charges.txt,report.tsv} + explore_summary.tsv" % (NLAM - 1))


if __name__ == "__main__":
    main()
