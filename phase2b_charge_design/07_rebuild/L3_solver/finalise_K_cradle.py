#!/usr/bin/env python3
"""Layer 3b (cradle) -- certified K catalytic charges ON TOP of the fixed counter-charge cradle."""
import os, sys
import numpy as np
import highspy

HERE   = os.path.dirname(os.path.abspath(__file__))
GRID   = os.path.join(HERE, "..", "L2_grid", "grid_final.tsv")
INP    = os.path.join(HERE, "..", "inputs")
REACT  = "reactant.xyz"
CRADLE = os.path.join(HERE, "design_cradle2_coords.xyz")
K      = int(sys.argv[1]) if len(sys.argv) > 1 else 6
QMAX   = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
PRESELECT = 2000
CLEAR_OF_CRADLE = 2.0


def load_grid():
    rows = open(GRID).read().splitlines()[1:]
    xyz, dv, md = [], [], []
    for ln in rows:
        p = ln.split("\t")
        xyz.append((float(p[1]), float(p[2]), float(p[3])))
        md.append(float(p[4])); dv.append(float(p[8]))
    return np.array(xyz), np.array(dv), np.array(md)


def load_cradle():
    L = open(CRADLE).read().splitlines()
    n = int(L[0].split()[0])
    return np.array([[float(v) for v in l.split()[1:4]] for l in L[2 + 24:2 + n]])


def solve_hardK(dV, K, qmax):
    N = dV.size
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("mip_rel_gap", 1e-2)
    h.setOptionValue("time_limit", 120.0)
    lower = np.concatenate([np.full(N, -qmax), np.zeros(N)])
    upper = np.concatenate([np.full(N,  qmax), np.ones(N)])
    cost  = np.concatenate([dV, np.zeros(N)]).astype(np.float64)
    h.addVars(2 * N, lower, upper)
    h.changeColsCost(2 * N, np.arange(2 * N, dtype=np.int32), cost)
    integ = np.concatenate([np.zeros(N), np.ones(N)]).astype(np.int32)
    h.changeColsIntegrality(2 * N, np.arange(2 * N, dtype=np.int32), integ)
    for i in range(N):
        h.addRow(-1e30, 0.0, 2, np.array([i, N + i], np.int32), np.array([1.0, -qmax]))
        h.addRow(0.0, 1e30, 2, np.array([i, N + i], np.int32), np.array([1.0,  qmax]))
    h.addRow(float(K), float(K), N, (np.arange(N) + N).astype(np.int32), np.ones(N))
    h.addRow(0.0, 0.0, N, np.arange(N, dtype=np.int32), np.ones(N))
    h.run()
    x = np.array(h.getSolution().col_value)
    info = h.getInfo()
    return x[:N], float(info.objective_function_value), getattr(info, "mip_gap", 0.0)


def main():
    xyz, dv, md = load_grid()
    cradle = load_cradle()
    clear = np.array([np.linalg.norm(cradle - p, axis=1).min() >= CLEAR_OF_CRADLE for p in xyz])
    xyz, dv, md = xyz[clear], dv[clear], md[clear]
    order = np.argsort(-np.abs(dv))[:PRESELECT]
    xs, ds, ms = xyz[order], dv[order], md[order]
    print("grid %d (clear of cradle) -> MILP over top %d by |dV| | K_catalytic=%d | QMAX=%.1f"
          % (len(dv), len(order), K, QMAX))
    q, obj, gap = solve_hardK(ds, K, QMAX)
    act = np.where(np.abs(q) >= 1e-4)[0]
    print("certified MILP: obj=%.6e | gap=%.2e | catalytic active=%d (K=%d) | net_q_cat=%+.2e"
          % (obj, gap, len(act), K, q.sum()))
    sub_lines = open(os.path.join(INP, REACT)).read().splitlines()[2:26]
    tag = "cradleK%dq%s" % (K, ("%.2f" % QMAX).replace(".", "p"))
    ncat = len(act)
    with open(os.path.join(HERE, "design_%s_coords.xyz" % tag), "w") as f:
        f.write("%d\ncradle(2x+1) + %d catalytic (certified, gap %.1e)\n" % (24 + 2 + ncat, ncat, gap))
        for l in sub_lines:
            f.write(l.rstrip() + "\n")
        for p in cradle:
            f.write("C   %14.8f %14.8f %14.8f\n" % (p[0], p[1], p[2]))
        for i in act:
            f.write("C   %14.8f %14.8f %14.8f\n" % (xs[i, 0], xs[i, 1], xs[i, 2]))
    with open(os.path.join(HERE, "design_%s_charges.txt" % tag), "w") as f:
        vals = ["0.0000"] * 24 + ["1.0000", "1.0000"] + ["%.4f" % q[i] for i in act]
        f.write(",".join(vals))
    with open(os.path.join(HERE, "design_%s_report.tsv" % tag), "w") as f:
        f.write("global_idx\tmin_dist\tdV\tq\n")
        for i in act:
            f.write("%d\t%.3f\t%+.6f\t%+.4f\n" % (order[i], ms[i], ds[i], q[i]))
    v, c = np.unique(np.floor(ms[act]).astype(int), return_counts=True)
    print("catalytic placement (min-dist): " + " ".join("%dA:%d" % (a, b) for a, b in zip(v, c)))
    print("wrote design_%s_* (24 substrate + 2 counter + %d catalytic)" % (tag, ncat))


if __name__ == "__main__":
    main()
