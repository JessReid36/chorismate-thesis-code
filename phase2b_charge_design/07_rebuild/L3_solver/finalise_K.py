#!/usr/bin/env python3
"""Layer 3b -- certified hard-K FINALISER (buildable catalyst; embodies Claim 1).

Selects EXACTLY K charge sites from the full L2 standoff grid minimising the Sokalski differential-
stabilisation objective, with a provable optimality certificate (MILP, gap -> 0). Certified discrete
counterpart to GOCAT's bounded-size embeddings: they find a small charge set by GA; we find the
provably-optimal K-charge set by MILP.

  vars : q_i in [-QMAX,QMAX], z_i in {0,1}
  min  : sum_i q_i dV_i
  s.t. : -QMAX z_i <= q_i <= QMAX z_i ;  sum z_i = K ;  sum q_i = 0
"""
import os, sys
import numpy as np
import highspy

HERE   = os.path.dirname(os.path.abspath(__file__))
GRID   = os.path.join(HERE, "..", "L2_grid", "grid_final.tsv")
INP    = os.path.join(HERE, "..", "inputs")
REACT  = "reactant.xyz"
QMAX   = 1.0
K      = int(sys.argv[1]) if len(sys.argv) > 1 else 20
CARD_NEUTRAL = True
PRESELECT = 2000


def load_grid():
    rows = open(GRID).read().splitlines()[1:]
    xyz, dv, md = [], [], []
    for ln in rows:
        p = ln.split("\t")
        xyz.append((float(p[1]), float(p[2]), float(p[3])))
        md.append(float(p[4])); dv.append(float(p[8]))
    return np.array(xyz), np.array(dv), np.array(md)


def solve_hardK(dV, K, qmax, neutral):
    N = dV.size
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("mip_rel_gap", 0.0)
    h.setOptionValue("mip_abs_gap", 0.0)
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
    if neutral:
        h.addRow(0.0, 0.0, N, np.arange(N, dtype=np.int32), np.ones(N))
    h.run()
    st = h.getModelStatus()
    if st != highspy.HighsModelStatus.kOptimal:
        raise RuntimeError("HiGHS MILP: %s" % h.modelStatusToString(st))
    x = np.array(h.getSolution().col_value)
    info = h.getInfo()
    return x[:N], float(info.objective_function_value), getattr(info, "mip_gap", 0.0)


def main():
    xyz, dv, md = load_grid()
    order = np.argsort(-np.abs(dv))[:PRESELECT]
    xs, ds, ms = xyz[order], dv[order], md[order]
    print("full grid %d -> MILP over top %d by |dV| | K=%d | QMAX=%.1f | neutral=%s"
          % (len(dv), len(order), K, QMAX, CARD_NEUTRAL))
    q, obj, gap = solve_hardK(ds, K, QMAX, CARD_NEUTRAL)
    act = np.where(np.abs(q) >= 1e-4)[0]
    print("certified MILP: obj=%.6e | gap=%.2e | active=%d (K=%d) | net_q=%+.2e"
          % (obj, gap, len(act), K, q.sum()))
    sub_lines = open(os.path.join(INP, REACT)).read().splitlines()[2:26]
    tag = "K%d" % K
    with open(os.path.join(HERE, "design_%s_coords.xyz" % tag), "w") as f:
        f.write("%d\ncertified K=%d design (gap %.1e)\n" % (24 + len(act), K, gap))
        for l in sub_lines:
            f.write(l.rstrip() + "\n")
        for i in act:
            f.write("C   %14.8f %14.8f %14.8f\n" % (xs[i, 0], xs[i, 1], xs[i, 2]))
    with open(os.path.join(HERE, "design_%s_charges.txt" % tag), "w") as f:
        f.write(",".join(["0.0000"] * 24 + ["%.4f" % q[i] for i in act]))
    with open(os.path.join(HERE, "design_%s_report.tsv" % tag), "w") as f:
        f.write("global_idx\tmin_dist\tdV\tq\n")
        for i in act:
            f.write("%d\t%.3f\t%+.6f\t%+.4f\n" % (order[i], ms[i], ds[i], q[i]))
    v, c = np.unique(np.floor(ms[act]).astype(int), return_counts=True)
    print("charge placement (min-dist): " + " ".join("%dA:%d" % (a, b) for a, b in zip(v, c)))
    print("wrote design_%s_{coords.xyz,charges.txt,report.tsv}" % tag)


if __name__ == "__main__":
    main()
