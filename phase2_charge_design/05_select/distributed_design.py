#!/usr/bin/env python3
# Distributed charge design: minimise the largest single-charge contribution, at fixed catalytic
# effect.
#
# The previous objective maximised barrier lowering, which always concentrates charge into the
# strongest lobe of the field and adds lowering without limit until the predicted barrier goes
# negative. It cannot produce a distributed arrangement at any K, because nothing in it values
# distribution.
#
# Here the catalytic effect is a constraint rather than the objective, held inside a band: at
# least TARGET so the design is worth comparing to the enzyme, and no more than BARE - B_MIN so
# the predicted barrier stays positive by construction rather than being caught afterwards.
# The objective is then the maximum contribution of any single charge. Reaching a fixed target
# with more charges means each carries less, so it falls with K and the optimiser prefers spread
# arrangements. That is a realizability criterion, not a cosmetic one: a design where one group
# carries half the effect fails when that group moves, which is the fragility Arg90 shows.
#
# Placing a passenger charge at a near-zero-field site forces the remaining charges to work
# harder, which raises the maximum, so the formulation excludes passengers without needing a
# separate constraint.

import numpy as np
from mip import Model, xsum, BINARY, CONTINUOUS, MINIMIZE, OptimizationStatus

HARTREE_TO_KCAL = 627.5094740631
MIN_SEP = 2.5
N_CHARGED_POCKET = 6
K_MAX = N_CHARGED_POCKET + 2
SOLVE_LIMIT = 900

# bare reaction, verified from the committed Stage-1 single points
BARE_BARRIER = (-836.340751 - (-836.365884)) * HARTREE_TO_KCAL

# catalytic band. TARGET is the enzyme's experimental electrostatic effect (24.5 - 15.4, Kast
# 1996); B_MIN leaves headroom for the polarisation correction, which ran 12-23% at K=2 and
# always deepens the lowering.
TARGET = 9.1
B_MIN = 5.0
MAX_LOWERING = BARE_BARRIER - B_MIN

GRID = "dv_grid.tsv"
OUT_NPZ = "sol_distributed.npz"
OUT_TSV = "distributed_summary.tsv"
PC_PREFIX = "design_K"
COL_X, COL_SHELL, COL_DV = 1, 4, 7


def load_grid(path):
    xyz, shell, dv = [], [], []
    for line in open(path):
        p = line.split()
        if len(p) <= COL_DV:
            continue
        try:
            xyz.append([float(p[COL_X]), float(p[COL_X + 1]), float(p[COL_X + 2])])
            shell.append(int(float(p[COL_SHELL])))
            dv.append(float(p[COL_DV]))
        except ValueError:
            continue
    return np.array(xyz), np.array(shell), np.array(dv) * HARTREE_TO_KCAL


def exclusion_pairs(xyz, cutoff):
    pairs = []
    cut2 = cutoff * cutoff
    for i in range(len(xyz)):
        dx = xyz[i + 1:, 0] - xyz[i, 0]
        dy = xyz[i + 1:, 1] - xyz[i, 1]
        dz = xyz[i + 1:, 2] - xyz[i, 2]
        d2 = dx * dx + dy * dy + dz * dz
        for off in np.nonzero(d2 < cut2)[0]:
            pairs.append((i, i + 1 + int(off)))
    return pairs


def solve(dv, pairs, K):
    n = len(dv)
    m = Model(sense=MINIMIZE)
    m.verbose = 0
    pos = [m.add_var(var_type=BINARY) for _ in range(n)]
    neg = [m.add_var(var_type=BINARY) for _ in range(n)]
    t = m.add_var(var_type=CONTINUOUS, lb=0.0)

    m.objective = t
    for i in range(n):
        m += pos[i] + neg[i] <= 1
        # t bounds the contribution of any placed charge; sign is free so the magnitude is |dv|
        m += t >= abs(dv[i]) * (pos[i] + neg[i])
    m += xsum(pos[i] + neg[i] for i in range(n)) == K
    total = xsum(dv[i] * (pos[i] - neg[i]) for i in range(n))
    m += total <= -TARGET
    m += total >= -MAX_LOWERING
    for (i, j) in pairs:
        m += pos[i] + neg[i] + pos[j] + neg[j] <= 1

    st = m.optimize(max_seconds=SOLVE_LIMIT)
    if st not in (OptimizationStatus.OPTIMAL, OptimizationStatus.FEASIBLE):
        return {"status": str(st).split(".")[-1]}
    sites, charges = [], []
    for i in range(n):
        if pos[i].x >= 0.99:
            sites.append(i); charges.append(1)
        elif neg[i].x >= 0.99:
            sites.append(i); charges.append(-1)
    contrib = [charges[k] * dv[sites[k]] for k in range(len(sites))]
    return {"status": str(st).split(".")[-1], "t": t.x, "gap": m.gap,
            "sites": sites, "charges": charges, "contrib": contrib,
            "total": float(sum(contrib))}


def write_pc(path, xyz, sites, charges):
    fh = open(path, "w")
    fh.write("%d\n" % len(sites))
    for s, q in zip(sites, charges):
        fh.write("%8.4f %12.6f %12.6f %12.6f\n" % (q, xyz[s][0], xyz[s][1], xyz[s][2]))
    fh.close()


def main():
    xyz, shell, dv = load_grid(GRID)
    pairs = exclusion_pairs(xyz, MIN_SEP)
    print("grid %d points, contribution range %+.3f to %+.3f kcal/mol per unit charge"
          % (len(dv), dv.min(), dv.max()))
    print("bare barrier %.3f; catalytic band %.1f to %.3f kcal/mol lowering"
          % (BARE_BARRIER, TARGET, MAX_LOWERING))
    print("so the predicted barrier is held between %.3f and %.1f\n"
          % (BARE_BARRIER - MAX_LOWERING, BARE_BARRIER - TARGET))

    results, rows = {}, []
    print("  K   max|contrib|  ideal   min|contrib|  spread   total    barrier  gap     net  shells")
    for K in range(2, K_MAX + 1):
        r = solve(dv, pairs, K)
        if "t" not in r:
            print("  %2d   %s" % (K, r["status"]))
            rows.append((K, r["status"], "", "", "", "", "", "", "", ""))
            continue
        mags = [abs(c) for c in r["contrib"]]
        ideal = TARGET / K
        spread = max(mags) / min(mags) if min(mags) > 0 else float("inf")
        barrier = BARE_BARRIER + r["total"]
        net = int(sum(r["charges"]))
        shells = "".join(str(int(shell[i])) for i in r["sites"])
        print("  %2d     %7.3f   %6.3f      %7.3f   %5.2f  %+8.3f  %+7.2f  %.4f  %+d  %s"
              % (K, r["t"], ideal, min(mags), spread, r["total"], barrier, r["gap"], net, shells))
        key = "K%d" % K
        results[key + "_t"] = r["t"]
        results[key + "_total"] = r["total"]
        results[key + "_barrier"] = barrier
        results[key + "_sites"] = np.array(r["sites"])
        results[key + "_charges"] = np.array(r["charges"])
        results[key + "_xyz"] = xyz[np.array(r["sites"], dtype=int)]
        results[key + "_contrib"] = np.array(r["contrib"])
        results[key + "_gap"] = r["gap"]
        rows.append((K, r["status"], "%.4f" % r["t"], "%.4f" % ideal, "%.4f" % min(mags),
                     "%.3f" % spread, "%.4f" % r["total"], "%.4f" % barrier,
                     "%.4f" % r["gap"], "%+d" % net))
        write_pc("%s%d.pc" % (PC_PREFIX, K), xyz, r["sites"], r["charges"])

    np.savez(OUT_NPZ, **results)
    fh = open(OUT_TSV, "w")
    fh.write("K\tstatus\tmax_contrib\tideal\tmin_contrib\tspread\ttotal_ddE\tbarrier\tgap\tnet\n")
    for row in rows:
        fh.write("\t".join(str(v) for v in row) + "\n")
    fh.close()
    print("\nwrote %s, %s and one %sN.pc per feasible K" % (OUT_NPZ, OUT_TSV, PC_PREFIX))


if __name__ == "__main__":
    main()
