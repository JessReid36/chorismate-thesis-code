#!/usr/bin/env python3
# Charge-count sweep under progressive standoff restriction.
#
# Every design from the unrestricted sweep placed all of its charges on the innermost grid shell.
# A point charge that close to the substrate density is the standard electron spill-out risk once
# the density is allowed to respond, so those designs cannot be carried into a self-consistent
# geometry optimisation as they stand. This sweep re-solves with the inner shells progressively
# excluded, to measure what a defensible standoff costs.
#
# It also tests a second thing. The frozen-density objective is a first-order expansion with no
# lower bound, so it happily predicts negative barriers once the perturbation is large enough.
# Pushing charges outward weakens the perturbation, which should pull the higher-K designs back
# toward the regime where the prediction means anything. The implied-barrier column is that test:
# a positive barrier is necessary for a design to be reportable, though not sufficient.

import numpy as np
from mip import Model, xsum, BINARY, MINIMIZE, OptimizationStatus

HARTREE_TO_KCAL = 627.5094740631
MIN_SEP = 2.5                 # A; below this two designed charges would occupy the same space
N_CHARGED_POCKET = 6          # Arg7, Lys60', Arg63', Arg90, Arg116 cationic; Glu78 anionic
K_MAX = N_CHARGED_POCKET + 2
SOLVE_LIMIT = 900

# bare reaction, from the committed Stage-1 single points (B3LYP-D3BJ/def2-SVP, CPCM eps=4)
E_REACTANT_HA = -836.365884
E_TS_HA = -836.340751
BARE_BARRIER = (E_TS_HA - E_REACTANT_HA) * HARTREE_TO_KCAL

GRID = "dv_grid.tsv"
OUT_NPZ = "sol_shell_restricted.npz"
OUT_TSV = "shell_restricted_summary.tsv"

COL_X, COL_SHELL, COL_DV = 1, 4, 7


def load_grid(path):
    xyz, shell, dv = [], [], []
    fh = open(path)
    for line in fh:
        parts = line.split()
        if len(parts) <= COL_DV:
            continue
        try:
            x = float(parts[COL_X])
            y = float(parts[COL_X + 1])
            z = float(parts[COL_X + 2])
            s = int(float(parts[COL_SHELL]))
            d = float(parts[COL_DV])
        except ValueError:
            continue
        xyz.append((x, y, z))
        shell.append(s)
        dv.append(d)
    fh.close()
    return np.array(xyz), np.array(shell), np.array(dv)


def exclusion_pairs(xyz, allowed, cutoff):
    # exclusion pairs among the allowed sites only, indexed into the allowed subset
    pairs = []
    cut2 = cutoff * cutoff
    sub = xyz[allowed]
    for i in range(len(sub)):
        dx = sub[i + 1:, 0] - sub[i, 0]
        dy = sub[i + 1:, 1] - sub[i, 1]
        dz = sub[i + 1:, 2] - sub[i, 2]
        d2 = dx * dx + dy * dy + dz * dz
        for off in np.nonzero(d2 < cut2)[0]:
            pairs.append((i, i + 1 + int(off)))
    return pairs


def solve(dv_sub, pairs, K):
    n = len(dv_sub)
    if K > n:
        return None
    m = Model(sense=MINIMIZE)
    m.verbose = 0
    pos = [m.add_var(var_type=BINARY) for _ in range(n)]
    neg = [m.add_var(var_type=BINARY) for _ in range(n)]
    m.objective = xsum(dv_sub[i] * (pos[i] - neg[i]) for i in range(n))
    for i in range(n):
        m += pos[i] + neg[i] <= 1
    m += xsum(pos[i] + neg[i] for i in range(n)) == K
    for (i, j) in pairs:
        m += pos[i] + neg[i] + pos[j] + neg[j] <= 1
    status = m.optimize(max_seconds=SOLVE_LIMIT)
    if status not in (OptimizationStatus.OPTIMAL, OptimizationStatus.FEASIBLE):
        return None
    local, charges = [], []
    for i in range(n):
        if pos[i].x >= 0.99:
            local.append(i); charges.append(1)
        elif neg[i].x >= 0.99:
            local.append(i); charges.append(-1)
    return {"ddE": m.objective_value * HARTREE_TO_KCAL, "gap": m.gap,
            "local": local, "charges": charges}


def main():
    xyz, shell, dv = load_grid(GRID)
    shells = sorted(set(int(s) for s in shell))
    print("grid %d points, Dv %+.6f to %+.6f Ha" % (len(dv), dv.min(), dv.max()))
    print("bare barrier %.3f kcal/mol\n" % BARE_BARRIER)

    print("field strength by shell:")
    for s in shells:
        sel = shell == s
        print("  shell %d  %3d points   Dv %+.6f to %+.6f   max|Dv| %.6f"
              % (s, sel.sum(), dv[sel].min(), dv[sel].max(), np.abs(dv[sel]).max()))
    print("")

    results = {}
    rows = []
    for smin in shells:
        allowed = np.nonzero(shell >= smin)[0]
        dv_sub = dv[allowed]
        pairs = exclusion_pairs(xyz, allowed, MIN_SEP)
        print("shells >= %d   %d sites available, %d exclusion pairs" % (smin, len(allowed), len(pairs)))
        print("   K   ddE        barrier    gap     net  shells")
        for K in range(2, K_MAX + 1):
            r = solve(dv_sub, pairs, K)
            if r is None:
                print("  %2d   infeasible" % K)
                rows.append((smin, K, "INFEASIBLE", "", "", "", "", ""))
                continue
            glob = allowed[np.array(r["local"], dtype=int)]
            barrier = BARE_BARRIER + r["ddE"]
            net = int(np.sum(r["charges"]))
            used = "".join(str(int(shell[i])) for i in glob)
            flag = "" if barrier > 0 else "   <- unphysical"
            print("  %2d  %+8.3f  %+8.2f  %.4f  %+d  %s%s"
                  % (K, r["ddE"], barrier, r["gap"], net, used, flag))
            key = "s%d_K%d" % (smin, K)
            results[key + "_ddE"] = r["ddE"]
            results[key + "_barrier"] = barrier
            results[key + "_sites"] = glob
            results[key + "_charges"] = np.array(r["charges"])
            results[key + "_xyz"] = xyz[glob]
            results[key + "_gap"] = r["gap"]
            rows.append((smin, K, "OK", "%.4f" % r["ddE"], "%.4f" % barrier,
                         "%.4f" % r["gap"], "%+d" % net, used))
        print("")

    print("designs leaving a positive barrier (the reportable set):")
    any_ok = False
    for smin in shells:
        for K in range(2, K_MAX + 1):
            b = results.get("s%d_K%d_barrier" % (smin, K))
            if b is not None and b > 0:
                any_ok = True
                print("  shells>=%d  K=%d   ddE %+8.3f   barrier %+6.2f"
                      % (smin, K, results["s%d_K%d_ddE" % (smin, K)], b))
    if not any_ok:
        print("  none")

    np.savez(OUT_NPZ, **results)
    fh = open(OUT_TSV, "w")
    fh.write("min_shell\tK\tstatus\tddE_kcal\timplied_barrier\tgap\tnet_charge\tshells_used\n")
    for row in rows:
        fh.write("\t".join(str(v) for v in row) + "\n")
    fh.close()
    print("\nwrote %s and %s" % (OUT_NPZ, OUT_TSV))


if __name__ == "__main__":
    main()
