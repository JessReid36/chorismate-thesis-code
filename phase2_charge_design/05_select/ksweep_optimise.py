#!/usr/bin/env python3
# Charge-count sweep of the certified charge-placement optimiser.
#
# Two curves are reported from a single solve per (K, constraint), because they are related by
# an exact scaling. For charges of equal magnitude m with free sign, the objective is
# m * sum_i s_i * dv_i with s_i in {-1,+1}; m > 0 is constant for fixed K, so the optimal
# configuration is independent of m and only the objective scales:
#
#   unit-charge curve : every placed charge has |q| = 1        -> ddE_unit(K)
#   budget curve      : total |charge| held at Q, so |q| = Q/K -> ddE_budget(K) = (Q/K)*ddE_unit(K)
#
# These answer different questions and must not be mixed. The unit-charge curve is the design
# series (each charge = one ionised group, directly comparable to the enzyme's charged residues).
# The budget curve is the dilution series: the same total charge spread more thinly, which is the
# only formulation under which adding charges can worsen the result.
#
# The solve is done at unit magnitude and scaled afterwards rather than solving the budget problem
# directly: at K=8 the budget coefficients fall to ~1e-4 and CBC returns a slightly sub-optimal
# solution within its tolerance.
#
# Net-free is the design constraint. Net-neutral is retained only as a reconciliation control
# against the committed run, and is infeasible for odd K since an odd number of +-1 charges
# cannot sum to zero.

import numpy as np
from mip import Model, xsum, BINARY, MINIMIZE, OptimizationStatus

HARTREE_TO_KCAL = 627.5094740631
MIN_SEP = 2.5                 # A; below this two designed charges would occupy the same space
Q_BUDGET = 2.0                # e; matches the total charge of the committed continuous ceiling
N_CHARGED_POCKET = 6          # Arg7, Lys60', Arg63', Arg90, Arg116 cationic; Glu78 anionic
K_MAX = N_CHARGED_POCKET + 2
SOLVE_LIMIT = 900             # s per solve; hitting this voids that solve's certificate

GRID = "dv_grid.tsv"
OUT_NPZ = "sol_ksweep.npz"
OUT_TSV = "ksweep_summary.tsv"

COL_X, COL_SHELL, COL_DV = 1, 4, 7

# committed budget-curve values the reconstruction must reproduce before new K are trusted.
# K2 is exact; K6/K8 are back-derived from recorded percentages and compared loosely.
REFERENCE = {2: (-10.888, 0.002), 6: (-10.11, 0.15), 8: (-9.76, 0.15)}


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
            continue                       # header line
        xyz.append((x, y, z))
        shell.append(s)
        dv.append(d)
    fh.close()
    return np.array(xyz), np.array(shell), np.array(dv)


def exclusion_pairs(xyz, cutoff):
    # site pairs closer than the minimum separation; each becomes a mutual-exclusion constraint
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


def solve_unit(dv, pairs, K, neutral):
    n = len(dv)
    m = Model(sense=MINIMIZE)
    m.verbose = 0
    pos = [m.add_var(var_type=BINARY) for _ in range(n)]
    neg = [m.add_var(var_type=BINARY) for _ in range(n)]

    # ddE = sum_i q_i * dv_i in Hartree; minimising it maximises the barrier lowering
    m.objective = xsum(dv[i] * (pos[i] - neg[i]) for i in range(n))

    for i in range(n):
        m += pos[i] + neg[i] <= 1
    m += xsum(pos[i] + neg[i] for i in range(n)) == K
    if neutral:
        m += xsum(pos[i] - neg[i] for i in range(n)) == 0
    for (i, j) in pairs:
        m += pos[i] + neg[i] + pos[j] + neg[j] <= 1

    status = m.optimize(max_seconds=SOLVE_LIMIT)
    if status not in (OptimizationStatus.OPTIMAL, OptimizationStatus.FEASIBLE):
        return None

    sites, charges = [], []
    for i in range(n):
        if pos[i].x >= 0.99:
            sites.append(i)
            charges.append(1)
        elif neg[i].x >= 0.99:
            sites.append(i)
            charges.append(-1)
    return {"ddE": m.objective_value * HARTREE_TO_KCAL, "gap": m.gap,
            "sites": sites, "charges": charges,
            "status": str(status).split(".")[-1]}


def main():
    xyz, shell, dv = load_grid(GRID)
    pairs = exclusion_pairs(xyz, MIN_SEP)
    print("grid %d points, Dv %+.6f to %+.6f Ha; %d exclusion pairs below %.1f A"
          % (len(dv), dv.min(), dv.max(), len(pairs), MIN_SEP))
    print("sweeping K = 2..%d (pocket has %d charged residues)\n" % (K_MAX, N_CHARGED_POCKET))

    results = {}
    rows = []
    budget_base = {}

    for neutral in (False, True):
        label = "neutral" if neutral else "free"
        print("%s" % ("net-neutral (reconciliation control)" if neutral else "net-free (design series)"))
        print("  K   ddE_unit   marginal   per_chg   ddE_budget   %%base   gap    net  shells")
        prev = None
        for K in range(2, K_MAX + 1):
            if neutral and K % 2 == 1:
                rows.append((K, label, "INFEASIBLE_ODD", "", "", "", "", "", "", ""))
                continue
            r = solve_unit(dv, pairs, K, neutral)
            if r is None:
                print("  %-3d %s" % (K, "no solution"))
                rows.append((K, label, "NO_SOLUTION", "", "", "", "", "", "", ""))
                continue

            budget = (Q_BUDGET / K) * r["ddE"]
            if K not in budget_base:
                budget_base.setdefault(label, budget)
            base = budget_base.setdefault(label, budget)
            marginal = "" if prev is None else "%+8.3f" % (r["ddE"] - prev)
            prev = r["ddE"]
            net = int(np.sum(r["charges"]))
            shells = "".join(str(shell[i]) for i in r["sites"])

            print("  %-3d %+9.3f  %8s  %+8.3f  %+10.3f  %5.1f  %.4f  %+d  %s"
                  % (K, r["ddE"], marginal, r["ddE"] / K, budget,
                     100.0 * budget / base, r["gap"], net, shells))

            key = "K%d_%s" % (K, label)
            results[key + "_ddE_unit"] = r["ddE"]
            results[key + "_ddE_budget"] = budget
            results[key + "_sites"] = np.array(r["sites"])
            results[key + "_charges"] = np.array(r["charges"])
            results[key + "_xyz"] = xyz[np.array(r["sites"], dtype=int)]
            results[key + "_gap"] = r["gap"]
            rows.append((K, label, r["status"], "%.4f" % r["ddE"], "%.4f" % (r["ddE"] / K),
                         "%.4f" % budget, "%.1f" % (100.0 * budget / base),
                         "%.4f" % r["gap"], "%+d" % net, shells))
        print("")

    print("reconciliation of the budget curve against the committed net-neutral run:")
    matched = 0
    checked = 0
    for K in sorted(REFERENCE):
        target, tol = REFERENCE[K]
        got = results.get("K%d_neutral_ddE_budget" % K)
        if got is None:
            continue
        checked += 1
        if abs(got - target) <= tol:
            matched += 1
            print("  K=%d  %+8.4f  matches committed %+.3f" % (K, got, target))
        else:
            print("  K=%d  %+8.4f  differs from committed %+.3f by %+.4f"
                  % (K, got, target, got - target))
    if checked and matched == checked:
        print("  the committed sweep was the fixed-budget formulation; new K values are comparable")
    else:
        print("  formulation not confirmed; do not report new K against the committed numbers")

    np.savez(OUT_NPZ, **results)
    fh = open(OUT_TSV, "w")
    fh.write("K\tconstraint\tstatus\tddE_unit\tddE_per_charge\tddE_budget\t"
             "pct_of_K2_budget\tgap\tnet_charge\tshells\n")
    for row in rows:
        fh.write("\t".join(str(v) for v in row) + "\n")
    fh.close()
    print("\nwrote %s and %s" % (OUT_NPZ, OUT_TSV))


if __name__ == "__main__":
    main()
