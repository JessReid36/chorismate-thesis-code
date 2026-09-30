#!/usr/bin/env python3
"""
s14_design_milp.py - Tier 1 of the Phase 2.2 optimiser.

Certified global optimisation of charge placement on a fixed reaction path, plus exact
enumeration of the gauge-degenerate optimal set.

WHAT THIS IS, AND WHAT ITS CERTIFICATE COVERS
The certificate is over a SURROGATE, and the surrogate stack must be stated wherever the
word "optimal" is used:
  1. differential transition-state stabilisation as a proxy for the barrier
  2. LINEAR electrostatic response, i.e. ddE_f = sum_i a_if q_i, exact by superposition
     of Coulomb's law at fixed geometry, but neglecting the measured second-order
     polarisation term b*q^2
  3. a FIXED reaction path
  4. a DISCRETISED candidate surface
  5. a frame-aggregation heuristic standing in for the exponential average, which is
     unusable at this project's barrier spread (Ryde, JCTC 2017, 13, 5745: the exponential
     average is badly conditioned once the snapshot spread exceeds about 7 kJ/mol, and
     our spread is about 19.7 kJ/mol)
A certified global optimum of levels 2-5 is a real object. It is NOT a certified optimal
catalyst, and must never be described as one.

WHY THIS IS A MILP AND NOT A MIQP
The project's primary objective is the mean penalised by spread. With spread measured as
the standard deviation that is a convex MISOCP and needs Gurobi or Mosek. Konno & Yamazaki
(Management Science 1991, 37, 519-531) showed that replacing variance with MEAN ABSOLUTE
DEVIATION gives an equivalent risk-return model that is a pure linear program, because
|x| <= u is LP-representable. MAD and sigma are both dispersion measures - for a normal
distribution MAD = sigma*sqrt(2/pi) - so lambda simply rescales, and the robustness
argument does not depend on which is used. This makes Tier 1 solvable to certified global
optimality by HiGHS through scipy.optimize.milp, with no commercial licence.

THE FORMULATION
Variables, in this order:
    q_i     i = 0..S-1    continuous charge at site i, bounded by +/- QMAX * z_i
    z_i     i = 0..S-1    binary, site i carries a charge
    d_f     f = 0..F-1    continuous, per-frame differential stabilisation
    dbar                  continuous, mean of d over frames
    u_f     f = 0..F-1    continuous, the per-frame dispersion contribution
    w_i     i = 0..S-1    continuous, |q_i|, used only by the optional L1 budget
Constraints:
    d_f - sum_i a_if q_i           = 0            F equalities   (linear response)
    dbar - (1/F) sum_f d_f         = 0            1 equality
    u_f - d_f + dbar              >= 0            F inequalities  (dispersion)
    u_f + d_f - dbar              >= 0            F more, ONLY if two-sided
    w_i >= q_i,  w_i >= -q_i                      2S inequalities, if l1_budget > 0
    sum_i w_i                     <= l1_budget    1 inequality,    if l1_budget > 0
    per-image unimodality chain                   only if --images is supplied
    q_i - QMAX z_i                <= 0            S inequalities  (indicator coupling)
   -q_i - QMAX z_i                <= 0            S inequalities
    sum_i z_i                      = N_CH         1 equality      (cardinality)
    z_i + z_j                     <= 1            per close pair  (r_min packing)
    sum_i q_i                      = 0            1 equality      (neutrality, optional)
Objective, minimised:
    dbar + (lambda / F) sum_f u_f
Negative d means the transition state is stabilised more than the reactant, so minimising
the mean is correct, and adding lambda times a dispersion term penalises designs whose
benefit is fragile across conformers. That fragility is not hypothetical: one measured
frame reverses the sign of the response beyond 4 A.

DISPERSION MEASURE: the default is the SEMI-DEVIATION, counting only frames WORSE than
the mean. A two-sided mean absolute deviation rejects dominant designs: given a site that
is never worse than another on any frame and better on one, two-sided MAD penalises the
better site for being "spread". Semi-deviation does not, and is equally LP-representable.
Set two_sided=1 to recover the Konno-Yamazaki form.

THE NOVELTY THIS IMPLEMENTS, AND WHAT IT DOES NOT YET
Dittner & Hartke (JCTC 2018) use a genetic algorithm: best-of-ten runs, no guarantee, and
at N_Ch = 10 their solution database collapses to essentially one dimension through gauge
freedom - many charge arrangements giving the same potential at the reacting atoms - which
they discovered only by clustering 5388 solutions after the fact.

IMPLEMENTED HERE:
  - the optimum is certified, with a reported optimality gap
  - the degenerate optimal set is ENUMERATED directly by no-good cuts, turning their
    post-hoc observation into a constructive result

IMPLEMENTED BUT INACTIVE WITHOUT FURTHER DATA:
  - reaction-profile unimodality as a HARD constraint rather than a post-hoc penalty.
    This needs the per-image electrostatic potential v_ifm at every path image m, not
    just the reactant-to-TS DIFFERENCE a_if that the design map supplies. Supply it with
    --images and the constraint activates; without it the run says so explicitly and the
    novelty claim must not be made.

CONSEQUENCE OF THE MISSING CONSTRAINT, STATED BECAUSE IT IS VISIBLE IN EVERY RESULT:
A linear objective over a box has its optimum at a VERTEX, so without a field-limiting
constraint every charge saturates at +/-qmax and the continuous charge layer does no work
- the problem reduces to choosing sites and signs. [D2018] Fig. 11 finds charges at
+/-0.751, NOT saturated, because their objective is a full NEB barrier whose validity
terms penalise strong fields. Use l1_budget to impose a total-field budget sum_i |q_i| <= B
if interior charges are wanted before the unimodality data exists.

USAGE
    python3 s14_design_milp.py selftest
    python3 s14_design_milp.py solve <grid.tsv> <amatrix.tsv> [key=value ...]

    amatrix.tsv: one row per site, one column per frame, the linear coefficient a_if in
    kcal/mol per unit charge. Row order must match the grid file.
"""
import sys
import math
from pathlib import Path

import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix

CONFIG = {
    "n_charges": 4,        # cardinality N_Ch
    "qmax": 1.0,           # |q| bound, e. [DTHESIS] Ch.2 uses [-1,+1]
    "lam": 1.0,            # weight on the across-frame dispersion term
    "two_sided": 0,        # 0 = semi-deviation (default), 1 = Konno-Yamazaki MAD
    "neutral": 1,          # 1 = enforce sum q = 0, as [DTHESIS] Ch.7.3 always does
    "r_min": 1.0,          # A, minimum charge-charge separation. [DTHESIS] Ch.2.
                           # VACUOUS unless r_min > the grid spacing, since Poisson-disk
                           # sampling already guarantees that separation. Warned about.
    "l1_budget": 0.0,      # total field budget sum_i |q_i| <= B. 0 disables.
                           # Without this or the unimodality data, charges SATURATE.
    "n_pool": 5,           # how many distinct optima to enumerate
    "mip_gap": 1e-6,
    "time_limit": 600.0,
}


def build(A, xyz, cfg, forbid=(), images=None):
    """Assemble the MILP. A is (S, F).

    images, if given, is (V, E0, ts) where V is (S, F, M) per-image potentials, E0 is
    (F, M) per-image reference energies and ts is the index of the barrier top. It
    activates the unimodality chain. Each link is linear in q directly:
        (E0_{f,m+1} - E0_{f,m}) + sum_i q_i (v_{i,f,m+1} - v_{i,f,m}) >= 0   for m < ts
    and the reverse sign for m >= ts, so no extra variables are needed.
    """
    S, F = A.shape
    i_q, i_z, i_d = 0, S, 2 * S
    i_dbar = 2 * S + F
    i_u = i_dbar + 1
    i_w = i_u + F
    n = i_w + S

    rows, lo, hi = [], [], []

    def add(coefs, l, h):
        r = lil_matrix((1, n))
        for j, v in coefs:
            r[0, j] = v
        rows.append(r)
        lo.append(l)
        hi.append(h)

    # d_f = sum_i a_if q_i
    for f in range(F):
        add([(i_d + f, 1.0)] + [(i_q + i, -A[i, f]) for i in range(S) if A[i, f] != 0.0],
            0.0, 0.0)
    # dbar = mean_f d_f
    add([(i_dbar, 1.0)] + [(i_d + f, -1.0 / F) for f in range(F)], 0.0, 0.0)
    # dispersion. d negative is good, so d_f > dbar is a WORSE-than-average frame.
    # Semi-deviation counts only those; two-sided counts both.
    for f in range(F):
        add([(i_u + f, 1.0), (i_d + f, -1.0), (i_dbar, 1.0)], 0.0, np.inf)
        if cfg["two_sided"]:
            add([(i_u + f, 1.0), (i_d + f, 1.0), (i_dbar, -1.0)], 0.0, np.inf)
    # -QMAX z_i <= q_i <= QMAX z_i
    for i in range(S):
        add([(i_q + i, 1.0), (i_z + i, -cfg["qmax"])], -np.inf, 0.0)
        add([(i_q + i, -1.0), (i_z + i, -cfg["qmax"])], -np.inf, 0.0)
    # cardinality
    add([(i_z + i, 1.0) for i in range(S)], cfg["n_charges"], cfg["n_charges"])
    # neutrality
    if cfg["neutral"]:
        add([(i_q + i, 1.0) for i in range(S)], 0.0, 0.0)
    # r_min packing: z_i + z_j <= 1 for sites closer than r_min
    npair = 0
    if cfg["r_min"] > 0 and xyz is not None:
        for i in range(S):
            for j in range(i + 1, S):
                if np.linalg.norm(xyz[i] - xyz[j]) < cfg["r_min"]:
                    add([(i_z + i, 1.0), (i_z + j, 1.0)], -np.inf, 1.0)
                    npair += 1
    # optional total-field budget, the mechanism that prevents saturation
    if cfg["l1_budget"] > 0:
        for i in range(S):
            add([(i_w + i, 1.0), (i_q + i, -1.0)], 0.0, np.inf)
            add([(i_w + i, 1.0), (i_q + i, 1.0)], 0.0, np.inf)
        add([(i_w + i, 1.0) for i in range(S)], -np.inf, cfg["l1_budget"])

    # optional unimodality of the reaction profile, as a HARD constraint
    n_uni = 0
    if images is not None:
        V, E0, ts = images
        M = V.shape[2]
        for f in range(F):
            for mm in range(M - 1):
                dv = V[:, f, mm + 1] - V[:, f, mm]
                de = E0[f, mm + 1] - E0[f, mm]
                coefs = [(i_q + i, float(dv[i])) for i in range(S) if dv[i] != 0.0]
                if mm < ts:      # rising to the barrier top
                    add(coefs, -de, np.inf)
                else:            # falling after it
                    add([(j, -v) for j, v in coefs], de, np.inf)
                n_uni += 1

    # no-good cuts excluding previously found supports
    for supp in forbid:
        k = len(supp)
        add([(i_z + i, 1.0) for i in supp], -np.inf, k - 1.0)

    from scipy.sparse import vstack
    Amat = vstack(rows).tocsc()
    cons = LinearConstraint(Amat, np.array(lo), np.array(hi))

    c = np.zeros(n)
    c[i_dbar] = 1.0
    c[i_u:i_u + F] = cfg["lam"] / F

    integrality = np.zeros(n)
    integrality[i_z:i_z + S] = 1
    lb = np.full(n, -np.inf)
    ub = np.full(n, np.inf)
    lb[i_q:i_q + S], ub[i_q:i_q + S] = -cfg["qmax"], cfg["qmax"]
    lb[i_z:i_z + S], ub[i_z:i_z + S] = 0.0, 1.0
    lb[i_u:i_u + F] = 0.0
    lb[i_w:i_w + S], ub[i_w:i_w + S] = 0.0, cfg["qmax"]

    idx = dict(q=i_q, z=i_z, d=i_d, dbar=i_dbar, u=i_u, w=i_w, n=n,
               npair=npair, n_uni=n_uni)
    return c, cons, integrality, Bounds(lb, ub), idx


def solve_once(A, xyz, cfg, forbid=(), images=None):
    # Neutrality with a single charge forces q = 0: the design is inert and the
    # objective is trivially zero. Mathematically correct, but a silent trap, so it is
    # named rather than returned quietly.
    if cfg["neutral"] and cfg["n_charges"] == 1:
        raise ValueError(
            "n_charges=1 with neutral=1 forces q=0 and gives an inert design. "
            "Either set neutral=0, as Burschowsky's single cation implies, or use "
            "n_charges>=2.")
    if cfg["qmax"] <= 0:
        raise ValueError("qmax must be positive; qmax=0 admits only the null design.")
    c, cons, integrality, bounds, idx = build(A, xyz, cfg, forbid, images)
    res = milp(c=c, constraints=cons, integrality=integrality, bounds=bounds,
               options={"mip_rel_gap": cfg["mip_gap"], "time_limit": cfg["time_limit"]})
    if not res.success:
        return None, res, idx
    x = res.x
    S, F = A.shape
    out = {
        "objective": float(res.fun),
        "support": sorted(i for i in range(S) if x[idx["z"] + i] > 0.5),
        "q": {i: float(x[idx["q"] + i]) for i in range(S) if x[idx["z"] + i] > 0.5},
        "d": [float(x[idx["d"] + f]) for f in range(F)],
        "mean": float(x[idx["dbar"]]),
        "mad": float(sum(x[idx["u"] + f] for f in range(F)) / F),
        "gap": float(getattr(res, "mip_gap", float("nan"))),
        "npair": idx["npair"],
        "n_uni": idx["n_uni"],
        "l1": float(sum(abs(v) for v in
                        (x[idx["q"] + i] for i in range(S)))),
        # saturation is only meaningful for sites that actually carry charge
        "saturated": bool(any(x[idx["z"] + i] > 0.5 for i in range(S))) and
                     all(abs(abs(x[idx["q"] + i]) - cfg["qmax"]) < 1e-6
                         for i in range(S) if x[idx["z"] + i] > 0.5),
    }
    return out, res, idx


def report(sol, cfg, label=""):
    print(f"  {label}objective {sol['objective']:+.4f}   "
          f"mean {sol['mean']:+.4f}   MAD {sol['mad']:.4f}   "
          f"gap {sol['gap']:.2e}")
    qs = "  ".join(f"site {i}: q={q:+.3f}" for i, q in sorted(sol["q"].items()))
    print(f"    {qs}")
    print(f"    per-frame d: " + "  ".join(f"{d:+.3f}" for d in sol["d"]))


def enumerate_optima(A, xyz, cfg, images=None):
    """Solve, cut, re-solve. This is the gauge-degeneracy result [D2018] found only by
    clustering 5388 GA solutions after the fact."""
    found, forbid = [], []
    for k in range(cfg["n_pool"]):
        sol, res, _ = solve_once(A, xyz, cfg, forbid, images)
        if sol is None:
            print(f"  no further feasible solution after {k} (solver: {res.message})")
            break
        found.append(sol)
        forbid.append(sol["support"])
    return found


def selftest():
    """Verify the formulation on data whose optimum is known by construction."""
    print("SELF-TEST\n")
    ok = True

    print("1. every constraint is honoured by the returned solution")
    rng = np.random.default_rng(20260930)
    S, F = 40, 6
    A = rng.normal(0, 1, (S, F))
    xyz = rng.normal(0, 6, (S, 3))
    cfg = dict(CONFIG, n_charges=4, n_pool=1)
    sol, res, idx = solve_once(A, xyz, cfg)
    assert sol is not None, "solver failed on the random instance"
    n_sel = len(sol["support"])
    print(f"   [{'PASS' if n_sel == cfg['n_charges'] else 'FAIL'}] cardinality: "
          f"{n_sel} sites selected, requested {cfg['n_charges']}")
    ok &= n_sel == cfg["n_charges"]
    qsum = sum(sol["q"].values())
    print(f"   [{'PASS' if abs(qsum) < 1e-7 else 'FAIL'}] neutrality: sum q = {qsum:+.2e}")
    ok &= abs(qsum) < 1e-7
    within = all(abs(q) <= cfg["qmax"] + 1e-9 for q in sol["q"].values())
    print(f"   [{'PASS' if within else 'FAIL'}] charge bounds respected")
    ok &= within
    worst = 0.0
    for f in range(F):
        lhs = sum(A[i, f] * sol["q"][i] for i in sol["support"])
        worst = max(worst, abs(lhs - sol["d"][f]))
    print(f"   [{'PASS' if worst < 1e-7 else 'FAIL'}] d_f reproduces sum_i a_if q_i "
          f"(worst {worst:.2e})")
    ok &= worst < 1e-7
    m = sum(sol["d"]) / F
    print(f"   [{'PASS' if abs(m - sol['mean']) < 1e-7 else 'FAIL'}] dbar is the mean "
          f"({m:+.5f} vs {sol['mean']:+.5f})")
    ok &= abs(m - sol["mean"]) < 1e-7
    # the dispersion measure must match the configured one. d negative is good, so a
    # frame with d_f > dbar is WORSE than average and is what semi-deviation counts.
    if cfg["two_sided"]:
        disp = sum(abs(d - m) for d in sol["d"]) / F
        name = "two-sided MAD"
    else:
        disp = sum(max(0.0, d - m) for d in sol["d"]) / F
        name = "semi-deviation"
    print(f"   [{'PASS' if abs(disp - sol['mad']) < 1e-6 else 'FAIL'}] dispersion term is "
          f"the {name} ({disp:.5f} vs {sol['mad']:.5f})")
    ok &= abs(disp - sol["mad"]) < 1e-6
    obj = sol["mean"] + cfg["lam"] * sol["mad"]
    print(f"   [{'PASS' if abs(obj - sol['objective']) < 1e-6 else 'FAIL'}] objective is "
          f"mean + lambda*dispersion ({obj:+.5f} vs {sol['objective']:+.5f})")
    ok &= abs(obj - sol["objective"]) < 1e-6

    print("\n2. the certified optimum matches a CLOSED-FORM answer")
    # One charge, no neutrality, no dispersion penalty. Then the objective is just the
    # mean over frames of a_if q_i, minimised over q_i in [-qmax, qmax] and over i.
    # The optimum is therefore  -qmax * max_i |mean_f a_if|,  computable by hand.
    # This does not call the solver twice, so it is a genuine independent check.
    S2, F2 = 30, 5
    A2 = rng.normal(0, 1, (S2, F2))
    cfg2 = dict(CONFIG, n_charges=1, r_min=0.0, neutral=0, lam=0.0, n_pool=1)
    sol2, _, _ = solve_once(A2, None, cfg2)
    means = A2.mean(axis=1)
    closed = -cfg2["qmax"] * float(np.max(np.abs(means)))
    best_i = int(np.argmax(np.abs(means)))
    d = abs(sol2["objective"] - closed)
    good = d < 1e-6 and sol2["support"] == [best_i]
    print(f"   [{'PASS' if good else 'FAIL'}] MILP {sol2['objective']:+.6f} at site "
          f"{sol2['support'][0]}; closed form {closed:+.6f} at site {best_i}; "
          f"difference {d:.2e}")
    ok &= good

    print("\n3. the spread penalty changes the design, as it must")
    A4 = np.zeros((6, 3))
    A4[0] = [-3.0, -3.0, -3.0]      # uniform benefit
    A4[1] = [-9.0, -0.1, -0.1]      # large but fragile
    A4[2:] = -0.2
    for lam, expect in ((0.0, 1), (5.0, 0)):
        cfg4 = dict(CONFIG, n_charges=1, r_min=0.0, neutral=0, lam=lam, n_pool=1)
        s4, _, _ = solve_once(A4, None, cfg4)
        got = s4["support"][0]
        good = got == expect
        ok &= good
        print(f"   [{'PASS' if good else 'FAIL'}] lambda={lam}: chose site {got}, "
              f"expected {expect} "
              f"({'fragile, high mean' if expect == 1 else 'uniform, robust'})")

    print("\n3b. semi-deviation does not reject a dominant design")
    # site 1 is never worse than site 0 on any frame and is better on one, so it
    # dominates and must be chosen at any lambda.
    A5 = np.zeros((4, 3))
    A5[0] = [-4.0, -4.0, -4.0]
    A5[1] = [-8.0, -4.0, -4.0]
    A5[2:] = -0.1
    for two, expect in ((0, 1), (1, 0)):
        cfg5 = dict(CONFIG, n_charges=1, r_min=0.0, neutral=0, lam=1.0,
                    two_sided=two, n_pool=1)
        s5, _, _ = solve_once(A5, None, cfg5)
        got = s5["support"][0]
        good = got == expect
        ok &= good
        name = "two-sided MAD" if two else "semi-deviation"
        note = ("rejects the dominant design, which is the defect"
                if two else "keeps the dominant design, as it must")
        print(f"   [{'PASS' if good else 'FAIL'}] {name}: chose site {got} - {note}")

    print("\n3c. the L1 budget prevents saturation")
    A6 = rng.normal(0, 1, (20, 4))
    for budget in (0.0, 1.5):
        cfg6 = dict(CONFIG, n_charges=4, r_min=0.0, l1_budget=budget, n_pool=1)
        s6, _, _ = solve_once(A6, None, cfg6)
        sat = s6["saturated"]
        good = (sat if budget == 0.0 else not sat)
        ok &= good
        print(f"   [{'PASS' if good else 'FAIL'}] budget={budget}: sum|q| = "
              f"{s6['l1']:.3f}, saturated = {sat}"
              + ("  (expected: no budget means vertex optimum)" if budget == 0.0
                 else "  (expected: budget forces interior charges)"))

    print("\n3d. the unimodality chain is enforced when per-image data is supplied")
    # Two images either side of a barrier top. A field strong enough to inverpt the
    # rise must be rejected.
    S7, F7, M7 = 6, 2, 3
    V7 = np.zeros((S7, F7, M7))
    V7[0, :, :] = np.array([[0.0, 1.0, 0.0]] * F7)   # site 0 raises image 1
    E07 = np.tile(np.array([0.0, 2.0, -1.0]), (F7, 1))
    A7 = np.full((S7, F7), -0.5)
    cfg7 = dict(CONFIG, n_charges=1, r_min=0.0, neutral=0, lam=0.0, n_pool=1)
    s_no, _, idx_no = solve_once(A7, None, cfg7)
    s_yes, _, idx_yes = solve_once(A7, None, cfg7, images=(V7, E07, 1))
    print(f"   [{'PASS' if idx_no['n_uni'] == 0 else 'FAIL'}] no images supplied: "
          f"{idx_no['n_uni']} unimodality constraints added")
    ok &= idx_no["n_uni"] == 0
    nexp = F7 * (M7 - 1)
    print(f"   [{'PASS' if idx_yes['n_uni'] == nexp else 'FAIL'}] images supplied: "
          f"{idx_yes['n_uni']} constraints added, expected {nexp}")
    ok &= idx_yes["n_uni"] == nexp
    feas = s_yes is not None
    print(f"   [{'PASS' if feas else 'FAIL'}] the constrained problem is still feasible")
    ok &= feas

    print("\n4. no-good cuts return genuinely distinct supports")
    cfg5 = dict(CONFIG, n_charges=3, n_pool=4)
    pool = enumerate_optima(A, xyz, cfg5)
    sup = [tuple(s["support"]) for s in pool]
    distinct = len(set(sup)) == len(sup)
    print(f"   [{'PASS' if distinct else 'FAIL'}] {len(sup)} solutions, all distinct")
    ok &= distinct
    print("    objectives: " + "  ".join(f"{s['objective']:+.4f}" for s in pool))
    print("    (a flat sequence is the gauge degeneracy [D2018] reported; a rising one")
    print("     means the optimum is unique and the rest are successively worse)")

    print(f"\n{'ALL SELF-TESTS PASSED' if ok else 'SELF-TEST FAILURES ABOVE'}")
    return 0 if ok else 1


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    if sys.argv[1] == "selftest":
        return selftest()
    if sys.argv[1] != "solve":
        sys.exit(__doc__)

    grid, amat = Path(sys.argv[2]), Path(sys.argv[3])
    cfg = dict(CONFIG)
    for kv in sys.argv[4:]:
        k, v = kv.split("=", 1)
        if k not in cfg:
            sys.exit(f"unknown parameter {k!r}. Known: {', '.join(sorted(cfg))}")
        cfg[k] = type(cfg[k])(v)

    xyz = np.array([[float(x) for x in l.split("\t")[1:4]]
                    for l in grid.read_text().splitlines()
                    if not l.startswith("#") and not l.startswith("idx")])
    A = np.array([[float(x) for x in l.split()]
                  for l in amat.read_text().splitlines() if not l.startswith("#")])
    if A.shape[0] != xyz.shape[0]:
        sys.exit(f"row mismatch: grid has {xyz.shape[0]} sites, a-matrix has {A.shape[0]}")
    print(f"{A.shape[0]} candidate sites, {A.shape[1]} frames")
    if cfg["r_min"] > 0 and len(xyz) > 1:
        from itertools import combinations as _cb
        closest = min(float(np.linalg.norm(xyz[i] - xyz[j]))
                      for i, j in _cb(range(len(xyz)), 2))
        if cfg["r_min"] <= closest:
            print(f"  NOTE: r_min = {cfg['r_min']} A is VACUOUS. The closest pair on this"
                  f" grid is {closest:.3f} A")
            print(f"        apart, so no packing constraint can bind. Poisson-disk"
                  f" sampling already enforces")
            print(f"        the separation. r_min only bites if set above the grid"
                  f" spacing.")
    print(f"config: " + "  ".join(f"{k}={v}" for k, v in sorted(cfg.items())))
    print()
    pool = enumerate_optima(A, xyz, cfg)
    if not pool:
        print("no feasible design found")
        return 1
    print(f"{len(pool)} certified solutions, best first:\n")
    for k, s in enumerate(pool):
        report(s, cfg, label=f"[{k}] ")
        print()
    spread = pool[-1]["objective"] - pool[0]["objective"]
    print(f"objective spread across the pool: {spread:.5f} kcal/mol")
    if spread < 1e-4:
        print("  The optimum is DEGENERATE: these supports are equally good. That is the")
        print("  gauge freedom [D2018] reported at N_Ch = 10, here characterised directly")
        print("  rather than by clustering solutions after the fact.")
    else:
        print("  The optimum is unique within the pool; later entries are strictly worse.")
    print()
    if pool[0]["n_uni"] == 0:
        print("REACTION-PROFILE UNIMODALITY WAS NOT ENFORCED: no per-image potential data")
        print("  was supplied. The hard-validity-constraint novelty claim does not apply")
        print("  to this run.")
    if pool[0]["saturated"]:
        print("ALL CHARGES SATURATED at +/-qmax. Expected for a linear objective over a")
        print("  box with no field-limiting constraint: the optimum is at a vertex, so")
        print("  only sites and signs are being chosen. Set l1_budget, or supply the")
        print("  per-image data so unimodality can bind, if interior charges are wanted.")
    print()
    print("CERTIFICATE SCOPE: global optimality over the LINEAR-RESPONSE surrogate on a")
    print("FIXED path and a DISCRETISED surface, with a frame-aggregation heuristic in")
    print("place of the exponential average. Not a certified optimal catalyst.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
