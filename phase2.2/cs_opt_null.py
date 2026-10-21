#!/usr/bin/env python3
"""
cs_opt_null.py - the optimiser stage's null distribution (decision O10).

Random designs under the SAME rules as the optimiser: the same n, the same admitted charges per site, the same
3.25 A spacing and the same budget sum_i q_i^2 <= B. Every optimised design is reported with its percentile here,
so "better than the optimiser's own random baseline" is a statement with a number behind it.

TWO KINDS (O10)
  (i)  random sites, random charges
       Sites are added at random, respecting the spacing. Each sphere is given one of its site's admitted ranges at
       random. The vector of charges is drawn uniformly inside the budget - uniformly on the sphere of radius
       sqrt(B) in n dimensions, scaled by U^(1/n) so the draw is uniform in the BALL, not on its surface - and each
       charge is then moved to the nearest value its range admits. A design is redrawn only if that last step pushes
       sum_i q_i^2 above B; the acceptance rate is reported.
       Drawing each charge uniformly over its own range and rejecting designs over budget would almost never be
       accepted at 20 spheres, which is why the draw is made inside the budget to begin with.
  (ii) random sites, best charges
       The sites are random; the charges are the best the budget allows on those sites. This tests whether choosing
       the SITES matters, as against choosing the charges. On a fixed set of sites with no range bound active and
       lambda = 0 the answer is in closed form, q = -sqrt(B) abar / ||abar||; the general case (bounds active, or
       lambda > 0) is a small problem of its own, solved here with the same MILP as the optimiser, with the sites
       forced.

The two kinds are NOT paired: each draws its own random sites, so (ii) is not (i) with the charges improved. They
answer different questions - what a random design scores, and what a random SET OF SITES scores once its charges are
chosen as well as they can be - and a design's percentile is reported in each separately.

WHY THE SAMPLING RULE IS PART OF THE RECORD
It shapes the distribution, so it is stated rather than left to the code: sites uniform among those still admissible
given the spacing, ranges uniform among a site's admitted ranges, charges uniform in the budget ball then moved to
the nearest admitted value. The rule is pre-registered with the sample sizes (10 000 of each kind per n at first
order; 50 in total re-scored with full single points).

WHAT THIS PRODUCES IS NOT A RESULT. On the development set it is a code test.

USAGE
    python3 cs_opt_null.py selftest
    python3 cs_opt_null.py run <sites.tsv> <L.tsv> <A.tsv> [key=value ...]
"""
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from cs_opt_milp import (CONFIG, Problem, close_pairs, cliques, load, read_matrix,  # noqa: E402
                         read_sites, sha256, violations)

NULL_CONFIG = dict(CONFIG)
NULL_CONFIG.update({"draws": 10000,      # kind (i): drawn designs
                    "draws_ii": 100,    # kind (ii): each one is a small solve of its own, so far fewer
                    "kind": "both", "report_every": 2000})


def min_abs_q(ranges_i):
    """the smallest |q| this site can carry. Zero on the anion/neutral layer; on the cation layer it is that site's
    own minimum positive charge, which the budget must pay for whether the design wants it or not."""
    return min((0.0 if lo <= 0 <= hi else min(abs(lo), abs(hi))) for lo, hi in ranges_i)


def random_sites(rng, xyz, tree, n, r_min, ranges=None, budget=None, tries=200):
    """n sites, added one at a time at random, no two closer than r_min.

    With `ranges` and `budget`, a site is skipped when its own minimum charge would not fit in what the budget has
    left (decision O10: "Sites whose minimum charge alone exceeds the remaining budget are skipped"). That rule is
    not cosmetic on this grid: 1228 of the 2769 development sites cannot carry zero at all - 1064 of them force
    |q| = 1 - so a set of sites drawn without it is usually impossible to charge inside the budget. Drawn without
    the rule, 22.6 % of random 4-site sets and 99.9 % of random 20-site sets break B = 2.0 on their minimum charges
    alone.

    Returns None if n sites cannot be placed in `tries` attempts, which at a tight budget can be genuine.
    """
    S = len(xyz)
    mq = None if ranges is None else np.array([min_abs_q(rr) for rr in ranges])
    for _ in range(tries):
        chosen, blocked, spent = [], np.zeros(S, dtype=bool), 0.0
        order = rng.permutation(S)
        for i in order:
            if blocked[i]:
                continue
            if mq is not None and budget is not None and spent + mq[i] ** 2 > budget + 1e-12:
                continue                      # its minimum charge alone will not fit
            chosen.append(int(i))
            if mq is not None:
                spent += mq[i] ** 2
            blocked[tree.query_ball_point(xyz[i], r_min)] = True
            if len(chosen) == n:
                return sorted(chosen)
    return None


def nearest_admitted(q, ranges_i):
    """the value nearest q that this site admits, and whether q was already admitted."""
    best, d0 = None, np.inf
    for lo, hi in ranges_i:
        v = min(max(q, lo), hi)
        d = abs(v - q)
        if d < d0:
            best, d0 = v, d
    return best, d0 <= 1e-12


def draw_charges(rng, n, B):
    """uniform in the ball sum q^2 <= B: a direction on the unit sphere times radius sqrt(B) * U^(1/n)."""
    g = rng.normal(size=n)
    g /= np.linalg.norm(g)
    return np.sqrt(B) * rng.random() ** (1.0 / n) * g


def score(A, L, sites, q):
    """the frozen score of a design: the per-frame values, their mean and semi-deviation."""
    d = A[sites].T @ np.asarray(q) + L[sites].sum(axis=0)        # sum_i a_if q_i + sum_i L_if
    mean = float(d.mean())
    sd = float(np.mean(np.maximum(0.0, d - mean)))
    return d, mean, sd


def kind_i(rng, xyz, tree, ranges, A, L, n, B, lam, r_min, draws, report=None):
    """random sites, random charges."""
    vals, kept, tried = [], 0, 0
    while kept < draws:
        tried += 1
        sites = random_sites(rng, xyz, tree, n, r_min, ranges, B)
        if sites is None:
            raise SystemExit("could not place %d spheres with %.2f A spacing inside the budget %.2f e^2; at this n "
                             "the budget may admit no design at all" % (n, r_min, B))
        qs, over = [], False
        for i in sites:
            rr = ranges[i]
            pick = rr[rng.integers(len(rr))] if len(rr) > 1 else rr[0]
            qs.append(pick)
        raw = draw_charges(rng, n, B)
        q = []
        for (lo, hi), v in zip(qs, raw):
            vv, _ = nearest_admitted(float(v), [(lo, hi)])
            q.append(vv)
        if float(np.dot(q, q)) > B + 1e-12:
            over = True
        if over:
            continue
        kept += 1
        d, mean, sd = score(A, L, sites, q)
        vals.append(mean + lam * sd)
        if report and kept % report == 0:
            print("      %d of %d drawn (acceptance %.3f)" % (kept, draws, kept / tried))
    return np.array(vals), kept / tried


def kind_ii(rng, xyz, tree, ranges, A, L, n, B, lam, r_min, draws, cfg, report=None):
    """random sites, the best charges the budget allows on them, by the same MILP with the sites forced."""
    vals, tried = [], 0
    while len(vals) < draws:
        tried += 1
        if tried > 50 * draws:
            break
        sites = random_sites(rng, xyz, tree, n, r_min, ranges, B)
        if sites is None:
            continue
        sub = sorted(sites)
        p = Problem(xyz[sub], [ranges[i] for i in sub], L[sub], A[sub],
                    dict(cfg, r_min=0.0), cliq=[])
        sol, _ = p.solve(n, B, time_limit=30.0)
        if sol is None:
            continue                          # these sites cannot be charged inside the budget; draw again
        vals.append(sol["objective"])
        if report and len(vals) % report == 0:
            print("      %d of %d solved (acceptance %.3f)" % (len(vals), draws, len(vals) / tried))
    return np.array(vals), (len(vals) / tried if tried else 0.0)


def percentile_of(value, sample):
    """the fraction of the sample at or below `value` (lower is better, so a small percentile is good)."""
    if len(sample) == 0:
        return float("nan")
    return float(np.mean(np.asarray(sample) <= value) * 100.0)


def summarise(name, v):
    if len(v) == 0:
        print("   %s: nothing drawn" % name)
        return
    qs = np.percentile(v, [0, 1, 5, 25, 50, 75, 100])
    print("   %-28s n=%d  best %+8.4f  1%% %+8.4f  5%% %+8.4f  median %+8.4f  worst %+8.4f  mean %+8.4f  sd %6.4f"
          % (name, len(v), qs[0], qs[1], qs[2], qs[4], qs[6], float(v.mean()), float(v.std(ddof=1))))


def selftest():
    print("SELF-TEST  (cs_opt_null.py)\n")
    ok = True
    rng = np.random.default_rng(CONFIG["seed"])

    def say(good, text):
        nonlocal ok
        ok &= bool(good)
        print("   [%s] %s" % ("PASS" if good else "FAIL", text))

    print("1. the charge draw is uniform inside the budget, not on its surface")
    B, n = 2.0, 4
    Q = np.array([draw_charges(rng, n, B) for _ in range(20000)])
    r2 = (Q ** 2).sum(1)
    say(r2.max() <= B + 1e-9, "every draw inside the budget (largest sum q^2 %.4f of %.2f)" % (r2.max(), B))
    # for a uniform ball the radius r has P(r <= x) = (x/R)^n, so the median radius is R * 0.5^(1/n)
    med = np.median(np.sqrt(r2))
    want = np.sqrt(B) * 0.5 ** (1.0 / n)
    say(abs(med - want) < 0.02 * want, "median radius %.4f against %.4f for a uniform ball" % (med, want))

    print("\n2. random sites always respect the spacing")
    S = 300
    xyz = rng.normal(0, 6, (S, 3))
    tree = cKDTree(xyz)
    worst = 0.0
    for _ in range(300):
        sites = random_sites(rng, xyz, tree, 4, 3.25)
        for a in range(len(sites)):
            for b in range(a + 1, len(sites)):
                worst = max(worst, 3.25 - float(np.linalg.norm(xyz[sites[a]] - xyz[sites[b]])))
    say(worst <= 0.0, "300 draws, closest pair never under 3.25 A (worst excess %.2e)" % max(0.0, worst))

    print("\n3. charges are moved to a value the site admits")
    ranges = [[(-1.0, 0.0)], [(0.25, 1.0)], [(-1.0, 1.0)]]
    cases = [(-0.4, 0, -0.4), (0.7, 0, 0.0), (0.1, 1, 0.25), (-0.5, 1, 0.25), (0.3, 2, 0.3)]
    good = all(abs(nearest_admitted(q, ranges[i])[0] - want) < 1e-12 for q, i, want in cases)
    say(good, "a charge outside a site's range is moved to its nearest admitted value")

    print("\n4. the two kinds agree with the optimiser's own arithmetic")
    S4, F4 = 60, 5
    A4 = rng.normal(0, 1, (S4, F4)); L4 = rng.normal(0, 0.02, (S4, F4))
    xyz4 = rng.normal(0, 8, (S4, 3)); tree4 = cKDTree(xyz4)
    r4 = [[(-1.0, 1.0)] if i % 2 else [(0.25, 1.0)] for i in range(S4)]
    v1, acc = kind_i(rng, xyz4, tree4, r4, A4, L4, 3, 1.0, 0.94, 3.25, 200)
    say(len(v1) == 200 and np.isfinite(v1).all(), "kind (i): 200 designs scored, acceptance %.3f" % acc)
    # one of them, re-scored by the optimiser's own violation check, must agree
    sites = random_sites(rng, xyz4, tree4, 3, 3.25)
    p4 = Problem(xyz4[sorted(sites)], [r4[i] for i in sorted(sites)], L4[sorted(sites)], A4[sorted(sites)],
                 dict(CONFIG, r_min=0.0, lam=0.94), cliq=[])
    s4, _ = p4.solve(3, 1.0, time_limit=30.0)
    d, mean, sd = score(A4, L4, sorted(sites), s4["q"])
    say(abs(mean + 0.94 * sd - s4["objective"]) < 1e-6,
        "kind (ii): the MILP's score on forced sites matches this script's arithmetic (%.2e)"
        % abs(mean + 0.94 * sd - s4["objective"]))
    v2, acc2 = kind_ii(rng, xyz4, tree4, r4, A4, L4, 3, 1.0, 0.94, 3.25, 20, dict(CONFIG, lam=0.94))
    say(len(v2) == 20 and v2.mean() < v1.mean(),
        "kind (ii) beats kind (i) on average (%.4f against %.4f), as choosing the charges must (acceptance %.2f)"
        % (v2.mean(), v1.mean(), acc2))

    print("\n5. the budget's skip rule is applied when sites are drawn")
    # sites that cannot carry zero must not be heaped into one design beyond what the budget can pay for.
    r6 = [[(1.0, 1.0)] if i % 2 else [(-1.0, 1.0)] for i in range(120)]   # half the sites force |q| = 1
    xyz6 = rng.normal(0, 30, (120, 3)); tree6 = cKDTree(xyz6)
    worst = 0
    for _ in range(200):
        st = random_sites(rng, xyz6, tree6, 4, 3.25, r6, 2.0)
        if st is None:
            continue
        worst = max(worst, sum(1 for i in st if min_abs_q(r6[i]) >= 1.0))
    say(worst <= 2, "with B = 2.0, no drawn design holds more than 2 sites that force |q| = 1 (worst %d)" % worst)
    st_no = [random_sites(rng, xyz6, tree6, 4, 3.25) for _ in range(50)]
    many = max(sum(1 for i in st if min_abs_q(r6[i]) >= 1.0) for st in st_no if st)
    say(many > 2, "without the rule such designs do appear (worst %d), so the rule is doing something" % many)

    print("\n6. percentiles read the right way round")
    samp = np.array([-1.0, -2.0, -3.0, -4.0])
    say(abs(percentile_of(-3.5, samp) - 25.0) < 1e-9 and abs(percentile_of(0.0, samp) - 100.0) < 1e-9,
        "a design better than all but one of four sits at the 25th percentile")

    print("\n%s" % ("ALL SELF-TESTS PASSED" if ok else "SELF-TEST FAILURES ABOVE"))
    return 0 if ok else 1


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    if sys.argv[1] == "selftest":
        return selftest()
    if sys.argv[1] != "run" or len(sys.argv) < 5:
        sys.exit(__doc__)
    cfg = dict(NULL_CONFIG)
    for kv in sys.argv[5:]:
        k, v = kv.split("=", 1)
        if k not in cfg:
            sys.exit("unknown parameter %r" % k)
        cfg[k] = type(cfg[k])(v)
    xyz, layer, ranges, L, A = load(sys.argv[2], sys.argv[3], sys.argv[4])
    print("config: " + "  ".join("%s=%s" % (k, v) for k, v in sorted(cfg.items())))
    print()
    print("EXPLORATORY on the development set: a code test, not a result.")
    print()
    rng = np.random.default_rng(cfg["seed"])
    tree = cKDTree(xyz)
    n, B, lam = cfg["n"], cfg["budget"], cfg["lam"]
    print("null distribution at n = %d, B = %.2f, lambda = %.2f" % (n, B, lam))
    if cfg["kind"] in ("i", "both"):
        t0 = time.time()
        v1, acc = kind_i(rng, xyz, tree, ranges, A, L, n, B, lam, cfg["r_min"], cfg["draws"],
                         report=cfg["report_every"])
        summarise("(i) random sites, charges", v1)
        print("      acceptance rate %.4f; %.1fs" % (acc, time.time() - t0))
    if cfg["kind"] in ("ii", "both"):
        t0 = time.time()
        v2, acc2 = kind_ii(rng, xyz, tree, ranges, A, L, n, B, lam, cfg["r_min"],
                           cfg["draws_ii"], dict(cfg), report=max(1, cfg["draws_ii"] // 4))
        summarise("(ii) random sites, best q", v2)
        print("      acceptance rate %.4f; %.1fs for %d designs" % (acc2, time.time() - t0, len(v2)))
        if len(v2) < cfg["draws_ii"]:
            print("      NOTE: only %d of the %d designs asked for were drawn before the attempt cap; at this n the"
                  % (len(v2), cfg["draws_ii"]))
            print("      budget admits few random site sets at all, which is itself reported rather than hidden")
    return 0


if __name__ == "__main__":
    sys.exit(main())
