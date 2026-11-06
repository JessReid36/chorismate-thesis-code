#!/usr/bin/env python3
"""DIAGNOSTIC 31: generate the inputs for the (n, q) envelope scan. Two arms.

A  the designs the optimiser actually produces, unconstrained: n in {1,2,3,5,10,20}, each sphere
   carrying the best charge its own site admits. Every one of those is +-1.00 e, so Arm A IS the
   unit-charge family, and Sigma|q| = n.
B  the magnitude dimension, on the sites that admit the whole ladder at their chosen sign:
   the same n values crossed with q in {0.25, 0.50, 0.75, 1.00}.

Nothing here invents a ranking or a charge. The sites come from the committed cs6_arrange ranking
(per_site_best over each site's own admissible ranges, then greedy_sites with the committed
spacing), and Arm A's charges are per_site_best's own. Validated in D31.2 against the committed
S6-2 and S6-3 arrangements, which this pipeline reproduces exactly.
"""
import argparse, hashlib, importlib.util, shutil, sys
from pathlib import Path

import numpy as np

NS = (1, 2, 3, 5, 10, 20)
LADDER = (0.25, 0.50, 0.75, 1.00)
KEYWORD = "! CAM-B3LYP D3BJ def2-SVPD def2/J RIJCOSX TightSCF"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def load_committed(code):
    """cs_opt_milp and cs6_arrange, by explicit path, from the code repository."""
    code = str(Path(code) / "phase2.2")
    if code not in sys.path: sys.path.insert(0, code)
    import cs_opt_milp as opt
    spec = importlib.util.spec_from_file_location("_arr", Path(code) / "cs6_arrange.py")
    arr = importlib.util.module_from_spec(spec); sys.modules["_arr"] = arr
    try: spec.loader.exec_module(arr)
    except SystemExit: pass          # its __main__ guard runs a CLI; the functions are what we want
    return opt, arr


def in_range(ranges_i, q):
    return any(lo - 1e-9 <= q <= hi + 1e-9 for lo, hi in ranges_i)


def job_inp(charged):
    L = [KEYWORD, "%maxcore 3000", "%pal nprocs 8 end", "%scf", "  MaxIter 300",
         "  STABPerform true", "end"]
    # STABPerform on every run: D30 showed it returns a bit-identical energy where the solution is
    # stable and the correct one where it is not, so it costs nothing and removes a whole failure mode.
    if charged: L += ['%pointcharges "pc.pc"']
    L += ["* xyzfile -2 1 substrate.xyz", ""]
    return "\n".join(L)


def write_pc(path, sites_xyz, q):
    """the committed format, read from a committed pc.pc before this was written:
    a count, then '  %9.6f' and three 10-decimal coordinates separated by six spaces."""
    L = ["%d" % len(q)]
    for qq, (x, y, z) in zip(q, sites_xyz):
        L.append(" %9.6f      %.10f      %.10f      %.10f" % (qq, x, y, z))
    Path(path).write_text("\n".join(L) + "\n")


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--code", required=True); a.add_argument("--stage6", required=True)
    a.add_argument("--geom", required=True, help="a committed control directory pair's parent, for substrate.xyz")
    a.add_argument("--out", required=True)
    a = a.parse_args()
    opt, arr = load_committed(a.code)
    S6, OUT = Path(a.stage6), Path(a.out)
    sites_p, L_p, A_p = S6 / "sites_f42436.tsv", S6 / "L_f42436.tsv", S6 / "A_f42436.tsv"
    print("inputs")
    for p in (sites_p, L_p, A_p): print("   %-22s sha256 %s" % (p.name, sha(p)))
    xyz, layer, ranges, L, A = opt.load(str(sites_p), str(L_p), str(A_p))
    if A.shape[1] != 1: raise SystemExit("A has %d frames; this diagnostic is one frame." % A.shape[1])
    cfg = dict(opt.CONFIG, budget=arr.BUDGET_DEFAULT)
    am, Lm = A.mean(1), L.mean(1)
    best = arr.per_site_best(ranges, Lm, am)
    vals = np.array([b[0] for b in best]); bq = np.array([b[1] for b in best])
    order = list(np.argsort(vals, kind="stable"))

    # Arm B's site pool: admits every ladder value at the sign per_site_best chose
    def admits_ladder(i):
        s = 1.0 if bq[i] > 0 else -1.0
        return all(in_range(ranges[i], s * q) for q in LADDER)
    pool = [i for i in order if admits_ladder(i)]
    print("   Arm B pool: %d of %d sites admit the whole ladder at their chosen sign" % (len(pool), len(xyz)))

    plan, problems = [], []
    selA = {n: arr.greedy_sites(order, xyz, n, cfg["r_min"]) for n in NS}
    selB = {n: arr.greedy_sites(np.array(pool), xyz, n, cfg["r_min"]) for n in NS}
    for sel, nm in ((selA, "A"), (selB, "B")):
        ks = sorted(sel)
        for i in range(len(ks) - 1):
            if sel[ks[i]] != sel[ks[i + 1]][:ks[i]]:
                problems.append("arm %s: the n=%d site set is not a prefix of n=%d" % (nm, ks[i], ks[i + 1]))
    for n in NS:
        s = selA[n]; q = [float(bq[i]) for i in s]
        plan.append(("A/n%02d" % n, s, q))
    for n in NS:
        s = selB[n]
        for qm in LADDER:
            q = [float(np.sign(bq[i]) * qm) for i in s]
            plan.append(("B/n%02d_q%.2f" % (n, qm), s, q))
    # every charge must sit inside its own site's committed range
    for name, s, q in plan:
        for i, qq in zip(s, q):
            if not in_range(ranges[i], qq):
                problems.append("%s: q %+0.4f is outside site %d's admitted ranges %s" % (name, qq, i, ranges[i]))
    if problems:
        for p in problems[:10]: print("   PROBLEM " + p)
        raise SystemExit("refusing to write inputs: %d problem(s)" % len(problems))

    GE = Path(a.geom)
    made = 0
    for name, s, q in plan + [("ctl", [], [])]:
        for st in ("R", "TS"):
            d = OUT / ("%s_%s" % (name, st)); d.mkdir(parents=True, exist_ok=True)
            shutil.copy2(GE / ("ctl_%s" % st) / "substrate.xyz", d / "substrate.xyz")
            if s: write_pc(d / "pc.pc", [xyz[i] for i in s], q)
            (d / "job.inp").write_text(job_inp(bool(s)))
            made += 1
    print("\n%d job directories written" % made)
    print("\n  condition          n   Sigma|q|   max|q|   Sigma q^2     score")
    for name, s, q in plan:
        print("  %-18s %-3d %8.4f %8.4f %10.4f %9.3f"
              % (name, len(s), sum(abs(x) for x in q), max(abs(x) for x in q),
                 sum(x * x for x in q), sum(am[i] * qq + Lm[i] for i, qq in zip(s, q))))
    print("\n  R and TS carry identical charges for every condition: %s"
          % all((OUT / ("%s_R" % n) / "pc.pc").read_bytes() == (OUT / ("%s_TS" % n) / "pc.pc").read_bytes()
                for n, s, _ in plan if s))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
