#!/usr/bin/env python3
"""
cs_opt_milp.py - the optimiser stage's MILP, on the two-layer grid.

The decisions this implements are CHARGE_MODEL_DECISIONS.md, "Update: the optimiser stage decided
(20 October 2026)", O1-O10, with the reasoning in phase2.2/OPTIMISER_DECISIONS_20261020.md. s14_design_milp.py is
the reference this replaces; what is new here is L, two admitted charge ranges per site, the spacing as cliques,
the sum-of-squares budget, and distinct-design enumeration.

WHAT IT OPTIMISES (O1, O5)
Per frame f the frozen score of a design is
    d_f = sum_i a_if q_i + sum_i L_if z_i
and the design's score is
    dbar + (lambda / F) sum_f u_f,    u_f >= d_f - dbar,  u_f >= 0
i.e. the mean over frames plus lambda times the semi-deviation, which counts only frames worse than the mean.
d negative means the transition state is stabilised more than the reactant, so the score is MINIMISED. Polarisation
is not in it (O1): the pool is re-scored afterwards with full single points.

WHAT IT MAY CHOOSE
Each site admits one or two charge ranges (Stage 5, with Amendment 2: q <= 0 on the anion/neutral layer only). Per
site i and admitted range r:
    y_ir  binary, "a sphere sits here with its charge in range r"; at most one r per site
    q_ir  continuous, lo_r * y_ir <= q_ir <= hi_r * y_ir, so q_ir = 0 unless y_ir = 1
    z_i = sum_r y_ir   (is there a sphere)        q_i = sum_r q_ir   (its charge)

THE RULES
    sum_i z_i = n                                     exactly n spheres (O6)
    sum_i s_i <= B                                    the budget (O3), s_i standing in for q_i^2
    sum_{i in K} z_i <= 1   for each clique K         the 3.25 A spacing
    one cut per earlier design                        distinct designs (O7)

THE BUDGET'S STAND-IN (O3)
HiGHS takes no quadratic rule, so q_i^2 is replaced by the upper envelope of 21 tangent lines
    s_i >= 2c q_i - c^2,   c = -1.0, -0.9, ..., +1.0
Tangents lie below q^2, so a design's true sum_i q_i^2 may exceed B by at most (0.1/2)^2 = 0.0025 e^2 per sphere
(0.05 e^2 at n = 20). The true value is computed and reported for every design.

Each tangent is written with its constant tied to the site's indicator, s_i >= 2c q_i - c^2 z_i, rather than
s_i >= 2c q_i - c^2. For z_i in {0, 1} the two are the same rule - identical at z_i = 1, and at z_i = 0 the charge is
already forced to zero - so the designs allowed are exactly those the decisions describe, and the 0.0025 e^2 bound is
unchanged (selftest 9 checks both on the same instance). The difference is in the solver's relaxations, where z_i may
be fractional: without the indicator a hundredth of a sphere buys a hundredth of a charge at a ten-thousandth of the
budget, so the relaxed bound collapses. On the development grid at n = 4 the plain form relaxes to -19.32 kcal/mol
and HiGHS finds no design at all within 120 s; with the indicator the bound is far tighter and the same solve
finishes (selftest 9 and the timing table). A tangent that cannot bind for a
site's admitted charges is not added, which does not change the feasible set: for c > 0 it can only bind where
hi_i > c/2, for c < 0 only where lo_i < c/2, and c = 0 never binds because s_i >= 0 already.

THE SPACING AS CLIQUES
One rule per close pair (z_i + z_j <= 1) is valid but weak. Sites all within 3.25 A of one another form a clique,
of which at most one may hold a sphere, and one row then carries what many pairwise rows carried. The cliques
forbid exactly the designs the pairwise rules forbid (selftest 6 checks this on a small instance); the solver's
relaxations are tighter, so it finishes sooner.

WHAT THIS PRODUCES IS NOT A RESULT
On the development set this is a code test: the sites, L and A are the development ones
(cs_stage5/exploratory/dev_sites_a2.tsv, dev_L.tsv, dev_A/A_dev.tsv), not Stage 5 production. No setting is chosen
by the barriers its designs reach: B comes from the electronic-screen calibration, lambda from the confidence
interval, the rest are fixed in the pre-registration.

USAGE
    python3 cs_opt_milp.py selftest
    python3 cs_opt_milp.py solve <sites.tsv> <L.tsv> <A.tsv> [key=value ...]
    python3 cs_opt_milp.py time  <sites.tsv> <L.tsv> <A.tsv> [key=value ...]
Python 3.8+, numpy, scipy (milp/HiGHS). Checked on python 3.12.3 + scipy 1.11.4 and python 3.13 + scipy 1.18.1.
"""
import hashlib
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import csc_matrix, coo_matrix
from scipy.spatial import cKDTree

CONFIG = {
    "n": 4,                 # spheres in the design (O6)
    "lam": 0.94,            # semi-deviation weight; O5 ties it to F (0.94 at F = 30)
    "budget": 2.0,          # B in e^2, sum_i q_i^2 <= B (O3). Production B comes from the screens.
    "r_min": 3.25,          # A, minimum sphere-sphere separation (Stage 5)
    "tangents": 21,         # tangent lines standing in for q^2 (O3); 21 = step 0.1 on [-1, +1]
    "pool": 10,             # distinct designs per n (O7)
    "distinct_r": 1.5,      # A, the distinct-design radius (O7)
    "mip_gap": 1e-6,
    "time_limit": 600.0,
    "seed": 20261020,
    "strategy": "auto",     # full | lp | auto (try the full search, fall back to the LP-guided one)
    "probe": 120.0,         # s, how long 'auto' gives the full search before falling back
    "grow": 1.0,            # A, how far around the relaxation's own sites the candidate set reaches
}
TANGENT_SLACK = 0.0025      # e^2 per sphere, the gap the stand-in may leave (O3)


# ---------------------------------------------------------------- inputs

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_sites(path):
    """site table -> xyz, layer, and the admitted charge ranges per site.

    A site admits q in [q_min_neg, 0] when q_min_neg is -1.00, and q in [q_min_pos, +1] when q_min_pos is a number
    ("none" = no positive charge). Where both hold and q_min_pos is 0.00 the two meet, so they are one interval
    [-1, +1]; otherwise they are two ranges and at most one may be chosen.
    """
    rows = [l.rstrip("\n").split("\t") for l in Path(path).read_text().splitlines() if not l.startswith("#")]
    head, rows = rows[0], rows[1:]
    c = {k: head.index(k) for k in ("site", "layer", "x", "y", "z", "q_min_neg", "q_min_pos")}
    xyz = np.array([[float(r[c["x"]]), float(r[c["y"]]), float(r[c["z"]])] for r in rows])
    layer = [r[c["layer"]] for r in rows]
    ranges = []
    for r in rows:
        neg = r[c["q_min_neg"]] == "-1.00"
        pos = None if r[c["q_min_pos"]] == "none" else float(r[c["q_min_pos"]])
        if neg and pos is not None and pos == 0.0:
            rr = [(-1.0, 1.0)]                       # one interval
        else:
            rr = ([(-1.0, 0.0)] if neg else []) + ([(pos, 1.0)] if pos is not None else [])
        if not rr:
            raise SystemExit("site %s admits no charge at all" % r[c["site"]])
        ranges.append(rr)
    return xyz, layer, ranges


def read_matrix(path, nsite, name):
    M = np.array([[float(v) for v in l.split()] for l in Path(path).read_text().splitlines()
                  if not l.startswith("#")])
    if M.shape[0] != nsite:
        raise SystemExit("%s has %d rows, the site table has %d" % (name, M.shape[0], nsite))
    return M


# ---------------------------------------------------------------- spacing cliques

def close_pairs(xyz, r_min):
    """every pair of sites closer than r_min, as a sorted list."""
    return sorted(cKDTree(xyz).query_pairs(r_min))


def cliques(xyz, r_min, pairs=None):
    """Greedy cliques covering every close pair, as the decisions specify:
    take the close pairs in site order; for each pair not yet covered start a clique from it; repeatedly add the
    common neighbour with the most neighbours among the remaining candidates, ties by site number. Any pair left
    uncovered becomes its own two-site clique, so coverage is exact by construction.
    """
    pairs = close_pairs(xyz, r_min) if pairs is None else pairs
    nb = {}
    for i, j in pairs:
        nb.setdefault(i, set()).add(j)
        nb.setdefault(j, set()).add(i)
    covered, out = set(), []
    for i, j in pairs:
        if (i, j) in covered:
            continue
        K = [i, j]
        cand = nb[i] & nb[j]
        while cand:
            # the candidate with the most neighbours among the remaining candidates; ties by site number
            k = min(cand, key=lambda t: (-len(nb[t] & cand), t))
            K.append(k)
            cand = cand & nb[k]
        K.sort()
        for a in range(len(K)):
            for b in range(a + 1, len(K)):
                covered.add((K[a], K[b]))
        out.append(K)
    missing = [p for p in pairs if p not in covered]
    out.extend([list(p) for p in missing])
    return out, len(missing)


# ---------------------------------------------------------------- the MILP

class Problem:
    """The variable layout and the rows that do not depend on n, B or the cuts.

    Variables, in this order:
        y_k   k = 0..R-1   binary, one per admitted range (R = sum_i len(ranges_i))
        q_k   k = 0..R-1   the charge in that range
        s_i   i = 0..S-1   the stand-in for q_i^2
        d_f   f = 0..F-1   the frame's frozen score
        dbar               their mean
        u_f   f = 0..F-1   how far frame f comes out above the mean
    """

    def __init__(self, xyz, ranges, L, A, cfg, cliq=None, tie=True):
        # tie=False writes the tangents as s_i >= 2c q_i - c^2, without the indicator. Same rule at z in {0,1};
        # kept only so selftest 9 can show the two agree and how differently they relax.
        self.xyz, self.ranges, self.L, self.A, self.cfg, self.tie = xyz, ranges, L, A, cfg, tie
        self.S, self.F = A.shape
        self.site_of = np.array([i for i, rr in enumerate(ranges) for _ in rr])
        self.lo = np.array([lo for rr in ranges for lo, hi in rr])
        self.hi = np.array([hi for rr in ranges for lo, hi in rr])
        self.R = len(self.lo)
        self.ks_of = [[] for _ in ranges]          # site -> its range-variable indices
        for k, i in enumerate(self.site_of):
            self.ks_of[int(i)].append(k)
        self.i_y, self.i_q = 0, self.R
        self.i_s = 2 * self.R
        self.i_d = self.i_s + self.S
        self.i_dbar = self.i_d + self.F
        self.i_u = self.i_dbar + 1
        self.nvar = self.i_u + self.F
        self.cliques = cliq if cliq is not None else cliques(xyz, cfg["r_min"])[0]
        self._static()

    def _static(self):
        """rows that never change: range coupling, one range per site, the tangents, d_f, dbar, u_f."""
        r, c, v, lo, hi = [], [], [], [], []
        row = 0

        def put(rr, cc, vv):
            r.extend([rr] * len(cc)); c.extend(cc); v.extend(vv)

        # q_k - hi_k y_k <= 0  and  lo_k y_k - q_k <= 0
        for k in range(self.R):
            put(row, [self.i_q + k, self.i_y + k], [1.0, -self.hi[k]]); lo.append(-np.inf); hi.append(0.0); row += 1
            put(row, [self.i_q + k, self.i_y + k], [-1.0, self.lo[k]]); lo.append(-np.inf); hi.append(0.0); row += 1
        # at most one range per site, only where a site has more than one
        for i, rr in enumerate(self.ranges):
            if len(rr) > 1:
                ks = self.ks_of[i]
                put(row, [self.i_y + k for k in ks], [1.0] * len(ks)); lo.append(-np.inf); hi.append(1.0); row += 1
        # tangents: s_i - sum_k(in i) 2c q_k >= -c^2, only where the tangent can bind
        T = np.linspace(-1.0, 1.0, self.cfg["tangents"])
        self.n_tangent = 0
        for i, rr in enumerate(self.ranges):
            ks = self.ks_of[i]
            s_lo, s_hi = min(self.lo[k] for k in ks), max(self.hi[k] for k in ks)
            for t in T:
                if t == 0.0 or (t > 0 and not s_hi > t / 2) or (t < 0 and not s_lo < t / 2):
                    continue
                # s_i - sum_k 2c q_k + c^2 sum_k y_k >= 0, the tangent with its constant tied to the indicator
                if self.tie:
                    put(row, [self.i_s + i] + [self.i_q + k for k in ks] + [self.i_y + k for k in ks],
                        [1.0] + [-2.0 * t] * len(ks) + [t * t] * len(ks))
                    lo.append(0.0)
                else:
                    put(row, [self.i_s + i] + [self.i_q + k for k in ks], [1.0] + [-2.0 * t] * len(ks))
                    lo.append(-t * t)
                hi.append(np.inf); row += 1
                self.n_tangent += 1
        # d_f - sum_i a_if q_i - sum_i L_if z_i = 0, written over ranges (q_k and y_k)
        for f in range(self.F):
            cc = [self.i_d + f]; vv = [1.0]
            for i in range(self.S):
                a, Lf = float(self.A[i, f]), float(self.L[i, f])
                for k in self.ks_of[i]:
                    if a != 0.0:
                        cc.append(self.i_q + k); vv.append(-a)
                    if Lf != 0.0:
                        cc.append(self.i_y + k); vv.append(-Lf)
            put(row, cc, vv); lo.append(0.0); hi.append(0.0); row += 1
        # dbar - mean_f d_f = 0
        put(row, [self.i_dbar] + [self.i_d + f for f in range(self.F)],
            [1.0] + [-1.0 / self.F] * self.F); lo.append(0.0); hi.append(0.0); row += 1
        # u_f - d_f + dbar >= 0
        for f in range(self.F):
            put(row, [self.i_u + f, self.i_d + f, self.i_dbar], [1.0, -1.0, 1.0])
            lo.append(0.0); hi.append(np.inf); row += 1
        self._sr, self._sc, self._sv = r, c, v
        self._slo, self._shi, self._srow = lo, hi, row

    def constraints(self, n, budget, cuts=()):
        """the static rows plus the ones that depend on n, B and the cuts."""
        r, c, v = list(self._sr), list(self._sc), list(self._sv)
        lo, hi = list(self._slo), list(self._shi)
        row = self._srow

        def put(rr, cc, vv):
            r.extend([rr] * len(cc)); c.extend(cc); v.extend(vv)

        # sum_i z_i = n, with z_i = sum_r y_ir
        put(row, [self.i_y + k for k in range(self.R)], [1.0] * self.R)
        lo.append(float(n)); hi.append(float(n)); row += 1
        # sum_i s_i <= B
        put(row, [self.i_s + i for i in range(self.S)], [1.0] * self.S)
        lo.append(-np.inf); hi.append(float(budget)); row += 1
        # spacing cliques
        for K in self.cliques:
            ks = [k for i in K for k in self.ks_of[i]]
            put(row, [self.i_y + k for k in ks], [1.0] * len(ks))
            lo.append(-np.inf); hi.append(1.0); row += 1
        # distinct-design cuts: at most n-1 spheres within distinct_r of an earlier design's spheres
        for near in cuts:
            ks = [k for i in near for k in self.ks_of[i]]
            put(row, [self.i_y + k for k in ks], [1.0] * len(ks))
            lo.append(-np.inf); hi.append(float(n - 1)); row += 1
        Am = coo_matrix((v, (r, c)), shape=(row, self.nvar)).tocsc()
        return LinearConstraint(Am, np.array(lo), np.array(hi)), row

    def bounds_objective(self):
        lb = np.full(self.nvar, -np.inf); ub = np.full(self.nvar, np.inf)
        lb[self.i_y:self.i_y + self.R], ub[self.i_y:self.i_y + self.R] = 0.0, 1.0
        lb[self.i_q:self.i_q + self.R] = np.minimum(self.lo, 0.0)
        ub[self.i_q:self.i_q + self.R] = np.maximum(self.hi, 0.0)
        lb[self.i_s:self.i_s + self.S], ub[self.i_s:self.i_s + self.S] = 0.0, 1.0
        lb[self.i_u:self.i_u + self.F] = 0.0
        integrality = np.zeros(self.nvar); integrality[self.i_y:self.i_y + self.R] = 1
        c = np.zeros(self.nvar)
        c[self.i_dbar] = 1.0
        c[self.i_u:self.i_u + self.F] = self.cfg["lam"] / self.F
        return Bounds(lb, ub), integrality, c

    def solve(self, n, budget, cuts=(), time_limit=None):
        cons, nrow = self.constraints(n, budget, cuts)
        bounds, integrality, c = self.bounds_objective()
        t0 = time.time()
        res = milp(c=c, constraints=cons, integrality=integrality, bounds=bounds,
                   options={"mip_rel_gap": self.cfg["mip_gap"],
                            "time_limit": float(self.cfg["time_limit"] if time_limit is None else time_limit)})
        wall = time.time() - t0
        if not res.success or res.x is None:
            return None, dict(message=res.message, wall=wall, nrow=nrow, nvar=self.nvar)
        return self.read(res, n, budget, wall, nrow), dict(message=res.message, wall=wall,
                                                           nrow=nrow, nvar=self.nvar)

    def read(self, res, n, budget, wall, nrow):
        x = res.x
        y = x[self.i_y:self.i_y + self.R]
        on = [k for k in range(self.R) if y[k] > 0.5]
        q = {}
        for k in on:
            q[int(self.site_of[k])] = float(x[self.i_q + k])
        d = x[self.i_d:self.i_d + self.F]
        sites = sorted(q)
        sol = dict(
            n=n, budget=float(budget), objective=float(res.fun), sites=sites,
            q=[q[i] for i in sites], mean=float(x[self.i_dbar]),
            semidev=float(np.sum(x[self.i_u:self.i_u + self.F]) / self.F),
            d=[float(t) for t in d],
            sumq2=float(sum(v * v for v in q.values())),              # the TRUE budget use
            s_sum=float(np.sum(x[self.i_s:self.i_s + self.S])),       # what the stand-in charged for it
            net=float(sum(q.values())), n_neutral=int(sum(1 for v in q.values() if abs(v) < 1e-9)),
            gap=float(getattr(res, "mip_gap", float("nan"))),
            dual_bound=float(getattr(res, "mip_dual_bound", float("nan"))),
            status=int(res.status), wall=wall, nrow=nrow, nvar=self.nvar,
        )
        return sol

    def near(self, sites):
        """sites within distinct_r of any sphere of a design (O7's cut)."""
        tree = cKDTree(self.xyz)
        out = set()
        for i in sites:
            out.update(tree.query_ball_point(self.xyz[i], self.cfg["distinct_r"]))
        return sorted(out)

    def relax(self, n, budget, cuts=(), time_limit=None):
        """the relaxation of the FULL problem: a valid lower bound on every design it admits, and the sites it uses."""
        cons, _ = self.constraints(n, budget, cuts)
        bounds, _, c = self.bounds_objective()
        t0 = time.time()
        res = milp(c=c, constraints=cons, integrality=np.zeros(self.nvar), bounds=bounds,
                   options={"time_limit": float(self.cfg["time_limit"] if time_limit is None else time_limit)})
        if res.x is None:
            return None, None, time.time() - t0
        y = res.x[self.i_y:self.i_y + self.R]
        used = sorted({int(self.site_of[k]) for k in range(self.R) if y[k] > 1e-6})
        return float(res.fun), used, time.time() - t0

    def restrict(self, sites):
        """the same problem on a subset of the sites. Its designs are designs of the full problem: the rules are the
        same and the induced cliques still cover every close pair inside the subset, since any pair was covered by a
        clique holding both."""
        sub = sorted(set(int(i) for i in sites))
        idx = {g: t for t, g in enumerate(sub)}
        cl = [[idx[i] for i in K if i in idx] for K in self.cliques]
        cl = [K for K in cl if len(K) > 1]
        q = Problem(self.xyz[sub], [self.ranges[i] for i in sub], self.L[sub], self.A[sub], self.cfg,
                    cliq=cl, tie=self.tie)
        q.parent_sites = sub
        return q

    def solve_lp_guided(self, n, budget, cuts=(), time_limit=None, grow=None):
        """Relax the full problem for a bound, then search only the sites that relaxation uses.

        The design returned is a design of the FULL problem - the restricted problem's rules are the full problem's,
        on fewer sites - and the bound is the full problem's, so the difference between them is a true optimality
        gap: no design anywhere on the grid can beat this one by more than that. This is how n beyond about 9 is
        reached at all: the relaxation is nearly exact here (at n = 6 it equals the optimum), but the solver cannot
        find the matching design among thousands of near-equivalent sites inside any usable time.
        """
        grow = self.cfg["grow"] if grow is None else grow
        LB, used, t_lp = self.relax(n, budget, cuts, time_limit)
        if LB is None:
            return None, dict(message="the relaxation itself did not solve", wall=t_lp)
        cand = set(used)
        if grow > 0:
            tree = cKDTree(self.xyz)
            for i in used:
                cand.update(tree.query_ball_point(self.xyz[i], grow))
        q = self.restrict(cand)
        back = {t: g for t, g in enumerate(q.parent_sites)}
        sub_cuts = [[t for t, g in back.items() if g in set(cc)] for cc in cuts]
        sol, info = q.solve(n, budget, sub_cuts, time_limit)
        info["wall"] = info.get("wall", 0.0) + t_lp
        if sol is None:
            return None, info
        sol["sites"] = [back[t] for t in sol["sites"]]
        sol["bound"] = LB
        sol["certified_gap"] = sol["objective"] - LB
        sol["candidates"] = len(q.parent_sites)
        sol["relax_sites"] = len(used)
        sol["wall"] = info["wall"]
        sol["strategy"] = "lp"
        return sol, info

    def best(self, n, budget, cuts=(), time_limit=None):
        """The strategy the settings ask for.

        'auto' gives the full search `probe` seconds and keeps its answer only if it closed the gap itself;
        otherwise it falls back to the LP-guided search, which always returns a design and a bound. The fallback is
        decided by what the solver does on THIS machine rather than by a chosen n, because how far the full search
        reaches depends on the HiGHS build: the same n = 7, 8 and 9 solved in 39, 102 and 194 s under scipy 1.18.1
        and returned nothing at all in 600 s under scipy 1.11.4.
        """
        st = self.cfg["strategy"]
        if st == "lp":
            return self.solve_lp_guided(n, budget, cuts, time_limit)
        if st == "auto":
            probe = min(float(self.cfg["probe"]), float(self.cfg["time_limit"] if time_limit is None else time_limit))
            sol, info = self.solve(n, budget, cuts, probe)
            if sol is not None:
                db = sol.get("dual_bound", float("nan"))
                sol["bound"] = db if np.isfinite(db) else float("-inf")
                sol["certified_gap"] = max(0.0, sol["objective"] - db) if np.isfinite(db) else float("inf")
                sol["strategy"] = "full"
            if sol is not None and sol["status"] == 0:
                return sol, info          # the full search finished: optimal to its gap tolerance (1e-6)
            out, info2 = self.solve_lp_guided(n, budget, cuts, time_limit)
            probe_wall = info.get("wall", 0.0)
            if out is None:
                # the LP-guided search found nothing; an unproven design from the probe is still a design, with the
                # solver's own bound and gap reported as they are
                if sol is not None:
                    sol["wall"] = sol["wall"] + info2.get("wall", 0.0)
                return sol, (info if sol is not None else info2)
            out["wall"] = out.get("wall", 0.0) + probe_wall
            if sol is not None:
                # both bounds are valid for the whole grid: the relaxation's, and whatever the probe's search had
                # proved when it stopped. The larger is the tighter. Keep whichever design scores better.
                bound = max(out["bound"], sol["bound"])
                best_ = sol if sol["objective"] < out["objective"] - 1e-12 else out
                best_["bound"] = bound
                best_["certified_gap"] = best_["objective"] - bound
                best_["wall"] = out["wall"]
                return best_, info2
            return out, info2
        sol, info = self.solve(n, budget, cuts, time_limit)
        if sol is not None:
            sol.setdefault("strategy", "full")
            # the solver's own dual bound is a valid bound on this problem, so objective - bound is the gap in
            # kcal/mol. It is 0 when the search closed, and the real distance when it stopped at the time limit.
            db = sol.get("dual_bound", float("nan"))
            sol.setdefault("bound", db)
            sol.setdefault("certified_gap", max(0.0, sol["objective"] - db) if np.isfinite(db) else float("nan"))
        return sol, info

    def pool(self, n, budget, k, time_limit=None):
        """the k best distinct designs at this n (O7)."""
        found, cuts = [], []
        for _ in range(k):
            sol, info = self.best(n, budget, cuts, time_limit)
            if sol is None:
                break
            found.append(sol)
            cuts.append(self.near(sol["sites"]))
        return found


# ---------------------------------------------------------------- checks on a solution

def violations(p, sol):
    """every rule, checked on the returned numbers rather than trusted."""
    out = {}
    out["cardinality"] = abs(len(sol["sites"]) - sol["n"])
    bad = 0.0
    for i, qv in zip(sol["sites"], sol["q"]):
        ok = any(lo - 1e-9 <= qv <= hi + 1e-9 for lo, hi in p.ranges[i])
        bad = max(bad, 0.0 if ok else 1.0)
    out["charge_range"] = bad
    d = 0.0
    for a in range(len(sol["sites"])):
        for b in range(a + 1, len(sol["sites"])):
            r = float(np.linalg.norm(p.xyz[sol["sites"][a]] - p.xyz[sol["sites"][b]]))
            d = max(d, max(0.0, p.cfg["r_min"] - r))
    out["spacing"] = d
    out["budget_true"] = max(0.0, sol["sumq2"] - sol["budget"])
    out["budget_slack_per_sphere"] = out["budget_true"] / max(1, sol["n"])
    q = np.zeros(p.S)
    for i, qv in zip(sol["sites"], sol["q"]):
        q[i] = qv
    z = np.zeros(p.S); z[sol["sites"]] = 1.0
    d_chk = p.A.T @ q + p.L.T @ z
    out["d_f"] = float(np.max(np.abs(d_chk - np.array(sol["d"]))))
    mean = float(np.mean(d_chk))
    out["mean"] = abs(mean - sol["mean"])
    sd = float(np.mean(np.maximum(0.0, d_chk - mean)))
    out["semidev"] = abs(sd - sol["semidev"])
    out["objective"] = abs(mean + p.cfg["lam"] * sd - sol["objective"])
    return out


# ---------------------------------------------------------------- self-tests

def _mk(xyz, ranges, L, A, **kw):
    return Problem(xyz, ranges, L, A, dict(CONFIG, **kw))


def selftest():
    """Checks on instances whose answers are known independently of the solver."""
    print("SELF-TEST  (cs_opt_milp.py)\n")
    ok = True
    rng = np.random.default_rng(CONFIG["seed"])

    def say(good, text):
        nonlocal ok
        ok &= bool(good)
        print("   [%s] %s" % ("PASS" if good else "FAIL", text))

    print("1. one sphere, closed form")
    # With n = 1, lambda = 0 and a budget that does not bind, the score is min over sites and admitted q of
    # mean_f (a_if q + L_if). Each site's best q is at a bound of its range, so the answer is computable directly.
    S, F = 25, 4
    A = rng.normal(0, 1, (S, F)); L = rng.normal(0, 0.01, (S, F))
    xyz = rng.normal(0, 40, (S, 3))                      # far apart: no spacing rule binds
    ranges = [[(-1.0, 1.0)] if i % 3 == 0 else ([(-1.0, 0.0)] if i % 3 == 1 else [(0.25, 1.0)]) for i in range(S)]
    p = _mk(xyz, ranges, L, A, n=1, lam=0.0, budget=1.0)
    sol, info = p.solve(1, 1.0)
    am, Lm = A.mean(1), L.mean(1)
    best, best_i, best_q = np.inf, None, None
    for i in range(S):
        for lo, hi in ranges[i]:
            for qv in (lo, hi):
                val = am[i] * qv + Lm[i]
                if val < best:
                    best, best_i, best_q = val, i, qv
    say(sol is not None and abs(sol["objective"] - best) < 1e-6 and sol["sites"] == [best_i],
        "MILP %+.6f at site %s q %+.2f; closed form %+.6f at site %d q %+.2f"
        % (sol["objective"], sol["sites"], sol["q"][0], best, best_i, best_q))

    print("\n2. every rule honoured, on a crowded instance with both layers")
    S, F = 120, 6
    A = rng.normal(0, 2, (S, F)); L = rng.normal(0, 0.02, (S, F))
    xyz = rng.normal(0, 5, (S, 3))                        # crowded: spacing binds
    ranges = [[(-1.0, 1.0)] if i % 4 == 0 else ([(-1.0, 0.0), (0.5, 1.0)] if i % 4 == 1 else
              ([(-1.0, 0.0)] if i % 4 == 2 else [(0.25, 1.0)])) for i in range(S)]
    p = _mk(xyz, ranges, L, A, n=4, lam=0.94, budget=2.0)
    sol, info = p.solve(4, 2.0)
    v = violations(p, sol)
    say(v["cardinality"] == 0, "exactly n spheres")
    say(v["charge_range"] == 0, "every charge inside a range its site admits")
    say(v["spacing"] == 0.0, "no two spheres closer than %.2f A" % p.cfg["r_min"])
    say(v["d_f"] < 1e-7, "d_f = sum_i a_if q_i + sum_i L_if z_i (worst %.1e)" % v["d_f"])
    say(v["mean"] < 1e-7 and v["semidev"] < 1e-7, "dbar is the mean and u is the semi-deviation")
    say(v["objective"] < 1e-6, "objective = mean + lambda x semi-deviation")

    print("\n3. the q^2 stand-in stays inside its bound")
    # The stand-in may let the true sum of squares exceed B, by at most 0.0025 e^2 per sphere.
    worst = 0.0
    for nn, bb in ((2, 0.5), (4, 1.0), (6, 2.0)):
        s2, _ = p.solve(nn, bb)
        if s2 is None:
            continue
        vv = violations(p, s2)
        worst = max(worst, vv["budget_slack_per_sphere"])
        print("      n=%d B=%.2f: true sum q^2 = %.4f, stand-in charged %.4f, over by %.4f e^2 (%.5f per sphere)"
              % (nn, bb, s2["sumq2"], s2["s_sum"], vv["budget_true"], vv["budget_slack_per_sphere"]))
    say(worst <= TANGENT_SLACK + 1e-9, "worst overshoot %.5f e^2 per sphere, bound %.4f" % (worst, TANGENT_SLACK))

    print("\n4. with lambda = 0 and no bound active, the charges are proportional to -a")
    # On a FIXED set of sites with q free in [-1, +1], minimising mean_f sum_i a_if q_i under sum q^2 <= B has
    # q = -sqrt(B) * abar / ||abar||, and the score -sqrt(B)*||abar|| + mean L. Checked against the MILP with the
    # sites forced by making every other site unavailable.
    k = 3
    S4, F4 = k, 5
    A4 = rng.normal(0, 1, (S4, F4)); L4 = np.zeros((S4, F4))
    xyz4 = rng.normal(0, 40, (S4, 3))
    ranges4 = [[(-1.0, 1.0)] for _ in range(S4)]
    B4 = 0.3                                              # small enough that no |q| reaches 1
    p4 = _mk(xyz4, ranges4, L4, A4, n=k, lam=0.0, budget=B4, tangents=201)
    s4, _ = p4.solve(k, B4)
    abar = A4.mean(1)
    qstar = -np.sqrt(B4) * abar / np.linalg.norm(abar)
    got = np.array(s4["q"])
    say(np.max(np.abs(got - qstar)) < 2e-2,
        "charges %s against the closed form %s (worst %.4f; 201 tangents)"
        % (np.round(got, 3).tolist(), np.round(qstar, 3).tolist(), float(np.max(np.abs(got - qstar)))))
    say(abs(s4["objective"] - (-np.sqrt(B4) * np.linalg.norm(abar))) < 2e-2,
        "score %+.4f against the closed form %+.4f" % (s4["objective"], -np.sqrt(B4) * np.linalg.norm(abar)))

    print("\n5. with a bound active, the optimality conditions hold directly")
    # Same, but B large enough that some |q| would exceed 1. Then at the optimum each site either sits at a bound
    # of its range, or has q_i = -a_i / (2 mu) for one common mu > 0 (the budget's multiplier). Checked by
    # recovering mu from an interior site and verifying every other site against it.
    B5 = 2.6
    p5 = _mk(xyz4, ranges4, L4, A4, n=k, lam=0.0, budget=B5, tangents=201)
    s5, _ = p5.solve(k, B5)
    qq, aa = np.array(s5["q"]), abar
    interior = [i for i in range(k) if abs(abs(qq[i]) - 1.0) > 1e-3]
    if interior:
        mu = -aa[interior[0]] / (2.0 * qq[interior[0]])
        resid = max(abs(qq[i] + aa[i] / (2.0 * mu)) for i in interior)
        at_bound_ok = all((qq[i] > 0) == (aa[i] < 0) for i in range(k) if i not in interior)
        say(mu > 0 and resid < 2e-2 and at_bound_ok,
            "mu = %.4f > 0; interior sites satisfy q = -a/(2mu) (worst %.4f); bounded sites take the sign that lowers the score"
            % (mu, resid))
    else:
        say(all((qq[i] > 0) == (aa[i] < 0) for i in range(k)),
            "every site at a bound, each with the sign that lowers the score")
    say(s5["sumq2"] <= B5 + k * TANGENT_SLACK + 1e-9,
        "the budget still holds: sum q^2 = %.4f against B = %.2f" % (s5["sumq2"], B5))

    print("\n6. the cliques forbid exactly what the pairwise rules forbid")
    S6 = 40
    xyz6 = rng.normal(0, 3.2, (S6, 3))
    A6 = rng.normal(0, 1, (S6, 6)); L6 = np.zeros((S6, 6))
    ranges6 = [[(-1.0, 1.0)] for _ in range(S6)]
    pairs = close_pairs(xyz6, 3.25)
    cl, missed = cliques(xyz6, 3.25, pairs)
    cover = set()
    for K in cl:
        for a in range(len(K)):
            for b in range(a + 1, len(K)):
                cover.add((K[a], K[b]))
    say(missed == 0 and all(pp in cover for pp in pairs),
        "%d close pairs all covered by %d cliques (%d pairs needed their own)" % (len(pairs), len(cl), missed))
    inside = all(float(np.linalg.norm(xyz6[K[a]] - xyz6[K[b]])) < 3.25
                 for K in cl for a in range(len(K)) for b in range(a + 1, len(K)))
    say(inside, "every pair inside a clique really is closer than 3.25 A, so no valid design is excluded")
    p6c = _mk(xyz6, ranges6, L6, A6, n=3, lam=0.0, budget=3.0)
    p6p = _mk(xyz6, ranges6, L6, A6, n=3, lam=0.0, budget=3.0, )
    p6p.cliques = [list(pp) for pp in pairs]              # the pairwise form
    a_c, _ = p6c.solve(3, 3.0)
    a_p, _ = p6p.solve(3, 3.0)
    say(abs(a_c["objective"] - a_p["objective"]) < 1e-6,
        "same optimum from cliques (%+.6f, %d rows) and pairs (%+.6f, %d rows)"
        % (a_c["objective"], a_c["nrow"], a_p["objective"], a_p["nrow"]))

    print("\n7. distinct-design cuts return designs that are not shifted copies")
    pool = p.pool(3, 1.5, 4)
    say(len(pool) >= 2, "%d designs found" % len(pool))
    sup = [tuple(s["sites"]) for s in pool]
    say(len(set(sup)) == len(sup), "all supports distinct")
    tree = cKDTree(p.xyz)
    far_enough = True
    for a in range(len(pool)):
        for b in range(a):
            near_b = set(pool[b]["sites"])
            for i in pool[b]["sites"]:
                near_b.update(tree.query_ball_point(p.xyz[i], p.cfg["distinct_r"]))
            if not any(i not in near_b for i in pool[a]["sites"]):
                far_enough = False
    say(far_enough, "each design places at least one sphere more than %.1f A from every sphere of each earlier design"
        % p.cfg["distinct_r"])
    say(all(pool[i]["objective"] <= pool[i + 1]["objective"] + 1e-9 for i in range(len(pool) - 1)),
        "the pool is in order, best first: " + "  ".join("%+.4f" % s["objective"] for s in pool))

    print("\n8a. with no electrostatic term, L alone chooses the sites")
    # A = 0 exactly, so the score is mean_f sum_i L_if z_i and the design must take the n most negative mean-L
    # sites. An independent check that L enters the score at all.
    S8, F8 = 30, 4
    L8 = rng.normal(0, 0.05, (S8, F8))
    xyz8 = rng.normal(0, 40, (S8, 3))
    ranges8 = [[(-1.0, 1.0)] for _ in range(S8)]
    p8 = _mk(xyz8, ranges8, L8, np.zeros((S8, F8)), n=3, lam=0.0, budget=1.0)
    s8, _ = p8.solve(3, 1.0)
    want = sorted(np.argsort(L8.mean(1))[:3].tolist())
    say(s8["sites"] == want, "sites %s, the three most negative mean L (expected %s)" % (s8["sites"], want))

    print("\n8b. B = 0 leaves the charges at the stand-in's own bound, not at zero")
    # The tangent envelope is zero for |q| <= 0.05 (half the 0.1 tangent step), so B = 0 does NOT force q = 0: it
    # allows |q| <= 0.05 per sphere, which is the stand-in's slack expressed as a charge. The closed form per site
    # is then mean L - 0.05 |mean a|, and that is what the design must minimise.
    A8 = rng.normal(0, 1, (S8, F8))
    p8b = _mk(xyz8, ranges8, L8, A8, n=3, lam=0.0, budget=0.0)
    s8b, _ = p8b.solve(3, 0.0)
    qcap = 0.5 * (2.0 / (CONFIG["tangents"] - 1))          # half the tangent step
    per = L8.mean(1) - qcap * np.abs(A8.mean(1))
    want_b = sorted(np.argsort(per)[:3].tolist())
    say(max(abs(q) for q in s8b["q"]) <= qcap + 1e-6,
        "every |q| <= %.3f e, the envelope's bound at B = 0 (largest %.4f)" % (qcap, max(abs(q) for q in s8b["q"])))
    say(s8b["sites"] == want_b and abs(s8b["objective"] - per[want_b].sum()) < 1e-6,
        "sites %s and score %+.6f against the closed form %s, %+.6f"
        % (s8b["sites"], s8b["objective"], want_b, per[want_b].sum()))
    say(s8b["sumq2"] <= 3 * TANGENT_SLACK + 1e-9,
        "the true sum q^2 is %.5f, within the 0.0025 e^2 per sphere the stand-in may leave" % s8b["sumq2"])

    print("\n9. the two tangent forms allow the same designs, and the tied one relaxes tighter")
    S9, F9 = 80, 5
    A9 = rng.normal(0, 1.5, (S9, F9)); L9 = rng.normal(0, 0.02, (S9, F9))
    xyz9 = rng.normal(0, 6, (S9, 3))
    r9 = [[(-1.0, 1.0)] if i % 3 == 0 else ([(-1.0, 0.0)] if i % 3 == 1 else [(0.25, 1.0)]) for i in range(S9)]
    cfg9 = dict(CONFIG, n=4, lam=0.94, budget=1.5)
    p_tie = Problem(xyz9, r9, L9, A9, cfg9, tie=True)
    p_plain = Problem(xyz9, r9, L9, A9, cfg9, tie=False)
    s_tie, _ = p_tie.solve(4, 1.5, time_limit=300.0)
    s_plain, _ = p_plain.solve(4, 1.5, time_limit=300.0)
    say(s_tie is not None and s_plain is not None
        and abs(s_tie["objective"] - s_plain["objective"]) < 1e-6,
        "same optimum: tied %+.6f, plain %+.6f" % (s_tie["objective"], s_plain["objective"]))
    rel = {}
    for name, pp in (("tied", p_tie), ("plain", p_plain)):
        cons, _ = pp.constraints(4, 1.5)
        bounds, _, cvec = pp.bounds_objective()
        res = milp(c=cvec, constraints=cons, integrality=np.zeros(pp.nvar), bounds=bounds,
                   options={"time_limit": 300.0})
        rel[name] = float(res.fun)
    say(rel["tied"] > rel["plain"] + 1e-6,
        "relaxed bound: tied %+.4f against plain %+.4f, so the tied form leaves the solver far less to search"
        % (rel["tied"], rel["plain"]))
    v9 = violations(p_tie, s_tie)
    say(v9["budget_slack_per_sphere"] <= TANGENT_SLACK + 1e-9,
        "the tied form keeps the same budget bound (over by %.5f e^2 per sphere)" % v9["budget_slack_per_sphere"])

    print("\n10. the LP-guided search returns a design of the FULL problem, with a true bound")
    S10, F10 = 200, 6
    A10 = rng.normal(0, 1.5, (S10, F10)); L10 = rng.normal(0, 0.02, (S10, F10))
    xyz10 = rng.normal(0, 7, (S10, 3))
    r10 = [[(-1.0, 1.0)] if i % 3 == 0 else ([(-1.0, 0.0)] if i % 3 == 1 else [(0.25, 1.0)]) for i in range(S10)]
    cfg10 = dict(CONFIG, n=5, lam=0.94, budget=1.5, grow=1.0)
    p10 = _mk(xyz10, r10, L10, A10, **{k: v for k, v in cfg10.items() if k in CONFIG})
    s_full, _ = p10.solve(5, 1.5, time_limit=300.0)
    s_lp, _ = p10.solve_lp_guided(5, 1.5, time_limit=300.0)
    say(s_lp is not None and s_full is not None, "both searches returned a design")
    v10 = violations(p10, s_lp)
    say(max(v10["cardinality"], v10["charge_range"], v10["spacing"], v10["d_f"], v10["objective"]) < 1e-6,
        "the restricted design obeys every rule of the full problem and its arithmetic re-checks")
    say(s_lp["objective"] >= s_full["objective"] - 1e-6,
        "it is no better than the full optimum (%+.4f against %+.4f), as a design on fewer sites must be"
        % (s_lp["objective"], s_full["objective"]))
    say(s_lp["bound"] <= s_full["objective"] + 1e-6,
        "the bound %+.4f really does bound the full optimum %+.4f" % (s_lp["bound"], s_full["objective"]))
    say(abs(s_lp["certified_gap"] - (s_lp["objective"] - s_lp["bound"])) < 1e-9
        and s_lp["objective"] - s_full["objective"] <= s_lp["certified_gap"] + 1e-9,
        "the reported gap %.4f is objective minus bound, and covers the real shortfall %.4f"
        % (s_lp["certified_gap"], s_lp["objective"] - s_full["objective"]))

    print("\n11. 'auto' keeps the full search when it finishes, and falls back correctly when it does not")
    p11 = _mk(xyz10, r10, L10, A10, **{k: v for k, v in dict(cfg10, strategy="auto", probe=300.0).items()
                                       if k in CONFIG})
    a_ok, _ = p11.best(5, 1.5, time_limit=300.0)
    say(a_ok is not None and a_ok["strategy"] == "full" and abs(a_ok["objective"] - s_full["objective"]) < 1e-6
        and a_ok["certified_gap"] <= 1e-5 * max(1.0, abs(a_ok["objective"])),
        "with time to finish, auto returns the full search's design (%+.4f, gap %.1e)"
        % (a_ok["objective"], a_ok["certified_gap"]))
    p11.cfg["probe"] = 0.0                       # a zero probe cannot finish: the fallback is forced
    a_fb, _ = p11.best(5, 1.5, time_limit=300.0)
    say(a_fb is not None and a_fb["strategy"] == "lp"
        and a_fb["bound"] <= s_full["objective"] + 1e-6 <= a_fb["objective"] + 2e-6
        and a_fb["objective"] - s_full["objective"] <= a_fb["certified_gap"] + 1e-9,
        "with no time, auto falls back to the LP-guided search; its bound and gap still hold (%+.4f, gap %.4f)"
        % (a_fb["objective"], a_fb["certified_gap"]))
    # the choice between an unproven probe design and a worse fallback design is pure logic, so test it directly
    probe_design = dict(s_full, status=1, dual_bound=s_full["objective"] - 0.30)
    worse = dict(s_lp, objective=s_full["objective"] + 0.05, bound=s_full["objective"] - 0.20, strategy="lp")
    p11.cfg["probe"] = 5.0
    real_solve, real_lp = p11.solve, p11.solve_lp_guided
    p11.solve = lambda *a, **k: (dict(probe_design), dict(wall=1.0))
    p11.solve_lp_guided = lambda *a, **k: (dict(worse), dict(wall=2.0))
    pick, _ = p11.best(5, 1.5, time_limit=300.0)
    p11.solve, p11.solve_lp_guided = real_solve, real_lp
    say(pick["strategy"] == "full" and abs(pick["objective"] - s_full["objective"]) < 1e-12
        and abs(pick["bound"] - (s_full["objective"] - 0.20)) < 1e-12
        and abs(pick["certified_gap"] - 0.20) < 1e-12,
        "an unproven probe design that beats the fallback is kept, with the tighter of the two bounds")

    print("\n%s" % ("ALL SELF-TESTS PASSED" if ok else "SELF-TEST FAILURES ABOVE"))
    return 0 if ok else 1


# ---------------------------------------------------------------- command line

def load(sites_p, L_p, A_p):
    xyz, layer, ranges = read_sites(sites_p)
    L = read_matrix(L_p, len(xyz), "L")
    A = read_matrix(A_p, len(xyz), "A")
    if L.shape[1] != A.shape[1]:
        raise SystemExit("L has %d frames, A has %d" % (L.shape[1], A.shape[1]))
    print("inputs")
    for p_ in (sites_p, L_p, A_p):
        print("   %s  sha256 %s" % (p_, sha256(p_)))
    nr = sum(len(r) for r in ranges)
    print("   %d sites (anion/neutral %d, cation %d), %d frames, %d yes/no choices"
          % (len(xyz), sum(1 for t in layer if t == "anion"), sum(1 for t in layer if t == "cation"),
             A.shape[1], nr))
    return xyz, layer, ranges, L, A


def describe(p, pairs):
    sizes = [len(K) for K in p.cliques]
    print("   %d close pairs under %.2f A; %d cliques (mean size %.1f, %d entries); %d tangent rows"
          % (len(pairs), p.cfg["r_min"], len(p.cliques), float(np.mean(sizes)),
             sum(len(K) for K in p.cliques), p.n_tangent))


def show(sol, p):
    v = violations(p, sol)
    print("   n=%-2d objective %+9.4f  (mean %+9.4f, semi-deviation %7.4f)  solver gap %.1e  %.1fs  [%s]"
          % (sol["n"], sol["objective"], sol["mean"], sol["semidev"], sol["gap"], sol["wall"],
             sol.get("strategy", "full")))
    if sol.get("strategy") == "lp":
        print("      the relaxation of the whole grid bounds every design at %+.4f and used %d sites; the search ran"
              % (sol["bound"], sol["relax_sites"]))
        print("      over %d candidates, so no design anywhere on the grid beats this one by more than %.4f kcal/mol"
              % (sol["candidates"], sol["certified_gap"]))
    print("      sites %s" % " ".join(str(i) for i in sol["sites"]))
    print("      q     %s" % " ".join("%+.3f" % q for q in sol["q"]))
    print("      net charge %+.3f, %d neutral sphere(s); true sum q^2 %.4f of B %.2f (over by %.5f e^2 per sphere)"
          % (sol["net"], sol["n_neutral"], sol["sumq2"], sol["budget"], v["budget_slack_per_sphere"]))
    worst = max(v["cardinality"], v["charge_range"], v["spacing"], v["d_f"], v["mean"], v["semidev"], v["objective"])
    print("      rules and arithmetic re-checked on the returned numbers: worst violation %.1e" % worst)
    return worst


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    mode = sys.argv[1]
    if mode == "selftest":
        return selftest()
    if mode not in ("solve", "time") or len(sys.argv) < 5:
        sys.exit(__doc__)
    cfg = dict(CONFIG)
    for kv in sys.argv[5:]:
        k, val = kv.split("=", 1)
        if k not in cfg:
            sys.exit("unknown parameter %r. Known: %s" % (k, ", ".join(sorted(cfg))))
        cfg[k] = type(cfg[k])(val)
    xyz, layer, ranges, L, A = load(sys.argv[2], sys.argv[3], sys.argv[4])
    print("config: " + "  ".join("%s=%s" % (k, v) for k, v in sorted(cfg.items())))
    print()
    t0 = time.time()
    pairs = close_pairs(xyz, cfg["r_min"])
    cl, missed = cliques(xyz, cfg["r_min"], pairs)
    p = Problem(xyz, ranges, L, A, cfg, cliq=cl)
    print("the problem (built in %.1fs)" % (time.time() - t0))
    describe(p, pairs)
    print("   %d variables" % p.nvar)
    print()
    print("EXPLORATORY on the development set: a code test, not a result. Nothing here is quoted as one.")
    print()
    if mode == "solve":
        print("the %d best distinct designs at n = %d, B = %.2f" % (cfg["pool"], cfg["n"], cfg["budget"]))
        pool = p.pool(cfg["n"], cfg["budget"], cfg["pool"])
        if not pool:
            print("   no feasible design")
            return 1
        worst = 0.0
        for k, s in enumerate(pool):
            print("  [%d]" % k)
            worst = max(worst, show(s, p))
        print("\n   objective spread across the pool: %.5f kcal/mol" % (pool[-1]["objective"] - pool[0]["objective"]))
        print("   worst rule or arithmetic violation anywhere in the pool: %.1e" % worst)
        return 0
    # mode == "time": the solve time and gap as n grows (decisions, section 4)
    print("time and gap against n, B = %.2f" % cfg["budget"])
    print("    n   objective        mean   semi-dev   solver    wall(s)  how    certified  check")
    for n in range(3, cfg["n"] + 1):
        s, info = p.best(n, cfg["budget"])
        if s is None:
            print("   %2d   no design: %s (%.1fs)" % (n, info["message"], info["wall"]))
            continue
        w = max(violations(p, s).get(k, 0.0) for k in ("cardinality", "charge_range", "spacing", "d_f", "objective"))
        cg = s.get("certified_gap", float("nan"))
        print("   %2d  %+10.4f  %+10.4f  %9.4f  %.1e  %8.1f  %-5s %8.4f  (worst violation %.1e)"
              % (n, s["objective"], s["mean"], s["semidev"], s["gap"], s["wall"],
                 s.get("strategy", "full"), cg, w))
    print()
    print("   'certified gap' is against the relaxation of the WHOLE grid, so it bounds how far the design can be")
    print("   from the best design this score and these rules admit anywhere on it. 'full' rows searched every site")
    print("   and closed that gap themselves; 'lp' rows searched the sites the relaxation used.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
