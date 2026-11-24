#!/usr/bin/env python3
"""Stage 6, DIAGNOSTIC 56: the search and its inputs. The committed optimiser (cs_opt_milp.py) places bare point charges of
at most 0.25 e on the committed grid of frame 42436; the best three distinct designs at each number of charges are then
measured at the working level and at CAM-B3LYP-D3BJ/def2-TZVPD, the dianion and the monoanion at R and TS.

THE SEARCH. Sites, admitted signs and the spacing are those of sites_f42436.tsv and Stage 5. The first-order coefficients
are the committed A_f42436_d33a.tsv (one frame, so the semi-deviation term of the score is zero). L is zero: the charges
are bare points. Every admitted range is cut to |q| <= CAP_Q; a site that then admits nothing is left out. The budget
is set where it cannot bind (n * CAP_Q^2). One rule is added to the committed program's rules: the charges of a design
sum to zero or more (H9 of D43). A design is taken only if the solver proves it optimal.

THE INPUTS. Every input is a committed input, unchanged: the working-level ones of D43, D31 and D40, the def2-TZVPD ones
of D44. The job script is D53's committed one with its names, its runs and its caps changed.

usage: cs6_d56_prepare.py <results repo> <code repo> <out dir> [<hpc dir>]     runs the search and writes the d56 directory
"""
import hashlib, importlib.util, os, sys

S6 = "phase2.2/cs_stage6"
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/d56"
P47 = "phase2.2/cs6_d47_prepare.py"
OPT = "phase2.2/cs_opt_milp.py"
CAP_Q = 0.25                              # e: the largest charge at a site (AMENDMENT 11)
NS = (3, 4, 5, 6, 8, 10, 12, 15, 20)      # the numbers of charges searched
POOL = 3                                  # distinct designs kept at each n, best first
R_MIN = 3.25                              # A: the spacing of Stage 5
DISTINCT_R = 1.5                          # A: decision O7
SOLVE_LIMIT = 3600.0                      # s: a search that is not proven optimal within this stops the program
LEVELS = (("SV", "CAM-B3LYP-D3BJ/def2-SVPD", 402, 1800), ("TZ", "CAM-B3LYP-D3BJ/def2-TZVPD", 682, 5000))   # key, name, basis functions, cap of one run in s
SRC = {("SV", "di"): "d43/B_n11_q0.25_di_R/job.inp", ("SV", "mo"): "d43/B_n11_q0.25_mo_R/job.inp", ("SV", "ctl"): "d31/ctl_R/job.inp", ("SV", "bare"): "d40/B_bare_R/job.inp",
       ("TZ", "di"): "d44/CAMTZ_B_n02_q0.50_di_R/job.inp", ("TZ", "mo"): "d44/CAMTZ_B_n02_q0.50_mo_R/job.inp", ("TZ", "ctl"): "d44/CAMTZ_C_ctl_R/job.inp", ("TZ", "bare"): "d44/CAMTZ_C_bare_R/job.inp"}
PBS53 = "d53/run_K3.pbs"
CAP53 = 10200
LIMIT_H = 24.0
STATES = ("R", "TS")


def stop(m):
    sys.stderr.write("STOP: %s\n" % m)
    sys.exit(1)


def sha_text(t):
    return hashlib.sha256(t.encode()).hexdigest()


def load(path, name):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def cut_ranges(ranges, cap=CAP_Q):
    """each admitted range cut to |q| <= cap; an empty range is dropped"""
    out = []
    for rr in ranges:
        mine = []
        for lo, hi in rr:
            lo2, hi2 = max(lo, -cap), min(hi, cap)
            if lo2 <= hi2 + 1e-12:
                mine.append((lo2, hi2))
        out.append(mine)
    return out


def search(s6, code, a_col=None, quiet=True):
    """the designs, best first within each n: dict(name, n, k, sites, q, score, net). a_col replaces the committed
    coefficients ONLY in rehearsals, with invented numbers."""
    import numpy as np
    import io, contextlib
    opt = load(os.path.join(code, OPT), "opt_for_d56")
    p47 = load(os.path.join(code, P47), "p47_for_d56")
    xyz, layer, ranges = opt.read_sites(os.path.join(s6, p47.SITES))
    a = np.array(p47.read_col(os.path.join(s6, p47.AMAT)) if a_col is None else a_col, dtype=float)
    if len(a) != len(xyz) or len(xyz) != p47.NSITE:
        stop("the site table and the coefficient table do not both hold %d rows." % p47.NSITE)
    cut = cut_ranges(ranges)
    keep = [i for i, rr in enumerate(cut) if rr]
    cfg = dict(opt.CONFIG, r_min=R_MIN, distinct_r=DISTINCT_R, lam=0.0, time_limit=SOLVE_LIMIT, mip_gap=1e-9)

    class Net(opt.Problem):
        """the committed problem with one more row: the charges sum to zero or more"""
        def constraints(self, n, budget, cuts=()):
            from scipy.optimize import LinearConstraint
            from scipy.sparse import csr_matrix, vstack
            con, nrow = opt.Problem.constraints(self, n, budget, cuts)
            row = csr_matrix(([1.0] * self.R, ([0] * self.R, [self.i_q + k for k in range(self.R)])), shape=(1, self.nvar))
            return LinearConstraint(vstack([con.A, row]).tocsc(), np.append(con.lb, 0.0), np.append(con.ub, np.inf)), nrow + 1

    sub_xyz = xyz[keep]
    A = a[keep].reshape(-1, 1)
    L = np.zeros_like(A)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf if quiet else sys.stdout):
        p = Net(sub_xyz, [cut[i] for i in keep], L, A, cfg)
    out = []
    for n in NS:
        cuts = []
        for k in range(1, POOL + 1):
            with contextlib.redirect_stdout(buf if quiet else sys.stdout):
                sol, info = p.solve(n, n * CAP_Q ** 2 + 1.0, cuts)
            if sol is None or sol["status"] != 0:
                stop("the search at n = %d, design %d, was not proven optimal: %s" % (n, k, info.get("message")))
            sites = [keep[i] for i in sol["sites"]]
            q = [round(v, 4) for v in sol["q"]]
            pairs = sorted((s, v) for s, v in zip(sites, q) if abs(v) > 1e-9)
            if len(pairs) != n:
                stop("the design %d at n = %d holds a charge of zero." % (k, n))
            if any(abs(v) > CAP_Q + 1e-9 for s, v in pairs) or sum(v for s, v in pairs) < -1e-6:
                stop("the design %d at n = %d breaks a rule of the search." % (k, n))
            for i in range(len(pairs)):
                for j in range(i):
                    if p47.dist(tuple(xyz[pairs[i][0]]), tuple(xyz[pairs[j][0]])) < R_MIN - 1e-9:
                        stop("two charges of design %d at n = %d are closer than %.2f A." % (k, n, R_MIN))
            out.append(dict(name="S_n%02d_k%d" % (n, k), n=n, k=k, sites=[s for s, v in pairs], q=[v for s, v in pairs],
                            xyz=[tuple(float(t) for t in xyz[s]) for s, v in pairs], score=float(sum(a[s] * v for s, v in pairs)), net=float(sum(v for s, v in pairs))))
            cuts.append(p.near(sol["sites"]))
    return out, dict(nsite=len(xyz), nkept=len(keep))


def runs_of(designs):
    """every run, in the order made within its job: dict(run, job, level, kind, design, state)"""
    out = []
    for lv, _, _, _ in LEVELS:
        for kind in ("ctl", "bare"):
            for st in STATES:
                out.append(dict(run="%s_C_%s_%s" % (lv, kind, st), job="C0", level=lv, kind=kind, design=None, state=st))
    for jn, n in enumerate(NS):
        mine = [d for d in designs if d["n"] == n]
        for lv, kinds in (("SV", ("di", "mo")), ("TZ", ("di",)), ("TZ", ("mo",))):      # the short runs first, the longest last
            for kind in kinds:
                for d in mine:
                    for st in STATES:
                        out.append(dict(run="%s_%s_%s_%s" % (lv, d["name"], kind, st), job="N%d" % (jn + 1), level=lv, kind=kind, design=d["name"], state=st))
    if len(set(r["run"] for r in out)) != len(out):
        stop("a run is listed twice.")
    return out


def jobs_of():
    return [("C0", "cm_d56_c0", "the eight controls without a charge, both levels")] + [
        ("N%d" % (j + 1), "cm_d56_n%d" % (j + 1), "the %d designs of %d charges, at both levels" % (POOL, n)) for j, n in enumerate(NS)]


def pbs_text(s6, job, runs, hpc):
    """D53's committed job script with its names and its runs changed, and the cap of a run set by its level"""
    jk, jn, what = job
    cap = dict((lv, c) for lv, _, _, c in LEVELS)
    src = open(os.path.join(s6, PBS53)).read().split("\n")
    k = [i for i, l in enumerate(src) if l.startswith("run ")]
    if not k or src[k[0] - 1] != "}":
        stop("%s is not the job script this program expects." % PBS53)
    head = src[:k[0]]
    mine = [r for r in runs if r["job"] == jk]
    worst = sum(cap[r["level"]] + 10 for r in mine) / 3600.0
    if worst >= LIMIT_H:
        stop("job %s could run for %.1f hours at its caps." % (jk, worst))
    cap53, capnew = 'timeout -k 10 %d "$ORCA/orca" job.inp' % CAP53, 'timeout -k 10 "$2" "$ORCA/orca" job.inp'
    if sum(cap53 in l for l in head) != 1:
        stop("%s does not hold the cap of %d s." % (PBS53, CAP53))
    out = []
    for l in head:
        l = l.replace(cap53, capnew).replace("cm_d53_k3", jn).replace("cs_stage6/d53", "cs_stage6/d56").replace("run_K3", "run_%s" % jk).replace(".running_K3", ".running_%s" % jk)
        if l.startswith("# Stage 6 DIAGNOSTIC 53"):
            l = "# Stage 6 DIAGNOSTIC 56: the search on the committed grid, each design measured at two levels."
        elif l.startswith("# Job K3 of"):
            l = "# Job %s of %d: %s. It makes %d runs in the order below (name, cap in s); a run that fails does not stop the ones after it." % (jk, len(jobs_of()), what, len(mine))
        out.append(l)
    if any("d53" in l or "K3" in l for l in out):
        stop("a line of D53's job script was not changed.")
    out += ["run %s %d" % (r["run"], cap[r["level"]]) for r in mine]
    out += ['echo "end=$(date)"',
            'n=0; for r in %s; do grep -q \'ORCA TERMINATED NORMALLY\' "$W/$r/job.out" 2>/dev/null && n=$((n+1)); done' % " ".join(r["run"] for r in mine),
            'echo "finished: $n of %d"' % len(mine), ""]
    return "\n".join(out)


def read_designs(path, s6, code):
    """the designs of a committed DESIGNS.tsv, each checked against the rules of the search"""
    opt = load(os.path.join(code, OPT), "opt_for_d56r")
    p47 = load(os.path.join(code, P47), "p47_for_d56r")
    xyz, layer, ranges = opt.read_sites(os.path.join(s6, p47.SITES))
    cut = cut_ranges(ranges)
    a = p47.read_col(os.path.join(s6, p47.AMAT))
    L = [l.split("\t") for l in open(path).read().splitlines()]
    out = []
    for r in L[1:]:
        d = dict(zip(L[0], r))
        sites = [int(v) for v in d["sites"].split(",")]; q = [float(v) for v in d["charges_e"].split(",")]
        n, k = int(d["n"]), int(d["rank_at_n"])
        bad = []
        if d["design"] != "S_n%02d_k%d" % (n, k) or len(sites) != n or len(q) != n or len(set(sites)) != n or n not in NS or not 1 <= k <= POOL:
            bad.append("its name, its number of charges or its rank")
        if any(not any(lo - 1e-9 <= v <= hi + 1e-9 for lo, hi in cut[s]) or abs(v) < 1e-9 for s, v in zip(sites, q)):
            bad.append("a charge its site does not admit, or one above %.2f e" % CAP_Q)
        if sum(q) < -1e-6:
            bad.append("charges that sum to less than zero")
        if any(p47.dist(tuple(xyz[sites[i]]), tuple(xyz[sites[j]])) < R_MIN - 1e-9 for i in range(n) for j in range(i)):
            bad.append("two charges closer than %.2f A" % R_MIN)
        if abs(sum(a[s] * v for s, v in zip(sites, q)) - float(d["first_order_score_kcal"])) > 1e-5:
            bad.append("a first-order score that is not the committed coefficients times its charges")
        if bad:
            stop("the design %s of DESIGNS.tsv has %s." % (d["design"], "; ".join(bad)))
        out.append(dict(name=d["design"], n=n, k=k, sites=sites, q=q, xyz=[tuple(float(t) for t in xyz[s]) for s in sites], score=float(d["first_order_score_kcal"]), net=float(sum(q))))
    if [(d["n"], d["k"]) for d in out] != [(n, k) for n in NS for k in range(1, POOL + 1)]:
        stop("DESIGNS.tsv does not hold %d designs at each of the %d numbers of charges, in order." % (POOL, len(NS)))
    for n in NS:
        sc = [d["score"] for d in out if d["n"] == n]
        if any(sc[i] > sc[i + 1] + 1e-9 for i in range(len(sc) - 1)):
            stop("the designs of n = %d are not in the order of their scores." % n)
    return out, dict(nsite=len(xyz), nkept=sum(1 for rr in cut if rr))


def build(s6, code, hpc=HPC, a_col=None, designs_from=None):
    """(files of d56 as {relative path: text}, designs, runs). designs_from: a committed DESIGNS.tsv to form the inputs
    again from, without searching"""
    p47 = load(os.path.join(code, P47), "p47_for_d56b")
    for f in [p47.GEOM["R"], p47.GEOM["TS"], p47.SITES, p47.AMAT, PBS53] + sorted(SRC.values()):
        if not os.path.isfile(os.path.join(s6, f)):
            stop("%s is missing." % f)
    designs, info = search(s6, code, a_col) if designs_from is None else read_designs(designs_from, s6, code)
    runs = runs_of(designs)
    F = {}
    inp = dict((k, open(os.path.join(s6, v)).read()) for k, v in SRC.items())
    for lv, name, nb, cap in LEVELS:
        want = p47.LINE["CAMTZ" if lv == "TZ" else "CAMSV"]
        if any(inp[(lv, kind)].split("\n")[0] != want for kind in ("di", "mo", "ctl", "bare")):
            stop("a committed input does not carry the method line of %s." % name)
    geo = dict((st, open(os.path.join(s6, p47.GEOM[st])).read()) for st in STATES)
    dmap = dict((d["name"], d) for d in designs)
    for d in designs:
        d["pc"] = p47.pc_text(list(zip(d["q"], d["xyz"])))
    for r in runs:
        F["%s/job.inp" % r["run"]] = inp[(r["level"], r["kind"])]
        F["%s/substrate.xyz" % r["run"]] = geo[r["state"]]
        if r["design"]:
            F["%s/pc.pc" % r["run"]] = dmap[r["design"]]["pc"]
    F["DESIGNS.tsv"] = "design\tn\trank_at_n\tfirst_order_score_kcal\tsum_q_e\tpositive_e\tsites\tcharges_e\tpc_sha256\n" + "".join(
        "%s\t%d\t%d\t%.6f\t%+.4f\t%.4f\t%s\t%s\t%s\n" % (d["name"], d["n"], d["k"], d["score"], d["net"], sum(v for v in d["q"] if v > 0), ",".join(str(s) for s in d["sites"]),
                                                    ",".join("%+.4f" % v for v in d["q"]), sha_text(d["pc"])) for d in designs)
    nbf = dict((lv, nb) for lv, _, nb, _ in LEVELS)
    F["RUNS.tsv"] = "run\tjob\tlevel\tkind\tdesign\tstate\tbasis_functions\n" + "".join(
        "%s\t%s\t%s\t%s\t%s\t%s\t%d\n" % (r["run"], r["job"], r["level"], r["kind"], r["design"] or "-", r["state"], nbf[r["level"]]) for r in runs)
    F["SEARCH.txt"] = "".join("%s = %s\n" % kv for kv in (("largest charge at a site, e", CAP_Q), ("numbers of charges", ", ".join(str(n) for n in NS)), ("designs kept at each", POOL),
                              ("spacing, A", R_MIN), ("distinct-design radius, A", DISTINCT_R), ("sites of the grid", info["nsite"]), ("sites that admit a charge of this size", info["nkept"]),
                              ("the charges of a design sum to", "zero or more"), ("every design", "proven optimal by the solver")))
    for job in jobs_of():
        F["run_%s.pbs" % job[0]] = pbs_text(s6, job, runs, hpc)
    srcs = sorted(set([p47.GEOM["R"], p47.GEOM["TS"], p47.SITES, p47.AMAT, PBS53] + list(SRC.values())))
    F["SOURCES.sha256"] = "".join("%s  %s\n" % (hashlib.sha256(open(os.path.join(s6, f), "rb").read()).hexdigest(), f) for f in srcs)
    F["INPUT_CHECKSUMS.sha256"] = "".join("%s  %s\n" % (sha_text(F[k]), k) for k in sorted(F))
    return F, designs, runs


def main():
    if len(sys.argv) not in (4, 5):
        stop("usage: cs6_d56_prepare.py <results repo> <code repo> <out dir> [<hpc dir>]")
    res, code, out = (os.path.abspath(x) for x in sys.argv[1:4])
    if os.path.exists(out):
        stop("%s exists already; this program overwrites nothing." % out)
    F, designs, runs = build(os.path.join(res, S6), code, sys.argv[4] if len(sys.argv) == 5 else HPC)
    for k in sorted(F):
        p = os.path.join(out, k)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as f:
            f.write(F[k])
    print("d56: %d designs, %d runs in %d jobs, %d files" % (len(designs), len(runs), len(jobs_of()), len(F) + 0))


if __name__ == "__main__":
    main()
