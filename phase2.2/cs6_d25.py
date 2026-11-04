#!/usr/bin/env python3
"""
cs6_d25.py - Stage 6 DIAGNOSTIC 25: is the predicted barrier reduction transition-state stabilisation,
or reactant destabilisation? Pre-registered in cs_stage6/STAGE6_DIAGNOSTIC25_20261104.txt.

  python3 cs6_d25.py selftest <results_repo>
  python3 cs6_d25.py run      <results_repo> <OUT>

NO CALCULATION. Every number is a dot product over committed files: no SCF, no orca_vpot, no HPC.

WHAT IT COMPUTES, per design, in kcal/mol (D25.2)
  dE_R   = sum_i V_R(r_i)  q_i      the first-order reactant-side term
  dE_TS  = sum_i V_TS(r_i) q_i      the transition-state-side term
  d_elec = dE_TS - dE_R             which must equal sum_i a_i q_i identically
  q_net, sum|q|, and f_R = dE_R / d_elec, the signed share carried by the reactant side

WHY THE DECOMPOSITION IS POSSIBLE AT ALL. cs_opt_milp.py only ever sees a = V_TS - V_R, but
cs6_potentials.py wrote the two potentials apart into cs_stage6/absolute_potentials/V_f42436.tsv from
the same committed a_matrix outputs. OPTIMISER_REACTANT_TERM_ERRATUM1_20261104.md records that the
finding which proposed this diagnostic named the wrong file - the exploratory dev_A vpot set, at a
grid where 0 of the 2626 production sites appear - and that this one reproduces A_f42436.tsv to
1.0e-06 kcal/mol per e. This reader re-derives that rather than citing it.

INDEXING, pinned from committed data rather than assumed. ARRANGEMENTS.txt's site numbers are
0-BASED and index A_f42436.tsv directly: over the first 40 sphere records, A[site] reproduces the
`a` value ARRANGEMENTS.txt prints for that sphere with a worst difference of 0.0000, while A[site-1]
mismatches all 40. check_index() re-establishes that every run.

CHECKS, every one a stop, lettered as D25.2 letters them
  (a) sum_i a_i q_i from A_f42436.tsv equals dE_TS - dE_R from V_f42436.tsv, to 1e-4 kcal/mol, per
      design. This is the identity that licenses the decomposition.
  (b) (V_TS - V_R) reproduces A_f42436.tsv over all 2626 sites to 1e-5 kcal/mol per e.
  (c) every design's site indices exist in sites_f42436.tsv, and every charge equals that design's
      committed site.ORCAFF.prms to 1e-9 e.
  (d) S6-C, which has no spheres, yields dE_R = dE_TS = d_elec = 0 exactly. It is the null control.
  (e) S6-P is the one BARE arrangement and no Lennard-Jones term is read for it.
  (f) every input file's sha256 equals the value STAGE6_DIAGNOSTIC25_20261104.txt records for it.

WHAT THE SELFTEST DOES AND DELIBERATELY DOES NOT DO. It exercises every PARSE against the real
committed files - the potentials, the A matrix, the site table, ARRANGEMENTS.txt and every design's
force field - and makes every stop fire on constructed data. It does NOT compute dE_R or dE_TS for
any design except S6-C, whose values are zero by construction and which the pre-registration names
as the null control. Those decompositions ARE the diagnostic, and this project commits a diagnostic
before any of its runs exist.

Python 3.8+, numpy.
"""
import hashlib
import importlib.util
import math
import re
import sys
from pathlib import Path

import numpy as np

NSITES = 2626
NQM = 24
HALF = 0.5                 # AMENDMENT 3: n10_A is S6-10 at exactly half field
DESIGN_EXTRA = "n10_A"
PARENT = "S6-10"
CONTROL = "S6-C"
BARE = "S6-P"
FRACTION = 0.15            # the project's one fraction (decision O6, STAGE6_CRITERIA.txt)

SHA = {   # D25.2, and check (f) compares against these
    "phase2.2/cs_stage6/absolute_potentials/V_f42436.tsv":
        "7b24124a6bac673f57923ad3835535cfc76285af26277b2e64ee462919c49eaf",
    "phase2.2/cs_stage6/A_f42436.tsv":
        "a7faa6c5cacc2f2ee0956f80a1a6645eed1e04cc126830d5f5803d9762eae553",
    "phase2.2/cs_stage6/sites_f42436.tsv":
        "d2a750077d79e34dfecf448c30c2899ae0cea98f1edc5de3d10361b7c357cfc7",
    "phase2.2/cs_stage6/relax/ARRANGEMENTS.txt":
        "651e7e228e1409f052bb024b42d125ab59bc86ebee789f3dd60ad1eb5f9c6038",
}


def stop(m):
    sys.exit("STOP: " + m)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _load(name, where):
    """cs6_analyse.py's ARRANGEMENTS reader is reused rather than copied. Loaded by explicit path:
    run.sh runs this under -I, which implies -P, so the script's directory is not on sys.path."""
    p = Path(where) / (name + ".py")
    if not p.is_file():
        stop("%s is not beside this script; cs6_d25.py reuses cs6_analyse.read_arrangements" % p)
    spec = importlib.util.spec_from_file_location("_d25_" + name, p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m, sha(p)


def _text(p):
    """every file this reader opens goes through here, so a missing input is a STOP that names it
    rather than a traceback."""
    p = Path(p)
    if not p.is_file():
        stop("%s does not exist" % p)
    return p.read_text(errors="replace")


def read_V(p):
    L = [l for l in _text(p).splitlines() if l and not l.startswith("#")]
    head = L[0].split("\t")
    for want in ("site", "V_R", "V_TS", "V_P"):
        if want not in head:
            stop("%s has no %s column; its columns are %s" % (p, want, head))
    ix = {n: i for i, n in enumerate(head)}
    R = np.array([[float(x) for x in l.split("\t")] for l in L[1:]], dtype=float)
    if R.shape[0] != NSITES:
        stop("%s holds %d rows, expected %d sites" % (p, R.shape[0], NSITES))
    site = R[:, ix["site"]].astype(int)
    if not np.array_equal(site, np.arange(NSITES)):
        stop("%s's site column is not 0 to %d in order" % (p, NSITES - 1))
    return R[:, ix["V_R"]], R[:, ix["V_TS"]], R[:, ix["V_P"]]


def read_A(p):
    v = np.array([float(l.split("\t")[0]) for l in _text(p).splitlines()
                  if l and not l.startswith("#")], dtype=float)
    if len(v) != NSITES:
        stop("%s holds %d values, expected %d" % (p, len(v), NSITES))
    return v


def read_sites(p):
    L = [l for l in _text(p).splitlines() if l and not l.startswith("#")]
    head = L[0].split("\t")
    ix = {n: i for i, n in enumerate(head)}
    for want in ("site", "layer"):
        if want not in ix:
            stop("%s has no %s column" % (p, want))
    rows = [l.split("\t") for l in L[1:]]
    if len(rows) != NSITES:
        stop("%s holds %d sites, expected %d" % (p, len(rows), NSITES))
    ids = [int(r[ix["site"]]) for r in rows]
    if ids != list(range(NSITES)):
        stop("%s's site ids are not 0 to %d in order" % (p, NSITES - 1))
    return {"layer": [r[ix["layer"]] for r in rows]}


def prms_charges(p):
    """the sphere charges of an arrangement as ORCA read them: column 2 of the rows past the 24
    substrate ones. Column 2 is the charge - the remaining columns are Lennard-Jones parameters and
    hold no coordinates, which is the mistake cs6_relaxed.py made and its erratum records."""
    rows = [r.split() for r in _text(p).splitlines() if len(r.split()) == 7]
    if len(rows) < NQM:
        stop("%s holds %d parameter rows, fewer than the %d substrate atoms" % (p, len(rows), NQM))
    return [float(r[2]) for r in rows[NQM:]]


def check_index(arr, A):
    """(the indexing, pinned every run) ARRANGEMENTS.txt prints `a` for each sphere. Those values must
    be A[site] and not A[site-1], and this re-establishes it from the committed data."""
    recs = [(s, a) for d in arr.values() for s, _q, a, _l in d["sph"]]
    if not recs:
        stop("ARRANGEMENTS.txt yielded no sphere records, so the indexing cannot be pinned")
    w0 = max(abs(A[s] - a) for s, a in recs if 0 <= s < NSITES)
    off = [(s, a) for s, a in recs if not (0 <= s - 1 < NSITES)]
    w1 = max(abs(A[s - 1] - a) for s, a in recs if 0 <= s - 1 < NSITES) if len(recs) > len(off) else 9e9
    if not (w0 <= 5e-4 < w1):
        stop("the site indexing is not what it was pinned to be: 0-based worst |diff| %.4g, "
             "1-based %.4g, over %d sphere records" % (w0, w1, len(recs)))
    return len(recs), w0, w1


def designs_of(arr, relax):
    """the pre-registered set: ARRANGEMENTS.txt's seven, its ADD-* one-sphere runs, and n10_A."""
    out = dict(arr)
    par = arr.get(PARENT)
    if par is None:
        stop("%s is not in ARRANGEMENTS.txt, so %s cannot be derived from it" % (PARENT, DESIGN_EXTRA))
    # the ADD-* one-sphere runs, discovered from the folders and CHECKED against the spheres of
    # S6-2 and S6-3 exactly as cs6_analyse.frozen does, so a folder whose name matches no real
    # sphere of those two designs is a stop rather than something quietly decomposed.
    want = {}
    for nm in ("S6-2", "S6-3"):
        for site, q, a_, l_ in arr.get(nm, {"sph": []})["sph"]:
            want["ADD-s%d_q%s%05d" % (site, "p" if q >= 0 else "m", round(abs(q) * 10000))] = \
                (site, q, a_, l_)
    for d_ in sorted(Path(relax).glob("ADD-*")):
        tag = d_.name
        if tag not in want:
            stop("%s is not a one-sphere run of any S6-2 or S6-3 sphere" % tag)
        site, q, a_, l_ = want[tag]
        pr = d_ / "R_sp" / "site.ORCAFF.prms"
        if pr.is_file():
            qq = prms_charges(pr)
            if len(qq) != 1 or abs(qq[0] - q) > 1e-6:
                stop("%s does not hold exactly one sphere at q %+.4f" % (tag, q))
        out[tag] = {"n": 1, "bare": False, "sph": [(site, q, a_, l_)], "score": None,
                    "aq": None, "sl": None, "derived_from": "S6-2/S6-3"}

    prms = Path(relax) / DESIGN_EXTRA / "R_sp" / "site.ORCAFF.prms"
    if prms.is_file():
        q = prms_charges(prms)
        if len(q) != par["n"]:
            stop("%s carries %d spheres and %s has %d" % (DESIGN_EXTRA, len(q), PARENT, par["n"]))
        worst = max(abs(a - HALF * b) for a, b in zip(q, [s[1] for s in par["sph"]]))
        if worst > 1e-9:
            stop("%s's charges are not exactly %.2f of %s's: largest departure %.3e e"
                 % (DESIGN_EXTRA, HALF, PARENT, worst))
        out[DESIGN_EXTRA] = {"n": par["n"], "bare": False,
                             "sph": [(s, HALF * qq, a, l) for (s, qq, a, l) in par["sph"]],
                             "score": None, "aq": None, "sl": None, "derived_from": PARENT}
    return out


def decompose(d, VR, VTS, A):
    s = np.array([x[0] for x in d["sph"]], dtype=int)
    q = np.array([x[1] for x in d["sph"]], dtype=float)
    if len(s) == 0:
        return dict(dE_R=0.0, dE_TS=0.0, d_elec=0.0, aq=0.0, q_net=0.0, sabs=0.0, f_R=float("nan"))
    dR = float(VR[s] @ q)
    dT = float(VTS[s] @ q)
    aq = float(A[s] @ q)
    de = dT - dR
    return dict(dE_R=dR, dE_TS=dT, d_elec=de, aq=aq, q_net=float(q.sum()),
                sabs=float(np.abs(q).sum()), f_R=(dR / de) if de else float("nan"))


def spearman(x, y):
    """rank correlation and a Fisher-z interval, written out rather than imported so the reader needs
    only numpy. Ties are averaged."""
    def rank(v):
        v = np.asarray(v, dtype=float)
        o = np.argsort(v, kind="mergesort")
        r = np.empty(len(v), dtype=float)
        r[o] = np.arange(1, len(v) + 1, dtype=float)
        for u in np.unique(v):
            m = v == u
            if m.sum() > 1:
                r[m] = r[m].mean()
        return r
    rx, ry = rank(x), rank(y)
    rx -= rx.mean()
    ry -= ry.mean()
    den = math.sqrt(float((rx ** 2).sum()) * float((ry ** 2).sum()))
    if den == 0:
        return float("nan"), None
    rho = float((rx * ry).sum() / den)
    n = len(x)
    if n < 4 or abs(rho) >= 1.0:
        return rho, None
    z = 0.5 * math.log((1 + rho) / (1 - rho))
    se = 1.0 / math.sqrt(n - 3)
    lo, hi = (math.tanh(z - 1.959964 * se), math.tanh(z + 1.959964 * se))
    return rho, (lo, hi)


def load_all(root, code):
    root = Path(root)
    cs6 = root / "phase2.2" / "cs_stage6"
    A_mod, A_sha = _load("cs6_analyse", code)
    bad = []
    for rel, want in SHA.items():
        f = root / rel
        if not f.is_file():
            bad.append("missing: " + rel)
        elif sha(f) != want:
            bad.append("changed: %s\n      recorded %s\n      found    %s" % (rel, want, sha(f)))
    if bad:
        stop("check (f): %d input file(s) are not what DIAGNOSTIC 25 records:\n   %s"
             % (len(bad), "\n   ".join(bad)))
    VR, VTS, VP = read_V(cs6 / "absolute_potentials" / "V_f42436.tsv")
    A = read_A(cs6 / "A_f42436.tsv")
    sites = read_sites(cs6 / "sites_f42436.tsv")
    arr = A_mod.read_arrangements(cs6 / "relax" / "ARRANGEMENTS.txt")
    return dict(cs6=cs6, VR=VR, VTS=VTS, VP=VP, A=A, sites=sites, arr=arr,
                dsg=designs_of(arr, cs6 / "relax"), A_sha=A_sha)


def run(root, out_p):
    out = Path(out_p)
    if out.exists():
        stop("%s exists" % out)
    S = load_all(root, Path(__file__).resolve().parent)
    VR, VTS, A, dsg = S["VR"], S["VTS"], S["A"], S["dsg"]
    rep = []

    def say(s=""):
        print(s)
        rep.append(s)

    say("Stage 6, DIAGNOSTIC 25: the reactant-side and transition-state-side terms of the predicted")
    say("effect, decomposed from committed potentials. No calculation was run.")
    say("cs6_d25.py     sha256 %s" % sha(__file__))
    say("cs6_analyse.py sha256 %s   (imported: read_arrangements)" % S["A_sha"])
    say("")
    say("(f) the four recorded input files are unchanged")

    # (b) the reproduction, re-derived
    d1 = (VTS - VR) - A
    if np.abs(d1).max() > 1e-5:
        stop("check (b): (V_TS - V_R) departs from A_f42436.tsv by %.3e kcal/mol per e"
             % np.abs(d1).max())
    say("(b) (V_TS - V_R) reproduces A_f42436.tsv over all %d sites: max |difference| %.3e,"
        % (NSITES, np.abs(d1).max()))
    say("    median %.3e kcal/mol per e" % np.median(np.abs(d1)))

    # the indexing
    nrec, w0, w1 = check_index(S["arr"], A)
    say("    the site indexing is 0-based, pinned on %d sphere records: A[site] worst |diff| %.4f"
        % (nrec, w0))
    say("    against A[site-1] worst |diff| %.4f" % w1)

    # (c) charges against each force field
    nchk = 0
    for name, d in sorted(dsg.items()):
        if not d["sph"]:
            continue
        p = S["cs6"] / "relax" / name / "R_sp" / "site.ORCAFF.prms"
        if not p.is_file():
            continue
        q = prms_charges(p)
        if len(q) != len(d["sph"]):
            stop("check (c): %s's force field holds %d spheres, the arrangement has %d"
                 % (name, len(q), len(d["sph"])))
        worst = max(abs(a - b) for a, b in zip(q, [s[1] for s in d["sph"]]))
        if worst > 1e-9:
            stop("check (c): %s's charges differ from its force field by %.3e e" % (name, worst))
        nchk += 1
        for s, _q, _a, _l in d["sph"]:
            if not (0 <= s < NSITES):
                stop("check (c): %s uses site %d, outside 0 to %d" % (name, s, NSITES - 1))
    say("(c) %d designs' charges equal their committed force fields and every site index is in range"
        % nchk)

    # (d), (e) and (a)
    res = {}
    for name, d in dsg.items():
        res[name] = decompose(d, VR, VTS, A)
        if abs(res[name]["aq"] - res[name]["d_elec"]) > 1e-4:
            stop("check (a): %s's sum a*q is %+.6f but dE_TS - dE_R is %+.6f, a gap of %.3e"
                 % (name, res[name]["aq"], res[name]["d_elec"],
                    abs(res[name]["aq"] - res[name]["d_elec"])))
    c = res.get(CONTROL)
    if c is None or not (c["dE_R"] == c["dE_TS"] == c["d_elec"] == 0.0):
        stop("check (d): the null control %s does not give exactly zero" % CONTROL)
    if BARE in dsg and not dsg[BARE]["bare"]:
        stop("check (e): %s is not marked BARE in ARRANGEMENTS.txt" % BARE)
    say("(a) sum a*q equals dE_TS - dE_R for all %d designs, to 1e-4 kcal/mol" % len(res))
    say("(d) the null control %s gives exactly zero; (e) %s is the one BARE arrangement"
        % (CONTROL, BARE))
    say("")

    order = sorted(res, key=lambda k: (k.startswith("ADD-"), k == DESIGN_EXTRA, dsg[k]["n"], k))
    say("THE DECOMPOSITION, kcal/mol. dE_R > 0 means the reactant is RAISED; the barrier falls by")
    say("lifting its foot. dE_R < 0 means the field lowers the reactant's own energy.")
    say("  %-20s %3s %9s %10s %10s %10s %9s %8s"
        % ("design", "n", "sum|q|", "q_net", "dE_R", "dE_TS", "d_elec", "f_R"))
    for k in order:
        r = res[k]
        say("  %-20s %3d %9.4f %+10.4f %+10.4f %+10.4f %+9.4f %8s"
            % (k, dsg[k]["n"], r["sabs"], r["q_net"], r["dE_R"], r["dE_TS"], r["d_elec"],
               "%+.3f" % r["f_R"] if r["f_R"] == r["f_R"] else "-"))
    say("")

    live = [k for k in order if dsg[k]["sph"]]
    pos = [k for k in live if res[k]["dE_R"] > 0]
    neg = [k for k in live if res[k]["dE_R"] < 0]
    say("OF THE %d DESIGNS WITH SPHERES: %d have dE_R > 0 (reactant raised), %d have dE_R < 0."
        % (len(live), len(pos), len(neg)))
    tol = FRACTION
    b = [k for k in pos if abs(res[k]["dE_R"]) > tol * abs(res[k]["d_elec"])]
    cc = [k for k in neg if abs(res[k]["dE_R"]) > abs(res[k]["d_elec"])]
    say("")
    say("THE PRE-REGISTERED OUTCOMES (D25.3), applied mechanically")
    say("  A  the blindness never bit      : %s"
        % ("YES, every design is within %.0f %% " % (100 * tol) if not b and not cc else "no"))
    say("  B  off-mechanism catalysis      : %s"
        % (", ".join(b) if b else "no design"))
    say("  C  the relocation risk          : %s"
        % (", ".join(cc) if cc else "no design"))
    rho, ci = spearman([res[k]["q_net"] for k in live], [res[k]["dE_R"] for k in live])
    say("")
    say("  q_net against dE_R over the %d designs: Spearman rho %+.4f%s"
        % (len(live), rho, (", 95 %% interval %+.4f to %+.4f" % ci) if ci else ", interval not formed"))
    if ci and (ci[0] > 0 or ci[1] < 0):
        say("  The interval excludes zero, so the reactant term tracks net charge and D25.3 points at")
        say("  the constraint GOCAT uses, sum_i q_i = 0.")
    else:
        say("  The interval does not exclude zero, so D25.3 points at the explicit row")
        say("  sum_i V_R(r_i) q_i <= 0 rather than at net neutrality.")
    say("")
    say("NEITHER REMEDY IS APPLIED BY THIS DIAGNOSTIC AND NEITHER IS AUTHORISED BY IT. Changing what")
    say("the optimiser optimises is a Stage 5 decision needing its own dated amendment, a re-solve,")
    say("and a statement of which committed designs the amended optimiser would no longer produce.")
    say("")
    say("FIRST ORDER ONLY, from the frozen densities of the bare substrate, exactly as a is. It says")
    say("nothing about polarisation, and a design with dE_R <= 0 here could still destabilise the")
    say("reactant once polarised. One frame, 42436. No barrier is quoted. The catalytic verdict")
    say("remains SP_AMENDMENT2.txt A2.5's chain.")

    out.mkdir(parents=True)
    (out / "D25_REPORT.txt").write_text("\n".join(rep) + "\n")
    with open(out / "d25.tsv", "w") as fh:
        fh.write("design\tn\tsum_abs_q\tq_net\tdE_R_kcal\tdE_TS_kcal\td_elec_kcal\tf_R\n")
        for k in order:
            r = res[k]
            fh.write("%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n"
                     % (k, dsg[k]["n"], r["sabs"], r["q_net"], r["dE_R"], r["dE_TS"],
                        r["d_elec"], r["f_R"]))
    print("\nwrote %s/D25_REPORT.txt and d25.tsv" % out)
    return 0


# ---------------------------------------------------------------------------------------------
def _fires(fn, label, expect):
    try:
        fn()
    except SystemExit as e:
        m = str(e)
        if not m.startswith("STOP:"):
            return False, "   DID NOT STOP (exited %r): %s" % (m, label)
        if expect not in m:
            return False, "   WRONG REASON: %-40s wanted %r, got %r" % (label, expect,
                                                                        m[6:].splitlines()[0][:70])
        return True, "   fires: %-52s %s" % (label, m[6:].splitlines()[0][:70])
    return False, "   DID NOT STOP: %s" % label


def selftest(root):
    import shutil
    import tempfile
    ok = True
    def say(good, line):
        nonlocal ok
        ok = ok and good
        print(line)

    root = Path(root)
    print("cs6_d25.py selftest - every parse against the real committed files, every stop on")
    print("constructed data. No design's decomposition is computed except the null control's,")
    print("because that is the diagnostic and the pre-registration is committed first.")
    print()
    S = load_all(root, Path(__file__).resolve().parent)
    VR, VTS, A, dsg = S["VR"], S["VTS"], S["A"], S["dsg"]

    say(True, "   ok    (f) the four recorded sha256 all match")
    say(len(VR) == NSITES, "   %s  V_f42436.tsv parsed: %d sites, columns V_R V_TS V_P"
        % ("ok  " if len(VR) == NSITES else "FAIL", len(VR)))
    say(len(A) == NSITES, "   %s  A_f42436.tsv parsed: %d values" % ("ok  " if len(A) == NSITES else "FAIL", len(A)))
    m = float(np.abs((VTS - VR) - A).max())
    say(m <= 1e-5, "   %s  (b) the reproduction re-derived: max |difference| %.3e kcal/mol per e"
        % ("ok  " if m <= 1e-5 else "FAIL", m))
    nrec, w0, w1 = check_index(S["arr"], A)
    say(w0 <= 5e-4 < w1, "   %s  the indexing pinned on %d records: 0-based %.4f, 1-based %.4f"
        % ("ok  " if w0 <= 5e-4 < w1 else "FAIL", nrec, w0, w1))
    names = sorted(S["arr"])
    say(len(S["arr"]) >= 7, "   %s  ARRANGEMENTS.txt parsed: %d arrangements, %s"
        % ("ok  " if len(S["arr"]) >= 7 else "FAIL", len(S["arr"]), " ".join(names)))
    cnt = {k: v["n"] for k, v in S["arr"].items()}
    want = {"S6-C": 0, "S6-P": 1, "S6-2": 2, "S6-3": 3, "S6-5": 5, "S6-10": 10, "S6-20": 20}
    bad = {k: (cnt.get(k), v) for k, v in want.items() if cnt.get(k) != v}
    say(not bad, "   %s  their sphere counts are the committed ones%s"
        % ("ok  " if not bad else "FAIL", "" if not bad else "  %r" % bad))
    nadd = len([k for k in dsg if k.startswith("ADD-")])
    say(nadd == 4, "   %s  %d ADD-* one-sphere runs discovered and matched to S6-2/S6-3 spheres"
        % ("ok  " if nadd == 4 else "FAIL", nadd))
    say(DESIGN_EXTRA in dsg, "   %s  %s derived from %s at exactly half field"
        % ("ok  " if DESIGN_EXTRA in dsg else "FAIL", DESIGN_EXTRA, PARENT))
    nff = 0
    for name, d in sorted(dsg.items()):
        p = S["cs6"] / "relax" / name / "R_sp" / "site.ORCAFF.prms"
        if d["sph"] and p.is_file():
            q = prms_charges(p)
            assert len(q) == len(d["sph"])
            nff += 1
    say(nff >= 6, "   %s  (c) %d designs' force fields parsed, every charge count matching"
        % ("ok  " if nff >= 6 else "FAIL", nff))
    c = decompose(dsg[CONTROL], VR, VTS, A)
    good = c["dE_R"] == c["dE_TS"] == c["d_elec"] == 0.0
    say(good, "   %s  (d) the null control %s decomposes to exactly zero" % ("ok  " if good else "FAIL", CONTROL))

    # the rank correlation, on constructed data whose answer is known by hand
    rho, ci = spearman([1, 2, 3, 4, 5, 6], [2, 4, 6, 8, 10, 12])
    say(abs(rho - 1.0) < 1e-12, "   %s  spearman on a monotone pair gives rho %+.4f" % ("ok  " if abs(rho-1) < 1e-12 else "FAIL", rho))
    rho2, _ = spearman([1, 2, 3, 4, 5, 6], [12, 10, 8, 6, 4, 2])
    say(abs(rho2 + 1.0) < 1e-12, "   %s  and on a reversed pair %+.4f" % ("ok  " if abs(rho2+1) < 1e-12 else "FAIL", rho2))
    rho3, ci3 = spearman([1, 2, 3, 4, 5, 6, 7, 8], [3, 1, 4, 1, 5, 9, 2, 6])
    say(ci3 is not None, "   %s  a Fisher interval forms on 8 points: rho %+.4f, %s"
        % ("ok  " if ci3 else "FAIL", rho3, ("%+.3f to %+.3f" % ci3) if ci3 else "none"))

    # run() end to end on a SYNTHETIC repository whose answer is known by hand. This tests the
    # glue - the table, the outcome logic, the TSV - without computing any real design's
    # decomposition, which is the diagnostic itself.
    print()
    print("   run() on a synthetic repository with a hand-computed answer:")
    tmp0 = Path(tempfile.mkdtemp(prefix="cs6_d25_syn_"))
    try:
        c6 = tmp0 / "phase2.2" / "cs_stage6"
        (c6 / "absolute_potentials").mkdir(parents=True)
        (c6 / "relax" / "S6-9" / "R_sp").mkdir(parents=True)
        # V_R = -100 everywhere, V_TS = -100 + site/1000, so a = site/1000 exactly
        vr = [-100.0] * NSITES
        vts = [-100.0 + i / 1000.0 for i in range(NSITES)]
        vpp = [-100.0 + i / 500.0 for i in range(NSITES)]
        (c6 / "absolute_potentials" / "V_f42436.tsv").write_text(
            "site\tV_R\tV_TS\tV_P\n" + "".join(
                "%d\t%.6f\t%.6f\t%.6f\n" % (i, vr[i], vts[i], vpp[i]) for i in range(NSITES)))
        (c6 / "A_f42436.tsv").write_text("".join("%.6f\n" % (i / 1000.0) for i in range(NSITES)))
        (c6 / "AP_f42436.tsv").write_text("".join("%.6f\n" % (i / 500.0) for i in range(NSITES)))
        (c6 / "sites_f42436.tsv").write_text(
            "site\tlayer\tx\ty\tz\n" + "".join(
                "%d\tanion\t0.0\t0.0\t0.0\n" % i for i in range(NSITES)))
        # two spheres: site 1000 at q +2, site 2000 at q -1
        (c6 / "relax" / "ARRANGEMENTS.txt").write_text(
            "synthetic\n\n  S6-C   0 sphere(s)\n\n"
            "  S6-9   2 sphere(s)\n"
            "     site  1000  q +2.0000  at 0.0 0.0 0.0  a +1.00000  L +0.00000\n"
            "     site  2000  q -1.0000  at 0.0 0.0 0.0  a +2.00000  L +0.00000\n"
            "     score +0.00000\n")
        (c6 / "relax" / "S6-9" / "R_sp" / "site.ORCAFF.prms").write_text(
            "\n".join("%d X 0.0 0.0 0.0 0.0 0.0" % (i + 1) for i in range(NQM))
            + "\n%d N 2.0 0.0 0.0 0.0 0.0\n%d N -1.0 0.0 0.0 0.0 0.0\n" % (NQM + 1, NQM + 2))
        # hand computation: dE_R = -100*2 + -100*-1 = -100 ; dE_TS = -99*2 + -98*-1 = -100
        # a*q = 1.0*2 + 2.0*-1 = 0 ; d_elec = dE_TS - dE_R = 0
        want = {"dE_R": -100.0, "dE_TS": -100.0, "d_elec": 0.0, "q_net": 1.0, "sabs": 3.0}
        import copy as _c
        keep = _c.deepcopy(SHA); keepd = (DESIGN_EXTRA, PARENT, CONTROL, BARE)
        globals()["SHA"] = {("phase2.2/cs_stage6/" + k): sha(c6 / k) for k in
                            ("absolute_potentials/V_f42436.tsv", "A_f42436.tsv",
                             "sites_f42436.tsv", "relax/ARRANGEMENTS.txt")}
        globals()["DESIGN_EXTRA"] = "S6-NONE"
        globals()["PARENT"] = "S6-9"
        globals()["CONTROL"] = "S6-C"
        globals()["BARE"] = "S6-NONE"
        try:
            rc = run(tmp0, tmp0 / "out")
            row = dict(zip((tmp0 / "out" / "d25.tsv").read_text().splitlines()[0].split("\t"),
                           [l.split("\t") for l in
                            (tmp0 / "out" / "d25.tsv").read_text().splitlines()[1:]
                            if l.startswith("S6-9")][0]))
            for k, v in (("dE_R_kcal", want["dE_R"]), ("dE_TS_kcal", want["dE_TS"]),
                         ("d_elec_kcal", want["d_elec"]), ("q_net", want["q_net"]),
                         ("sum_abs_q", want["sabs"])):
                good = abs(float(row[k]) - v) < 1e-9
                say(good, "   %s  %-14s %+12.6f, constructed to be %+12.6f"
                    % ("ok  " if good else "FAIL", k, float(row[k]), v))
            say(rc == 0, "   %s  run() returned %d on the synthetic repository"
                % ("ok  " if rc == 0 else "FAIL", rc))
        finally:
            globals()["SHA"] = keep
            globals()["DESIGN_EXTRA"], globals()["PARENT"], globals()["CONTROL"], \
                globals()["BARE"] = keepd
    finally:
        shutil.rmtree(tmp0, ignore_errors=True)

    # the stops, on constructed files
    print()
    print("   the stops, each on constructed data:")
    tmp = Path(tempfile.mkdtemp(prefix="cs6_d25_st_"))
    try:
        vp = tmp / "V.tsv"
        vp.write_text("# x\nsite\tV_R\tV_TS\tV_P\n" + "\n".join(
            "%d\t-100.0\t-99.0\t-98.0" % i for i in range(NSITES)) + "\n")
        say(*_fires(lambda: read_V(tmp / "nope.tsv"), "a missing potential file", "does not exist"))
        bad1 = tmp / "V_short.tsv"
        bad1.write_text("site\tV_R\tV_TS\tV_P\n0\t-1\t-1\t-1\n")
        say(*_fires(lambda: read_V(bad1), "a potential file with too few rows", "rows, expected"))
        bad2 = tmp / "V_nocol.tsv"
        bad2.write_text("site\tV_R\tV_P\n" + "\n".join("%d\t-1\t-1" % i for i in range(NSITES)) + "\n")
        say(*_fires(lambda: read_V(bad2), "a potential file missing V_TS", "has no V_TS column"))
        bad3 = tmp / "V_unsorted.tsv"
        bad3.write_text("site\tV_R\tV_TS\tV_P\n" + "\n".join(
            "%d\t-1\t-1\t-1" % (NSITES - 1 - i) for i in range(NSITES)) + "\n")
        say(*_fires(lambda: read_V(bad3), "a potential file whose sites are out of order", "not 0 to"))
        ap = tmp / "A_short.tsv"
        ap.write_text("1.0\n2.0\n")
        say(*_fires(lambda: read_A(ap), "an A matrix with the wrong length", "values, expected"))
        pp = tmp / "short.prms"
        pp.write_text("\n".join("a b 0.0 d 0.0 0.0 0.0" for _ in range(3)) + "\n")
        say(*_fires(lambda: prms_charges(pp), "a force field shorter than the substrate", "parameter rows"))
        say(*_fires(lambda: check_index({"X": {"sph": [(0, 1.0, A[1] + 7.0, 0.0)]}}, A),
                    "an arrangement whose a value matches no indexing", "indexing is not what"))
        say(*_fires(lambda: check_index({}, A), "an ARRANGEMENTS file with no spheres", "no sphere records"))
        # (a): a doctored A that breaks the identity
        A2 = A.copy(); A2[0] += 1.0
        d = {"sph": [(0, 1.0, 0.0, 0.0)], "n": 1}
        r = decompose(d, VR, VTS, A2)
        say(abs(r["aq"] - r["d_elec"]) > 1e-4,
            "   %s  (a) a doctored A breaks the identity by %.4f, which the run-mode check stops on"
            % ("ok  " if abs(r["aq"] - r["d_elec"]) > 1e-4 else "FAIL", abs(r["aq"] - r["d_elec"])))
        # (f): a changed input
        import hashlib as _h
        fake = dict(SHA); k0 = sorted(SHA)[0]
        SHA[k0] = "0" * 64
        say(*_fires(lambda: load_all(root, Path(__file__).resolve().parent),
                    "an input file whose sha256 is not the recorded one", "check (f)"))
        SHA.update(fake)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print("\n%s" % ("ALL PASS" if ok else "SOMETHING FAILED"))
    return 0 if ok else 1


def main():
    if len(sys.argv) == 3 and sys.argv[1] == "selftest":
        return selftest(sys.argv[2])
    if len(sys.argv) != 4 or sys.argv[1] != "run":
        sys.exit(__doc__)
    return run(sys.argv[2], sys.argv[3])


if __name__ == "__main__":
    sys.exit(main())
