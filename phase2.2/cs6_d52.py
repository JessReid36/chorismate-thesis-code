#!/usr/bin/env python3
"""cs6_d52.py - the reader of Stage 6 DIAGNOSTIC 52: charges together, one electron first. Committed before any
job of D52 is submitted.

One electron in the field of charges of +0.25 e at positive sites of the design B_n11_q0.25, with the functions of
D50 (sets X2 and X3) at all six positive sites. B is minus the energy of the one electron, in eV: ORCA's final
energy less the nuclear repulsion it prints, which is expected to be zero. THE LIMIT of a system is the lowest level
in a large basis, computed beforehand in a second program (cs6_d52_expected.py); for one centre it is 13.6057 q^2.
Each run is read with the committed reader of outputs (cs6_d40.py), as D48's to D51's readers read theirs.

  ORCA AGAINST THE SECOND PROGRAM   ORCA's Hartree-Fock B within 0.001 eV of the second program's in the same
                           functions. Hartree-Fock is exact for one electron and uses no grid, so this tests that
                           the two programs were given the same system.
  THE FUNCTIONS            Hartree-Fock within 0.05 eV of the limit.
  THE GRID                 the working functional's B on the two highest rungs measured (IntAcc 40 and 20, or 20
                           and 12) agree to 0.005 eV: CONVERGED. THE LOWEST RUNG THAT SERVES is the first whose B,
                           and that of every rung above it, lies within 0.005 eV of the top one.
  THE FUNCTIONS FOR THE FUNCTIONAL   where the grid is converged in both sets, B in X2 within 0.05 eV of B in X3.

RULE A names the functions and the grid for the bare charges of the dianion: set X2 if it SERVES at all six systems,
set X3 if X2 does not and X3 does. RULE B names a second level: the candidate that lies nearest the limit at all six
systems, if it lies nearer than the working functional. RULE C states the limit of D51's pre-registration: three
times the working functional's largest distance from the limit, and not less than three times D51's 0.0631 eV.

The inputs are not taken on trust: every file of d52 is formed again by cs6_d52_prepare.py from the committed files
and compared, every output must print the committed input beside it, and the committed files this leans on must be
the committed ones. A gate that fails stops the reader before anything is written. D52 adopts nothing.

usage: python3 cs6_d52.py --res <results repository> --code <code repository> --prereg <the pre-registration>
                          --result <the result to write> --tsv <the table of runs to write>
"""
import argparse
import hashlib
import importlib.util
import os
import re
import subprocess
import sys
import textwrap

S = "phase2.2/cs_stage6"
PREREG = "STAGE6_DIAGNOSTIC52_20261123.txt"
A10 = "STAGE6_AMENDMENT10_20261123.txt"
D51PRE = "STAGE6_DIAGNOSTIC51_20261123.txt"
D51RES = "STAGE6_DIAGNOSTIC51_RESULT_20261123.txt"
D50PRE = "STAGE6_DIAGNOSTIC50_20261122.txt"
EV = 27.211386245988
TOL_ASKED = 0.0015     # an IntAcc ORCA prints against the one asked for
# quoted sentences: each must be in the document it is quoted from and in the pre-registration given
Q_A = "Charges together are measured at 0.25 e with functions at each positive charge"                                       # AMENDMENT 10, A10.4
Q_B = "a run that takes electrons with a gain in energy not larger than a few times that error is not settled at this level alone."   # D51's pre-registration
Q_C = "THE FUNCTIONS FOR THE CAPTURE TEST: SET X2, with auxiliary treatment A."                                              # D51's result
Q_D = "From 0.25 e up, where the grid is converged, the functional lies at most +0.0631 eV from exact (0.40 e)."             # D51's result
Q_E = "The next diagnostic then puts charges together first"                                                                 # D50's pre-registration
GUARDED = ["d37", "d43", "d44", "d50", "d51"]
DOCS = ["d51.tsv", "d51_rule.tsv", A10, D50PRE, D51PRE, D51RES]
W = 96
NAMES = {"HF": "Hartree-Fock", "CAM": "CAM-B3LYP", "WXV": "wB97X-V", "WMV": "wB97M-V", "RDS": "revDSD-PBEP86-D4/2021"}


def stop(m):
    sys.stderr.write("STOP  %s\n" % m)
    raise SystemExit(1)


def load(path, name):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def echoed(text):
    m = re.search(r"INPUT FILE\n=+\nNAME = [^\n]*\n((?:\|\s*\d+>[^\n]*\n)+)", text)
    if not m:
        return None
    ls = [re.sub(r"^\|\s*\d+> ?", "", l) for l in m.group(1).split("\n") if l]
    while ls and (not ls[-1].strip() or "****END OF INPUT****" in ls[-1]):
        ls.pop()
    return ls


def nbf(text):
    nb = sorted(set(int(v) for v in re.findall(r"Number of basis functions\s*\.+\s*(\d+)", text)))
    return nb[0] if len(nb) == 1 else None


def count(text):
    """the electrons ORCA's grid of the functional integrates: the last N(Total) it prints"""
    v = re.findall(r"(?m)^N\(Total\)\s*:\s*(-?\d+\.\d+)", text)
    return float(v[-1]) if v else None


def enuc(text):
    """the nuclear repulsion ORCA prints, hartree: every value it prints must be the same one"""
    v = sorted(set(float(x) for x in re.findall(r"Nuclear Repulsion\s+ENuc\s*\.+\s*(-?\d+\.\d+)", text)))
    return v[0] if len(v) == 1 else None


def grid(text):
    """(IntAcc and number of points of the grid of the functional), as ORCA prints them"""
    xc, pts = None, None
    k = text.find("DFT GRID GENERATION")
    kx = text.find("COSX GRID GENERATION")
    if k >= 0:
        part = text[k:kx] if kx > k else text[k:]
        m = re.search(r"General Integration Accuracy\s+IntAcc\s*\.+\s*(\d+\.\d+)", part)
        p = re.search(r"Total number of grid points\s*\.+\s*(\d+)", part)
        xc, pts = (float(m.group(1)) if m else None), (int(p.group(1)) if p else None)
    return xc, pts


def f(v, fmt="%+.4f", none="-"):
    return none if v is None else fmt % v


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    a.res, a.code, a.prereg = (os.path.abspath(x) for x in (a.res, a.code, a.prereg))
    S6 = os.path.join(a.res, S)
    D = os.path.join(S6, "d52")
    paths = dict((k, os.path.join(a.code, "phase2.2", fn)) for k, fn in (
        ("d40", "cs6_d40.py"), ("d33b", "cs6_d33b.py"), ("p48", "cs6_d48_prepare.py"), ("p50", "cs6_d50_prepare.py"), ("p52", "cs6_d52_prepare.py"),
        ("exp", "cs6_d52_expected.py")))
    for k, p in paths.items():
        if not os.path.isfile(p):
            stop("the committed script %s is missing." % p)
    for fn in [a.prereg, os.path.join(D, "RUNS.tsv")] + [os.path.join(S6, x) for x in DOCS]:
        if not os.path.isfile(fn):
            stop("%s is missing." % fn)
    tsv2, tsv3 = os.path.join(os.path.dirname(a.tsv), "d52_systems.tsv"), os.path.join(os.path.dirname(a.tsv), "d52_rules.tsv")
    for fn in (a.result, a.tsv, tsv2, tsv3):
        if os.path.exists(fn):
            stop("%s exists already; this reader overwrites nothing." % fn)
    d40, d33b = load(paths["d40"], "d52_cs6_d40"), load(paths["d33b"], "d52_cs6_d33b")
    p52 = load(paths["p52"], "d52_cs6_d52_prepare")
    gates = []

    def git(*args):
        try:
            q = subprocess.run(["git", "-C", a.res] + list(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            return q.returncode, q.stdout.strip()
        except OSError:
            return 127, ""

    # ------------------------------------------------------------ G52.3: the committed files this leans on
    guard = ["%s/%s" % (S, g) for g in GUARDED] + ["%s/%s" % (S, x) for x in DOCS]
    rc, st = git("status", "--porcelain", "--untracked-files=all", "--", *guard)
    if rc != 0:
        gates.append("G52.3 could not be checked: git did not run")
    elif st:
        gates.append("G52.3 a file under %s, or a committed table or document this reader names, is not the committed one" % ", ".join(GUARDED))
    rc, pcm = git("log", "--diff-filter=A", "--format=%H", "-1", "--", "%s/%s" % (S, PREREG))
    if rc != 0 or not pcm:
        gates.append("G52.3 the commit that added %s could not be found" % PREREG)
    else:
        rc1, b1 = git("rev-parse", "%s:%s/%s" % (pcm, S, PREREG))
        rc2, b2 = git("hash-object", a.prereg)
        if rc1 != 0 or rc2 != 0 or b1 != b2:
            gates.append("G52.3 the pre-registration given is not the committed %s" % PREREG)

    # ------------------------------------------------------------ G52.1 and G52.2: formed again from the committed files
    try:
        F, runs, p50, p48, ex, geo = p52.build(S6, a.code)
    except Exception as e:
        sys.stdout.write("\n".join("FAILED  " + g for g in gates + ["G52.1 the files of d52 could not be formed again from the committed files: %s: %s" % (type(e).__name__, e)]) + "\n")
        stop("%d gate(s) failed. Nothing was written." % (len(gates) + 1))
    psrc = os.path.join(D, "SOURCES.sha256")
    if not os.path.isfile(psrc) or open(psrc).read() != F["SOURCES.sha256"]:
        gates.append("G52.2 d52/SOURCES.sha256 is not the list formed again from the committed files")
    nsrc = 0
    for l in F["SOURCES.sha256"].splitlines():
        h, fn = l.split(None, 1)
        nsrc += 1
        p = os.path.join(S6, fn)
        if not os.path.isfile(p) or sha(p) != h:
            gates.append("G52.2 %s is not the file d52/SOURCES.sha256 lists" % fn)
    nsum = 0
    pic = os.path.join(D, "INPUT_CHECKSUMS.sha256")
    if not os.path.isfile(pic):
        gates.append("G52.1 d52/INPUT_CHECKSUMS.sha256 is missing")
    else:
        if open(pic).read() != F["INPUT_CHECKSUMS.sha256"]:
            gates.append("G52.1 d52/INPUT_CHECKSUMS.sha256 is not the list formed again here")
        for l in open(pic).read().splitlines():
            if l.strip():
                h, fn = l.split(None, 1)
                nsum += 1
                p = os.path.join(D, fn)
                if not os.path.isfile(p) or sha(p) != h:
                    gates.append("G52.1 %s is not the file listed in INPUT_CHECKSUMS.sha256" % fn)
    for k in sorted(F):
        p = os.path.join(D, k)
        if not os.path.isfile(p) or open(p).read() != F[k]:
            gates.append("G52.1 d52/%s is not the file formed again from the committed files" % k)

    # ------------------------------------------------------------ G52.4: the quoted sentences
    def has(path, sent):
        return " ".join(sent.split()) in " ".join(open(path, errors="replace").read().split())
    for doc, q, what in ((A10, Q_A, "step 2a of AMENDMENT 10"), (D51PRE, Q_B, "D51's sentence on a gain of a few times the error"),
                         (D51RES, Q_C, "D51's rule for the functions"), (D51RES, Q_D, "D51's largest distance of the functional from exact"),
                         (D50PRE, Q_E, "D50's sentence on charges together")):
        if not has(os.path.join(S6, doc), q):
            gates.append("G52.4 %s is not in %s" % (what, doc))
        if not has(a.prereg, q):
            gates.append("G52.4 the words \"%s\" are not in the pre-registration given" % (q[:60] + " ..."))
    if ("%+.4f eV" % p52.D51_ERR) not in Q_D:
        gates.append("G52.4 the number cs6_d52_prepare.py holds for D51's largest distance is not the one in D51's result")

    # ------------------------------------------------------------ reading
    RUNG = dict(p52.RUNGS)
    nsite = len(p52.SITES)

    def read_one(d):
        r = d40.read(d, "pc", 1, d33b)
        r["status"], r["B"], r["count"], r["nb"], r["xc"], r["pts"], r["enuc"] = "measured", None, None, None, None, None, None
        if not r.get("exists"):
            r["status"] = "NOT MEASURED"
            return r
        r["nb"] = nbf(r["text"])
        r["xc"], r["pts"] = grid(r["text"])
        r["enuc"] = enuc(r["text"])
        if r["E"] is None:
            r["status"] = "DID NOT CONVERGE" if "SCF NOT CONVERGED" in r["text"] else "NOT MEASURED"
            return r
        if r["enuc"] is None:
            r["status"], r["why"] = "NOT MEASURED", ["prints no single nuclear repulsion"]
            return r
        r["B"], r["count"] = -(r["E"] - r["enuc"]) * EV, count(r["text"])
        return r

    RUN = {}
    for r in runs:
        d = os.path.join(D, r["run"])
        x = read_one(d)
        x.update(r)
        RUN[r["run"]] = x
        if not x.get("exists"):
            continue
        e = echoed(x["text"])
        want = open(os.path.join(d, "job.inp")).read().split("\n") if os.path.isfile(os.path.join(d, "job.inp")) else None
        while want and not want[-1].strip():
            want.pop()
        if e is None:
            if x.get("normal"):
                gates.append("G52.5 %s/job.out ends normally and does not print the input it ran" % r["run"])
        elif want is not None and e != want:
            gates.append("G52.5 the input printed at the top of %s/job.out is not the committed %s/job.inp" % (r["run"], r["run"]))
        if x["status"] == "measured" and x["nb"] != nsite * p50.NBF[r["set"]]:
            gates.append("G52.6 %s reports %s basis functions, not %d" % (r["run"], x["nb"] if x["nb"] is not None else "no single count of", nsite * p50.NBF[r["set"]]))
        # a rung is what it says only if ORCA prints the IntAcc that was asked for
        if x["status"] == "measured" and r["rung"] in p52.INTACC:
            want_acc = p52.INTACC[r["rung"]]
            if x["xc"] is None or abs(x["xc"] - want_acc) > TOL_ASKED:
                x["status"], x["asked"] = "NOT AS ASKED", "IntAcc %.1f was asked for the grid of the functional; ORCA prints %s" % (want_acc, f(x["xc"], "%.3f", "none"))

    if gates:
        sys.stdout.write("\n".join("FAILED  " + g for g in gates) + "\n")
        stop("%d gate(s) failed. Nothing was written." % len(gates))

    sets, systems = list(p52.SETS), [k for k, _ in p52.SYSTEMS]
    sites_of = dict(p52.SYSTEMS)
    LAD = list(p52.LADDER)
    funcs = [k for k, _ in p52.FUNCS]
    cands = list(p52.CANDIDATES)
    TOP = p52.TOP

    def X(level, name, rung, k):
        return RUN.get(p52.hf_name(name, k) if level == "HF" else p52.fn_name(level, name, rung, k))

    def B(level, name, rung, k):
        x = X(level, name, rung, k)
        return x["B"] if x is not None and x["status"] == "measured" else None

    def CNT(level, name, rung, k):
        x = X(level, name, rung, k)
        return x["count"] if x is not None and x["status"] == "measured" else None

    # ---- R52.2 and R52.3: the two programs, and the functions
    TWO, FUN = {}, {}
    for name in sets:
        for k in systems:
            b = B("HF", name, None, k)
            d2 = None if b is None else b - ex.HF[(name, k)]
            TWO[(name, k)] = (d2, "NOT MEASURED" if b is None else ("AGREES" if abs(d2) <= p52.TOL_TWO else "DOES NOT AGREE"))
            d3 = None if b is None else ex.LIMIT[k] - b
            FUN[(name, k)] = (d3, "NOT MEASURED" if b is None else ("ADEQUATE" if abs(d3) <= p52.TOL_ADEQ else "NOT ADEQUATE"))

    # ---- R52.4: the grid of the working functional
    GR = {}
    for name in sets:
        for k in systems:
            top = next((rg for rg in ("i40", "i20") if B("CAM", name, rg, k) is not None), None)
            nxt = None
            if top is not None:
                below = LAD[:LAD.index(top)]
                nxt = below[-1] if below and B("CAM", name, below[-1], k) is not None else None
            g = dict(top=top, nxt=nxt, d=None, word="not decided", low=None)
            if top is not None and nxt is not None:
                g["d"] = B("CAM", name, top, k) - B("CAM", name, nxt, k)
                g["word"] = "CONVERGED" if abs(g["d"]) <= p52.TOL_GRID else "NOT CONVERGED"
                if g["word"] == "CONVERGED":
                    upto = LAD[:LAD.index(top) + 1]
                    for i, rg in enumerate(upto):
                        vals = [B("CAM", name, r2, k) for r2 in upto[i:]]
                        if all(v is not None and abs(v - B("CAM", name, top, k)) <= p52.TOL_GRID for v in vals):
                            g["low"] = rg
                            break
            GR[(name, k)] = g

    # ---- R52.5: set X2 against set X3 with the working functional
    CROSS = {}
    for k in systems:
        g2, g3 = GR[("X2", k)], GR[("X3", k)]
        if g2["word"] == "CONVERGED" and g3["word"] == "CONVERGED":
            dv = B("CAM", "X2", g2["top"], k) - B("CAM", "X3", g3["top"], k)
            CROSS[k] = (dv, "HOLDS" if abs(dv) <= p52.TOL_ADEQ else "DOES NOT HOLD")
        else:
            CROSS[k] = (None, "not decided")

    # ---- R52.6: RULE A
    SERVE = {}
    for name in sets:
        for k in systems:
            why, bad, open_ = [], False, False
            tw, fu, g = TWO[(name, k)][1], FUN[(name, k)][1], GR[(name, k)]
            if tw == "DOES NOT AGREE":
                bad = True
                why.append("ORCA's Hartree-Fock is %+.4f eV from the second program's" % TWO[(name, k)][0])
            elif tw == "NOT MEASURED":
                open_ = True
                why.append("the Hartree-Fock run is not measured")
            if fu == "NOT ADEQUATE":
                bad = True
                why.append("the functions: Hartree-Fock is %.4f eV from the limit" % abs(FUN[(name, k)][0]))
            if g["word"] == "NOT CONVERGED":
                bad = True
                why.append("the grid: %s minus %s is %+.4f eV" % (g["top"], g["nxt"], g["d"]))
            elif g["word"] != "CONVERGED":
                open_ = True
                why.append("the grid is not decided")
            if name == "X2":
                cw = CROSS[k][1]
                if cw == "DOES NOT HOLD":
                    bad = True
                    why.append("the functional holds %+.4f eV in X2 against X3" % CROSS[k][0])
                elif cw != "HOLDS":
                    open_ = True
                    why.append("X2 against X3 is not decided")
            SERVE[(name, k)] = ("DOES NOT SERVE" if bad else ("not decided" if open_ else "SERVES"), why)

    def verdict(name):
        ws = [SERVE[(name, k)][0] for k in systems]
        return "SERVES" if all(w == "SERVES" for w in ws) else ("DOES NOT SERVE" if any(w == "DOES NOT SERVE" for w in ws) else "NOT DECIDED")
    v2, v3 = verdict("X2"), verdict("X3")
    if v2 == "SERVES":
        chosen = "X2"
    elif v2 == "NOT DECIDED":
        chosen = "NOT DECIDED"
    else:
        chosen = "X3" if v3 == "SERVES" else ("NOT DECIDED" if v3 == "NOT DECIDED" else "NONE")
    grid_rung = None
    if chosen in sets:
        grid_rung = max((GR[(chosen, k)]["low"] for k in systems), key=LAD.index)

    # ---- R52.7: the candidates
    CG, CX = {}, {}
    for c in cands:
        for k in systems:
            lo, hi = B(c, "X2", "i12", k), B(c, "X2", TOP, k)
            CG[(c, k)] = (None, "not decided") if None in (lo, hi) else (hi - lo, "CONVERGED" if abs(hi - lo) <= p52.TOL_GRID else "NOT CONVERGED")
            b2, b3 = B(c, "X2", TOP, k), B(c, "X3", TOP, k)
            CX[(c, k)] = (None, "not decided") if None in (b2, b3) else (b2 - b3, "HOLDS" if abs(b2 - b3) <= p52.TOL_ADEQ else "DOES NOT HOLD")

    # ---- R52.8: every functional against the limit, in set X3 at the rung TOP
    ERR = dict(((fn, k), None if B(fn, "X3", TOP, k) is None else B(fn, "X3", TOP, k) - ex.LIMIT[k]) for fn in funcs for k in systems)

    def worst(fn):
        vs = [(ERR[(fn, k)], k) for k in systems]
        if any(v is None for v, _ in vs):
            return None
        return max(vs, key=lambda t: abs(t[0]))
    WORST = dict((fn, worst(fn)) for fn in funcs)

    # ---- R52.9: RULE B and RULE C
    complete = [c for c in cands if WORST[c] is not None]
    if WORST["CAM"] is None:
        second, second_why = "NOT DECIDED", "the working functional is not measured at all six systems in set X3 at %s" % TOP
    elif not complete:
        second, second_why = "NONE", "no candidate is measured at all six systems in set X3 at %s" % TOP
    else:
        best = min(complete, key=lambda c: (abs(WORST[c][0]), cands.index(c)))
        if abs(WORST[best][0]) < abs(WORST["CAM"][0]):
            second, second_why = best, "its largest distance from the limit is %.4f eV, the working functional's %.4f eV" % (abs(WORST[best][0]), abs(WORST["CAM"][0]))
        else:
            second, second_why = "NONE", "no candidate lies nearer the limit than the working functional (%.4f eV); the nearest is %s at %.4f eV" % (
                abs(WORST["CAM"][0]), NAMES[best], abs(WORST[best][0]))
    cam_far = None if WORST["CAM"] is None else abs(WORST["CAM"][0])
    limit_L = p52.FEW * max(p52.D51_ERR, cam_far if cam_far is not None else 0.0)

    # ---- R52.10: the grid given to hydrogen alone
    fs, fk = p52.FORM
    xg, xz, xi = X("CAM", fs, "g2", fk), X("CAM", fs, "z1", fk), X("CAM", fs, "i12", fk)
    if any(x is None or x["status"] != "measured" for x in (xg, xz, xi)):
        zword = "NOT DECIDED: one of the three runs is not measured"
    elif abs(xi["E"] - xg["E"]) <= p52.TOL_EH:
        zword = "NOT DECIDED: IntAcc 12 and the default grid give the same energy here"
    elif abs(xz["E"] - xi["E"]) <= p52.TOL_EH:
        zword = "THE SAME AS IntAcc 12 FOR EVERY ATOM: ORCA gives a ghost atom the grid asked for hydrogen"
    elif abs(xz["E"] - xg["E"]) <= p52.TOL_EH:
        zword = "THE SAME AS THE DEFAULT GRID: the form changed nothing"
    else:
        zword = "NEITHER: the form gives another energy than IntAcc 12 for every atom and than the default grid"

    # ------------------------------------------------------------ writing
    L = []
    bar = "=" * W

    def head(t):
        L.extend(["", bar, t, bar, ""])

    def para(t, indent=""):
        L.extend(textwrap.wrap(t, W, initial_indent=indent, subsequent_indent=indent))

    rc, hc = git("rev-parse", "HEAD")
    L += ["STAGE 6 - DIAGNOSTIC 52 - RESULT", "Charges together, one electron first", bar, "",
          "Read by cs6_d52.py under %s, which was committed before any job was submitted." % PREREG,
          "   results repository at %s; the pre-registration was added by %s" % (hc[:12] if rc == 0 else "?", pcm[:12])]
    for label, fn in (("reader", "cs6_d52.py"), ("inputs formed again", "cs6_d52_prepare.py"), ("what Claude expects", "cs6_d52_expected.py"),
                      ("D49's input by", "cs6_d50_prepare.py"), ("its lines from", "cs6_d48_prepare.py"), ("outputs read by", "cs6_d40.py"),
                      ("stability by", "cs6_d33b.py")):
        L.append("   %-20s %-22s sha256 %s" % (label, fn, sha(os.path.join(a.code, "phase2.2", fn))))
    L.append("")
    para("One electron, charges of +0.25 e at positive sites of B_n11_q0.25, no substrate and no negative charge. The functions of D50 "
         "(def2-TZVPD hydrogen with an s and a p function added at each of 0.035, 0.012 [X2] and also 0.004 [X3]) stand at all six "
         "positive sites in every run. B is minus the energy of the one electron, eV. THE LIMIT of a system was computed beforehand in "
         "a second program and is not a measurement. Tolerances: %.3f eV between the two programs, %.2f eV for the functions, %.3f eV "
         "for the grid. D52 adopts nothing: it names no design and quotes no ddE." % (p52.TOL_TWO, p52.TOL_ADEQ, p52.TOL_GRID))
    L += ["",
          "   system   charges of +0.25 e at                    the limit, eV"]
    for k in systems:
        L.append("   %-8s %-40s %.4f" % (k, ", ".join(str(s) for s in sites_of[k]), ex.LIMIT[k]))
    L += ["", "THE GATES, as they held:",
          "   G52.1  all %d files of d52/INPUT_CHECKSUMS.sha256 are the ones listed, and every file of d52 is" % nsum,
          "          the one formed again here from the committed files.",
          "   G52.2  the %d files of d52/SOURCES.sha256 are the ones listed." % nsrc,
          "   G52.3  no file under %s, and no committed table or document named here," % ", ".join(GUARDED),
          "          differs from the committed one; the pre-registration given is the committed one.",
          "   G52.4  the sentences quoted here are in the documents they are quoted from.",
          "   G52.5  every output that exists prints, as the input it ran, the committed input beside it.",
          "   G52.6  every measured run reports the number of basis functions its input asks for."]

    head("R52.1  THE RUNS")
    L.append("   what is run                                     measured   IntAcc ORCA prints")
    groups = [("HF", None, "Hartree-Fock")]
    groups += [("CAM", rg, "CAM-B3LYP, %s" % {"g2": "ORCA's default grid", "z1": "IntAcc 12 for hydrogen alone"}.get(rg, "IntAcc %.0f" % p52.INTACC.get(rg, 0))) for rg in ("g2", "z1", "i12", "i20", "i40")]
    groups += [(c, rg, "%s, IntAcc %.0f" % (NAMES[c], p52.INTACC[rg])) for c in cands for rg in ("i12", "i20")]
    for lv, rg, what in groups:
        xs = [x for x in RUN.values() if x["level"] == lv and (rg is None or x["rung"] == rg)]
        acc = [] if lv == "HF" else sorted(set("%.3f" % x["xc"] for x in xs if x["xc"] is not None))
        L.append("   %-47s %2d of %2d   %s" % (what, sum(1 for x in xs if x["status"] == "measured"), len(xs), ", ".join(acc) if acc else "-"))
    bad = [x for x in RUN.values() if x["status"] != "measured"]
    L.append("")
    if not bad:
        L.append("   Every one of the %d runs is measured." % len(RUN))
    else:
        L.append("   Measured: %d of %d. The others, with ORCA's words:" % (len(RUN) - len(bad), len(RUN)))
    for x in sorted(bad, key=lambda v: v["run"]):
        L.append("   %s: %s: %s" % (x["run"], x["status"], x.get("asked") or "; ".join(x.get("why") or ["-"])))
    nz = [abs(x["enuc"]) for x in RUN.values() if x["status"] == "measured" and x["enuc"] is not None]
    L += ["", "   The nuclear repulsion ORCA prints, over the measured runs: at most %s hartree." % ("%.1e" % max(nz) if nz else "not formed")]

    head("R52.2  ORCA AGAINST THE SECOND PROGRAM: HARTREE-FOCK IN THE SAME FUNCTIONS")
    para("Hartree-Fock is exact for one electron in the functions it is given and uses no grid. 'expected' is the lowest level in "
         "the same functions, computed beforehand in PySCF. AGREES: within %.3f eV." % p52.TOL_TWO)
    L += ["",
          "   set  system   B, ORCA    expected   ORCA minus expected   word"]
    for name in sets:
        for k in systems:
            L.append("   %-4s %-8s %-10s %+.4f    %-21s %s" % (name, k, f(B("HF", name, None, k)), ex.HF[(name, k)], f(TWO[(name, k)][0], "%+.5f"), TWO[(name, k)][1]))
    L += ["", "   The two programs agree at %d of %d (set, system)." % (sum(1 for v in TWO.values() if v[1] == "AGREES"), len(TWO))]

    head("R52.3  THE FUNCTIONS: HARTREE-FOCK AGAINST THE LIMIT")
    L += ["ADEQUATE: Hartree-Fock within %.2f eV of the limit." % p52.TOL_ADEQ, "",
          "   system   the limit   X2: B, HF  limit-HF   word           X3: B, HF  limit-HF   word"]
    for k in systems:
        L.append(("   %-8s %-11s " % (k, "%.4f" % ex.LIMIT[k]) + "".join("%-10s %-10s %-14s " % (
            f(B("HF", n, None, k)), f(FUN[(n, k)][0]), FUN[(n, k)][1]) for n in sets)).rstrip())

    head("R52.4  THE GRID OF THE WORKING FUNCTIONAL")
    para("CAM-B3LYP with auxiliary treatment A on IntAcc 12, 20 and 40. CONVERGED if the two highest rungs measured agree to %.3f eV. "
         "THE LOWEST RUNG THAT SERVES: the first whose B, and that of every rung above it, lies within %.3f eV of the top one. "
         "'expected' is CAM-B3LYP in the same functions on a fine grid with exact integrals, computed beforehand in PySCF; it is not "
         "a measurement. 'count' is the electrons the grid integrates (N(Total) in ORCA's output); it should be 1." % (p52.TOL_GRID, p52.TOL_GRID))
    L += ["", "   set  system   B: i12     i20        i40        top-next   the grid        serves from   expected"]
    for name in sets:
        for k in systems:
            g = GR[(name, k)]
            L.append("   %-4s %-8s %-10s %-10s %-10s %-10s %-15s %-13s %+.4f" % (
                name, k, f(B("CAM", name, "i12", k)), f(B("CAM", name, "i20", k)), f(B("CAM", name, "i40", k)),
                f(g["d"]), g["word"], g["low"] or "-", ex.CAM[(name, k)][0]))
    L += ["", "   set  system   count: i12   i20       i40"]
    for name in sets:
        for k in systems:
            L.append("   %-4s %-8s        %-8s %-9s %s" % (name, k, f(CNT("CAM", name, "i12", k), "%.4f"), f(CNT("CAM", name, "i20", k), "%.4f"),
                                                       f(CNT("CAM", name, "i40", k), "%.4f")))
    nconv = sum(1 for v in GR.values() if v["word"] == "CONVERGED")
    L += ["", "   The grid is CONVERGED at %d of %d (set, system), NOT CONVERGED at %d, not decided at %d." % (
        nconv, len(GR), sum(1 for v in GR.values() if v["word"] == "NOT CONVERGED"), sum(1 for v in GR.values() if v["word"] == "not decided"))]

    head("R52.5  THE FUNCTIONS FOR THE FUNCTIONAL: SET X2 AGAINST SET X3")
    para("Where the grid is CONVERGED in both sets: what CAM-B3LYP holds in X2 minus what it holds in X3, each on its top rung. "
         "HOLDS if they agree to %.2f eV." % p52.TOL_ADEQ)
    L += ["", "   system   B in X2    B in X3    X2 minus X3   word"]
    for k in systems:
        g2, g3 = GR[("X2", k)], GR[("X3", k)]
        L.append("   %-8s %-10s %-10s %-13s %s" % (k, f(None if g2["top"] is None else B("CAM", "X2", g2["top"], k)),
                                                f(None if g3["top"] is None else B("CAM", "X3", g3["top"], k)), f(CROSS[k][0]), CROSS[k][1]))

    head("R52.6  RULE A: THE FUNCTIONS AND THE GRID FOR THE BARE CHARGES OF THE DIANION")
    para("A set SERVES at a system if the two programs agree there (R52.2), the functions are ADEQUATE (R52.3) and the grid is "
         "CONVERGED (R52.4); set X2 must also hold what set X3 holds (R52.5). The rule takes set X2 if it SERVES at all six systems, "
         "and set X3 if X2 does not serve somewhere and X3 serves at all six.")
    L.append("")
    for name in sets:
        for k in systems:
            w, why = SERVE[(name, k)]
            L.append("   %-4s %-8s %-15s %s" % (name, k, w, "; ".join(why)))
    L.append("")
    if chosen in sets:
        L.append("   RULE A GIVES: SET %s, with auxiliary treatment A, and IntAcc %.0f for the ghost atoms (rung %s:" % (chosen, p52.INTACC[grid_rung], grid_rung))
        L.append("   the lowest that serves at every system).")
    elif chosen == "NONE":
        L.append("   RULE A GIVES: NONE. Neither set serves at all six systems.")
    else:
        L.append("   RULE A GIVES: NOT DECIDED. At a system the rule could not be applied.")

    head("R52.7  THE CANDIDATES FOR A SECOND LEVEL: THEIR GRID AND THEIR FUNCTIONS")
    para("Each candidate in set X2 on IntAcc 12 and 20, and in set X3 on IntAcc 20. The grid: IntAcc 20 minus IntAcc 12 in X2, "
         "CONVERGED within %.3f eV. The functions: X2 minus X3 at IntAcc 20, HOLDS within %.2f eV. The electron counts are in "
         "d52_systems.tsv." % (p52.TOL_GRID, p52.TOL_ADEQ))
    L.append("")
    for c in cands:
        xs = [x for x in RUN.values() if x["level"] == c]
        L.append("   %s: measured %d of %d runs." % (NAMES[c], sum(1 for x in xs if x["status"] == "measured"), len(xs)))
        L.append("      system   B, X2 i12  X2 i20     X3 i20     i20-i12    the grid        X2-X3      the functions")
        for k in systems:
            L.append("      %-8s %-10s %-10s %-10s %-10s %-15s %-10s %s" % (
                k, f(B(c, "X2", "i12", k)), f(B(c, "X2", TOP, k)), f(B(c, "X3", TOP, k)), f(CG[(c, k)][0]), CG[(c, k)][1], f(CX[(c, k)][0]), CX[(c, k)][1]))
        L.append("")
    L.pop()

    head("R52.8  EVERY FUNCTIONAL AGAINST THE LIMIT")
    para("B in set X3 at IntAcc 20, minus the limit, eV. Above zero: the functional holds the electron at the charges more strongly "
         "than it is held. 'expected' is the same difference for CAM-B3LYP on a fine grid with exact integrals, computed beforehand "
         "in PySCF. D51 found for ONE charge alone, from 0.25 e up, at most %+.4f eV." % p52.D51_ERR)
    L += ["", "   system   the limit   " + "".join("%-10s" % fn for fn in funcs) + "expected (CAM)   CAM at its top rung"]
    for k in systems:
        g3 = GR[("X3", k)]
        L.append("   %-8s %-11s " % (k, "%.4f" % ex.LIMIT[k]) + "".join("%-10s" % f(ERR[(fn, k)]) for fn in funcs) + "%+.4f          %s" % (
            ex.CAM[("X3", k)][0] - ex.LIMIT[k], f(None if g3["top"] is None else B("CAM", "X3", g3["top"], k) - ex.LIMIT[k])))
    L.append("")
    for fn in funcs:
        L.append("   %-24s %s" % (NAMES[fn] + ":", "not formed: not measured at all six systems" if WORST[fn] is None else
                                  "the largest distance from the limit is %+.4f eV (%s)" % (WORST[fn][0], WORST[fn][1])))

    head("R52.9  RULE B AND RULE C")
    para("RULE B. The second level is the candidate whose largest distance from the limit (R52.8) is smallest, if it is measured "
         "at all six systems and lies nearer the limit than the working functional.")
    L.append("")
    if second in cands:
        L.append("   RULE B GIVES: %s. %s." % (NAMES[second], second_why[0].upper() + second_why[1:]))
        L.append("   Its grid at IntAcc 12 is CONVERGED at %d of %d systems, and it holds in X2 what it holds in X3 at %d." % (
            sum(1 for k in systems if CG[(second, k)][1] == "CONVERGED"), len(systems), sum(1 for k in systems if CX[(second, k)][1] == "HOLDS")))
    else:
        L.append("   RULE B GIVES: %s. %s." % (second, second_why[0].upper() + second_why[1:]))
    L.append("")
    para("RULE C. D51's pre-registration: \"%s\" The limit L is %d times the larger of the working functional's largest distance from "
         "the limit over the six systems and D51's %.4f eV for one charge." % (Q_B, p52.FEW, p52.D51_ERR))
    L += ["",
          "   The working functional's largest distance over the six systems: %s." % ("not formed" if cam_far is None else "%.4f eV" % cam_far),
          "   RULE C GIVES: L = %.4f eV.%s" % (limit_L, "" if cam_far is not None else " It rests on D51's number alone.")]

    head("R52.10  THE GRID GIVEN TO HYDROGEN ALONE")
    para("At (%s, %s): CAM-B3LYP on ORCA's default grid (g2), with IntAcc 12 for every atom (i12), and with IntAcc 12 asked for "
         "hydrogen alone by SpecialGridAtoms 1 (z1). In these runs every atom is a ghost hydrogen, so z1 should be i12 if ORCA counts "
         "a ghost atom as hydrogen. The same: to %.0e hartree." % (fs, fk, p52.TOL_EH))
    L += ["",
          "   rung   status         E, hartree        count     grid points   IntAcc ORCA prints"]
    for rg, x in (("g2", xg), ("i12", xi), ("z1", xz)):
        L.append("   %-6s %-14s %-17s %-9s %-13s %s" % (rg, "-" if x is None else x["status"], f(None if x is None or x["status"] != "measured" else x["E"], "%.10f"),
                                                        f(None if x is None else x["count"], "%.6f"), "-" if x is None or x["pts"] is None else str(x["pts"]),
                                                        f(None if x is None else x["xc"], "%.3f")))
    L += ["", "   THE FORM: %s." % zword]

    head("WHAT THIS CANNOT SAY (section 6 of the pre-registration, in short; its own words govern)")
    L += ["   One electron in the field of bare charges says what the functions, the grid and a functional do",
          "   for one electron. It does not say what a run of the dianion does: there the electron that may",
          "   leave is one of 118, and the negative charges and the substrate are present.",
          "   The limit and the expected numbers were computed by Claude in a second program. The limit is the",
          "   reference of the rules; that it is converged rests on its repetition in a larger basis.",
          "   A candidate that ORCA did not run for one electron is not thereby unfit for the dianion, and one",
          "   that lies near the limit here is not thereby right there.",
          "   D50's and D51's results stand as committed. D52 adopts nothing."]
    summary = ("Runs measured: %d of %d. The two programs agree at %d of %d. RULE A: %s. RULE B: %s. RULE C: L = %.4f eV. "
               "The working functional's largest distance from the limit over the six systems: %s. The grid for hydrogen alone: %s." % (
                   len(RUN) - len(bad), len(RUN), sum(1 for v in TWO.values() if v[1] == "AGREES"), len(TWO),
                   ("set %s with treatment A at IntAcc %.0f" % (chosen, p52.INTACC[grid_rung])) if chosen in sets else chosen,
                   NAMES[second] if second in cands else second, limit_L, "not formed" if WORST["CAM"] is None else "%+.4f eV (%s)" % WORST["CAM"],
                   zword.split(":")[0]))
    L += ["", summary]

    T1 = ["run\tlevel\tset\trung\tsystem\tstatus\tE_Eh\tnuclear_repulsion_Eh\tB_eV\telectron_count\tbasis_functions\tIntAcc_printed\tgrid_points\twhy"]
    for r in runs:
        x = RUN[r["run"]]
        ok = x["status"] == "measured"
        T1.append("\t".join([r["run"], r["level"], r["set"], r["rung"], r["system"], x["status"], f(x["E"] if ok else None, "%.10f", ""),
                             f(x["enuc"], "%.10f", ""), f(x["B"] if ok else None, "%.6f", ""), f(x["count"] if ok else None, "%.6f", ""),
                             "" if x["nb"] is None else str(x["nb"]), f(x["xc"], "%.3f", ""), "" if x["pts"] is None else str(x["pts"]),
                             "" if ok else (x.get("asked") or "; ".join(x.get("why") or []))]))
    T2 = ["level\tset\tsystem\tlimit_eV\tB_HF_eV\tB_i12_eV\tB_i20_eV\tB_i40_eV\tcount_i12\tcount_i20\tcount_i40\tB_X3_i20_minus_limit_eV\texpected_HF_eV\texpected_CAM_eV"]
    for name in sets:
        for k in systems:
            T2.append("\t".join(["HF", name, k, "%.6f" % ex.LIMIT[k], f(B("HF", name, None, k), "%.6f", ""), "", "", "", "", "", "", "", "%.6f" % ex.HF[(name, k)], ""]))
    for fn in funcs:
        for name in sets:
            for k in systems:
                T2.append("\t".join([fn, name, k, "%.6f" % ex.LIMIT[k], ""] + [f(B(fn, name, rg, k), "%.6f", "") for rg in LAD] +
                                    [f(CNT(fn, name, rg, k), "%.6f", "") for rg in LAD] + [f(ERR[(fn, k)], "%.6f", "") if name == "X3" else "", "",
                                     ("%.6f" % ex.CAM[(name, k)][0]) if fn == "CAM" else ""]))
    T3 = ["what\tvalue"]
    T3 += ["rule_A_set\t%s" % chosen, "rule_A_rung\t%s" % (grid_rung or ""), "rule_A_IntAcc\t%s" % (("%.1f" % p52.INTACC[grid_rung]) if grid_rung else ""),
           "rule_B_second_level\t%s" % (NAMES[second] if second in cands else second), "rule_C_L_eV\t%.6f" % limit_L,
           "working_functional_largest_distance_eV\t%s" % f(None if WORST["CAM"] is None else WORST["CAM"][0], "%.6f", ""),
           "grid_for_hydrogen_alone\t%s" % zword]
    for fn in funcs:
        T3.append("largest_distance_from_the_limit_%s_eV\t%s" % (fn, f(None if WORST[fn] is None else WORST[fn][0], "%.6f", "")))
    for name in sets:
        for k in systems:
            T3.append("serves_%s_%s\t%s" % (name, k, SERVE[(name, k)][0]))
    L = [l.rstrip() for l in L]
    for path, lines in ((a.result, L), (a.tsv, T1), (tsv2, T2), (tsv3, T3)):
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")
    sys.stdout.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
