#!/usr/bin/env python3
"""cs6_d51.py - the reader of Stage 6 DIAGNOSTIC 51: the calibration of D50 done again, with the integration grid
measured. Committed before any job of D51 is submitted.

One electron at a bare charge q with the functions of D50 (sets X2 and X3). B is minus the energy of the one
electron, in eV; exact, in a complete basis, is 13.605693122994 q^2. Each run is read with the committed reader of
outputs (cs6_d40.py), as D48's, D49's and D50's readers read theirs. Three things that D50's rule held together are
read apart:

  THE FUNCTIONS            Hartree-Fock, which is exact for one electron and uses no grid, within 0.05 eV of exact.
  THE GRID                 the functional's value on the two highest rungs of the ladder that are measured (IntAcc 40
                           and 20, or 20 and 12) agree to 0.005 eV: CONVERGED. The electron count ORCA prints is set
                           beside it. THE LOWEST RUNG THAT SERVES is the first whose value, and that of every rung
                           above it, lies within 0.005 eV of the top one.
  THE FUNCTIONS FOR THE FUNCTIONAL   where the grid is converged in both sets, the functional's value in set X2
                           within 0.05 eV of its value in set X3.

THE RULE: set X2 is ADEQUATE at a size if all three hold there, and is the set for the capture test if it is ADEQUATE
at every size from 0.25 e up. What the functional itself lies from exact is reported and is not part of the rule.
The rule is applied only if the Hartree-Fock runs and rung g2, which are D50's inputs byte for byte, give D50's
committed energies again.

The inputs are not taken on trust: every file of d51 is formed again by cs6_d51_prepare.py from the committed files
and compared, every output must print the committed input beside it, and the committed files this leans on must be
the committed ones. A gate that fails stops the reader before anything is written. D51 adopts nothing.

usage: python3 cs6_d51.py --res <results repository> --code <code repository> --prereg <the pre-registration>
                          --result <the result to write> --tsv <the table of runs to write>
"""
import argparse
import hashlib
import importlib.util
import os
import re
import subprocess
import sys

S = "phase2.2/cs_stage6"
PREREG = "STAGE6_DIAGNOSTIC51_20261123.txt"
D50PRE = "STAGE6_DIAGNOSTIC50_20261122.txt"
D50RES = "STAGE6_DIAGNOSTIC50_RESULT_20261123.txt"
EV = 27.211386245988
RYD = EV / 2.0
TOL_ASKED = 0.0015     # an IntAcc ORCA prints against the one asked for
# quoted sentences: each must be in the document it is quoted from and in the pre-registration given
Q_A = "THE FUNCTIONS FOR THE CAPTURE TEST: NONE. No set has a treatment that is adequate at every size from 0.25 e up."     # D50's result
Q_B = "THE RULE GIVES NONE: there are no functions for the test."                                                             # D50, section 5
Q_C = "A run that DID NOT TAKE with these functions is a test of capture at 0.25 e only if R50.3"                            # D50, D50.3
GUARDED = ["d48", "d49", "d50"]
DOCS = ["d50.tsv", "d50_calibration.tsv", "d50_dianion.tsv", D50PRE, D50RES]
W = 96


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


def grids(text):
    """(IntAcc and number of points of the grid of the functional, the IntAcc of each COSX grid), as ORCA prints them"""
    xc, pts, cosx = None, None, []
    k = text.find("DFT GRID GENERATION")
    kx = text.find("COSX GRID GENERATION")
    if k >= 0:
        part = text[k:kx] if kx > k else text[k:]
        m = re.search(r"General Integration Accuracy\s+IntAcc\s*\.+\s*(\d+\.\d+)", part)
        p = re.search(r"Total number of grid points\s*\.+\s*(\d+)", part)
        xc, pts = (float(m.group(1)) if m else None), (int(p.group(1)) if p else None)
    if kx >= 0:
        for i in (1, 2, 3):                                  # the three grids ORCA calls GRIDX 1, 2 and 3
            m = re.search(r"GRIDX %d\b.*?General Integration Accuracy\s+IntAcc\s*\.+\s*(\d+\.\d+)" % i, text[kx:], re.S)
            if m:
                cosx.append(float(m.group(1)))
    return xc, pts, cosx


def f(v, fmt="%+.4f", none="-"):
    return none if v is None else fmt % v


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    a.res, a.code, a.prereg = (os.path.abspath(x) for x in (a.res, a.code, a.prereg))
    S6 = os.path.join(a.res, S)
    D = os.path.join(S6, "d51")
    paths = dict((k, os.path.join(a.code, "phase2.2", fn)) for k, fn in (
        ("d40", "cs6_d40.py"), ("d33b", "cs6_d33b.py"), ("p48", "cs6_d48_prepare.py"), ("p50", "cs6_d50_prepare.py"), ("p51", "cs6_d51_prepare.py"),
        ("exp", "cs6_d51_expected.py")))
    for k, p in paths.items():
        if not os.path.isfile(p):
            stop("the committed script %s is missing." % p)
    for fn in [a.prereg, os.path.join(D, "RUNS.tsv")] + [os.path.join(S6, x) for x in DOCS]:
        if not os.path.isfile(fn):
            stop("%s is missing." % fn)
    tsv2, tsv3 = os.path.join(os.path.dirname(a.tsv), "d51_ladder.tsv"), os.path.join(os.path.dirname(a.tsv), "d51_rule.tsv")
    for fn in (a.result, a.tsv, tsv2, tsv3):
        if os.path.exists(fn):
            stop("%s exists already; this reader overwrites nothing." % fn)
    d40, d33b = load(paths["d40"], "d51_cs6_d40"), load(paths["d33b"], "d51_cs6_d33b")
    p51 = load(paths["p51"], "d51_cs6_d51_prepare")
    gates = []

    def git(*args):
        try:
            q = subprocess.run(["git", "-C", a.res] + list(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            return q.returncode, q.stdout.strip()
        except OSError:
            return 127, ""

    # ------------------------------------------------------------ G51.3: the committed files this leans on
    guard = ["%s/%s" % (S, g) for g in GUARDED] + ["%s/%s" % (S, x) for x in DOCS]
    rc, st = git("status", "--porcelain", "--untracked-files=all", "--", *guard)
    if rc != 0:
        gates.append("G51.3 could not be checked: git did not run")
    elif st:
        gates.append("G51.3 a file under %s, or a committed table or document this reader names, is not the committed one" % ", ".join(GUARDED))
    rc, pcm = git("log", "--diff-filter=A", "--format=%H", "-1", "--", "%s/%s" % (S, PREREG))
    if rc != 0 or not pcm:
        gates.append("G51.3 the commit that added %s could not be found" % PREREG)
    else:
        rc1, b1 = git("rev-parse", "%s:%s/%s" % (pcm, S, PREREG))
        rc2, b2 = git("hash-object", a.prereg)
        if rc1 != 0 or rc2 != 0 or b1 != b2:
            gates.append("G51.3 the pre-registration given is not the committed %s" % PREREG)

    # ------------------------------------------------------------ G51.1 and G51.2: formed again from the committed files
    try:
        F, runs, p50, p48, ex = p51.build(S6, a.code)
    except Exception as e:
        sys.stdout.write("\n".join("FAILED  " + g for g in gates + ["G51.1 the files of d51 could not be formed again from the committed files: %s: %s" % (type(e).__name__, e)]) + "\n")
        stop("%d gate(s) failed. Nothing was written." % (len(gates) + 1))
    psrc = os.path.join(D, "SOURCES.sha256")
    if not os.path.isfile(psrc) or open(psrc).read() != F["SOURCES.sha256"]:
        gates.append("G51.2 d51/SOURCES.sha256 is not the list formed again from the committed files")
    nsrc = 0
    for l in F["SOURCES.sha256"].splitlines():
        h, fn = l.split(None, 1)
        nsrc += 1
        p = os.path.join(S6, fn)
        if not os.path.isfile(p) or sha(p) != h:
            gates.append("G51.2 %s is not the file d51/SOURCES.sha256 lists" % fn)
    nsum = 0
    pic = os.path.join(D, "INPUT_CHECKSUMS.sha256")
    if not os.path.isfile(pic):
        gates.append("G51.1 d51/INPUT_CHECKSUMS.sha256 is missing")
    else:
        if open(pic).read() != F["INPUT_CHECKSUMS.sha256"]:
            gates.append("G51.1 d51/INPUT_CHECKSUMS.sha256 is not the list formed again here")
        for l in open(pic).read().splitlines():
            if l.strip():
                h, fn = l.split(None, 1)
                nsum += 1
                p = os.path.join(D, fn)
                if not os.path.isfile(p) or sha(p) != h:
                    gates.append("G51.1 %s is not the file listed in INPUT_CHECKSUMS.sha256" % fn)
    for k in sorted(F):
        p = os.path.join(D, k)
        if not os.path.isfile(p) or open(p).read() != F[k]:
            gates.append("G51.1 d51/%s is not the file formed again from the committed files" % k)

    # ------------------------------------------------------------ G51.4: the quoted sentences
    def has(path, sent):
        return " ".join(sent.split()) in " ".join(open(path, errors="replace").read().split())
    for doc, q, what in ((D50RES, Q_A, "D50's word NONE"), (D50PRE, Q_B, "D50's sentence on the rule giving NONE"), (D50PRE, Q_C, "D50.3")):
        if not has(os.path.join(S6, doc), q):
            gates.append("G51.4 %s is not in %s" % (what, doc))
        if not has(a.prereg, q):
            gates.append("G51.4 the words \"%s\" are not in the pre-registration given" % (q[:60] + " ..."))

    # ------------------------------------------------------------ reading
    RUNG = dict((r[0], r[1:]) for r in p51.RUNGS)

    def read_one(d):
        r = d40.read(d, "pc", 1, d33b)
        r["status"], r["B"], r["count"], r["nb"], r["xc"], r["pts"], r["cosx"] = "measured", None, None, None, None, None, []
        if not r.get("exists"):
            r["status"] = "NOT MEASURED"
            return r
        r["nb"] = nbf(r["text"])
        r["xc"], r["pts"], r["cosx"] = grids(r["text"])
        if r["E"] is None:
            r["status"] = "DID NOT CONVERGE" if "SCF NOT CONVERGED" in r["text"] else "NOT MEASURED"
            return r
        r["B"], r["count"] = -r["E"] * EV, count(r["text"])
        return r

    # ---- G51.7: D50's committed outputs
    OLD = {}
    for r in runs:
        if r["d50"]:
            OLD[r["run"]] = read_one(os.path.join(S6, "d50", r["d50"]))
            if OLD[r["run"]]["status"] != "measured":
                gates.append("G51.7 D50's committed output d50/%s is not readable" % r["d50"])

    # ---- every run of D51, and G51.5, G51.6
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
                gates.append("G51.5 %s/job.out ends normally and does not print the input it ran" % r["run"])
        elif want is not None and e != want:
            gates.append("G51.5 the input printed at the top of %s/job.out is not the committed %s/job.inp" % (r["run"], r["run"]))
        if x["status"] == "measured" and x["nb"] != p50.NBF[r["set"]]:
            gates.append("G51.6 %s reports %s basis functions, not %d" % (r["run"], x["nb"] if x["nb"] is not None else "no single count of", p50.NBF[r["set"]]))
        # a rung is what it says only if ORCA prints the IntAcc that was asked for
        if x["status"] == "measured" and r["rung"]:
            kind, word, xc, cosx = RUNG[r["rung"]]
            if xc is not None and (x["xc"] is None or abs(x["xc"] - xc) > TOL_ASKED):
                x["status"], x["asked"] = "NOT AS ASKED", "IntAcc %.1f was asked for the grid of the functional; ORCA prints %s" % (xc, f(x["xc"], "%.3f", "none"))
            elif cosx is not None and (len(x["cosx"]) != 3 or any(abs(v - cosx) > TOL_ASKED for v in x["cosx"])):
                x["status"], x["asked"] = "NOT AS ASKED", "IntAcc %.1f was asked for the three COSX grids; ORCA prints %s" % (
                    cosx, ", ".join("%.3f" % v for v in x["cosx"]) or "none")

    if gates:
        sys.stdout.write("\n".join("FAILED  " + g for g in gates) + "\n")
        stop("%d gate(s) failed. Nothing was written." % len(gates))

    sets, sizes, LAD = sorted(p50.ADD), list(p50.SIZES), list(p51.LADDER)
    rungs = [r[0] for r in p51.RUNGS]

    def B(name, rung, k):
        x = RUN[p51.cam_name(name, rung, k)]
        return x["B"] if x["status"] == "measured" else None

    def CNT(name, rung, k):
        x = RUN[p51.cam_name(name, rung, k)]
        return x["count"] if x["status"] == "measured" else None

    def HF(name, k):
        x = RUN[p51.hf_name(name, k)]
        return x["B"] if x["status"] == "measured" else None

    # ---- R51.2: the cluster gives D50's numbers again
    REP = []
    for r in runs:
        if r["d50"]:
            x, o = RUN[r["run"]], OLD[r["run"]]
            de = None if x["status"] != "measured" else x["E"] - o["E"]
            dc = None if (x["status"] != "measured" or x["count"] is None or o["count"] is None) else x["count"] - o["count"]
            REP.append((r, de, dc))
    nrep = sum(1 for r, de, dc in REP if de is not None and abs(de) <= p51.TOL_EH)
    reproduced = nrep == len(REP)

    # ---- R51.4: the grid at each size
    GR = {}
    for name in sets:
        for k in sizes:
            top = next((rg for rg in ("i40", "i20") if B(name, rg, k) is not None), None)
            nxt = None
            if top is not None:
                below = LAD[:LAD.index(top)]
                nxt = below[-1] if below and B(name, below[-1], k) is not None else None
            g = dict(top=top, nxt=nxt, d=None, word="not decided", low=None)
            if top is not None and nxt is not None:
                g["d"] = B(name, top, k) - B(name, nxt, k)
                g["word"] = "CONVERGED" if abs(g["d"]) <= p51.TOL_GRID else "NOT CONVERGED"
                if g["word"] == "CONVERGED":
                    upto = LAD[:LAD.index(top) + 1]
                    for i, rg in enumerate(upto):
                        vals = [B(name, r2, k) for r2 in upto[i:]]
                        if all(v is not None and abs(v - B(name, top, k)) <= p51.TOL_GRID for v in vals):
                            g["low"] = rg
                            break
            GR[(name, k)] = g

    # ---- R51.5 to R51.7
    RULE = {}
    for k in sizes:
        q = k / 100.0
        exact = RYD * q * q
        z = dict(exact=exact, hf={}, hfw={}, fd=None, fw="not decided", word="not decided", why=[])
        for name in sets:
            z["hf"][name] = HF(name, k)
            z["hfw"][name] = "NOT MEASURED" if z["hf"][name] is None else ("ADEQUATE" if abs(exact - z["hf"][name]) <= p51.TOL_ADEQ else "NOT ADEQUATE")
        g2, g3 = GR[("X2", k)], GR[("X3", k)]
        if g2["word"] == "CONVERGED" and g3["word"] == "CONVERGED":
            z["fd"] = B("X2", g2["top"], k) - B("X3", g3["top"], k)
            z["fw"] = "HOLDS" if abs(z["fd"]) <= p51.TOL_ADEQ else "DOES NOT HOLD"
        if not reproduced:
            z["why"].append("D50's energies were not given again (R51.2)")
        elif z["hfw"]["X2"] == "NOT MEASURED":
            z["why"].append("the Hartree-Fock run is not measured")
        elif z["hfw"]["X2"] == "NOT ADEQUATE":
            z["word"] = "NOT ADEQUATE"
            z["why"].append("the functions: Hartree-Fock is %.4f eV from exact" % abs(exact - z["hf"]["X2"]))
        elif z["fw"] == "not decided":
            z["why"].append("the grid is %s in X2 and %s in X3" % (g2["word"], g3["word"]))
        elif z["fw"] == "DOES NOT HOLD":
            z["word"] = "NOT ADEQUATE"
            z["why"].append("the functional holds %+.4f eV in X2 against X3" % z["fd"])
        else:
            z["word"] = "ADEQUATE"
        RULE[k] = z
    up = [k for k in sizes if k >= p51.FROM]
    if any(RULE[k]["word"] == "NOT ADEQUATE" for k in up):
        chosen = "NONE"
    elif any(RULE[k]["word"] != "ADEQUATE" for k in up):
        chosen = "NOT DECIDED"
    else:
        chosen = "X2"
    lows = [GR[("X2", k)]["low"] for k in up]
    grid_rung = None if any(v is None for v in lows) else max(lows, key=LAD.index)

    # ---- what was expected, beside what was measured
    worst_c, worst_b = None, None
    for name in sets:
        for k in sizes:
            for rg in LAD:
                c = CNT(name, rg, k)
                if c is not None:
                    dv = min((c - ex.MODEL[(n0, name, rg, k)][0] for n0 in ex.MODELS), key=abs)
                    if worst_c is None or abs(dv) > abs(worst_c[0]):
                        worst_c = (dv, name, k, rg)
            g = GR[(name, k)]
            if g["word"] == "CONVERGED":
                dv = B(name, g["top"], k) - ex.FULL[(name, k)][1]
                if worst_b is None or abs(dv) > abs(worst_b[0]):
                    worst_b = (dv, name, k)

    # ------------------------------------------------------------ writing
    L = []
    bar = "=" * W

    def head(t):
        L.extend(["", bar, t, bar, ""])

    rc, hc = git("rev-parse", "HEAD")
    L += ["STAGE 6 - DIAGNOSTIC 51 - RESULT", "The calibration of D50 done again, with the integration grid measured", bar, "",
          "Read by cs6_d51.py under %s, which was committed before any job was submitted." % PREREG,
          "   results repository at %s; the pre-registration was added by %s" % (hc[:12] if rc == 0 else "?", pcm[:12])]
    for label, fn in (("reader", "cs6_d51.py"), ("inputs formed again", "cs6_d51_prepare.py"), ("what Claude expects", "cs6_d51_expected.py"),
                      ("D50's inputs by", "cs6_d50_prepare.py"), ("its lines from", "cs6_d48_prepare.py"), ("outputs read by", "cs6_d40.py"),
                      ("stability by", "cs6_d33b.py")):
        L.append("   %-20s %-22s sha256 %s" % (label, fn, sha(os.path.join(a.code, "phase2.2", fn))))
    L += ["", "One electron, a bare charge q, the functions of D50: def2-TZVPD hydrogen with an s and a p function",
          "added at each of 0.035, 0.012 (set X2) and also 0.004 (set X3). B is minus the energy of the one",
          "electron, eV. Tolerances: %.2f eV for the functions, %.3f eV for the grid, %.0e hartree for giving" % (p51.TOL_ADEQ, p51.TOL_GRID, p51.TOL_EH),
          "D50's energies again. D51 adopts nothing: it opens and closes no size, no distance and no design.", "",
          "THE GATES, as they held:",
          "   G51.1  all %d files of d51/INPUT_CHECKSUMS.sha256 are the ones listed, and every file of d51 is" % nsum,
          "          the one formed again here from the committed files.",
          "   G51.2  the %d files of d51/SOURCES.sha256 are the ones listed." % nsrc,
          "   G51.3  no file under %s, and no committed table or document named here, differs" % ", ".join(GUARDED),
          "          from the committed one; the pre-registration given is the committed one.",
          "   G51.4  the sentences quoted here are in the documents they are quoted from.",
          "   G51.5  every output that exists prints, as the input it ran, the committed input beside it.",
          "   G51.6  every measured run reports the number of basis functions its input asks for.",
          "   G51.7  D50's %d committed outputs set beside the runs are readable." % len(OLD)]

    head("R51.1  THE RUNS")
    L.append("   what is run                               measured   IntAcc ORCA prints   grid points (X2, 0.25 e)")
    rows = [("HF", "Hartree-Fock, D50's input", [RUN[p51.hf_name(n, k)] for n in sets for k in sizes], None)]
    for rg in rungs:
        kind, word, xc, cosx = RUNG[rg]
        what = {"g2": "D50's input (default grid)", "g3": "DefGrid3", "c40": "IntAcc 40, COSX grids 40 too", "n2": "NORI, default grid", "n40": "NORI, IntAcc 40"}.get(rg, "IntAcc %.0f" % (xc or 0))
        rows.append((rg, what, [RUN[p51.cam_name(n, rg, k)] for n in sets for k in sizes], RUN[p51.cam_name("X2", rg, 25)]))
    for rg, what, xs, x25 in rows:
        acc = [] if rg == "HF" else sorted(set("%.3f" % x["xc"] for x in xs if x["xc"] is not None))
        L.append("   %-4s %-36s %2d of %2d   %-19s  %s" % (rg, what, sum(1 for x in xs if x["status"] == "measured"), len(xs), ", ".join(acc) if acc else "-",
                                                         "-" if x25 is None or x25["pts"] is None else str(x25["pts"])))
    bad = [x for x in RUN.values() if x["status"] != "measured"]
    L.append("")
    if not bad:
        L.append("   Every one of the %d runs is measured." % len(RUN))
    for x in sorted(bad, key=lambda v: v["run"]):
        L.append("   %s: %s: %s" % (x["run"], x["status"], x.get("asked") or "; ".join(x.get("why") or ["-"])))

    head("R51.2  DOES THE CLUSTER GIVE D50'S NUMBERS AGAIN?")
    L += ["The Hartree-Fock runs and rung g2 are D50's committed inputs, byte for byte. Each energy is set",
          "against the committed energy of the same input in d50, and each electron count against D50's.", ""]
    des = [abs(de) for r, de, dc in REP if de is not None]
    dcs = [abs(dc) for r, de, dc in REP if dc is not None]
    L.append("   energies the same to %.0e hartree: %d of %d; the largest difference is %s hartree." % (p51.TOL_EH, nrep, len(REP), "%.2e" % max(des) if des else "not formed"))
    L.append("   electron counts of rung g2 against D50's: formed at %d of %d; the largest difference is %s." % (len(dcs), sum(1 for r, de, dc in REP if r["rung"]), "%.2e" % max(dcs) if dcs else "not formed"))
    for r, de, dc in REP:
        if de is None or abs(de) > p51.TOL_EH:
            L.append("   NOT THE SAME: %s: %s" % (r["run"], "not measured" if de is None else "%+.3e hartree from d50/%s" % (de, r["d50"])))
    L += ["", "   " + ("D50'S ENERGIES ARE GIVEN AGAIN: the rule is applied." if reproduced else "D50'S ENERGIES ARE NOT GIVEN AGAIN: the ladder is printed and the rule is NOT applied.")]

    head("R51.3  THE LADDER: WHAT THE GRID OF THE FUNCTIONAL COUNTS AND WHAT THE FUNCTIONAL HOLDS")
    L += ["Rungs g2 to i40: the same run with more radial points in the grid of the functional. 'count' is",
          "the electrons that grid integrates (N(Total) in ORCA's output); it should be 1. 'model' is what",
          "Claude's two models of ORCA's radial grid gave beforehand (cs6_d51_expected.py). 'full' is what",
          "the functional holds in that set on a grid that counts the whole electron, computed by Claude in",
          "a second program. None of these expectations is a measurement.", ""]
    for name in sets:
        L.append("SET %s, the electron count" % name)
        L.append("   q, e   " + "".join("%-9s" % rg for rg in LAD) + "  model at g2      model at i40")
        for k in sizes:
            L.append("   %.2f   " % (k / 100.0) + "".join("%-9s" % f(CNT(name, rg, k), "%.4f") for rg in LAD) +
                     "  %.4f/%.4f    %.4f/%.4f" % tuple(ex.MODEL[(n0, name, rg, k)][0] for rg in ("g2", "i40") for n0 in ex.MODELS))
        L += ["", "SET %s, B with the functional, eV" % name]
        L.append("   q, e   " + "".join("%-10s" % rg for rg in LAD) + " full      top minus full")
        for k in sizes:
            g = GR[(name, k)]
            L.append("   %.2f   " % (k / 100.0) + "".join("%-10s" % f(B(name, rg, k)) for rg in LAD) + " %+.4f   %s" % (
                ex.FULL[(name, k)][1], f(None if g["top"] is None else B(name, g["top"], k) - ex.FULL[(name, k)][1])))
        L.append("")
    L.append("   The largest difference between a measured count and the nearer of the two models: %s." % (
        "not formed" if worst_c is None else "%+.4f (%s, %.2f e, %s)" % (worst_c[0], worst_c[1], worst_c[2] / 100.0, worst_c[3])))
    L.append("   The largest difference between B on a converged grid and 'full': %s." % (
        "not formed" if worst_b is None else "%+.4f eV (%s, %.2f e)" % (worst_b[0], worst_b[1], worst_b[2] / 100.0)))

    head("R51.4  THE GRID AT EACH SIZE")
    L += ["CONVERGED if the two highest rungs that are measured agree to %.3f eV. THE LOWEST RUNG THAT SERVES:" % p51.TOL_GRID,
          "the first whose B, and that of every rung above it, lies within %.3f eV of the top one." % p51.TOL_GRID, ""]
    L.append("   set  q, e   top   next  top minus next, eV   the grid         lowest rung that serves   count at top")
    for name in sets:
        for k in sizes:
            g = GR[(name, k)]
            L.append("   %-4s %.2f   %-5s %-5s %-21s %-16s %-25s %s" % (name, k / 100.0, g["top"] or "-", g["nxt"] or "-", f(g["d"]), g["word"], g["low"] or "-",
                                                                     f(None if g["top"] is None else CNT(name, g["top"], k), "%.4f")))
    nconv = sum(1 for v in GR.values() if v["word"] == "CONVERGED")
    L += ["", "   The grid is CONVERGED at %d of %d (set, size), NOT CONVERGED at %d, not decided at %d." % (
        nconv, len(GR), sum(1 for v in GR.values() if v["word"] == "NOT CONVERGED"), sum(1 for v in GR.values() if v["word"] == "not decided"))]

    head("R51.5  THE FUNCTIONS: HARTREE-FOCK AGAINST EXACT")
    L += ["Hartree-Fock is exact for one electron and uses no grid, so this reads the functions alone.", "",
          "   q, e   exact      X2: B, HF  exact-HF   word           X3: B, HF  exact-HF   word"]
    for k in sizes:
        z = RULE[k]
        L.append("   %.2f   %-10s " % (k / 100.0, "%.4f" % z["exact"]) + "".join("%-10s %-10s %-14s " % (
            f(z["hf"][n]), f(None if z["hf"][n] is None else z["exact"] - z["hf"][n]), z["hfw"][n]) for n in sets).rstrip())

    head("R51.6  THE FUNCTIONS FOR THE FUNCTIONAL: SET X2 AGAINST SET X3")
    L += ["Where the grid is CONVERGED in both sets: what the functional holds in X2 minus what it holds in",
          "X3, each on its top rung. HOLDS if they agree to %.2f eV." % p51.TOL_ADEQ, "",
          "   q, e   B in X2    B in X3    X2 minus X3   word"]
    for k in sizes:
        z, g2, g3 = RULE[k], GR[("X2", k)], GR[("X3", k)]
        L.append("   %.2f   %-10s %-10s %-13s %s" % (k / 100.0, f(None if g2["top"] is None else B("X2", g2["top"], k)), f(None if g3["top"] is None else B("X3", g3["top"], k)),
                                                   f(z["fd"]), z["fw"]))

    head("R51.7  THE RULE")
    L += ["Set X2 is ADEQUATE at a size if the functions hold (R51.5), the grid is converged in both sets",
          "(R51.4) and the functional holds in X2 what it holds in X3 (R51.6). It is the set for the capture",
          "test if it is ADEQUATE at every size from %.2f e up." % (p51.FROM / 100.0), ""]
    for k in sizes:
        z = RULE[k]
        L.append("   %.2f e   %-13s %s" % (k / 100.0, z["word"], "; ".join(z["why"])))
    L.append("")
    if chosen == "X2":
        L.append("   THE FUNCTIONS FOR THE CAPTURE TEST: SET X2, with auxiliary treatment A.")
    elif chosen == "NONE":
        L.append("   THE FUNCTIONS FOR THE CAPTURE TEST: NONE. Set X2 is NOT ADEQUATE at a size from %.2f e up." % (p51.FROM / 100.0))
    else:
        L.append("   THE FUNCTIONS FOR THE CAPTURE TEST: NOT DECIDED. At a size from %.2f e up the rule could not be applied." % (p51.FROM / 100.0))
    if grid_rung is None:
        L.append("   THE GRID FOR A GHOST ATOM: NONE OF THE LADDER. At a size from %.2f e up in X2 no rung serves." % (p51.FROM / 100.0))
    else:
        L.append("   THE GRID FOR A GHOST ATOM: rung %s (IntAcc %.3f): the lowest that serves at every size from %.2f e up in X2." % (
            grid_rung, p51.INTACC[grid_rung], p51.FROM / 100.0))

    head("R51.8  THE FUNCTIONAL AGAINST EXACT: REPORTED, NOT PART OF THE RULE")
    L += ["What the functional holds in set X3 minus exact: with no RI and no COSX at IntAcc 40 (rung n40),",
          "and on the top rung of the ladder. 'computed' is the same difference in a complete basis on a",
          "full grid, computed by Claude in a second program. The grid word is X3's, from R51.4.", "",
          "   q, e   exact      n40 minus exact   top minus exact   computed   the grid"]
    fx = []
    for k in sizes:
        z, g3 = RULE[k], GR[("X3", k)]
        dn = None if B("X3", "n40", k) is None else B("X3", "n40", k) - z["exact"]
        dt = None if g3["top"] is None else B("X3", g3["top"], k) - z["exact"]
        if g3["word"] == "CONVERGED" and k >= p51.FROM and (dn if dn is not None else dt) is not None:
            fx.append((dn if dn is not None else dt, k))
        L.append("   %.2f   %-10s %-17s %-17s %+.4f    %s" % (k / 100.0, "%.4f" % z["exact"], f(dn), f(dt), ex.LIMIT[k][1] - z["exact"], g3["word"]))
    L.append("")
    if fx:
        w = max(fx, key=lambda t: abs(t[0]))
        L.append("   From %.2f e up, where the grid is converged, the functional lies at most %+.4f eV from exact (%.2f e)." % (p51.FROM / 100.0, w[0], w[1] / 100.0))
    else:
        L.append("   From %.2f e up the grid is converged at no size, so no such difference is formed." % (p51.FROM / 100.0))

    head("R51.9  WHAT RI AND COSX ADD")
    L += ["i40 minus n40: RI-J with AutoAux and COSX on its default grids, against exact integrals, both at",
          "IntAcc 40. c40 minus i40: the COSX grids raised to IntAcc 40. g2 minus n2: the same on the default grid.", "",
          "   set  q, e   i40 minus n40   c40 minus i40   g2 minus n2"]
    big = {"i40-n40": [], "c40-i40": [], "g2-n2": []}
    for name in sets:
        for k in sizes:
            d1 = None if None in (B(name, "i40", k), B(name, "n40", k)) else B(name, "i40", k) - B(name, "n40", k)
            d2 = None if None in (B(name, "c40", k), B(name, "i40", k)) else B(name, "c40", k) - B(name, "i40", k)
            d3 = None if None in (B(name, "g2", k), B(name, "n2", k)) else B(name, "g2", k) - B(name, "n2", k)
            for key, v in (("i40-n40", d1), ("c40-i40", d2), ("g2-n2", d3)):
                if v is not None:
                    big[key].append(abs(v))
            L.append("   %-4s %.2f   %-15s %-15s %s" % (name, k / 100.0, f(d1, "%+.5f"), f(d2, "%+.5f"), f(d3, "%+.5f")))
    L += ["", "   The largest: " + "; ".join("%s %s" % (key.replace("-", " minus "), ("%.5f eV" % max(v)) if v else "not formed") for key, v in big.items()) + "."]

    head("WHAT THIS CANNOT SAY (section 6 of the pre-registration, in short; its own words govern)")
    L += ["   One electron and a bare charge say what the functions, the grid and the functional do for one",
          "   electron. They do not say what a run of the dianion does.",
          "   The grid of a lone ghost atom is not the grid of a molecule, where the grids of the atoms",
          "   overlap. A rung that serves here is the setting to try on the ghost atoms of the dianion; that",
          "   it serves there is for the electron count of each run to show.",
          "   What Claude computed in a second program is an expectation and not part of the record.",
          "   D50's result stands as committed. D51 adopts nothing."]
    summary = ("Runs measured: %d of %d. D50's energies given again: %d of %d. The grid converged at %d of %d (set, size). "
               "The functions for the capture test, by the rule: %s. The grid for a ghost atom: %s." % (
                   sum(1 for x in RUN.values() if x["status"] == "measured"), len(RUN), nrep, len(REP), nconv, len(GR),
                   {"X2": "SET X2 with treatment A", "NONE": "NONE", "NOT DECIDED": "NOT DECIDED"}[chosen],
                   "none of the ladder" if grid_rung is None else "rung %s (IntAcc %.3f)" % (grid_rung, p51.INTACC[grid_rung])))
    L += ["", summary]

    T1 = ["run\tlevel\tset\trung\tsize_e\tstatus\tE_Eh\tB_eV\telectron_count\tbasis_functions\tIntAcc_printed\tgrid_points\tIntAcc_COSX_printed\tE_minus_d50_Eh"]
    rep = dict((r["run"], de) for r, de, dc in REP)
    for r in runs:
        x = RUN[r["run"]]
        T1.append("\t".join([r["run"], r["level"], r["set"], r["rung"], "%.2f" % (r["size"] / 100.0), x["status"], f(x["E"] if x["status"] == "measured" else None, "%.10f", ""),
                             f(x["B"] if x["status"] == "measured" else None, "%.6f", ""), f(x["count"], "%.6f", ""), "" if x["nb"] is None else str(x["nb"]),
                             f(x["xc"], "%.3f", ""), "" if x["pts"] is None else str(x["pts"]), ",".join("%.3f" % v for v in x["cosx"]),
                             f(rep.get(r["run"]), "%.3e", "")]))
    T2 = ["set\tsize_e\t" + "\t".join("count_%s\tB_%s_eV" % (rg, rg) for rg in rungs) + "\ttop\tnext\ttop_minus_next_eV\tgrid\tlowest_rung_that_serves\tB_full_expected_eV"]
    for name in sets:
        for k in sizes:
            g = GR[(name, k)]
            T2.append("\t".join([name, "%.2f" % (k / 100.0)] + [v for rg in rungs for v in (f(CNT(name, rg, k), "%.6f", ""), f(B(name, rg, k), "%.6f", ""))] +
                                [g["top"] or "", g["nxt"] or "", f(g["d"], "%.6f", ""), g["word"], g["low"] or "", "%.6f" % ex.FULL[(name, k)][1]]))
    T3 = ["size_e\texact_eV\tB_HF_X2_eV\tword_HF_X2\tB_HF_X3_eV\tword_HF_X3\tX2_minus_X3_functional_eV\tfunctional_holds\tX2_by_the_rule\tn40_X3_minus_exact_eV\tcomputed_limit_minus_exact_eV\tset_by_the_rule\tgrid_rung_by_the_rule"]
    for k in sizes:
        z = RULE[k]
        T3.append("\t".join(["%.2f" % (k / 100.0), "%.6f" % z["exact"], f(z["hf"]["X2"], "%.6f", ""), z["hfw"]["X2"], f(z["hf"]["X3"], "%.6f", ""), z["hfw"]["X3"],
                             f(z["fd"], "%.6f", ""), z["fw"], z["word"], f(None if B("X3", "n40", k) is None else B("X3", "n40", k) - z["exact"], "%.6f", ""),
                             "%.6f" % (ex.LIMIT[k][1] - z["exact"]), chosen, grid_rung or ""]))
    L = [l.rstrip() for l in L]
    for path, lines in ((a.result, L), (a.tsv, T1), (tsv2, T2), (tsv3, T3)):
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")
    sys.stdout.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
