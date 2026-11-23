#!/usr/bin/env python3
"""Stage 6, DIAGNOSTIC 53: charges of 0.25 e together, on the dianion, at two levels. THE READER.

Reads the runs of phase2.2/cs_stage6/d53 under the pre-registration STAGE6_DIAGNOSTIC53_20261123.txt and writes the
result and three tables. It was committed, with the pre-registration and the inputs, before any job was submitted.

Each run of the dianion is read as D48 and D50 read theirs: by the committed reader of outputs (cs6_d40.py), the
committed reading of stability (cs6_d33b.py) and the committed classifier (cs6_audit.py). The functions that form the
legs and the words (legs, word1, took, moved) are D50's, line for line, lifted out of its main() so that they can be
run on D50's own committed outputs and compared with D50's committed table.

Nothing is written unless every gate holds: the inputs are formed again from the committed files and compared, every
output must print the committed input beside it, and the committed files this leans on must be the committed ones.
D53 adopts nothing: it names no design and forms no ddE.

usage: python3 cs6_d53.py --res <results repository> --code <code repository> --prereg <the pre-registration>
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
PREREG = "STAGE6_DIAGNOSTIC53_20261123.txt"
D52PRE = "STAGE6_DIAGNOSTIC52_20261123.txt"
D52RES = "STAGE6_DIAGNOSTIC52_RESULT_20261123.txt"
D51PRE = "STAGE6_DIAGNOSTIC51_20261123.txt"
D50PRE = "STAGE6_DIAGNOSTIC50_20261122.txt"
EV = 27.211386245988
NEL = 118
COUNT_TOL = 0.02       # e: how far the electrons a grid integrates may lie from 118
L_PRINTED = "0.3127"   # eV: RULE C of D52, as its result prints it
L2_PRINTED = "0.1151"  # eV: three times wB97X-V's largest distance from the limit in D52
# quoted sentences: each must be in the document it is quoted from and in the pre-registration given
Q_A = "a run that takes electrons with a gain in energy not larger than a few times that error is not settled at this level alone."     # D51
Q_B = "TOOK ELECTRONS the SCF converged and either a leg fires or the run LEFT THE ATOMS"                                              # D50
Q_C = "RULE C GIVES: L = 0.3127 eV."                                                                                                   # D52's result
Q_D = "RULE B GIVES: wB97X-V."                                                                                                         # D52's result
Q_E = "with Hartree-Fock beside both as a bracket and not as a judge."                                                                 # D52
GUARDED = ["d37", "d44", "d50", "d52"]
DOCS = ["d52_rules.tsv", "d50_dianion.tsv", D52PRE, D52RES, D51PRE, D50PRE, "STAGE6_AMENDMENT10_20261123.txt"]
D50_DESIGN = "B_n11_q0.25"


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


def loewdin(text):
    """the charges of the last LOEWDIN ATOMIC CHARGES block, in the order printed (D48's reader)"""
    k = text.rfind("LOEWDIN ATOMIC CHARGES")
    q = []
    if k >= 0:
        for l in text[k:].split("\n")[2:]:
            m = re.match(r"\s*\d+\s*[A-Za-z]{1,2}\s*:\s*(-?\d+\.\d+)", l)
            if not m:
                break
            q.append(float(m.group(1)))
    return q


def nbf(text):
    nb = sorted(set(int(v) for v in re.findall(r"Number of basis functions\s*\.+\s*(\d+)", text)))
    return nb[0] if len(nb) == 1 else None


def count(text):
    """the electrons the grid integrates: the last N(Total) ORCA prints; None where it prints none"""
    v = re.findall(r"N\(Total\)\s*:\s*(-?\d+\.\d+)", text)
    return float(v[-1]) if v else None


def dist(p, q):
    return ((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2) ** 0.5


def read_di(d40, d33b, audit, d, kind):
    """one run of the dianion (D50's read_di)"""
    r = d40.read(d, kind, NEL, d33b)
    r["status"], r["c"], r["lw"], r["nbf"], r["count"] = "measured", None, None, None, None
    r.setdefault("stab", dict(ran=False, first=None, last=None))
    r.setdefault("trah", False)
    r.setdefault("lshift", False)
    if not r.get("exists"):
        r["status"] = "NOT MEASURED"
        return r
    r["nbf"] = nbf(r["text"])
    r["count"] = count(r["text"])
    if r["E"] is None:
        r["status"] = "DID NOT CONVERGE" if "SCF NOT CONVERGED" in r["text"] else "NOT MEASURED"
        return r
    c = audit.read(r["text"])
    lw = loewdin(r["text"])
    if c is None or c.get("maxlw") is None or c.get("dip") is None or c.get("nb") is None or len(lw) < 24:
        r["status"], r["E"] = "NOT MEASURED", None
        r["why"] = ["the committed classifier does not read a Loewdin charge of each atom, a dipole and a basis-function count from it"]
        return r
    r["c"], r["lw"] = c, lw
    return r


def moved(lw, lw0, atoms, charges):
    """the atom with the most negative change, that change, its distance to the nearest positive charge, the sum over the 24 (D50's)"""
    dq = [lw[i] - lw0[i] for i in range(24)]
    i = min(range(24), key=lambda j: dq[j])
    pos = [xyz for q, xyz in charges if q > 0]
    dn = min(dist(atoms[i][1:], p) for p in pos) if pos else None
    return "%s%d" % (atoms[i][0], i + 1), dq[i], dn, sum(lw[:24]) - sum(lw0[:24])


def legs(x, c0, atoms, charges, T):
    """the flags and holds of one dianion run against its control (D50's reader, which is D48's, line for line)"""
    D_LW, R_DIP, GAP_ABS, R_GAP = T
    z = dict(flags=[], held=[], rise=None, ratio=None, gap=None, gap_read=False, mv=None)
    if x["status"] == "DID NOT CONVERGE":
        z["flags"].append("did not converge")
        return z
    if x["status"] != "measured":
        return z
    g, kk = x["c"], c0["c"]
    z["rise"] = g["maxlw"] - kk["maxlw"]
    z["ratio"] = g["dip"] / kk["dip"] if kk["dip"] else None
    if z["rise"] >= D_LW:
        z["flags"].append("relocation: Loewdin rise %+.4f e" % z["rise"])
    if z["ratio"] is None:
        z["held"].append("the dipole ratio could not be formed")
    elif z["ratio"] >= R_DIP:
        z["flags"].append("dipole %.2fx control" % z["ratio"])
    z["gap"] = g.get("gap")
    if x["trah"] or x["lshift"] or c0["trah"] or c0["lshift"] or g.get("gap") is None or kk.get("gap") is None:
        z["held"].append("gap withheld")
    else:
        z["gap_read"] = True
        if g["gap"] < GAP_ABS and g["gap"] < R_GAP * kk["gap"]:
            z["flags"].append("gap %.4f eV, %.2f of control" % (g["gap"], g["gap"] / kk["gap"]))
    sb = x["stab"]
    if not sb["ran"] or sb["first"] is None or sb["last"] is None:
        z["held"].append("stability analysis not read")
    elif sb["last"] == "UNSTABLE":
        z["flags"].append("unstable at the last analysis")
    z["mv"] = moved(x["lw"], c0["lw"], atoms, charges)
    return z


NOZ = dict(flags=[], held=[], rise=None, ratio=None, gap=None, gap_read=False, mv=None)


def word1(x):
    if x["noread"]:
        return "NOT READ"
    if x["z"]["flags"]:
        return "FLAGGED"
    if x["status"] != "measured":
        return "NOT MEASURED"
    return "HELD" if x["z"]["held"] else "NOT FLAGGED"


def took(x, D_LW):
    mv = x["z"]["mv"]
    if x["noread"] or x["status"] != "measured" or mv is None:
        return "not decided", None
    if mv[3] >= D_LW or x["z"]["flags"]:
        return "TOOK ELECTRONS", mv[3]
    return ("DID NOT TAKE" if x["word"] == "NOT FLAGGED" else "not decided"), mv[3]


def control_holds(x):
    """D50's word for a control with functions read against the plain control: the three legs that do not divide the density"""
    hard = [f for f in x["z"]["flags"] if not f.startswith("relocation")]
    if x["noread"] or (x["status"] != "measured" and not hard):
        return "not decided"
    if hard:
        return "DOES NOT HOLD"
    return "not decided" if x["z"]["held"] else "HOLDS"


def stable(x):
    return x["status"] == "measured" and bool(x["stab"].get("ran")) and x["stab"].get("last") == "STABLE"


def together(c, w, second):
    """R53.6: (the word, its mark, why it is NOT DECIDED) for one (arrangement, state), from its cell at each level"""
    if not second:
        return "NOT DECIDED", "", "the second level is not read"
    for lv, k in (("CAM", c), ("WXV", w)):
        if k["did"] == "not decided":
            return "NOT DECIDED", "", "%s: %s" % (lv, k["why"] or "not decided")
        if k["holds"] != "HOLDS":
            return "NOT DECIDED", "", "%s: its control %s %s" % (lv, k["control"], "DOES NOT HOLD" if k["holds"] == "DOES NOT HOLD" else "is not decided")
        if not k["count_ok"]:
            return "NOT DECIDED", "", "%s: a count does not hold" % lv
    ct, wt = c["did"] == "TOOK ELECTRONS", w["did"] == "TOOK ELECTRONS"
    if ct and wt:
        if w["gain"] is None:
            return "NOT DECIDED", "", "WXV: the gain is not formed (%s)" % w["nogain"]
        return "TAKES ELECTRONS AT BOTH LEVELS", "%s L2 at WXV" % w["bar"], ""
    if ct:
        if c["gain"] is None:
            return "NOT DECIDED", "", "CAM: the gain is not formed (%s)" % c["nogain"]
        return "TAKEN AT THE WORKING FUNCTIONAL ONLY", "%s L" % c["bar"], ""
    if wt:
        return "TAKEN AT THE SECOND LEVEL ONLY", "", ""
    return "HOLDS AT BOTH LEVELS", "", ""


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    a.res, a.code, a.prereg = (os.path.abspath(x) for x in (a.res, a.code, a.prereg))
    S6 = os.path.join(a.res, S)
    D = os.path.join(S6, "d53")
    paths = dict((k, os.path.join(a.code, "phase2.2", fn)) for k, fn in (
        ("d40", "cs6_d40.py"), ("d33b", "cs6_d33b.py"), ("audit", "cs6_audit.py"), ("p53", "cs6_d53_prepare.py")))
    for k, p in paths.items():
        if not os.path.isfile(p):
            stop("the committed script %s is missing." % p)
    for fn in [a.prereg, os.path.join(D, "RUNS.tsv")] + [os.path.join(S6, x) for x in DOCS]:
        if not os.path.isfile(fn):
            stop("%s is missing." % fn)
    tsv2, tsv3 = os.path.join(os.path.dirname(a.tsv), "d53_words.tsv"), os.path.join(os.path.dirname(a.tsv), "d53_controls.tsv")
    for fn in (a.result, a.tsv, tsv2, tsv3):
        if os.path.exists(fn):
            stop("%s exists already; this reader overwrites nothing." % fn)
    d40, d33b, audit = load(paths["d40"], "d53_cs6_d40"), load(paths["d33b"], "d53_cs6_d33b"), load(paths["audit"], "d53_cs6_audit")
    p53 = load(paths["p53"], "d53_cs6_d53_prepare")
    T = (audit.D_LW, audit.R_DIP, audit.GAP_ABS, audit.R_GAP)
    D_LW = audit.D_LW
    gates = []

    def git(*args):
        try:
            q = subprocess.run(["git", "-C", a.res] + list(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            return q.returncode, q.stdout.strip()
        except OSError:
            return 127, ""

    # ------------------------------------------------------------ G53.3: the committed files this leans on
    guard = ["%s/%s" % (S, g) for g in GUARDED] + ["%s/%s" % (S, x) for x in DOCS]
    rc, st = git("status", "--porcelain", "--untracked-files=all", "--", *guard)
    if rc != 0:
        gates.append("G53.3 could not be checked: git did not run")
    elif st:
        gates.append("G53.3 a file under %s, or a committed table or document this reader names, is not the committed one" % ", ".join(GUARDED))
    rc, pcm = git("log", "--diff-filter=A", "--format=%H", "-1", "--", "%s/%s" % (S, PREREG))
    if rc != 0 or not pcm:
        gates.append("G53.3 the commit that added %s could not be found" % PREREG)
    else:
        rc1, b1 = git("rev-parse", "%s:%s/%s" % (pcm, S, PREREG))
        rc2, b2 = git("hash-object", a.prereg)
        if rc1 != 0 or rc2 != 0 or b1 != b2:
            gates.append("G53.3 the pre-registration given is not the committed %s" % PREREG)

    # ------------------------------------------------------------ G53.1 and G53.2: formed again from the committed files
    try:
        F, runs = p53.build(S6, a.code)
    except Exception as e:
        sys.stdout.write("\n".join("FAILED  " + g for g in gates + ["G53.1 the files of d53 could not be formed again from the committed files: %s: %s" % (type(e).__name__, e)]) + "\n")
        stop("%d gate(s) failed. Nothing was written." % (len(gates) + 1))
    psrc = os.path.join(D, "SOURCES.sha256")
    if not os.path.isfile(psrc) or open(psrc).read() != F["SOURCES.sha256"]:
        gates.append("G53.2 d53/SOURCES.sha256 is not the list formed again from the committed files")
    nsrc = 0
    for l in F["SOURCES.sha256"].splitlines():
        h, fn = l.split(None, 1)
        nsrc += 1
        p = os.path.join(S6, fn)
        if not os.path.isfile(p) or sha(p) != h:
            gates.append("G53.2 %s is not the file d53/SOURCES.sha256 lists" % fn)
    nsum = 0
    pic = os.path.join(D, "INPUT_CHECKSUMS.sha256")
    if not os.path.isfile(pic):
        gates.append("G53.1 d53/INPUT_CHECKSUMS.sha256 is missing")
    else:
        if open(pic).read() != F["INPUT_CHECKSUMS.sha256"]:
            gates.append("G53.1 d53/INPUT_CHECKSUMS.sha256 is not the list formed again here")
        for l in open(pic).read().splitlines():
            if l.strip():
                h, fn = l.split(None, 1)
                nsum += 1
                p = os.path.join(D, fn)
                if not os.path.isfile(p) or sha(p) != h:
                    gates.append("G53.1 d53/%s is not the file the checksum list holds" % fn)
    for k in sorted(F):
        p = os.path.join(D, k)
        if not os.path.isfile(p) or open(p).read() != F[k]:
            gates.append("G53.1 d53/%s is not the file formed again from the committed files" % k)
    for r in runs:
        for fn in ("pc.pc", "substrate.xyz"):
            if os.path.exists(os.path.join(D, r["run"], fn)) and ("%s/%s" % (r["run"], fn)) not in F:
                gates.append("G53.1 %s holds a file %s that its input does not read" % (r["run"], fn))

    # ------------------------------------------------------------ G53.4: the quoted sentences, and the two bars
    def has(path, sent):
        return " ".join(sent.split()) in " ".join(open(path, errors="replace").read().split())
    for doc, q, what in ((D51PRE, Q_A, "D51's sentence on a gain within the functional's error"), (D50PRE, Q_B, "the reading TOOK ELECTRONS"),
                         (D52RES, Q_C, "RULE C of D52"), (D52RES, Q_D, "RULE B of D52"), (D52PRE, Q_E, "D52's sentence on Hartree-Fock")):
        if not has(os.path.join(S6, doc), q):
            gates.append("G53.4 %s is not in %s" % (what, doc))
        if not has(a.prereg, q):
            gates.append("G53.4 the words \"%s\" are not in the pre-registration given" % (q[:60] + " ..."))
    rules = dict(l.split("\t")[:2] for l in open(os.path.join(S6, "d52_rules.tsv")).read().split("\n")[1:] if "\t" in l)
    try:
        L = float(rules["rule_C_L_eV"])
        L2 = 3.0 * float(rules["largest_distance_from_the_limit_WXV_eV"])
    except (KeyError, ValueError):
        L = L2 = None
    if L is None or "%.4f" % L != L_PRINTED or "%.4f" % L2 != L2_PRINTED or rules.get("rule_B_second_level") != "wB97X-V" \
            or rules.get("rule_A_set") != "X2" or rules.get("rule_A_IntAcc") != "12.0":
        gates.append("G53.4 d52_rules.tsv does not hold the rules the pre-registration names (X2, IntAcc 12, wB97X-V, L %s eV, L2 %s eV)" % (L_PRINTED, L2_PRINTED))

    # ------------------------------------------------------------ reading every run, and G53.5, G53.6
    geom = dict((stt, p53.geometry(S6, stt)) for stt in p53.STATES)
    arr = dict(p53.ARR)
    RUN = {}
    for r in runs:
        d = os.path.join(D, r["run"])
        x = read_di(d40, d33b, audit, d, "pc" if r["kind"] in ("F", "N") else "plain")
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
                gates.append("G53.5 %s/job.out ends normally and does not print the input it ran" % r["run"])
        elif want is not None and e != want:
            gates.append("G53.5 the input printed at the top of %s/job.out is not the committed %s/job.inp" % (r["run"], r["run"]))
        if x["status"] == "measured" and x["nbf"] != r["nb"]:
            gates.append("G53.6 %s reports %s basis functions; its input asks for %d" % (r["run"], x["nbf"], r["nb"]))

    if gates:
        sys.stdout.write("\n".join("FAILED  " + g for g in gates) + "\n")
        stop("%d gate(s) failed. Nothing was written." % len(gates))

    # ------------------------------------------------------------ the count
    for x in RUN.values():
        if x["level"] == "HF":
            x["count_ok"], x["count_d"] = True, None
        elif x["status"] != "measured" or x["count"] is None:
            x["count_ok"], x["count_d"] = False, None
        else:
            x["count_d"] = x["count"] - NEL
            x["count_ok"] = abs(x["count_d"]) <= COUNT_TOL

    # ------------------------------------------------------------ the plain controls, the controls with functions, the partners
    for x in RUN.values():
        x["noread"], x["z"], x["word"], x["took"], x["how"], x["gain"], x["bar"], x["ghost"], x["holds"], x["nogain"] = None, dict(NOZ), "", ("not decided", None), "", None, "", None, "", ""
    for x in RUN.values():
        if x["kind"] == "GN":
            x["word"] = "STABLE" if stable(x) else ("NOT MEASURED" if x["status"] != "measured" else "NOT STABLE")
    for x in RUN.values():
        if x["kind"] not in ("GF", "N"):
            continue
        c0 = RUN[x["control"]]
        if not stable(c0):
            x["noread"] = "its plain control %s is %s" % (x["control"], c0["word"])
        x["z"] = legs(x, c0, geom[x["state"]][0], p53.charges_of(S6, arr[x["arr"]], x["state"]) if x["kind"] == "N" else [], T) if not x["noread"] else dict(NOZ)
        x["word"] = word1(x)
        if x["kind"] == "GF":
            x["holds"] = control_holds(x)
            mv = x["z"]["mv"]
            x["left"] = mv[3] if mv else None
            x["ghost"] = sum(x["lw"][24:30]) if x["status"] == "measured" and x["lw"] and len(x["lw"]) >= 30 else None
            x["dE"] = (x["E"] - c0["E"]) * EV if (x["status"] == "measured" and c0["status"] == "measured") else None
    # ------------------------------------------------------------ each arrangement with the functions
    for x in RUN.values():
        if x["kind"] != "F":
            continue
        c0 = RUN[x["control"]]
        why = None
        if c0["status"] != "measured":
            why = "its control %s is %s" % (x["control"], c0["status"])
        elif c0["noread"]:
            why = "its control %s is not read: %s" % (x["control"], c0["noread"])
        elif not c0["stab"].get("ran") or c0["stab"].get("last") is None:
            why = "the stability analysis of its control %s was not read" % x["control"]
        elif c0["stab"]["last"] != "STABLE":
            why = "its control %s is UNSTABLE at the last analysis" % x["control"]
        x["noread"] = why
        ch = p53.charges_of(S6, arr[x["arr"]], x["state"])
        x["z"] = legs(x, c0, geom[x["state"]][0], ch, T) if not x["noread"] else dict(NOZ)
        x["word"] = word1(x)
        x["took"] = took(x, D_LW)
        x["how"] = ("ONTO THE FUNCTIONS" if x["took"][1] >= D_LW else "A LEG ONLY") if x["took"][0] == "TOOK ELECTRONS" else ""
        if x["status"] == "measured" and x["lw"] and c0.get("lw") and len(x["lw"]) >= 30 and len(c0["lw"]) >= 30 and not x["noread"]:
            x["ghost"] = sum(x["lw"][24:30]) - sum(c0["lw"][24:30])
        if x["partner"] != "-":
            n = RUN[x["partner"]]
            gn = RUN[n["control"]]
            four = (x, c0, n, gn)
            if all(v["status"] == "measured" for v in four) and all(v["count_ok"] for v in four):
                x["gain"] = ((x["E"] - c0["E"]) - (n["E"] - gn["E"])) * EV
            else:
                x["nogain"] = "; ".join("%s %s" % (v["run"], v["status"] if v["status"] != "measured" else "has a count that does not hold") for v in four
                                        if v["status"] != "measured" or not v["count_ok"])
        if x["took"][0] == "TOOK ELECTRONS" and x["level"] in ("CAM", "WXV"):
            bar = L if x["level"] == "CAM" else L2
            x["bar"] = "not formed" if x["gain"] is None else ("ABOVE" if -x["gain"] > bar else "WITHIN")

    # ------------------------------------------------------------ R53.2 and R53.6
    wctl = [RUN[p53.name_of("WXV", k, "-", stt)] for stt in p53.STATES for k in ("GN", "GF")]
    second_why = ["%s is %s" % (x["run"], x["status"] if x["status"] != "measured" else ("without a stability analysis that was read" if not x["stab"].get("ran") or x["stab"].get("last") is None else "UNSTABLE at its last analysis"))
                  for x in wctl if not stable(x)]
    second = not second_why

    def cell(level, an, stt):
        x = RUN[p53.name_of(level, "F", an, stt)]
        c0 = RUN[x["control"]]
        why = x["noread"] or (None if x["status"] == "measured" else ("the run DID NOT CONVERGE" if x["status"] == "DID NOT CONVERGE" else "the run is NOT MEASURED")) or (", ".join(x["z"]["held"]) if x["took"][0] == "not decided" else None)
        return dict(did=x["took"][0], left=x["took"][1], gain=x["gain"], bar=x["bar"], count_ok=x["count_ok"] and c0["count_ok"], why=why, how=x["how"], x=x,
                    holds=c0["holds"], control=c0["run"], nogain=x["nogain"])
    WORDS = []
    for an, sites in p53.ARR:
        for stt in p53.STATES:
            c, w, h = cell("CAM", an, stt), cell("WXV", an, stt), cell("HF", an, stt)
            word, mark, reason = together(c, w, second)
            WORDS.append(dict(arr=an, state=stt, c=c, w=w, h=h, word=word, mark=mark, reason=reason))

    # ------------------------------------------------------------ the result
    rep = []
    say = rep.append

    def f(v, fmt="%+.4f", none="-"):
        return none if v is None else fmt % v

    def para(text, lead="", width=96):
        for l in textwrap.wrap(" ".join(text.split()), width=width - len(lead), break_long_words=False, break_on_hyphens=False):
            say(lead + l)

    def title(t):
        say("")
        say("=" * 96)
        say(t)
        say("=" * 96)
        say("")
    rc, head = git("rev-parse", "HEAD")
    say("STAGE 6 - DIAGNOSTIC 53 - RESULT")
    say("Charges of 0.25 e together, on the dianion, at two levels")
    say("=" * 96)
    say("")
    para("Read by cs6_d53.py under %s, which was committed before any job was submitted." % PREREG)
    say("   results repository at %s; the pre-registration was added by %s" % (head[:12] if rc == 0 else "?", pcm[:12]))
    for lab, k in (("reader", None), ("inputs formed again", "p53"), ("outputs read by", "d40"), ("stability by", "d33b"), ("the four legs by", "audit")):
        p = os.path.abspath(__file__) if k is None else paths[k]
        say("   %-20s %s" % (lab, os.path.basename(p)))
        say("   %-20s sha256 %s" % ("", sha(p)))
    say("")
    para("""The dianion of frame 42436 at R and TS, frozen, in vacuo, def2-TZVPD, with bare charges of 0.25 e at sites of
B_n11_q0.25. F: the charges, and the functions of D50 (set X2, treatment A) at all six positive sites. N: the same charges and
no functions. GF: the functions and no charge. GN: neither. CAM is CAM-B3LYP-D3BJ, WXV is wB97X-V, HF is Hartree-Fock. In every
run of CAM and WXV hydrogen has IntAcc 12. LEFT THE ATOMS: the Loewdin charges of the 24 atoms sum to %.2f e or more above the
control's. The gain: (F minus GF) minus (N minus GN), eV; below zero, the functions lower the energy more when the charges are
there. L = %.4f eV (D52, RULE C); L2 = %.4f eV. D53 adopts nothing: it names no design and forms no ddE.""" % (D_LW, L, L2))
    say("")
    say("   arrangement   charges of 0.25 e at (minus sign: -0.25 e)")
    for an, sites in p53.ARR:
        order = [s for s in p53.PC_SITES if s in sites]
        say("   %-13s %s" % (an, ", ".join(("%d" if s in p53.POS else "-%d") % s for s in order)))
    say("")
    say("THE GATES, as they held:")

    def item(label, text, ind=10):
        for i, l in enumerate(textwrap.wrap(" ".join(text.split()), width=96 - ind, break_long_words=False, break_on_hyphens=False)):
            say((("   " + label).ljust(ind) if i == 0 else " " * ind) + l)
    item("G53.1", "all %d files of d53/INPUT_CHECKSUMS.sha256 are the ones listed, and every file of d53 is the one formed again here from the committed files." % nsum)
    item("G53.2", "the %d files of d53/SOURCES.sha256 are the ones listed." % nsrc)
    item("G53.3", "no file under %s or %s, and no committed table or document named here, differs from the committed one; the pre-registration given is the committed one." % (", ".join(GUARDED[:-1]), GUARDED[-1]))
    item("G53.4", "the sentences quoted here are in the documents they are quoted from, and d52_rules.tsv holds the rules named.")
    item("G53.5", "every output that exists prints, as the input it ran, the committed input beside it.")
    item("G53.6", "every measured run reports the number of basis functions its input asks for.")

    # ---- R53.1
    title("R53.1  THE RUNS, AND THE COUNT")
    say("   level   kind   measured")
    for lv, _ in p53.LEVELS:
        for kind in p53.KINDS[lv]:
            xs = [x for x in RUN.values() if x["level"] == lv and x["kind"] == kind]
            say("   %-7s %-6s %2d of %2d" % (lv, kind, sum(1 for x in xs if x["status"] == "measured"), len(xs)))
    say("")
    bad = [RUN[r["run"]] for r in runs if RUN[r["run"]]["status"] != "measured"]
    if bad:
        say("   Measured: %d of %d. The others, with ORCA's words:" % (len(runs) - len(bad), len(runs)))
        groups = []
        for x in bad:
            key = (x["status"], "; ".join(x.get("why") or ["no output"]))
            g = next((g for g in groups if g[0] == key), None)
            if g is None:
                groups.append((key, [x["run"]]))
            else:
                g[1].append(x["run"])
        for (status, words), names in groups:
            para("%s (%d run%s): %s" % (status, len(names), "" if len(names) == 1 else "s", words), "   ")
            para(", ".join(names), "      ")
    else:
        say("   Every one of the %d runs is measured." % len(runs))
    say("")
    para("THE COUNT, the electrons the grid integrates (N(Total)). It HOLDS within %.2f e of %d. Hartree-Fock uses no grid of this kind and has no count." % (COUNT_TOL, NEL), "   ")
    cs = [x for x in RUN.values() if x["level"] != "HF" and x["status"] == "measured"]
    off = [x for x in cs if not x["count_ok"]]
    if cs and not off:
        far = max(cs, key=lambda x: abs(x["count_d"]))
        say("   The count holds in every measured run of CAM and WXV (%d)." % len(cs))
        say("   The farthest from %d: %s, %+.4f e." % (NEL, far["run"], far["count_d"]))
    else:
        say("   The count holds in %d of %d measured runs of CAM and WXV. It does not hold in:" % (len(cs) - len(off), len(cs)))
        if not cs:
            say("   (no run of CAM or WXV is measured)")
        for x in off:
            say("   %s: %s" % (x["run"], "no count printed" if x["count_d"] is None else "%+.4f e from %d" % (x["count_d"], NEL)))

    # ---- R53.2
    title("R53.2  CAN THE SECOND LEVEL BE READ")
    para("The four controls of wB97X-V, each with the stability analysis. This record had not run wB97X-V with it before.")
    say("")
    say("   control          status         first analysis   last analysis")
    for x in wctl:
        say("   %-16s %-14s %-16s %s" % (x["run"], x["status"], x["stab"].get("first") or "-", x["stab"].get("last") or "-"))
    say("")
    if second:
        say("   THE SECOND LEVEL IS READ.")
    else:
        say("   THE SECOND LEVEL IS NOT READ:")
        for w_ in second_why:
            para(w_ + ".", "   ")

    # ---- R53.3
    title("R53.3  THE CONTROLS")
    para("""The control with functions (GF) against the plain control (GN) of its level and state, by the dipole, the gap and the
stability. 'missing': the Loewdin charge the 24 atoms hold less than in the plain control, e. 'on the functions': the Loewdin
charge of the six ghost atoms, e (below zero: electrons). 'energy': GF minus GN, eV.""")
    say("")
    say("   level  state  plain control   with functions   missing    on the functions   energy")
    CTL = []
    for lv, _ in p53.LEVELS:
        for stt in p53.STATES:
            gn, gf = RUN[p53.name_of(lv, "GN", "-", stt)], RUN[p53.name_of(lv, "GF", "-", stt)]
            say("   %-6s %-6s %-15s %-16s %-10s %-18s %s" % (lv, stt, gn["word"], gf["holds"] if not gf["noread"] else "NOT READ", f(gf.get("left")), f(gf["ghost"]), f(gf.get("dE"))))
            CTL.append((lv, stt, gn, gf))
    notes = ["%s %s" % (gf["run"], "; ".join(gf["z"]["flags"] + gf["z"]["held"])) for lv, stt, gn, gf in CTL if gf["z"]["flags"] or gf["z"]["held"]]
    if notes:
        say("")
        for n in notes:
            para(n, "   ")

    # ---- R53.4
    title("R53.4  EACH ARRANGEMENT AT EACH LEVEL AND STATE, AND THE BAR (R53.5)")
    para("""Each run with the charges and the functions (F) against its control with the functions (GF), as D48 and D50 read
theirs: "%s". 'left': the charge that left the 24 atoms, e. 'functions': the change of the Loewdin charge of the six ghost
atoms, e. 'gain': eV. The bar is read only for a run that took electrons: WITHIN if it gained no more than L (CAM) or L2
(WXV), ABOVE if more.""" % Q_B)
    for lv, _ in p53.LEVELS:
        say("")
        say("   %s (%s)" % (lv, p53.LEVEL_NAME[lv]))
        say("   arr.  state  the run         how                 left      functions  gain      the bar")
        for an, sites in p53.ARR:
            for stt in p53.STATES:
                x = RUN[p53.name_of(lv, "F", an, stt)]
                say(("   %-5s %-6s %-15s %-19s %-9s %-10s %-9s %s" % (an, stt, x["took"][0] if not x["noread"] else "NOT READ", x["how"] or "-", f(x["took"][1]), f(x["ghost"]),
                     f(x["gain"]), x["bar"] or "-")).rstrip())
        why = ["%s: %s" % (x["run"], x["noread"] or "; ".join(x["z"]["flags"] + x["z"]["held"])) for x in RUN.values()
               if x["level"] == lv and x["kind"] == "F" and (x["noread"] or x["z"]["flags"] or x["z"]["held"])]
        for wl in why:
            para(wl, "      ")
    pn = [x for x in RUN.values() if x["kind"] == "N" and x["status"] == "measured" and x["word"] != "NOT FLAGGED"]
    say("")
    if pn:
        para("The partners without functions (N) that are not NOT FLAGGED against their plain control: " + "; ".join("%s %s (%s)" % (x["run"], x["word"], "; ".join(x["z"]["flags"] + x["z"]["held"]) or x["noread"]) for x in pn) + ".", "   ")
    else:
        para("Every measured partner without functions (N) is NOT FLAGGED against its plain control.", "   ")

    # ---- R53.6
    title("R53.6  THE TWO LEVELS TOGETHER")
    para("""D51's pre-registration: "%s" The words are those of R53.6 of the pre-registration, and NOT DECIDED is taken
first: the second level not read; or, at CAM or at WXV, the F run not decided, its control not HOLDING, or the count of F or of
its control not holding; or the gain not formed where the word needs its bar. TAKES ELECTRONS AT BOTH LEVELS is marked by the
gain at WXV against L2; TAKEN AT THE WORKING FUNCTIONAL ONLY by the gain at CAM against L.""" % Q_A)
    say("")
    say("   arr.  state  CAM           WXV           the two levels together")
    short = {"TOOK ELECTRONS": "TOOK", "DID NOT TAKE": "DID NOT TAKE", "not decided": "not decided"}
    for wd in WORDS:
        say("   %-5s %-6s %-13s %-13s %s%s" % (wd["arr"], wd["state"], short[wd["c"]["did"]], short[wd["w"]["did"]], wd["word"], (" (%s)" % wd["mark"]) if wd["mark"] else ""))
    nd = [wd for wd in WORDS if wd["word"] == "NOT DECIDED"]
    if nd:
        say("")
        say("   NOT DECIDED, because:")
        seen = []
        for wd in nd:
            if wd["reason"] == "the second level is not read":
                if wd["reason"] not in seen:
                    seen.append(wd["reason"])
                    para("every one: the second level is not read (R53.2).", "   ")
            else:
                para("%s at %s: %s." % (wd["arr"], wd["state"], wd["reason"]), "   ")
    say("")
    order = ["HOLDS AT BOTH LEVELS", "TAKES ELECTRONS AT BOTH LEVELS", "TAKEN AT THE WORKING FUNCTIONAL ONLY", "TAKEN AT THE SECOND LEVEL ONLY", "NOT DECIDED"]
    tally = "; ".join("%s %d" % (o, sum(1 for wd in WORDS if wd["word"] == o)) for o in order)
    para("Of the %d (arrangement, state): %s." % (len(WORDS), tally), "   ")
    both = [wd for wd in WORDS if wd["word"] == "TAKES ELECTRONS AT BOTH LEVELS"]
    if both:
        n1 = min(len(arr[wd["arr"]]) for wd in both)
        para("The fewest charges in an arrangement that takes electrons at both levels: %d (%s)." % (n1, ", ".join("%s at %s" % (wd["arr"], wd["state"]) for wd in both if len(arr[wd["arr"]]) == n1)), "   ")
    dd = [wd for wd in WORDS if wd["arr"] == "d11"]
    design = "; ".join("at %s, %s%s" % (wd["state"], wd["word"], (" (%s)" % wd["mark"]) if wd["mark"] else "") for wd in dd)
    para("The design (d11): %s." % design, "   ")

    # ---- R53.7
    title("R53.7  HARTREE-FOCK, THE BRACKET")
    para("""D52's pre-registration runs the second level beside the working functional "%s" What each run with the charges and
the functions did at Hartree-Fock against its own control. It changes no word of R53.6.""" % Q_E)
    say("")
    say("   arr.  state  the run         how                 left      functions")
    for wd in WORDS:
        x = wd["h"]["x"]
        say("   %-5s %-6s %-15s %-19s %-9s %s" % (wd["arr"], wd["state"], x["took"][0] if not x["noread"] else "NOT READ", x["how"] or "-", f(x["took"][1]), f(x["ghost"])))

    # ---- R53.8
    title("R53.8  THE DESIGN AGAINST D50")
    para("""D50 ran the design with the same functions on ORCA's default grid; here hydrogen, and with it each ghost atom, has
IntAcc 12. From the committed d50_dianion.tsv (treatment A) and from R53.4 (CAM, d11):""")
    say("")
    d50 = [dict(zip(l0.split("\t"), l.split("\t"))) for l0 in [open(os.path.join(S6, "d50_dianion.tsv")).read().split("\n")[0]]
           for l in open(os.path.join(S6, "d50_dianion.tsv")).read().split("\n")[1:] if l.strip()]
    say("   state  D50: what it did    left       gain       here: what it did   left       gain")
    for stt in p53.STATES:
        o = [r for r in d50 if r.get("condition") == D50_DESIGN and r.get("aux") == "A" and r.get("state") == stt]
        x = RUN[p53.name_of("CAM", "F", "d11", stt)]
        o = o[0] if len(o) == 1 else {}
        say("   %-6s %-19s %-10s %-10s %-19s %-10s %s" % (stt, o.get("what_the_run_did", "-"), o.get("sum_24_minus_control_e", "-")[:7], o.get("gain_eV", "-")[:7],
            x["took"][0] if not x["noread"] else "NOT READ", f(x["took"][1]), f(x["gain"])))

    title("WHAT THIS CANNOT SAY (section 6 of the pre-registration, in short; its own words govern)")
    for t in ("One frame, frozen geometries, in vacuo, bare charges of one size, seven arrangements at sites of one design. The fifth leg (the monoanion) is not run.",
              "D52 measured how strongly each functional holds one whole electron at the charges. How firmly each has the dianion hold its own electron is the other side of a capture; D39's coupled cluster gives it by two routes that lie 0.7 eV apart. Two levels that agree are not thereby right. L and L2 are bars for one whole electron.",
              "A Loewdin charge of a diffuse function divides the density; it does not count electrons that moved. Each run is read against a control with the same functions, and the gain is given beside it.",
              "D50's, D51's and D52's results stand as committed. D53 adopts nothing."):
        para(t, "   ")
    say("")
    nm = sum(1 for x in RUN.values() if x["status"] == "measured")
    para("Runs measured: %d of %d. %s Of the %d (arrangement, state): %s. The design: %s." % (
        nm, len(runs), "THE SECOND LEVEL IS READ." if second else "THE SECOND LEVEL IS NOT READ.", len(WORDS), tally,
        design))

    # ------------------------------------------------------------ the tables
    def tv(v, fmt="%.6f"):
        return "" if v is None else fmt % v
    cols = ["run", "job", "level", "kind", "arrangement", "state", "status", "E_Eh", "count", "count_minus_118", "count_holds", "basis_functions",
            "stability_first", "stability_last", "lowdin_rise_e", "dipole_ratio", "gap_eV", "flags", "held", "sum_24_minus_control_e",
            "ghosts_minus_control_e", "word", "what_the_run_did", "how_it_took", "gain_eV", "the_bar", "not_read_because"]
    out = ["\t".join(cols)]
    for r in runs:
        x = RUN[r["run"]]
        z = x["z"]
        isf = x["kind"] == "F"
        out.append("\t".join([x["run"], x["job"], x["level"], x["kind"], x["arr"], x["state"], x["status"], tv(x["E"], "%.10f"), tv(x["count"], "%.6f"), tv(x["count_d"], "%+.6f"),
                              "" if x["level"] == "HF" else ("yes" if x["count_ok"] else "no"), "" if x["nbf"] is None else str(x["nbf"]),
                              x["stab"].get("first") or "" if x.get("exists") else "", x["stab"].get("last") or "" if x.get("exists") else "",
                              tv(z["rise"], "%+.6f"), tv(z["ratio"], "%.4f"), tv(z["gap"], "%.4f"), "; ".join(z["flags"]), "; ".join(z["held"]),
                              tv(z["mv"][3] if z["mv"] else None, "%+.6f"), tv(x["ghost"], "%+.6f"), x["word"],
                              (x["took"][0] if isf else ""), (x["how"] if isf else ""), tv(x["gain"], "%+.6f"), x["bar"], x["noread"] or ""]))
    w2 = ["\t".join(["arrangement", "state", "charges", "CAM_did", "CAM_left_e", "CAM_gain_eV", "CAM_bar", "WXV_did", "WXV_left_e", "WXV_gain_eV", "WXV_bar",
                     "the_two_levels_together", "mark", "not_decided_because", "HF_did", "HF_left_e", "L_eV", "L2_eV"])]
    for wd in WORDS:
        c, w, h = wd["c"], wd["w"], wd["h"]
        w2.append("\t".join([wd["arr"], wd["state"], str(len(arr[wd["arr"]])), c["did"], tv(c["left"], "%+.6f"), tv(c["gain"], "%+.6f"), c["bar"], w["did"], tv(w["left"], "%+.6f"),
                             tv(w["gain"], "%+.6f"), w["bar"], wd["word"], wd["mark"], wd["reason"], h["did"], tv(h["left"], "%+.6f"), "%.6f" % L, "%.6f" % L2]))
    w3 = ["\t".join(["level", "state", "plain_control", "plain_E_Eh", "with_functions", "with_functions_E_Eh", "missing_from_the_24_atoms_e", "on_the_functions_e",
                     "GF_minus_GN_eV", "flags", "held", "count_GN_minus_118", "count_GF_minus_118", "the_second_level_is_read"])]
    for lv, stt, gn, gf in CTL:
        w3.append("\t".join([lv, stt, gn["word"], tv(gn["E"], "%.10f"), gf["holds"] if not gf["noread"] else "NOT READ", tv(gf["E"], "%.10f"), tv(gf.get("left"), "%+.6f"),
                             tv(gf["ghost"], "%+.6f"), tv(gf.get("dE"), "%+.6f"), "; ".join(gf["z"]["flags"]), "; ".join(gf["z"]["held"]),
                             tv(gn["count_d"], "%+.6f"), tv(gf["count_d"], "%+.6f"), "yes" if second else "no"]))
    text = "\n".join(l.rstrip() for l in rep) + "\n"
    for p, t in ((a.result, text), (a.tsv, "\n".join(out) + "\n"), (tsv2, "\n".join(w2) + "\n"), (tsv3, "\n".join(w3) + "\n")):
        with open(p, "w") as fh:
            fh.write(t)
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
