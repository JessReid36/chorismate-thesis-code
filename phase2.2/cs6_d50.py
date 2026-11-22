#!/usr/bin/env python3
"""cs6_d50.py - the reader of Stage 6 DIAGNOSTIC 50: functions at a charge that can hold an electron at 0.25 e.
Committed before any job of D50 is submitted.

PART C, one electron. Each run is read with the committed reader of outputs (cs6_d40.py), as D48's and D49's readers
read their calibrations. B is minus the energy of the one electron, in eV; exact, in a complete basis, is
13.605693122994 q^2. A set is ADEQUATE at a size, with one treatment of the auxiliary basis, if the Hartree-Fock run
and the run with the functional are both measured and both hold within 0.05 eV of exact. THE RULE (R50.3) takes, for a
set, the first treatment in the order A, E, J that is ADEQUATE at every size from 0.25 e up.

PART G, the dianion. Every run is read as D48's reader read its runs with functions at the charge, by the committed
classifier (cs6_audit.py) and the committed reading of stability (cs6_d33b.py): the legs, the charge that left the 24
atoms, TOOK ELECTRONS or DID NOT TAKE, and ONTO THE FUNCTIONS or A LEG ONLY. A run with a charge is read against the
control that holds the same functions and no charge. That control is read against the plain committed control by the
legs that do not divide the density among atoms (the dipole, the gap, the stability): THE CONTROL HOLDS if none of them
fires and none is withheld. The Loewdin charge its functions are given is reported beside that word and is not part
of it.

The inputs are not taken on trust: every file of d50 is formed again by cs6_d50_prepare.py from the committed files
and compared, every output must print the committed input beside it, and the committed files this leans on must be
the committed ones. A gate that fails stops the reader before anything is written. D50 adopts nothing.

usage: python3 cs6_d50.py --res <results repository> --code <code repository> --prereg <the pre-registration>
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
PREREG = "STAGE6_DIAGNOSTIC50_20261122.txt"
D49PRE = "STAGE6_DIAGNOSTIC49_20261121.txt"
D49RES = "STAGE6_DIAGNOSTIC49_RESULT_20261122.txt"
D47PRE = "STAGE6_DIAGNOSTIC47_20261121.txt"
EV = 27.211386245988
RYD = EV / 2.0
H2K = 627.5094740631
TOL_EH = 1e-6          # hartree: a job's plain control against the committed one
TOL_HF = 0.005         # eV: a measured HF value against the one computed for the charge at the centre, before a note is made
FROM = 25              # hundredths of e: the rule of R50.3 asks for every size from this one up
NEL_DI = 118
# quoted sentences: each must be in the document it is quoted from and in the pre-registration given
Q_A = "B at or below zero: CAM at 0.25 e. These functions do not hold an electron at a charge of that size."      # D49's result
Q_B = "B_n11_q0.25: NOT A TEST OF CAPTURE or A WEAK TEST: its run with functions showed nothing, or little, beyond A9.4."   # D49, section 5
Q_C = "TOOK ELECTRONS the SCF converged and either a leg fires or the run LEFT THE ATOMS"                           # D47, as D48 quotes it
GUARDED = ["d31", "d33a", "d37", "d40", "d43", "d44", "d47", "d48", "d49"]
DOCS = ["d44.tsv", "d47.tsv", "d48.tsv", "d49.tsv", "d49_calibration.tsv", "d49_sprime.tsv", D49PRE, D49RES, D47PRE,
        "STAGE6_DIAGNOSTIC48_20261121.txt", "STAGE6_DIAGNOSTIC48_RESULT_20261121.txt"]


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


def counts(text):
    """(orbital basis functions, auxiliary functions for J) as ORCA prints them; None where it prints none or two values"""
    nb = sorted(set(int(v) for v in re.findall(r"Number of basis functions\s*\.+\s*(\d+)", text)))
    na = sorted(set(int(v) for v in re.findall(r"# of basis functions in Aux-J\s*\.+\s*(\d+)", text)))
    return (nb[0] if len(nb) == 1 else None), (na[0] if len(na) == 1 else None)


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    a.res, a.code, a.prereg = (os.path.abspath(x) for x in (a.res, a.code, a.prereg))
    S6 = os.path.join(a.res, S)
    D = os.path.join(S6, "d50")
    paths = dict((k, os.path.join(a.code, "phase2.2", f)) for k, f in (
        ("d40", "cs6_d40.py"), ("d33b", "cs6_d33b.py"), ("audit", "cs6_audit.py"), ("p48", "cs6_d48_prepare.py"), ("p50", "cs6_d50_prepare.py")))
    for k, p in paths.items():
        if not os.path.isfile(p):
            stop("the committed script %s is missing." % p)
    for f in [a.prereg, os.path.join(D, "RUNS.tsv")] + [os.path.join(S6, x) for x in DOCS]:
        if not os.path.isfile(f):
            stop("%s is missing." % f)
    tsv2, tsv3 = os.path.join(os.path.dirname(a.tsv), "d50_calibration.tsv"), os.path.join(os.path.dirname(a.tsv), "d50_dianion.tsv")
    for f in (a.result, a.tsv, tsv2, tsv3):
        if os.path.exists(f):
            stop("%s exists already; this reader overwrites nothing." % f)
    d40, d33b, audit = load(paths["d40"], "d50_cs6_d40"), load(paths["d33b"], "d50_cs6_d33b"), load(paths["audit"], "d50_cs6_audit")
    p50 = load(paths["p50"], "d50_cs6_d50_prepare")
    D_LW, R_DIP, GAP_ABS, R_GAP = audit.D_LW, audit.R_DIP, audit.GAP_ABS, audit.R_GAP
    gates = []

    def git(*args):
        try:
            q = subprocess.run(["git", "-C", a.res] + list(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            return q.returncode, q.stdout.strip()
        except OSError:
            return 127, ""

    # ------------------------------------------------------------ G50.3: the committed files this leans on
    guard = ["%s/%s" % (S, g) for g in GUARDED] + ["%s/%s" % (S, f) for f in DOCS]
    rc, st = git("status", "--porcelain", "--untracked-files=all", "--", *guard)
    if rc != 0:
        gates.append("G50.3 could not be checked: git did not run")
    elif st:
        gates.append("G50.3 a file under %s, or a committed table or document this reader names, is not the committed one" % ", ".join(GUARDED))
    rc, pcm = git("log", "--diff-filter=A", "--format=%H", "-1", "--", "%s/%s" % (S, PREREG))
    if rc != 0 or not pcm:
        gates.append("G50.3 the commit that added %s could not be found" % PREREG)
    else:
        rc1, b1 = git("rev-parse", "%s:%s/%s" % (pcm, S, PREREG))
        rc2, b2 = git("hash-object", a.prereg)
        if rc1 != 0 or rc2 != 0 or b1 != b2:
            gates.append("G50.3 the pre-registration given is not the committed %s" % PREREG)

    # ------------------------------------------------------------ G50.2 and G50.7: formed again from the committed files
    try:
        F, runs, sel, p48 = p50.build(S6, a.code)
    except Exception as e:
        sys.stdout.write("\n".join("FAILED  " + g for g in gates + ["G50.2 the files of d50 could not be formed again from the committed files: %s: %s" % (type(e).__name__, e)]) + "\n")
        stop("%d gate(s) failed. Nothing was written." % (len(gates) + 1))
    psrc = os.path.join(D, "SOURCES.sha256")
    if not os.path.isfile(psrc) or open(psrc).read() != F["SOURCES.sha256"]:
        gates.append("G50.2 d50/SOURCES.sha256 is not the list formed again from the committed files")
    nsrc = 0
    for l in F["SOURCES.sha256"].splitlines():
        h, f = l.split(None, 1)
        nsrc += 1
        p = os.path.join(S6, f)
        if not os.path.isfile(p) or sha(p) != h:
            gates.append("G50.2 %s is not the file d50/SOURCES.sha256 lists" % f)

    # ------------------------------------------------------------ G50.1: the inputs
    nsum = 0
    pic = os.path.join(D, "INPUT_CHECKSUMS.sha256")
    if not os.path.isfile(pic):
        gates.append("G50.1 d50/INPUT_CHECKSUMS.sha256 is missing")
    else:
        if open(pic).read() != F["INPUT_CHECKSUMS.sha256"]:
            gates.append("G50.1 d50/INPUT_CHECKSUMS.sha256 is not the list formed again here")
        for l in open(pic).read().splitlines():
            if l.strip():
                h, f = l.split(None, 1)
                nsum += 1
                p = os.path.join(D, f)
                if not os.path.isfile(p) or sha(p) != h:
                    gates.append("G50.1 %s is not the file listed in INPUT_CHECKSUMS.sha256" % f)
    for k in sorted(F):
        p = os.path.join(D, k)
        if not os.path.isfile(p) or open(p).read() != F[k]:
            gates.append("G50.1 d50/%s is not the file formed again from the committed files" % k)
    for r in runs:
        for fn in ("pc.pc", "substrate.xyz"):
            if os.path.exists(os.path.join(D, r["run"], fn)) and ("%s/%s" % (r["run"], fn)) not in F:
                gates.append("G50.1 %s holds a file %s that its input does not read" % (r["run"], fn))

    # ------------------------------------------------------------ G50.4: the quoted sentences
    def has(path, sent):
        return " ".join(sent.split()) in " ".join(open(path, errors="replace").read().split())
    for doc, q, what in ((D49RES, Q_A, "D49's finding at 0.25 e"), (D49PRE, Q_B, "D49's sentence on the test of B_n11_q0.25"),
                         (D47PRE, Q_C, "D47's reading TOOK ELECTRONS")):
        if not has(os.path.join(S6, doc), q):
            gates.append("G50.4 %s is not in %s" % (what, doc))
        if not has(a.prereg, q):
            gates.append("G50.4 the words \"%s\" are not in the pre-registration given" % (q[:60] + " ..."))

    # ------------------------------------------------------------ reading
    def read_di(d, kind, nel=NEL_DI):
        r = d40.read(d, kind, nel, d33b)
        r["status"], r["c"], r["lw"], r["nb"], r["naux"] = "measured", None, None, None, None
        if not r.get("exists"):
            r["status"] = "NOT MEASURED"
            return r
        r["nb"], r["naux"] = counts(r["text"])
        if r["E"] is None:
            r["status"] = "DID NOT CONVERGE" if "SCF NOT CONVERGED" in r["text"] else "NOT MEASURED"
            return r
        if nel != NEL_DI:
            return r
        c = audit.read(r["text"])
        lw = loewdin(r["text"])
        if c is None or c.get("maxlw") is None or c.get("dip") is None or c.get("nb") is None or len(lw) < 24:
            r["status"], r["E"] = "NOT MEASURED", None
            r["why"] = ["the committed classifier does not read a Loewdin charge of each atom, a dipole and a basis-function count from it"]
            return r
        r["c"], r["lw"] = c, lw
        return r

    # ---- G50.8: the committed outputs this reader leans on
    CT, G48, P0, C48 = {}, {}, {}, {}
    cnames = [x["name"] for x in sel["conds"]]
    plain_of = {p50.A_NAME: "d47/CAMTZ_%s_di_%%s" % p50.A_NAME, p50.B_NAME: "d44/CAMTZ_%s_di_%%s" % p50.B_NAME}
    for stt in p48.STATES:
        CT[stt] = read_di(os.path.join(S6, "%s_%s" % (p48.CTL, stt)), "plain")
        if CT[stt]["status"] != "measured" or not CT[stt]["stab"].get("ran") or CT[stt]["stab"].get("last") != "STABLE":
            gates.append("G50.8 the committed control %s_%s is not a measured dianion whose last stability analysis is STABLE" % (p48.CTL, stt))
        elif CT[stt]["nb"] != p50.NB_PLAIN:
            gates.append("G50.8 the committed control %s_%s reports %s basis functions, not %d" % (p48.CTL, stt, CT[stt]["nb"], p50.NB_PLAIN))
        for g in sel["gctl"]:
            G48[(g["name"], stt)] = read_di(os.path.join(S6, "d48", "CAMTZG_%s_di_%s" % (g["name"], stt)), "plain")
        for n in cnames:
            P0[(n, stt)] = read_di(os.path.join(S6, plain_of[n] % stt), "pc")
        C48[(p50.B_NAME, stt)] = read_di(os.path.join(S6, "d48", "CAMTZG_%s_di_%s" % (p50.B_NAME, stt)), "pc")
    for k, v in list(G48.items()) + list(P0.items()) + list(C48.items()):
        if v["status"] != "measured":
            gates.append("G50.8 a committed output this reader sets beside, of %s at %s, is not readable" % k)

    # ---- every run of D50, and G50.5, G50.6
    RUN = {}
    for r in runs:
        d = os.path.join(D, r["run"])
        x = read_di(d, "pc" if r["kind"] in ("one", "cond") else "plain", 1 if r["kind"] == "one" else NEL_DI)
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
                gates.append("G50.5 %s/job.out ends normally and does not print the input it ran" % r["run"])
        elif want is not None and e != want:
            gates.append("G50.5 the input printed at the top of %s/job.out is not the committed %s/job.inp" % (r["run"], r["run"]))
        if x["status"] == "measured" and x["nb"] != p50.nb_of(sel, r):
            gates.append("G50.6 %s reports %s basis functions, not %d" % (r["run"], x["nb"] if x["nb"] is not None else "no single count of", p50.nb_of(sel, r)))

    if gates:
        sys.stdout.write("\n".join("FAILED  " + g for g in gates) + "\n")
        stop("%d gate(s) failed. Nothing was written." % len(gates))

    # ------------------------------------------------------------ PART C
    sets, sizes = sorted(p50.ADD), list(p50.SIZES)
    CAL = {}
    for name in sets:
        for k in sizes:
            q = k / 100.0
            hf = RUN[p50.one_name("HF", name, None, k)]
            c = dict(set=name, k=k, q=q, exact=RYD * q * q, hf=-hf["E"] * EV if hf["status"] == "measured" else None, exp=-p50.EXPECTED[name][k] * EV, cam={}, word={})
            for t in p50.AUX:
                x = RUN[p50.one_name("CAM", name, t, k)]
                c["cam"][t] = -x["E"] * EV if x["status"] == "measured" else None
                if c["hf"] is None or c["cam"][t] is None:
                    c["word"][t] = "NOT MEASURED"
                else:
                    c["word"][t] = "ADEQUATE" if (abs(c["exact"] - c["hf"]) <= p50.TOL_ADEQ and abs(c["exact"] - c["cam"][t]) <= p50.TOL_ADEQ) else "NOT ADEQUATE"
            CAL[(name, k)] = c
    RULE, DOWN = {}, {}
    for name in sets:
        RULE[name] = None
        for t in p50.AUX:
            ok = [k for k in sizes if CAL[(name, k)]["word"][t] == "ADEQUATE"]
            low = None
            for k in reversed(sizes):                       # the smallest size from which every size up is ADEQUATE
                if k in ok:
                    low = k
                else:
                    break
            DOWN[(name, t)] = low
            if RULE[name] is None and all(k in ok for k in sizes if k >= FROM):
                RULE[name] = t
    chosen = next(((name, RULE[name]) for name in sets if RULE[name] is not None), None)

    # ------------------------------------------------------------ PART G
    G = sel["G"]
    cm = dict((x["name"], x) for x in sel["conds"])
    gm = dict((g["name"], g) for g in sel["gctl"])

    def moved(lw, lw0, state, charges):
        """the atom with the most negative change, that change, its distance to the nearest positive charge, the sum over the 24"""
        dq = [lw[i] - lw0[i] for i in range(24)]
        i = min(range(24), key=lambda j: dq[j])
        pos = [xyz for q, xyz in charges if q > 0]
        dn = min(p48.dist(G[state][i][1:], p) for p in pos) if pos else None
        return p48.ATOMS[i], dq[i], dn, sum(lw[:24]) - sum(lw0[:24])

    def legs(x, c0, charges):
        """the flags and holds of one dianion run against its control (D48's reader, line for line)"""
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
        z["mv"] = moved(x["lw"], c0["lw"], x["state"], charges)
        return z

    NOZ = dict(flags=[], held=[], rise=None, ratio=None, gap=None, gap_read=False, mv=None)
    JUSE = {}
    for jk, pbs in p50.JOBS:
        u = None
        if jk != "K1":
            c = RUN["%s_C_ctl_R" % jk]
            if c["status"] != "measured":
                u = "its control %s_C_ctl_R is %s" % (jk, c["status"])
            elif not c["stab"].get("ran") or c["stab"].get("last") != "STABLE":
                u = "the last stability analysis of its control %s_C_ctl_R is not STABLE" % jk
            elif abs(c["E"] - CT["R"]["E"]) > TOL_EH:
                u = "its control %s_C_ctl_R does not reproduce the committed energy (differ %.2e Eh)" % (jk, abs(c["E"] - CT["R"]["E"]))
        JUSE[jk] = u

    def word1(x):
        if x["noread"]:
            return "NOT READ"
        if x["z"]["flags"]:
            return "FLAGGED"
        if x["status"] != "measured":
            return "NOT MEASURED"
        return "HELD" if x["z"]["held"] else "NOT FLAGGED"

    def took(x):
        mv = x["z"]["mv"]
        if x["noread"] or x["status"] != "measured" or mv is None:
            return "not decided", None
        if mv[3] >= D_LW or x["z"]["flags"]:
            return "TOOK ELECTRONS", mv[3]
        return ("DID NOT TAKE" if x["word"] == "NOT FLAGGED" else "not decided"), mv[3]

    for r in runs:
        x = RUN[r["run"]]
        if r["kind"] == "gctl":
            c0 = CT[r["state"]]
            x["noread"] = JUSE[r["job"]]
            x["z"] = legs(x, c0, []) if not x["noread"] else dict(NOZ)
            x["word"] = word1(x)
            mv = x["z"]["mv"]
            x["left"] = mv[3] if mv else None
            hard = [f for f in x["z"]["flags"] if not f.startswith("relocation")]
            if x["noread"] or (x["status"] != "measured" and not hard):
                x["holds"] = "not decided"
            elif hard:
                x["holds"] = "DOES NOT HOLD"
            else:
                x["holds"] = "not decided" if x["z"]["held"] else "HOLDS"
            x["share"] = bool(x["holds"] == "HOLDS" and (x["left"] >= D_LW or len(hard) != len(x["z"]["flags"])))
            ng = len(gm[r["gname"]]["ghosts"])
            x["ghost"] = x["lw"][24:24 + ng] if x["status"] == "measured" and x["lw"] and len(x["lw"]) >= 24 + ng else None
            o = G48[(r["gname"], r["state"])]
            x["dE_plain"] = (x["E"] - c0["E"]) * EV if x["status"] == "measured" else None
            x["dE_48"] = (x["E"] - o["E"]) * EV if x["status"] == "measured" else None
    for r in runs:
        x = RUN[r["run"]]
        if r["kind"] != "cond":
            continue
        cd = cm[r["cond"]]
        c0 = RUN[r["ref"]]
        why = None
        if c0["status"] != "measured":
            why = "its control %s is %s" % (r["ref"], c0["status"])
        elif not c0["stab"].get("ran") or c0["stab"].get("last") is None:
            why = "the stability analysis of its control %s was not read" % r["ref"]
        elif c0["stab"]["last"] != "STABLE":
            why = "its control %s is UNSTABLE at the last analysis" % r["ref"]
        x["noread"] = JUSE[r["job"]] or why
        x["z"] = legs(x, c0, cd["charges"]) if not x["noread"] else dict(NOZ)
        x["word"] = word1(x)
        x["took"] = took(x)
        x["how"] = ("ONTO THE FUNCTIONS" if x["took"][1] >= D_LW else "A LEG ONLY") if x["took"][0] == "TOOK ELECTRONS" else ""
        x["S"], x["pred"] = cd["S"][r["state"]], cd["pred"][r["state"]]
        did = x["took"][0]
        x["verdict"] = "none made" if x["pred"] == "none" else ("not decided" if did == "not decided" else ("borne out" if did == x["pred"] else "NOT BORNE OUT"))
        ng = len(cd["ghosts"])
        x["ghost"] = x["lw"][24:24 + ng] if x["status"] == "measured" and x["lw"] and len(x["lw"]) >= 24 + ng else None
        x["gmost"] = None
        if x["ghost"] is not None and not x["noread"] and c0.get("lw") and len(c0["lw"]) >= 24 + ng:
            dg = [x["lw"][24 + k] - c0["lw"][24 + k] for k in range(ng)]
            k = min(range(ng), key=lambda j: dg[j])
            x["gmost"] = ([i for q, i in cd["pairs"] if q > 0][k], dg[k])
        p0 = P0[(r["cond"], r["state"])]
        ok = x["status"] == "measured" and c0["status"] == "measured"
        x["gain"] = ((x["E"] - c0["E"]) - (p0["E"] - CT[r["state"]]["E"])) * EV if ok else None
        o = C48.get((r["cond"], r["state"]))
        x["gain48"] = ((o["E"] - G48[(cd["gctl"], r["state"])]["E"]) - (p0["E"] - CT[r["state"]]["E"])) * EV if o else None

    def dde(cond, aux):
        """(ddE with the functions, ddE of the committed runs without them), kcal/mol; None where a run is not measured"""
        x = dict((stt, RUN[p50.xname("cond", cond, aux, stt)]) for stt in p48.STATES)
        c = dict((stt, RUN[p50.xname("gctl", cm[cond]["gctl"], aux, stt)]) for stt in p48.STATES)
        plain = ((P0[(cond, "TS")]["E"] - P0[(cond, "R")]["E"]) - (CT["TS"]["E"] - CT["R"]["E"])) * H2K
        if any(v["status"] != "measured" or v.get("noread") for v in list(x.values()) + list(c.values())):
            return None, plain
        return ((x["TS"]["E"] - x["R"]["E"]) - (c["TS"]["E"] - c["R"]["E"])) * H2K, plain

    # ------------------------------------------------------------ the result
    rep = []
    say = rep.append

    def fs(v, fmt="%.6f", none="-"):
        return none if v is None else fmt % v

    def wrapsay(text, lead="   ", width=96, hang=""):
        for i, l in enumerate(textwrap.wrap(text, width - len(lead))):
            say(lead + (hang if i else "") + l)
    rc, head = git("rev-parse", "HEAD")
    say("STAGE 6 - DIAGNOSTIC 50 - RESULT")
    say("Functions at a charge that can hold an electron at 0.25 e")
    say("=" * 96)
    say("")
    say("Read by cs6_d50.py under %s, which was committed before any job was submitted." % PREREG)
    say("   results repository at %s; the pre-registration was added by %s" % (head[:12] if rc == 0 else "?", pcm[:12]))
    say("   reader               cs6_d50.py           sha256 %s" % sha(os.path.abspath(__file__)))
    say("   inputs formed again  cs6_d50_prepare.py   sha256 %s" % sha(paths["p50"]))
    say("   conditions and S     cs6_d48_prepare.py   sha256 %s" % sha(paths["p48"]))
    say("   outputs read by      cs6_d40.py           sha256 %s" % sha(paths["d40"]))
    say("   thresholds and legs  cs6_audit.py         sha256 %s" % sha(paths["audit"]))
    say("   stability by         cs6_d33b.py          sha256 %s" % sha(paths["d33b"]))
    say("")
    wrapsay("Part C: one electron, a bare charge q, the def2-TZVPD functions of hydrogen with diffuse functions added. Part G: the "
            "dianion of frame 42436 at CAM-B3LYP-D3BJ/def2-TZVPD with set %s at each positive charge. Legs at D_LW %.2f e, R_DIP %.1f, "
            "GAP_ABS %.1f eV, R_GAP %.1f, stability on the last analysis. D50 adopts nothing: it opens and closes no size, no "
            "distance and no design." % (p50.G_SET, D_LW, R_DIP, GAP_ABS, R_GAP), lead="")
    say("")
    say("THE GATES, as they held:")
    say("   G50.1  all %d files of d50/INPUT_CHECKSUMS.sha256 are the ones listed, and every file of d50 is the" % nsum)
    say("          one formed again here from the committed files.")
    say("   G50.2  the %d files of d50/SOURCES.sha256 are the ones listed." % nsrc)
    say("   G50.3  no file under %s," % ", ".join(GUARDED))
    say("          and no committed table or document named here, differs from the committed one; the")
    say("          pre-registration given is the committed one.")
    say("   G50.4  the sentences quoted here are in the documents they are quoted from.")
    say("   G50.5  every output that exists prints, as the input it ran, the committed input beside it.")
    say("   G50.6  every measured run reports the number of basis functions its input asks for.")
    say("   G50.7  S of the two conditions, formed again here, is the S that d47.tsv and d48.tsv hold.")
    say("   G50.8  the committed controls and the committed outputs set beside the runs are readable.")
    say("")

    def title(t):
        say("=" * 96)
        say(t)
        say("=" * 96)
        say("")

    # ---- R50.1
    title("R50.1  THE RUNS OF ONE ELECTRON")
    say("   %-20s %-4s %-5s %-5s %-17s %-15s %-6s %s" % ("run", "set", "level", "q, e", "status", "energy, Eh", "basis", "auxiliary"))
    ones = [r for r in runs if r["kind"] == "one"]
    for r in ones:
        x = RUN[r["run"]]
        say("   %-20s %-4s %-5s %.2f  %-17s %-15s %-6s %s" % (r["run"], r["set"], r["level"].replace("ONE_", "") + (r["aux"] or ""), r["size"] / 100.0, x["status"],
                                                          fs(x["E"], "%.10f"), fs(x["nb"], "%d"), fs(x["naux"], "%d")))
    say("")
    nm = sum(1 for r in ones if RUN[r["run"]]["status"] == "measured")
    say("   Measured: %d of %d." % (nm, len(ones)))
    for lab, pick in (("HF", lambda r: r["level"] == "ONE_HF"),) + tuple(("CAM with treatment %s" % t, (lambda r, t=t: r["level"] == "ONE_CAM" and r["aux"] == t)) for t in p50.AUX):
        sub = [r for r in ones if pick(r)]
        say("      %-22s %d of %d" % (lab + ":", sum(1 for r in sub if RUN[r["run"]]["status"] == "measured"), len(sub)))
    notes, why = [], {}
    for r in ones:
        x = RUN[r["run"]]
        if x["status"] != "measured":
            why.setdefault("%s%s" % (x["status"], (": " + "; ".join(x.get("why") or [])) if x.get("why") else ""), []).append(r["run"])
        elif r["level"] == "ONE_CAM":
            want = p50.naux(r["set"], r["aux"])
            if (want is not None and x["naux"] != want) or (want is None and (x["naux"] is None or x["naux"] <= p50.NAUX_J)):
                notes.append("%s: reports %s auxiliary functions; %s" % (r["run"], fs(x["naux"], "%d", "no count of"),
                                                                          ("%d were expected" % want) if want is not None else "more than %d were expected of AutoAux" % p50.NAUX_J))
    for w, rs in why.items():                               # runs that failed in one way are named together
        notes.insert(0, "%d run(s) %s. They are: %s" % (len(rs), w, ", ".join(rs)))
    for n in notes:
        wrapsay(n, hang="   ")
    if notes:
        say("")

    # ---- R50.2
    title("R50.2  WHAT THE SETS HOLD AT A CHARGE q")
    wrapsay("B is minus the energy of the one electron, eV. Exact, in a complete basis: 13.6057 q^2. 'computed' is the Hartree-Fock value "
            "Claude computed for the pre-registration, with the charge at the centre of the functions; the runs have it 0.01 A away. "
            "A, E, J: the functional with the three treatments of the ghost atom's auxiliary basis. A word under each: ADEQUATE if "
            "the HF run and that run both hold within %.2f eV of exact." % p50.TOL_ADEQ, lead="")
    say("")
    for name in sets:
        say("SET %s: def2-TZVPD hydrogen and an s and a p function at each of %s (%d functions)." % (
            name, ", ".join("%g" % e for e in p50.ADD[name]), p50.NBF[name]))
        say("   q, e   exact     B, HF      computed   B, CAM A   B, CAM E   B, CAM J   exact-HF   A            E            J")
        for k in sizes:
            c = CAL[(name, k)]
            say("   %.2f   %-9s %-10s %-10s %-10s %-10s %-10s %-10s %-12s %-12s %s" % (
                c["q"], "%.4f" % c["exact"], fs(c["hf"], "%+.4f"), "%+.4f" % c["exp"], fs(c["cam"]["A"], "%+.4f"), fs(c["cam"]["E"], "%+.4f"), fs(c["cam"]["J"], "%+.4f"),
                fs(None if c["hf"] is None else c["exact"] - c["hf"], "%+.4f"), c["word"]["A"], c["word"]["E"], c["word"]["J"]))
        far = [(c["q"], c["hf"] - c["exp"]) for c in (CAL[(name, k)] for k in sizes) if c["hf"] is not None and abs(c["hf"] - c["exp"]) > TOL_HF]
        got = [abs(c["hf"] - c["exp"]) for c in (CAL[(name, k)] for k in sizes) if c["hf"] is not None]
        if got:
            say("   HF measured against computed: formed at %d of %d sizes; the largest difference is %.5f eV." % (len(got), len(sizes), max(got)))
        else:
            say("   HF measured against computed: not formed, no HF run of this set is measured.")
        if far:
            wrapsay("MORE THAN %.3f eV from the computed value at %d size(s): %s. The functions ORCA used are then not the ones "
                    "the computed values are for, or the 0.01 A matters more than D49 found." % (TOL_HF, len(far), ", ".join("%.2f e (%+.4f)" % f for f in far)), hang="   ")
        for t in p50.AUX:
            d = [(c["q"], c["cam"][t] - c["hf"]) for c in (CAL[(name, k)] for k in sizes) if c["hf"] is not None and c["cam"][t] is not None]
            if d:
                big = max(d, key=lambda v: abs(v[1]))
                say("   CAM %s minus HF: formed at %d of %d sizes; the largest is %+.4f eV (%.2f e)." % (t, len(d), len(sizes), big[1], big[0]))
            else:
                say("   CAM %s minus HF: not formed." % t)
        say("")

    # ---- R50.3
    title("R50.3  THE FUNCTIONS FOR THE CAPTURE TEST, BY THE RULE")
    wrapsay("The rule, fixed before the runs: for a set, the first treatment in the order A, E, J that is ADEQUATE at every size from "
            "%.2f e up. The set is %s if the rule gives it a treatment, otherwise %s." % (FROM / 100.0, sets[0], " then ".join(sets[1:])), lead="")
    say("")
    for name in sets:
        say("   %s: %s" % (name, ("treatment %s" % RULE[name]) if RULE[name] else "NO TREATMENT is adequate at every size from %.2f e up" % (FROM / 100.0)))
        for t in p50.AUX:
            low = DOWN[(name, t)]
            bad = [k for k in sizes if k >= FROM and CAL[(name, k)]["word"][t] != "ADEQUATE"]
            say("      %s  adequate without a break from %s up%s" % (t, ("%.2f e" % (low / 100.0)) if low is not None else "no size",
                                                                 ("; not at %s" % ", ".join("%.2f e (%s)" % (k / 100.0, CAL[(name, k)]["word"][t]) for k in bad)) if bad else ""))
    say("")
    if chosen:
        say("   THE FUNCTIONS FOR THE CAPTURE TEST: set %s with treatment %s." % chosen)
    else:
        say("   THE FUNCTIONS FOR THE CAPTURE TEST: NONE. No set has a treatment that is adequate at every size from %.2f e up." % (FROM / 100.0))
    say("")

    # ---- R50.4
    title("R50.4  THE DIANION WITH THE FUNCTIONS AND NO CHARGE")
    wrapsay("Each control holds set %s at the ghost centres of its condition and no charge, and is read against the plain committed "
            "control of its state (%s). THE CONTROL HOLDS if the dipole, the gap and the stability legs do not fire and none is "
            "withheld; DOES NOT HOLD if one fires. 'left' is the electron charge missing from the 24 atoms' Loewdin charges against "
            "that control and 'on the functions' the sum of the Loewdin charges ORCA gives the ghost atoms: a way of dividing the "
            "density, reported beside the word and not part of it. dE: the energy against the plain control, and against D48's "
            "control with the def2-TZVPD functions of hydrogen at the same centres, eV." % (p50.G_SET, p48.CTL), lead="")
    say("")
    for jk, pbs in p50.JOBS:
        if jk == "K1":
            continue
        c = RUN["%s_C_ctl_R" % jk]
        say("   job %s: its plain control %s" % (jk, ("reproduces the committed energy (%.2e Eh apart) and is STABLE" % abs(c["E"] - CT["R"]["E"])) if not JUSE[jk] else "FAILS: %s. Its runs are NOT READ." % JUSE[jk]))
    say("")
    say("   %-30s %-2s %-17s %-6s %-6s %-9s %-10s %-8s %-8s %-9s %-9s %-13s %s" % (
        "run", "st", "status", "basis", "aux", "left, e", "on functs", "gap, eV", "HOMO, eV", "dE plain", "dE D48", "stability", "the control"))
    gruns = [r for r in runs if r["kind"] == "gctl"]
    tail = []
    for r in gruns:
        x = RUN[r["run"]]
        sb = x.get("stab") or {}
        say("   %-30s %-2s %-17s %-6s %-6s %-9s %-10s %-8s %-8s %-9s %-9s %-13s %s" % (
            r["run"], r["state"], x["status"], fs(x["nb"], "%d"), fs(x["naux"], "%d"), fs(x["left"], "%+.4f"), fs(sum(x["ghost"]) if x["ghost"] else None, "%+.4f"),
            fs(x["z"]["gap"], "%.4f"), fs(x["c"]["homo"] if x["c"] else None, "%+.3f"), fs(x["dE_plain"], "%+.4f"), fs(x["dE_48"], "%+.4f"),
            "NOT FOUND" if not sb.get("ran") else ("NOT READ" if sb.get("first") is None else "%s/%s" % (sb["first"], sb["last"])), x["holds"]))
        for w in x["z"]["flags"]:
            tail.append("%s: %s" % (r["run"], w))
        for w in x["z"]["held"]:
            tail.append("%s: held: %s" % (r["run"], w))
        if x["noread"]:
            tail.append("%s: NOT READ: %s" % (r["run"], x["noread"]))
        elif x["status"] != "measured":
            tail.append("%s: %s%s" % (r["run"], x["status"], (": " + "; ".join(x.get("why") or [])) if x.get("why") else ""))
    say("")
    for n in tail:
        wrapsay(n, hang="   ")
    if tail:
        say("")
    nh = sum(1 for r in gruns if RUN[r["run"]]["holds"] == "HOLDS")
    nd = sum(1 for r in gruns if RUN[r["run"]]["holds"] == "DOES NOT HOLD")
    say("   The control HOLDS in %d of %d runs, DOES NOT HOLD in %d, is not decided in %d." % (nh, len(gruns), nd, len(gruns) - nh - nd))
    shr = [r["run"] for r in gruns if RUN[r["run"]].get("share")]
    if shr:
        wrapsay("THE LOEWDIN CHARGES GIVE THE FUNCTIONS %.2f e OR MORE while the control holds, in %d run(s): %s. No leg shows the "
                "electrons moved; the division of the density has." % (D_LW, len(shr), ", ".join(shr)), hang="   ")
    else:
        say("   In no control that holds do the Loewdin charges give the functions %.2f e or more." % D_LW)
    say("")

    # ---- R50.5
    title("R50.5  THE DIANION WITH THE FUNCTIONS AND THE CHARGES")
    wrapsay("Each run is read against the control that holds the same functions and no charge, as D48 read its runs. S is D47's "
            "estimate, formed from the committed maps; it uses 13.6057 q^2 and predicts DID NOT TAKE at +0.5 eV or more. 'gain' is "
            "the energy the functions add when the charges are there beyond what they add without them, eV: (run minus its "
            "control) minus (the committed run of the same charges without functions minus the plain control).", lead="")
    say("")
    say("   %-34s %-2s %-17s %-9s %-9s %-8s %-11s %-14s %-18s %-7s %-13s %s" % (
        "run", "st", "status", "left, e", "gain, eV", "gap, eV", "four legs", "the run", "how", "S, eV", "S predicts", "prediction"))
    cruns = [r for r in runs if r["kind"] == "cond"]
    tail = []
    for r in cruns:
        x = RUN[r["run"]]
        say("   %-34s %-2s %-17s %-9s %-9s %-8s %-11s %-14s %-18s %-7s %-13s %s" % (
            r["run"], r["state"], x["status"], fs(x["took"][1], "%+.4f"), fs(x["gain"], "%+.4f"), fs(x["z"]["gap"], "%.4f"), x["word"], x["took"][0], x["how"] or "-",
            "%+.2f" % x["S"], x["pred"], x["verdict"]))
        for w in x["z"]["flags"]:
            tail.append("%s: %s" % (r["run"], w))
        for w in x["z"]["held"]:
            tail.append("%s: held: %s" % (r["run"], w))
        if x["noread"]:
            tail.append("%s: NOT READ: %s" % (r["run"], x["noread"]))
        elif x["status"] != "measured":
            tail.append("%s: %s%s" % (r["run"], x["status"], (": " + "; ".join(x.get("why") or [])) if x.get("why") else ""))
        if x["took"][0] == "TOOK ELECTRONS" and x["gmost"]:
            tail.append("%s: the ghost centre that gained most is at site %d: %+.4f e against its control" % (r["run"], x["gmost"][0], x["gmost"][1]))
        if RUN[r["ref"]].get("holds") == "DOES NOT HOLD":
            tail.append("%s: its control %s DOES NOT HOLD; the run is read against it all the same, by D48's rule" % (r["run"], r["ref"]))
    say("")
    for n in tail:
        wrapsay(n, hang="   ")
    if tail:
        say("")
    for n in cnames:
        say("   %s: %s" % (n, cm[n]["where"]))
        for t in p50.G_AUX:
            w, plain = dde(n, t)
            say("      ddE with the functions, treatment %s: %s; without functions (committed runs): %+.4f kcal/mol; difference %s" % (
                t, fs(w, "%+.4f kcal/mol", "not formed"), plain, fs(None if w is None else w - plain, "%+.4f", "not formed")))
        if (n, "R") in C48:
            say("      gain with D48's def2-TZVPD functions at the same centres (committed runs): R %+.4f eV, TS %+.4f eV" % (
                RUN[p50.xname("cond", n, p50.G_AUX[0], "R")]["gain48"], RUN[p50.xname("cond", n, p50.G_AUX[0], "TS")]["gain48"]))
    say("")
    did = [RUN[r["run"]]["took"][0] for r in cruns]
    ver = [RUN[r["run"]]["verdict"] for r in cruns]
    say("   Runs with charges: %d; took electrons %d, did not take %d, not decided %d." % (len(cruns), did.count("TOOK ELECTRONS"), did.count("DID NOT TAKE"), did.count("not decided")))
    say("   Predictions of S: borne out %d, NOT borne out %d, none made %d, not decided %d." % (ver.count("borne out"), ver.count("NOT BORNE OUT"), ver.count("none made"), ver.count("not decided")))
    for r in cruns:
        x = RUN[r["run"]]
        if x["verdict"] == "NOT BORNE OUT":
            say("      NOT BORNE OUT: %s: S %+.2f eV predicted %s; the run %s%s." % (r["run"], x["S"], x["pred"], x["took"][0], ", %s" % x["how"] if x["how"] else ""))
    say("")

    # ---- R50.6
    REF = p50.G_AUX[-1]
    others = [u for u in p50.G_AUX if u != REF]
    title("R50.6  THE TREATMENTS OF THE AUXILIARY BASIS, ON THE DIANION")
    wrapsay("The same run with each treatment: its energy minus that of treatment %s, eV, and what each run did." % REF, lead="")
    say("")
    groups = []
    for r in runs:
        if r["kind"] in ("gctl", "cond") and r["aux"] == p50.G_AUX[0]:
            xs = dict((u, RUN[p50.xname(r["kind"], r["cond"], u, r["state"])]) for u in p50.G_AUX)
            wd = dict((u, xs[u]["holds"] if r["kind"] == "gctl" else xs[u]["took"][0]) for u in p50.G_AUX)
            dE = dict((u, (xs[u]["E"] - xs[REF]["E"]) * EV if (xs[u]["status"] == "measured" and xs[REF]["status"] == "measured") else None) for u in others)
            groups.append((r, wd, dE))
            say("   %-24s %-2s  %s   %s" % (("control " if r["kind"] == "gctl" else "") + r["cond"], r["state"],
                                           "  ".join("%s minus %s: %-10s" % (u, REF, fs(dE[u], "%+.5f", "not formed")) for u in others),
                                           "  ".join("%s: %-14s" % (u, wd[u]) for u in p50.G_AUX).rstrip()))
    dd = [abs(v) for r, wd, dE in groups for v in dE.values() if v is not None]
    full = [g for g in groups if all(w != "not decided" for w in g[1].values())]
    diff = [g[0] for g in full if len(set(g[1].values())) > 1]
    say("")
    say("   Every treatment is decided in %d of %d runs%s." % (len(full), len(groups), ("; the largest energy difference from %s is %.5f eV" % (REF, max(dd))) if dd else ""))
    if not full:
        say("   No run is decided with every treatment, so the treatments are not compared.")
    elif not diff:
        say("   Where every treatment is decided, they give the same word.")
    else:
        say("   THE TREATMENTS GIVE DIFFERENT WORDS in %d run(s): %s." % (len(diff), ", ".join("%s%s at %s" % ("control " if r["kind"] == "gctl" else "", r["cond"], r["state"]) for r in diff)))
    say("")

    title("WHAT THIS CANNOT SAY (section 6 of the pre-registration, in short; its own words govern)")
    for t in ("One electron and a bare charge say what the functions hold. They do not say what a run of the dianion does.",
              "Part G is one site and one design, at R and TS, one frame, frozen geometries, in vacuo, bare point charges. It is "
              "not the test of the inner sites, of the chemical groups or of charges placed together; that is the next diagnostic.",
              "A run that DID NOT TAKE with these functions is a test of capture at 0.25 e only if R50.3 finds the set adequate with "
              "the treatment that run used.",
              "S is an estimate. D48 found that S does not separate what its runs did, and D50 adopts nothing: no size, no distance "
              "and no design is opened or closed here.",
              "The Loewdin charge of a diffuse function is a way of dividing the density, not a count of electrons that moved; the "
              "control with the same functions is what each run is read against."):
        wrapsay(t, hang="")
    say("")
    tookn = [r["run"] for r in cruns if RUN[r["run"]]["took"][0] == "TOOK ELECTRONS"]
    last = ("One electron measured: %d of %d. The functions for the capture test, by the rule: %s. The control with the functions and no charge: holds in %d of %d runs, "
            "does not hold in %d%s. Runs with charges: took electrons %d of %d%s, did not take %d, not decided %d. Predictions of S: borne out %d, not borne out %d, "
            "not decided %d." % (nm, len(ones), ("set %s with treatment %s" % chosen) if chosen else "NONE", nh, len(gruns), nd,
                                 (" (the Loewdin charges give the functions %.2f e or more in %d that hold)" % (D_LW, len(shr))) if shr else "", len(tookn), len(cruns), (" (%s)" % ", ".join(tookn)) if tookn else "", did.count("DID NOT TAKE"), did.count("not decided"),
                                 ver.count("borne out"), ver.count("NOT BORNE OUT"), ver.count("not decided")))
    say(last)

    # ------------------------------------------------------------ the tables
    T1 = ["run\tjob\tlevel\tkind\tcondition\tstate\tset\taux\tsize_e\tstatus\tE_Eh\tbasis_functions\tauxiliary_functions"]
    for r in runs:
        x = RUN[r["run"]]
        T1.append("\t".join([r["run"], r["job"], r["level"], r["kind"], r["cond"], r["state"], r["set"], r["aux"], "" if r["size"] is None else "%.2f" % (r["size"] / 100.0),
                             x["status"], fs(x["E"], "%.10f", ""), fs(x["nb"], "%d", ""), fs(x["naux"], "%d", "")]))
    T2 = ["set\tsize_e\texact_eV\tB_HF_eV\tB_HF_computed_eV\tB_CAM_A_eV\tB_CAM_E_eV\tB_CAM_J_eV\tword_A\tword_E\tword_J\ttreatment_by_the_rule"]
    for name in sets:
        for k in sizes:
            c = CAL[(name, k)]
            T2.append("\t".join([name, "%.2f" % c["q"], "%.6f" % c["exact"], fs(c["hf"], "%.6f", ""), "%.6f" % c["exp"], fs(c["cam"]["A"], "%.6f", ""), fs(c["cam"]["E"], "%.6f", ""),
                                 fs(c["cam"]["J"], "%.6f", ""), c["word"]["A"], c["word"]["E"], c["word"]["J"], RULE[name] or ""]))
    T3 = ["run\tjob\tkind\tcondition\taux\tstate\tstatus\tE_Eh\tcontrol\tE_control_Eh\tbasis_functions\tauxiliary_functions\tlowdin_rise_e\tdipole_ratio\tgap_eV\tgap_read\tHOMO_eV\t"
          "stab_first\tstab_last\tS2\tTRAH\tlevel_shift\tflags\theld\tmost_negative_atom\tmost_negative_change_e\tsum_24_minus_control_e\tghost_rows_e\tword\tthe_control\t"
          "what_the_run_did\thow_it_took\tS_eV\tS_predicts\tprediction\tgain_eV\tdE_vs_plain_control_eV\tdE_vs_D48_control_eV\tghost_most_site\tghost_most_change_e"]
    for r in runs:
        if r["kind"] == "one":
            continue
        x = RUN[r["run"]]
        if r["kind"] == "ctl":
            c0n, c0e = "%s_R (committed)" % p48.CTL, CT["R"]["E"]
            z, word = dict(NOZ), ("NOT READ" if JUSE[r["job"]] else "reproduces the committed control")
        elif r["kind"] == "gctl":
            c0n, c0e, z, word = "%s_%s (committed)" % (p48.CTL, r["state"]), CT[r["state"]]["E"], x["z"], x["word"]
        else:
            c0n, c0e, z, word = r["ref"], RUN[r["ref"]].get("E"), x["z"], x["word"]
        sb = x.get("stab") or {}
        mv = z["mv"]
        T3.append("\t".join([
            r["run"], r["job"], r["kind"], r["cond"], r["aux"], r["state"], x["status"], fs(x["E"], "%.8f", ""), c0n, fs(c0e, "%.8f", ""), fs(x["nb"], "%d", ""), fs(x["naux"], "%d", ""),
            fs(z["rise"], "%+.4f", ""), fs(z["ratio"], "%.4f", ""), fs(z["gap"], "%.4f", ""), "yes" if z["gap_read"] else "no", fs(x["c"]["homo"] if x.get("c") else None, "%+.4f", ""),
            str(sb.get("first") or ""), str(sb.get("last") or ""), fs(x.get("s2"), "%.6f", ""), "yes" if x.get("trah") else "no", "yes" if x.get("lshift") else "no",
            "; ".join(z["flags"]), "; ".join(z["held"] + ([x["noread"]] if x.get("noread") else [])), mv[0] if mv else "", fs(mv and mv[1], "%+.4f", ""), fs(mv and mv[3], "%+.4f", ""),
            " ".join("%+.4f" % v for v in (x.get("ghost") or [])), word, x.get("holds", ""), x["took"][0] if r["kind"] == "cond" else "", x.get("how", ""),
            ("%+.4f" % x["S"]) if r["kind"] == "cond" else "", x.get("pred", ""), x.get("verdict", ""), fs(x.get("gain"), "%+.6f", ""), fs(x.get("dE_plain"), "%+.6f", ""),
            fs(x.get("dE_48"), "%+.6f", ""), str(x["gmost"][0]) if x.get("gmost") else "", fs(x["gmost"][1] if x.get("gmost") else None, "%+.4f", "")]))
    for p, t in ((a.result, rep), (a.tsv, T1), (tsv2, T2), (tsv3, T3)):
        with open(p, "w") as f:
            f.write("\n".join(t) + "\n")
    sys.stdout.write("\n".join(rep) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
