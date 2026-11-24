#!/usr/bin/env python3
"""Stage 6, DIAGNOSTIC 56: the reader. The designs of the search, each read at the working level and at
CAM-B3LYP-D3BJ/def2-TZVPD by the four legs, the fifth leg (AMENDMENT 6) and the tolerance of D44, and named by A9.4 (i).

Committed before any job is submitted. It reads each output with the committed functions of cs6_d40.py, cs6_d33b.py,
cs6_audit.py and cs6_d53.py. Before it reads any output of D56 it applies the same functions to the committed outputs of
D43 and D44 and stops unless they give the committed words and numbers of d43.tsv and d44.tsv. It overwrites nothing and
stops, writing nothing, if a gate fails.

usage: cs6_d56.py --res <results repo> --code <code repo> --prereg <pre-registration> --result <file> --tsv <file>
"""
import argparse, csv, os, subprocess, sys, textwrap

S = "phase2.2/cs_stage6"
PREREG = "STAGE6_DIAGNOSTIC56_20261124.txt"
AMEND11 = "STAGE6_AMENDMENT11_20261124.txt"
AMEND9 = "STAGE6_AMENDMENT9_20261120.txt"
AMEND6 = "STAGE6_AMENDMENT6_20261117.txt"
Q_J = "okay yes please build that, but ensure you incorperate those tests correctly"                  # Jessica, 24 November 2026
Q_A9 = "only if it passes the five legs at the working level AND, at CAM-B3LYP-D3BJ/def2-TZVPD, is RELIABLE on four legs, passes the fifth leg, and has a ddE that is WITHIN the tolerance of its ddE at the working level."
Q_A6 = "A design PASSES THE FIFTH LEG if it is BOUND at both states."
GUARDED = ["d31", "d40", "d43", "d44", "d53"]
DOCS = ["sites_f42436.tsv", "A_f42436_d33a.tsv", "d43.tsv", "d44.tsv", AMEND6, AMEND9, AMEND11]
H2K = 627.5094740631
EV = 27.211386245988
TOL_EH = 1e-6          # Eh: how closely a control must give the committed energy of the same input (D47's)
TOL_REL = 0.15         # D44's tolerance: 15 % of the working-level ddE (STAGE6_CRITERIA.txt)
FLOOR = 3.2            # kcal/mol: below this working-level ddE the agreement is NOT ASSESSABLE (STAGE6_CRITERIA.txt)
COUNT_TOL = 0.02       # e: D53's bound on the electrons the grid integrates; reported, not a leg
NEL = {"di": 118, "ctl": 118, "mo": 117, "bare": 117}
COMMITTED_CTL = {("SV", "ctl", "R"): "d31/ctl_R", ("SV", "ctl", "TS"): "d31/ctl_TS", ("SV", "bare", "R"): "d40/B_bare_R", ("SV", "bare", "TS"): "d40/B_bare_TS",
                 ("TZ", "ctl", "R"): "d44/CAMTZ_C_ctl_R", ("TZ", "ctl", "TS"): "d44/CAMTZ_C_ctl_TS", ("TZ", "bare", "R"): "d44/CAMTZ_C_bare_R", ("TZ", "bare", "TS"): "d44/CAMTZ_C_bare_TS"}
RECORD = ("B_n02_q0.50", "B_n02_q0.75", "B_n02_q1.00", "B_n03_q0.75", "B_n05_q0.50", "B_n11_q0.25")     # D43's six designs


def stop(m):
    sys.stderr.write("STOP: %s\n" % m)
    sys.exit(1)


def read_run(M, d, kind):
    """one run: D53's read_di for a dianion; a monoanion is read for its energy and its stability alone"""
    nel = NEL[kind]
    r = M["d40"].read(d, "pc" if kind in ("di", "mo") else "plain", nel, M["d33b"])
    r["status"], r["c"], r["lw"], r["nbf"], r["count"], r["kind"] = "measured", None, None, None, None, kind
    r.setdefault("stab", dict(ran=False, first=None, last=None))
    r.setdefault("trah", False)
    r.setdefault("lshift", False)
    if not r.get("exists"):
        r["status"] = "NOT MEASURED"
        return r
    r["nbf"], r["count"] = M["r53"].nbf(r["text"]), M["r53"].count(r["text"])
    if r["E"] is None:
        r["status"] = "DID NOT CONVERGE" if "SCF NOT CONVERGED" in r["text"] else "NOT MEASURED"
        return r
    if nel != 118:
        return r
    c, lw = M["audit"].read(r["text"]), M["r53"].loewdin(r["text"])
    if c is None or c.get("maxlw") is None or c.get("dip") is None or c.get("nb") is None or len(lw) < 24:
        r["status"], r["E"] = "NOT MEASURED", None
        r["why"] = ["the committed classifier does not read a Loewdin charge of each atom, a dipole and a basis-function count from it"]
        return r
    r["c"], r["lw"] = c, lw
    return r


def four_legs(M, di, ctl, geom, charges):
    """(the word of D43, {state: legs}) for a design at one level: di and ctl are {state: run}"""
    T = (M["audit"].D_LW, M["audit"].R_DIP, M["audit"].GAP_ABS, M["audit"].R_GAP)
    z = {}
    for st in ("R", "TS"):
        z[st] = M["r53"].legs(di[st], ctl[st], geom[st], charges, T) if ctl[st]["status"] == "measured" and di[st]["status"] in ("measured", "DID NOT CONVERGE") else dict(M["r53"].NOZ)
    if any(z[st]["flags"] for st in z):
        return "UNRELIABLE", z
    if any(di[st]["status"] != "measured" or ctl[st]["status"] != "measured" for st in z):
        return "NOT CLASSIFIED", z
    if any(z[st]["held"] for st in z):
        return "INCONCLUSIVE", z
    return "RELIABLE", z


def fifth_leg(di, mo):
    """(PASSED or NOT PASSED: why, {state: (VDE in eV or None, word)}) by AMENDMENT 6"""
    w = {}
    for st in ("R", "TS"):
        if di[st]["status"] != "measured" or mo[st]["status"] != "measured":
            w[st] = (None, "NOT MEASURED")
            continue
        v = (mo[st]["E"] - di[st]["E"]) * EV
        w[st] = (v, "METASTABLE" if v <= 0 else ("BOUND" if mo[st]["stab"]["last"] == "STABLE" else "ABOVE ZERO, MONOANION NOT STABLE"))
    if all(w[st][1] == "BOUND" for st in w):
        return "PASSED", w
    return "NOT PASSED: " + ", ".join(sorted(set(w[st][1] for st in w if w[st][1] != "BOUND"))), w


def dde_of(di, ctl):
    if any(di[st]["status"] != "measured" or ctl[st]["status"] != "measured" for st in ("R", "TS")):
        return None
    return ((di["TS"]["E"] - ctl["TS"]["E"]) - (di["R"]["E"] - ctl["R"]["E"])) * H2K


def agreement(dde, work):
    """D44's word for a ddE against the working level's"""
    if dde is None or work is None:
        return "NOT FORMED"
    if abs(work) < FLOOR:
        return "NOT ASSESSABLE"
    return "WITHIN" if abs(dde - work) < TOL_REL * abs(work) else "OUTSIDE"


def at_level(four, fifth, fw, agree):
    """D44's word at a level"""
    if four == "RELIABLE" and fifth == "PASSED" and agree == "WITHIN":
        return "PASSES"
    if four == "UNRELIABLE" or any(fw[st][1] == "METASTABLE" for st in fw) or agree == "OUTSIDE":
        return "FAILS"
    return "OPEN"


def pc_charges(d):
    L = [l.split() for l in open(os.path.join(d, "pc.pc")).read().splitlines()[1:] if l.strip()]
    return [(float(f[0]), (float(f[1]), float(f[2]), float(f[3]))) for f in L]


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    a.res, a.code, a.prereg = (os.path.abspath(x) for x in (a.res, a.code, a.prereg))
    S6 = os.path.join(a.res, S)
    D = os.path.join(S6, "d56")
    paths = dict((k, os.path.join(a.code, "phase2.2", fn)) for k, fn in (
        ("d40", "cs6_d40.py"), ("d33b", "cs6_d33b.py"), ("audit", "cs6_audit.py"), ("r53", "cs6_d53.py"), ("p47", "cs6_d47_prepare.py"), ("p56", "cs6_d56_prepare.py"), ("opt", "cs_opt_milp.py")))
    for k, p in paths.items():
        if not os.path.isfile(p):
            stop("the committed script %s is missing." % p)
    for fn in [a.prereg, os.path.join(D, "RUNS.tsv"), os.path.join(D, "DESIGNS.tsv")] + [os.path.join(S6, x) for x in DOCS]:
        if not os.path.isfile(fn):
            stop("%s is missing." % fn)
    tsvs = dict((k, os.path.join(os.path.dirname(a.tsv), "d56_%s.tsv" % k)) for k in ("designs", "survivors"))
    for fn in [a.result, a.tsv] + list(tsvs.values()):
        if os.path.exists(fn):
            stop("%s exists already; this reader overwrites nothing." % fn)
    import importlib.util

    def load(path, name):
        sp = importlib.util.spec_from_file_location(name, path)
        m = importlib.util.module_from_spec(sp)
        sp.loader.exec_module(m)
        return m
    M = dict((k, load(paths[k], "d56_" + k)) for k in ("d40", "d33b", "audit", "r53", "p47", "p56"))
    R, p47, p56 = M["r53"], M["p47"], M["p56"]
    gates = []

    def git(*args):
        try:
            q = subprocess.run(["git", "-C", a.res] + list(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            return q.returncode, q.stdout.strip()
        except OSError:
            return 127, ""

    # ------------------------------------------------------------ G56.3: the committed files this leans on
    guard = ["%s/%s" % (S, g) for g in GUARDED] + ["%s/%s" % (S, x) for x in DOCS] + ["%s/d56/DESIGNS.tsv" % S, "%s/d56/INPUT_CHECKSUMS.sha256" % S]
    rc, st = git("status", "--porcelain", "--untracked-files=all", "--", *guard)
    if rc != 0:
        gates.append("G56.3 could not be checked: git did not run")
    elif st:
        gates.append("G56.3 a file under d31, d40, d43, d44 or d53, a committed table or document this reader names, or the committed list of designs or of inputs of d56, is not the committed one")
    rc, pcm = git("log", "--diff-filter=A", "--format=%H", "-1", "--", "%s/%s" % (S, PREREG))
    rc_d, dcm = git("log", "--diff-filter=A", "--format=%H", "-1", "--", "%s/d56/DESIGNS.tsv" % S)
    if rc != 0 or not pcm:
        gates.append("G56.3 the commit that added %s could not be found" % PREREG)
    else:
        rc1, b1 = git("rev-parse", "%s:%s/%s" % (pcm, S, PREREG))
        rc2, b2 = git("hash-object", a.prereg)
        if rc1 != 0 or rc2 != 0 or b1 != b2:
            gates.append("G56.3 the pre-registration given is not the committed %s" % PREREG)
        if rc_d != 0 or not dcm:
            gates.append("G56.3 the commit that added d56/DESIGNS.tsv could not be found")
        elif dcm == pcm or git("merge-base", "--is-ancestor", pcm, dcm)[0] != 0:
            gates.append("G56.3 the designs were not committed after the pre-registration")

    # ------------------------------------------------------------ G56.1 and G56.2: formed again from the committed designs
    try:
        F, designs, runs = p56.build(S6, a.code, designs_from=os.path.join(D, "DESIGNS.tsv"))
    except SystemExit:
        stop("G56.1 the designs of d56/DESIGNS.tsv do not obey the rules of the search (the line above says which). Nothing was written.")
    except Exception as e:
        sys.stdout.write("\n".join("FAILED  " + g for g in gates + ["G56.1 the files of d56 could not be formed again: %s: %s" % (type(e).__name__, e)]) + "\n")
        stop("%d gate(s) failed. Nothing was written." % (len(gates) + 1))
    nsrc = nsum = nbad = 0
    psrc = os.path.join(D, "SOURCES.sha256")
    if not os.path.isfile(psrc) or open(psrc).read() != F["SOURCES.sha256"]:
        gates.append("G56.2 d56/SOURCES.sha256 is not the list formed again from the committed files")
    for l in F["SOURCES.sha256"].splitlines():
        h, fn = l.split(None, 1)
        nsrc += 1
        if not os.path.isfile(os.path.join(S6, fn)) or R.sha(os.path.join(S6, fn)) != h:
            gates.append("G56.2 %s is not the file d56/SOURCES.sha256 lists" % fn)
    pic = os.path.join(D, "INPUT_CHECKSUMS.sha256")
    if not os.path.isfile(pic) or open(pic).read() != F["INPUT_CHECKSUMS.sha256"]:
        gates.append("G56.1 d56/INPUT_CHECKSUMS.sha256 is not the list formed again here")
    for k in sorted(F):
        p = os.path.join(D, k)
        nsum += 1
        if not os.path.isfile(p) or open(p).read() != F[k]:
            nbad += 1
            if nbad <= 5:
                gates.append("G56.1 d56/%s is not the file formed again from the committed files and the committed designs" % k)
    if nbad > 5:
        gates.append("G56.1 and %d more files of d56 are not the files formed again" % (nbad - 5))

    # ------------------------------------------------------------ G56.4: the quoted sentences
    def has(path, sent):
        return " ".join(sent.split()) in " ".join(open(path, errors="replace").read().split())
    if not has(os.path.join(S6, AMEND9), Q_A9):
        gates.append("G56.4 the words of A9.4 (i) are not in %s" % AMEND9)
    if not has(os.path.join(S6, AMEND6), Q_A6):
        gates.append("G56.4 the words of A6.2 are not in %s" % AMEND6)
    for doc in (os.path.join(S6, AMEND11), a.prereg):
        if not has(doc, Q_J):
            gates.append("G56.4 Jessica's words of 24 November are not in %s" % os.path.basename(doc))

    geom = dict((st, p47.read_geom(os.path.join(S6, p47.GEOM[st]))) for st in ("R", "TS"))

    # ------------------------------------------------------------ G56.7: the same functions on the committed record
    rec = []
    d43 = dict(((r["condition"], r["state"]), r) for r in csv.DictReader(open(os.path.join(S6, "d43.tsv")), delimiter="\t"))
    d44 = dict(((r["design"], r["level_key"], r["state"]), r) for r in csv.DictReader(open(os.path.join(S6, "d44.tsv")), delimiter="\t"))
    cc = dict((k, read_run(M, os.path.join(S6, v), k[1])) for k, v in COMMITTED_CTL.items())
    nrec = 0
    for des in RECORD:
        work = None
        for lv, pre, tab in (("SV", "d43/", d43), ("TZ", "d44/CAMTZ_", d44)):
            di = dict((st, read_run(M, os.path.join(S6, "%s%s_di_%s" % (pre, des, st)), "di")) for st in ("R", "TS"))
            mo = dict((st, read_run(M, os.path.join(S6, "%s%s_mo_%s" % (pre, des, st)), "mo")) for st in ("R", "TS"))
            ctl = dict((st, cc[(lv, "ctl", st)]) for st in ("R", "TS"))
            four, z = four_legs(M, di, ctl, geom, pc_charges(os.path.join(S6, "%s%s_di_R" % (pre, des))))
            fifth, fw = fifth_leg(di, mo)
            dde = dde_of(di, ctl)
            row = tab[(des, "R")] if lv == "SV" else tab[(des, "CAMTZ", "R")]
            got = [four, fifth, "" if dde is None else "%+.4f" % dde]
            want = [row["four_legs"], row["fifth_leg"], row["ddE_kcal"]]
            if lv == "SV":
                work = dde
            else:
                ag = agreement(dde, work)
                got += [ag, at_level(four, fifth, fw, ag)]
                want += [row["agreement"], row["at_this_level"]]
            nrec += len(got)
            for g, w, what in zip(got, want, ("the four legs", "the fifth leg", "ddE", "the agreement", "the word at the level")):
                same = g == w or (what == "the fifth leg" and g.split(":")[0] == w.split(":")[0]) or (what == "ddE" and g and w and abs(float(g) - float(w)) <= 1.5e-4)
                if not same:
                    rec.append("%s at %s, %s: these functions give '%s', the committed table '%s'" % (des, "the working level" if lv == "SV" else "CAMTZ", what, g, w))
    for r_ in rec[:6]:
        gates.append("G56.7 " + r_)
    if len(rec) > 6:
        gates.append("G56.7 and %d more words or numbers of D43 and D44 are not reproduced" % (len(rec) - 6))

    # ------------------------------------------------------------ reading every run, and G56.5, G56.6
    nbf = dict((lv, nb) for lv, _, nb, _ in p56.LEVELS)
    RUN = {}
    for r in runs:
        d = os.path.join(D, r["run"])
        x = read_run(M, d, r["kind"])
        x.update(r)
        x["count_ok"] = x["status"] == "measured" and x["count"] is not None and abs(x["count"] - NEL[r["kind"]]) <= COUNT_TOL
        RUN[r["run"]] = x
        if not x.get("exists"):
            continue
        e = R.echoed(x["text"])
        want = open(os.path.join(d, "job.inp")).read().split("\n") if os.path.isfile(os.path.join(d, "job.inp")) else None
        while want and not want[-1].strip():
            want.pop()
        if e is None:
            if x.get("normal"):
                gates.append("G56.5 %s/job.out ends normally and does not print the input it ran" % r["run"])
        elif want is not None and e != want:
            gates.append("G56.5 the input printed at the top of %s/job.out is not the committed %s/job.inp" % (r["run"], r["run"]))
        if x["status"] == "measured" and x["nbf"] != nbf[r["level"]]:
            gates.append("G56.6 %s reports %s basis functions; its level has %d on the substrate" % (r["run"], x["nbf"], nbf[r["level"]]))
    if gates:
        sys.stdout.write("\n".join("FAILED  " + g for g in gates) + "\n")
        stop("%d gate(s) failed. Nothing was written." % len(gates))

    # ------------------------------------------------------------ the readings
    LV = [lv for lv, _, _, _ in p56.LEVELS]
    lname = dict((lv, nm) for lv, nm, _, _ in p56.LEVELS)
    C = dict(((lv, k, st), RUN["%s_C_%s_%s" % (lv, k, st)]) for lv in LV for k in ("ctl", "bare") for st in ("R", "TS"))
    DES = []
    for d_ in designs:
        o = dict(d_)
        charges = list(zip(d_["q"], d_["xyz"]))
        for lv in LV:
            di = dict((st, RUN["%s_%s_di_%s" % (lv, d_["name"], st)]) for st in ("R", "TS"))
            mo = dict((st, RUN["%s_%s_mo_%s" % (lv, d_["name"], st)]) for st in ("R", "TS"))
            ctl = dict((st, C[(lv, "ctl", st)]) for st in ("R", "TS"))
            four, z = four_legs(M, di, ctl, geom, charges)
            fifth, fw = fifth_leg(di, mo)
            o[lv] = dict(four=four, z=z, fifth=fifth, fw=fw, dde=dde_of(di, ctl), di=di, mo=mo)
        o["agree"] = agreement(o["TZ"]["dde"], o["SV"]["dde"])
        o["at_tz"] = at_level(o["TZ"]["four"], o["TZ"]["fifth"], o["TZ"]["fw"], o["agree"])
        o["work"] = "PASSES THE FIVE LEGS" if o["SV"]["four"] == "RELIABLE" and o["SV"]["fifth"] == "PASSED" else (
            "FAILS" if o["SV"]["four"] == "UNRELIABLE" or any(o["SV"]["fw"][st][1] == "METASTABLE" for st in ("R", "TS")) else "OPEN")
        o["word"] = "SURVIVOR" if o["work"] == "PASSES THE FIVE LEGS" and o["at_tz"] == "PASSES" else ("NOT A SURVIVOR: FAILS" if "FAILS" in (o["work"], o["at_tz"]) else "NOT A SURVIVOR: OPEN")
        DES.append(o)
    BEST = {}
    for n in p56.NS:
        s_ = [o for o in DES if o["n"] == n and o["word"] == "SURVIVOR"]
        BEST[n] = min(s_, key=lambda o: o["k"]) if s_ else None
    surv = sorted((o for o in DES if o["word"] == "SURVIVOR"), key=lambda o: o["SV"]["dde"])

    # ------------------------------------------------------------ writing
    out = []
    say = out.append

    def para(text, lead=""):
        for l in textwrap.wrap(" ".join(text.split()), 96 - len(lead)):
            say(lead + l)

    def head(t):
        say(""); say("=" * 96); say(t); say("=" * 96); say("")
    f = lambda v, fmt="%+.4f": "-" if v is None else fmt % v
    rc, hd = git("rev-parse", "--short=12", "HEAD")
    say("STAGE 6 - DIAGNOSTIC 56 - RESULT")
    say("The search on the committed grid: designs of at most %.2f e a site, each measured at two levels" % p56.CAP_Q)
    say("=" * 96)
    say("")
    para("Read by cs6_d56.py under %s, which was committed before the search was run and before any job was submitted." % PREREG)
    say("   results repository at %s" % hd)
    say("   the pre-registration was added by %s, the designs by %s" % (pcm[:12], dcm[:12]))
    for label, k in (("reader", "cs6_d56.py"), ("search and inputs", "cs6_d56_prepare.py"), ("the optimiser", "cs_opt_milp.py"), ("reading functions", "cs6_d53.py"),
                     ("outputs read by", "cs6_d40.py"), ("stability by", "cs6_d33b.py"), ("the four legs by", "cs6_audit.py")):
        say("   %-20s %s" % (label, k))
        say("   %-20s sha256 %s" % ("", R.sha(os.path.join(a.code, "phase2.2", k))))
    say("")
    para("""Bare point charges of at most %.2f e, no basis function at a charge, frame 42436, frozen geometries, in vacuo. Working level:
%s. Second level: %s. A design is a SURVIVOR OF THE SEARCH, and nothing more, by A9.4 (i). D56 adopts no level of theory and names no
design a result.""" % (p56.CAP_Q, lname["SV"], lname["TZ"]))
    say("")
    say("THE GATES, as they held:")
    for g, t in (("G56.1", "all %d files of d56 are the ones formed again here from the committed files and the committed designs, and every design obeys the rules of the search." % nsum),
                 ("G56.2", "the %d files of d56/SOURCES.sha256 are the ones listed." % nsrc),
                 ("G56.3", "no guarded file differs from the committed one; the pre-registration given is the committed one, and the designs were committed after it."),
                 ("G56.4", "the words quoted from AMENDMENT 6, AMENDMENT 9 and Jessica are in the documents they are quoted from."),
                 ("G56.5", "every output that exists prints, as the input it ran, the committed input beside it."),
                 ("G56.6", "every measured run reports the basis functions of its level (%d and %d)." % (nbf["SV"], nbf["TZ"])),
                 ("G56.7", "these reading functions, applied to the committed outputs of D43 and D44, give all %d committed words and numbers of the six designs at the working level and at CAMTZ." % nrec)):
        ls = textwrap.wrap(t, 96 - 10)
        say("   %-6s %s" % (g, ls[0]))
        for l in ls[1:]:
            say("          " + l)

    head("R56.1  THE RUNS AND THE CONTROLS")
    for lv in LV:
        for kind, label in (("ctl", "controls, dianion"), ("bare", "controls, monoanion"), ("di", "dianion with a design"), ("mo", "monoanion with a design")):
            mine = [RUN[r["run"]] for r in runs if r["kind"] == kind and r["level"] == lv]
            say("   %-2s %-24s measured %3d of %3d" % (lv, label, sum(1 for x in mine if x["status"] == "measured"), len(mine)))
    nm = sum(1 for r in runs if RUN[r["run"]]["status"] == "measured")
    say("")
    bad = [RUN[r["run"]] for r in runs if RUN[r["run"]]["status"] != "measured"]
    if bad:
        say("   Not measured (%d):" % len(bad))
        for x in bad[:40]:
            para("%s: %s%s" % (x["run"], x["status"], (": " + "; ".join(x.get("why") or [])) if x.get("why") else ""), "      ")
        if len(bad) > 40:
            say("      ... and %d more; d56.tsv lists every run." % (len(bad) - 40))
    else:
        say("   Every one of the %d runs is measured." % len(runs))
    off = [RUN[r["run"]] for r in runs if RUN[r["run"]]["status"] == "measured" and not RUN[r["run"]]["count_ok"]]
    para("Runs whose count of electrons is off by more than %.2f e (reported, not a leg): %d%s." % (COUNT_TOL, len(off), (": " + ", ".join(x["run"] for x in off[:12])) if off else ""), "   ")
    say("")
    say("   control          E, Eh             against the committed run of the same input")
    for lv in LV:
        for kind in ("ctl", "bare"):
            for stt in ("R", "TS"):
                x, ref = C[(lv, kind, stt)], cc[(lv, kind, stt)]
                cmp_ = "-" if x["status"] != "measured" or ref.get("E") is None else (("reproduces it" if abs(x["E"] - ref["E"]) <= TOL_EH else "DIFFERS") + " (%+.1e Eh)" % (x["E"] - ref["E"]))
                say("   %-16s %-17s %s" % (x["run"], f(x["E"], "%.8f"), cmp_))
    say("")
    for lv in LV:
        cb = dde_of(dict((st, C[(lv, "ctl", st)]) for st in ("R", "TS")), dict((st, dict(status="measured", E=0.0)) for st in ("R", "TS")))
        vd = [((C[(lv, "bare", st)]["E"] - C[(lv, "ctl", st)]["E"]) * EV if C[(lv, "bare", st)]["status"] == C[(lv, "ctl", st)]["status"] == "measured" else None) for st in ("R", "TS")]
        say("   %s: barrier without charges %s kcal/mol" % (lv, f(cb)))
        say("       binding of the last electron without charges %s eV at R, %s eV at TS" % (f(vd[0]), f(vd[1])))

    head("R56.2  EACH DESIGN AT EACH LEVEL")
    para("""'four legs': RELIABLE, UNRELIABLE, INCONCLUSIVE or NOT CLASSIFIED, as D43 words them. 'fifth': the dianion keeps its last electron
at both states (AMENDMENT 6). 'ddE': the change of the barrier, kcal/mol; negative lowers it. 'score': the first-order value the search
ranked by.""")
    for n in p56.NS:
        say("")
        say("   %d charges" % n)
        say("   design      lv  four legs       fifth leg                      ddE        score      net e")
        for o in [o for o in DES if o["n"] == n]:
            for lv in LV:
                say(("   %-11s %-3s %-15s %-30s %-10s %-10s %s" % (o["name"] if lv == "SV" else "", lv, o[lv]["four"], o[lv]["fifth"][:30], f(o[lv]["dde"]), f(o["score"]) if lv == "SV" else "", ("%+.2f" % o["net"]) if lv == "SV" else "")).rstrip())
        for o in [o for o in DES if o["n"] == n]:
            for lv in LV:
                for stt in ("R", "TS"):
                    z = o[lv]["z"][stt]
                    if z["flags"] or z["held"]:
                        para("%s: %s." % (o[lv]["di"][stt]["run"], "; ".join(["fired: " + t for t in z["flags"]] + ["withheld: " + t for t in z["held"]])), "      ")

    head("R56.3  THE WORD OF EACH DESIGN, BY A9.4 (i)")
    para("""A design is a SURVIVOR only if it passes the five legs at the working level AND, at the second level, is RELIABLE, passes the
fifth leg and has a ddE WITHIN %d %% of its working-level ddE. NOT ASSESSABLE: the working-level ddE is below %.1f kcal/mol.""" % (round(100 * TOL_REL), FLOOR))
    say("")
    say("   design      working level         ddE TZ - SV  agreement       at TZ   the design")
    for o in DES:
        df = None if o["TZ"]["dde"] is None or o["SV"]["dde"] is None else o["TZ"]["dde"] - o["SV"]["dde"]
        say("   %-11s %-21s %-12s %-15s %-7s %s" % (o["name"], o["work"], f(df), o["agree"], o["at_tz"], o["word"]))

    head("R56.4  THE SURVIVORS")
    if not surv:
        say("   NO DESIGN IS A SURVIVOR.")
    else:
        say("   Ordered by the working-level ddE, the most lowering first.")
        say("")
        say("   design      charges  net e   positive e  ddE SV     ddE TZ     nearest atom, A")
        for o in surv:
            dmin = min(min(p47.dist(x, g[1:]) for g in geom[st]) for x in o["xyz"] for st in ("R", "TS"))
            say("   %-11s %-8d %+.2f   %-11.2f %-10s %-10s %.2f" % (o["name"], o["n"], o["net"], sum(v for v in o["q"] if v > 0), f(o["SV"]["dde"]), f(o["TZ"]["dde"]), dmin))
    say("")
    say("   At each number of charges, the best-ranked design that is a SURVIVOR:")
    for n in p56.NS:
        say("   %2d charges: %s" % (n, BEST[n]["name"] if BEST[n] else "none of the %d" % p56.POOL))

    head("WHAT THIS CANNOT SAY (section 6 of the pre-registration, in short; its own words govern)")
    for t in ("One frame, frozen geometries, in vacuo. No barrier here is a prediction of a physical barrier.",
              "The charges are bare points and stand for no residue. No basis function stands at a charge and no wall is in any run: by AMENDMENT 11 the test with functions and the helium wall are run on the survivors afterwards and reported beside them.",
              "A SURVIVOR is not a result. The cross-check of AMENDMENT 7 and the reference of AMENDMENT 8 have not been run on it.",
              "Three designs at each number of charges, ranked at first order. A design the search did not keep is not measured.",
              "D56 adopts nothing."):
        para(t, "   ")
    say("")
    cnt = {}
    for o in DES:
        cnt[o["word"]] = cnt.get(o["word"], 0) + 1
    summ = "; ".join("%s: %d" % (k, cnt[k]) for k in ("SURVIVOR", "NOT A SURVIVOR: FAILS", "NOT A SURVIVOR: OPEN") if k in cnt)
    para("Runs measured: %d of %d. The %d designs: %s." % (nm, len(runs), len(DES), summ))

    res = "\n".join(out) + "\n"
    wide = [l for l in res.split("\n") if len(l) > 96]
    if wide:
        stop("a line of the result is %d characters wide. Nothing was written." % max(len(l) for l in wide))
    g = lambda v, fmt="%+.6f": "" if v is None else fmt % v
    T1 = ["run\tjob\tlevel\tkind\tdesign\tstate\tstatus\tE_Eh\tcount_minus_electrons\tbasis_functions\tstability_first\tstability_last\tTRAH\tlevel_shift_set"]
    for r in runs:
        x = RUN[r["run"]]
        T1.append("\t".join([r["run"], r["job"], r["level"], r["kind"], r["design"] or "-", r["state"], x["status"], g(x["E"], "%.10f"), g(x["count"] - NEL[r["kind"]]) if x["count"] is not None else "",
                             str(x["nbf"] or ""), x["stab"]["first"] or "", x["stab"]["last"] or "", "yes" if x["trah"] else "no", "yes" if x["lshift"] else "no"]))
    T2 = ["design\tn\trank_at_n\tsum_q_e\tfirst_order_score_kcal\tlevel\tfour_legs\tflags_R\tflags_TS\twithheld_R\twithheld_TS\tlowdin_rise_R_e\tlowdin_rise_TS_e\tdipole_ratio_R\tdipole_ratio_TS\tgap_R_eV\tgap_TS_eV\tVDE_R_eV\tVDE_TS_eV\tbinding_R\tbinding_TS\tfifth_leg\tddE_kcal"]
    for o in DES:
        for lv in LV:
            v, z = o[lv], o[lv]["z"]
            T2.append("\t".join([o["name"], str(o["n"]), str(o["k"]), "%+.4f" % o["net"], "%.6f" % o["score"], lv, v["four"], "; ".join(z["R"]["flags"]), "; ".join(z["TS"]["flags"]), "; ".join(z["R"]["held"]), "; ".join(z["TS"]["held"]),
                                 g(z["R"]["rise"]), g(z["TS"]["rise"]), g(z["R"]["ratio"], "%.4f"), g(z["TS"]["ratio"], "%.4f"), g(z["R"]["gap"], "%.4f"), g(z["TS"]["gap"], "%.4f"),
                                 g(v["fw"]["R"][0]), g(v["fw"]["TS"][0]), v["fw"]["R"][1], v["fw"]["TS"][1], v["fifth"], g(v["dde"])]))
    T3 = ["design\tn\trank_at_n\tworking_level\tddE_SV_kcal\tddE_TZ_kcal\tagreement\tat_TZ\tthe_design\tbest_survivor_at_its_n"]
    for o in DES:
        T3.append("\t".join([o["name"], str(o["n"]), str(o["k"]), o["work"], g(o["SV"]["dde"]), g(o["TZ"]["dde"]), o["agree"], o["at_tz"], o["word"], "yes" if BEST[o["n"]] is o else "no"]))
    for fn, t in ((a.result, res), (a.tsv, "\n".join(T1) + "\n"), (tsvs["designs"], "\n".join(T2) + "\n"), (tsvs["survivors"], "\n".join(T3) + "\n")):
        with open(fn, "w") as fh:
            fh.write(t)
    sys.stdout.write(res)


if __name__ == "__main__":
    main()
