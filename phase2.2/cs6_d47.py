#!/usr/bin/env python3
"""DIAGNOSTIC 47's reader. A preliminary: when does a bare positive charge take electrons from the substrate?

Every run is a dianion single point on frame 42436. A run with a charge is read against a charge-free control of
the same level and state: at CAMSV and CAMTZ the COMMITTED control (d31/ctl, d44/CAMTZ_C_ctl); at CAMTZG the
ghost control of this diagnostic with the same ghost centre; at CAMQZ the control of this diagnostic.

THE LEGS, on each run, by the committed classifier cs6_audit.py at its committed thresholds, with the rules D33.3
and AMENDMENT 5 settled and cs6_d33b.py holds: the Loewdin rise, the dipole ratio, the gap (withheld in a run that
engaged TRAH or a level shift, or whose control did) and the LAST stability analysis. A dianion output that says
the SCF did not converge is a flag (D31.3). The fifth leg is not read: no monoanion with a charge is run.

WHERE THE CHARGE WENT, on each run: the Loewdin charge of each of the 24 substrate atoms minus the control's; the
atom with the most negative change, that change, and its distance to the nearest positive point charge; and the
sum over the 24 atoms, which at CAMTZG says how much sits on the functions at the charge.

THE MAPS: orca_vpot at every committed site and every rung, on the charge-free dianion and monoanion of CAMSV and
CAMTZ at R and TS. From them, by formulas fixed in the pre-registration, section 2: the first-order coefficient
a = V_TS - V_R of the dianion; v = V_monoanion - V_dianion, the first-order change of the detachment energy per unit
charge; and for a positive charge q_i among charges q_j the estimate
    S_i = VDE_bare + sum_j q_j v_j - EV [ V_monoanion(i) + sum_(j not i) q_j / r_ij ] - RYD q_i^2      (eV)
whose sign says whether an electron bound to the bare charge in a complete basis would lie below the electron the
dianion holds. S is an estimate. It is compared here with what the runs did; it decides nothing.

The inputs are not taken on trust: the committed cs6_d47_prepare.py is imported by path, its select() is run
again on the committed files, and every table, input, point-charge file and the points file are formed again and
compared with the files that ran. The input ORCA echoes at the top of each output must be that file.

usage: cs6_d47.py --res <results-repo> --code <code-repo> --prereg <file> --result <file> --tsv <file>
"""
import argparse, importlib.util, os, re, sys, textwrap

EV = 27.211386245988
RYD = EV / 2.0
H2K = 627.5094740631
S = "phase2.2/cs_stage6"
NEL_MONO, NEL_DI = 117, 118
TOL_EH = 1e-6
FRAC, FLOOR = 0.15, 3.2
S_BAND = 0.5                    # eV: S at or above +S_BAND predicts no flag, at or below -S_BAND predicts a flag
PREREG = "STAGE6_DIAGNOSTIC47_20261121.txt"
AMEND5 = "STAGE6_AMENDMENT5_20261109.txt"
AMEND9 = "STAGE6_AMENDMENT9_20261120.txt"
CRIT = "STAGE6_CRITERIA.txt"
CRIT_A = "LOSS MUST BE BELOW 15 % OF |ddE(frozen)|."
CRIT_C = "if |ddE(frozen)| is smaller than 3.2 kcal/mol"
A9_A = ("A design is a SURVIVOR OF THE SEARCH, and nothing more, only if it passes the five legs at the working level AND, at "
        "CAM-B3LYP-D3BJ/def2-TZVPD, is RELIABLE on four legs, passes the fifth leg, and has a ddE that is WITHIN the tolerance of "
        "its ddE at the working level.")
D44RES = "STAGE6_DIAGNOSTIC44_RESULT_20261120.txt"
D44_A = "Of the 6 designs. Pass the cross-check: 0. Do not pass: 4. Incomplete: 2."
GUARDED = ("d31", "d33a", "d37", "d40", "d43", "d44")
STATES = ("R", "TS")


def stop(m):
    sys.stderr.write("STOP  %s\n" % m)
    raise SystemExit(1)


def load(path, name):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def echoed(text):
    m = re.search(r"INPUT FILE\n=+\nNAME = [^\n]*\n((?:\|\s*\d+>[^\n]*\n)+)", text)
    if not m:
        return None
    ls = [re.sub(r"^\|\s*\d+> ?", "", l) for l in m.group(1).split("\n") if l]
    while ls and (not ls[-1].strip() or "****END OF INPUT****" in ls[-1]):
        ls.pop()
    return ls


def loewdin(text):
    """the charges of the last LOEWDIN ATOMIC CHARGES block, in the order printed"""
    k = text.rfind("LOEWDIN ATOMIC CHARGES")
    q = []
    if k >= 0:
        for l in text[k:].split("\n")[2:]:
            m = re.match(r"\s*\d+\s*[A-Za-z]{1,2}\s*:\s*(-?\d+\.\d+)", l)
            if not m:
                break
            q.append(float(m.group(1)))
    return q


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    a.res, a.code, a.prereg = (os.path.abspath(x) for x in (a.res, a.code, a.prereg))
    S6 = os.path.join(a.res, S)
    D = os.path.join(S6, "d47")
    paths = dict((k, os.path.join(a.code, "phase2.2", f)) for k, f in (
        ("d40", "cs6_d40.py"), ("d33b", "cs6_d33b.py"), ("audit", "cs6_audit.py"), ("prep", "cs6_d47_prepare.py")))
    for p in paths.values():
        if not os.path.isfile(p):
            stop("%s is missing; this reader imports it by path." % p)
    d40 = load(paths["d40"], "d47_cs6_d40")
    d33b = load(paths["d33b"], "d47_cs6_d33b")
    audit = load(paths["audit"], "d47_cs6_audit")
    prep = load(paths["prep"], "d47_cs6_d47_prepare")
    sha, git = d40.sha, d40.git
    D_LW, R_DIP, GAP_ABS, R_GAP = audit.D_LW, audit.R_DIP, audit.GAP_ABS, audit.R_GAP
    BOHR, NS = prep.BOHR, prep.NSITE

    rep, gates = [], []
    say = rep.append
    say("STAGE 6, DIAGNOSTIC 47 - RESULT")
    say("A preliminary: when does a bare positive charge take electrons from the substrate - its size, its")
    say("place, its distance, and the reach of the basis set?")
    say("")
    say("Read by phase2.2/cs6_d47.py                 sha256 %s" % sha(__file__))
    say("   every output by cs6_d40.py read           sha256 %s" % sha(paths["d40"]))
    say("   stability by cs6_d33b.py stability        sha256 %s" % sha(paths["d33b"]))
    say("   thresholds and legs by cs6_audit.py       sha256 %s" % sha(paths["audit"]))
    say("   the inputs by cs6_d47_prepare.py          sha256 %s" % sha(paths["prep"]))
    say("   pre-registration %-23s sha256 %s" % (os.path.basename(a.prereg), sha(a.prereg)))
    say("")
    say("Frame 42436, bare point charges, dianion. Legs at D_LW %.2f e, R_DIP %.1f, GAP_ABS %.1f eV, R_GAP %.1f," % (D_LW, R_DIP, GAP_ABS, R_GAP))
    say("stability on the last analysis (AMENDMENT 5). %.12f eV per Eh, %.10f kcal/mol per Eh," % (EV, H2K))
    say("%.10f A per bohr. S is an estimate (section 2 of the pre-registration) and decides nothing." % BOHR)
    say("")

    def has(path, sent):
        return os.path.isfile(path) and sent in " ".join(open(path).read().split())

    # ------------------------------------------------------------ G47.2
    sel = None
    try:
        sel = prep.select(S6)
        runs = prep.runs_of(sel)
    except prep.Refused as e:
        gates.append("G47.2 the diagnostic could not be formed again from the committed files: %s" % e)
        runs = []
    except Exception as e:
        gates.append("G47.2 the diagnostic could not be formed again from the committed files: %s: %s" % (type(e).__name__, e))
        runs = []
    conds = sel["conds"] if sel else []
    cmap = dict((x["name"], x) for x in conds)
    if sel and runs:
        for fn, want in (("CONDITIONS.tsv", prep.conditions_text(sel)), ("RUNS.tsv", prep.runs_text(sel)), ("LEVELS.tsv", prep.levels_text(sel)),
                         ("points.xyz", prep.points_text(sel))):
            p = os.path.join(D, fn)
            if not os.path.isfile(p) or open(p).read() != want:
                gates.append("G47.2 d47/%s is not what the committed files give" % fn)
        for job in prep.JOBS:
            p = os.path.join(D, "run_%s.pbs" % job[0])
            if not os.path.isfile(p) or open(p).read() != prep.pbs_text(job, sel, runs, prep.HPC):
                gates.append("G47.2 d47/run_%s.pbs is not the job script of %s" % (job[0], job[0]))

    # ------------------------------------------------------------ G47.3
    psrc = os.path.join(D, "SOURCES.sha256")
    srcs = []
    if not os.path.isfile(psrc):
        gates.append("G47.3 d47/SOURCES.sha256 is missing")
    else:
        for l in open(psrc).read().splitlines():
            if l.strip():
                h, f = l.split(None, 1)
                srcs.append(f)
                p = os.path.join(S6, f)
                if not os.path.isfile(p) or sha(p) != h:
                    gates.append("G47.3 %s is not the file d47/SOURCES.sha256 lists" % f)
        if sel and srcs != prep.sources(sel):
            gates.append("G47.3 d47/SOURCES.sha256 does not list the files this diagnostic reads")
    spath = [os.path.normpath("%s/%s" % (S, x)) for x in srcs]
    rc, _ = git(a.res, "diff", "--quiet", "HEAD", "--", *(["%s/%s" % (S, x) for x in GUARDED] + spath))
    if rc is None:
        gates.append("G47.3 could not be checked: git did not run")
    elif rc != 0:
        gates.append("G47.3 a committed file under %s, or one this diagnostic reads, differs from HEAD" % ", ".join(GUARDED))
    rc, pcm = git(a.res, "log", "--diff-filter=A", "--format=%H", "-1", "--", "%s/%s" % (S, PREREG))
    if rc != 0 or not pcm:
        gates.append("G47.3 the commit that added %s could not be found" % PREREG)
        pcm = ""
    else:
        rc1, b1 = git(a.res, "rev-parse", "%s:%s/%s" % (pcm, S, PREREG))
        rc2, b2 = git(a.res, "hash-object", a.prereg)
        if rc1 != 0 or rc2 != 0 or b1 != b2:
            gates.append("G47.3 the pre-registration given is not the committed %s" % PREREG)
        for f, sp in zip(srcs, spath):
            rc1, b1 = git(a.res, "rev-parse", "%s:%s" % (pcm, sp))
            rc2, b2 = git(a.res, "hash-object", os.path.join(S6, f))
            if rc1 != 0 or rc2 != 0 or b1 != b2:
                gates.append("G47.3 %s is not the file at the commit the pre-registration is committed on" % f)

    # ------------------------------------------------------------ G47.5
    for q in (CRIT_A, CRIT_C):
        if not has(os.path.join(S6, CRIT), q):
            gates.append("G47.5 the words \"%s\" are not in %s" % (q, CRIT))
    if "%d %%" % round(100 * FRAC) not in CRIT_A or "%.1f kcal/mol" % FLOOR not in CRIT_C:
        gates.append("G47.5 the tolerance and the floor this reader uses are not the ones its quotations of %s state" % CRIT)
    if not has(os.path.join(S6, AMEND9), A9_A):
        gates.append("G47.5 the sentence of A9.4 is not in %s" % AMEND9)
    if not has(os.path.join(S6, D44RES), D44_A):
        gates.append("G47.5 the words \"%s\" are not in %s" % (D44_A, D44RES))
    for q in (CRIT_A, A9_A, D44_A):
        if not has(a.prereg, q):
            gates.append("G47.5 the words \"%s\" are not in the pre-registration given" % (q[:60] + (" ..." if len(q) > 60 else "")))

    # ------------------------------------------------------------ G47.1
    pic = os.path.join(D, "INPUT_CHECKSUMS.sha256")
    nsum = 0
    if not os.path.isfile(pic):
        gates.append("G47.1 d47/INPUT_CHECKSUMS.sha256 is missing")
    else:
        for l in open(pic).read().splitlines():
            if l.strip():
                h, f = l.split(None, 1)
                nsum += 1
                p = os.path.join(D, f)
                if not os.path.isfile(p) or sha(p) != h:
                    gates.append("G47.1 d47/%s is not the file INPUT_CHECKSUMS.sha256 lists" % f)
        if sel and runs and nsum != prep.nsum(sel, runs):
            gates.append("G47.1 INPUT_CHECKSUMS.sha256 lists %d files, not %d" % (nsum, prep.nsum(sel, runs)))
    for r in runs:
        d = os.path.join(D, r["run"])
        x = cmap.get(r["cond"])
        try:
            want = prep.input_text(S6, sel, r["level"], r["kind"], x, r["state"], r["gname"])
        except Exception as e:
            want = None
            gates.append("G47.1 the input of %s could not be formed again: %s" % (r["run"], e))
        pi = os.path.join(d, "job.inp")
        if want is not None and (not os.path.isfile(pi) or open(pi).read() != want):
            gates.append("G47.1 %s/job.inp is not the input the committed files give" % r["run"])
        px = os.path.join(d, "substrate.xyz")
        if r["level"] == "CAMTZG":
            if os.path.exists(px):
                gates.append("G47.1 %s holds a geometry file beside its written-out geometry" % r["run"])
        elif not os.path.isfile(px) or sha(px) != sha(os.path.join(S6, prep.GEOM[r["state"]])):
            gates.append("G47.1 %s/substrate.xyz is not %s" % (r["run"], prep.GEOM[r["state"]]))
        pp = os.path.join(d, "pc.pc")
        if x is not None:
            if not os.path.isfile(pp) or open(pp).read() != x["pc"]:
                gates.append("G47.1 %s/pc.pc does not hold the charges of %s" % (r["run"], x["name"]))
        elif os.path.exists(pp):
            gates.append("G47.1 the control %s holds a point-charge file" % r["run"])

    # ------------------------------------------------------------ reading
    def read_di(d, kind, nel=NEL_DI):
        r = d40.read(d, kind, nel, d33b)
        r["status"], r["c"], r["lw"] = "measured", None, None
        if not r.get("exists"):
            r["status"] = "NOT MEASURED"
            return r
        if r["E"] is None:
            r["status"] = "DID NOT CONVERGE" if "SCF NOT CONVERGED" in r["text"] else "NOT MEASURED"
            return r
        if nel == NEL_MONO:
            return r
        c = audit.read(r["text"])
        lw = loewdin(r["text"])
        if c is None or c.get("maxlw") is None or c.get("dip") is None or c.get("nb") is None or len(lw) < 24:
            r["status"], r["E"] = "NOT MEASURED", None
            r["why"] = ["the committed classifier does not read a Loewdin charge of each atom, a dipole and a basis-function count from it"]
            return r
        r["c"], r["lw"] = c, lw
        return r

    def check_echo(run, r):
        if not r.get("exists"):
            return
        e = echoed(r["text"])
        pin = os.path.join(D, run, "job.inp")
        if not os.path.isfile(pin):
            return
        want = open(pin).read().split("\n")
        while want and not want[-1].strip():
            want.pop()
        if e is None:
            if r.get("normal"):
                gates.append("G47.6 %s/job.out ends normally and does not print the input it ran" % run)
        elif e != want:
            gates.append("G47.6 the input printed at the top of %s/job.out is not the committed %s/job.inp" % (run, run))

    # ---- G47.4: the committed controls and outputs this reader leans on
    CT, BM = {}, {}
    lmap = dict((l["key"], l) for l in prep.LEVELS)
    for k in ("CAMSV", "CAMTZ"):
        for st in STATES:
            c = read_di(os.path.join(S6, "%s_%s" % (lmap[k]["ctl"], st)), "plain")
            b = read_di(os.path.join(S6, "%s_%s" % (lmap[k]["bare"], st)), "plain", NEL_MONO)
            CT[(k, st)], BM[(k, st)] = c, b
            if c["status"] != "measured" or not c["stab"].get("ran") or c["stab"].get("last") != "STABLE":
                gates.append("G47.4 the committed control %s_%s is not a measured dianion whose last stability analysis is STABLE" % (lmap[k]["ctl"], st))
            if b["status"] != "measured":
                gates.append("G47.4 the committed monoanion control %s_%s is not readable" % (lmap[k]["bare"], st))
    D44 = {}
    for d in prep.D44_DESIGNS:
        for st in STATES:
            r = read_di(os.path.join(S6, "d44", "CAMTZ_%s_di_%s" % (d, st)), "pc")
            D44[(d, st)] = r
    D43E = {}
    for x in conds:
        if x["part"] == "A" and abs(x["q"] - 1.0) < 1e-9:
            for st in STATES:
                D43E[(x["name"], st)] = read_di(os.path.join(S6, "d43", "%s_di_%s" % (x["name"], st)), "pc")
                if D43E[(x["name"], st)]["status"] != "measured":
                    gates.append("G47.4 D43's committed output of %s at %s is not readable" % (x["name"], st))

    def vread(p, pts):
        """the potentials of one vpot.out, Eh per e, in point order; a string where it cannot be read"""
        if not os.path.isfile(p):
            return "no vpot.out"
        L = [l.split() for l in open(p).read().splitlines() if l.strip()]
        if not L or len(L[0]) != 1 or not L[0][0].isdigit() or int(L[0][0]) != len(pts) or len(L) != len(pts) + 1:
            return "its count is not that of the points"
        out = []
        for row, ref in zip(L[1:], pts):
            try:
                v = [float(t) for t in row]
            except ValueError:
                return "a line is not four numbers"
            if len(v) != 4 or any(abs(v[j] - ref[j] / BOHR) > 1e-6 for j in range(3)):
                return "its points are not those of points.xyz"
            out.append(v[3])
        return out

    V33 = {}
    if sel:
        for st in STATES:
            V33[st] = vread(os.path.join(S6, prep.D33A_VPOT[st]), sel["points"][:NS])
            if isinstance(V33[st], str):
                gates.append("G47.4 the committed %s cannot be read: %s" % (prep.D33A_VPOT[st], V33[st]))

    # ---- every run of D47
    RUN = {}
    for r in runs:
        d = os.path.join(D, r["run"])
        x = read_di(d, "pc" if r["kind"] == "cond" else "plain", NEL_MONO if r["kind"] == "bare" else NEL_DI)
        check_echo(r["run"], x)
        x.update(r)
        RUN[r["run"]] = x

    def ref_of(x):
        """the control a run is read against: (the control's reading, its name, why it cannot be used or None)"""
        if x["ref"]:
            c = RUN[x["ref"]]
            why = None
            if c["status"] != "measured":
                why = "its control %s is %s" % (x["ref"], c["status"])
            elif not c["stab"].get("ran") or c["stab"].get("last") is None:
                why = "the stability analysis of its control %s was not read" % x["ref"]
            elif c["stab"]["last"] != "STABLE":
                why = "its control %s is UNSTABLE at the last analysis" % x["ref"]
            return c, x["ref"], why
        return CT[(x["level"], x["state"])], "%s_%s" % (lmap[x["level"]]["ctl"], x["state"]), None

    # ---- whether a job lets its runs be read: its plain control gives back the committed energy
    JUSE = {}
    for jk, lk, pbs in prep.JOBS:
        u = None
        if jk != "Q1":
            c = RUN.get("%s_C_ctl_R" % jk)
            ck = "CAMSV" if jk == "S1" else "CAMTZ"
            if c is None or c["status"] != "measured":
                u = "its control %s_C_ctl_R is %s" % (jk, c["status"] if c else "not a run")
            elif not c["stab"].get("ran") or c["stab"].get("last") != "STABLE":
                u = "the last stability analysis of its control %s_C_ctl_R is not STABLE" % jk
            elif abs(c["E"] - CT[(ck, "R")]["E"]) > TOL_EH:
                u = "its control %s_C_ctl_R does not reproduce the committed energy (differ %.2e Eh)" % (jk, abs(c["E"] - CT[(ck, "R")]["E"]))
        JUSE[jk] = u

    # ---- G47.7
    for r in runs:
        x = RUN[r["run"]]
        if x["status"] != "measured" or r["kind"] == "bare":
            continue
        if r["kind"] == "cond":
            c0, cn, why = ref_of(x)
        else:
            c0, cn = (CT[("CAMSV" if r["level"] == "CAMSV" else "CAMTZ", r["state"])], "the committed control") if r["level"] in ("CAMSV", "CAMTZ") else (None, None)
        if c0 is not None and c0.get("c") and x["c"]["nb"] != c0["c"]["nb"]:
            gates.append("G47.7 %s reports %s basis functions against the %s of %s" % (r["run"], x["c"]["nb"], c0["c"]["nb"], cn))

    if gates:
        say("GATES")
        for g in gates:
            say("   FAILED  %s" % g)
        say("")
        say("A FAILURE OF G47.1 TO G47.7 STOPS THE DIAGNOSTIC. Nothing else is readable.")
        print("\n".join(rep))
        stop("%d gate(s) failed. Nothing is written." % len(gates))

    # ------------------------------------------------------------ helpers
    G = sel["G"]

    def sbcell(sb):
        return "NOT FOUND" if not sb.get("ran") else ("NOT READ" if sb.get("first") is None else "%s/%s" % (sb["first"], sb["last"]))

    def fs(v, fmt="%.6f", none="-"):
        return none if v is None else fmt % v

    def wrapsay(text, lead="   ", width=92):
        for l in textwrap.wrap(text, width):
            say("%s%s" % (lead, l))

    def moved(lw, lw0, state, charges):
        """the atom with the most negative change, that change, its distance to the nearest positive charge, the sum over the 24"""
        dq = [lw[i] - lw0[i] for i in range(24)]
        i = min(range(24), key=lambda j: dq[j])
        pos = [xyz for q, xyz in charges if q > 0]
        dn = min(prep.dist(G[state][i][1:], p) for p in pos) if pos else None
        return prep.ATOMS[i], dq[i], dn, sum(lw[:24]) - sum(lw0[:24])

    def legs(x, c0, charges):
        """the flags and holds of one dianion run against its control"""
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

    def word1(x):
        """the word of one run"""
        if x["noread"]:
            return "NOT READ"
        if x["z"]["flags"]:
            return "FLAGGED"
        if x["status"] != "measured":
            return "NOT MEASURED"
        return "HELD" if x["z"]["held"] else "NOT FLAGGED"

    for r in runs:
        x = RUN[r["run"]]
        if r["kind"] != "cond":
            continue
        c0, cn, why = ref_of(x)
        x["refname"], x["noread"] = cn, JUSE[r["job"]] or why
        x["z"] = legs(x, c0, cmap[r["cond"]]["charges"]) if not x["noread"] else dict(flags=[], held=[], rise=None, ratio=None, gap=None, gap_read=False, mv=None)
        x["E0"] = c0.get("E")
        x["word"] = word1(x)

    def took(x):
        """AT CAMTZG (section 2): (outcome, charge missing from the 24 atoms against the control or None)"""
        mv = x["z"]["mv"]
        if x["noread"] or x["status"] != "measured" or mv is None:
            return "not decided", None
        if mv[3] >= D_LW or x["z"]["flags"]:
            return "TOOK ELECTRONS", mv[3]
        return ("DID NOT TAKE" if x["word"] == "NOT FLAGGED" else "not decided"), mv[3]

    for r in runs:
        if r["kind"] == "cond":
            RUN[r["run"]]["took"] = took(RUN[r["run"]]) if r["level"] == "CAMTZG" else (None, None)

    def cw(name, level, st):
        r = RUN.get(prep.rname(level, name, st))
        return r["word"] if r else None

    def two(name, level):
        """the four-leg word of a condition run at both states at a level"""
        w = [cw(name, level, st) for st in STATES]
        if "FLAGGED" in w:
            return "UNRELIABLE"
        if "NOT MEASURED" in w or "NOT READ" in w or None in w:
            return "NOT CLASSIFIED"
        return "INCONCLUSIVE" if "HELD" in w else "RELIABLE"

    def dde(name, level):
        e = [RUN.get(prep.rname(level, name, st)) for st in STATES]
        if any(r is None or r["status"] != "measured" or r["noread"] or r["E0"] is None for r in e):
            return None
        return ((e[1]["E"] - e[1]["E0"]) - (e[0]["E"] - e[0]["E0"])) * H2K

    def agree(new, ref):
        if new is None or ref is None:
            return "NOT FORMED"
        if abs(ref) < FLOOR:
            return "NOT ASSESSABLE"
        return "WITHIN" if abs(new - ref) < FRAC * abs(ref) else "OUTSIDE"

    # ------------------------------------------------------------ the maps
    pts = sel["points"]
    MAP, mapwhy = {}, []
    for jk, k in (("S1", "CAMSV"), ("T1", "CAMTZ")):
        for sp, tag in (("di", "C_ctl"), ("mo", "C_bare")):
            for st in STATES:
                run = "%s_%s_%s" % (jk, tag, st)
                x = RUN[run]
                ref = (CT if sp == "di" else BM)[(k, st)]["E"]
                v = vread(os.path.join(D, run, "vpot.out"), pts)
                why = None
                if x["status"] != "measured":
                    why = "the run is %s" % x["status"]
                elif abs(x["E"] - ref) > TOL_EH:
                    why = "the run does not reproduce the committed energy (differ %.2e Eh)" % abs(x["E"] - ref)
                elif sp == "di" and x["stab"].get("last") != "STABLE":
                    why = "the run's last stability analysis is not STABLE"
                elif isinstance(v, str):
                    why = v
                if why:
                    mapwhy.append("%s: no map: %s" % (run, why))
                else:
                    MAP[(k, sp, st)] = v
    BARE = dict(((k, st), (BM[(k, st)]["E"] - CT[(k, st)]["E"]) * EV) for k in ("CAMSV", "CAMTZ") for st in STATES)

    def have(k, st):
        return (k, "di", st) in MAP and (k, "mo", st) in MAP

    def vcoef(k, st, p):
        return (MAP[(k, "mo", st)][p] - MAP[(k, "di", st)][p]) * EV

    def s_of(k, st, charges, rows):
        """the smallest S over the positive charges of an arrangement, with the index of that charge"""
        if not have(k, st) or any(p is None for p in rows):
            return None
        vde = BARE[(k, st)] + sum(q * vcoef(k, st, p) for (q, xyz), p in zip(charges, rows))
        best = None
        for i, ((q, xyz), p) in enumerate(zip(charges, rows)):
            if q <= 0:
                continue
            other = sum(qj * BOHR / prep.dist(xyz, xj) for j, ((qj, xj), pj) in enumerate(zip(charges, rows)) if j != i)
            s = vde - EV * (MAP[(k, "mo", st)][p] + other) - RYD * q * q
            if best is None or s < best[0]:
                best = (s, i)
        return best

    def qstar(k, st, p):
        """the size at which S of one positive charge alone at point p is zero"""
        v = vcoef(k, st, p)
        c = BARE[(k, st)] - EV * MAP[(k, "mo", st)][p]
        disc = v * v + 4 * RYD * c
        return None if disc < 0 else (v + disc ** 0.5) / (2 * RYD)

    def rows_of(x):
        return x["ids"] if x["part"] == "Q" else [x["pt"]]

    def s_cell(x, k, st):
        b = s_of(k, st, x["charges"], rows_of(x))
        return None if b is None else b[0]

    def predict(s):
        return "none (no map)" if s is None else ("DID NOT TAKE" if s >= S_BAND else ("TOOK ELECTRONS" if s <= -S_BAND else "none (S within %.1f eV of zero)" % S_BAND))

    # ------------------------------------------------------------ the gates, as they held
    nb = {}
    for k in ("CAMSV", "CAMTZ"):
        nb[k] = CT[(k, "R")]["c"]["nb"]
    say("=" * 96)
    say("G47.1 TO G47.7  THE STOPPING GATES")
    say("=" * 96)
    say("")
    say("   G47.1  all %d files of INPUT_CHECKSUMS.sha256 are the ones listed; every input, geometry and" % nsum)
    say("          point-charge file of the %d runs is what the committed files give." % len(runs))
    say("   G47.2  the %d conditions and %d runs, formed again here from the committed files, are those of" % (len(conds), len(runs)))
    say("          CONDITIONS.tsv, RUNS.tsv and LEVELS.tsv; points.xyz is the %d committed sites, each line" % NS)
    say("          that of d33a/points.xyz, and the %d rungs; the %d job scripts are what those give." % (len(pts) - NS, len(prep.JOBS)))
    say("   G47.3  the %d files of SOURCES.sha256 are the ones listed; no committed file under %s" % (len(srcs), ", ".join(GUARDED)))
    say("          differs from HEAD; the pre-registration and every source are the files at the commit the")
    say("          pre-registration was committed on, %s." % pcm[:8])
    say("   G47.4  the committed controls of CAMSV and CAMTZ are readable and each dianion control's last")
    say("          stability analysis is STABLE; the committed maps of d33a and D43's outputs of the")
    say("          condition this diagnostic shares with it are readable.")
    say("   G47.5  the sentences quoted here are in the documents they are quoted from, and the tolerance of")
    say("          %d %% and its floor of %.1f kcal/mol are %s's." % (round(100 * FRAC), FLOOR, CRIT))
    say("   G47.6  every output that exists prints, as the input it ran, the committed input beside it.")
    say("   G47.7  every dianion run that was read reports the basis functions of its control (CAMSV %s," % nb["CAMSV"])
    say("          CAMTZ %s; at CAMTZG and CAMQZ, those of its own control)." % nb["CAMTZ"])
    say("")
    say("   %d of the %d runs are measured." % (sum(1 for x in RUN.values() if x["status"] == "measured"), len(RUN)))
    say("")

    # ------------------------------------------------------------ R47.1
    say("=" * 96)
    say("R47.1  THE CONTROLS")
    say("=" * 96)
    say("")
    say("   run                    level    status            E, Eh             committed, Eh     differ, Eh   stability      basis   map")
    for r in runs:
        if r["kind"] == "cond":
            continue
        x = RUN[r["run"]]
        k = r["level"]
        ref = None
        if k in ("CAMSV", "CAMTZ"):
            ref = (BM if r["kind"] == "bare" else CT)[(k, r["state"])]["E"]
        mp = "-"
        if r["vpot"]:
            mp = "read" if (k, "mo" if r["kind"] == "bare" else "di", r["state"]) in MAP else "NO MAP"
        say("   %-21s  %-7s  %-16s  %-16s  %-16s  %-11s  %-13s  %-6s  %s" % (
            r["run"], k, x["status"], fs(x.get("E"), "%.8f"), fs(ref, "%.8f", "none: new here"),
            "-" if (x.get("E") is None or ref is None) else "%+.2e" % (x["E"] - ref), sbcell(x["stab"]) if x.get("exists") else "-",
            x["c"]["nb"] if x.get("c") else "-", mp))
    say("")
    for jk, lk, pbs in prep.JOBS:
        if JUSE[jk]:
            wrapsay("JOB %s: %s. Section 2: the runs of that job are then NOT READ." % (jk, JUSE[jk]))
    if not any(JUSE.values()):
        say("   Every job's plain control gave back the committed energy within %.0e Eh (Q1 has none to give back)." % TOL_EH)
    for w in mapwhy:
        wrapsay(w)
    say("")
    for (k, st) in (("CAMSV", "R"), ("CAMSV", "TS"), ("CAMTZ", "R"), ("CAMTZ", "TS")):
        say("   bare detachment energy at %s, %s (committed): %+.4f eV" % (k, st, BARE[(k, st)]))
    say("")
    for g in sel["gctl"]:
        for st in g["states"]:
            x = RUN["CAMTZG_%s_di_%s" % (g["name"], st)]
            if x["status"] == "measured":
                say("   the ghost centre alone, %s at %s: charge-free energy minus the committed CAMTZ control's %+.6f Eh;" % (g["name"], st, x["E"] - CT[("CAMTZ", st)]["E"]))
                say("      sum of the Loewdin charges of the 24 atoms %+.4f (the committed control's %+.4f)" % (sum(x["lw"][:24]), sum(CT[("CAMTZ", st)]["lw"][:24])))
    say("")
    if D43E:
        w = []
        for (cn, st), r0 in sorted(D43E.items()):
            r = RUN.get(prep.rname("CAMSV", cn, st))
            if r and r.get("E") is not None and r0.get("E") is not None:
                w.append((abs(r["E"] - r0["E"]), "%s at %s" % (cn, st)))
        if w:
            say("   D43 ran %s at the working level from the same input and point-charge file. Largest difference" % sorted(D43E)[0][0])
            say("   between this diagnostic's energy and D43's committed one: %.2e Eh (%s)." % max(w))
            say("")

    # ------------------------------------------------------------ R47.2
    say("=" * 96)
    say("R47.2  EVERY RUN WITH A CHARGE: THE LEGS, AND WHERE THE CHARGE WENT")
    say("=" * 96)
    say("")
    say("rise: largest Loewdin charge minus the control's, e. ratio: dipole over the control's. 'most negative")
    say("change': the substrate atom whose Loewdin charge fell most against the control, by how much, and its")
    say("distance to the nearest positive point charge. 'on 24': the sum of the 24 atoms' charges minus the")
    say("control's; at CAMTZG, what is missing from the atoms sits on the functions at the charge.")
    say("")
    say("   run                               status            rise      ratio    gap      stab             most negative change        on 24     word")
    tail = []
    for r in runs:
        if r["kind"] != "cond":
            continue
        x = RUN[r["run"]]
        z = x["z"]
        mv = z["mv"]
        say("   %-32s  %-16s  %-8s  %-7s  %-7s  %-15s  %-26s  %-8s  %s" % (
            r["run"], x["status"], fs(z["rise"], "%+.4f"), fs(z["ratio"], "%.3f"), fs(z["gap"], "%.4f") if z["gap_read"] else "-",
            sbcell(x["stab"]) if x.get("exists") else "-",
            "-" if mv is None else "%-3s %+.3f e, %s" % (mv[0], mv[1], "-" if mv[2] is None else "%.2f A" % mv[2]),
            "-" if mv is None else "%+.3f" % mv[3], x["word"]))
        for f in z["flags"]:
            tail.append("%s: FLAG: %s" % (r["run"], f))
        if x["took"][1] is not None and x["took"][1] >= D_LW:
            tail.append("%s: LEFT THE ATOMS: the 24 atoms' Loewdin charges sum to %+.3f e above the control's (D_LW %.2f e)" % (r["run"], x["took"][1], D_LW))
        for h in z["held"]:
            tail.append("%s: held: %s" % (r["run"], h))
        if x["noread"]:
            tail.append("%s: NOT READ: %s" % (r["run"], x["noread"]))
        elif x["status"] not in ("measured", "DID NOT CONVERGE") and x.get("why"):
            tail.append("%s: %s: %s" % (r["run"], x["status"], "; ".join(x["why"])[:300]))
    say("")
    for l in tail:
        for i, w in enumerate(textwrap.wrap(l, 88)):
            say("   %s%s" % ("NOTE  " if i == 0 else "      ", w))
    if tail:
        say("")
    say("FLAGGED: a leg fires, or the SCF did not converge. NOT FLAGGED: no leg fires, the gap was readable")
    say("and the stability analysis was read. HELD: no leg fires, and the gap was withheld or an analysis was")
    say("not read. NOT READ: the run's job or its control cannot be used. NOT MEASURED: no readable output.")
    say("LEFT THE ATOMS, at CAMTZG only: the 24 atoms' Loewdin charges sum to D_LW or more above the control's.")
    say("It is not a leg and changes no word; R47.6 uses it.")
    say("")

    # ------------------------------------------------------------ R47.3
    def srow(x, k="CAMTZ"):
        return "  ".join("S %s %s" % (st, fs(s_cell(x, k, st), "%+.2f")) for st in STATES)

    say("=" * 96)
    say("R47.3  PART A: THE SIZE LADDER AT SITE %d" % prep.LADDER_SITE)
    say("=" * 96)
    say("")
    say("   size   CAMSV R / TS                four legs       CAMTZ R / TS                four legs       ddE CAMSV   ddE CAMTZ   first order  agreement       S at CAMTZ, eV")
    A = [x for x in conds if x["part"] == "A"]
    for x in A:
        ws, wt = dde(x["name"], "CAMSV"), dde(x["name"], "CAMTZ")
        say("   %.2f   %-12s/ %-12s  %-14s  %-12s/ %-12s  %-14s  %-10s  %-10s  %+9.4f    %-14s  %s" % (
            x["q"], cw(x["name"], "CAMSV", "R"), cw(x["name"], "CAMSV", "TS"), two(x["name"], "CAMSV"),
            cw(x["name"], "CAMTZ", "R"), cw(x["name"], "CAMTZ", "TS"), two(x["name"], "CAMTZ"),
            fs(ws, "%+.4f"), fs(wt, "%+.4f"), x["first"], agree(wt, ws), srow(x)))
    say("")
    okA = [x for x in A if two(x["name"], "CAMTZ") == "RELIABLE"]
    badA = [x for x in A if two(x["name"], "CAMTZ") == "UNRELIABLE"]
    upto, ok = None, True
    for x in A:
        if ok and two(x["name"], "CAMTZ") == "RELIABLE":
            upto = x["q"]
        else:
            ok = False
    say("   RELIABLE at CAMTZ UP TO: %s" % ("%.2f e (that size and every smaller one tested)" % upto if upto else "no size: the smallest size tested is not RELIABLE at CAMTZ"))
    say("   smallest size UNRELIABLE at CAMTZ: %s" % ("%.2f e" % min(x["q"] for x in badA) if badA else "none"))
    say("")
    say("ddE in kcal/mol against the committed controls; 'first order' is a_i q from the committed coefficient.")
    say("Agreement: CAMTZ within %d %% of CAMSV where CAMSV's ddE is %.1f kcal/mol or more in size." % (round(100 * FRAC), FLOOR))
    say("")

    # ------------------------------------------------------------ R47.4, R47.5
    def onestate(title, part, note):
        say("=" * 96)
        say(title)
        say("=" * 96)
        say("")
        say("   condition              where                                     q, e    nearest atom     CAMSV at R     CAMTZ at R     S at CAMTZ, eV")
        for x in [y for y in conds if y["part"] == part]:
            say("   %-21s  %-40s  %+.2f   %-3s %.3f A      %-13s  %-13s  %s" % (
                x["name"], x["where"], x["q"], x["atom"], x["d"], cw(x["name"], "CAMSV", "R"), cw(x["name"], "CAMTZ", "R"), srow(x)))
        say("")
        wrapsay(note, lead="")
        say("")

    onestate("R47.4  PART B: OTHER SURROUNDINGS, AT THE REACTANT", "B",
             "Each site is the outer site with the most negative coefficient among those whose nearest atom is of the class named. Whether a positive charge of the size that is flagged at site %d is flagged here too says whether that site is a special case." % prep.LADDER_SITE)
    onestate("R47.5  PART C: THE WALK OUT, AT THE REACTANT", "C",
             "The same charge on the line from %s through site %d, at larger clearances. NOT FLAGGED at a clearance at CAMTZ says only that def2-TZVPD describes no relocation there; R47.6 asks the same with functions at the charge." % (sel["line"]["atom"], prep.LADDER_SITE))

    # ------------------------------------------------------------ R47.6
    say("=" * 96)
    say("R47.6  PART D1: FUNCTIONS AT THE CHARGE (CAMTZG), AND WHAT S SAID BEFOREHAND")
    say("=" * 96)
    say("")
    say("   run                               S, eV    S predicted                       word          on 24     what the run did   prediction       same point at CAMTZ")
    npred = dict(borne=0, notborne=0, none=0, open=0)
    ntook = [0, 0]
    for r in runs:
        if r["kind"] != "cond" or r["level"] != "CAMTZG":
            continue
        x = RUN[r["run"]]
        c = cmap[r["cond"]]
        s = s_cell(c, "CAMTZ", r["state"])
        p = predict(s)
        did, left = x["took"]
        ntook[1] += 1
        ntook[0] += 1 if did == "TOOK ELECTRONS" else 0
        if p.startswith("none"):
            v = "none made"
            npred["none"] += 1
        elif did == "not decided":
            v = "not decided"
            npred["open"] += 1
        elif did == p:
            v = "borne out"
            npred["borne"] += 1
        else:
            v = "NOT BORNE OUT"
            npred["notborne"] += 1
        say("   %-32s  %-7s  %-32s  %-12s  %-8s  %-17s  %-15s  %s" % (r["run"], fs(s, "%+.2f"), p, x["word"], fs(left, "%+.3f"), did, v, cw(c["name"], "CAMTZ", r["state"]) or "-"))
    say("")
    say("   took electrons at CAMTZG: %d of %d runs." % (ntook[0], ntook[1]))
    say("   predictions: borne out %d, NOT borne out %d, none made %d, not decided %d." % (npred["borne"], npred["notborne"], npred["none"], npred["open"]))
    say("")
    say("TOOK ELECTRONS: a leg fires on a converged SCF, or the run LEFT THE ATOMS ('on 24' at D_LW %.2f e or more)." % D_LW)
    say("DID NOT TAKE: the run is NOT FLAGGED and did not leave the atoms. Otherwise not decided.")
    say("S at or above %+.1f eV predicts DID NOT TAKE, at or below %+.1f eV predicts TOOK ELECTRONS, between them nothing." % (S_BAND, -S_BAND))
    say("S is formed from the maps of CAMTZ, which hold no functions at the charge.")
    say("")

    # ------------------------------------------------------------ R47.7
    say("=" * 96)
    say("R47.7  PART D2: THE TWO DESIGNS THAT PASS IN D44'S LARGER BASIS SET, ONE BASIS STEP UP (CAMQZ)")
    say("=" * 96)
    say("")
    t44 = dict(((r["design"], r["state"]), r) for r in prep.table(os.path.join(S6, "d44.tsv")) if r.get("level_key") == "CAMTZ")

    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    say("   design        CAMQZ R / TS               four legs       ddE working   ddE CAMTZ (D44)   ddE CAMQZ    CAMQZ against working   against CAMTZ    S at CAMTZ, eV")
    Qw = {}
    for x in [y for y in conds if y["part"] == "Q"]:
        dq = dde(x["name"], "CAMQZ")
        w, t = num(t44[(x["name"], "R")].get("ddE_working_kcal")), num(t44[(x["name"], "R")].get("ddE_kcal"))
        Qw[x["name"]] = two(x["name"], "CAMQZ")
        say("   %-12s  %-12s/ %-12s  %-14s  %-12s  %-16s  %-11s  %-22s  %-15s  %s" % (
            x["name"], cw(x["name"], "CAMQZ", "R"), cw(x["name"], "CAMQZ", "TS"), Qw[x["name"]], fs(w, "%+.4f"), fs(t, "%+.4f"), fs(dq, "%+.4f"),
            agree(dq, w), agree(dq, t), srow(x)))
    say("")
    say("ddE at the working level and at CAMTZ are the committed ones of d44.tsv; at CAMQZ it is formed here")
    say("against this diagnostic's own CAMQZ controls. The fifth leg is not read at CAMQZ.")
    say("")

    # ------------------------------------------------------------ R47.8
    say("=" * 96)
    say("R47.8  THE MAPS")
    say("=" * 96)
    say("")
    amat = prep.read_col(os.path.join(S6, prep.AMAT))
    if ("CAMSV", "di", "R") in MAP and ("CAMSV", "di", "TS") in MAP:
        d33 = max(abs(MAP[("CAMSV", "di", st)][i] - V33[st][i]) for st in STATES for i in range(NS))
        aS = [(MAP[("CAMSV", "di", "TS")][i] - MAP[("CAMSV", "di", "R")][i]) * H2K for i in range(NS)]
        say("   The dianion's map at CAMSV against the committed map of d33a, over the %d sites and both states:" % NS)
        say("   largest difference %.2e Eh per e. The coefficient a formed from it against the committed" % d33)
        say("   A_f42436_d33a.tsv: largest difference %.2e kcal/mol per e." % max(abs(aS[i] - amat[i]) for i in range(NS)))
    else:
        say("   The dianion's map at CAMSV is not read at both states; the check against d33a is not formed.")
    say("")
    aT = None
    if ("CAMTZ", "di", "R") in MAP and ("CAMTZ", "di", "TS") in MAP:
        aT = [(MAP[("CAMTZ", "di", "TS")][i] - MAP[("CAMTZ", "di", "R")][i]) * H2K for i in range(NS)]
        df = [aT[i] - amat[i] for i in range(NS)]
        big = [i for i in range(NS) if abs(amat[i]) >= 1.0]
        within = sum(1 for i in big if abs(df[i]) < FRAC * abs(amat[i]))
        flip = sum(1 for i in big if aT[i] * amat[i] < 0)
        say("   THE FIRST-ORDER COEFFICIENT IN THE LARGER BASIS SET. a = V_TS - V_R of the charge-free dianion,")
        say("   kcal/mol per e, at CAMTZ against the committed one at the working level, over the %d sites:" % NS)
        say("      largest difference %.4f, root mean square %.4f" % (max(abs(v) for v in df), (sum(v * v for v in df) / NS) ** 0.5))
        say("      of the %d sites whose committed |a| is 1 or more: %d within %d %%, %d of opposite sign" % (len(big), within, round(100 * FRAC), flip))
        say("      the ten sites of largest committed |a|:")
        say("         site    a working    a CAMTZ     differ")
        for i in sorted(range(NS), key=lambda j: (-abs(amat[j]), j))[:10]:
            say("         %-5d   %+8.4f     %+8.4f    %+7.4f" % (i, amat[i], aT[i], df[i]))
    else:
        say("   The dianion's map at CAMTZ is not read at both states; its coefficients are not formed.")
    say("")
    adm = [s["i"] for s in sel["sites"] if s["qp"] != "none"]
    for k in ("CAMSV", "CAMTZ"):
        for st in STATES:
            if not have(k, st):
                say("   %s, %s: the dianion and monoanion maps are not both read; v and S are not formed." % (k, st))
                continue
            vs = sorted(vcoef(k, st, i) for i in range(NS))
            qs = sorted(q for q in (qstar(k, st, i) for i in adm) if q is not None)
            say("   %s, %s: v = V_monoanion - V_dianion over the %d sites: %.3f to %.3f eV per e (median %.3f)." % (k, st, NS, vs[0], vs[-1], vs[NS // 2]))
            say("      the size at which S of ONE positive charge alone is zero, over the %d sites that admit a positive" % len(adm))
            say("      charge: %.2f to %.2f e (median %.2f); at site %d: %s" % (qs[0], qs[-1], qs[len(qs) // 2], prep.LADDER_SITE, fs(qstar(k, st, prep.LADDER_SITE), "%.2f e")))
    say("")
    say("   THE SIX DESIGNS OF D43, as D44 committed them at CAMTZ, beside S formed here from the maps of CAMTZ")
    say("   and the same reading of where the charge went, made on D44's committed outputs.")
    say("")
    say("   design        st   S, eV    at the charge of     most negative change in D44's output         four legs at CAMTZ (D44)")
    d43c = dict((r["condition"], r) for r in prep.table(os.path.join(S6, "d43/CONDITIONS.tsv")))
    nsort = [0, 0]
    for d in prep.D44_DESIGNS:
        txt = open(os.path.join(S6, "d44/CAMTZ_%s_di_R/pc.pc" % d)).read()
        ch = [(float(f[0]), (float(f[1]), float(f[2]), float(f[3]))) for f in (l.split() for l in txt.split("\n")[1:] if l.strip())]
        ids = [int(v) for v in d43c[d]["sites"].split(",")]
        for st in STATES:
            b = s_of("CAMTZ", st, ch, ids)
            r = D44[(d, st)]
            mv = moved(r["lw"], CT[("CAMTZ", st)]["lw"], st, ch) if r["status"] == "measured" and r.get("lw") else None
            w = t44[(d, st)].get("four_legs")
            say("   %-12s  %-2s   %-7s  %-19s  %-43s  %s" % (
                d, st, fs(b and b[0], "%+.2f"), "-" if b is None else "%+.2f e, site %d" % (ch[b[1]][0], ids[b[1]]),
                "-" if mv is None else "%s %+.3f e, %.2f A from a positive charge" % (mv[0], mv[1], mv[2]), w))
        bR, bT = s_of("CAMTZ", "R", ch, ids), s_of("CAMTZ", "TS", ch, ids)
        if bR and bT and w in ("RELIABLE", "UNRELIABLE"):
            smin = min(bR[0], bT[0])
            if abs(smin) >= S_BAND:
                nsort[0 if (smin > 0) == (w == "RELIABLE") else 1] += 1
    say("")
    say("   Of the six, the sign of the smaller S (where it is %.1f eV or more from zero) goes with D44's word in %d and against it in %d." % (S_BAND, nsort[0], nsort[1]))
    say("   This is a comparison after the fact with results already committed. It tests nothing.")
    say("")

    # ------------------------------------------------------------ what this cannot say
    say("=" * 96)
    say("WHAT THIS CANNOT SAY (section 7 of the pre-registration, in short; its own words govern)")
    say("=" * 96)
    say("")
    for t in (
        "It is a preliminary. It fixes no distance, no size and no rule, adopts no level, rebuilds no grid and runs no search. Section 6 of the pre-registration says what each outcome sends the next step to.",
        "S is an estimate: one electron, the bare charge treated as a nucleus in a complete basis, potentials to first order, nothing relaxed. Where it is borne out here, that is at the few points run.",
        "Parts B and C and most of part D1 are run at the reactant only. A run NOT FLAGGED at the reactant says nothing of the transition state.",
        "The fifth leg is not read. No monoanion with a charge is run.",
        "The ghost centre carries the def2-TZVPD functions of hydrogen. They describe an electron bound to a charge of 0.75 or 1.00 e well, one bound to 0.50 e less well, and one bound to 0.25 e hardly at all; no run at CAMTZG is made at 0.25 e.",
        "CAMQZ is one step. A design that holds there is not thereby shown to hold in every basis set.",
        "One frame, frozen geometries, in vacuo, bare point charges that stand for no chemical group.",
    ):
        wrapsay(t)
        say("")
    nflag = dict((k, sum(1 for r in runs if r["kind"] == "cond" and r["level"] == k and RUN[r["run"]]["word"] == "FLAGGED")) for k in ("CAMSV", "CAMTZ", "CAMTZG", "CAMQZ"))
    ntot = dict((k, sum(1 for r in runs if r["kind"] == "cond" and r["level"] == k)) for k in nflag)
    say("Runs with a charge FLAGGED: CAMSV %d of %d, CAMTZ %d of %d, CAMTZG %d of %d, CAMQZ %d of %d. Ladder at site %d: RELIABLE at CAMTZ up to %s. Took electrons at CAMTZG: %d of %d. Predictions from S at CAMTZG: borne out %d, not borne out %d, none made %d, not decided %d. At CAMQZ: %s." % (
        nflag["CAMSV"], ntot["CAMSV"], nflag["CAMTZ"], ntot["CAMTZ"], nflag["CAMTZG"], ntot["CAMTZG"], nflag["CAMQZ"], ntot["CAMQZ"], prep.LADDER_SITE,
        "%.2f e" % upto if upto else "no size", ntook[0], ntook[1], npred["borne"], npred["notborne"], npred["none"], npred["open"], ", ".join("%s %s" % (n, w) for n, w in sorted(Qw.items()))))

    # ------------------------------------------------------------ the tables
    L = ["run\tjob\tlevel\tkind\tcondition\tpart\tstate\tstatus\tE_Eh\tcontrol\tE_control_Eh\tlowdin_rise_e\tdipole_ratio\tgap_eV\tgap_read\tstab_first\tstab_last\tS2\tTRAH\tlevel_shift\t"
         "flags\theld\tmost_negative_atom\tmost_negative_change_e\tA_to_nearest_positive_charge\tsum_24_minus_control_e\tword\tS_CAMTZ_eV\tat_CAMTZG\tbasis_functions"]

    def sc(r, key):
        sb = r.get("stab") or {}
        return sb.get(key) or ("NOT FOUND" if not sb.get("ran") else "NOT READ")

    for r in runs:
        x = RUN[r["run"]]
        f = [r["run"], r["job"], r["level"], r["kind"], r["cond"] or "", cmap[r["cond"]]["part"] if r["cond"] else "", r["state"], x["status"], fs(x.get("E"), "%.8f", "")]
        if r["kind"] == "cond":
            z, mv = x["z"], x["z"]["mv"]
            f += [x["refname"], fs(x["E0"], "%.8f", ""), fs(z["rise"], "%+.4f", ""), fs(z["ratio"], "%.4f", ""), fs(z["gap"], "%.4f", ""), "yes" if z["gap_read"] else "no"]
        else:
            z, mv = None, None
            f += ["", "", "", "", "", ""]
        ex = x.get("exists")
        f += [sc(x, "first") if ex else "", sc(x, "last") if ex else "", fs(x.get("s2"), "%.6f", ""), ("yes" if x.get("trah") else "no") if ex else "", ("yes" if x.get("lshift") else "no") if ex else ""]
        if z is not None:
            f += ["; ".join(z["flags"]), "; ".join(z["held"] + ([x["noread"]] if x["noread"] else [])), mv[0] if mv else "", fs(mv and mv[1], "%+.4f", ""),
                  fs(mv and mv[2], "%.4f", ""), fs(mv and mv[3], "%+.4f", ""), x["word"], fs(s_cell(cmap[r["cond"]], "CAMTZ", r["state"]), "%+.4f", ""), x["took"][0] or ""]
        else:
            f += [""] * 9
        f.append(str(x["c"]["nb"]) if x.get("c") else "")
        L.append("\t".join(f))
    M = ["point\tkind\tx\ty\tz\tnearest_atom\tnearest_A\ta_committed_working\ta_CAMSV_here\ta_CAMTZ\t" + "\t".join(
        "%s_%s_%s" % (q, k, st) for k in ("CAMSV", "CAMTZ") for st in STATES for q in ("V_dianion_Eh", "V_monoanion_Eh", "v_eV_per_e", "q_at_S_zero_e"))]
    for i, p in enumerate(pts):
        site = sel["sites"][i] if i < NS else None
        f = [str(i), "site" if site else "rung", "%.6f" % p[0], "%.6f" % p[1], "%.6f" % p[2]]
        if site:
            f += [site["atom"], "%.4f" % site["d"], "%+.6f" % amat[i]]
        else:
            c = [x for x in conds if x["part"] == "C" and x["pt"] == i][0]
            f += [c["atom"], "%.4f" % c["d"], ""]
        for k in ("CAMSV", "CAMTZ"):
            ok2 = (k, "di", "R") in MAP and (k, "di", "TS") in MAP
            f.append("%+.6f" % ((MAP[(k, "di", "TS")][i] - MAP[(k, "di", "R")][i]) * H2K) if ok2 else "")
        for k in ("CAMSV", "CAMTZ"):
            for st in STATES:
                vd, vm = MAP.get((k, "di", st)), MAP.get((k, "mo", st))
                f += ["%.10f" % vd[i] if vd else "", "%.10f" % vm[i] if vm else "", "%+.5f" % vcoef(k, st, i) if (vd and vm) else "",
                      fs(qstar(k, st, i) if (vd and vm) else None, "%.4f", "")]
        M.append("\t".join(f))
    C = ["run\tjob\tlevel\tkind\tstate\tstatus\tE_Eh\tcommitted_E_Eh\tdiffer_Eh\tstab_first\tstab_last\tS2\tbasis_functions\tmap"]
    for r in runs:
        if r["kind"] == "cond":
            continue
        x = RUN[r["run"]]
        k = r["level"]
        ref = (BM if r["kind"] == "bare" else CT)[(k, r["state"])]["E"] if k in ("CAMSV", "CAMTZ") else None
        ok2 = x.get("E") is not None
        C.append("\t".join([r["run"], r["job"], k, r["kind"], r["state"], x["status"], "%.8f" % x["E"] if ok2 else "", fs(ref, "%.8f", ""),
                            "%+.3e" % (x["E"] - ref) if (ok2 and ref is not None) else "", sc(x, "first") if x.get("exists") else "", sc(x, "last") if x.get("exists") else "",
                            fs(x.get("s2"), "%.6f", ""), str(x["c"]["nb"]) if x.get("c") else "",
                            ("read" if (k, "mo" if r["kind"] == "bare" else "di", r["state"]) in MAP else "no map") if r["vpot"] else ""]))
    rep[:] = [l.rstrip() for l in rep]
    base = os.path.splitext(a.tsv)[0]
    open(a.tsv, "w").write("\n".join(L) + "\n")
    open(base + "_maps.tsv", "w").write("\n".join(M) + "\n")
    open(base + "_controls.tsv", "w").write("\n".join(C) + "\n")
    open(a.result, "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
