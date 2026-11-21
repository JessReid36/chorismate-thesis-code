#!/usr/bin/env python3
"""DIAGNOSTIC 48's reader. Is the estimate S of DIAGNOSTIC 47 a rule? A test on more sites and on arrangements, with
basis functions at every positive charge.

Every run with a charge is a dianion single point on frame 42436 at CAM-B3LYP-D3BJ/def2-TZVPD that carries, at each
positive charge, the def2-TZVPD functions of a hydrogen atom and no atom (a ghost atom). It is read against the
charge-free control of this diagnostic that carries the same ghost atoms.

THE LEGS, on each run, by the committed classifier cs6_audit.py at its committed thresholds, with the rules D33.3
and AMENDMENT 5 settled and cs6_d33b.py holds: the Loewdin rise, the dipole ratio, the gap (withheld in a run that
engaged TRAH or a level shift, or whose control did) and the LAST stability analysis. A dianion output that says
the SCF did not converge is a flag (D31.3). The fifth leg is not read: no monoanion with a charge is run.

WHAT THE RUN DID, by the reading D47's pre-registration made for runs with functions at the charge:
   LEFT THE ATOMS   the 24 atoms' Loewdin charges sum to D_LW or more above the control's
   TOOK ELECTRONS   the SCF converged and either a leg fires or the run LEFT THE ATOMS
   DID NOT TAKE     the run is NOT FLAGGED and did not leave the atoms
   not decided      anything else
HOW IT TOOK THEM, a mark D48 adds beside that reading: ONTO THE FUNCTIONS if the run TOOK ELECTRONS and LEFT THE
ATOMS; A LEG ONLY if a leg fired and it did not leave the atoms. S estimates the first only, so the counts and the
separation are given by D47's reading and again with the runs A LEG ONLY set aside.

S, formed BEFORE the runs by cs6_d48_prepare.py from the committed maps of D47 and written into CONDITIONS.tsv:
    S_i = VDE_bare + sum_j q_j v_j - EV [ V_monoanion(i) + sum_(j not i) q_j / r_ij ] - RYD q_i^2      (eV)
and S of an arrangement is the smallest S_i over its positive charges. S at or above +S_BAND predicts DID NOT TAKE,
at or below -S_BAND predicts TOOK ELECTRONS, between them nothing. This reader forms S again on its own and holds
each prediction against what the run did.

THE CALIBRATION: one electron, the functions of the ghost atom and a bare charge q. B(q) is minus its energy. S' is
S with RYD q_i^2 replaced by B(q_i) of the level of the runs: what S becomes when its last term is what these
functions and this functional give. S' predicts nothing; it is set beside S.

THE OTHER HALF, on results already committed: the first-order change of the detachment energy, sum_j q_j v_j,
against the change D43 and D44 measured, for the conditions their four legs find RELIABLE. And, for the designs
D44 measured so at def2-TZVPD, S_m: S with the measured detachment energy in place of the first-order one. S_m
predicts nothing.

PART V, THE DETACHMENT ENERGY OF THE PAIRS, MEASURED: for each pair of part B, the dianion and the monoanion with the
pair's charges and no functions at them, from D44's committed inputs. VDE = E(monoanion) - E(dianion) (AMENDMENT 6),
worded by D42's committed word_of(); the dianion is read on the four legs against the committed control. Where that
dianion is NOT FLAGGED, both energies are measured and the monoanion's last stability analysis is STABLE, S_m of the
pair is S with this VDE in place of the first-order one.

The inputs are not taken on trust: the committed cs6_d48_prepare.py is imported by path, its select() is run
again on the committed files, and every table, input and point-charge file are formed again and compared with the
files that ran. The input ORCA echoes at the top of each output must be that file.

usage: cs6_d48.py --res <results-repo> --code <code-repo> --prereg <file> --result <file> --tsv <file>
"""
import argparse, importlib.util, os, re, sys, textwrap

EV = 27.211386245988
RYD = EV / 2.0
S = "phase2.2/cs_stage6"
NEL_MONO, NEL_DI = 117, 118
TOL_EH = 1e-6
PREREG = "STAGE6_DIAGNOSTIC48_20261121.txt"
AMEND9 = "STAGE6_AMENDMENT9_20261120.txt"
D47PRE = "STAGE6_DIAGNOSTIC47_20261121.txt"
D47RES = "STAGE6_DIAGNOSTIC47_RESULT_20261121.txt"
A9_A = ("A design is a SURVIVOR OF THE SEARCH, and nothing more, only if it passes the five legs at the working level AND, at "
        "CAM-B3LYP-D3BJ/def2-TZVPD, is RELIABLE on four legs, passes the fifth leg, and has a ddE that is WITHIN the tolerance of "
        "its ddE at the working level.")
D47_A = ("Every prediction made is borne out: S may be put forward, in a pre-registration of its own, as a rule to be tested on "
         "more sites and on arrangements. Nothing is built on it before that.")
D47_B = "TOOK ELECTRONS the SCF converged and either a leg fires or the run LEFT THE ATOMS"
D47_C = "predictions: borne out 4, NOT borne out 0, none made 0, not decided 0."
GUARDED = ("d31", "d33a", "d37", "d40", "d43", "d44", "d47")
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
    D = os.path.join(S6, "d48")
    paths = dict((k, os.path.join(a.code, "phase2.2", f)) for k, f in (
        ("d40", "cs6_d40.py"), ("d33b", "cs6_d33b.py"), ("d42", "cs6_d42.py"), ("audit", "cs6_audit.py"), ("prep", "cs6_d48_prepare.py")))
    for p in paths.values():
        if not os.path.isfile(p):
            stop("%s is missing; this reader imports it by path." % p)
    d40 = load(paths["d40"], "d48_cs6_d40")
    d33b = load(paths["d33b"], "d48_cs6_d33b")
    audit = load(paths["audit"], "d48_cs6_audit")
    d42 = load(paths["d42"], "d48_cs6_d42")
    word_of, NOTSTABLE = d42.word_of, d42.NOTSTABLE
    prep = load(paths["prep"], "d48_cs6_d48_prepare")
    sha, git = d40.sha, d40.git
    D_LW, R_DIP, GAP_ABS, R_GAP = audit.D_LW, audit.R_DIP, audit.GAP_ABS, audit.R_GAP
    BOHR, S_BAND = prep.BOHR, prep.S_BAND
    if abs(prep.EV - EV) > 0 or abs(prep.RYD - RYD) > 0:
        stop("this reader and cs6_d48_prepare.py do not hold the same constants.")

    rep, gates = [], []
    say = rep.append
    say("STAGE 6, DIAGNOSTIC 48 - RESULT")
    say("Is the estimate S of DIAGNOSTIC 47 a rule? A test on more sites and on arrangements, with basis")
    say("functions at every positive charge.")
    say("")
    say("Read by phase2.2/cs6_d48.py                 sha256 %s" % sha(__file__))
    say("   every output by cs6_d40.py read           sha256 %s" % sha(paths["d40"]))
    say("   stability by cs6_d33b.py stability        sha256 %s" % sha(paths["d33b"]))
    say("   thresholds and legs by cs6_audit.py       sha256 %s" % sha(paths["audit"]))
    say("   the words of binding by cs6_d42.py        sha256 %s" % sha(paths["d42"]))
    say("   the inputs and S by cs6_d48_prepare.py    sha256 %s" % sha(paths["prep"]))
    say("   pre-registration %-23s sha256 %s" % (os.path.basename(a.prereg), sha(a.prereg)))
    say("")
    say("Frame 42436, bare point charges, dianion, CAM-B3LYP-D3BJ/def2-TZVPD. Legs at D_LW %.2f e, R_DIP %.1f," % (D_LW, R_DIP))
    say("GAP_ABS %.1f eV, R_GAP %.1f, stability on the last analysis (AMENDMENT 5). %.12f eV per Eh," % (GAP_ABS, R_GAP, EV))
    say("%.10f A per bohr. S predicts DID NOT TAKE at %+.1f eV or more and TOOK ELECTRONS at %+.1f eV or less." % (BOHR, S_BAND, -S_BAND))
    say("")

    def has(path, sent):
        return os.path.isfile(path) and sent in " ".join(open(path).read().split())

    # ------------------------------------------------------------ G48.2
    sel = None
    try:
        sel = prep.select(S6)
        runs = prep.runs_of(sel)
    except prep.Refused as e:
        gates.append("G48.2 the diagnostic could not be formed again from the committed files: %s" % e)
        runs = []
    except Exception as e:
        gates.append("G48.2 the diagnostic could not be formed again from the committed files: %s: %s" % (type(e).__name__, e))
        runs = []
    conds = sel["conds"] if sel else []
    cmap = dict((x["name"], x) for x in conds)
    if sel and runs:
        for fn, want in (("CONDITIONS.tsv", prep.conditions_text(sel)), ("SITES_A.tsv", prep.sites_text(sel)), ("RUNS.tsv", prep.runs_text(sel)),
                         ("LEVELS.tsv", prep.levels_text(sel))):
            p = os.path.join(D, fn)
            if not os.path.isfile(p) or open(p).read() != want:
                gates.append("G48.2 d48/%s is not what the committed files give" % fn)
        for job in prep.JOBS:
            p = os.path.join(D, "run_%s.pbs" % job[0])
            if not os.path.isfile(p) or open(p).read() != prep.pbs_text(job, sel, runs, prep.HPC):
                gates.append("G48.2 d48/run_%s.pbs is not the job script of %s" % (job[0], job[0]))

    # ------------------------------------------------------------ G48.3
    psrc = os.path.join(D, "SOURCES.sha256")
    srcs = []
    if not os.path.isfile(psrc):
        gates.append("G48.3 d48/SOURCES.sha256 is missing")
    else:
        for l in open(psrc).read().splitlines():
            if l.strip():
                h, f = l.split(None, 1)
                srcs.append(f)
                p = os.path.join(S6, f)
                if not os.path.isfile(p) or sha(p) != h:
                    gates.append("G48.3 %s is not the file d48/SOURCES.sha256 lists" % f)
        if sel and srcs != prep.sources(sel):
            gates.append("G48.3 d48/SOURCES.sha256 does not list the files this diagnostic reads")
    spath = [os.path.normpath("%s/%s" % (S, x)) for x in srcs]
    rc, _ = git(a.res, "diff", "--quiet", "HEAD", "--", *(["%s/%s" % (S, x) for x in GUARDED] + spath))
    if rc is None:
        gates.append("G48.3 could not be checked: git did not run")
    elif rc != 0:
        gates.append("G48.3 a committed file under %s, or one this diagnostic reads, differs from HEAD" % ", ".join(GUARDED))
    rc, pcm = git(a.res, "log", "--diff-filter=A", "--format=%H", "-1", "--", "%s/%s" % (S, PREREG))
    if rc != 0 or not pcm:
        gates.append("G48.3 the commit that added %s could not be found" % PREREG)
        pcm = ""
    else:
        rc1, b1 = git(a.res, "rev-parse", "%s:%s/%s" % (pcm, S, PREREG))
        rc2, b2 = git(a.res, "hash-object", a.prereg)
        if rc1 != 0 or rc2 != 0 or b1 != b2:
            gates.append("G48.3 the pre-registration given is not the committed %s" % PREREG)
        for f, sp in zip(srcs, spath):
            rc1, b1 = git(a.res, "rev-parse", "%s:%s" % (pcm, sp))
            rc2, b2 = git(a.res, "hash-object", os.path.join(S6, f))
            if rc1 != 0 or rc2 != 0 or b1 != b2:
                gates.append("G48.3 %s is not the file at the commit the pre-registration is committed on" % f)

    # ------------------------------------------------------------ G48.5
    for doc, q, what in ((AMEND9, A9_A, "the sentence of A9.4"), (D47PRE, D47_A, "the sentence of D47 on what a borne-out S may be put forward as"),
                         (D47PRE, D47_B, "D47's reading TOOK ELECTRONS"), (D47RES, D47_C, "D47's count of predictions")):
        if not has(os.path.join(S6, doc), q):
            gates.append("G48.5 %s is not in %s" % (what, doc))
    for q in (A9_A, D47_A):
        if not has(a.prereg, q):
            gates.append("G48.5 the words \"%s\" are not in the pre-registration given" % (q[:60] + " ..."))

    # ------------------------------------------------------------ G48.1
    pic = os.path.join(D, "INPUT_CHECKSUMS.sha256")
    nsum = 0
    if not os.path.isfile(pic):
        gates.append("G48.1 d48/INPUT_CHECKSUMS.sha256 is missing")
    else:
        for l in open(pic).read().splitlines():
            if l.strip():
                h, f = l.split(None, 1)
                nsum += 1
                p = os.path.join(D, f)
                if not os.path.isfile(p) or sha(p) != h:
                    gates.append("G48.1 d48/%s is not the file INPUT_CHECKSUMS.sha256 lists" % f)
        if sel and runs and nsum != prep.nsum(sel, runs):
            gates.append("G48.1 INPUT_CHECKSUMS.sha256 lists %d files, not %d" % (nsum, prep.nsum(sel, runs)))
    for r in runs:
        d = os.path.join(D, r["run"])
        x = cmap.get(r["cond"])
        try:
            want = prep.input_text(S6, sel, r["kind"], x, r["state"], r["gname"], r["one"])
        except Exception as e:
            want = None
            gates.append("G48.1 the input of %s could not be formed again: %s" % (r["run"], e))
        pi = os.path.join(d, "job.inp")
        if want is not None and (not os.path.isfile(pi) or open(pi).read() != want):
            gates.append("G48.1 %s/job.inp is not the input the committed files give" % r["run"])
        px = os.path.join(d, "substrate.xyz")
        if r["kind"] not in prep.WITH_XYZ:
            if os.path.exists(px):
                gates.append("G48.1 %s holds a geometry file beside its written-out geometry" % r["run"])
        elif not os.path.isfile(px) or sha(px) != sha(os.path.join(S6, prep.GEOM[r["state"]])):
            gates.append("G48.1 %s/substrate.xyz is not %s" % (r["run"], prep.GEOM[r["state"]]))
        pp = os.path.join(d, "pc.pc")
        wantpc = x["pc"] if x is not None else (prep.one_pc(r["size"]) if r["kind"] == "one" else None)
        if wantpc is not None:
            if not os.path.isfile(pp) or open(pp).read() != wantpc:
                gates.append("G48.1 %s/pc.pc does not hold the charges of %s" % (r["run"], r["cond"] or r["run"]))
        elif os.path.exists(pp):
            gates.append("G48.1 the control %s holds a point-charge file" % r["run"])

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
                gates.append("G48.6 %s/job.out ends normally and does not print the input it ran" % run)
        elif e != want:
            gates.append("G48.6 the input printed at the top of %s/job.out is not the committed %s/job.inp" % (run, run))

    # ---- G48.4: the committed controls and outputs this reader leans on
    CT, BM, G47, P47 = {}, {}, {}, {}
    for st in STATES:
        CT[st] = read_di(os.path.join(S6, "%s_%s" % (prep.CTL, st)), "plain")
        BM[st] = read_di(os.path.join(S6, "%s_%s" % (prep.BARE, st)), "plain", NEL_MONO)
        G47[st] = read_di(os.path.join(S6, prep.D47_GCTL % st), "plain")
        P47[st] = read_di(os.path.join(S6, prep.D47_P050 % st), "pc")
        if CT[st]["status"] != "measured" or not CT[st]["stab"].get("ran") or CT[st]["stab"].get("last") != "STABLE":
            gates.append("G48.4 the committed control %s_%s is not a measured dianion whose last stability analysis is STABLE" % (prep.CTL, st))
        if BM[st]["status"] != "measured":
            gates.append("G48.4 the committed monoanion control %s_%s is not readable" % (prep.BARE, st))
        if G47[st]["status"] != "measured" or P47[st]["status"] != "measured":
            gates.append("G48.4 D47's committed runs with functions at site %d at %s are not readable" % (prep.LADDER_SITE, st))
    t43 = t44 = None
    try:
        t43 = prep.table(os.path.join(S6, "d43.tsv"))
        t44 = prep.table(os.path.join(S6, "d44.tsv"))
        for t, cols, nm in ((t43, ("condition", "state", "sites", "charges", "word", "four_legs", "dVDE_vs_bare_eV", "VDE_eV"), "d43.tsv"),
                            (t44, ("design", "state", "level_key", "word", "four_legs", "dVDE_vs_bare_eV", "VDE_eV"), "d44.tsv")):
            if not t or any(k not in t[0] for k in cols):
                gates.append("G48.4 %s does not hold the columns this reader uses" % nm)
    except Exception as e:
        gates.append("G48.4 d43.tsv or d44.tsv cannot be read: %s" % e)

    # ---- every run of D48
    RUN = {}
    for r in runs:
        d = os.path.join(D, r["run"])
        if r["kind"] == "one":
            x = read_di(d, "pc", 1)
        elif r["kind"] == "pmo":
            x = read_di(d, "pc", NEL_MONO)
        else:
            x = read_di(d, "pc" if r["kind"] in ("cond", "pdi") else "plain")
        check_echo(r["run"], x)
        x.update(r)
        RUN[r["run"]] = x

    def ref_of(x):
        """the control a run with a charge is read against: (its reading, its name, why it cannot be used or None)"""
        c = RUN[x["ref"]]
        why = None
        if c["status"] != "measured":
            why = "its control %s is %s" % (x["ref"], c["status"])
        elif not c["stab"].get("ran") or c["stab"].get("last") is None:
            why = "the stability analysis of its control %s was not read" % x["ref"]
        elif c["stab"]["last"] != "STABLE":
            why = "its control %s is UNSTABLE at the last analysis" % x["ref"]
        return c, x["ref"], why

    # ---- whether a job lets its runs be read: its plain control gives back the committed energy
    JUSE = {}
    for jk, pbs in prep.JOBS:
        u = None
        if jk != "E1":
            c = RUN.get("%s_C_ctl_R" % jk)
            if c is None or c["status"] != "measured":
                u = "its control %s_C_ctl_R is %s" % (jk, c["status"] if c else "not a run")
            elif not c["stab"].get("ran") or c["stab"].get("last") != "STABLE":
                u = "the last stability analysis of its control %s_C_ctl_R is not STABLE" % jk
            elif abs(c["E"] - CT["R"]["E"]) > TOL_EH:
                u = "its control %s_C_ctl_R does not reproduce the committed energy (differ %.2e Eh)" % (jk, abs(c["E"] - CT["R"]["E"]))
        JUSE[jk] = u

    # ---- G48.7
    for r in runs:
        x = RUN[r["run"]]
        if x["status"] != "measured" or r["kind"] == "one" or not x.get("c"):
            continue
        if r["kind"] == "cond":
            c0 = RUN[r["ref"]]
            if c0.get("c") and x["c"]["nb"] != c0["c"]["nb"]:
                gates.append("G48.7 %s reports %s basis functions against the %s of its control %s" % (r["run"], x["c"]["nb"], c0["c"]["nb"], r["ref"]))
        elif r["kind"] in ("ctl", "pdi"):
            if CT["R"].get("c") and x["c"]["nb"] != CT["R"]["c"]["nb"]:
                gates.append("G48.7 %s reports %s basis functions against the %s of the committed control" % (r["run"], x["c"]["nb"], CT["R"]["c"]["nb"]))
        elif CT["R"].get("c") and x["c"]["nb"] <= CT["R"]["c"]["nb"]:
            gates.append("G48.7 the ghost control %s reports %s basis functions, no more than the %s of the committed control" % (r["run"], x["c"]["nb"], CT["R"]["c"]["nb"]))

    if gates:
        say("GATES")
        for g in gates:
            say("   FAILED  %s" % g)
        say("")
        say("A FAILURE OF G48.1 TO G48.7 STOPS THE DIAGNOSTIC. Nothing else is readable.")
        print("\n".join(rep))
        stop("%d gate(s) failed. Nothing is written." % len(gates))

    # ------------------------------------------------------------ helpers
    G = sel["G"]
    M = sel["M"]

    def sbcell(sb):
        return "NOT FOUND" if not sb.get("ran") else ("NOT READ" if sb.get("first") is None else "%s/%s" % (sb["first"], sb["last"]))

    def fs(v, fmt="%.6f", none="-"):
        return none if v is None else fmt % v

    def wrapsay(text, lead="   ", width=92, hang=""):
        for n, l in enumerate(textwrap.wrap(text, width - len(hang))):
            say("%s%s%s" % (lead, hang if n else "", l))

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
        if x["noread"]:
            return "NOT READ"
        if x["z"]["flags"]:
            return "FLAGGED"
        if x["status"] != "measured":
            return "NOT MEASURED"
        return "HELD" if x["z"]["held"] else "NOT FLAGGED"

    def took(x):
        """(what the run did, the charge missing from the 24 atoms against the control or None)"""
        mv = x["z"]["mv"]
        if x["noread"] or x["status"] != "measured" or mv is None:
            return "not decided", None
        if mv[3] >= D_LW or x["z"]["flags"]:
            return "TOOK ELECTRONS", mv[3]
        return ("DID NOT TAKE" if x["word"] == "NOT FLAGGED" else "not decided"), mv[3]

    def s_again(x, st, bind=None):
        """S formed here, apart from cs6_d48_prepare.py: (smallest S_i, the site of that charge, the first-order detachment energy).
        With bind, RYD q^2 is replaced by bind(q): that is S'. None where a charge has no bind."""
        mo, di = M[("mo", st)], M[("di", st)]
        vde = (BM[st]["E"] - CT[st]["E"]) * EV + sum(q * (mo[p] - di[p]) * EV for (q, xyz), p in zip(x["charges"], x["rows"]))
        best = None
        for i, ((q, xyz), p) in enumerate(zip(x["charges"], x["rows"])):
            if q <= 0:
                continue
            other = sum(qj * BOHR / prep.dist(xyz, xj) for j, (qj, xj) in enumerate(x["charges"]) if j != i)
            last = RYD * q * q if bind is None else bind(q)
            if last is None:
                return None, None, vde
            s = vde - EV * (mo[p] + other) - last
            if best is None or s < best[0]:
                best = (s, x["pairs"][i][1])
        return best[0], best[1], vde

    # ---- the calibration
    ONE = {}
    for r in runs:
        if r["kind"] == "one":
            x = RUN[r["run"]]
            ONE[(r["one"], r["size"])] = -x["E"] * EV if x["status"] == "measured" else None

    def bind_cam(q):
        return ONE.get(("CAM", round(100 * q)))

    for r in runs:
        x = RUN[r["run"]]
        if r["kind"] != "cond":
            continue
        cd = cmap[r["cond"]]
        c0, cn, why = ref_of(x)
        x["refname"], x["noread"] = cn, JUSE[r["job"]] or why
        x["z"] = legs(x, c0, cd["charges"]) if not x["noread"] else dict(flags=[], held=[], rise=None, ratio=None, gap=None, gap_read=False, mv=None)
        x["E0"] = c0.get("E")
        x["word"] = word1(x)
        x["took"] = took(x)
        x["how"] = ("ONTO THE FUNCTIONS" if x["took"][1] >= D_LW else "A LEG ONLY") if x["took"][0] == "TOOK ELECTRONS" else ""
        s, site, vde = s_again(cd, r["state"])
        if abs(s - cd["S"][r["state"]]) > 1e-9 or site != cd["names"][r["state"]] or abs(vde - cd["vde"][r["state"]]) > 1e-9:
            stop("S of %s at %s formed here is not the one cs6_d48_prepare.py forms." % (r["cond"], r["state"]))
        x["S"], x["Ssite"], x["vde1"] = s, site, vde
        x["Sp"] = s_again(cd, r["state"], bind_cam)[0]
        x["pred"] = prep.predict(s)
        did = x["took"][0]
        x["verdict"] = "none made" if x["pred"] == "none" else ("not decided" if did == "not decided" else ("borne out" if did == x["pred"] else "NOT BORNE OUT"))
        # the ghost centre that gained most
        ng = len(cd["ghosts"])
        x["gmost"] = None
        if x["status"] == "measured" and not x["noread"] and len(x["lw"]) >= 24 + ng and len(c0["lw"]) >= 24 + ng:
            dg = [x["lw"][24 + k] - c0["lw"][24 + k] for k in range(ng)]
            k = min(range(ng), key=lambda j: dg[j])
            x["gmost"] = ([i for q, i in cd["pairs"] if q > 0][k], dg[k])

    def crun(name, st):
        return RUN[prep.rname(name, st)]

    def first_order(key, st, pairs):
        return sum(q * (M[(key[1], st)][i] - M[(key[0], st)][i]) * EV for q, i in pairs)

    MEAS = {}
    for row in t44:
        if row["level_key"] == "CAMTZ" and row["design"] in sel["designs"] and row.get("four_legs") == "RELIABLE":
            try:
                MEAS[(row["design"], row["state"])] = float(row["dVDE_vs_bare_eV"])
            except ValueError:
                pass
    # ---- part V: the detachment energy of the pairs, measured
    PV = {}
    for k, j in sel["pairsB"]:
        nm = "P_s%04d_p%03d_s%04d_m%03d" % (prep.LADDER_SITE, k, j, k)
        for st in STATES:
            di, mo = RUN[prep.vname(nm, st, "di")], RUN[prep.vname(nm, st, "mo")]
            di["noread"] = JUSE["V1"]
            di["z"] = legs(di, CT[st], cmap[nm]["charges"]) if not di["noread"] else dict(flags=[], held=[], rise=None, ratio=None, gap=None, gap_read=False, mv=None)
            di["word"] = word1(di)
            v = None
            if not JUSE["V1"] and di["status"] == "measured" and mo["status"] == "measured":
                v = (mo["E"] - di["E"]) * EV
            sb = mo.get("stab") or {}
            PV[(nm, st)] = dict(di=di, mo=mo, v=v, w5=word_of(v, sb), vde1=cmap[nm]["vde"][st],
                                ok=(v is not None and di["word"] == "NOT FLAGGED" and bool(sb.get("ran")) and sb.get("last") == "STABLE"))
    for r in runs:
        if r["kind"] == "cond":
            x = RUN[r["run"]]
            m = MEAS.get((r["cond"], r["state"]))
            x["Sm"] = None if m is None else x["S"] - (first_order(prep.TZ, r["state"], cmap[r["cond"]]["pairs"]) - m)
            pv = PV.get((r["cond"], r["state"]))
            if pv is not None:
                if m is not None:
                    stop("%s has a detachment energy measured twice." % r["cond"])
                x["Sm"] = x["S"] - (x["vde1"] - pv["v"]) if pv["ok"] else None

    # ------------------------------------------------------------ the gates, as they held
    say("=" * 96)
    say("G48.1 TO G48.7  THE STOPPING GATES")
    say("=" * 96)
    say("")
    say("   G48.1  all %d files of INPUT_CHECKSUMS.sha256 are the ones listed; every input, geometry and" % nsum)
    say("          point-charge file of the %d runs is what the committed files give." % len(runs))
    say("   G48.2  the %d conditions and %d runs, with S and what it predicts, formed again here from the" % (len(conds), len(runs)))
    say("          committed files, are those of CONDITIONS.tsv, SITES_A.tsv, RUNS.tsv and LEVELS.tsv; the")
    say("          %d job scripts are what those give." % len(prep.JOBS))
    say("   G48.3  the %d files of SOURCES.sha256 are the ones listed; no committed file under %s," % (len(srcs), ", ".join(GUARDED[:4])))
    say("          %s or %s differs from HEAD; the pre-registration and every source are the files at the" % (", ".join(GUARDED[4:-1]), GUARDED[-1]))
    say("          commit the pre-registration was committed on, %s." % pcm[:8])
    say("   G48.4  the committed controls of def2-TZVPD are readable and each dianion control's last stability")
    say("          analysis is STABLE; D47's committed runs with functions at site %d, and d43.tsv and" % prep.LADDER_SITE)
    say("          d44.tsv, are readable.")
    say("   G48.5  the sentences quoted here are in the documents they are quoted from.")
    say("   G48.6  every output that exists prints, as the input it ran, the committed input beside it.")
    say("   G48.7  every run with functions at its charges reports the basis functions of its control; every")
    say("          plain control and every dianion of part V the %s of the committed one; every ghost control" % CT["R"]["c"]["nb"])
    say("          more.")
    say("")
    nrun = [x for x in RUN.values() if x["kind"] != "one"]
    say("   %d of the %d runs on the substrate are measured, and %d of the %d runs of the calibration." % (
        sum(1 for x in nrun if x["status"] == "measured"), len(nrun), sum(1 for x in RUN.values() if x["kind"] == "one" and x["status"] == "measured"),
        sum(1 for x in RUN.values() if x["kind"] == "one")))
    say("")

    # ------------------------------------------------------------ R48.1
    say("=" * 96)
    say("R48.1  THE CONTROLS")
    say("=" * 96)
    say("")
    say("   run                          job  status            E, Eh             minus committed control, Eh   stability      basis   sum on 24")
    for r in runs:
        if r["kind"] not in ("ctl", "gctl"):
            continue
        x = RUN[r["run"]]
        ref = CT[r["state"]]["E"]
        say("   %-27s  %-3s  %-16s  %-16s  %-28s  %-13s  %-6s  %s" % (
            r["run"], r["job"], x["status"], fs(x.get("E"), "%.8f"), "-" if x.get("E") is None else "%+.2e" % (x["E"] - ref),
            sbcell(x["stab"]) if x.get("exists") else "-", x["c"]["nb"] if x.get("c") else "-", fs(sum(x["lw"][:24]) if x.get("lw") else None, "%+.4f")))
    say("")
    for jk, pbs in prep.JOBS:
        if JUSE[jk]:
            wrapsay("JOB %s: %s. Section 2: the runs of that job are then NOT READ." % (jk, JUSE[jk]))
    if not any(JUSE.values()):
        say("   Every job's plain control gave back the committed energy within %.0e Eh (E1 has none)." % TOL_EH)
    say("   The committed control's Loewdin charges of the 24 atoms sum to %+.4f at R and %+.4f at TS." % (sum(CT["R"]["lw"][:24]), sum(CT["TS"]["lw"][:24])))
    say("   Bare detachment energy (committed): %+.4f eV at R, %+.4f eV at TS." % ((BM["R"]["E"] - CT["R"]["E"]) * EV, (BM["TS"]["E"] - CT["TS"]["E"]) * EV))
    say("")

    # ------------------------------------------------------------ R48.2
    say("=" * 96)
    say("R48.2  EVERY RUN WITH A CHARGE: THE LEGS, WHERE THE CHARGE WENT, AND WHAT THE RUN DID")
    say("=" * 96)
    say("")
    say("rise: largest Loewdin charge minus the control's, e. ratio: dipole over the control's. 'on 24': the")
    say("sum of the 24 atoms' charges minus the control's; what is missing from the atoms sits on the")
    say("functions at the charges. 'ghost': the positive charge whose functions gained most, by site, and")
    say("the change of their Loewdin charge.")
    say("")
    say("   run                                    status            rise      ratio    gap      stab             on 24     ghost            word          what the run did")
    tail = []
    for r in runs:
        if r["kind"] != "cond":
            continue
        x = RUN[r["run"]]
        z, mv = x["z"], x["z"]["mv"]
        gm = "-" if x["gmost"] is None else "s%d %+.3f" % x["gmost"]
        say("   %-37s  %-16s  %-8s  %-7s  %-7s  %-15s  %-8s  %-15s  %-12s  %s" % (
            r["run"], x["status"], fs(z["rise"], "%+.4f"), fs(z["ratio"], "%.3f"), fs(z["gap"], "%.4f") if z["gap_read"] else "-", sbcell(x["stab"]) if x.get("exists") else "-",
            fs(mv and mv[3], "%+.3f"), gm, x["word"], x["took"][0]))
        for f in z["flags"]:
            tail.append("%s: FLAG: %s" % (r["run"], f))
        if x["took"][1] is not None and x["took"][1] >= D_LW:
            tail.append("%s: LEFT THE ATOMS: the 24 atoms' Loewdin charges sum to %+.3f e above the control's (D_LW %.2f e)" % (r["run"], x["took"][1], D_LW))
        if x["how"] == "A LEG ONLY":
            tail.append("%s: A LEG ONLY: a leg fired and the run did not leave the atoms (%+.3f e)" % (r["run"], x["took"][1]))
        for h in z["held"]:
            tail.append("%s: held: %s" % (r["run"], h))
        if x["noread"]:
            tail.append("%s: NOT READ: %s" % (r["run"], x["noread"]))
        elif x["status"] not in ("measured", "DID NOT CONVERGE"):
            tail.append("%s: %s: %s" % (r["run"], x["status"], "; ".join(x.get("why") or ["no readable output"])))
    say("")
    for w in tail:
        wrapsay("NOTE  " + w, "   ", 92, "      ")
    say("")

    # ------------------------------------------------------------ R48.6 first formed, printed later
    def cell(x):
        return {"TOOK ELECTRONS": "TOOK", "DID NOT TAKE": "did not", "not decided": "not decided"}[x["took"][0]]

    def pcell(p):
        return {"TOOK ELECTRONS": "TOOK", "DID NOT TAKE": "did not", "none": "nothing"}[p]

    # ------------------------------------------------------------ R48.3
    say("=" * 96)
    say("R48.3  PART A: ONE POSITIVE CHARGE, SITE BY SITE")
    say("=" * 96)
    say("")
    say("S' is S with its last term as the calibration gives it (R48.6); it predicts nothing.")
    for asite in sel["asites"]:
        say("")
        say("   SITE %d (%s): nearest atom %s at %s, %.3f A. S of one charge is zero at %.3f e (R), %.3f e (TS)." % (
            asite["site"], asite["why"], asite["atom"], asite["kst"], asite["d"], asite["zero"]["R"], asite["zero"]["TS"]))
        say("   size   S R      S TS     S' R     S' TS    S predicted R / TS       on 24 R / TS       the run did R / TS            prediction R / TS")
        seq = dict((st, []) for st in STATES)
        for k in asite["sizes"]:
            nm = "A_s%04d_p%03d" % (asite["site"], k)
            xr, xt = crun(nm, "R"), crun(nm, "TS")
            say("   %.2f   %-7s  %-7s  %-7s  %-7s  %-8s / %-11s  %-7s / %-7s  %-12s / %-13s  %s / %s" % (
                k / 100.0, "%+.2f" % xr["S"], "%+.2f" % xt["S"], fs(xr["Sp"], "%+.2f"), fs(xt["Sp"], "%+.2f"), pcell(xr["pred"]), pcell(xt["pred"]),
                fs(xr["took"][1], "%+.3f"), fs(xt["took"][1], "%+.3f"), cell(xr), cell(xt), xr["verdict"], xt["verdict"]))
            for st, x in (("R", xr), ("TS", xt)):
                seq[st].append((k, x["took"][0], x["S"]))
        for st in STATES:
            o = seq[st]
            if any(w == "not decided" for k, w, s in o):
                say("   at %-2s: NOT ORDERED: a run is not decided." % st)
                continue
            no = [k for k, w, s in o if w == "DID NOT TAKE"]
            yes = [k for k, w, s in o if w == "TOOK ELECTRONS"]
            if no and yes and min(yes) < max(no):
                say("   at %-2s: NOT ORDERED: a smaller charge took electrons and a larger one did not." % st)
            elif not yes:
                say("   at %-2s: took electrons at no size run (largest %.2f e, S %+.2f)." % (st, max(no) / 100.0, dict((k, s) for k, w, s in o)[max(no)]))
            elif not no:
                say("   at %-2s: took electrons at every size run (smallest %.2f e, S %+.2f)." % (st, min(yes) / 100.0, dict((k, s) for k, w, s in o)[min(yes)]))
            else:
                sd = dict((k, s) for k, w, s in o)
                say("   at %-2s: did not take up to %.2f e (S %+.2f) and took from %.2f e (S %+.2f). S is zero at %.3f e." % (
                    st, max(no) / 100.0, sd[max(no)], min(yes) / 100.0, sd[min(yes)], asite["zero"][st]))
    say("")

    # ------------------------------------------------------------ R48.4
    say("=" * 96)
    say("R48.4  PART B: ARRANGEMENTS")
    say("=" * 96)
    say("")
    say("S of an arrangement is the smallest S_i over its positive charges; 'S names' is the site of that")
    say("charge. VDE1 is the first-order detachment energy with the charges present, eV. S_m is S with the")
    say("detachment energy measured at def2-TZVPD in place of VDE1: by D44 for the designs its four legs find")
    say("RELIABLE, by part V of this diagnostic for the pairs (R48.9). It predicts nothing.")
    say("")
    say("   condition                  st   S, eV    S', eV   S_m, eV  S names   VDE1     S predicted      on 24     ghost            the run did       prediction")
    for x0 in conds:
        if x0["part"] != "B":
            continue
        for st in STATES:
            x = crun(x0["name"], st)
            say("   %-25s  %-2s   %-7s  %-7s  %-7s  %-8s  %-7s  %-15s  %-8s  %-15s  %-16s  %s" % (
                x0["name"], st, "%+.2f" % x["S"], fs(x["Sp"], "%+.2f"), fs(x["Sm"], "%+.2f"), "s%d" % x["Ssite"], "%+.2f" % x["vde1"],
                {"none": "nothing"}.get(x["pred"], x["pred"]), fs(x["took"][1], "%+.3f"), "-" if x["gmost"] is None else "s%d %+.3f" % x["gmost"], x["took"][0], x["verdict"]))
    say("")
    for x0 in conds:
        if x0["part"] == "B":
            wrapsay("%s: %s." % (x0["name"], x0["where"]))
    say("")
    for k, j in sel["pairsB"]:
        nm = "P_s%04d_p%03d_s%04d_m%03d" % (prep.LADDER_SITE, k, j, k)
        al = "A_s%04d_p%03d" % (prep.LADDER_SITE, k)
        wrapsay("+%.2f e at site %d: alone, the run %s at R and %s at TS; with -%.2f e at site %d, it %s at R and %s at TS." % (
            k / 100.0, prep.LADDER_SITE, cell(crun(al, "R")), cell(crun(al, "TS")), k / 100.0, j, cell(crun(nm, "R")), cell(crun(nm, "TS"))))
    say("")

    # ------------------------------------------------------------ R48.5
    say("=" * 96)
    say("R48.5  THE PREDICTIONS, AND WHETHER S SEPARATES WHAT THE RUNS DID")
    say("=" * 96)
    say("")
    cr = [RUN[r["run"]] for r in runs if r["kind"] == "cond"]

    def counts_and_separation(rows):
        npred = dict((k, sum(1 for x in rows if x["verdict"] == k)) for k in ("borne out", "NOT BORNE OUT", "none made", "not decided"))
        say("   predictions: borne out %d, NOT borne out %d, none made %d, not decided %d." % (npred["borne out"], npred["NOT BORNE OUT"], npred["none made"], npred["not decided"]))
        for x in rows:
            if x["verdict"] == "NOT BORNE OUT":
                say("      NOT BORNE OUT: %s: S %+.2f eV predicted %s; the run %s%s." % (x["run"], x["S"], x["pred"], x["took"][0], ", %s" % x["how"] if x["how"] else ""))
                if x["pred"] == "DID NOT TAKE":
                    x["accounts"] = "NOT KNOWN" if x["Sm"] is None else ("ACCOUNTS FOR IT" if x["Sm"] < S_BAND else "DOES NOT ACCOUNT FOR IT")
                    say("         the first-order detachment energy: %s%s." % (x["accounts"], " (no S_m)" if x["Sm"] is None else ": S_m is %+.2f eV" % x["Sm"]))
        nt = sum(1 for x in rows if x["took"][0] == "TOOK ELECTRONS")
        nn = sum(1 for x in rows if x["took"][0] == "DID NOT TAKE")
        say("   of the %d runs with a charge: took electrons %d, did not take %d, not decided %d." % (len(rows), nt, nn, len(rows) - nt - nn))
        say("")
        sepw = {}
        for key, lab in (("S", "S"), ("Sp", "S'")):
            yes = [(x[key], x["run"]) for x in rows if x["took"][0] == "TOOK ELECTRONS" and x[key] is not None]
            no = [(x[key], x["run"]) for x in rows if x["took"][0] == "DID NOT TAKE" and x[key] is not None]
            if not yes or not no:
                sepw[key] = "NOT FORMED"
                say("   %s: the separation is NOT FORMED: %s." % (lab, "S' is not formed" if (key == "Sp" and not any(x["Sp"] is not None for x in rows)) else "the runs did not do both things"))
                say("")
                continue
            hi, lo = max(yes), min(no)
            say("   %s: the largest at which a run took electrons is %+.2f eV (%s);" % (lab, hi[0], hi[1]))
            say("   %s  the smallest at which a run did not take is %+.2f eV (%s)." % (" " * len(lab), lo[0], lo[1]))
            if hi[0] < lo[0]:
                sepw[key] = "SEPARATES"
                say("   %s SEPARATES what the runs did: any threshold between %+.2f and %+.2f eV gives every one of them." % (lab, hi[0], lo[0]))
            else:
                sepw[key] = "DOES NOT SEPARATE"
                inside = sorted([(v, n, "TOOK") for v, n in yes if v >= lo[0]] + [(v, n, "did not") for v, n in no if v <= hi[0]])
                say("   %s DOES NOT SEPARATE what the runs did. Between %+.2f and %+.2f eV lie:" % (lab, lo[0], hi[0]))
                for v, n, w in inside:
                    say("      %+.2f  %-8s %s" % (v, w, n))
            say("")
        return npred, sepw, nt, nn

    say("BY D47'S READING, every run with a charge:")
    npred, sepw, ntook, nnot = counts_and_separation(cr)
    legonly = [x for x in cr if x["how"] == "A LEG ONLY"]
    say("HOW THE RUNS TOOK THEM: of the %d that took electrons, %d ONTO THE FUNCTIONS and %d A LEG ONLY." % (ntook, ntook - len(legonly), len(legonly)))
    for x in legonly:
        wrapsay("A LEG ONLY: %s: %s; the 24 atoms' charges sum to %+.3f e against the control's; S %+.2f eV predicted %s." % (
            x["run"], "; ".join(x["z"]["flags"]), x["took"][1], x["S"], {"none": "nothing"}.get(x["pred"], x["pred"])), "   ", 92, "   ")
    say("")
    if legonly:
        say("WITH THE %d RUN(S) A LEG ONLY SET ASIDE, as section 6 of the pre-registration reads its outcomes:" % len(legonly))
        npred2, sepw2, ntook2, nnot2 = counts_and_separation([x for x in cr if x["how"] != "A LEG ONLY"])
    else:
        say("No run is A LEG ONLY, so nothing is set aside and the counts above are the ones section 6 of the")
        say("pre-registration reads.")
        say("")
        npred2, sepw2 = npred, sepw
    say("S' is formed only where the calibration is measured at each size the arrangement holds.")
    say("")

    # ------------------------------------------------------------ R48.6
    say("=" * 96)
    say("R48.6  PART C: THE CALIBRATION - WHAT THE FUNCTIONS AT THE CHARGE CAN BIND")
    say("=" * 96)
    say("")
    say("One electron, the def2-TZVPD functions of hydrogen %.2f A from a bare charge q. B is minus its" % prep.GHOST_OFF)
    say("energy, eV. Exact, in a complete basis: 13.6057 q^2. HF is exact for one electron but for the")
    say("functions; CAM is the functional of the runs, without its dispersion term.")
    say("")
    say("   q, e    13.6057 q^2   B, HF         B, CAM        exact minus HF   CAM minus HF   exact minus CAM (what S' adds to S)")
    for k in sel["sizes"]:
        q = k / 100.0
        ex, bh, bc = RYD * q * q, ONE.get(("HF", k)), ONE.get(("CAM", k))
        say("   %.2f    %-11s   %-12s  %-12s  %-15s  %-13s  %s" % (q, "%.4f" % ex, fs(bh, "%.4f", "NOT MEASURED"), fs(bc, "%.4f", "NOT MEASURED"),
                                                              fs(None if bh is None else ex - bh, "%+.4f"), fs(None if (bh is None or bc is None) else bc - bh, "%+.4f"),
                                                              fs(None if bc is None else ex - bc, "%+.4f")))
    say("")
    for r in runs:
        x = RUN[r["run"]]
        if r["kind"] == "one" and x["status"] != "measured":
            wrapsay("NOTE  %s: %s: %s" % (r["run"], x["status"], "; ".join(x.get("why") or ["no readable output"])))
    say("")

    # ------------------------------------------------------------ R48.7
    say("=" * 96)
    say("R48.7  WHAT D47 RAN AND D48 RAN AGAIN")
    say("=" * 96)
    say("")
    for st in STATES:
        g = RUN["CAMTZG_G_s%04d_di_%s" % (prep.LADDER_SITE, st)]
        p = crun("A_s%04d_p050" % prep.LADDER_SITE, st)
        say("   %s: the ghost control of site %d: this energy minus D47's committed one %s Eh; +0.50 e there: %s Eh." % (
            st, prep.LADDER_SITE, "-" if g.get("E") is None else "%+.2e" % (g["E"] - G47[st]["E"]), "-" if p.get("E") is None else "%+.2e" % (p["E"] - P47[st]["E"])))
    say("")
    say("   The inputs and point-charge files of those four runs are D47's, byte for byte.")
    say("")

    # ------------------------------------------------------------ R48.8
    say("=" * 96)
    say("R48.8  THE OTHER HALF, ON RESULTS ALREADY COMMITTED: THE FIRST-ORDER DETACHMENT ENERGY")
    say("=" * 96)
    say("")
    say("The change of the detachment energy that the charges of a condition make, measured by D43 at the")
    say("working level and by D44 at def2-TZVPD (their column dVDE_vs_bare_eV), against sum_j q_j v_j from the")
    say("committed maps of D47 at the same level. Only conditions whose four legs D43 or D44 finds RELIABLE,")
    say("at a state where the value was measured. No run of D48 enters. eV.")
    say("")
    VR = []
    for row in t43:
        try:
            pairs = [(float(q), int(i)) for q, i in zip(row["charges"].split(","), row["sites"].split(","))]
            meas = float(row["dVDE_vs_bare_eV"])
        except ValueError:
            continue
        if row.get("four_legs") != "RELIABLE":
            continue
        VR.append(("working", row["condition"], row["state"], len(pairs), sum(q for q, i in pairs), meas, first_order(prep.SV, row["state"], pairs), float(row["VDE_eV"]), row["word"]))
    for row in t44:
        if row["level_key"] != "CAMTZ" or row["design"] not in sel["designs"] or row.get("four_legs") != "RELIABLE":
            continue
        try:
            meas = float(row["dVDE_vs_bare_eV"])
        except ValueError:
            continue
        pairs = sel["designs"][row["design"]]
        VR.append(("def2-TZVPD", row["design"], row["state"], len(pairs), sum(q for q, i in pairs), meas, first_order(prep.TZ, row["state"], pairs), float(row["VDE_eV"]), row["word"]))
    say("   level        condition       st   charges   sum q, e   measured   first order   first order minus measured   measured VDE   the electron (fifth leg)")
    for lv, cn, st, nq, sq, meas, fo, vde, w in VR:
        say("   %-11s  %-14s  %-2s   %-7d   %+.2f      %+.4f    %+.4f       %+.4f                      %+.4f        %s" % (lv, cn, st, nq, sq, meas, fo, fo - meas, vde, w))
    say("")
    for lv in ("working", "def2-TZVPD"):
        for lab, pick in (("one charge", lambda n: n == 1), ("several charges", lambda n: n > 1)):
            d = [fo - meas for l, cn, st, nq, sq, meas, fo, vde, w in VR if l == lv and pick(nq)]
            if d:
                say("   %-11s %-16s %2d rows; first order minus measured: largest in size %+.4f, root mean square %.4f." % (
                    lv + ",", lab + ":", len(d), max(d, key=abs), (sum(v * v for v in d) / len(d)) ** 0.5))
            else:
                say("   %-11s %-16s no row." % (lv + ",", lab + ":"))
    say("")
    say("This compares an estimate with results committed before this diagnostic was written. It tests")
    say("nothing new and decides nothing.")
    say("")

    # ------------------------------------------------------------ R48.9
    say("=" * 96)
    say("R48.9  PART V: THE DETACHMENT ENERGY OF THE PAIRS, MEASURED")
    say("=" * 96)
    say("")
    say("The dianion and the monoanion with the pair's charges and no functions at them, from D44's committed")
    say("inputs. VDE = E(monoanion) - E(dianion), eV; positive means the dianion holds its last electron")
    say("(AMENDMENT 6, worded by D42's committed word_of()). The dianion is read on the four legs against the")
    say("committed control. VDE1 is the first-order value that S uses. S_m is formed where the dianion is NOT")
    say("FLAGGED, both energies are measured and the monoanion's last stability analysis is STABLE.")
    say("")
    say("   pair                       st   dianion        word          monoanion      <S**2>    stability        VDE       VDE1     VDE1 minus VDE   the electron   S, eV    S_m, eV")
    tail = []
    for (nm, st), pv in PV.items():
        di, mo = pv["di"], pv["mo"]
        say("   %-25s  %-2s   %-13s  %-12s  %-13s  %-8s  %-15s  %-8s  %-7s  %-15s  %-13s  %-7s  %s" % (
            nm, st, di["status"], di["word"], mo["status"], fs(mo.get("s2"), "%.4f"), sbcell(mo["stab"]) if mo.get("exists") else "-", fs(pv["v"], "%+.4f"), "%+.2f" % pv["vde1"],
            fs(None if pv["v"] is None else pv["vde1"] - pv["v"], "%+.4f"), "NOT MEASURED" if pv["w5"] is None else ("ABOVE ZERO *" if pv["w5"] == NOTSTABLE else pv["w5"]),
            "%+.2f" % crun(nm, st)["S"], fs(crun(nm, st)["Sm"], "%+.2f")))
        for f in di["z"]["flags"]:
            tail.append("%s: FLAG: %s" % (di["run"], f))
        for h in di["z"]["held"]:
            tail.append("%s: held: %s" % (di["run"], h))
        if di["noread"]:
            tail.append("%s and %s: NOT READ: %s" % (di["run"], mo["run"], di["noread"]))
        else:
            for x in (di, mo):
                if x["status"] not in ("measured", "DID NOT CONVERGE"):
                    tail.append("%s: %s: %s" % (x["run"], x["status"], "; ".join(x.get("why") or ["no readable output"])))
        if pv["w5"] == NOTSTABLE:
            tail.append("%s: * %s (AMENDMENT 6)" % (mo["run"], NOTSTABLE))
    say("")
    for w in tail:
        wrapsay("NOTE  " + w, "   ", 92, "      ")
    nbound = sum(1 for pv in PV.values() if pv["w5"] == "BOUND")
    say("   Of the %d pair states: the detachment energy is measured at %d; BOUND at %d; S_m formed at %d." % (
        len(PV), sum(1 for pv in PV.values() if pv["v"] is not None), nbound, sum(1 for pv in PV.values() if pv["ok"])))
    say("")

    # ------------------------------------------------------------ what this cannot say
    say("=" * 96)
    say("WHAT THIS CANNOT SAY (section 7 of the pre-registration, in short; its own words govern)")
    say("=" * 96)
    say("")
    for para in (
        "It adopts no rule, changes no word of D44, and opens and closes no size, no distance and no design. Section 6 of the pre-registration says what each outcome sends the next step to.",
        "S is an estimate: one electron, the bare charge treated as a nucleus in a complete basis, potentials to first order, nothing relaxed. Where it is borne out here, that is at the sites, sizes and arrangements run.",
        "The functions at a charge are the def2-TZVPD functions of hydrogen. They bind an electron to a small charge less well than a complete basis would (R48.6). A run that DID NOT TAKE shows that these functions describe no electron leaving; it does not show that none would in a complete basis.",
        "The fifth leg is read for the two pairs only (R48.9), and only without functions at the charges. R48.8 is a comparison on committed results.",
        "Every charge is at an outer site, 3.50 A or more from every atom. Nothing is said of sites closer in.",
        "One frame, frozen geometries, in vacuo, bare point charges that stand for no chemical group.",
    ):
        wrapsay(para)
        say("")
    line = ("Runs with a charge: %d; took electrons %d, did not take %d, not decided %d. Predictions from S: borne out %d, not borne out %d, none made %d, "
            "not decided %d. S %s what the runs did.%s Calibration measured: %d of %d. Detachment energy of the pairs measured: %d of %d." % (
                len(cr), ntook, nnot, len(cr) - ntook - nnot, npred["borne out"], npred["NOT BORNE OUT"], npred["none made"], npred["not decided"],
                {"SEPARATES": "separates", "DOES NOT SEPARATE": "does not separate", "NOT FORMED": "is not held against"}[sepw["S"]],
                "" if not legonly else " With %d run(s) A LEG ONLY set aside: borne out %d, not borne out %d, none made %d, not decided %d; S %s what the runs did." % (
                    len(legonly), npred2["borne out"], npred2["NOT BORNE OUT"], npred2["none made"], npred2["not decided"],
                    {"SEPARATES": "separates", "DOES NOT SEPARATE": "does not separate", "NOT FORMED": "is not held against"}[sepw2["S"]]),
                sum(1 for v in ONE.values() if v is not None), len(ONE), sum(1 for pv in PV.values() if pv["v"] is not None), len(PV)))
    say(line)

    # ------------------------------------------------------------ the tables
    L = ["run\tjob\tcondition\tpart\tstate\tstatus\tE_Eh\tcontrol\tE_control_Eh\tlowdin_rise_e\tdipole_ratio\tgap_eV\tgap_read\tstab_first\tstab_last\tS2\tTRAH\tlevel_shift\t"
         "flags\theld\tmost_negative_atom\tmost_negative_change_e\tsum_24_minus_control_e\tword\tS_eV\tS_prime_eV\tS_m_eV\tS_names_site\tfirst_order_VDE_eV\tS_predicts\twhat_the_run_did\thow_it_took\t"
         "prediction\tghost_most_site\tghost_most_change_e\tbasis_functions"]

    def sc(r, key):
        sb = r.get("stab") or {}
        return sb.get(key) or ("NOT FOUND" if not sb.get("ran") else "NOT READ")

    for r in runs:
        if r["kind"] != "cond":
            continue
        x = RUN[r["run"]]
        z, mv = x["z"], x["z"]["mv"]
        ex = x.get("exists")
        L.append("\t".join([
            r["run"], r["job"], r["cond"], cmap[r["cond"]]["part"], r["state"], x["status"], fs(x.get("E"), "%.8f", ""), x["refname"], fs(x["E0"], "%.8f", ""),
            fs(z["rise"], "%+.4f", ""), fs(z["ratio"], "%.4f", ""), fs(z["gap"], "%.4f", ""), "yes" if z["gap_read"] else "no",
            sc(x, "first") if ex else "", sc(x, "last") if ex else "", fs(x.get("s2"), "%.6f", ""), ("yes" if x.get("trah") else "no") if ex else "", ("yes" if x.get("lshift") else "no") if ex else "",
            "; ".join(z["flags"]), "; ".join(z["held"] + ([x["noread"]] if x["noread"] else [])), mv[0] if mv else "", fs(mv and mv[1], "%+.4f", ""), fs(mv and mv[3], "%+.4f", ""),
            x["word"], "%+.4f" % x["S"], fs(x["Sp"], "%+.4f", ""), fs(x["Sm"], "%+.4f", ""), str(x["Ssite"]), "%+.4f" % x["vde1"], x["pred"], x["took"][0], x["how"], x["verdict"],
            str(x["gmost"][0]) if x["gmost"] else "", fs(x["gmost"][1] if x["gmost"] else None, "%+.4f", ""), str(x["c"]["nb"]) if x.get("c") else ""]))
    C = ["run\tjob\tlevel\tkind\tstate\tsize_e\tstatus\tE_Eh\tminus_committed_control_Eh\tstab_first\tstab_last\tbasis_functions\tsum_24_e\tB_eV"]
    for r in runs:
        if r["kind"] in ("cond", "pdi", "pmo"):
            continue
        x = RUN[r["run"]]
        ok2 = x.get("E") is not None
        one = r["kind"] == "one"
        C.append("\t".join([r["run"], r["job"], r["level"], r["kind"], "" if one else r["state"], "%.2f" % (r["size"] / 100.0) if one else "", x["status"],
                            "%.8f" % x["E"] if ok2 else "", "%+.3e" % (x["E"] - CT[r["state"]]["E"]) if (ok2 and not one) else "",
                            sc(x, "first") if (x.get("exists") and not one) else "", sc(x, "last") if (x.get("exists") and not one) else "",
                            str(x["c"]["nb"]) if x.get("c") else "", fs(sum(x["lw"][:24]) if x.get("lw") else None, "%+.4f", ""),
                            fs(ONE.get((r["one"], r["size"])) if one else None, "%.5f", "")]))
    V = ["level\tcondition\tstate\tcharges\tsum_q_e\tmeasured_dVDE_eV\tfirst_order_dVDE_eV\tfirst_order_minus_measured_eV\tmeasured_VDE_eV\tfifth_leg_word"]
    for lv, cn, st, nq, sq, meas, fo, vde, w in VR:
        V.append("\t".join([lv, cn, st, str(nq), "%+.4f" % sq, "%+.4f" % meas, "%+.4f" % fo, "%+.4f" % (fo - meas), "%+.4f" % vde, w]))
    P = ["condition\tstate\tdi_run\tdi_status\tE_di_Eh\tE_control_Eh\tlowdin_rise_e\tdipole_ratio\tgap_eV\tgap_read\tdi_stab_first\tdi_stab_last\tflags\theld\tword\t"
         "mo_run\tmo_status\tE_mo_Eh\tS2\tmo_stab_first\tmo_stab_last\tVDE_eV\tfifth_leg_word\tfirst_order_VDE_eV\tfirst_order_minus_measured_eV\tS_eV\tS_m_eV"]
    for (nm, st), pv in PV.items():
        di, mo, z = pv["di"], pv["mo"], pv["di"]["z"]
        P.append("\t".join([
            nm, st, di["run"], di["status"], fs(di.get("E"), "%.8f", ""), "%.8f" % CT[st]["E"], fs(z["rise"], "%+.4f", ""), fs(z["ratio"], "%.4f", ""), fs(z["gap"], "%.4f", ""),
            "yes" if z["gap_read"] else "no", sc(di, "first") if di.get("exists") else "", sc(di, "last") if di.get("exists") else "", "; ".join(z["flags"]),
            "; ".join(z["held"] + ([di["noread"]] if di["noread"] else [])), di["word"], mo["run"], mo["status"], fs(mo.get("E"), "%.8f", ""), fs(mo.get("s2"), "%.6f", ""),
            sc(mo, "first") if mo.get("exists") else "", sc(mo, "last") if mo.get("exists") else "", fs(pv["v"], "%+.4f", ""), pv["w5"] or "", "%+.4f" % pv["vde1"],
            fs(None if pv["v"] is None else pv["vde1"] - pv["v"], "%+.4f", ""), "%+.4f" % crun(nm, st)["S"], fs(crun(nm, st)["Sm"], "%+.4f", "")]))
    rep[:] = [l.rstrip() for l in rep]
    base = os.path.splitext(a.tsv)[0]
    open(base + "_pairs.tsv", "w").write("\n".join(P) + "\n")
    open(a.tsv, "w").write("\n".join(L) + "\n")
    open(base + "_controls.tsv", "w").write("\n".join(C) + "\n")
    open(base + "_vde.tsv", "w").write("\n".join(V) + "\n")
    open(a.result, "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
