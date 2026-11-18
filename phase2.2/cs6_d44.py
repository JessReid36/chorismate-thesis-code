#!/usr/bin/env python3
"""DIAGNOSTIC 44's reader: the cross-check of AMENDMENT 7 on the six survivors of D43.

Each of the six designs is one committed point-charge file of D43, run at each of the five
levels of A7.3 at the reactant and at the transition state. Each level has controls with no
charges. Four things are read at each level.

THE FOUR LEGS, on the dianion, against that level's own charge-free control at the same state.
The thresholds are those of the committed classifier cs6_audit.py, read from that file and not
typed in, and its read() gives the gap, the largest Loewdin charge, the dipole and the number
of basis functions. The rules are the ones D33.3 and AMENDMENT 5 settled and cs6_d33b.py holds:
stability on the LAST analysis; the gap withheld in a run that engaged TRAH or a level shift.

  AT THE HYBRID LEVEL, CAMTZ, the dianion run (di) carries the stability analysis and the four
  legs are read on it exactly as DIAGNOSTIC 43 read them at the working level.

  AT A DOUBLE HYBRID the dianion is run twice. di has no analysis: its energy is on the
  restricted solution, the one the record's comparisons of functionals used, and the Loewdin
  rise, the dipole ratio and the gap are read on it. ds is the same dianion with the analysis:
  the search for a lower self-consistent solution. The stability leg is read on ds. If its first
  analysis is STABLE, nothing lower was found. If it is UNSTABLE and ORCA re-converged to a
  solution whose last analysis is STABLE, the three legs are read again on that solution: a
  solution with the charge moved is a flag, and one that keeps the charge in place is reported
  and is not a flag. A last analysis that is UNSTABLE is a flag. The energy of ds is never used
  for a ddE or a detachment energy.

THE FIFTH LEG, by AMENDMENT 6. VDE = E(monoanion) - E(dianion) at the same state with the same
point charges, positive means bound, at 27.211386245988 eV per Eh. The dianion energy is that
of di. The words are D42's, by its committed word_of().

ddE, on the dianion (di), against the same level's control, at 627.5094740631 kcal/mol per Eh.

THE AGREEMENT of that ddE with the design's ddE at the working level, which is recomputed here
from D43's committed outputs and must be the number d43.tsv prints. The tolerance is the one
constant the record carries for method and basis, 15 per cent, with its floor of 3.2 kcal/mol,
both checked against STAGE6_CRITERIA.txt before they are used.

Reading one output. The committed cs6_d40.py is imported by path and its read() is used as it
stands, for dianions and monoanions alike: a run is read only if it terminates normally,
converges its SCF, reports the right number of electrons and says it read as many point
charges as its pc.pc holds; the energy is the last FINAL SINGLE POINT ENERGY; the stability
analysis is read by cs6_d33b.py stability() with ORCA's own verdict sentences as the fallback,
and is never guessed. A dianion output (di) that says the SCF did not converge is a flag (D31.3).

The inputs are not taken on trust. The committed cs6_d44_prepare.py is imported by path; its
select() is run again here on the committed files, and every input is built again and compared
with the file that ran. The input ORCA echoes at the top of each output must be that file.

Nothing is typed in except committed sentences and constants, each named with its source and
each checked against that source before it is used.

usage: cs6_d44.py --res <results-repo> --code <code-repo> --prereg <file> --result <file> --tsv <file>
"""
import argparse, importlib.util, os, re, sys, textwrap

EV = 27.211386245988            # eV per Eh, as D40, D42 and D43 use it
H2K = 627.5094740631
S = "phase2.2/cs_stage6"
NEL_MONO, NEL_DI = 117, 118
TOL_EH = 1e-6                   # a control reproduces the committed charge-free energy within this
FRAC, FLOOR = 0.15, 3.2         # STAGE6_CRITERIA.txt; checked against its sentences below
PREREG = "STAGE6_DIAGNOSTIC44_20261118.txt"
AMEND6 = "STAGE6_AMENDMENT6_20261117.txt"
AMEND7 = "STAGE6_AMENDMENT7_20261117.txt"
CRIT = "STAGE6_CRITERIA.txt"
CRIT_A = "LOSS MUST BE BELOW 15 % OF |ddE(frozen)|."
CRIT_B = ("15 % is not a new constant: it is the fraction decision O6 already requires for the method-and-basis "
          "uncertainty")
CRIT_C = "if |ddE(frozen)| is smaller than 3.2 kcal/mol"
A7_A = "All five are run for every survivor and all five are reported, whichever way they come out. None is dropped."
A7_B = "the tolerance within which a survivor's ddE must agree across the levels"
A7_C = "WHAT A SURVIVOR THAT PASSES THE CROSS-CHECK MAY BE CALLED. That is a dated document of its own and a decision for Jessica."
A7_D = "Agreement of three functionals is not that measurement."
P2_FILE = "vde_probe/PROBE2_RESULT.txt"
P2_SENT = "Any ddE measured at such a condition is not a barrier effect and must not be reported as one."
GUARDED = ("d31", "d33a", "d37", "d37a1", "d40", "d43")
STATES = ("R", "TS")
WORKING = "CAM-B3LYP-D3BJ/def2-SVPD"


def stop(m):
    sys.stderr.write("STOP  %s\n" % m)
    raise SystemExit(1)


def load(path, name):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def echoed(text):
    """The input ORCA prints at the top of its output, as lines; None if there is none."""
    m = re.search(r"INPUT FILE\n=+\nNAME = [^\n]*\n((?:\|\s*\d+>[^\n]*\n)+)", text)
    if not m:
        return None
    ls = [re.sub(r"^\|\s*\d+> ?", "", l) for l in m.group(1).split("\n") if l]
    while ls and (not ls[-1].strip() or "****END OF INPUT****" in ls[-1]):
        ls.pop()
    return ls


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    a.res, a.code, a.prereg = (os.path.abspath(x) for x in (a.res, a.code, a.prereg))
    S6 = os.path.join(a.res, S)
    D = os.path.join(S6, "d44")
    paths = dict((k, os.path.join(a.code, "phase2.2", f)) for k, f in (
        ("d40", "cs6_d40.py"), ("d33b", "cs6_d33b.py"), ("d42", "cs6_d42.py"), ("audit", "cs6_audit.py"), ("prep", "cs6_d44_prepare.py")))
    for p in paths.values():
        if not os.path.isfile(p):
            stop("%s is missing; this reader imports it by path." % p)
    d40 = load(paths["d40"], "d44_cs6_d40")
    d33b = load(paths["d33b"], "d44_cs6_d33b")
    d42 = load(paths["d42"], "d44_cs6_d42")
    audit = load(paths["audit"], "d44_cs6_audit")
    prep = load(paths["prep"], "d44_cs6_d44_prepare")
    sha, git = d40.sha, d40.git
    word_of, NOTSTABLE = d42.word_of, d42.NOTSTABLE
    D_LW, R_DIP, GAP_ABS, R_GAP = audit.D_LW, audit.R_DIP, audit.GAP_ABS, audit.R_GAP

    rep, gates = [], []
    say = rep.append
    say("STAGE 6, DIAGNOSTIC 44 - RESULT")
    say("The cross-check of AMENDMENT 7: the six survivors of D43 at five further levels. Do they")
    say("pass the four legs there, does the dianion stay bound, and does the barrier change agree?")
    say("")
    say("Read by phase2.2/cs6_d44.py                 sha256 %s" % sha(__file__))
    say("   every output by cs6_d40.py read           sha256 %s" % sha(paths["d40"]))
    say("   stability by cs6_d33b.py stability        sha256 %s" % sha(paths["d33b"]))
    say("   thresholds and legs by cs6_audit.py       sha256 %s" % sha(paths["audit"]))
    say("   the words of binding by cs6_d42.py        sha256 %s" % sha(paths["d42"]))
    say("   the inputs by cs6_d44_prepare.py          sha256 %s" % sha(paths["prep"]))
    say("   pre-registration %-23s sha256 %s" % (os.path.basename(a.prereg), sha(a.prereg)))
    for f in (AMEND6, AMEND7, CRIT):
        p = os.path.join(S6, f)
        say("   %-40s sha256 %s" % (f, sha(p) if os.path.isfile(p) else "MISSING"))
    say("")
    say("Frame 42436, bare point charges. Four legs at D_LW %.2f e, R_DIP %.1f, GAP_ABS %.1f eV," % (D_LW, R_DIP, GAP_ABS))
    say("R_GAP %.1f, against the charge-free control of the same level and state; stability on the" % R_GAP)
    say("last analysis (AMENDMENT 5). Fifth leg by AMENDMENT 6: VDE = E(monoanion) - E(dianion),")
    say("same state, same charges, %.12f eV per Eh. Agreement: a ddE within %d %% of the" % (EV, round(100 * FRAC)))
    say("design's ddE at the working level, %s." % WORKING)
    say("")

    def has(path, sent):
        return os.path.isfile(path) and sent in " ".join(open(path).read().split())

    # ------------------------------------------------------------ G44.2: the levels and the designs, again
    sel = None
    try:
        sel = prep.select(S6)
    except prep.Refused as e:
        gates.append("G44.2 the levels and the designs could not be read again from the committed files: %s" % e)
    except Exception as e:                                   # a file that does not parse is a gate, not a crash
        gates.append("G44.2 the levels and the designs could not be read again from the committed files: %s: %s" % (type(e).__name__, e))
    levels = sel["levels"] if sel else []
    designs = sel["designs"] if sel else []
    runs = []
    if sel:
        try:
            runs = prep.runs_of(sel)
        except prep.Refused as e:
            gates.append("G44.2 %s" % e)
        for fn, want in (("CONDITIONS.tsv", prep.conditions_text(designs)), ("LEVELS.tsv", prep.levels_text(levels))):
            p = os.path.join(D, fn)
            if not os.path.isfile(p) or open(p).read() != want:
                gates.append("G44.2 d44/%s is not what the committed files give" % fn)
        for lv in levels:
            p = os.path.join(D, "run_%s.pbs" % lv["key"])
            if not os.path.isfile(p) or open(p).read() != prep.pbs_text(lv, runs, prep.HPC):
                gates.append("G44.2 d44/run_%s.pbs is not the job script of %s" % (lv["key"], lv["name"]))
    lmap = dict((l["key"], l) for l in levels)
    dmap = dict((d["name"], d) for d in designs)

    # ------------------------------------------------------------ G44.3: the committed files
    psrc = os.path.join(D, "SOURCES.sha256")
    srcs = []
    if not os.path.isfile(psrc):
        gates.append("G44.3 d44/SOURCES.sha256 is missing")
    else:
        for l in open(psrc).read().splitlines():
            if l.strip():
                h, f = l.split(None, 1)
                srcs.append(f)
                p = os.path.join(S6, f)
                if not os.path.isfile(p) or sha(p) != h:
                    gates.append("G44.3 %s is not the file d44/SOURCES.sha256 lists" % f)
        if sel and srcs != prep.sources([d["name"] for d in designs]):
            gates.append("G44.3 d44/SOURCES.sha256 does not list the files this diagnostic reads")
    rc, _ = git(a.res, "diff", "--quiet", "HEAD", "--", *(["%s/%s" % (S, x) for x in GUARDED] + ["%s/%s" % (S, x) for x in srcs]))
    if rc is None:
        gates.append("G44.3 could not be checked: git did not run")
    elif rc != 0:
        gates.append("G44.3 a committed file under %s, or one this diagnostic reads, differs from HEAD" % ", ".join(GUARDED))
    rc, pcm = git(a.res, "log", "--diff-filter=A", "--format=%H", "-1", "--", "%s/%s" % (S, PREREG))
    if rc != 0 or not pcm:
        gates.append("G44.3 the commit that added %s could not be found" % PREREG)
        pcm = ""
    else:
        rc1, b1 = git(a.res, "rev-parse", "%s:%s/%s" % (pcm, S, PREREG))
        rc2, b2 = git(a.res, "hash-object", a.prereg)
        if rc1 != 0 or rc2 != 0 or b1 != b2:
            gates.append("G44.3 the pre-registration given is not the committed %s" % PREREG)
        for f in srcs:
            rc1, b1 = git(a.res, "rev-parse", "%s:%s/%s" % (pcm, S, f))
            rc2, b2 = git(a.res, "hash-object", os.path.join(S6, f))
            if rc1 != 0 or rc2 != 0 or b1 != b2:
                gates.append("G44.3 %s is not the file at the commit the pre-registration is committed on" % f)

    # ------------------------------------------------------------ G44.5: the sentences and constants this result uses
    pcrit = os.path.join(S6, CRIT)
    for q in (CRIT_A, CRIT_B, CRIT_C):
        if not has(pcrit, q):
            gates.append("G44.5 the words \"%s\" are not in %s" % (q, CRIT))
    if "%d %%" % round(100 * FRAC) not in CRIT_A or "%.1f kcal/mol" % FLOOR not in CRIT_C:
        gates.append("G44.5 the tolerance and the floor this reader uses are not the ones its quotations of %s state" % CRIT)
    for q in (A7_A, A7_B, A7_C, A7_D):
        if not has(os.path.join(S6, AMEND7), q):
            gates.append("G44.5 the words \"%s\" are not in %s" % (q, AMEND7))
    if not has(os.path.join(S6, P2_FILE), P2_SENT):
        gates.append("G44.5 Probe 2's sentence is not in %s" % P2_FILE)
    for q in (CRIT_A, A7_A, A7_C):
        if not has(a.prereg, q):
            gates.append("G44.5 the words \"%s\" are not in the pre-registration given" % q)

    # ------------------------------------------------------------ G44.1: every committed input
    pic = os.path.join(D, "INPUT_CHECKSUMS.sha256")
    nsum = 0
    if not os.path.isfile(pic):
        gates.append("G44.1 d44/INPUT_CHECKSUMS.sha256 is missing")
    else:
        for l in open(pic).read().splitlines():
            if l.strip():
                h, f = l.split(None, 1)
                nsum += 1
                p = os.path.join(D, f)
                if not os.path.isfile(p) or sha(p) != h:
                    gates.append("G44.1 d44/%s is not the file INPUT_CHECKSUMS.sha256 lists" % f)
        if nsum != prep.NSUM:
            gates.append("G44.1 INPUT_CHECKSUMS.sha256 lists %d files, not %d" % (nsum, prep.NSUM))
    nfile = 0
    for run, lk, dn, kind, st, cm, stab in runs:
        d = os.path.join(D, run)
        lv = lmap[lk]
        x, y = os.path.join(d, "substrate.xyz"), os.path.join(S6, prep.GEOM[st])
        if not (os.path.isfile(x) and os.path.isfile(y)) or sha(x) != sha(y):
            gates.append("G44.1 %s/substrate.xyz is not %s" % (run, prep.GEOM[st]))
        nfile += 1
        pi = os.path.join(d, "job.inp")
        src = os.path.join(S6, "%s_%s" % (lv["src"], st), "job.inp")
        if not os.path.isfile(pi) or open(pi).read() != prep.build_input(lv, cm, dn is not None, stab):
            gates.append("G44.1 %s/job.inp is not the input of %s built again from the committed files" % (run, lv["name"]))
        elif open(pi).read().split("\n")[0] != open(src).read().split("\n")[0]:
            gates.append("G44.1 the method line of %s is not line 1 of %s_%s/job.inp" % (run, lv["src"], st))
        elif dn is None and not stab and sha(pi) != sha(src):
            gates.append("G44.1 the control %s is not %s_%s/job.inp byte for byte" % (run, lv["src"], st))
        nfile += 1
        if dn is not None:
            pp = os.path.join(d, "pc.pc")
            if not os.path.isfile(pp) or sha(pp) != dmap[dn]["pcsha"] or sha(pp) != sha(os.path.join(S6, dmap[dn]["pc"])):
                gates.append("G44.1 %s/pc.pc is not D43's committed point-charge file of %s" % (run, dn))
            nfile += 1
        elif os.path.exists(os.path.join(d, "pc.pc")):
            gates.append("G44.1 the control %s holds a point-charge file" % run)

    # ------------------------------------------------------------ reading one output
    def read_di(d, kind):
        """A dianion output: cs6_d40.read, then the committed classifier's read on the same text."""
        r = d40.read(d, kind, NEL_DI, d33b)
        r["status"], r["c"] = "measured", None
        if not r.get("exists"):
            r["status"] = "NOT MEASURED"
            return r
        if r["E"] is None:
            r["status"] = "DID NOT CONVERGE" if "SCF NOT CONVERGED" in r["text"] else "NOT MEASURED"
            return r
        c = audit.read(r["text"])
        if c is None or c.get("maxlw") is None or c.get("dip") is None or c.get("nb") is None:
            r["status"], r["E"] = "NOT MEASURED", None
            r["why"] = ["the committed classifier does not read a Loewdin charge, a dipole and a basis-function count from it"]
            return r
        r["c"] = c
        return r

    def read_mo(d, kind):
        r = d40.read(d, kind, NEL_MONO, d33b)
        r["status"] = "measured" if r.get("E") is not None else "NOT MEASURED"
        return r

    def check_echo(run, r):
        """G44.6: an output that exists must be the output of the committed input beside it."""
        if not r.get("exists"):
            return
        e = echoed(r["text"])
        pin = os.path.join(D, run, "job.inp")
        if not os.path.isfile(pin):
            return                                           # G44.1 has already failed on it
        want = open(pin).read().split("\n")
        while want and not want[-1].strip():
            want.pop()
        if e is None:
            if r.get("normal"):
                gates.append("G44.6 %s/job.out ends normally and does not print the input it ran" % run)
        elif e != want:
            gates.append("G44.6 the input printed at the top of %s/job.out is not the committed %s/job.inp" % (run, run))

    # ------------------------------------------------------------ G44.4: the working level, from D43's committed outputs
    W = {}
    CTw, BMw = {}, {}
    for st in STATES:
        CTw[st] = read_di(os.path.join(S6, "d31", "ctl_%s" % st), "plain")
        BMw[st] = read_mo(os.path.join(S6, "d40", "B_bare_%s" % st), "plain")
        if CTw[st]["E"] is None or BMw[st]["E"] is None:
            gates.append("G44.4 the committed controls of the working level at %s are not readable" % st)
    t43 = {}
    p43 = os.path.join(S6, prep.D43_TSV)
    if os.path.isfile(p43):
        for r in prep.table(p43):
            t43[(r.get("condition"), r.get("state"))] = r
    if not any(g.startswith("G44.4") for g in gates):
        for d in designs:
            e, v = {}, {}
            for st in STATES:
                di = read_di(os.path.join(S6, "d43", "%s_di_%s" % (d["name"], st)), "pc")
                mo = read_mo(os.path.join(S6, "d43", "%s_mo_%s" % (d["name"], st)), "pc")
                if di["E"] is None or mo["E"] is None:
                    gates.append("G44.4 D43's committed output of %s at %s is not readable" % (d["name"], st))
                    continue
                e[st] = di["E"]
                v[st] = (mo["E"] - di["E"]) * EV
                if "%+.4f" % v[st] != t43.get((d["name"], st), {}).get("VDE_eV"):
                    gates.append("G44.4 the working-level detachment energy of %s at %s, recomputed from D43's outputs, is not the %s eV d43.tsv prints"
                                 % (d["name"], st, t43.get((d["name"], st), {}).get("VDE_eV")))
            if len(e) == 2:
                dde = ((e["TS"] - CTw["TS"]["E"]) - (e["R"] - CTw["R"]["E"])) * H2K
                if any("%+.4f" % dde != t43.get((d["name"], st), {}).get("ddE_kcal") for st in STATES):
                    gates.append("G44.4 the working-level ddE of %s, recomputed from D43's outputs, is not the %s kcal/mol d43.tsv prints"
                                 % (d["name"], t43.get((d["name"], "R"), {}).get("ddE_kcal")))
                W[d["name"]] = dict(dde=dde, v=v)
    barw = None if any(CTw[st]["E"] is None for st in STATES) else (CTw["TS"]["E"] - CTw["R"]["E"]) * H2K
    barew = dict((st, None if (BMw[st]["E"] is None or CTw[st]["E"] is None) else (BMw[st]["E"] - CTw[st]["E"]) * EV) for st in STATES)

    # ------------------------------------------------------------ the committed charge-free points of each level (D37, D37 Amendment 1)
    OLD = {}
    for lv in levels:
        for st in STATES:
            r = read_di(os.path.join(S6, "%s_%s" % (lv["src"], st)), "plain")
            OLD[(lv["key"], st)] = r
            if r["E"] is None:
                gates.append("G44.4 the committed charge-free output %s_%s is not readable: %s" % (lv["src"], st, "; ".join(r.get("why") or [])))

    # ------------------------------------------------------------ read every run of D44
    CT, CS, BM, X = {}, {}, {}, {}
    for lv in levels:
        k = lv["key"]
        for st in STATES:
            run = prep.cname(k, "di", st)
            CT[(k, st)] = read_di(os.path.join(D, run), "plain")
            check_echo(run, CT[(k, st)])
            if lv["dh"]:
                run = prep.cname(k, "ds", st)
                CS[(k, st)] = read_di(os.path.join(D, run), "plain")
                check_echo(run, CS[(k, st)])
            run = prep.cname(k, "mo", st)
            BM[(k, st)] = read_mo(os.path.join(D, run), "plain")
            check_echo(run, BM[(k, st)])
        for d in designs:
            for st in STATES:
                x = dict(rdi=prep.rname(k, d["name"], "di", st), rmo=prep.rname(k, d["name"], "mo", st), rds=None, ds=None)
                x["di"] = read_di(os.path.join(D, x["rdi"]), "pc")
                x["mo"] = read_mo(os.path.join(D, x["rmo"]), "pc")
                check_echo(x["rdi"], x["di"])
                check_echo(x["rmo"], x["mo"])
                if lv["dh"]:
                    x["rds"] = prep.rname(k, d["name"], "ds", st)
                    x["ds"] = read_di(os.path.join(D, x["rds"]), "pc")
                    check_echo(x["rds"], x["ds"])
                X[(k, d["name"], st)] = x
        # G44.7: the basis of the level
        for st in STATES:
            old = OLD.get((k, st))
            if old is None or old["c"] is None:
                continue
            every = [(prep.cname(k, "di", st), CT[(k, st)])] + [(X[(k, d["name"], st)]["rdi"], X[(k, d["name"], st)]["di"]) for d in designs]
            if lv["dh"]:
                every += [(prep.cname(k, "ds", st), CS[(k, st)])] + [(X[(k, d["name"], st)]["rds"], X[(k, d["name"], st)]["ds"]) for d in designs]
            for nm, r in every:
                if r["status"] == "measured" and r["c"]["nb"] != old["c"]["nb"]:
                    gates.append("G44.7 %s reports %s basis functions against the %s of %s_%s" % (nm, r["c"]["nb"], old["c"]["nb"], lv["src"], st))

    if gates:
        say("GATES")
        for g in gates:
            say("   FAILED  %s" % g)
        say("")
        say("A FAILURE OF G44.1 TO G44.7 STOPS THE DIAGNOSTIC. Nothing else is readable.")
        say("The reader stops in the same way where a committed file it must quote or read is not")
        say("as the pre-registration names it.")
        print("\n".join(rep))
        stop("%d gate(s) failed. Nothing is written." % len(gates))

    # ------------------------------------------------------------ helpers
    def sbcell(sb):
        return "NOT FOUND" if not sb.get("ran") else ("NOT READ" if sb.get("first") is None else "%s/%s" % (sb["first"], sb["last"]))

    def fs(v, fmt="%.6f", none="-"):
        return none if v is None else fmt % v

    def wrapsay(text, lead="   ", width=76):
        for l in textwrap.wrap(text, width):
            say("%s%s" % (lead, l))

    def notes(tail):
        for l in tail:
            for i, w in enumerate(textwrap.wrap(l, 80)):
                say("   %s%s" % ("NOTE  " if i == 0 else "      ", w))

    def negs(sb):
        return "" if sb.get("nneg") is None else " (%d negative, lowest %+.6f Eh)" % (sb["nneg"], sb["lo"])

    def three(r, c0, flags, held, what=""):
        """The Loewdin rise, the dipole ratio and the gap of one dianion solution against a control's.
        Appends to flags and held; returns (rise, ratio, gap, gap_read)."""
        g, kk = r["c"], c0["c"]
        rise = g["maxlw"] - kk["maxlw"]
        ratio = g["dip"] / kk["dip"] if kk["dip"] else None
        if rise >= D_LW:
            flags.append("%srelocation: Loewdin rise %+.4f e" % (what, rise))
        if ratio is None:
            held.append("%sthe dipole ratio could not be formed" % what)
        elif ratio >= R_DIP:
            flags.append("%sdipole %.2fx control" % (what, ratio))
        gap, gap_read = g.get("gap"), False
        if r["trah"] or r["lshift"] or c0["trah"] or c0["lshift"] or g.get("gap") is None or kk.get("gap") is None:
            held.append("%sgap withheld" % what)
        else:
            gap_read = True
            if g["gap"] < GAP_ABS and g["gap"] < R_GAP * kk["gap"]:
                flags.append("%sgap %.4f eV, %.2f of control" % (what, g["gap"], g["gap"] / kk["gap"]))
        return rise, ratio, gap, gap_read

    # ------------------------------------------------------------ whether the control of a level can be used
    USE = {}
    for lv in levels:
        k = lv["key"]
        for st in STATES:
            c = CT[(k, st)]
            u = None
            if c["status"] != "measured":
                u = "its control at %s is %s" % (st, c["status"])
            elif lv["dh"]:
                if abs(c["E"] - OLD[(k, st)]["E"]) > TOL_EH:
                    u = "its control at %s does not reproduce the committed energy of the same input (differ %.2e Eh)" % (st, abs(c["E"] - OLD[(k, st)]["E"]))
            elif not c["stab"].get("ran") or c["stab"].get("last") is None:
                u = "the stability analysis of its control at %s was not read" % st
            elif c["stab"]["last"] != "STABLE":
                u = "its control at %s is UNSTABLE at the last analysis" % st
            USE[(k, st)] = u

    # ------------------------------------------------------------ the legs of every design at every level
    for lv in levels:
        k = lv["key"]
        for d in designs:
            for st in STATES:
                x = X[(k, d["name"], st)]
                di, mo, c0 = x["di"], x["mo"], CT[(k, st)]
                x.update(flags=[], held=[], notes=[], rise=None, ratio=None, gap=None, gap_read=False, noctl=USE[(k, st)],
                         srise=None, sratio=None, sgap=None, sdiff=None, search="")
                if di["status"] == "DID NOT CONVERGE":
                    x["flags"].append("did not converge")
                elif di["status"] == "measured" and x["noctl"] is None:
                    x["rise"], x["ratio"], x["gap"], x["gap_read"] = three(di, c0, x["flags"], x["held"])
                if not lv["dh"]:
                    if di["status"] == "measured":
                        sb = di["stab"]
                        if not sb["ran"] or sb["first"] is None or sb["last"] is None:
                            x["held"].append("stability analysis not read")
                        elif sb["last"] == "UNSTABLE":
                            x["flags"].append("unstable at the last analysis")
                        elif sb["first"] == "UNSTABLE":
                            x["notes"].append("the first solution was unstable%s and was corrected; the last analysis is STABLE (AMENDMENT 5: not a flag)" % negs(sb))
                else:
                    ds = x["ds"]
                    if ds["status"] != "measured":
                        x["search"] = "NOT READ"
                        x["held"].append("the stability search %s" % ("did not converge" if ds["status"] == "DID NOT CONVERGE" else "is NOT MEASURED"))
                    else:
                        sb = ds["stab"]
                        if di["E"] is not None:
                            x["sdiff"] = ds["E"] - di["E"]
                        if not sb["ran"] or sb["first"] is None or sb["last"] is None:
                            x["search"] = "NOT READ"
                            x["held"].append("stability analysis not read")
                        elif sb["last"] == "UNSTABLE":
                            x["search"] = "UNSTABLE AT THE END"
                            x["flags"].append("unstable at the last analysis of the stability search")
                        elif sb["first"] == "STABLE":
                            x["search"] = "NOTHING LOWER"
                            if x["sdiff"] is not None and abs(x["sdiff"]) > TOL_EH:
                                x["held"].append("the stability search did not start from the solution whose energy is used (differ %.2e Eh)" % abs(x["sdiff"]))
                        elif x["noctl"] is None:
                            f2, h2 = [], []
                            x["srise"], x["sratio"], x["sgap"], _ = three(ds, c0, f2, h2, "the solution the search ended on: ")
                            x["flags"] += f2
                            x["held"] += h2
                            x["search"] = "CHARGE MOVED" if f2 else "LOWER, CHARGE IN PLACE"
                            if not f2:
                                x["notes"].append("the restricted solution is unstable%s; ORCA re-converged to a solution whose last analysis is STABLE and which keeps the charge in place (rise %+.4f e, dipole %sx, <S**2> %s); not a flag, and the energy used is the restricted one" % (
                                    negs(sb), x["srise"], fs(x["sratio"], "%.2f"), fs(ds.get("s2"))))
                        else:
                            x["search"] = "LOWER, NOT READ"
                x["v"] = None if (di["E"] is None or mo["E"] is None) else (mo["E"] - di["E"]) * EV
                x["word"] = word_of(x["v"], mo.get("stab") or {})

    def verdict(k, d):
        xs = [X[(k, d["name"], st)] for st in STATES]
        if any(x["flags"] for x in xs):
            return "UNRELIABLE"
        if any(x["di"]["status"] != "measured" or x["noctl"] for x in xs):
            return "NOT CLASSIFIED"
        if any(x["held"] for x in xs):
            return "INCONCLUSIVE"
        return "RELIABLE"

    def fifth(k, d):
        w = [X[(k, d["name"], st)]["word"] for st in STATES]
        if w == ["BOUND", "BOUND"]:
            return "PASSED"
        if "METASTABLE" in w:
            return "NOT PASSED: METASTABLE"
        if None in w:
            return "NOT PASSED: NOT MEASURED"
        return "NOT PASSED: NOT STABLE"

    def cellw(k, d, st):
        w = X[(k, d["name"], st)]["word"]
        return "NOT MEASURED" if w is None else ("ABOVE ZERO *" if w == NOTSTABLE else w)

    Rz = {}
    for lv in levels:
        k = lv["key"]
        for d in designs:
            z = dict(verdict=verdict(k, d), fifth=fifth(k, d), dde=None, dev=None, frac=None)
            e = dict((st, X[(k, d["name"], st)]["di"]["E"]) for st in STATES)
            w = W[d["name"]]["dde"]
            if None in e.values() or any(USE[(k, st)] for st in STATES):
                z["agree"] = "NOT FORMED"
            else:
                z["dde"] = ((e["TS"] - CT[(k, "TS")]["E"]) - (e["R"] - CT[(k, "R")]["E"])) * H2K
                z["dev"] = z["dde"] - w
                z["frac"] = abs(z["dev"]) / abs(w)
                if abs(w) < FLOOR:
                    z["agree"] = "NOT ASSESSABLE"
                elif abs(z["dev"]) < FRAC * abs(w):
                    z["agree"] = "WITHIN"
                else:
                    z["agree"] = "OUTSIDE"
            if z["verdict"] == "UNRELIABLE" or z["fifth"] == "NOT PASSED: METASTABLE" or z["agree"] == "OUTSIDE":
                z["at"] = "FAILS"
            elif z["verdict"] == "RELIABLE" and z["fifth"] == "PASSED" and z["agree"] == "WITHIN":
                z["at"] = "PASSES"
            else:
                z["at"] = "OPEN"
            Rz[(k, d["name"])] = z
    for d in designs:
        ats = [Rz[(lv["key"], d["name"])]["at"] for lv in levels]
        d["cross"] = ("DOES NOT PASS THE CROSS-CHECK" if "FAILS" in ats else
                      "PASSES THE CROSS-CHECK" if ats == ["PASSES"] * len(levels) else "THE CROSS-CHECK IS INCOMPLETE")

    # ------------------------------------------------------------ the gates, as they held
    say("=" * 78)
    say("G44.1 TO G44.7  THE STOPPING GATES")
    say("=" * 78)
    say("")
    say("   G44.1  all %d files of INPUT_CHECKSUMS.sha256 are the ones listed. Of them, %d" % (nsum, nfile))
    say("          geometry, input and point-charge files: every geometry is d33a's; every input,")
    say("          built again here from the committed files, is the one that ran; every control")
    say("          without the stability analysis is the committed charge-free input of its level")
    say("          byte for byte; every point-charge file is D43's, byte for byte.")
    say("   G44.2  the five levels are those of A7.3, each method line read from the committed")
    say("          input AMENDMENT 7 names for it; the %d designs are the survivors d43.tsv lists;" % len(designs))
    say("          CONDITIONS.tsv, LEVELS.tsv and the five job scripts are what those give.")
    say("   G44.3  the %d files of SOURCES.sha256 are the ones listed; no committed file under" % len(srcs))
    say("          %s or %s differs from HEAD; the pre-registration and" % (", ".join(GUARDED[:-1]), GUARDED[-1]))
    say("          every one of those files are the files at the commit the pre-registration was")
    say("          committed on, %s." % pcm[:8])
    say("   G44.4  the working level, recomputed from D43's committed outputs, gives the ddE and")
    say("          the detachment energies d43.tsv prints for each design; the committed")
    say("          charge-free outputs of the five levels are readable.")
    say("   G44.5  the sentences quoted here are in the documents they are quoted from, and the")
    say("          tolerance of %d %% and its floor of %.1f kcal/mol are %s's." % (round(100 * FRAC), FLOOR, CRIT))
    say("   G44.6  every output that exists prints, as the input it ran, the committed input beside it.")
    say("   G44.7  every dianion run that was read reports the basis functions of its level:")
    say("          %s." % ", ".join("%s %s" % (lv["key"], OLD[(lv["key"], "R")]["c"]["nb"]) for lv in levels))
    say("")

    # ------------------------------------------------------------ R44.1
    say("=" * 78)
    say("R44.1  THE CONTROLS OF EACH LEVEL")
    say("=" * 78)
    say("")
    say("   C_ctl: the charge-free dianion, beside the committed charge-free energy of the same")
    say("   level, which D37 and D37 Amendment 1 ran without the stability analysis. At a double")
    say("   hybrid C_ctl is that committed input byte for byte and must REPRODUCE its energy within")
    say("   %.0e Eh for the level to be read. At CAMTZ it carries the analysis, so the comparison is" % TOL_EH)
    say("   a reading and its stability is the condition.")
    say("")
    say("   level  st   E control, Eh     committed, Eh      differ     stability          TRAH  word")
    tail = []
    for lv in levels:
        k = lv["key"]
        for st in STATES:
            c, o = CT[(k, st)], OLD[(k, st)]
            if c["status"] != "measured":
                say("   %-5s  %-2s   %s (G44.8)" % (k, st, c["status"]))
                for w in c.get("why") or []:
                    say("        %s" % w)
            else:
                dif = abs(c["E"] - o["E"])
                say("   %-5s  %-2s   %.8f   %.8f   %.2e   %-17s  %-4s  %s" % (
                    k, st, c["E"], o["E"], dif, "no analysis" if lv["dh"] else sbcell(c["stab"]), "yes" if c["trah"] else "no",
                    "REPRODUCES" if dif <= TOL_EH else "DIFFERS"))
                if not lv["dh"] and c["stab"].get("first") == "UNSTABLE":
                    tail.append("%s %s: the control's first solution was unstable%s; its last analysis is %s. <S**2> of the solution it ended on: %s." % (
                        k, st, negs(c["stab"]), c["stab"].get("last"), fs(c.get("s2"))))
            if USE[(k, st)]:
                tail.append("%s %s: THE CONTROL CANNOT BE USED: %s. The legs that need it and the ddE of this level are not formed." % (k, st, USE[(k, st)]))
    say("")
    notes(tail)
    if tail:
        say("")
    say("   C_ctlS: the stability search on the charge-free dianion, at the double hybrids.")
    say("   differ: its final energy minus C_ctl's. rise, dip: of the solution it ended on, against")
    say("   C_ctl's. The words are those of R44.2.")
    say("")
    say("   level  st   stability          <S**2>     differ, Eh    rise     dip    the search")
    tail = []
    for lv in levels:
        k = lv["key"]
        if not lv["dh"]:
            continue
        for st in STATES:
            s, c = CS[(k, st)], CT[(k, st)]
            if s["status"] != "measured":
                say("   %-5s  %-2s   %s (G44.8)" % (k, st, s["status"]))
                for w in s.get("why") or []:
                    say("        %s" % w)
                continue
            sb = s["stab"]
            ok = c["status"] == "measured"
            if not sb["ran"] or sb["first"] is None or sb["last"] is None:
                word = "NOT READ"
            elif sb["last"] == "UNSTABLE":
                word = "UNSTABLE AT THE END"
            elif sb["first"] == "STABLE":
                word = "NOTHING LOWER"
            elif ok:
                f2 = []
                three(s, c, f2, [])
                word = "CHARGE MOVED" if f2 else "LOWER, CHARGE IN PLACE"
            else:
                word = "LOWER, NOT READ"
            say("   %-5s  %-2s   %-17s  %-8s  %s  %s  %s   %s" % (
                k, st, sbcell(sb), fs(s.get("s2")), fs(s["E"] - c["E"] if ok else None, "%+.2e", "        -"),
                fs(s["c"]["maxlw"] - c["c"]["maxlw"] if ok else None, "%+.4f", "      -"),
                fs(s["c"]["dip"] / c["c"]["dip"] if ok and c["c"]["dip"] else None, "%5.2f", "    -"), word))
            if sb.get("first") == "UNSTABLE":
                tail.append("%s %s: the restricted solution of the charge-free dianion is unstable%s at this level (F5). The control's energy stays the restricted one, which is the committed one." % (k, st, negs(sb)))
    say("")
    notes(tail)
    if tail:
        say("")
    say("   C_bare: the bare monoanion, which gives the bare substrate's detachment energy.")
    say("")
    say("   level   barrier, kcal/mol   committed    differ      bare VDE at R, eV   at TS, eV   monoanion stability R, TS")
    say("   %-6s  %9.4f           -            -           %s            %s     (the working level, from d31 and d40)" % (
        "work", barw, fs(barew["R"], "%+.4f"), fs(barew["TS"], "%+.4f")))
    BARE = {}
    for lv in levels:
        k = lv["key"]
        bo = (OLD[(k, "TS")]["E"] - OLD[(k, "R")]["E"]) * H2K
        bn = None if any(CT[(k, st)]["E"] is None for st in STATES) else (CT[(k, "TS")]["E"] - CT[(k, "R")]["E"]) * H2K
        for st in STATES:
            BARE[(k, st)] = None if (BM[(k, st)]["E"] is None or CT[(k, st)]["E"] is None) else (BM[(k, st)]["E"] - CT[(k, st)]["E"]) * EV
        say("   %-6s  %s           %9.4f    %s   %s            %s     %s, %s" % (
            k, fs(bn, "%9.4f", "NOT FORMED"), bo, fs(None if bn is None else abs(bn - bo), "%.2e", "       -"),
            fs(BARE[(k, "R")], "%+.4f", "      -"), fs(BARE[(k, "TS")], "%+.4f", "      -"),
            sbcell(BM[(k, "R")].get("stab") or {}) if BM[(k, "R")]["E"] is not None else "NOT MEASURED",
            sbcell(BM[(k, "TS")].get("stab") or {}) if BM[(k, "TS")]["E"] is not None else "NOT MEASURED"))
    say("")
    for lv in levels:
        for st in STATES:
            m = BM[(lv["key"], st)]
            if m["E"] is None:
                wrapsay("NOTE  %s: NOT MEASURED. %s" % (prep.cname(lv["key"], "mo", st), "; ".join(m.get("why") or [])))
    say("   The bare detachment energies are reported, not interpreted; <S**2> of each bare")
    say("   monoanion is in d44.tsv.")
    say("")

    # ------------------------------------------------------------ R44.2
    say("=" * 78)
    say("R44.2  THE FOUR LEGS AT EACH LEVEL, ON THE DIANION")
    say("=" * 78)
    say("")
    say("   rise: the largest Loewdin charge minus the level's control's, e; flagged at %.2f." % D_LW)
    say("   dip: the dipole over the control's; flagged at %.1f." % R_DIP)
    say("   gap: eV; flagged below %.1f and below %.1f of the control's; 'withheld' where TRAH or" % (GAP_ABS, R_GAP))
    say("        a level shift was engaged in the run or in its control.")
    say("   read: the number of point charges ORCA says it read.")
    say("   At CAMTZ, stability is the dianion's first analysis, then its last.")
    say("   At a double hybrid, rise, dip and gap are of di, the restricted solution whose energy")
    say("   is used, and the search is ds: its first and last analysis, <S**2> of the solution it")
    say("   ended on, that solution's rise and dip against the control, and its final energy minus")
    say("   di's. NOTHING LOWER: the first analysis is STABLE. LOWER, CHARGE IN PLACE: it found a")
    say("   lower solution that passes the three legs; not a flag. CHARGE MOVED: it found one that")
    say("   does not; a flag. UNSTABLE AT THE END: a flag.")
    say("")
    for lv in levels:
        k = lv["key"]
        say("   %s  %s%s" % (k, lv["name"], "" if all(USE[(k, st)] is None for st in STATES) else "   THE CONTROL CANNOT BE USED (R44.1)"))
        gc = [CT[(k, st)]["c"].get("gap") if CT[(k, st)]["status"] == "measured" else None for st in STATES]
        say("   the control's gap: %s eV at R, %s eV at TS" % (fs(gc[0], "%.4f"), fs(gc[1], "%.4f")))
        if lv["dh"]:
            say("   design        st    rise    dip    gap        TRAH shift read | the search: stability   <S**2>    rise     dip   differ, Eh  word | flags")
        else:
            say("   design        st    rise    dip    gap        stability          TRAH shift read  flags")
        tail = []
        for d in designs:
            for st in STATES:
                x = X[(k, d["name"], st)]
                di = x["di"]
                if di["status"] != "measured":
                    say("   %-13s %-2s  %s" % (d["name"], st, "DID NOT CONVERGE: a flag (D31.3)" if di["status"] == "DID NOT CONVERGE" else "NOT MEASURED (G44.8)"))
                    for w in di.get("why") or []:
                        say("        %s" % w)
                else:
                    mm = re.search(r"Reading point charge file\s*\.*\s*ok \((\d+) point charges?\)", di["text"])
                    x["npc"] = int(mm.group(1)) if mm else None
                    head = "   %-13s %-2s  %s  %s  %-9s" % (
                        d["name"], st, fs(x["rise"], "%+.4f", "      -"), fs(x["ratio"], "%5.2f", "    -"),
                        ("%.4f" % x["gap"]) if x["gap_read"] else ("-" if x["noctl"] else "withheld"))
                    tsr = "%-4s %-5s %-4s" % ("yes" if di["trah"] else "no", "yes" if di["lshift"] else "no", "-" if x["npc"] is None else str(x["npc"]))
                    if lv["dh"]:
                        ds = x["ds"]
                        if ds["status"] == "measured":
                            cell = "%-17s  %-8s  %s  %s  %s  %s" % (
                                sbcell(ds["stab"]), fs(ds.get("s2")), fs(x["srise"], "%+.4f", "      -"), fs(x["sratio"], "%5.2f", "    -"),
                                fs(x["sdiff"], "%+.2e", "        -"), x["search"])
                        else:
                            cell = "%s" % ("DID NOT CONVERGE" if ds["status"] == "DID NOT CONVERGE" else "NOT MEASURED (G44.8)")
                        say("%s  %s | %s | %s" % (head, tsr, cell, "; ".join(x["flags"]) or "none"))
                    else:
                        say("%s  %-17s  %s  %s" % (head, sbcell(di["stab"]), tsr, "; ".join(x["flags"]) or "none"))
                for n_ in x["notes"]:
                    tail.append("%s %s: %s." % (d["name"], st, n_))
                for h_ in x["held"]:
                    if h_ != "gap withheld":
                        tail.append("%s %s: %s." % (d["name"], st, h_))
                r_ = x["ds"] if lv["dh"] else di
                if r_ is not None and r_["status"] == "measured":
                    sb = r_["stab"]
                    if not sb["ran"]:
                        tail.append("%s %s: no stability analysis was found in %s, though its input asks for one. Lines of the output that mention stability: %s" % (
                            d["name"], st, x["rds"] if lv["dh"] else x["rdi"], " | ".join(sb["quote"]) or "none"))
                    elif sb["first"] is None:
                        tail.append("%s %s: %s. ORCA's own lines: %s" % (d["name"], st, sb["how"], " | ".join(sb["quote"]) or "(none found)"))
                    elif sb["how"] != "cs6_d33b.stability":
                        tail.append("%s %s: stability read from %s." % (d["name"], st, sb["how"]))
                elif r_ is not None and lv["dh"]:
                    for w in r_.get("why") or []:
                        tail.append("%s %s: the stability search %s." % (d["name"], st, w))
        notes(tail)
        say("")
    say("   design         %s" % "  ".join("%-14s" % lv["key"] for lv in levels))
    for d in designs:
        say("   %-13s  %s" % (d["name"], "  ".join("%-14s" % Rz[(lv["key"], d["name"])]["verdict"] for lv in levels)))
    say("")
    say("   RELIABLE: no flag at either state, with the gap readable and the stability analysis")
    say("   read at both. UNRELIABLE: a flag at either state. INCONCLUSIVE: no flag, and a gap")
    say("   withheld, or a stability analysis or search not read, or a search that did not start")
    say("   from the solution whose energy is used. NOT CLASSIFIED: no flag, and a state NOT")
    say("   MEASURED or a control of the level that cannot be used.")
    say("")

    # ------------------------------------------------------------ R44.3
    say("=" * 78)
    say("R44.3  THE FIFTH LEG AT EACH LEVEL: THE DETACHMENT ENERGY WITH THE CHARGES PRESENT")
    say("=" * 78)
    say("")
    say("   BOUND: VDE above zero and the monoanion's last stability analysis STABLE.")
    say("   METASTABLE: VDE at or below zero. change: VDE minus the bare substrate's at the same")
    say("   level and state (R44.1). Stability is the monoanion's, first analysis then last. The")
    say("   dianion energy is di's.")
    say("")
    anyns = False
    for lv in levels:
        k = lv["key"]
        say("   %s  %s" % (k, lv["name"]))
        say("   design        st    VDE eV    change   <S**2>    stability          TRAH shift read  word")
        tail = []
        for d in designs:
            for st in STATES:
                x = X[(k, d["name"], st)]
                m, di = x["mo"], x["di"]
                if x["v"] is None:
                    say("   %-13s %-2s  NOT MEASURED (G44.8, G44.9)" % (d["name"], st))
                    if m["E"] is None:
                        for w in m.get("why") or []:
                            say("        the monoanion %s" % w)
                    if di["E"] is None:
                        say("        the dianion of this state is not read (R44.2)")
                    continue
                sb = m["stab"]
                mm = re.search(r"Reading point charge file\s*\.*\s*ok \((\d+) point charges?\)", m["text"])
                x["npc_mo"] = int(mm.group(1)) if mm else None
                b0 = BARE[(k, st)]
                say("   %-13s %-2s  %+8.4f  %s  %s  %-17s  %-4s %-5s %-4s  %s" % (
                    d["name"], st, x["v"], fs(None if b0 is None else x["v"] - b0, "%+8.4f", "       -"), fs(m["s2"]), sbcell(sb),
                    "yes" if m["trah"] else "no", "yes" if m["lshift"] else "no", "-" if x["npc_mo"] is None else str(x["npc_mo"]), x["word"]))
                if not sb["ran"]:
                    tail.append("%s %s: no stability analysis was found in the monoanion's output, though its input asks for one. Lines of the output that mention stability: %s" % (
                        d["name"], st, " | ".join(sb["quote"]) or "none"))
                elif sb["first"] is None:
                    tail.append("%s %s: %s. ORCA's own lines: %s" % (d["name"], st, sb["how"], " | ".join(sb["quote"]) or "(none found)"))
                else:
                    if sb["how"] != "cs6_d33b.stability":
                        tail.append("%s %s: the monoanion's stability read from %s." % (d["name"], st, sb["how"]))
                    if sb["first"] == "UNSTABLE":
                        tail.append("%s %s: the monoanion's first solution was unstable%s; %s." % (
                            d["name"], st, negs(sb), "1 analysis was made" if sb["n"] == 1 else "%d analyses were made" % sb["n"]))
                if x["word"] == NOTSTABLE:
                    anyns = True
                    tail.append("%s %s: the value is above zero and the monoanion's last stability analysis is not STABLE, so by AMENDMENT 6 it is not counted as BOUND." % (d["name"], st))
        notes(tail)
        say("")
    say("   design         level   R               TS              the fifth leg")
    for d in designs:
        for lv in levels:
            k = lv["key"]
            say("   %-13s  %-5s   %-14s  %-14s  %s" % (d["name"], k, cellw(k, d, "R"), cellw(k, d, "TS"), Rz[(k, d["name"])]["fifth"]))
    if anyns:
        say("   * %s (AMENDMENT 6)." % NOTSTABLE)
    say("")
    say("   The change is reported, not interpreted. D40.4: where a detachment energy is negative,")
    say("   its sign is the reading and its size is not. <S**2> is reported and not corrected; the")
    say("   ideal value is 0.75. Both energies of every row are in d44.tsv.")
    say("")

    # ------------------------------------------------------------ R44.4
    say("=" * 78)
    say("R44.4  THE BARRIER CHANGE AT EACH LEVEL, BESIDE THE WORKING LEVEL")
    say("=" * 78)
    say("")
    say("   ddE = (E(TS) - E(TS, control)) - (E(R) - E(R, control)) on the dianion, kcal/mol, with")
    say("   the controls of the same level. working: the design's ddE at %s," % WORKING)
    say("   from D43's committed outputs. differ: ddE minus working. frac: |differ| over |working|.")
    say("")
    say("   The tolerance, from %s:" % CRIT)
    for l in textwrap.wrap('"%s" "%s"' % (CRIT_A, CRIT_B), 72):
        say("      %s" % l)
    say("   WITHIN: frac below %.2f. OUTSIDE: frac at or above %.2f, in either direction." % (FRAC, FRAC))
    say("   NOT ASSESSABLE: |working| below %.1f kcal/mol, the floor of the same document." % FLOOR)
    say("   NOT FORMED: a dianion state of the design is not read, or the level's control cannot")
    say("   be used.")
    say("")
    say("   design         level      ddE     working    differ    frac   agreement        four legs       the fifth leg")
    for d in designs:
        w = W[d["name"]]["dde"]
        for lv in levels:
            k = lv["key"]
            z = Rz[(k, d["name"])]
            say("   %-13s  %-5s  %s  %+9.4f  %s  %s   %-15s  %-14s  %s" % (
                d["name"], k, fs(z["dde"], "%+9.4f", "        -"), w, fs(z["dev"], "%+8.4f", "       -"), fs(z["frac"], "%.3f", "    -"),
                z["agree"], z["verdict"], z["fifth"]))
    say("")
    say("   A ddE is printed for a level at which the design is UNRELIABLE or not BOUND, as D33's")
    say("   and D43's results did; the words beside it govern. Where a state is METASTABLE at a")
    say("   level, Probe 2's sentence applies to the ddE of that level:")
    for l in textwrap.wrap('"%s"' % P2_SENT, 72):
        say("      %s" % l)
    say("")

    # ------------------------------------------------------------ R44.5
    say("=" * 78)
    say("R44.5  THE CROSS-CHECK OF EACH DESIGN")
    say("=" * 78)
    say("")
    say("   At a level a design PASSES if it is RELIABLE on four legs, passes the fifth and its ddE")
    say("   is WITHIN. It FAILS if it is UNRELIABLE, or METASTABLE at a state, or its ddE is")
    say("   OUTSIDE. Anything else is OPEN: not established either way.")
    say("   A design PASSES THE CROSS-CHECK if it passes at all five levels. It DOES NOT PASS if")
    say("   it fails at any. Otherwise THE CROSS-CHECK IS INCOMPLETE for it.")
    say("")
    say("   design         %s   the cross-check" % "  ".join("%-6s" % lv["key"] for lv in levels))
    for d in designs:
        say("   %-13s  %s   %s" % (d["name"], "  ".join("%-6s" % Rz[(lv["key"], d["name"])]["at"] for lv in levels), d["cross"]))
    say("")
    for d in designs:
        why = []
        for lv in levels:
            z = Rz[(lv["key"], d["name"])]
            if z["at"] == "PASSES":
                continue
            r = []
            if z["verdict"] != "RELIABLE":
                r.append("four legs %s" % z["verdict"])
            if z["fifth"] != "PASSED":
                r.append("the fifth leg %s" % z["fifth"])
            if z["agree"] != "WITHIN":
                r.append("ddE %s%s" % (z["agree"], "" if z["dde"] is None else " (%+.4f against %+.4f)" % (z["dde"], W[d["name"]]["dde"])))
            why.append("%s %s: %s" % (lv["key"], z["at"], ", ".join(r)))
        if why:
            wrapsay("%s. %s." % (d["name"], "; ".join(why)))
    npass = sum(1 for d in designs if d["cross"] == "PASSES THE CROSS-CHECK")
    nfail = sum(1 for d in designs if d["cross"] == "DOES NOT PASS THE CROSS-CHECK")
    say("")
    say("   Of the %d designs. Pass the cross-check: %d. Do not pass: %d. Incomplete: %d." % (
        len(designs), npass, nfail, len(designs) - npass - nfail))
    say("")
    say("   AMENDMENT 7, A7.4:")
    for l in textwrap.wrap('"%s"' % A7_C, 72):
        say("      %s" % l)
    say("   D44 names nothing.")
    say("")
    for lv in levels:
        k = lv["key"]
        every = [X[(k, d["name"], st)][t] for d in designs for st in STATES for t in (("di", "ds", "mo") if lv["dh"] else ("di", "mo"))]
        withs = [X[(k, d["name"], st)][t] for d in designs for st in STATES for t in ("ds", "mo")] if lv["dh"] else []
        if not any(r["status"] == "measured" for r in every):
            say("   F4 APPLIES AT %s: none of its %d charged runs is measured. The level did not run, and" % (k, len(every)))
            say("   that is the primary result for it. A7.3: \"None is dropped.\"")
            say("")
        elif withs and not any(r["status"] == "measured" for r in withs):
            say("   F4 APPLIES AT %s: none of its %d charged runs with the stability analysis (ds, mo) is" % (k, len(withs)))
            say("   measured, so the stability leg and the fifth leg are not established there. ORCA's")
            say("   own words for each are above. A7.3: \"None is dropped.\"")
            say("")

    # ------------------------------------------------------------ R44.6
    say("=" * 78)
    say("R44.6  THE SIX LEVELS SIDE BY SIDE")
    say("=" * 78)
    say("")
    say("   ddE, kcal/mol. work is the working level. range: the smallest and the largest of the")
    say("   six, and their difference; it is formed only where all five levels give a ddE.")
    say("")
    say("   design           work     %s    range" % "    ".join("%-6s" % lv["key"] for lv in levels))
    for d in designs:
        cells, vals = [], [W[d["name"]]["dde"]]
        for lv in levels:
            z = Rz[(lv["key"], d["name"])]
            if z["dde"] is None:
                cells.append("%-9s" % "    -")
            else:
                cells.append("%+9.4f" % z["dde"])
                vals.append(z["dde"])
        rng = "%+.4f to %+.4f, %.4f" % (min(vals), max(vals), max(vals) - min(vals)) if len(vals) == len(levels) + 1 else "not formed"
        say("   %-13s %+9.4f  %s   %s" % (d["name"], W[d["name"]]["dde"], " ".join(cells), rng))
    say("")
    say("   The smaller detachment energy of the two states, eV.")
    say("")
    say("   design           work     %s" % "    ".join("%-6s" % lv["key"] for lv in levels))
    for d in designs:
        cells = []
        for lv in levels:
            vs = [X[(lv["key"], d["name"], st)]["v"] for st in STATES]
            cells.append("%-9s" % "    -" if None in vs else "%+9.4f" % min(vs))
        say("   %-13s %+9.4f  %s" % (d["name"], min(W[d["name"]]["v"].values()), " ".join(cells)))
    say("")
    say("   The larger basis minus the smaller, same functional, and the double hybrid minus")
    say("   CAM-B3LYP, same basis: differences of ddE, kcal/mol, formed where both exist.")
    say("")
    pairs = (("CAMTZ - work", "CAMTZ", None), ("PWTZ - PWSV", "PWTZ", "PWSV"), ("RDTZ - RDSV", "RDTZ", "RDSV"),
             ("PWSV - work", "PWSV", None), ("RDSV - work", "RDSV", None), ("PWTZ - CAMTZ", "PWTZ", "CAMTZ"), ("RDTZ - CAMTZ", "RDTZ", "CAMTZ"))
    say("   design         %s" % "  ".join("%-12s" % p[0] for p in pairs))
    for d in designs:
        cells = []
        for label, k1, k2 in pairs:
            v1 = Rz[(k1, d["name"])]["dde"]
            v2 = W[d["name"]]["dde"] if k2 is None else Rz[(k2, d["name"])]["dde"]
            cells.append("%-12s" % ("-" if None in (v1, v2) else "%+.4f" % (v1 - v2)))
        say("   %-13s  %s" % (d["name"], "  ".join(cells)))
    say("")
    say("   Reported, not interpreted. Nothing is fitted, extrapolated or averaged. A7.4:")
    say("      \"%s\"" % A7_D)
    say("")

    # ------------------------------------------------------------ limits
    say("=" * 78)
    say("WHAT THIS CANNOT SAY (section 7 of the pre-registration, in short; its own words govern)")
    say("=" * 78)
    say("")
    for l in ("D44.1  It adopts no production level of theory and does not say which benchmark is nearer",
              "       the truth. It names no design: what a survivor that passes may be called is a dated",
              "       document and a decision for Jessica (A7.4).",
              "D44.2  One frame, frozen geometries, in vacuo, no continuum. No barrier here is a prediction",
              "       of a physical barrier.",
              "D44.3  The charges are bare point charges and stand for no residue. No wall is in any run.",
              "D44.4  Three functionals in two basis sets are not independent tests. The working level and",
              "       CAMTZ share a functional. Whether a functional's error on the barrier cancels in ddE",
              "       is not measured by their agreement.",
              "D44.5  The tolerance is the record's one constant for method and basis, fixed on 21 October",
              "       2026 for another use. It was applied here with D43's six numbers known.",
              "D44.6  At a double hybrid the energy is on the restricted solution whether or not the",
              "       search finds that solution stable, the three legs read a self-consistent density and",
              "       not one relaxed with the perturbative term, and the monoanion is an unrestricted",
              "       reference whose spin contamination is reported and not corrected.",
              "D44.7  The six designs share two sites and were placed by a greedy rule on one grid.",
              "D44.8  The detachment energies are vertical and in vacuo; where one is negative its sign is",
              "       the reading, not its size.",
              "D44.9  It does not re-measure r_safe, does not take route (a) and does not adopt the wall.",
              "D44.10 Nothing here is shaped: the tolerance is symmetric and every level is reported",
              "       whichever way it comes out."):
        say("   %s" % l)
    say("")
    allq = [(lv, d["name"], st) for lv in levels for d in designs for st in STATES]
    ndi = sum(1 for lv, n, st in allq if X[(lv["key"], n, st)]["di"]["status"] == "measured")
    nnc = sum(1 for lv, n, st in allq if X[(lv["key"], n, st)]["di"]["status"] == "DID NOT CONVERGE")
    nmo = sum(1 for lv, n, st in allq if X[(lv["key"], n, st)]["mo"]["E"] is not None)
    nds = sum(1 for lv, n, st in allq if lv["dh"] and X[(lv["key"], n, st)]["ds"]["status"] == "measured")
    tds = sum(1 for lv, n, st in allq if lv["dh"])
    nct = sum(1 for lv in levels for st in STATES if CT[(lv["key"], st)]["status"] == "measured") \
        + sum(1 for lv in levels for st in STATES if BM[(lv["key"], st)]["E"] is not None) \
        + sum(1 for lv in levels for st in STATES if lv["dh"] and CS[(lv["key"], st)]["status"] == "measured")
    tct = sum(6 if lv["dh"] else 4 for lv in levels)
    say("%d of the %d dianion runs (di), %d of the %d stability searches (ds), %d of the %d monoanion runs and %d of the %d controls are measured.%s%s" % (
        ndi, len(allq), nds, tds, nmo, len(allq), nct, tct,
        " %d dianion run%s did not converge, which is a flag." % (nnc, "" if nnc == 1 else "s") if nnc else "",
        "" if (ndi + nnc == len(allq) and nmo == len(allq) and nds == tds and nct == tct) else " The rest are NOT MEASURED, each listed above with ORCA's own words."))

    # ------------------------------------------------------------ the table
    L = ["design\tlevel_key\tlevel\tdouble_hybrid\tstate\tn\tsum_q_e\t"
         "di_run\tdi_status\tE_di_Eh\tE_control_Eh\tlowdin_rise_e\tdipole_ratio\tgap_eV\tgap_read\tdi_stab_first\tdi_stab_last\tdi_S2\tdi_TRAH\tdi_level_shift\t"
         "di_charges_read\t"
         "ds_run\tds_status\tE_ds_Eh\tE_ds_minus_E_di_Eh\tds_stab_first\tds_stab_last\tds_S2\tds_lowdin_rise_e\tds_dipole_ratio\tds_gap_eV\tsearch\tflags\theld\t"
         "mo_run\tmo_status\tE_mo_Eh\tS2\tmo_stab_first\tmo_stab_last\tmo_TRAH\tmo_level_shift\tmo_charges_read\tVDE_eV\tword\tbare_VDE_eV\tdVDE_vs_bare_eV\t"
         "four_legs\tfifth_leg\tddE_kcal\tddE_working_kcal\tddE_minus_working_kcal\tfraction_of_working\tagreement\tat_this_level\tcross_check"]

    def sc(r, key):
        sb = r.get("stab") or {}
        return sb.get(key) or ("NOT FOUND" if not sb.get("ran") else "NOT READ")

    for d in designs:
        for lv in levels:
            k = lv["key"]
            z = Rz[(k, d["name"])]
            for st in STATES:
                x = X[(k, d["name"], st)]
                di, m, c0, ds = x["di"], x["mo"], CT[(k, st)], x["ds"]
                f = [d["name"], k, lv["name"], "yes" if lv["dh"] else "no", st, "%d" % d["n"], d["sum"], x["rdi"], di["status"]]
                if di["status"] == "measured":
                    f += ["%.8f" % di["E"], "" if c0["E"] is None else "%.8f" % c0["E"], "" if x["rise"] is None else "%+.4f" % x["rise"],
                          "" if x["ratio"] is None else "%.4f" % x["ratio"], "" if x["gap"] is None else "%.4f" % x["gap"],
                          "yes" if x["gap_read"] else "no", "" if lv["dh"] else sc(di, "first"), "" if lv["dh"] else sc(di, "last"),
                          "" if di.get("s2") is None else "%.6f" % di["s2"],
                          "yes" if di["trah"] else "no", "yes" if di["lshift"] else "no", "" if x.get("npc") is None else str(x["npc"])]
                else:
                    f += ["", "" if c0["E"] is None else "%.8f" % c0["E"]] + [""] * 10
                if not lv["dh"]:
                    f += [""] * 11
                elif ds["status"] != "measured":
                    f += [x["rds"], ds["status"]] + [""] * 8 + [x["search"]]
                else:
                    f += [x["rds"], ds["status"], "%.8f" % ds["E"], "" if x["sdiff"] is None else "%+.3e" % x["sdiff"], sc(ds, "first"), sc(ds, "last"),
                          "" if ds.get("s2") is None else "%.6f" % ds["s2"], "" if x["srise"] is None else "%+.4f" % x["srise"],
                          "" if x["sratio"] is None else "%.4f" % x["sratio"], "" if x["sgap"] is None else "%.4f" % x["sgap"], x["search"]]
                f += ["; ".join(x["flags"]), "; ".join(x["held"]), x["rmo"], m["status"]]
                if m["E"] is not None:
                    f += ["%.8f" % m["E"], "" if m["s2"] is None else "%.6f" % m["s2"], sc(m, "first"), sc(m, "last"),
                          "yes" if m["trah"] else "no", "yes" if m["lshift"] else "no", "" if x.get("npc_mo") is None else str(x["npc_mo"])]
                else:
                    f += [""] * 7
                b0 = BARE[(k, st)]
                if x["v"] is not None:
                    f += ["%+.4f" % x["v"], x["word"], "" if b0 is None else "%+.4f" % b0, "" if b0 is None else "%+.4f" % (x["v"] - b0)]
                else:
                    f += ["", "", "" if b0 is None else "%+.4f" % b0, ""]
                f += [z["verdict"], z["fifth"], "" if z["dde"] is None else "%+.4f" % z["dde"], "%+.4f" % W[d["name"]]["dde"],
                      "" if z["dev"] is None else "%+.4f" % z["dev"], "" if z["frac"] is None else "%.4f" % z["frac"], z["agree"], z["at"], d["cross"]]
                L.append("\t".join(f))
    C = ["control\tlevel_key\tstate\tkind\tstatus\tE_Eh\tcommitted_charge_free_E_Eh\tstab_first\tstab_last\tS2\tTRAH\tlevel_shift\tcan_be_used"]
    for lv in levels:
        k = lv["key"]
        for st in STATES:
            rows = [("di", CT[(k, st)])] + ([("ds", CS[(k, st)])] if lv["dh"] else []) + [("mo", BM[(k, st)])]
            for kind, r in rows:
                ok = r["status"] == "measured"
                stabbed = ok and not (lv["dh"] and kind == "di")
                C.append("\t".join([prep.cname(k, kind, st), k, st, kind, r["status"], "%.8f" % r["E"] if ok else "",
                                    "%.8f" % OLD[(k, st)]["E"] if kind != "mo" else "", sc(r, "first") if stabbed else "", sc(r, "last") if stabbed else "",
                                    "" if r.get("s2") is None else "%.6f" % r["s2"], ("yes" if r["trah"] else "no") if ok else "",
                                    ("yes" if r["lshift"] else "no") if ok else "",
                                    ("no: %s" % USE[(k, st)] if USE[(k, st)] else "yes") if kind == "di" else ""]))
    rep[:] = [l.rstrip() for l in rep]
    open(a.tsv, "w").write("\n".join(L) + "\n")
    open(os.path.splitext(a.tsv)[0] + "_controls.tsv", "w").write("\n".join(C) + "\n")
    open(a.result, "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main())
