#!/usr/bin/env python3
"""DIAGNOSTIC 43's reader: bare charges at or beyond r_safe, on the committed grid.

Every condition is one point-charge file run four times: the dianion and the monoanion, at the
reactant and at the transition state. Three things are read from them.

THE FOUR LEGS, on the dianion, against d31's committed charge-free control at the same state.
The thresholds are those of the committed classifier cs6_audit.py, read from that file and not
typed in, and its read() gives the gap, the largest Loewdin charge, the dipole and the number
of basis functions. The rules are the ones D33.3 and AMENDMENT 5 settled and cs6_d33b.py holds:
stability on the LAST analysis; the gap withheld in a run that engaged TRAH or a level shift.

THE FIFTH LEG, by AMENDMENT 6. VDE = E(monoanion) - E(dianion) at the same state with the same
point charges, positive means bound, at 27.211386245988 eV per Eh. The words are D42's, by its
committed word_of(): BOUND needs a value above zero AND a monoanion whose last stability
analysis is STABLE; METASTABLE is a value at or below zero.

ddE, on the dianion, against the same control, at 627.5094740631 kcal/mol per Eh.

Reading one output. The committed cs6_d40.py is imported by path and its read() is used as it
stands, for dianions and monoanions alike: a run is read only if it terminates normally,
converges its SCF, reports the right number of electrons and says it read as many point
charges as its pc.pc holds; the energy is the last FINAL SINGLE POINT ENERGY; the stability
analysis is read by cs6_d33b.py stability() with ORCA's own verdict sentences as the fallback,
and is never guessed. A dianion output that says the SCF did not converge is a flag (D31.3).

The selection is not taken on trust. The committed cs6_d43_prepare.py is imported by path and
its select() is run again here on the committed files; what it gives must be what
d43/CONDITIONS.tsv lists and what the point-charge files hold. That part needs the standard
library only.

Nothing is typed in except committed reference values, each named with its source and each
checked against that source before it is used.

usage: cs6_d43.py --res <results-repo> --code <code-repo> --prereg <file> --result <file> --tsv <file>
"""
import argparse, importlib.util, math, os, re, sys, textwrap

EV = 27.211386245988            # eV per Eh, as D40 and D42 use it
H2K = 627.5094740631
S = "phase2.2/cs_stage6"
NEL_MONO, NEL_DI = 117, 118
CTL_DI_E = -836.31486034        # G43.4: d31/ctl_R, which is D28's CAM-B3LYP-D3BJ/def2-SVPD reactant (cs6_d40.py CTL_E)
CTL_MO_E = -836.26897688        # G43.4: D40's committed monoanion B_bare_R, d40.tsv
TOL_EH = 1e-6
BARE = {"R": "+1.2486", "TS": "+0.7681"}   # D40's bare detachment energies, d40.tsv rows B_bare_R and B_bare_TS
D31_BARRIER, TOL_K = 26.8248, 5e-5      # G43.5: the charge-free barrier of d31/ctl_{R,TS}
CTL_SHA = {"R": "2f5148d27e932aad2280f0ea7960bdf32fafdf5d769d9efe8911fef0ccadf3f2",     # as D33.2 pins them
           "TS": "5a650d9ef9d41ff957af9018aa65b9d597dc584148aab78144b976b4a6ba1f80"}
P2_FILE = "vde_probe/PROBE2_RESULT.txt"
P2_SENT = "Any ddE measured at such a condition is not a barrier effect and must not be reported as one."
PREREG = "STAGE6_DIAGNOSTIC43_20261117.txt"
AMEND6 = "STAGE6_AMENDMENT6_20261117.txt"
D11P = "STAGE6_DIAGNOSTIC11_20261027.txt"
D11_Q = "one charge at one site is not a design"
H3_A = "Outside r_safe a bare charge is reliable at every magnitude up to one."
H3_B = "a unit charge at 3.50 A or beyond fires any leg of the criterion at CAM-B3LYP-D3BJ/def2-SVPD."
H2_A = "Relocation is onto the positive charges: an arrangement with no positive charge does not relocate, at any magnitude up to one."
H2_B = "any arrangement with no positive charge fires a leg of the criterion at CAM-B3LYP-D3BJ/def2-SVPD with the stability analysis"
H9_A = "a RELIABLE arrangement whose charges sum to zero or more keeps it bound"
H9_B = ("a RELIABLE arrangement whose charges sum to zero or more that is METASTABLE at either state, "
        "at another site, outside r_safe, or among designs an optimiser places")
GUARDED = ("d31", "d33a", "d40")
STATES = ("R", "TS")
NCOND, NRUN, NSUM = 18, 74, 224


def stop(m):
    sys.stderr.write("STOP  %s\n" % m)
    raise SystemExit(1)


def load(path, name):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    a.res, a.code, a.prereg = (os.path.abspath(x) for x in (a.res, a.code, a.prereg))
    S6 = os.path.join(a.res, S)
    D = os.path.join(S6, "d43")
    paths = dict((k, os.path.join(a.code, "phase2.2", f)) for k, f in (
        ("d40", "cs6_d40.py"), ("d33b", "cs6_d33b.py"), ("d42", "cs6_d42.py"), ("audit", "cs6_audit.py"), ("prep", "cs6_d43_prepare.py")))
    for p in paths.values():
        if not os.path.isfile(p):
            stop("%s is missing; this reader imports it by path." % p)
    d40 = load(paths["d40"], "d43_cs6_d40")
    d33b = load(paths["d33b"], "d43_cs6_d33b")
    d42 = load(paths["d42"], "d43_cs6_d42")
    audit = load(paths["audit"], "d43_cs6_audit")
    prep = load(paths["prep"], "d43_cs6_d43_prepare")
    sha, git = d40.sha, d40.git
    word_of, NOTSTABLE = d42.word_of, d42.NOTSTABLE
    D_LW, R_DIP, GAP_ABS, R_GAP = audit.D_LW, audit.R_DIP, audit.GAP_ABS, audit.R_GAP
    pam = os.path.join(S6, AMEND6)

    rep, gates = [], []
    say = rep.append
    say("STAGE 6, DIAGNOSTIC 43 - RESULT")
    say("Bare charges at or beyond r_safe, on the committed grid: do they pass the four legs,")
    say("does the dianion stay bound, and what do designs placed there do to the barrier?")
    say("")
    say("Read by phase2.2/cs6_d43.py                 sha256 %s" % sha(__file__))
    say("   every output by cs6_d40.py read           sha256 %s" % sha(paths["d40"]))
    say("   stability by cs6_d33b.py stability        sha256 %s" % sha(paths["d33b"]))
    say("   thresholds and legs by cs6_audit.py       sha256 %s" % sha(paths["audit"]))
    say("   the words of binding by cs6_d42.py        sha256 %s" % sha(paths["d42"]))
    say("   the selection by cs6_d43_prepare.py       sha256 %s" % sha(paths["prep"]))
    say("   pre-registration %-23s sha256 %s" % (os.path.basename(a.prereg), sha(a.prereg)))
    say("   %-40s sha256 %s" % (AMEND6, sha(pam) if os.path.isfile(pam) else "MISSING"))
    say("")
    say("Frame 42436, CAM-B3LYP-D3BJ/def2-SVPD with the stability analysis, bare point charges.")
    say("Four legs at D_LW %.2f e, R_DIP %.1f, GAP_ABS %.1f eV, R_GAP %.1f, against d31/ctl_R and" % (D_LW, R_DIP, GAP_ABS, R_GAP))
    say("d31/ctl_TS; stability on the last analysis (AMENDMENT 5). Fifth leg by AMENDMENT 6:")
    say("VDE = E(monoanion) - E(dianion), same state, same charges, %.12f eV per Eh." % EV)
    say("")

    # ------------------------------------------------------------ G43.6: the committed files
    psrc = os.path.join(D, "SOURCES.sha256")
    srcs = []
    if not os.path.isfile(psrc):
        gates.append("G43.6 d43/SOURCES.sha256 is missing")
    else:
        for l in open(psrc).read().splitlines():
            if l.strip():
                h, f = l.split(None, 1)
                srcs.append(f)
                p = os.path.join(S6, f)
                if not os.path.isfile(p) or sha(p) != h:
                    gates.append("G43.6 %s is not the file d43/SOURCES.sha256 lists" % f)
        if srcs != list(prep.SOURCES):
            gates.append("G43.6 d43/SOURCES.sha256 does not list the files the selection reads")
    rc, _ = git(a.res, "diff", "--quiet", "HEAD", "--", *(["%s/%s" % (S, x) for x in GUARDED] + ["%s/%s" % (S, x) for x in prep.SOURCES]))
    if rc is None:
        gates.append("G43.6 could not be checked: git did not run")
    elif rc != 0:
        gates.append("G43.6 a committed file under %s, or one the selection reads, differs from HEAD" % ", ".join(GUARDED))
    committed = ["d31/ctl_R", "d31/ctl_TS", "d40/B_bare_R", "d40/B_bare_TS"]
    rc, pcm = git(a.res, "log", "--diff-filter=A", "--format=%H", "-1", "--", "%s/%s" % (S, PREREG))
    if rc != 0 or not pcm:
        gates.append("G43.6 the commit that added %s could not be found" % PREREG)
        pcm = ""
    else:
        rc1, b1 = git(a.res, "rev-parse", "%s:%s/%s" % (pcm, S, PREREG))
        rc2, b2 = git(a.res, "hash-object", a.prereg)
        if rc1 != 0 or rc2 != 0 or b1 != b2:
            gates.append("G43.6 the pre-registration given is not the committed %s" % PREREG)
        rc1, b1 = git(a.res, "rev-parse", "%s:%s/%s" % (pcm, S, AMEND6))
        rc2, b2 = git(a.res, "hash-object", pam)
        if rc1 != 0 or rc2 != 0 or b1 != b2:
            gates.append("G43.6 %s is not the file at the commit the pre-registration is committed on" % AMEND6)
        for s in committed:
            rc1, b1 = git(a.res, "rev-parse", "%s:%s/%s/job.out" % (pcm, S, s))
            rc2, b2 = git(a.res, "hash-object", os.path.join(S6, s, "job.out"))
            if rc1 != 0 or rc2 != 0 or b1 != b2:
                gates.append("G43.6 %s/job.out is not the file at the commit the pre-registration is committed on" % s)

    # ------------------------------------------------------------ the sentences this result quotes
    def has(path, sent):
        return os.path.isfile(path) and sent in " ".join(open(path).read().split())

    if not has(os.path.join(S6, P2_FILE), P2_SENT):
        gates.append("R43.4 Probe 2's sentence is not in %s" % P2_FILE)
    if not has(os.path.join(S6, D11P), D11_Q):
        gates.append("R43.8 the words quoted from D11.4, \"%s\", are not in %s" % (D11_Q, D11P))
    for q in (H3_A, H3_B, H2_A, H2_B, H9_A, H9_B):
        if not has(a.prereg, q):
            gates.append("the words of a hypothesis, \"%s\", are not in the pre-registration given" % q)

    # ------------------------------------------------------------ G43.2: the selection, repeated
    sel, cond = None, []
    try:
        sel = prep.select(S6)
    except prep.Refused as e:
        gates.append("G43.2 the selection could not be repeated from the committed files: %s" % e)
    except Exception as e:                                   # a file that does not parse is a gate, not a crash
        gates.append("G43.2 the selection could not be repeated from the committed files: %s: %s" % (type(e).__name__, e))
    ctsv = []
    pc_ = os.path.join(D, "CONDITIONS.tsv")
    if not os.path.isfile(pc_):
        gates.append("G43.2 d43/CONDITIONS.tsv is missing")
    else:
        ls = open(pc_).read().splitlines()
        hdr = ls[0].split("\t") if ls else []
        ctsv = [dict(zip(hdr, l.split("\t"))) for l in ls[1:] if l.strip()]
    if sel is not None and ctsv:
        if [c["condition"] for c in ctsv if "condition" in c] != [c[0] for c in sel["cond"]]:
            gates.append("G43.2 CONDITIONS.tsv does not list the %d conditions the selection gives, in its order" % len(sel["cond"]))
        else:
            for row, (name, part, why, ss, qs) in zip(ctsv, sel["cond"]):
                try:
                    fo = math.fsum(sel["a"][i] * v for i, v in zip(ss, qs))
                    Ls = math.fsum(sel["L"][i] for i in ss)
                    s = sum(qs)
                    ok = (row["part"] == part and row["rule"] == why and int(row["n"]) == len(ss)
                          and row["sites"] == ",".join(str(i) for i in ss)
                          and row["charges"] == ",".join("%+.4f" % v for v in qs)
                          and abs(float(row["q_e"]) - abs(qs[0])) < 1e-9
                          and abs(float(row["sum_q_e"]) - s) < 5e-5 and abs(float(row["sum_abs_q_e"]) - sum(abs(v) for v in qs)) < 5e-5
                          and row["class"] == ("H" if s > -5e-7 else "N")
                          and int(row["n_positive"]) == sum(1 for v in qs if v > 0) and int(row["n_negative"]) == sum(1 for v in qs if v < 0)
                          and abs(float(row["first_order_kcal"]) - fo) < 2e-6 and abs(float(row["L_kcal"]) - Ls) < 2e-6
                          and abs(float(row["score_kcal"]) - (fo + Ls)) < 2e-6
                          and abs(float(row["nearest_A_R"]) - min(sel["near"]["R"][i][0] for i in ss)) < 1e-4
                          and abs(float(row["nearest_A_TS"]) - min(sel["near"]["TS"][i][0] for i in ss)) < 1e-4)
                except (KeyError, ValueError):
                    ok = False
                if not ok:
                    gates.append("G43.2 CONDITIONS.tsv does not describe %s as the selection gives it" % name)
                    continue
                cond.append(dict(name=name, part=part, why=why, sites=ss, q=qs, n=len(ss), mag=abs(qs[0]), sum=s,
                                 total=sum(abs(v) for v in qs), cls=row["class"], fo=fo, L=Ls, score=fo + Ls,
                                 nearR=min(sel["near"]["R"][i][0] for i in ss), nearTS=min(sel["near"]["TS"][i][0] for i in ss)))
    if not gates and len(cond) != NCOND:
        gates.append("G43.2 %d conditions were read, not %d" % (len(cond), NCOND))
    pst = os.path.join(D, "SITES.tsv")
    if not os.path.isfile(pst):
        gates.append("G43.2 d43/SITES.tsv is missing")
    elif sel is not None:
        ls = open(pst).read().splitlines()
        hdr = ls[0].split("\t") if ls else []
        srows = [dict(zip(hdr, l.split("\t"))) for l in ls[1:] if l.strip()]
        used = sorted(set(i for c in sel["cond"] for i in c[3]))
        try:
            ok = [int(r["site"]) for r in srows] == used
            for r in srows if ok else []:
                i = int(r["site"])
                ok = ok and (r["layer"] == sel["layer"][i] and abs(float(r["a_kcal_per_e"]) - sel["a"][i]) < 2e-6 and abs(float(r["L_kcal"]) - sel["L"][i]) < 2e-6
                             and max(abs(float(r[k]) - sel["xyz"][i][j]) for j, k in enumerate(("x", "y", "z"))) < 1e-9
                             and abs(float(r["nearest_A_R"]) - sel["near"]["R"][i][0]) < 1e-4 and r["nearest_atom_R"] == sel["near"]["R"][i][1]
                             and abs(float(r["nearest_A_TS"]) - sel["near"]["TS"][i][0]) < 1e-4 and r["nearest_atom_TS"] == sel["near"]["TS"][i][1]
                             and r["used_in"] == ",".join(c[0] for c in sel["cond"] if i in c[3]))
        except (KeyError, ValueError):
            ok = False
        if not ok:
            gates.append("G43.2 SITES.tsv does not describe the sites the selection uses")

    # ------------------------------------------------------------ G43.1: every committed input
    pic = os.path.join(D, "INPUT_CHECKSUMS.sha256")
    nsum = 0
    if not os.path.isfile(pic):
        gates.append("G43.1 d43/INPUT_CHECKSUMS.sha256 is missing")
    else:
        for l in open(pic).read().splitlines():
            if l.strip():
                h, f = l.split(None, 1)
                nsum += 1
                p = os.path.join(D, f)
                if not os.path.isfile(p) or sha(p) != h:
                    gates.append("G43.1 d43/%s is not the file INPUT_CHECKSUMS.sha256 lists" % f)
        if nsum != NSUM:
            gates.append("G43.1 INPUT_CHECKSUMS.sha256 lists %d files, not %d" % (nsum, NSUM))

    # ------------------------------------------------------------ every run; G43.1 and G43.3
    runs = []
    if len(cond) == NCOND:
        try:
            runs = prep.runs_of([(c["name"], c["part"], c["why"], c["sites"], c["q"]) for c in cond])
        except prep.Refused as e:
            gates.append("G43.2 %s" % e)
    cmap = dict((c["name"], c) for c in cond)
    tpl = os.path.join(S6, prep.TEMPLATE)
    tl = open(tpl).read().split("\n") if os.path.isfile(tpl) else []
    atoms = {}
    for st in STATES:
        try:
            atoms[st] = prep.read_geom(os.path.join(S6, prep.GEOM[st]))
        except Exception:
            gates.append("G43.1 the committed geometry %s could not be read" % prep.GEOM[st])
    nfile, ncharge = 0, 0
    for run, src, st, sp_, cm in runs:
        d = os.path.join(D, run)
        x, y = os.path.join(d, "substrate.xyz"), os.path.join(S6, prep.GEOM[st])
        if not (os.path.isfile(x) and os.path.isfile(y)) or sha(x) != sha(y):
            gates.append("G43.1 %s/substrate.xyz is not %s" % (run, prep.GEOM[st]))
        nfile += 1
        pi = os.path.join(d, "job.inp")
        if not os.path.isfile(pi):
            gates.append("G43.1 %s has no input" % run)
            continue
        nfile += 1
        if sp_ is None:
            ps = os.path.join(S6, src, "job.inp")
            if not os.path.isfile(ps) or sha(pi) != sha(ps):
                gates.append("G43.1 the control input %s is not byte-identical to %s/job.inp" % (run, src))
            continue
        li = open(pi).read().split("\n")
        if sp_ == "di":
            if sha(pi) != sha(tpl):
                gates.append("G43.1 %s/job.inp is not byte-identical to %s" % (run, prep.TEMPLATE))
        else:
            diff = [i for i in range(min(len(li), len(tl))) if li[i] != tl[i]]
            if len(li) != len(tl) or len(diff) != 1 \
                    or tl[diff[0]] != "* xyzfile -2 1 substrate.xyz" or li[diff[0]] != "* xyzfile -1 2 substrate.xyz":
                gates.append("G43.1 %s does not differ from %s in the charge-and-multiplicity line alone" % (run, prep.TEMPLATE))
        pp = os.path.join(d, "pc.pc")
        first = os.path.join(D, prep.rname(src, "di", "R"), "pc.pc")
        if not (os.path.isfile(pp) and os.path.isfile(first)) or sha(pp) != sha(first):
            gates.append("G43.1 %s/pc.pc is not byte-identical to the other point-charge files of %s" % (run, src))
            continue
        nfile += 1
        c = cmap[src]
        try:
            pc = prep.read_pc(pp)
        except Exception:
            gates.append("G43.2 %s/pc.pc could not be read" % run)
            continue
        if len(pc) != c["n"]:
            gates.append("G43.2 %s/pc.pc does not hold %d charges" % (run, c["n"]))
            continue
        for (q, px, py, pz), i, v in zip(pc, c["sites"], c["q"]):
            sx = sel["xyz"][i]
            if abs(q - v) > 1e-9 or max(abs(px - sx[0]), abs(py - sx[1]), abs(pz - sx[2])) > 1e-9:
                gates.append("G43.2 %s/pc.pc does not hold site %d at %+.2f e" % (run, i, v))
            elif not prep.admits(sel["ranges"], i, q):
                gates.append("G43.2 %s/pc.pc holds %+.2f e at site %d, which does not admit it" % (run, q, i))
            if len(atoms) == 2 and min(prep.nearest((px, py, pz), atoms[s2])[0] for s2 in STATES) < prep.R_SAFE:
                gates.append("G43.3 %s/pc.pc holds a charge inside %.2f A" % (run, prep.R_SAFE))
            ncharge += 1
        if c["part"] == "B" and sp_ == "di" and st == "R":
            qs = [r[0] for r in pc]
            if (sum(qs) < -5e-7 or sum(1 for v in qs if v < 0) > len(qs) // 2 or sum(abs(v) for v in qs) > prep.CONV + 1e-9
                    or len(set(round(abs(v), 9) for v in qs)) != 1
                    or any(prep.dist(pc[i][1:], pc[j][1:]) < prep.R_MIN for i in range(len(pc)) for j in range(i))):
                gates.append("G43.3 %s does not obey the rules of part B" % src)

    # ------------------------------------------------------------ the committed references
    ref = {}
    p40t = os.path.join(S6, "d40.tsv")
    if os.path.isfile(p40t):
        l40 = open(p40t).read().splitlines()
        h40 = l40[0].split("\t")
        for l in l40[1:]:
            r = dict(zip(h40, l.split("\t")))
            ref[r.get("run")] = r
    if ref.get("B_bare_R", {}).get("E_mono_Eh") != "%.8f" % CTL_MO_E:
        gates.append("G43.4 the committed source of %.8f Eh, d40.tsv, does not print it for B_bare_R" % CTL_MO_E)
    if ref.get("B_bare_R", {}).get("E_di_Eh") != "%.8f" % CTL_DI_E:
        gates.append("G43.4 the committed source of %.8f Eh, d40.tsv, does not print it as the dianion of B_bare_R" % CTL_DI_E)
    for st in STATES:
        if ref.get("B_bare_%s" % st, {}).get("VDE_eV") != BARE[st]:
            gates.append("G43.5 the committed source of the bare substrate's %s eV at %s, d40.tsv, does not print it" % (BARE[st], st))
        p = os.path.join(S6, "d31", "ctl_%s" % st, "job.out")
        if not os.path.isfile(p) or sha(p) != CTL_SHA[st]:
            gates.append("G43.5 d31/ctl_%s/job.out is not the file D33 pinned" % st)

    # ------------------------------------------------------------ read everything
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

    CT = {}
    for st in STATES:
        CT[st] = read_di(os.path.join(S6, "d31", "ctl_%s" % st), "plain")
        if CT[st]["E"] is None:
            gates.append("G43.5 the committed control d31/ctl_%s is not readable: %s" % (st, "; ".join(CT[st]["why"])))
        elif CT[st]["c"].get("gap") is None or CT[st]["trah"] or CT[st]["lshift"] or CT[st]["stab"].get("last") != "STABLE":
            gates.append("G43.5 the committed control d31/ctl_%s does not carry a readable gap and a STABLE last analysis" % st)
    BM = {}
    for st in STATES:
        BM[st] = d40.read(os.path.join(S6, "d40", "B_bare_%s" % st), "plain", NEL_MONO, d33b)
        if BM[st]["E"] is None:
            gates.append("G43.5 the committed monoanion d40/B_bare_%s is not readable: %s" % (st, "; ".join(BM[st]["why"])))
    c_di = read_di(os.path.join(D, prep.CTL_DI[0]), "plain")
    c_mo = d40.read(os.path.join(D, prep.CTL_MO[0]), "plain", NEL_MONO, d33b)
    for nm, r, e0 in ((prep.CTL_DI[0], c_di, CTL_DI_E), (prep.CTL_MO[0], c_mo, CTL_MO_E)):
        if r["E"] is None:
            gates.append("G43.4 the control %s is not readable: %s" % (nm, "; ".join(r["why"])))
        elif abs(r["E"] - e0) > TOL_EH:
            gates.append("G43.4 the control %s gives %.8f Eh, which differs from the committed %.8f by %.2e, over the %.0e limit"
                         % (nm, r["E"], e0, abs(r["E"] - e0), TOL_EH))
    bare, barrier = {}, None
    if all(CT[st]["E"] is not None and BM[st]["E"] is not None for st in STATES):
        for st in STATES:
            bare[st] = (BM[st]["E"] - CT[st]["E"]) * EV
            if "%+.4f" % bare[st] != BARE[st]:
                gates.append("G43.5 the bare substrate's detachment energy at %s, recomputed from the committed outputs, is %+.4f eV and not D40's %s"
                             % (st, bare[st], BARE[st]))
        barrier = (CT["TS"]["E"] - CT["R"]["E"]) * H2K
        if abs(barrier - D31_BARRIER) > TOL_K:
            gates.append("G43.5 the charge-free barrier from d31/ctl_R and d31/ctl_TS is %.4f kcal/mol, not %.4f within %.0e"
                         % (barrier, D31_BARRIER, TOL_K))
        if c_di["E"] is not None and CT["R"]["E"] is not None and "%.8f" % CT["R"]["E"] != "%.8f" % CTL_DI_E:
            gates.append("G43.4 d31/ctl_R/job.out does not print %.8f Eh" % CTL_DI_E)

    X = {}
    for c in cond:
        for st in STATES:
            di = read_di(os.path.join(D, prep.rname(c["name"], "di", st)), "pc")
            mo = d40.read(os.path.join(D, prep.rname(c["name"], "mo", st)), "pc", NEL_MONO, d33b)
            x = dict(di=di, mo=mo, flags=[], held=[], notes=[], rise=None, ratio=None, gap=None, gap_read=False)
            if di["status"] == "DID NOT CONVERGE":
                x["flags"].append("did not converge")
            elif di["status"] == "measured" and CT[st]["c"] is not None:
                g, k = di["c"], CT[st]["c"]
                if g["nb"] != k["nb"]:
                    gates.append("G43.7 %s reports %s basis functions against the control's %s" % (prep.rname(c["name"], "di", st), g["nb"], k["nb"]))
                x["rise"] = g["maxlw"] - k["maxlw"]
                x["ratio"] = g["dip"] / k["dip"] if k["dip"] else None
                if x["rise"] >= D_LW:
                    x["flags"].append("relocation: Loewdin rise %+.4f e" % x["rise"])
                if x["ratio"] is None:
                    x["held"].append("the dipole ratio could not be formed")
                elif x["ratio"] >= R_DIP:
                    x["flags"].append("dipole %.2fx control" % x["ratio"])
                x["gap"] = g.get("gap")
                if di["trah"] or di["lshift"] or g.get("gap") is None or k.get("gap") is None:
                    x["held"].append("gap withheld")
                else:
                    x["gap_read"] = True
                    if g["gap"] < GAP_ABS and g["gap"] < R_GAP * k["gap"]:
                        x["flags"].append("gap %.4f eV, %.2f of control" % (g["gap"], g["gap"] / k["gap"]))
                sb = di["stab"]
                if not sb["ran"] or sb["first"] is None or sb["last"] is None:
                    x["held"].append("stability analysis not read")
                elif sb["last"] == "UNSTABLE":
                    x["flags"].append("unstable at the last analysis")
                elif sb["first"] == "UNSTABLE":
                    x["notes"].append("the first solution was unstable%s and was corrected; the last analysis is STABLE (AMENDMENT 5: not a flag)" % (
                        "" if sb["nneg"] is None else " (%d negative, lowest %+.6f Eh)" % (sb["nneg"], sb["lo"])))
            x["v"] = None if (di["E"] is None or mo["E"] is None) else (mo["E"] - di["E"]) * EV
            x["word"] = word_of(x["v"], mo.get("stab") or {})
            X[(c["name"], st)] = x

    if gates:
        say("GATES")
        for g in gates:
            say("   FAILED  %s" % g)
        say("")
        say("A FAILURE OF G43.1 TO G43.7 STOPS THE DIAGNOSTIC. Nothing else is readable.")
        say("The reader stops in the same way where a committed file it must quote or read is not")
        say("as the pre-registration names it.")
        print("\n".join(rep))
        stop("%d gate(s) failed. Nothing is written." % len(gates))

    # ------------------------------------------------------------ the words of each condition
    def verdict(c):
        xs = [X[(c["name"], st)] for st in STATES]
        if any(x["flags"] for x in xs):
            return "UNRELIABLE"
        if any(x["di"]["status"] != "measured" for x in xs):
            return "NOT CLASSIFIED"
        if any(x["held"] for x in xs):
            return "INCONCLUSIVE"
        return "RELIABLE"

    def fifth(c):
        w = [X[(c["name"], st)]["word"] for st in STATES]
        if w == ["BOUND", "BOUND"]:
            return "PASSED"
        if "METASTABLE" in w:
            return "NOT PASSED: METASTABLE"
        if None in w:
            return "NOT PASSED: NOT MEASURED"
        return "NOT PASSED: NOT STABLE"

    def probe2(c):
        w = [X[(c["name"], st)]["word"] for st in STATES]
        if "METASTABLE" in w:
            return "APPLIES"
        if w == ["BOUND", "BOUND"]:
            return "does not reach"
        if None in w:
            return "not said: %s NOT MEASURED" % " and ".join(st for st in STATES if X[(c["name"], st)]["word"] is None)
        return "not said: monoanion not stable"

    def cellw(c, st):
        w = X[(c["name"], st)]["word"]
        return "NOT MEASURED" if w is None else ("ABOVE ZERO *" if w == NOTSTABLE else w)

    def fs(v, fmt="%.6f", none="-"):
        return none if v is None else fmt % v

    def wrapsay(text, lead="   ", width=76):
        for l in textwrap.wrap(text, width):
            say("%s%s" % (lead, l))

    for c in cond:
        c["verdict"], c["fifth"], c["probe2"] = verdict(c), fifth(c), probe2(c)
        e = dict((st, X[(c["name"], st)]["di"]["E"]) for st in STATES)
        c["dde"] = None if None in e.values() else ((e["TS"] - CT["TS"]["E"]) - (e["R"] - CT["R"]["E"])) * H2K
    PB = [c for c in cond if c["part"] == "B"]
    PA = [c for c in cond if c["part"] == "A"]
    PAp = [c for c in PA if c["sum"] > 0]
    PAm = [c for c in PA if c["sum"] < 0]
    if len(PB) != 6 or len(PAp) != 6 or len(PAm) != 6:
        stop("the reader's own lists hold %d designs, %d positive and %d negative single charges, not 6 of each" % (len(PB), len(PAp), len(PAm)))

    say("=" * 78)
    say("G43.1 TO G43.7  THE STOPPING GATES")
    say("=" * 78)
    say("")
    say("   G43.1  all %d files of INPUT_CHECKSUMS.sha256 are the ones listed. Of them, %d" % (nsum, nfile))
    say("          geometry, input and point-charge files: every geometry is d33a's, every")
    say("          dianion input is d31/B/n01_q0.25_R/job.inp byte for byte, every monoanion input")
    say("          differs from it in the charge-and-multiplicity line alone, the two control")
    say("          inputs are their committed sources, and the four point-charge files of a")
    say("          condition are identical.")
    say("   G43.2  the selection, repeated here from the committed files, gives the %d conditions" % NCOND)
    say("          of CONDITIONS.tsv and the sites of SITES.tsv: D33 stage B's six committed designs")
    say("          come back from the rule,")
    say("          D12's sites at or beyond %.2f A are %s, the sites of most" % (prep.R_SAFE, ", ".join(str(i) for i in prep.D12_OUTER)))
    say("          negative and most positive a_i are %d and %d, and the six designs are the" % (prep.ANCHOR_POS, prep.ANCHOR_NEG))
    say("          rule's. Every point-charge file holds the sites and charges listed.")
    say("   G43.3  all %d charges of the %d point-charge files lie at least %.2f A from every" % (ncharge, 4 * NCOND, prep.R_SAFE))
    say("          substrate atom at both states; every design of part B obeys its rules.")
    say("   G43.6  the %d files of SOURCES.sha256 are the ones listed; no committed file under" % len(srcs))
    say("          %s or %s differs from HEAD; the pre-registration, %s" % (", ".join(GUARDED[:-1]), GUARDED[-1], AMEND6))
    say("          and the %d committed outputs read are the files at the commit the" % len(committed))
    say("          pre-registration was committed on, %s." % pcm[:8])
    say("   G43.7  every dianion run that was read reports the control's %s basis functions." % CT["R"]["c"]["nb"])
    say("   G43.4 and G43.5 are the controls, under R43.1.")
    say("")

    # ------------------------------------------------------------ R43.1
    say("=" * 78)
    say("R43.1  THE CONTROLS")
    say("=" * 78)
    say("")
    say("   G43.4  C_ctl_R, D31's committed charge-free dianion input, run again:")
    say("          %.8f Eh   committed %.8f   differ %.2e   limit %.0e   PASS" % (c_di["E"], CTL_DI_E, abs(c_di["E"] - CTL_DI_E), TOL_EH))
    say("          C_bare_R, D40's committed monoanion input for the bare substrate, run again:")
    say("          %.8f Eh   committed %.8f   differ %.2e   limit %.0e   PASS" % (c_mo["E"], CTL_MO_E, abs(c_mo["E"] - CTL_MO_E), TOL_EH))
    say("   G43.5  d31/ctl_R and d31/ctl_TS are the files D33 pinned. The charge-free barrier:")
    say("          %.4f kcal/mol   D31 %.4f   differ %.2e   limit %.0e   PASS" % (barrier, D31_BARRIER, abs(barrier - D31_BARRIER), TOL_K))
    say("          The bare substrate's detachment energy, from the committed outputs:")
    say("          %+.4f eV at the reactant and %+.4f eV at the transition state, as D40 has them." % (bare["R"], bare["TS"]))
    say("")

    # ------------------------------------------------------------ R43.2
    say("=" * 78)
    say("R43.2  THE FOUR LEGS, ON THE DIANION")
    say("=" * 78)
    say("")
    say("   rise: the largest Loewdin charge minus the control's, e; flagged at %.2f." % D_LW)
    say("   dip: the dipole over the control's; flagged at %.1f." % R_DIP)
    say("   gap: eV; flagged below %.1f and below %.1f of the control's (%.4f at R, %.4f at TS);" % (
        GAP_ABS, R_GAP, CT["R"]["c"]["gap"], CT["TS"]["c"]["gap"]))
    say("        'withheld' where TRAH or a level shift was engaged.")
    say("   stability: the dianion's first analysis, then its last. read: the number of point")
    say("   charges ORCA says it read.")
    say("")
    say("   condition     st    rise    dip    gap        stability          TRAH shift read  flags")
    tail = []
    for group, title in ((PB, "part B, the designs"), (PAp, "part A, +1.00 e"), (PAm, "part A, -1.00 e")):
        say("   %s" % title)
        for c in group:
            for st in STATES:
                x = X[(c["name"], st)]
                di = x["di"]
                if di["status"] != "measured":
                    say("   %-13s %-2s  %s" % (c["name"], st, "DID NOT CONVERGE: a flag (D31.3)" if di["status"] == "DID NOT CONVERGE" else "NOT MEASURED (G43.8)"))
                    for w in di["why"]:
                        say("        %s" % w)
                    continue
                sb = di["stab"]
                cell = "NOT FOUND" if not sb["ran"] else ("NOT READ" if sb["first"] is None else "%s/%s" % (sb["first"], sb["last"]))
                mm = re.search(r"Reading point charge file\s*\.*\s*ok \((\d+) point charges?\)", di["text"])
                x["npc"] = int(mm.group(1)) if mm else None
                say("   %-13s %-2s  %+.4f  %5.2f  %-9s  %-17s  %-4s %-5s %-4s  %s" % (
                    c["name"], st, x["rise"], x["ratio"] if x["ratio"] is not None else float("nan"),
                    ("%.4f" % x["gap"]) if x["gap_read"] else "withheld", cell, "yes" if di["trah"] else "no",
                    "yes" if di["lshift"] else "no", "-" if x["npc"] is None else str(x["npc"]), "; ".join(x["flags"]) or "none"))
                for n_ in x["notes"]:
                    tail.append("%s %s: %s." % (c["name"], st, n_))
                if not sb["ran"]:
                    tail.append("%s %s: no stability analysis was found in the output, though its input asks for one. Lines of the output that mention stability: %s" % (
                        c["name"], st, " | ".join(sb["quote"]) or "none"))
                elif sb["first"] is None:
                    tail.append("%s %s: %s. ORCA's own lines: %s" % (c["name"], st, sb["how"], " | ".join(sb["quote"]) or "(none found)"))
                elif sb["how"] != "cs6_d33b.stability":
                    tail.append("%s %s: stability read from %s." % (c["name"], st, sb["how"]))
    say("")
    for l in tail:
        for i, w in enumerate(textwrap.wrap(l, 80)):
            say("   %s%s" % ("NOTE  " if i == 0 else "      ", w))
    if tail:
        say("")
    say("   condition      sum e   total e  nearest atom A   on four legs")
    for c in cond:
        say("   %-13s  %+.2f    %.2f     %.3f / %.3f    %s" % (c["name"], c["sum"], c["total"], c["nearR"], c["nearTS"], c["verdict"]))
    say("")
    say("   RELIABLE: no flag at either state, with the gap readable and the stability analysis")
    say("   read at both. UNRELIABLE: a flag at either state. INCONCLUSIVE: no flag, and a gap")
    say("   withheld or a stability analysis not read. NOT CLASSIFIED: no flag, and a state NOT")
    say("   MEASURED. nearest atom: the smallest distance from a charge to a substrate atom, at")
    say("   the reactant and at the transition state.")
    say("")

    # ------------------------------------------------------------ R43.3
    say("=" * 78)
    say("R43.3  THE FIFTH LEG: THE DETACHMENT ENERGY WITH THE CHARGES PRESENT")
    say("=" * 78)
    say("")
    say("   BOUND: VDE above zero and the monoanion's last stability analysis STABLE.")
    say("   METASTABLE: VDE at or below zero. change: VDE minus the bare substrate's at the same")
    say("   state, %s eV at R and %s eV at TS. Stability is the monoanion's, first analysis" % (BARE["R"], BARE["TS"]))
    say("   then last.")
    say("")
    say("   condition     st    VDE eV    change   <S**2>    stability          TRAH shift read  word")
    tail = []
    for group, title in ((PB, "part B, the designs"), (PAp, "part A, +1.00 e"), (PAm, "part A, -1.00 e")):
        say("   %s" % title)
        for c in group:
            for st in STATES:
                x = X[(c["name"], st)]
                m, di = x["mo"], x["di"]
                if x["v"] is None:
                    say("   %-13s %-2s  NOT MEASURED (G43.8, G43.9)" % (c["name"], st))
                    if m["E"] is None:
                        for w in m["why"]:
                            say("        the monoanion %s" % w)
                    if di["E"] is None:
                        say("        the dianion of this state is not read (R43.2)")
                    continue
                sb = m["stab"]
                cell = "NOT FOUND" if not sb["ran"] else ("NOT READ" if sb["first"] is None else "%s/%s" % (sb["first"], sb["last"]))
                mm = re.search(r"Reading point charge file\s*\.*\s*ok \((\d+) point charges?\)", m["text"])
                x["npc_mo"] = int(mm.group(1)) if mm else None
                say("   %-13s %-2s  %+8.4f  %+8.4f  %s  %-17s  %-4s %-5s %-4s  %s" % (
                    c["name"], st, x["v"], x["v"] - bare[st], fs(m["s2"]), cell, "yes" if m["trah"] else "no",
                    "yes" if m["lshift"] else "no", "-" if x["npc_mo"] is None else str(x["npc_mo"]), x["word"]))
                if not sb["ran"]:
                    tail.append("%s %s: no stability analysis was found in the monoanion's output, though its input asks for one. Lines of the output that mention stability: %s" % (
                        c["name"], st, " | ".join(sb["quote"]) or "none"))
                elif sb["first"] is None:
                    tail.append("%s %s: %s. ORCA's own lines: %s" % (c["name"], st, sb["how"], " | ".join(sb["quote"]) or "(none found)"))
                else:
                    if sb["how"] != "cs6_d33b.stability":
                        tail.append("%s %s: the monoanion's stability read from %s." % (c["name"], st, sb["how"]))
                    if sb["first"] == "UNSTABLE":
                        tail.append("%s %s: the monoanion's first solution was unstable%s; %s." % (
                            c["name"], st, "" if sb["nneg"] is None else " (%d negative, lowest %+.6f Eh)" % (sb["nneg"], sb["lo"]),
                            "1 analysis was made" if sb["n"] == 1 else "%d analyses were made" % sb["n"]))
                if x["word"] == NOTSTABLE:
                    tail.append("%s %s: the value is above zero and the monoanion's last stability analysis is not STABLE, so by AMENDMENT 6 it is not counted as BOUND." % (c["name"], st))
    say("")
    for l in tail:
        for i, w in enumerate(textwrap.wrap(l, 80)):
            say("   %s%s" % ("NOTE  " if i == 0 else "      ", w))
    if tail:
        say("")
    say("   condition      sum e    R             TS            the fifth leg")
    for c in cond:
        say("   %-13s  %+.2f    %-12s  %-12s  %s" % (c["name"], c["sum"], cellw(c, "R"), cellw(c, "TS"), c["fifth"]))
    if any(X[(c["name"], st)]["word"] == NOTSTABLE for c in cond for st in STATES):
        say("   * %s (AMENDMENT 6)." % NOTSTABLE)
    say("")
    say("   The change is reported, not interpreted. Each is formed from the energies and then")
    say("   rounded. D40.4: where a detachment energy is negative, its sign is the reading and its")
    say("   size is not. <S**2> is reported and not corrected; the ideal value is 0.75. Both")
    say("   energies of every row are in d43.tsv.")
    say("")

    # ------------------------------------------------------------ R43.4
    say("=" * 78)
    say("R43.4  THE BARRIER CHANGE BESIDE THE TWO TESTS")
    say("=" * 78)
    say("")
    say("   ddE = (E(TS) - E(TS, control)) - (E(R) - E(R, control)) on the dianion, kcal/mol.")
    say("   first: the first-order value, the sum of a_i q_i with the matrix of D33 stage A.")
    say("   score: first plus the Lennard-Jones term the committed objective carries; these bare")
    say("   charges carry none in the calculation.")
    say("")
    say("   Probe 2, %s:" % P2_FILE)
    for l in textwrap.wrap('"%s"' % P2_SENT, 72):
        say("      %s" % l)
    say("   A lowering is a difference between a transition state and a reactant, so the")
    say("   sentence applies to a condition's ddE if either of its states is METASTABLE (the")
    say("   rule of D40's R40.7).")
    say("")
    say("   condition          ddE     first     score  ddE-first  ratio  four legs       Probe 2's sentence")
    for c in cond:
        if c["dde"] is None:
            say("   %-13s  NOT FORMED: a dianion state is not read (G43.9)   first %+8.4f   %-14s  %s" % (
                c["name"], c["fo"], c["verdict"], c["probe2"]))
            continue
        say("   %-13s %+9.4f %+9.4f %+9.4f  %+8.4f  %5.2f  %-14s  %s" % (
            c["name"], c["dde"], c["fo"], c["score"], c["dde"] - c["fo"], c["dde"] / c["fo"], c["verdict"], c["probe2"]))
    say("")
    say("   APPLIES: a state is METASTABLE, so Probe 2's sentence applies to this ddE.")
    say("   does not reach: both states are BOUND, so this objection does not reach this ddE at")
    say("   this level and basis. ratio: ddE over first. Where the first-order value is small,")
    say("   as at sites 4 and 1043, the ratio says little and the difference is the reading.")
    say("   A ddE is printed for an UNRELIABLE condition, as D33's result did; the word beside")
    say("   it governs.")
    say("")

    def three(group, n, label):
        """FALSIFIED / NOT FALSIFIED BY THESE n / INCOMPLETE on a group read on four legs."""
        bad = [c for c in group if c["verdict"] == "UNRELIABLE"]
        open_ = [c for c in group if c["verdict"] in ("INCONCLUSIVE", "NOT CLASSIFIED")]
        if bad:
            say("   %s IS FALSIFIED." % label)
            for c in bad:
                for st in STATES:
                    x = X[(c["name"], st)]
                    if x["flags"]:
                        wrapsay("%s at %s: %s." % (c["name"], st, "; ".join(x["flags"])), "      ", 72)
            if open_:
                wrapsay("Not established, and not needed for the falsification: " + "; ".join("%s, %s" % (c["name"], c["verdict"]) for c in open_) + ".")
        elif open_:
            say("   THE TEST OF %s IS INCOMPLETE." % label)
            say("   None of the %d is UNRELIABLE, and these are not established:" % n)
            for c in open_:
                say("      %s: %s" % (c["name"], c["verdict"]))
        else:
            say("   %s IS NOT FALSIFIED BY THESE %d CONDITIONS." % (label, n))
        return bad, open_

    # ------------------------------------------------------------ R43.5
    say("=" * 78)
    say("R43.5  THE TEST OF H3, ON THE 12 SINGLE CHARGES OF PART A")
    say("=" * 78)
    say("")
    say("   H3, as section 1 of the pre-registration quotes it from the consolidated write-up:")
    for l in textwrap.wrap('"%s" To be given up if "%s"' % (H3_A, H3_B), 72):
        say("      %s" % l)
    say("")
    say("   site    a_i       nearest atom A    +1.00 e        -1.00 e")
    for cp, cm_ in zip(PAp, PAm):
        i = cp["sites"][0]
        say("   %-5d  %+.4f   %.3f / %.3f    %-13s  %s" % (i, sel["a"][i], cp["nearR"], cp["nearTS"], cp["verdict"], cm_["verdict"]))
    say("")
    say("   RELIABLE on four legs: %d of the 6 at +1.00 e and %d of the 6 at -1.00 e." % (
        sum(1 for c in PAp if c["verdict"] == "RELIABLE"), sum(1 for c in PAm if c["verdict"] == "RELIABLE")))
    say("")
    bad3, open3 = three(PA, 12, "H3")
    if not bad3 and not open3:
        say("   D43.4: not falsified is not established. Six sites of 616, at one magnitude, where H3")
        say("   speaks of every magnitude up to one.")
    say("")

    # ------------------------------------------------------------ R43.6
    say("=" * 78)
    say("R43.6  THE TEST OF H2, ON THE 6 SINGLE CHARGES OF -1.00 e, AND THE SIGN-SWAPPED PAIRS")
    say("=" * 78)
    say("")
    say("   H2, as section 1 of the pre-registration quotes it:")
    for l in textwrap.wrap('"%s" To be given up if "%s".' % (H2_A, H2_B), 72):
        say("      %s" % l)
    say("")

    def maxrise(c):
        r = [X[(c["name"], st)]["rise"] for st in STATES if X[(c["name"], st)]["rise"] is not None]
        return max(r) if r else None

    say("   site    +1.00 e: word, larger rise e      -1.00 e: word, larger rise e")
    for cp, cm_ in zip(PAp, PAm):
        say("   %-5d   %-14s %-9s          %-14s %s" % (cp["sites"][0], cp["verdict"], fs(maxrise(cp), "%+.4f"), cm_["verdict"], fs(maxrise(cm_), "%+.4f")))
    say("")
    bad2, open2 = three(PAm, 6, "H2")
    if not bad2 and not open2:
        say("   D43.4: not falsified is not established; these are six single charges.")
    say("   The pairs are reported. One site at two signs holds the position fixed and changes")
    say("   the sign alone, which the record did not have at this level.")
    say("")

    # ------------------------------------------------------------ R43.7
    say("=" * 78)
    say("R43.7  THE TEST OF H9, ON THE 12 CONDITIONS WHOSE CHARGES SUM TO ZERO OR MORE")
    say("=" * 78)
    say("")
    say("   H9, as section 1 of the pre-registration quotes it:")
    for l in textwrap.wrap('"%s". What would still falsify it: "%s".' % (H9_A, H9_B), 72):
        say("      %s" % l)
    say("   \"It\" is the dianion. H9 speaks of RELIABLE arrangements, so its test is on those of")
    say("   the 12 that are RELIABLE on four legs.")
    say("")
    HC = PB + PAp
    if any(c["cls"] != "H" for c in HC) or any(c["cls"] != "N" for c in PAm):
        stop("the classes of CONDITIONS.tsv are not H for the designs and the positive charges and N for the negative ones")
    say("   condition      sum e   on four legs    R             TS            in the test")
    for c in HC:
        say("   %-13s  %+.2f   %-14s  %-12s  %-12s  %s" % (c["name"], c["sum"], c["verdict"], cellw(c, "R"), cellw(c, "TS"),
                                                           "yes" if c["verdict"] == "RELIABLE" else "no: not RELIABLE"))
    if any(X[(c["name"], st)]["word"] == NOTSTABLE for c in HC for st in STATES):
        say("   * %s (AMENDMENT 6)." % NOTSTABLE)
    say("")
    T = [c for c in HC if c["verdict"] == "RELIABLE"]
    fals = [(c["name"], st) for c in T for st in STATES if X[(c["name"], st)]["word"] == "METASTABLE"]
    open9 = [(c["name"], st, "NOT MEASURED" if X[(c["name"], st)]["word"] is None else NOTSTABLE)
             for c in T for st in STATES if X[(c["name"], st)]["word"] in (None, NOTSTABLE)]
    say("   RELIABLE, and so in the test: %d of the 6 designs and %d of the 6 single charges." % (
        sum(1 for c in T if c["part"] == "B"), sum(1 for c in T if c["part"] == "A")))
    if not T:
        say("   H9 IS NOT TESTED. None of the 12 is RELIABLE on four legs.")
    elif fals:
        say("   H9 IS FALSIFIED.")
        wrapsay("RELIABLE on four legs and METASTABLE: " + "; ".join("%s at %s, %+.4f eV" % (n, st, X[(n, st)]["v"]) for n, st in fals) + ".")
        if open9:
            wrapsay("Not established, and not needed for the falsification: " + "; ".join("%s at %s, %s" % x for x in open9) + ".")
    elif open9:
        say("   THE TEST OF H9 IS INCOMPLETE.")
        say("   No RELIABLE condition is METASTABLE at either state, and these are not established:")
        for x in open9:
            say("      %s at %s: %s" % x)
    else:
        say("   H9 IS NOT FALSIFIED BY THESE %d CONDITIONS." % len(T))
        say("   D43.4: not falsified is not established. The six designs share two sites.")
    say("   The words of the conditions that are not RELIABLE, and of the 6 at -1.00 e, are in the")
    say("   table of R43.3. H9 neither predicts nor forbids either word for them.")
    say("")

    # ------------------------------------------------------------ R43.8
    say("=" * 78)
    say("R43.8  THE SURVIVORS OF THE SEARCH")
    say("=" * 78)
    say("")
    say("   A design of part B is a SURVIVOR if it is RELIABLE on four legs and passes the fifth.")
    say("   Every design has its total inside %.2f e and every charge at or beyond %.2f A (G43.3)." % (prep.CONV, prep.R_SAFE))
    say("")
    say("   design         n   sum e  total e   on four legs    the fifth leg            ddE kcal/mol")
    surv = []
    for c in PB:
        ok = c["verdict"] == "RELIABLE" and c["fifth"] == "PASSED"
        if ok:
            surv.append(c)
        say("   %-13s %2d   %+.2f   %.2f     %-14s  %-23s  %s%s" % (
            c["name"], c["n"], c["sum"], c["total"], c["verdict"], c["fifth"], fs(c["dde"], "%+9.4f", "NOT FORMED"), "   SURVIVOR" if ok else ""))
    say("")
    if surv:
        say("   THE SURVIVORS:")
        for c in surv:
            say("      %-13s  ddE %+9.4f kcal/mol   first order %+9.4f   %d charge%s, sum %+.2f e" % (
                c["name"], c["dde"], c["fo"], c["n"], "" if c["n"] == 1 else "s", c["sum"]))
        say("")
        say("   D43.1 and AMENDMENT 7: each ddE above is a working-level number and is quoted as no")
        say("   more than that. A survivor goes to the cross-check of AMENDMENT 7, under a")
        say("   pre-registration written before any cross-check output exists (F5).")
    else:
        say("   NO DESIGN OF PART B IS RELIABLE ON FOUR LEGS AND PASSES THE FIFTH.")
        rel = [c for c in PB if c["verdict"] == "RELIABLE"]
        if not rel:
            say("   None of the six is RELIABLE on four legs, which is the outcome section 6 of the")
            say("   pre-registration calls F3.")
        elif all(c["fifth"] == "NOT PASSED: METASTABLE" for c in rel):
            say("   %d of the six %s RELIABLE on four legs and %s METASTABLE at a state, which is the" % (
                len(rel), "is" if len(rel) == 1 else "are", "it is" if len(rel) == 1 else "each of them is"))
            say("   outcome section 6 of the pre-registration calls F4.")
        else:
            say("   %d of the six %s RELIABLE on four legs and none of them passes the fifth leg:" % (len(rel), "is" if len(rel) == 1 else "are"))
            for c in rel:
                say("      %-13s  %s" % (c["name"], c["fifth"]))
            say("   That is neither F3 nor F4 as section 6 of the pre-registration writes them: F4 is for")
            say("   RELIABLE designs that are all METASTABLE, and a fifth leg that is not established is")
            say("   not METASTABLE.")
    say("")
    singles = [c for c in PA if c["verdict"] == "RELIABLE" and c["fifth"] == "PASSED"]
    say("   Single charges of part A that are RELIABLE and pass the fifth leg, listed apart")
    say("   because D11.4 says \"%s\":" % D11_Q)
    if singles:
        for c in singles:
            say("      %-13s  ddE %s kcal/mol   first order %+9.4f" % (c["name"], fs(c["dde"], "%+9.4f", "NOT FORMED"), c["fo"]))
    else:
        say("      none.")
    say("")
    nbad = sum(1 for c in cond for st in STATES for r in (X[(c["name"], st)]["di"], X[(c["name"], st)]["mo"])
               if r.get("status", "measured" if r["E"] is not None else "NOT MEASURED") != "measured")
    if all(c["verdict"] == "UNRELIABLE" for c in cond):
        say("   F1 APPLIES: all 18 conditions are UNRELIABLE, the six with no positive charge among")
        say("   them. Nothing is concluded; the inputs and the reader are diagnosed first.")
        say("")
    if nbad > 18:
        say("   F2 APPLIES: %d of the 72 charged runs are NOT MEASURED or did not converge, which is" % nbad)
        say("   more than a quarter. That is the primary result.")
        say("")

    # ------------------------------------------------------------ R43.9
    say("=" * 78)
    say("R43.9  FIRST ORDER AT EACH SITE OF PART A")
    say("=" * 78)
    say("")
    say("   odd: half of ddE(+1.00 e) minus ddE(-1.00 e), the part of the response that changes")
    say("   sign with the charge, which is what a_i predicts. even: half their sum, the part that")
    say("   does not, which the first-order objective leaves out. kcal/mol.")
    say("")
    say("   site    a_i       ddE(+1.00)  ddE(-1.00)      odd      even    odd - a_i   a state METASTABLE")
    for cp, cm_ in zip(PAp, PAm):
        i = cp["sites"][0]
        met = [c["name"] for c in (cp, cm_) if c["probe2"] == "APPLIES"]
        if cp["dde"] is None or cm_["dde"] is None:
            say("   %-5d  %+.4f   %s  %s   NOT FORMED" % (i, sel["a"][i], fs(cp["dde"], "%+9.4f", "        -"), fs(cm_["dde"], "%+9.4f", "        -")))
            continue
        odd, even = (cp["dde"] - cm_["dde"]) / 2, (cp["dde"] + cm_["dde"]) / 2
        say("   %-5d  %+.4f   %+9.4f   %+9.4f   %+8.4f  %+8.4f   %+8.4f   %s" % (
            i, sel["a"][i], cp["dde"], cm_["dde"], odd, even, odd - sel["a"][i], ", ".join(met) if met else "neither"))
    say("")
    say("   Reported, not interpreted. Where a condition is named in the last column, Probe 2's")
    say("   sentence applies to its ddE, and so to both parts formed from it.")
    say("")

    # ------------------------------------------------------------ R43.10
    say("=" * 78)
    say("R43.10 ONE PAIR AT THREE MAGNITUDES, AND ITS TWO CHARGES ONE AT A TIME")
    say("=" * 78)
    say("")
    pairs = [c for c in PB if c["n"] == 2]
    say("   The pair: -q at site %d and +q at site %d. dVDE: the change in detachment energy" % (pairs[0]["sites"][0], pairs[0]["sites"][1]))
    say("   from the bare substrate.")
    say("")
    say("   design          q e    ddE kcal/mol   ddE / q    dVDE R eV   per e    dVDE TS eV   per e")

    def dv(c, st):
        v = X[(c["name"], st)]["v"]
        return None if v is None else v - bare[st]

    for c in pairs:
        q = c["mag"]
        say("   %-13s  %.2f   %s   %s   %s  %s   %s  %s" % (
            c["name"], q, fs(c["dde"], "%+10.4f", "         -"), fs(None if c["dde"] is None else c["dde"] / q, "%+8.4f", "       -"),
            fs(dv(c, "R"), "%+9.4f", "        -"), fs(None if dv(c, "R") is None else dv(c, "R") / q, "%+7.4f", "      -"),
            fs(dv(c, "TS"), "%+9.4f", "        -"), fs(None if dv(c, "TS") is None else dv(c, "TS") / q, "%+7.4f", "      -")))
    say("")
    p1 = [c for c in pairs if abs(c["mag"] - 1.0) < 1e-9][0]
    sm = cmap[prep.aname(p1["sites"][0], -1.0)]
    sp1 = cmap[prep.aname(p1["sites"][1], +1.0)]
    say("   The pair at 1.00 e beside its two charges one at a time, %s and %s:" % (sm["name"], sp1["name"]))
    say("                         the pair    sum of the two   pair minus sum")

    def addrow(label, a_, b_, c_, unit):
        if None in (a_, b_, c_):
            say("   %-20s  NOT FORMED: a state is not read" % label)
        else:
            say("   %-20s %+9.4f     %+9.4f       %+9.4f   %s" % (label, a_, b_ + c_, a_ - b_ - c_, unit))

    addrow("ddE", p1["dde"], sm["dde"], sp1["dde"], "kcal/mol")
    for st in STATES:
        addrow("dVDE at %s" % st, dv(p1, st), dv(sm, st), dv(sp1, st), "eV")
    say("")
    say("   Reported, not interpreted. Nothing is fitted. Where a state named here is METASTABLE")
    say("   (R43.3), its detachment energy has a sign and not a size, and its ddE is not a")
    say("   barrier effect.")
    say("")

    # ------------------------------------------------------------ limits
    say("=" * 78)
    say("WHAT THIS CANNOT SAY (section 7 of the pre-registration, in short; its own words govern)")
    say("=" * 78)
    say("")
    for l in ("D43.1  It makes no ddE quotable as more than a working-level number, withdraws nothing",
              "       and changes no verdict.",
              "D43.2  It does not re-measure r_safe. The 3.50 A is D11's, at another level, on one axis.",
              "D43.3  It does not take route (a) and does not adopt the wall. The grid is not rebuilt.",
              "D43.4  Not falsified is not established: six sites of 616 at one magnitude; six designs",
              "       that share two sites.",
              "D43.5  The six sites of part A are not a sample of the shell.",
              "D43.6  Part B is six designs by a greedy rule, not an optimisation.",
              "D43.7  One functional and one basis.",
              "D43.8  The detachment energies are vertical, in vacuo, unrestricted Kohn-Sham differences;",
              "       where one is negative its sign is the reading, not its size.",
              "D43.9  One frame, frozen geometries, in vacuo, no continuum.",
              "D43.10 The charges are bare point charges and stand for no residue. No wall is in any run.",
              "D43.11 It adopts no level of theory and says nothing about D39.",
              "D43.12 Nothing here is shaped: a result that falsifies a hypothesis is reported as plainly",
              "       as one that does not, by rules fixed before any output existed."):
        say("   %s" % l)
    say("")
    ndi = sum(1 for c in cond for st in STATES if X[(c["name"], st)]["di"]["status"] == "measured")
    nnc = sum(1 for c in cond for st in STATES if X[(c["name"], st)]["di"]["status"] == "DID NOT CONVERGE")
    nmo = sum(1 for c in cond for st in STATES if X[(c["name"], st)]["mo"]["E"] is not None)
    say("%d of the 36 dianion runs and %d of the 36 monoanion runs are measured.%s%s" % (
        ndi, nmo, " %d dianion run%s did not converge, which is a flag." % (nnc, "" if nnc == 1 else "s") if nnc else "",
        "" if ndi + nnc == 36 and nmo == 36 else " The rest are NOT MEASURED, each listed above with ORCA's own words."))

    # ------------------------------------------------------------ the table
    L = ["condition\tpart\tstate\tn\tq_e\tsum_q_e\tclass\tsites\tcharges\tnearest_A\t"
         "di_run\tdi_status\tE_di_Eh\tlowdin_rise_e\tdipole_ratio\tgap_eV\tgap_read\tdi_stab_first\tdi_stab_last\tdi_TRAH\tdi_level_shift\tdi_charges_read\tflags\t"
         "mo_run\tmo_status\tE_mo_Eh\tS2\tmo_stab_first\tmo_stab_last\tmo_TRAH\tmo_level_shift\tmo_charges_read\tVDE_eV\tword\tdVDE_vs_bare_eV\t"
         "four_legs\tfifth_leg\tddE_kcal\tfirst_order_kcal\tscore_kcal\tprobe2_sentence\tsurvivor"]

    def sbcell(sb, k):
        return sb.get(k) or ("NOT FOUND" if not sb.get("ran") else "NOT READ")

    for c in cond:
        for st in STATES:
            x = X[(c["name"], st)]
            di, m = x["di"], x["mo"]
            f = [c["name"], c["part"], st, "%d" % c["n"], "%.2f" % c["mag"], "%+.4f" % c["sum"], c["cls"],
                 ",".join(str(i) for i in c["sites"]), ",".join("%+.4f" % v for v in c["q"]), "%.4f" % (c["nearR"] if st == "R" else c["nearTS"]),
                 prep.rname(c["name"], "di", st), di["status"]]
            if di["status"] == "measured":
                dsb = di["stab"]
                f += ["%.8f" % di["E"], "%+.4f" % x["rise"], "" if x["ratio"] is None else "%.4f" % x["ratio"],
                      "" if x["gap"] is None else "%.4f" % x["gap"], "yes" if x["gap_read"] else "no", sbcell(dsb, "first"), sbcell(dsb, "last"),
                      "yes" if di["trah"] else "no", "yes" if di["lshift"] else "no", "" if x.get("npc") is None else str(x["npc"])]
            else:
                f += [""] * 10
            f += ["; ".join(x["flags"]), prep.rname(c["name"], "mo", st), "measured" if m["E"] is not None else "NOT MEASURED"]
            if m["E"] is not None:
                msb = m["stab"]
                f += ["%.8f" % m["E"], "" if m["s2"] is None else "%.6f" % m["s2"], sbcell(msb, "first"), sbcell(msb, "last"),
                      "yes" if m["trah"] else "no", "yes" if m["lshift"] else "no", "" if x.get("npc_mo") is None else str(x["npc_mo"])]
            else:
                f += [""] * 7
            if x["v"] is not None:
                f += ["%+.4f" % x["v"], x["word"], "%+.4f" % (x["v"] - bare[st])]
            else:
                f += ["", "", ""]
            f += [c["verdict"], c["fifth"], "" if c["dde"] is None else "%+.4f" % c["dde"], "%+.4f" % c["fo"], "%+.4f" % c["score"],
                  c["probe2"], "yes" if c in surv else "no"]
            L.append("\t".join(f))
    rep[:] = [l.rstrip() for l in rep]
    open(a.tsv, "w").write("\n".join(L) + "\n")
    open(a.result, "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main())
