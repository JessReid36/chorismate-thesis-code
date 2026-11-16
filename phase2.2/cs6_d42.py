#!/usr/bin/env python3
"""DIAGNOSTIC 42's reader: is the dianion bound with the charges of D31's RELIABLE conditions?

VDE = E(monoanion) - E(dianion), same geometry, same external charges, positive means bound:
D40's convention, which is Probe 2's, at 27.211386245988 eV per Eh. The dianion half of every
difference is the committed output of the D31 run the monoanion input was copied from. It is
not run again and it is read here by the same rule as the monoanion.

Reading one output. The committed cs6_d40.py is imported by path and its read() is used as it
stands: a run is read only if it terminates normally, converges its SCF, reports the right
number of electrons and says it read as many point charges as its pc.pc holds; the energy is
the last FINAL SINGLE POINT ENERGY; the stability analysis is read by cs6_d33b.py stability()
with ORCA's own verdict sentences as the fallback, and is never guessed. That reader was run
on D40's twenty open-shell outputs, twelve of them with a stability analysis.

The words, fixed in section 2 of the pre-registration:
   METASTABLE                         VDE at or below zero, whatever the stability analysis says
   BOUND                              VDE above zero AND the monoanion's last analysis is STABLE
   ABOVE ZERO, MONOANION NOT STABLE   VDE above zero and the last analysis is not STABLE, or
                                      was not found, or could not be read. Not counted as BOUND.

The 13 conditions are those of d42/CONDITIONS.tsv, which is committed with the inputs. They
are checked here against the table R31.2 of D31's committed result, through the committed
cs6_d42_prepare.py, and against the committed charge files.

Nothing is typed in except committed reference values, each named with its source and each
checked against that source before it is used.

usage: cs6_d42.py --res <results-repo> --code <code-repo> --prereg <file> --result <file> --tsv <file>
"""
import argparse, importlib.util, os, re, sys, textwrap

EV = 27.211386245988            # eV per Eh, vde_designs/cs11_design_vde.py, as D40 uses it
H2K = 627.5094740631
S = "phase2.2/cs_stage6"
NEL_MONO, NEL_DI = 117, 118
CTL_E, TOL_EH = -836.26897688, 1e-6     # G42.4: D40's committed monoanion B_bare_R, d40.tsv
BARE = {"R": "+1.2486", "TS": "+0.7681"}   # D40's bare detachment energies, d40.tsv rows B_bare_R and B_bare_TS
D31_BARRIER, TOL_K = 26.8248, 5e-5      # G42.5: the charge-free barrier of d31/ctl_{R,TS}
P2_FILE = "vde_probe/PROBE2_RESULT.txt"
P2_SENT = "Any ddE measured at such a condition is not a barrier effect and must not be reported as one."
PREREG = "STAGE6_DIAGNOSTIC42_20261115.txt"
D31R = "STAGE6_DIAGNOSTIC31_RESULT_20261107.txt"
D31P = "STAGE6_DIAGNOSTIC31_20261106.txt"       # D31.8, its limitations
D31_Q2 = "No barrier reported from this diagnostic is a prediction of a physical barrier"
D31_Q4 = "ARM B IS NOT THE OPTIMISER'S BEST"
D31_Q1 = "Every condition here above about 3 e inherits that doubt"
H9_A = "a RELIABLE arrangement whose charges sum to zero or more keeps it bound"      # as section 1 of the
H9_B = ("a RELIABLE arrangement whose charges sum to zero or more is METASTABLE at "   # pre-registration quotes H9
        "CAM-B3LYP-D3BJ/def2-SVPD at either state.")
GUARDED = ("d31", "d33a", "d40")
CONTROL = ("C_bare_R", "d40/B_bare_R")
STATES = ("R", "TS")
BIG = "n20_q0.25"                        # the one condition above the total at which D30 found def2-SVPD converged
NOTSTABLE = "ABOVE ZERO, MONOANION NOT STABLE"
NCOND = 13


def stop(m):
    sys.stderr.write("STOP  %s\n" % m)
    raise SystemExit(1)


def load(path, name):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def word_of(v, sb):
    """The word of section 2. The one asymmetry can only take a BOUND away."""
    if v is None:
        return None
    if v <= 0:
        return "METASTABLE"
    if sb.get("ran") and sb.get("last") == "STABLE":
        return "BOUND"
    return NOTSTABLE


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    a.res, a.code, a.prereg = (os.path.abspath(x) for x in (a.res, a.code, a.prereg))
    S6 = os.path.join(a.res, S)
    D = os.path.join(S6, "d42")
    p40 = os.path.join(a.code, "phase2.2", "cs6_d40.py")
    p33 = os.path.join(a.code, "phase2.2", "cs6_d33b.py")
    ppr = os.path.join(a.code, "phase2.2", "cs6_d42_prepare.py")
    for p in (p40, p33, ppr):
        if not os.path.isfile(p):
            stop("%s is missing; this reader imports it by path." % p)
    d40 = load(p40, "d42_cs6_d40")
    d33b = load(p33, "d42_cs6_d33b")
    prep = load(ppr, "d42_cs6_d42_prepare")
    sha, git = d40.sha, d40.git

    rep, gates = [], []
    say = rep.append
    say("STAGE 6, DIAGNOSTIC 42 - RESULT")
    say("Is the dianion bound with the charges of Diagnostic 31's RELIABLE conditions? The first")
    say("test of the hypothesis that binding is lost because the net charge is negative.")
    say("")
    say("Read by phase2.2/cs6_d42.py                 sha256 %s" % sha(__file__))
    say("   every output by cs6_d40.py read           sha256 %s" % sha(p40))
    say("   stability by cs6_d33b.py stability        sha256 %s" % sha(p33))
    say("   R31.2 parsed by cs6_d42_prepare.py        sha256 %s" % sha(ppr))
    say("   pre-registration %-23s sha256 %s" % (os.path.basename(a.prereg), sha(a.prereg)))
    say("")
    say("VDE = E(monoanion) - E(dianion) at the same geometry with the same external charges,")
    say("in eV at %.12f eV per Eh. Positive means bound. Frame 42436," % EV)
    say("CAM-B3LYP-D3BJ/def2-SVPD. The dianion energies are the committed outputs of the D31 runs")
    say("the monoanion inputs were copied from. They are not run again.")
    say("")

    # ------------------------------------------------------------ the conditions, G42.3
    cond = []
    pc_ = os.path.join(D, "CONDITIONS.tsv")
    if not os.path.isfile(pc_):
        stop("d42/CONDITIONS.tsv is missing.")
    ls = open(pc_).read().splitlines()
    hdr = ls[0].split("\t") if ls else []
    for l in ls[1:]:
        if l.strip():
            cond.append(dict(zip(hdr, l.split("\t"))))
    try:
        for c in cond:
            int(c["n"]), float(c["q_e"]), float(c["sum_q_e"]), c["condition"], c["class"]
    except (KeyError, ValueError):
        stop("d42/CONDITIONS.tsv does not hold the columns condition, n, q_e, sum_q_e and class as numbers and names (G42.3). Nothing is written.")
    p31 = os.path.join(S6, D31R)
    try:
        rel = prep.reliable_in_r312(open(p31).read()) if os.path.isfile(p31) else []
    except (SystemExit, Exception):
        rel = []
    mine = [(int(c["n"]), c["q_e"]) for c in cond]
    if len(cond) != NCOND or sorted(mine) != sorted(rel) or len(set(mine)) != len(mine):
        gates.append("G42.3 the %d conditions of CONDITIONS.tsv are not exactly the %d cells that R31.2 of %s calls RELIABLE"
                     % (len(cond), len(rel), D31R))
    for c in cond:
        name = c["condition"]
        if name != prep.cname(int(c["n"]), c["q_e"]):
            gates.append("G42.3 %s is not the name of n = %s, q = %s" % (name, c["n"], c["q_e"]))
        qs = {}
        for st in STATES:
            p = os.path.join(D, "%s_%s" % (name, st), "pc.pc")
            if not os.path.isfile(p):
                gates.append("G42.3 %s_%s/pc.pc is missing" % (name, st))
                continue
            try:
                qs[st] = prep.charges(p)
            except (SystemExit, Exception):
                gates.append("G42.3 %s_%s/pc.pc could not be read" % (name, st))
        if len(qs) == 2:
            if qs["R"] != qs["TS"]:
                gates.append("G42.3 %s does not carry the same charges at the two states" % name)
            s = sum(qs["R"])
            if (len(qs["R"]) != int(c["n"]) or abs(s - float(c["sum_q_e"])) > 5e-7
                    or c["class"] != ("H" if s > -5e-7 else "N")
                    or any(abs(abs(x) - float(c["q_e"])) > 1e-9 for x in qs["R"])):
                gates.append("G42.3 %s: CONDITIONS.tsv does not describe the charges its pc.pc holds" % name)
        c["sum"] = float(c["sum_q_e"])
        c["total"] = int(c["n"]) * float(c["q_e"])
    H = [c for c in cond if c["class"] == "H"]
    N = [c for c in cond if c["class"] == "N"]

    # ------------------------------------------------------------ every run and its source
    runs = [(CONTROL[0], CONTROL[1], "plain", "R")]
    for c in H + N:
        runs += [("%s_%s" % (c["condition"], st), "d31/B/%s_%s" % (c["condition"], st), "pc", st) for st in STATES]
    if len(runs) != 1 + 2 * NCOND and not gates:
        stop("the reader's own list holds %d runs, not %d" % (len(runs), 1 + 2 * NCOND))

    # ------------------------------------------------------------ G42.1 and G42.2
    nfile = 0
    for run, src, kind, st in runs:
        for f in (("substrate.xyz", "pc.pc") if kind == "pc" else ("substrate.xyz",)):
            x, y = os.path.join(D, run, f), os.path.join(S6, src, f)
            if not (os.path.isfile(x) and os.path.isfile(y)) or sha(x) != sha(y):
                gates.append("G42.1 %s/%s is not byte-identical to %s/%s" % (run, f, src, f))
            nfile += 1
        x, y = os.path.join(D, run, "substrate.xyz"), os.path.join(S6, "d33a", "sp_42436_%s" % st, "substrate.xyz")
        if not (os.path.isfile(x) and os.path.isfile(y)) or sha(x) != sha(y):
            gates.append("G42.1 %s/substrate.xyz is not d33a/sp_42436_%s/substrate.xyz" % (run, st))
        pi, ps = os.path.join(D, run, "job.inp"), os.path.join(S6, src, "job.inp")
        if not (os.path.isfile(pi) and os.path.isfile(ps)):
            gates.append("G42.2 %s or its source has no input" % run)
            continue
        if run == CONTROL[0]:
            if sha(pi) != sha(ps):
                gates.append("G42.2 the control input is not byte-identical to %s/job.inp" % src)
            continue
        li, lsrc = open(pi).read().split("\n"), open(ps).read().split("\n")
        diff = [i for i in range(min(len(li), len(lsrc))) if li[i] != lsrc[i]]
        if len(li) != len(lsrc) or len(diff) != 1 \
                or lsrc[diff[0]] != "* xyzfile -2 1 substrate.xyz" or li[diff[0]] != "* xyzfile -1 2 substrate.xyz":
            gates.append("G42.2 %s does not differ from %s/job.inp in the charge-and-multiplicity line alone" % (run, src))

    # ------------------------------------------------------------ G42.6
    rc, _ = git(a.res, "diff", "--quiet", "HEAD", "--", *["%s/%s" % (S, x) for x in GUARDED])
    if rc is None:
        gates.append("G42.6 could not be checked: git did not run")
    elif rc != 0:
        gates.append("G42.6 a committed file under %s differs from HEAD" % ", ".join(GUARDED))
    committed = sorted(set([s for _, s, _, _ in runs[1:]] + ["d31/ctl_R", "d31/ctl_TS", "d40/B_bare_R", "d40/B_bare_TS"]))
    rc, pcm = git(a.res, "log", "--diff-filter=A", "--format=%H", "-1", "--", "%s/%s" % (S, PREREG))
    if rc != 0 or not pcm:
        gates.append("G42.6 the commit that added %s could not be found" % PREREG)
        pcm = ""
    else:
        rc1, b1 = git(a.res, "rev-parse", "%s:%s/%s" % (pcm, S, PREREG))
        rc2, b2 = git(a.res, "hash-object", a.prereg)
        if rc1 != 0 or rc2 != 0 or b1 != b2:
            gates.append("G42.6 the pre-registration given is not the committed %s" % PREREG)
        for s in committed:
            rc1, b1 = git(a.res, "rev-parse", "%s:%s/%s/job.out" % (pcm, S, s))
            rc2, b2 = git(a.res, "hash-object", os.path.join(S6, s, "job.out"))
            if rc1 != 0 or rc2 != 0 or b1 != b2:
                gates.append("G42.6 %s/job.out is not the file at the commit the pre-registration is committed on" % s)

    # ------------------------------------------------------------ the committed references
    ref = {}
    p40t = os.path.join(S6, "d40.tsv")
    if os.path.isfile(p40t):
        l40 = open(p40t).read().splitlines()
        h40 = l40[0].split("\t")
        for l in l40[1:]:
            r = dict(zip(h40, l.split("\t")))
            ref[r.get("run")] = r
    if ref.get("B_bare_R", {}).get("E_mono_Eh") != "%.8f" % CTL_E:
        gates.append("G42.4 the committed source of %.8f Eh, d40.tsv, does not print it for B_bare_R" % CTL_E)
    for st in STATES:
        if ref.get("B_bare_%s" % st, {}).get("VDE_eV") != BARE[st]:
            gates.append("R42.5 the committed source of the bare substrate's %s eV at %s, d40.tsv, does not print it" % (BARE[st], st))
    pp = os.path.join(S6, P2_FILE)
    if not os.path.isfile(pp) or P2_SENT not in " ".join(open(pp).read().split()):
        gates.append("R42.6 Probe 2's sentence is not in %s" % P2_FILE)
    pq = os.path.join(S6, D31P)
    tq = " ".join(open(pq).read().split()) if os.path.isfile(pq) else ""
    for q in (D31_Q1, D31_Q2, D31_Q4):
        if q not in tq:
            gates.append("the sentence quoted from D31's limitations, \"%s\", is not in %s" % (q, D31P))
    tp = " ".join(open(a.prereg).read().split())
    for q in (H9_A, H9_B):
        if q not in tp:
            gates.append("the words of H9, \"%s\", are not in the pre-registration given" % q)

    # ------------------------------------------------------------ read everything
    DI = {}
    for s in committed:
        mono = s.startswith("d40/")
        r = d40.read(os.path.join(S6, s), "pc" if s.startswith("d31/B/") else "plain", NEL_MONO if mono else NEL_DI, d33b)
        if r["E"] is None:
            gates.append("the committed output %s is not readable by the rule of G42.7 at %d electrons: %s"
                         % (s, NEL_MONO if mono else NEL_DI, "; ".join(r["why"])))
        DI[s] = r
    MO = {}
    for run, src, kind, st in runs:
        MO[run] = d40.read(os.path.join(D, run), kind, NEL_MONO, d33b)

    # ------------------------------------------------------------ G42.4 and G42.5
    c0 = MO[CONTROL[0]]
    if c0["E"] is None:
        gates.append("G42.4 the control %s is not readable: %s" % (CONTROL[0], "; ".join(c0["why"])))
    elif abs(c0["E"] - CTL_E) > TOL_EH:
        gates.append("G42.4 the control gives %.8f Eh, which differs from D40's %.8f by %.2e, over the %.0e limit"
                     % (c0["E"], CTL_E, abs(c0["E"] - CTL_E), TOL_EH))
    bare, barrier = {}, None
    if all(DI[k]["E"] is not None for k in ("d31/ctl_R", "d31/ctl_TS", "d40/B_bare_R", "d40/B_bare_TS")):
        for st in STATES:
            bare[st] = (DI["d40/B_bare_%s" % st]["E"] - DI["d31/ctl_%s" % st]["E"]) * EV
            if "%+.4f" % bare[st] != BARE[st]:
                gates.append("R42.5 the bare substrate's detachment energy at %s, recomputed from the committed outputs, is %+.4f eV and not D40's %s"
                             % (st, bare[st], BARE[st]))
        barrier = (DI["d31/ctl_TS"]["E"] - DI["d31/ctl_R"]["E"]) * H2K
        if abs(barrier - D31_BARRIER) > TOL_K:
            gates.append("G42.5 the charge-free barrier from d31/ctl_R and d31/ctl_TS is %.4f kcal/mol, not %.4f within %.0e"
                         % (barrier, D31_BARRIER, TOL_K))

    if gates:
        say("GATES")
        for g in gates:
            say("   FAILED  %s" % g)
        say("")
        say("A FAILURE OF G42.1 TO G42.6 STOPS THE DIAGNOSTIC. Nothing else is readable.")
        say("The reader stops in the same way where a committed file it must quote or read is not")
        say("as the pre-registration names it.")
        print("\n".join(rep))
        stop("%d gate(s) failed. Nothing is written." % len(gates))

    say("=" * 78)
    say("G42.1 TO G42.6  THE STOPPING GATES")
    say("=" * 78)
    say("")
    say("   G42.1  %d geometry and point-charge files are byte-identical to the committed runs'" % nfile)
    say("          own, and every geometry to d33a/sp_42436_{R,TS}.")
    say("   G42.2  each of the %d monoanion inputs differs from its committed D31 source in the" % (2 * NCOND))
    say("          charge-and-multiplicity line alone; the control input is d40/B_bare_R/job.inp.")
    say("   G42.3  the %d conditions are exactly the cells R31.2 of D31's result calls RELIABLE," % NCOND)
    say("          and CONDITIONS.tsv describes the charges the committed pc.pc files hold.")
    say("   G42.6  no committed file under %s or %s differs from HEAD, and the %d" % (", ".join(GUARDED[:-1]), GUARDED[-1], len(committed)))
    say("          committed outputs read are the files at the commit the pre-registration was")
    say("          committed on, %s." % pcm[:8])
    say("   G42.4 and G42.5 are the two controls, under R42.1.")
    say("")

    # ------------------------------------------------------------ R42.1
    say("=" * 78)
    say("R42.1  THE CONTROLS")
    say("=" * 78)
    say("")
    say("   G42.4  C_bare_R, D40's committed monoanion input for the bare substrate, run again:")
    say("          %.8f Eh   D40 %.8f   differ %.2e   limit %.0e   PASS" % (c0["E"], CTL_E, abs(c0["E"] - CTL_E), TOL_EH))
    say("          against the committed dianion d31/ctl_R: VDE %+.4f eV   D40 %s"
        % ((c0["E"] - DI["d31/ctl_R"]["E"]) * EV, BARE["R"]))
    say("   G42.5  the charge-free barrier from d31/ctl_R and d31/ctl_TS:")
    say("          %.4f kcal/mol   D31 %.4f   differ %.2e   limit %.0e   PASS" % (barrier, D31_BARRIER, abs(barrier - D31_BARRIER), TOL_K))
    say("")

    # ------------------------------------------------------------ R42.2
    def fs(v, fmt="%.6f"):
        return "-" if v is None else fmt % v

    say("=" * 78)
    say("R42.2  THE 13 CONDITIONS AT BOTH STATES")
    say("=" * 78)
    say("")
    say("   BOUND: VDE above zero and the monoanion's last stability analysis STABLE.")
    say("   METASTABLE: VDE at or below zero. Stability is the monoanion's, first analysis then")
    say("   last. read: the number of point charges ORCA says it read.")
    say("")
    say("   condition  st   n  q     sum e    cl    VDE eV  <S**2>    stability          TRAH  shift read  word")
    B, tail = {}, []
    for c in H + N:
        name = c["condition"]
        for st in STATES:
            run, src = "%s_%s" % (name, st), "d31/B/%s_%s" % (name, st)
            m, d = MO[run], DI[src]
            v = None if m["E"] is None else (m["E"] - d["E"]) * EV
            sb = m.get("stab") or {}
            w = word_of(v, sb)
            npc = None
            if m.get("exists"):
                mm = re.search(r"Reading point charge file\s*\.*\s*ok \((\d+) point charges?\)", m["text"])
                npc = int(mm.group(1)) if mm else None
            B[(name, st)] = dict(v=v, word=w, m=m, d=d, src=src, run=run, npc=npc)
            if v is None:
                say("   %-10s %-2s  %2s  %s  %+.4f  %s   NOT MEASURED (G42.7)   E(dianion) %.8f"
                    % (name, st, c["n"], c["q_e"], c["sum"], c["class"], d["E"]))
                for x in m["why"]:
                    say("        %s" % x)
                continue
            cell = "NOT FOUND" if not sb["ran"] else ("NOT READ" if sb["first"] is None else "%s/%s" % (sb["first"], sb["last"]))
            say("   %-10s %-2s  %2s  %s  %+.4f  %s   %+8.4f  %s  %-17s  %-4s  %-5s %-4s  %s"
                % (name, st, c["n"], c["q_e"], c["sum"], c["class"], v, fs(m["s2"]), cell,
                   "yes" if m["trah"] else "no", "yes" if m["lshift"] else "no", "-" if npc is None else str(npc), w))
            if not sb["ran"]:
                tail.append("%s %s: no stability analysis was found in the output, though its input asks for one. Lines of the output that mention stability: %s" % (
                    name, st, " | ".join(sb["quote"]) or "none"))
            elif sb["first"] is None:
                tail.append("%s %s: %s. ORCA's own lines: %s" % (name, st, sb["how"], " | ".join(sb["quote"]) or "(none found)"))
            else:
                if sb["how"] != "cs6_d33b.stability":
                    tail.append("%s %s: read from %s." % (name, st, sb["how"]))
                if sb["first"] == "UNSTABLE":
                    tail.append("%s %s: the first solution was unstable%s; %s." % (
                        name, st, "" if sb["nneg"] is None else " (%d negative, lowest %+.6f Eh)" % (sb["nneg"], sb["lo"]),
                        "1 analysis was made" if sb["n"] == 1 else "%d analyses were made" % sb["n"]))
            if w == NOTSTABLE:
                tail.append("%s %s: the value is above zero and the monoanion's last stability analysis is not STABLE, so by section 2 of the pre-registration it is not counted as BOUND." % (name, st))
    say("")
    say("   cl: H, the charges sum to zero or more; N, they sum to less than zero.")
    say("   Both energies of every row are in d42.tsv.")
    for l in tail:
        for i, x in enumerate(textwrap.wrap(l, 80)):
            say("   %s%s" % ("NOTE  " if i == 0 else "      ", x))
    say("")
    say("   D42.4: where a detachment energy is negative, its sign is the reading and its size")
    say("   is not. <S**2> is reported and not corrected; the ideal value is 0.75.")
    say("")

    def keeps(c):
        """G42.8 first: a condition with a state NOT MEASURED is not classified. Then R42.3:
        yes, BOUND at both states; no, METASTABLE at either; open, neither."""
        w = [B[(c["condition"], st)]["word"] for st in STATES]
        if None in w:
            return "unclassified"
        if "METASTABLE" in w:
            return "no"
        if w == ["BOUND", "BOUND"]:
            return "yes"
        return "open"

    KEEPS = {"yes": "yes", "no": "NO", "open": "not established", "unclassified": "NOT CLASSIFIED (G42.8)"}

    def cellw(c, st):
        w = B[(c["condition"], st)]["word"]
        return "NOT MEASURED" if w is None else ("ABOVE ZERO *" if w == NOTSTABLE else w)

    # ------------------------------------------------------------ R42.3
    say("=" * 78)
    say("R42.3  THE TEST OF H9, ON THE 8 CONDITIONS WHOSE CHARGES SUM TO ZERO OR MORE")
    say("=" * 78)
    say("")
    say("   H9, as section 1 of the pre-registration quotes it from section 11 of the")
    say("   consolidated write-up:")
    for l in textwrap.wrap('"%s", to be given up if "%s"' % (H9_A, H9_B), 72):
        say("      %s" % l)
    say("   \"It\" is the dianion. A condition KEEPS THE DIANION BOUND if it is BOUND at both")
    say("   states.")
    say("")
    say("   condition    sum e     R             TS            keeps the dianion bound")
    for c in H:
        say("   %-10s  %+.4f   %-12s  %-12s  %s" % (c["condition"], c["sum"], cellw(c, "R"), cellw(c, "TS"), KEEPS[keeps(c)]))
    if any(B[(c["condition"], st)]["word"] == NOTSTABLE for c in H for st in STATES):
        say("   * %s (section 2 of the pre-registration)." % NOTSTABLE)
    say("")
    nyes = sum(1 for c in H if keeps(c) == "yes")
    small = [c for c in H if c["condition"] != BIG]
    say("   Of the %d conditions of this class, %d keep the dianion bound at both states." % (len(H), nyes))
    big = [c for c in H if c["condition"] == BIG]
    if big:
        say("   %s carries %.2f e in all, and D31's first limitation says of its conditions:" % (BIG, big[0]["total"]))
        say("   \"%s\"," % D31_Q1)
        say("   which is D30's doubt about def2-SVPD. Without %s: %d of %d." % (BIG, sum(1 for c in small if keeps(c) == "yes"), len(small)))
    say("")
    fals = [(c["condition"], st) for c in H for st in STATES if B[(c["condition"], st)]["word"] == "METASTABLE"]
    open_ = [(c["condition"], st, "NOT MEASURED" if B[(c["condition"], st)]["word"] is None else NOTSTABLE)
             for c in H for st in STATES if B[(c["condition"], st)]["word"] in (None, NOTSTABLE)]
    if fals:
        say("   H9 IS FALSIFIED.")
        for l in textwrap.wrap("METASTABLE in this class: " + "; ".join("%s at %s, %+.4f eV" % (n, st, B[(n, st)]["v"]) for n, st in fals) + ".", 76):
            say("   %s" % l)
        if set(n for n, _ in fals) == {BIG}:
            say("   The only condition that falsifies it is %s, which carries the limitation above." % BIG)
        if open_:
            for l in textwrap.wrap("Not established in this class, and not needed for the falsification: "
                                   + "; ".join("%s at %s, %s" % x for x in open_) + ".", 76):
                say("   %s" % l)
    elif open_:
        say("   THE TEST IS INCOMPLETE.")
        say("   No condition of this class is METASTABLE at either state, and these are not")
        say("   established:")
        for x in open_:
            say("      %s at %s: %s" % x)
    else:
        say("   H9 IS NOT FALSIFIED BY THESE %d CONDITIONS." % len(H))
        say("   D42.3: not falsified is not established. The class is %d conditions on one site" % len(H))
        say("   pool, and the sign of the sum is confounded with n and with site throughout.")
    say("")

    # ------------------------------------------------------------ R42.4
    say("=" * 78)
    say("R42.4  THE 5 CONDITIONS WHOSE CHARGES SUM TO LESS THAN ZERO")
    say("=" * 78)
    say("")
    say("   condition    sum e     R             TS")
    for c in N:
        say("   %-10s  %+.4f   %-12s  %-12s" % (c["condition"], c["sum"], cellw(c, "R"), cellw(c, "TS")))
    if any(B[(c["condition"], st)]["word"] == NOTSTABLE for c in N for st in STATES):
        say("   * %s (section 2 of the pre-registration)." % NOTSTABLE)
    say("")
    say("   Reported. H9 neither predicts nor forbids either word for these.")
    say("")

    # ------------------------------------------------------------ R42.5
    say("=" * 78)
    say("R42.5  CHANGE IN DETACHMENT ENERGY FROM THE BARE SUBSTRATE AT THE SAME STATE")
    say("=" * 78)
    say("")
    say("   The bare substrate: %s eV at the reactant and %s eV at the transition state (D40)." % (BARE["R"], BARE["TS"]))
    say("")
    say("   condition    sum e      n        R eV       TS eV")
    for c in sorted(H + N, key=lambda c: (-c["sum"], int(c["n"]), c["q_e"])):
        cells = []
        for st in STATES:
            v = B[(c["condition"], st)]["v"]
            cells.append("%+10.4f" % (v - bare[st]) if v is not None else "%10s" % "-")
        say("   %-10s  %+.4f   %2s   %s  %s" % (c["condition"], c["sum"], c["n"], cells[0], cells[1]))
    say("")
    say("   In order of the sum of the charges. Reported, not interpreted. Nothing is fitted.")
    say("   Each difference is formed from the energies and then rounded, so it can differ by")
    say("   0.0001 from the difference of two rounded values. A dash is a state NOT MEASURED.")
    say("")

    # ------------------------------------------------------------ R42.6
    say("=" * 78)
    say("R42.6  THE BARRIER CHANGE BESIDE THE TWO TESTS")
    say("=" * 78)
    say("")
    say("   ddE = (E(TS) - E(TS, control)) - (E(R) - E(R, control)), from D31's committed dianion")
    say("   outputs with d31/ctl_R and d31/ctl_TS as the control, at %.10f kcal/mol per" % H2K)
    say("   Eh. D31's result reports verdicts and no ddE; these are computed here.")
    say("")
    say("   Probe 2, %s:" % P2_FILE)
    for l in textwrap.wrap('"%s"' % P2_SENT, 72):
        say("      %s" % l)
    say("   A lowering is a difference between a transition state and a reactant, so the")
    say("   sentence applies to a condition's ddE if either of its states is METASTABLE (the")
    say("   rule of D40's R40.7).")
    say("")
    say("   condition  cl    sum e    ddE kcal/mol   R             TS            Probe 2's sentence")
    both = []
    for c in H + N:
        name = c["condition"]
        e = {st: DI["d31/B/%s_%s" % (name, st)]["E"] for st in STATES}
        dde = ((e["TS"] - DI["d31/ctl_TS"]["E"]) - (e["R"] - DI["d31/ctl_R"]["E"])) * H2K
        c["dde"] = dde
        w = [B[(name, st)]["word"] for st in STATES]
        if None in w:
            verdict = "NOT CLASSIFIED (G42.8): %s NOT MEASURED" % " and ".join(st for st in STATES if B[(name, st)]["word"] is None)
        elif "METASTABLE" in w:
            verdict = "APPLIES"
        elif w == ["BOUND", "BOUND"]:
            verdict = "does not reach"
            both.append(c)
        else:
            verdict = "neither is said: %s above zero, monoanion not stable" % " and ".join(
                st for st in STATES if B[(name, st)]["word"] == NOTSTABLE)
        say("   %-10s %s   %+.4f   %+11.4f    %-12s  %-12s  %s" % (name, c["class"], c["sum"], dde, cellw(c, "R"), cellw(c, "TS"), verdict))
    say("")
    say("   APPLIES: a state is METASTABLE, so Probe 2's sentence applies to this ddE.")
    say("   does not reach: both states are BOUND, so this objection does not reach this ddE at")
    say("   this level and basis.")
    say("")
    if both:
        say("   RELIABLE ON FOUR LEGS IN D31 AND BOUND AT BOTH STATES HERE:")
        for c in both:
            say("      %-10s  ddE %+9.4f kcal/mol   %s charge%s, sum %+.4f e%s"
                % (c["condition"], c["dde"], c["n"], "" if c["n"] == "1" else "s", c["sum"],
                   "   (%.2f e in all: D31's first limitation)" % c["total"] if c["condition"] == BIG else ""))
        say("")
        say("   D42.8: no ddE listed here is made quotable by D42. D31's own limitations stand:")
        say("   \"%s\"," % D31_Q2)
        say("   and \"%s\". D42.7: every one of these conditions sits" % D31_Q4)
        say("   inside r_safe.")
        say("   Whether any of these numbers may be quoted, and as what, is a decision for Jessica")
        say("   and a dated document of its own.")
    else:
        say("   NO CONDITION IS RELIABLE ON FOUR LEGS IN D31 AND BOUND AT BOTH STATES HERE.")
    say("")

    # ------------------------------------------------------------ limits
    say("=" * 78)
    say("WHAT THIS CANNOT SAY (section 6 of the pre-registration, in short; its own words govern)")
    say("=" * 78)
    say("")
    for l in ("D42.1  It withdraws nothing and changes no verdict.",
              "D42.2  It does not make binding a leg of the criterion, and it does not take the",
              "       decision on a neutrality constraint.",
              "D42.3  Not falsified is not established: 8 conditions on one site pool, with the",
              "       sign of the sum confounded with n and with site.",
              "D42.4  Unrestricted Kohn-Sham differences; spin contamination is not corrected; where",
              "       a detachment energy is negative its sign is the reading, not its size.",
              "D42.5  One basis, def2-SVPD; n20_q0.25 lies above the total at which D30 found that",
              "       basis converged.",
              "D42.6  Vertical, not adiabatic. One frame. In vacuo.",
              "D42.7  Every one of the 13 conditions sits inside r_safe.",
              "D42.8  A ddE listed under R42.6 is not made quotable by D42.",
              "D42.9  It adopts no level of theory and says nothing about D39.",
              "D42.10 Nothing here is shaped: a result that falsifies H9 is reported as plainly as",
              "       one that does not, by rules fixed before any output existed."):
        say("   %s" % l)
    say("")
    nm = sum(1 for b in B.values() if b["v"] is not None)
    say("%d of the %d monoanion runs are measured.%s" % (
        nm, 2 * NCOND, "" if nm == 2 * NCOND else " %d NOT MEASURED, each listed above with ORCA's own words." % (2 * NCOND - nm)))

    # ------------------------------------------------------------ the table
    L = ["condition\tstate\trun\tsource\tn\tq_e\tsum_q_e\tclass\tstatus\tE_mono_Eh\tE_di_Eh\tVDE_eV\tword\tS2\tstab_first\tstab_last\tTRAH\tlevel_shift\tcharges_read\tdVDE_vs_bare_eV\tddE_kcal"]
    for c in H + N:
        for st in STATES:
            b = B[(c["condition"], st)]
            m, d, v = b["m"], b["d"], b["v"]
            head = "%s\t%s\t%s\t%s\t%s\t%s\t%+.4f\t%s" % (c["condition"], st, b["run"], b["src"], c["n"], c["q_e"], c["sum"], c["class"])
            if v is None:
                L.append("%s\tNOT MEASURED\t\t%.8f\t\t\t\t\t\t\t\t\t\t%+.4f" % (head, d["E"], c["dde"]))
                continue
            sb = m.get("stab") or {}
            L.append("%s\tmeasured\t%.8f\t%.8f\t%+.4f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%+.4f\t%+.4f" % (
                head, m["E"], d["E"], v, b["word"], "" if m["s2"] is None else "%.6f" % m["s2"],
                sb.get("first") or ("NOT FOUND" if not sb.get("ran") else "NOT READ"),
                sb.get("last") or ("NOT FOUND" if not sb.get("ran") else "NOT READ"),
                "yes" if m["trah"] else "no", "yes" if m["lshift"] else "no", "" if b["npc"] is None else str(b["npc"]),
                v - bare[st], c["dde"]))
    rep[:] = [l.rstrip() for l in rep]
    open(a.tsv, "w").write("\n".join(L) + "\n")
    open(a.result, "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main())
