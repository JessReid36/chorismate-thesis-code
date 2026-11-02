#!/usr/bin/env python3
"""
cs6_d22.py - STAGE6_DIAGNOSTIC22_20261101.txt.

Do the walled and wall-at-zero legs of d_shield sit on stable solutions?

  python3 cs6_d22.py prepare <results_repo> <OUT> <HPC_DIR>
  python3 cs6_d22.py analyse <results_repo> <RUNS> <OUT_DIR>

NOTHING ALREADY COMMITTED IS RETYPED OR RE-RUN.
  cs6_transfer    head, H
  cs6_additivity  PBS
  cs13_qscale     FRAME_SHA, GEOM
  cs6_d19         read_stab                   (as corrected by D19 NOTE 1)
  cs6_d21         gate1, stab_head, fsp_last, gate_targets, CTL, D20_CTL_E0, D20_RESULT
                  and D21's bare verdicts, READ from its committed result, never re-run.
                  READING 4 takes the bare leg's energies from D21.gate_targets, so they are
                  the same numbers D21 itself gated against its own pre-registration.

ERRATUM 1 (2 November 2026) is implemented here: two of D22.2's nine conditions carry no
readable Dset, so READING 2 labels them separately and the ALL ARMS VERIFIED text names the
seven that do. No gate, tolerance, arithmetic or verdict rule is changed by it.

Every job.inp is COPIED VERBATIM from the committed run it measures, with exactly one
substitution - the production head swapped for the stability head - and gated to reduce back to
that run's committed input BYTE FOR BYTE. D21 rebuilt its inputs from parts; a walled input
carries 24 coordinates and up to 20 Ne> blocks of literal ECP parameters inline, so rebuilding
it would insert a transcription risk for no gain. D22.2 records that departure.
"""
import hashlib, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cs13_qscale as Q13
import cs6_transfer as TR
import cs6_additivity as AD
import cs6_d19 as D19
import cs6_d21 as D21

PREREG = "STAGE6_DIAGNOSTIC22_20261101.txt"
D21_RESULT = "d21/DIAGNOSTIC21_RESULT.txt"
BASIS = "def2-SVPD"
GATE_EH = 1e-6          # D22.2 GATE 2
BUDGET = 1.067          # D11.3, for READING 4
FRAME_TOL = 1e-8        # D22.2 GATE 3(c), angstrom
REP_REL, REP_ABS = 1e-5, 5e-8    # D22.3 READING 3

# D22.2 TABLE A: name, n, sum|q|, walled source, wall-at-zero source.
COND = (
    ("S6-2_D17",     2, 1.0777, "d17/def2-SVPD/S6-2_wall_M",       "d16/def2-SVPD/S6-2_wall0"),
    ("S6-2_D16",     2, 1.3000, "d16/def2-SVPD/S6-2_wall_F1",      "d16/def2-SVPD/S6-2_wall0"),
    ("S6-3_D16",     3, 1.3000, "d16/def2-SVPD/S6-3_wall_F1",      "d16/def2-SVPD/S6-3_wall0"),
    ("S6-5_D16",     5, 1.3000, "d16/def2-SVPD/S6-5_wall_F1",      "d16/def2-SVPD/S6-5_wall0"),
    ("S6-5_D17",     5, 1.6166, "d17/def2-SVPD/S6-5_wall_M",       "d16/def2-SVPD/S6-5_wall0"),
    ("S6-10_D16",   10, 1.3000, "d16/def2-SVPD/S6-10_wall_F1",     "d16/def2-SVPD/S6-10_wall0"),
    ("S6-10_D15",   10, 2.1900, "additivity/def2-SVPD/n10_wall_A", "additivity/def2-SVPD/n10_wall0"),
    ("S6-20_D16F1", 20, 1.3000, "d16/def2-SVPD/S6-20_wall_F1",     "d16/def2-SVPD/S6-20_wall0"),
    ("S6-20_D16F2", 20, 2.1900, "d16/def2-SVPD/S6-20_wall_F2",     "d16/def2-SVPD/S6-20_wall0"),
)

# D22.2 GATE 2: the committed FINAL SINGLE POINT ENERGY of every source, as D22.2 states it.
TARGETS = {
    "additivity/def2-SVPD/n10_wall0":   (-836.241243269, -836.213174266),
    "additivity/def2-SVPD/n10_wall_A":  (-836.095139076, -836.083416299),
    "d16/def2-SVPD/S6-10_wall0":        (-836.241243269, -836.213174266),
    "d16/def2-SVPD/S6-10_wall_F1":      (-836.165820348, -836.147139467),
    "d16/def2-SVPD/S6-20_wall0":        (-836.238159169, -836.209931251),
    "d16/def2-SVPD/S6-20_wall_F1":      (-836.207698852, -836.187162127),
    "d16/def2-SVPD/S6-20_wall_F2":      (-836.171935394, -836.157046998),
    "d16/def2-SVPD/S6-2_wall0":         (-836.247877118, -836.218666964),
    "d16/def2-SVPD/S6-2_wall_F1":       (-836.281479626, -836.266268010),
    "d16/def2-SVPD/S6-3_wall0":         (-836.247205444, -836.218280315),
    "d16/def2-SVPD/S6-3_wall_F1":       (-836.213409241, -836.197765944),
    "d16/def2-SVPD/S6-5_wall0":         (-836.242668794, -836.214156860),
    "d16/def2-SVPD/S6-5_wall_F1":       (-836.254917840, -836.238873106),
    "d17/def2-SVPD/S6-2_wall_M":        (-836.271028893, -836.253277715),
    "d17/def2-SVPD/S6-5_wall_M":        (-836.258061280, -836.245325259),
}

# The two sources D22.2 names as the replicate pair: same position set, byte-identical input.
REPLICATE = ("additivity/def2-SVPD/n10_wall0", "d16/def2-SVPD/S6-10_wall0")

# STAGE6_DIAGNOSTIC22_ERRATUM1_20261102.txt E22.1: two of D22.2's nine conditions carry NO
# readable Dset. Their own diagnostics report them NOT READABLE because the bare arm at R is
# STRAINED under D6.3 - an electronic-health failure, which is a DIFFERENT test from the one
# this diagnostic applies. Stable arms do not make their d_shield quotable (E22.2). The verdict
# RULE of D22.4 is unchanged; only the label and the printed text below change (E22.3).
NO_DSET = {
    "S6-2_D16": "DIAGNOSTIC 16 READING 1: NOT READABLE, bare R STRAINED, gap 0.567 eV",
    "S6-5_D17": "DIAGNOSTIC 17 READING 2: NOT READABLE, bare R STRAINED, gap 0.880 eV",
}


def stop(m):
    print("STOP: " + m); sys.stdout.flush(); sys.exit(1)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def tag(src):
    """A source's run directory name. Tree plus basename, which is unique across the three trees."""
    p = src.split("/")
    if len(p) != 3:
        stop("%r is not a <tree>/<basis>/<run> source path" % src)
    if p[1] != BASIS:
        stop("%r is not at %s" % (src, BASIS))
    return "%s_%s" % (p[0], p[2])


def arms():
    """Every (source, arm, [conditions served]) D22.2 requires, walled first then wall-at-zero."""
    walled, zero = [], {}
    for nm, n, f, w, w0 in COND:
        walled.append((w, "walled", [nm]))
        zero.setdefault(w0, []).append(nm)
    out = list(walled) + [(k, "wall0", v) for k, v in sorted(zero.items())]
    names = [tag(s) for s, _, _ in out]
    if len(set(names)) != len(names):
        stop("two sources collide on the same run directory name: %s" % names)
    if len(walled) != len(COND):
        stop("%d walled sources were assembled, expected the %d D22.2 conditions"
             % (len(walled), len(COND)))
    if len(out) * len(Q13.GEOM) != 30:
        stop("%d sources x %d states is %d runs, not the 30 D22.2 requires"
             % (len(out), len(Q13.GEOM), len(out) * len(Q13.GEOM)))
    return out


def frames(T6):
    """D22.2 GATE 3(c): the committed frame, verified by SHA before it is used as a reference."""
    geo = {}
    for st, fn in (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz")):
        p = T6 / "frame_42436" / fn
        if not p.is_file():
            stop("the committed frame is missing: %s" % p)
        if sha(p) != Q13.FRAME_SHA[fn]:
            stop("%s is not the committed frame" % p)
        L = p.read_text().splitlines()
        n = int(L[0].split()[0])
        geo[st] = [(x.split()[0], float(x.split()[1]), float(x.split()[2]), float(x.split()[3]))
                   for x in L[2:2 + n]]
    return geo


def gate_targets(T6):
    """D22.2 GATE 2: every target READ from the committed file and checked against D22.2."""
    out = {}
    for src, (eR, eTS) in sorted(TARGETS.items()):
        for st, want in (("R", eR), ("TS", eTS)):
            p = T6 / ("%s_%s" % (src, st)) / "job.out"
            if not p.is_file():
                stop("the committed source arm is missing: %s" % p)
            e = D21.fsp_last(p.read_text(errors="replace"))
            if abs(e - want) > 5e-9:
                stop("the committed energy for %s %s is %.9f, not the %.9f D22.2 states"
                     % (src, st, e, want))
            out[(src, st)] = e
    if len(out) != 2 * len(TARGETS):
        stop("%d targets were read, expected %d" % (len(out), 2 * len(TARGETS)))
    return out


def gate_selfcontained(src, st, text):
    """D22.2 GATE 3(b). job.inp alone must be the whole input, or copying it copies too little."""
    for bad in ("xyzfile", "%pointcharges"):
        if bad in text:
            stop("%s_%s references an external file (%r), so job.inp alone is not the complete\n"
                 "   input and copying it verbatim would run this against whatever file sat\n"
                 "   beside it. D22.2 GATE 3(b)." % (src, st, bad))


def gate_geometry(src, st, text, geo):
    """D22.2 GATE 3(c). d_shield compares bare and walled AT THE SAME GEOMETRY."""
    m = re.search(r"^\* xyz\s+-2\s+1\s*$", text, re.M)
    if not m:
        stop("%s_%s carries no '* xyz -2 1' coordinate block; D22.2 GATE 3(c) cannot check its\n"
             "   geometry against the committed frame." % (src, st))
    rows = []
    for line in text[m.end():].splitlines()[1:]:
        p = line.split()
        if not p or p[0].startswith("Ne>"):
            break
        rows.append((p[0], float(p[1]), float(p[2]), float(p[3])))
    ref = geo[st]
    if len(rows) != len(ref):
        stop("%s_%s has %d substrate atoms, the committed frame has %d. D22.2 GATE 3(c)."
             % (src, st, len(rows), len(ref)))
    worst, where = 0.0, None
    for i, (a, b) in enumerate(zip(rows, ref)):
        if a[0] != b[0]:
            stop("%s_%s atom %d is %s, the committed frame has %s. D22.2 GATE 3(c)."
                 % (src, st, i + 1, a[0], b[0]))
        for j in range(1, 4):
            if abs(a[j] - b[j]) > worst:
                worst, where = abs(a[j] - b[j]), i + 1
    if worst > FRAME_TOL:
        stop("%s_%s departs from the committed frame by %.2e A at atom %d, over the %.0e A\n"
             "   D22.2 GATE 3(c) allows. A walled arm on a different geometry from its bare arm\n"
             "   makes d_shield meaningless." % (src, st, worst, where, FRAME_TOL))
    return worst


def make_input(src, st, src_text):
    """The one substitution: the production head prefix swapped for the stability head."""
    h = TR.head(BASIS)
    if not src_text.startswith(h):
        stop("%s_%s does not begin with cs6_transfer.head(%r), so the one substitution D22.2\n"
             "   defines cannot be made without editing its body, which this does not do."
             % (src, st, BASIS))
    return D21.stab_head() + src_text[len(h):]


def gate_reduces(src, st, mine, src_text):
    """D22.2 GATE 3(a): reversing the substitution must give back the committed bytes."""
    s = D21.stab_head()
    if not mine.startswith(s):
        stop("%s_%s does not begin with the stability head" % (src, st))
    back = TR.head(BASIS) + mine[len(s):]
    if back != src_text:
        n = min(len(back), len(src_text))
        k = next((i for i in range(n) if back[i] != src_text[i]), n)
        stop("%s_%s does not reduce to its committed input BYTE FOR BYTE. First difference at\n"
             "   byte %d. The only difference between them may be STABPerform. D22.2 GATE 3(a)."
             % (src, st, k))


def rep1(t, old, new, what):
    if t.count(old) != 1:
        stop("the PBS template carries %r %d times, not once, so %s cannot be set without\n"
             "   retyping the template, which this does not do." % (old, t.count(old), what))
    return t.replace(old, new)


def d21_bare(T6):
    """READING 2 reuses DIAGNOSTIC 21's bare verdicts. They are READ, never re-run."""
    p = T6 / D21_RESULT
    if not p.is_file():
        stop("%s is not there; READING 2 reuses DIAGNOSTIC 21's bare verdicts and cannot be\n"
             "   assembled without the committed result." % p)
    t = p.read_text(errors="replace")
    m = re.search(r"READING 1,.*?\n(.*?)\n\s*VERDICT:\s*([^\n]*)", t, re.S)
    if not m:
        stop("could not find DIAGNOSTIC 21's READING 1 block in %s" % p)
    out = {}
    for nm, st, v in re.findall(r"^\s+(\S+)\s+(R|TS)\s+\d+\.\d+\s+(STABLE|UNSTABLE)\b",
                                m.group(1), re.M):
        out[(nm, st)] = v
    if len(out) != 18:
        stop("DIAGNOSTIC 21's READING 1 parsed to %d bare verdicts, not the 18 it reports.\n"
             "   This is a PARSER problem, not a result problem, and re-reading costs nothing."
             % len(out))
    missing = sorted("%s %s" % (nm, st) for nm, _, _, _, _ in COND for st in Q13.GEOM
                     if (nm, st) not in out)
    if missing:
        stop("DIAGNOSTIC 21 carries no bare verdict for: %s" % ", ".join(missing))
    if m.group(2).strip() != "ALL STABLE":
        print("   NOTE: DIAGNOSTIC 21's READING 1 verdict is %r, not ALL STABLE; READING 2 "
              "reports its per-arm verdicts as committed." % m.group(2).strip())
    return out


def prepare(root, out, hpc):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    out = Path(out)
    if out.exists():
        stop("%s already exists; this never overwrites" % out)
    if not (T6 / PREREG).is_file():
        stop("%s is not committed; nothing is prepared before its pre-registration exists"
             % (T6 / PREREG))

    D21.gate1(T6)                      # GATE 1, reused from D20 and not repeated
    gate_targets(T6)                   # GATE 2
    geo = frames(T6)                   # the reference for GATE 3(c)

    runs, worst = [], 0.0
    for src, arm, serves in arms():
        for st in Q13.GEOM:
            sp = T6 / ("%s_%s" % (src, st)) / "job.inp"
            if not sp.is_file():
                stop("the committed source input is missing: %s" % sp)
            src_text = sp.read_text()
            gate_selfcontained(src, st, src_text)
            worst = max(worst, gate_geometry(src, st, src_text, geo))
            mine = make_input(src, st, src_text)
            gate_reduces(src, st, mine, src_text)
            d = out / ("%s_%s" % (tag(src), st))
            d.mkdir(parents=True)
            (d / "job.inp").write_text(mine)
            runs.append("%s_%s" % (tag(src), st))
    if len(runs) != 30:
        stop("%d runs were built, expected the 30 D22.2 states" % len(runs))

    pbs = AD.PBS % dict(tag="22", hpc=hpc, basis=".", n=len(runs),
                        runs="\n".join('run "%s"' % r for r in runs))
    pbs = rep1(pbs, "cm_d1422", "cm_d22", "the job name")
    pbs = rep1(pbs, "DIAGNOSTIC 14 (cs6_additivity.py)", "DIAGNOSTIC 22 (cs6_d22.py)",
               "the comment line")
    pbs = rep1(pbs, "%d single points at .." % len(runs),
               "%d single points with STABPerform, at %s." % (len(runs), BASIS),
               "the description")
    pbs = rep1(pbs, "%s/./" % hpc, "%s/" % hpc, "the output path")
    pbs = rep1(pbs, "W=%s/." % hpc, "W=%s" % hpc, "the working directory")
    if re.search(r"^#PBS([ \t].*)?[ \t]-q([ \t]|$)", pbs, re.M):
        stop("the PBS script requests a queue; the allocation on hpc1 is automatic")
    (out / "run.pbs").write_text(pbs)

    print("DIAGNOSTIC 22 prepared at %s" % out)
    print("   %d runs, %s: the walled and wall-at-zero legs of d_shield" % (len(runs), BASIS))
    for nm, n, f, w, w0 in COND:
        print("      %-13s n=%-3d sum|q| %.4f e   walled %-34s wall0 %s" % (nm, n, f, w, w0))
    print("   gate 1     DIAGNOSTIC 20's %s control is committed and STABLE" % BASIS)
    print("   gate 2     all %d targets read from committed files agree with D22.2"
          % (2 * len(TARGETS)))
    print("   gate 3(a)  every input reduces to its source's committed input BYTE FOR BYTE")
    print("   gate 3(b)  no source input references an external geometry or charge file")
    print("   gate 3(c)  every inline geometry matches the committed frame; worst %.2e A" % worst)
    return 0


def analyse(root, runs, outdir):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    runs = Path(runs)
    H = TR.H

    D21.gate1(T6)
    T = gate_targets(T6)
    BARE = d21_bare(T6)

    S, SERVES = {}, {}
    for src, arm, serves in arms():
        SERVES[src] = (arm, serves)
        for st in Q13.GEOM:
            S[(src, st)] = D19.read_stab(runs / ("%s_%s" % (tag(src), st)) / "job.out")

    L = []; w = L.append
    w("Stage 6, DIAGNOSTIC 22 RESULT.")
    w("Read by cs6_d22.py against %s. Criteria as pre-registered." % PREREG)
    w("")
    w("GATE 1, the %s control, reused from DIAGNOSTIC 20 and not repeated" % BASIS)
    w("   ctl_R   STABLE  E(0) %+.8f      ctl_TS  STABLE  E(0) %+.8f"
      % (D21.D20_CTL_E0["R"], D21.D20_CTL_E0["TS"]))
    w("   GATE 1 passed: verified present and STABLE in the committed %s." % D21.D20_RESULT)
    w("")

    w("GATE 2, each run's FIRST converged energy against its own source's committed energy")
    w("   %-34s %-3s %19s %19s %11s" % ("source", "st", "first converged", "committed", "diff Eh"))
    g2, first_e = {}, {}
    for src, arm, serves in arms():
        for st in Q13.GEOM:
            an, tot, fsp = S[(src, st)]
            disp = fsp[-1] - tot[-1]
            first = tot[0] + disp
            ok = abs(first - T[(src, st)]) <= GATE_EH
            g2[(src, st)] = ok
            first_e[(src, st)] = first
            w("   %-34s %-3s %19.9f %19.9f %11.2e%s"
              % (src, st, first, T[(src, st)], first - T[(src, st)], "" if ok else "   FAILED"))
    nbad = sum(1 for v in g2.values() if not v)
    if nbad:
        w("   GATE 2 FAILED on %d run(s): STABPerform changed their SCF, so they are not the" % nbad)
        w("   production solution and READING 4 is withheld for every condition they serve.")
    else:
        w("   GATE 2 passed: the first solution IS the production solution, everywhere.")
    w("")

    w("READING 1, are the walled and wall-at-zero arms stable?")
    w("   %-6s %-34s %-3s %-9s %s" % ("arm", "source", "st", "verdict", "E(0..2)"))
    unstable = []
    for src, arm, serves in arms():
        for st in Q13.GEOM:
            an, tot, fsp = S[(src, st)]
            v, eig = an[0]
            if v == "UNSTABLE":
                unstable.append((src, arm, st))
            w("   %-6s %-34s %-3s %-9s %s" % (arm, src, st, v,
              "  ".join("%+.8f" % x for x in eig[:3])))
            if len(an) > 1:
                w("      ORCA reconverged; second analysis %s, E(0) %+.8f, total energy "
                  "%+.4f kcal/mol" % (an[-1][0], an[-1][1][0], (tot[-1] - tot[0]) * H))
    v1 = "ALL STABLE" if not unstable else "SOME UNSTABLE"
    w("   VERDICT: %s" % v1)
    if unstable:
        w("   unstable: %s" % ", ".join("%s %s %s" % (a, s, t) for s, a, t in unstable))
    w("")

    w("READING 2, the d_shield arm audit")
    w("   bare verdicts are DIAGNOSTIC 21's, read from its committed result and not re-run.")
    w("   the control leg is GATE 1 above, verified and not re-run.")
    w("   %-13s %-8s %-8s %-10s %-10s %-9s %-9s %s"
      % ("condition", "bare R", "bare TS", "walled R", "walled TS", "wall0 R", "wall0 TS", "verdict"))
    audit = {}
    for nm, n, f, wsrc, w0src in COND:
        cells = [BARE[(nm, "R")], BARE[(nm, "TS")],
                 S[(wsrc, "R")][0][0][0], S[(wsrc, "TS")][0][0][0],
                 S[(w0src, "R")][0][0][0], S[(w0src, "TS")][0][0][0]]
        ok = all(c == "STABLE" for c in cells)
        audit[nm] = ok
        if nm in NO_DSET:
            lbl = "ARMS STABLE, NO READABLE Dset" if ok else "NOT VERIFIED"
        else:
            lbl = "ARMS VERIFIED" if ok else "NOT VERIFIED"
        w("   %-13s %-8s %-8s %-10s %-10s %-9s %-9s %s"
          % (nm, cells[0], cells[1], cells[2], cells[3], cells[4], cells[5], lbl))
        if nm in NO_DSET:
            w("      %s" % NO_DSET[nm])
    v2 = "ALL ARMS VERIFIED" if all(audit.values()) else "ARMS NOT VERIFIED"
    w("   VERDICT: %s" % v2)
    read = [nm for nm, _, _, _, _ in COND if nm not in NO_DSET]
    if all(audit.values()):
        w("   The %d conditions that HAVE a readable d_shield now rest on verified-stable"
          % len(read))
        w("   solutions in all four legs - control, bare, walled and wall-at-zero:")
        w("      %s" % ", ".join(read))
        w("   The other %d have stable arms and NO readable Dset, per ERRATUM 1 E22.1. Their"
          % len(NO_DSET))
        w("   arms being minima does not make a quantity quotable that D6.3 withheld:")
        for nm in sorted(NO_DSET):
            w("      %-13s %s" % (nm, NO_DSET[nm]))
    w("")

    w("READING 3, the replicate: two committed runs with byte-identical input")
    a, b = REPLICATE
    w("   %s" % a)
    w("   %s" % b)
    agree = True
    for st in Q13.GEOM:
        va, ea = S[(a, st)][0][0]
        vb, eb = S[(b, st)][0][0]
        tol = max(REP_REL * max(abs(ea[0]), abs(eb[0])), REP_ABS)
        d = abs(ea[0] - eb[0])
        same = (va == vb) and (d <= tol)
        agree &= same
        w("   %-3s %-9s vs %-9s   E(0) %+.8f vs %+.8f   |diff| %.2e   tol %.2e   %s"
          % (st, va, vb, ea[0], eb[0], d, tol, "agree" if same else "DISAGREE"))
    v3 = "REPLICATE AGREES" if agree else "REPLICATE DISAGREES"
    w("   VERDICT: %s" % v3)
    if not agree:
        w("   D22.3 states the consequence in advance: this is a finding about how reproducibly a")
        w("   stability eigenvalue is computed - a property of the diagnostic - and not about the")
        w("   chemistry. The two runs are the same calculation.")
    w("")

    w("READING 4, what does any instability cost d_shield?")
    if not unstable:
        w("   No arm is unstable, so there is nothing to correct and nothing is reported here.")
        v4 = "NOTHING TO CORRECT"
    else:
        bad_src = {s for s, _, _ in unstable}
        BAREE = D21.gate_targets(T6)    # committed bare energies, full precision, D21-gated

        def leg(src, st, use_stable):
            an, tot, fsp = S[(src, st)]
            disp = fsp[-1] - tot[-1]
            return (tot[-1] if use_stable else tot[0]) + disp

        def ddE_bare(nm):
            return ((BAREE[(nm, "TS")] - D21.CTL["TS"])
                    - (BAREE[(nm, "R")] - D21.CTL["R"])) * H

        def ddE_wall(wsrc, w0src, use_stable):
            z = {st: leg(w0src, st, use_stable and w0src in bad_src) - D21.CTL[st]
                 for st in Q13.GEOM}
            d = {st: leg(wsrc, st, use_stable and wsrc in bad_src) - D21.CTL[st] - z[st]
                 for st in Q13.GEOM}
            return (d["TS"] - d["R"]) * H

        w("   d_shield = ddE(bare) - ddE(wall), formed exactly as D22.1 states it. The bare arm")
        w("   is NOT re-run: its committed energy is used unchanged and DIAGNOSTIC 21 verified")
        w("   it stable. Only the arm that came back unstable moves.")
        w("   %-13s %11s %11s %11s %11s %9s %s"
          % ("condition", "d_shield", "corrected", "change", "budget", "% budget", "moved by"))
        rows, withheld = [], []
        for nm, n, f, wsrc, w0src in COND:
            touched = [s for s in (wsrc, w0src) if s in bad_src]
            if not touched:
                continue
            if not all(g2[(s, st)] for s in (wsrc, w0src) for st in Q13.GEOM):
                withheld.append(nm); continue
            a0 = ddE_bare(nm) - ddE_wall(wsrc, w0src, False)
            a1 = ddE_bare(nm) - ddE_wall(wsrc, w0src, True)
            lbl = ", ".join("walled" if s == wsrc else "wall-at-zero" for s in touched)
            rows.append((nm, f, a0, a1, a1 - a0, lbl))
        for nm, f, a0, a1, dd, lbl in rows:
            w("   %-13s %11.4f %11.4f %11.4f %11.4f %9.1f %s"
              % (nm, a0, a1, dd, BUDGET, 100.0 * abs(dd) / BUDGET, lbl))
            if nm in NO_DSET:
                w("      ERRATUM 1 E22.3(c): this correction applies to a quantity that is "
                  "already")
                w("      NOT READABLE, so it bounds nothing this project quotes.")
        for nm in withheld:
            w("   %-13s withheld, GATE 2 failed on one of its arms" % nm)
        if not rows and not withheld:
            stop("READING 1 found an unstable arm but no condition is served by it; the COND "
                 "table and READING 1 disagree, which is a code fault, not a result.")
        over = [r for r in rows if abs(r[4]) > BUDGET]
        v4 = "CORRECTION MATERIAL" if over else "CORRECTION WITHIN BUDGET"
        if over:
            w("")
            w("   D22.4's stated consequence applies. These conditions' d_shield may not be quoted")
            w("   without the correction attached:")
            for nm, f, a0, a1, dd, _lbl in sorted(over):
                w("      %-13s sum|q| %.4f e   change %+.4f kcal/mol" % (nm, f, dd))
    w("   VERDICT: %s" % v4)
    w("")

    w("D22.5 stands whatever was returned. This is a restricted-to-restricted test: a STABLE")
    w("verdict does not make the system single-reference, and DIAGNOSTIC 19's READING 3 returned")
    w("DISAGREE, which stands. ORCA's restart is a bounded heuristic, so any correction above is")
    w("a LOWER BOUND. Nine conditions at %s and nothing else; sum|q| above 2.19 e is untouched." % BASIS)
    w("Stability is a necessary condition on the arithmetic, not evidence that the wall is the")
    w("right representation of Pauli repulsion. ERRATUM 7 is not revisited.")
    w("")
    w("VERDICTS: %s | %s | %s | %s" % (v1, v2, v3, v4))

    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    p = outdir / "DIAGNOSTIC22_RESULT.txt"
    p.write_text("\n".join(L) + "\n")
    print("\n".join(L))
    print()
    print("written to %s" % p)
    return 0


def main(argv):
    if len(argv) != 5:
        stop(__doc__.strip())
    if argv[1] == "prepare":
        return prepare(argv[2], argv[3], argv[4])
    if argv[1] == "analyse":
        return analyse(argv[2], argv[3], argv[4])
    stop("unknown mode %r" % argv[1])


if __name__ == "__main__":
    sys.exit(main(sys.argv))
