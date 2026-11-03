#!/usr/bin/env python3
"""
cs6_d23.py - STAGE6_DIAGNOSTIC23_20261102.txt.

Does the choice of pseudopotential parameter set reach d_shield?

  python3 cs6_d23.py prepare <results_repo> <OUT> <HPC_DIR>

NOTHING ALREADY COMMITTED IS RETYPED OR RE-RUN.
  cs6_transfer    head, H
  cs6_additivity  PBS
  cs13_qscale     FRAME_SHA, GEOM
  cs6_d21         gate1
  cs6_d22         tag, sha, rep1, gate_selfcontained, gate_geometry, frames, COND

THE ONE SUBSTITUTION. Every job.inp is COPIED VERBATIM from the committed run it measures,
with exactly one substitution: the Ne ecp-all parameter block swapped for the He or Ar block
of the same paper's Tables S1 and S3. The charge and the coordinates on each Ne> line are
untouched, the head, charge, multiplicity and geometry are untouched, and a gate requires that
reversing the substitution gives the committed bytes back BYTE FOR BYTE. This is DIAGNOSTIC 22's
construction, for the reason D22.2 records.

The substitution is two literal string replacements, not a reparse: the parameter block body is
byte-identical in every centre of every committed walled input (verified, GATE 4(b)), so the
whole operation is `body -> body_X` and `  Ne>  ->   X>  `, each exactly reversible.

THE CONTROL IS NOT SUBSTITUTED. It carries no ECP. Its two runs are copied verbatim and exist
to gate the allocation against the committed control at read time, as D23.3 GATE 2 states. The
control cancels exactly from both differences D23.4 READING 2 forms, so it enters no arithmetic.
"""
import hashlib, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cs13_qscale as Q13
import cs6_transfer as TR
import cs6_additivity as AD
import cs6_d21 as D21
import cs6_d22 as D22

PREREG = "STAGE6_DIAGNOSTIC23_20261102.txt"
SICHK = "STAGE6_DIAGNOSTIC23_SI_DIGIT_CHECK_20261102.txt"
BASIS = "def2-SVPD"
CTL = "d16/def2-SVPD/ctl"          # the charge-free control, copied verbatim

# D23.2: the seven readable conditions. Walled source, then wall-at-zero source.
# The wall-at-zero column names FIVE distinct position sets for seven conditions, which is
# what makes the per-set cost 24 rather than 28; GATE 4(c) proves they are distinct and
# GATE 4(d) proves the two S6-10 conditions really do share a position set.
COND = (
    ("S6-2_D17",     2, 1.0777, "d17/def2-SVPD/S6-2_wall_M",       "d16/def2-SVPD/S6-2_wall0"),
    ("S6-3_D16",     3, 1.3000, "d16/def2-SVPD/S6-3_wall_F1",      "d16/def2-SVPD/S6-3_wall0"),
    ("S6-5_D16",     5, 1.3000, "d16/def2-SVPD/S6-5_wall_F1",      "d16/def2-SVPD/S6-5_wall0"),
    ("S6-10_D16",   10, 1.3000, "d16/def2-SVPD/S6-10_wall_F1",     "d16/def2-SVPD/S6-10_wall0"),
    ("S6-10_D15",   10, 2.1900, "additivity/def2-SVPD/n10_wall_A", "d16/def2-SVPD/S6-10_wall0"),
    ("S6-20_D16F1", 20, 1.3000, "d16/def2-SVPD/S6-20_wall_F1",     "d16/def2-SVPD/S6-20_wall0"),
    ("S6-20_D16F2", 20, 2.1900, "d16/def2-SVPD/S6-20_wall_F2",     "d16/def2-SVPD/S6-20_wall0"),
)
# D22 referenced S6-10_D15's wall-at-zero arm to additivity/def2-SVPD/n10_wall0, which is a
# SECOND committed run at the SAME positions. D23 needs one run per position set, so it uses
# the d16 arm for both S6-10 conditions and GATE 4(d) checks that the two committed arms
# really do carry identical positions before that substitution of sources is allowed.
D22_N10_W0 = "additivity/def2-SVPD/n10_wall0"

# The committed Ne block, byte for byte. Not retyped: GATE 4(a) reads it out of the committed
# inputs and stops unless every centre of every source carries exactly this.
NE_BODY = ("  N_core 10\n  lmax f\n  s 1\n    1  2.0475   54.5100  2\n"
           "  p 1\n    1  0.44815   1.46500  2\n  d 1\n    1  0.49205  -0.83800  2\n"
           "  f 1\n    1  1.0000    0.0000  2\nend\n")
NE_PREFIX = "  Ne>  "

# D23.2, as the pre-registration prints them and as STAGE6_DIAGNOSTIC23_SI_DIGIT_CHECK
# verifies them against SI Tables S1 and S3.
SETS = {
    "he": dict(
        el="He", prefix="  He>  ",
        body=("  N_core 2\n  lmax d\n  s 1\n    1  0.50980   2.42000  2\n"
              "  p 1\n    1  0.49165  -0.43590  2\n  d 1\n    1  1.00000   0.00000  2\nend\n"),
        digits=["0.50980", "2.42000", "0.49165", "-0.43590", "1.00000", "0.00000"],
        ncore=2, lmax="d"),
    "ar": dict(
        el="Ar", prefix="  Ar>  ",
        body=("  N_core 18\n  lmax f\n  s 1\n    1  1.64100  275.00000  2\n"
              "  p 1\n    1  0.27330   1.90000  2\n  d 1\n    1  0.44000  -3.40000  2\n"
              "  f 1\n    1  1.00000   0.00000  2\nend\n"),
        digits=["1.64100", "275.00000", "0.27330", "1.90000", "0.44000", "-3.40000",
                "1.00000", "0.00000"],
        ncore=18, lmax="f"),
}


def stop(m):
    print("STOP: " + m)
    sys.stdout.flush()
    sys.exit(1)


def centres(text):
    """Every (Ne> line, block body) pair in a committed input, in order."""
    return re.findall(r"(^[ ]*\w+>[^\n]*NewECP\n)((?:.*?\n)*?end\n)", text, re.M)


def positions(text):
    """The coordinate triples of a committed input's ECP centres, as printed."""
    return tuple(re.findall(r"^[ ]*\w+>\s+-?[\d.]+\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)",
                            text, re.M))


def gate_ne_block(src, st, text):
    """GATE 4(a). Every centre must carry exactly the committed Ne block, or the one
    substitution below is not the one substitution it claims to be."""
    cs = centres(text)
    if not cs:
        stop("%s_%s carries no ECP centre, so there is nothing for D23 to substitute."
             % (src, st))
    for i, (line, body) in enumerate(cs):
        if body != NE_BODY:
            stop("centre %d of %s_%s does not carry the committed Ne block byte for byte, so a\n"
                 "   literal substitution would change something other than the parameters.\n"
                 "   expected %r\n   found    %r. D23 GATE 4(a)."
                 % (i + 1, src, st, NE_BODY[:48], body[:48]))
        if not line.startswith(NE_PREFIX):
            stop("centre %d of %s_%s does not begin %r, so swapping the element label would\n"
                 "   change its alignment. D23 GATE 4(a)." % (i + 1, src, st, NE_PREFIX))
    if text.count(NE_BODY) != len(cs) or text.count(NE_PREFIX) != len(cs):
        stop("%s_%s has %d centres but %d block bodies and %d %r prefixes; a literal\n"
             "   substitution would not be one to one. D23 GATE 4(a)."
             % (src, st, len(cs), text.count(NE_BODY), text.count(NE_PREFIX), NE_PREFIX))
    return len(cs)


def make_input(src, st, src_text, S):
    """The one substitution: the Ne parameter block swapped for this set's block."""
    out = src_text.replace(NE_BODY, S["body"]).replace(NE_PREFIX, S["prefix"])
    if NE_BODY in out or NE_PREFIX in out:
        stop("%s_%s still carries a Ne block after substitution. D23 GATE 4(a)." % (src, st))
    return out


def gate_reduces(src, st, mine, src_text, S):
    """GATE 1. Reversing the substitution must give back the committed bytes."""
    back = mine.replace(S["prefix"], NE_PREFIX).replace(S["body"], NE_BODY)
    if back != src_text:
        n = min(len(back), len(src_text))
        k = next((i for i in range(n) if back[i] != src_text[i]), n)
        stop("%s_%s at %s does not reduce to its committed input BYTE FOR BYTE. First\n"
             "   difference at byte %d. The only difference between them may be the\n"
             "   pseudopotential parameters. D23.3 GATE 1." % (src, st, S["el"], k))


def gate_body_digits(S):
    """GATE 4(g). The block body must actually contain every digit GATE 4(e) checks against
    the pre-registration. Without this the body and the digits list can drift apart and the
    inputs would carry parameters that nothing verified."""
    for d in S["digits"]:
        if d not in S["body"]:
            stop("%s's block body does not contain %s, which GATE 4(e) checks against the\n"
                 "   pre-registration. The body and the checked digits have drifted apart.\n"
                 "   D23 GATE 4(g)." % (S["el"], d))
    body_nums = set(re.findall(r"-?\d+\.\d+", S["body"]))
    listed = set(S["digits"])
    extra = {x for x in body_nums if x not in listed}
    if extra:
        stop("%s's block body carries %s, which GATE 4(e) does not check against the\n"
             "   pre-registration. Every parameter in the body must be checked.\n"
             "   D23 GATE 4(g)." % (S["el"], sorted(extra)))
    return len(listed)


def gate_si_digits(T6, S):
    """GATE 4(e). Every digit of this block appears in the committed SI digit check and in
    the committed pre-registration. The blocks are not retyped from the paper here."""
    for rel in (PREREG, SICHK):
        p = T6 / rel
        if not p.is_file():
            stop("%s is not committed; D23 prepares nothing before its pre-registration and\n"
                 "   its SI digit check exist." % rel)
        t = p.read_text()
        for d in S["digits"]:
            if d not in t:
                stop("the %s block's %s is not in %s, so this block is not the one that was\n"
                     "   pre-registered and checked against the paper. D23 GATE 4(e)."
                     % (S["el"], d, rel))
        if ("N_core %d" % S["ncore"]) not in t or ("lmax %s" % S["lmax"]) not in t:
            stop("%s's N_core or lmax is not in %s. D23 GATE 4(e)." % (S["el"], rel))
    m = re.search(r"primitives disagreeing with the SI page:\s*(\d+)", (T6 / SICHK).read_text())
    if not m or int(m.group(1)) != 0:
        stop("the committed SI digit check does not report zero disagreements, so the blocks\n"
             "   are not the published ones. D23 GATE 4(e).")


def arms():
    """Every (source, arm, serves) D23.2 requires: walled per condition, wall-at-zero per
    distinct position set."""
    walled = [(w, "walled", [nm]) for nm, n, f, w, w0 in COND]
    zero = {}
    for nm, n, f, w, w0 in COND:
        zero.setdefault(w0, []).append(nm)
    out = walled + [(k, "wall0", v) for k, v in sorted(zero.items())]
    if len(walled) != 7:
        stop("%d walled sources, expected the 7 readable conditions" % len(walled))
    if len(zero) != 5:
        stop("%d distinct wall-at-zero sources, but D23.2 names five position sets" % len(zero))
    names = [D22.tag(s) for s, _, _ in out]
    if len(set(names)) != len(names):
        stop("two sources collide on the same run directory name: %s" % names)
    return out


def prepare(root, out, hpc):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    out = Path(out)
    if out.exists():
        stop("%s already exists; this never overwrites" % out)
    for rel in (PREREG, SICHK):
        if not (T6 / rel).is_file():
            stop("%s is not committed; nothing is prepared before it exists" % (T6 / rel))

    D21.gate1(T6)                       # the committed def2-SVPD control is committed and STABLE
    geo = D22.frames(T6)                # the reference for the geometry gate
    for key in sorted(SETS):
        gate_body_digits(SETS[key])     # GATE 4(g), before 4(e) checks the same digits
        gate_si_digits(T6, SETS[key])   # GATE 4(e)

    # GATE 4(c)/(d): the five position sets are distinct, and the two S6-10 conditions really
    # do share one, which is what licenses running the d16 arm for both.
    seen, nset = {}, 0
    for src, arm, serves in arms():
        if arm != "wall0":
            continue
        t = (T6 / ("%s_R" % src) / "job.inp").read_text()
        pos = positions(t)
        if pos in seen:
            stop("%s and %s carry the same positions, so D23.2's count of five distinct\n"
                 "   wall-at-zero sets is wrong. D23 GATE 4(c)." % (src, seen[pos]))
        seen[pos] = src
        nset += 1
    if nset != 5:
        stop("%d wall-at-zero position sets, not the five D23.2 names. D23 GATE 4(c)." % nset)
    a = positions((T6 / "d16/def2-SVPD/S6-10_wall0_R" / "job.inp").read_text())
    b = positions((T6 / ("%s_R" % D22_N10_W0) / "job.inp").read_text())
    if a != b:
        stop("d16/def2-SVPD/S6-10_wall0 and %s do NOT carry identical positions, so one run\n"
             "   cannot serve both S6-10 conditions and D23.2's cost of 24 per set is wrong.\n"
             "   D23 GATE 4(d)." % D22_N10_W0)

    runs, worst, ncent = [], 0.0, 0
    for key in sorted(SETS):
        S = SETS[key]
        for src, arm, serves in arms():
            for st in Q13.GEOM:
                sp = T6 / ("%s_%s" % (src, st)) / "job.inp"
                if not sp.is_file():
                    stop("the committed source input is missing: %s" % sp)
                src_text = sp.read_text()
                D22.gate_selfcontained(src, st, src_text)
                worst = max(worst, D22.gate_geometry(src, st, src_text, geo))
                ncent += gate_ne_block(src, st, src_text)          # GATE 4(a)
                mine = make_input(src, st, src_text, S)
                gate_reduces(src, st, mine, src_text, S)           # GATE 1
                d = out / ("%s_%s_%s" % (key, D22.tag(src), st))
                d.mkdir(parents=True)
                (d / "job.inp").write_text(mine)
                runs.append(d.name)

    # The control, copied verbatim with no substitution. Unlike a walled arm it is NOT
    # self-contained: it names its geometry as an external xyzfile, so its companion file is
    # copied too and both are gated by sha256. GATE 4(f).
    for st in Q13.GEOM:
        srcdir = T6 / ("%s_%s" % (CTL, st))
        sp = srcdir / "job.inp"
        if not sp.is_file():
            stop("the committed control input is missing: %s" % sp)
        ct = sp.read_text()
        if "NewECP" in ct or re.search(r"^[ ]*\w+>", ct, re.M):
            stop("the control input carries an ECP centre; it is meant to be charge free.\n"
                 "   D23 GATE 4(f).")
        if "%pointcharges" in ct:
            stop("the control input references a point-charge file; it is meant to carry no\n"
                 "   charge at all. D23 GATE 4(f).")
        m = re.search(r"^\* xyzfile\s+-2\s+1\s+(\S+)\s*$", ct, re.M)
        if not m:
            stop("the control input does not name its geometry as '* xyzfile -2 1 <file>', so\n"
                 "   the companion file to copy cannot be identified. D23 GATE 4(f).")
        comp = srcdir / m.group(1)
        if not comp.is_file():
            stop("the control names %r but that file is not beside its committed input.\n"
                 "   D23 GATE 4(f)." % m.group(1))
        # its geometry must be the committed frame, checked on the xyz rather than inline
        rows = []
        for line in comp.read_text().splitlines()[2:]:
            q = line.split()
            if len(q) >= 4:
                rows.append((q[0], float(q[1]), float(q[2]), float(q[3])))
        ref = geo[st]
        if len(rows) != len(ref):
            stop("the control's %r has %d atoms, the committed frame has %d. D23 GATE 4(f)."
                 % (m.group(1), len(rows), len(ref)))
        for i, (aa, bb) in enumerate(zip(rows, ref)):
            if aa[0] != bb[0]:
                stop("the control's atom %d is %s, the committed frame has %s. D23 GATE 4(f)."
                     % (i + 1, aa[0], bb[0]))
            for j in range(1, 4):
                worst = max(worst, abs(aa[j] - bb[j]))
        d = out / ("ctl_%s" % st)
        d.mkdir(parents=True)
        for f in (sp, comp):
            (d / f.name).write_bytes(f.read_bytes())
            if hashlib.sha256((d / f.name).read_bytes()).hexdigest() != D22.sha(f):
                stop("the control copy of %s is not byte-identical to its committed file.\n"
                     "   D23 GATE 4(f)." % f.name)
        runs.append(d.name)

    if len(runs) != 50:
        stop("%d runs were built, expected the 50 D23.2 states" % len(runs))

    pbs = AD.PBS % dict(tag="23", hpc=hpc, basis=".", n=len(runs),
                        runs="\n".join('run "%s"' % r for r in runs))
    pbs = D22.rep1(pbs, "cm_d1423", "cm_d23", "the job name")
    pbs = D22.rep1(pbs, "DIAGNOSTIC 14 (cs6_additivity.py)", "DIAGNOSTIC 23 (cs6_d23.py)",
                   "the comment line")
    pbs = D22.rep1(pbs, "%d single points at .." % len(runs),
                   "%d single points at %s: the He and Ar parameter sets through both arms."
                   % (len(runs), BASIS), "the description")
    pbs = D22.rep1(pbs, "%s/./" % hpc, "%s/" % hpc, "the output path")
    pbs = D22.rep1(pbs, "W=%s/." % hpc, "W=%s" % hpc, "the working directory")
    if re.search(r"^#PBS([ \t].*)?[ \t]-q([ \t]|$)", pbs, re.M):
        stop("the PBS script requests a queue; the allocation on hpc1 is automatic")
    (out / "run.pbs").write_text(pbs)

    lines = []
    for d in sorted(out.iterdir()):
        if not d.is_dir():
            continue
        for f in sorted(d.iterdir()):
            lines.append("%s  %s/%s" % (hashlib.sha256(f.read_bytes()).hexdigest(), d.name, f.name))
    lines.append("%s  run.pbs" % hashlib.sha256((out / "run.pbs").read_bytes()).hexdigest())
    (out / "INPUT_CHECKSUMS.sha256").write_text("\n".join(lines) + "\n")

    print("DIAGNOSTIC 23 prepared at %s" % out)
    print("   %d runs at %s: %d walled, %d wall-at-zero, 2 charge-free control"
          % (len(runs), BASIS, 2 * 2 * 7, 2 * 2 * 5))
    for key in sorted(SETS):
        S = SETS[key]
        print("      %-3s N_core %-2d lmax %s   s/p%s channels from SI Table %s"
              % (S["el"], S["ncore"], S["lmax"], "/d" if S["lmax"] == "f" else "",
                 "S1" if S["el"] == "He" else "S3"))
    print("   gate 1     every input reduces to its committed input BYTE FOR BYTE")
    print("   gate 4(a)  all %d centres carried exactly the committed Ne block" % ncent)
    print("   gate 4(c)  5 distinct wall-at-zero position sets")
    print("   gate 4(d)  the two S6-10 conditions share one position set")
    print("   gate 4(e)  every block digit is in the committed pre-registration and SI check")
    print("   gate 4(f)  the control carries no charge and its geometry file is copied by sha")
    print("   gate 4(g)  every parameter in each block body is one GATE 4(e) checks")
    print("   geometry   every inline geometry matches the committed frame; worst %.2e A" % worst)
    return 0




# ======================================================================== READ
#
# D23.4 READINGS 1 to 5 and D23.5's criteria. Added after the fifty runs
# returned; nothing above this line was changed when it was added.

import cs6_audit as AUD
import cs6_negatives as NEG
import cs6_transfer as TRX

GATE_EH = 1e-6                 # D23.3 GATE 2, the figure cs6_d21.py already uses
SHIELD_ABS = 0.05              # D23.5, from cs11_standoff.py
DZ_FLOOR = 0.001               # D23.4 READING 3 withholds f below this
TOL_VERDICT = 0.0500           # the tolerance all four committed verdicts carry

# D23.4 READING 2's table of committed Ne values, re-derived from the committed
# energies at read time and gated against these printed figures.
NE_COMMITTED = {
    "S6-2_D17":     (-0.0026, -0.3468),
    "S6-3_D16":     (-0.1815, -0.4379),
    "S6-5_D16":     (-0.4407, -0.4311),
    "S6-10_D16":    (-0.7187, -0.2430),
    "S6-10_D15":    (-0.7187, -0.4140),
    "S6-20_D16F1":  (-0.6190, -0.1702),
    "S6-20_D16F2":  (-0.6190, -0.2919),
}
# the committed bare arms, needed for d_shield and for READING 4's excesses
BARE = {
    "S6-2_D17": "d17/def2-SVPD/S6-2_bare_M", "S6-3_D16": "d16/def2-SVPD/S6-3_bare_F1",
    "S6-5_D16": "d16/def2-SVPD/S6-5_bare_F1", "S6-10_D16": "d16/def2-SVPD/S6-10_bare_F1",
    "S6-10_D15": "additivity/def2-SVPD/n10_bare_A",
    "S6-20_D16F1": "d16/def2-SVPD/S6-20_bare_F1", "S6-20_D16F2": "d16/def2-SVPD/S6-20_bare_F2",
}
# READING 4's committed verdicts, with the rule each result file prints.
VERDICTS_COMMITTED = {
    "D16 READING 1": ("FALLS WITH n", 0.0068, "Dset monotone over n = 3, 5, 10, 20 at 1.3000 e"),
    "D16 READING 2": ("CONSISTENT", 0.0721, "Dset(S6-20 F2) - Dset(S6-10 A), tol 0.0500"),
    "D16 READING 3": ("COOPERATIVITY FALLS", 0.0070, "excess(n=20) - excess(n=10), tol 0.0500"),
    "D15 READING 1": ("NOT BOUNDED", 0.0119, "Dset 0.4140 - Dsum 0.3521, tol 0.0500"),
}
# READING 5's analogue, from the committed def2-SVP / def2-SVPD pair on arm B.
ANALOGUE = {"n10": (-0.5988, -0.3426), "n6neg": (-0.7785, -0.2777)}
# the single sites READING 4's excesses need, as DIAGNOSTIC 16 READING 3 lists them
S20_F2_SITES = [673, 590, 2152, 1695, 762, 615, 1089, 1234, 610, 691,
                396, 751, 1392, 475, 52, 888, 308, 832, 967, 688]
N10_A_SITES = [1089, 1234, 1695, 2152, 400, 592, 605, 650, 652, 756]
SITE_TREES = ["additivity/def2-SVPD", "negatives/def2-SVPD", "negatives_armA/def2-SVPD",
              "d16/def2-SVPD", "d17/def2-SVPD"]


def fsp(d):
    p = Path(d) / "job.out"
    if not p.is_file():
        return None
    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", p.read_text(errors="replace"))
    return float(m[-1]) if m else None


def ddE_pair(eR, eTS, cR, cTS):
    """The double difference, in kcal/mol. Both arms referenced to the SAME control."""
    return ((eTS - cTS) - (eR - cR)) * TR.H


def committed_ddE(T6, tree, stem):
    t = T6 / tree
    cR, cTS = fsp(t / "ctl_R"), fsp(t / "ctl_TS")
    eR, eTS = fsp(t / ("%s_R" % stem)), fsp(t / ("%s_TS" % stem))
    if None in (cR, cTS, eR, eTS):
        return None
    return ddE_pair(eR, eTS, cR, cTS)


def find_committed(T6, stem):
    for tr in SITE_TREES:
        if (T6 / tr / ("%s_R" % stem)).is_dir():
            v = committed_ddE(T6, tr, stem)
            if v is not None:
                return v
    return None


# ------------------------------------------------------- the arithmetic, pure

def deltas(newE, neE, bare):
    """D23.4 READING 2 and 3, as pure arithmetic on energies.

    newE and neE are {(arm, state): E} for one condition, arm in {"wall", "wall0"};
    bare is the committed ddE(bare). Returns dz, dd and f.

    Both differences are formed so that THE CONTROL CANCELS EXACTLY: each is a
    difference between the X arm and the Ne arm at the same two states, so no
    control energy enters. D23.3 GATE 2 bounds how far that differs from forming
    each ddE against its own control and subtracting, and READING 2 reports both.
    """
    dz = ddE_pair(newE[("wall0", "R")], newE[("wall0", "TS")],
                  neE[("wall0", "R")], neE[("wall0", "TS")])
    draw = ddE_pair(newE[("wall", "R")], newE[("wall", "TS")],
                    neE[("wall", "R")], neE[("wall", "TS")])
    # d_shield = ddE(bare) - ddE(wall,raw) + ddE(wall0); the bare arm is shared and cancels
    dd = -draw + dz
    f = (dd / dz) if abs(dz) > DZ_FLOOR else None
    return dz, draw, dd, f


def deltas_via_committed(newz, newd, nez, ned):
    """The same two differences formed the way READING 2 words them: each side's own
    ddE against its own control, then subtracted. Reported beside deltas() so the
    control's contribution is visible rather than assumed."""
    return newz - nez, newd - ned


# ------------------------------------------------------- classification

def mem_of(inp_text):
    """The (charge, position) of every ECP centre of an input, for A1.2."""
    import numpy as np
    out = []
    for m in re.finditer(r"^[ ]*\w+>\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)",
                         inp_text, re.M):
        q = float(m.group(1))
        X = np.array([float(m.group(2)), float(m.group(3)), float(m.group(4))])
        out.append((None, q, X))
    return out


def gate_centres_are_Q(run, text, ncent, ctl_nb, ctl_nel):
    """D23.3 GATE 3. The centre must still contribute no basis function and no electron, and
    ORCA must still label it Q in the Loewdin block. If it labels an He> or Ar> centre by its
    element instead, cs6_audit.read would fold that row into max|Loewdin| and the whole
    classification would be wrong, silently. This is the one assumption the build could not
    check, so it is checked here before any class is read."""
    nb = re.search(r"Number of basis functions\s+\.+\s+(\d+)", text)
    nel = re.search(r"Number of Electrons\s+NEL\s+\.+\s+(\d+)", text)
    if not nb or int(nb.group(1)) != ctl_nb:
        stop("%s reports %s basis functions, the control reports %d. The centre is no longer\n"
             "   basis-free. D23.3 GATE 3." % (run, nb.group(1) if nb else "none", ctl_nb))
    if not nel or int(nel.group(1)) != ctl_nel:
        stop("%s reports NEL %s, the control reports %d. The centre is no longer electron-free.\n"
             "   D23.3 GATE 3." % (run, nel.group(1) if nel else "none", ctl_nel))
    k = text.rfind("LOEWDIN ATOMIC CHARGES")
    if k < 0:
        stop("%s prints no Loewdin block, so D6.3 cannot be applied. D23.3 GATE 3." % run)
    syms = []
    for l in text[k:].splitlines()[2:]:
        mm = re.match(r"\s*\d+\s+([A-Za-z]{1,2})\s*:\s*(-?\d+\.\d+)", l)
        if not mm:
            break
        syms.append(mm.group(1))
    nq = sum(1 for s in syms if s == AUD.ECP_SYMBOL)
    if nq != ncent:
        stop("%s has %d ECP centres but %d rows labelled %r in its Loewdin block. cs6_audit\n"
             "   excludes rows by that symbol, so a mislabelled centre would enter max|Loewdin|\n"
             "   and change the class. D23.3 GATE 3." % (run, ncent, nq, AUD.ECP_SYMBOL))
    for bad in ("He", "Ar"):
        if bad in syms:
            stop("%s's Loewdin block carries a row labelled %r. ORCA has labelled the coreless\n"
                 "   centre by its element rather than %r, so max|Loewdin| would include it.\n"
                 "   D23.3 GATE 3." % (run, bad, AUD.ECP_SYMBOL))
    ch = re.search(r"Sum of atomic charges\s*:\s*(-?\d+\.\d+)", text)
    return float(ch.group(1)) if ch else None


def excesses(T6, dsets):
    """READING 4's two excesses: |Dset| minus |signed sum of the single-site d_shield|.
    The single-site arms are the committed ones and are NOT re-run, so only the two
    collective Dsets move with the parameter set."""
    def site_ds(sid, arms):
        b = w = None
        for a in arms:
            b = b if b is not None else find_committed(T6, "s%d_bare_%s" % (sid, a))
            w = w if w is not None else find_committed(T6, "s%d_wall_%s" % (sid, a))
        z = find_committed(T6, "s%d_wall0" % sid)
        if None in (b, w, z):
            stop("site %d has no readable committed single-site arms; READING 4's excess cannot\n"
                 "   be formed." % sid)
        return (b - w) + z
    d20 = [site_ds(s, ("F2", "B", "A")) for s in S20_F2_SITES]
    d10 = [site_ds(s, ("A", "B")) for s in N10_A_SITES]
    # TWO DIFFERENT QUANTITIES, and they are not interchangeable:
    #   signed sum  |sum_i d_i|   is what DIAGNOSTIC 16 READING 3's excess subtracts
    #   Dsum        sum_i |d_i|   is what DIAGNOSTIC 15 READING 1 compares Dset against
    # At n = 10 arm A they are 0.2499 and 0.3521. Using one for the other reproduces the
    # verdict and not the margin, which is why the null experiment checks both.
    s20, s10 = sum(d20), sum(d10)
    dsum10 = sum(abs(x) for x in d10)
    e20 = abs(dsets["S6-20_D16F2"]) - abs(s20)
    e10 = abs(dsets["S6-10_D15"]) - abs(s10)
    return e20, e10, s20, s10, dsum10


def recompute_verdicts(T6, dsets):
    """READING 4. Each committed verdict recomputed from the Dsets this parameter set gives.
    No tolerance is new: all four carry 0.0500."""
    out = {}
    # D16 READING 1: monotone over n = 3, 5, 10, 20 at the matched field 1.3000 e
    seq = [dsets["S6-3_D16"], dsets["S6-5_D16"], dsets["S6-10_D16"], dsets["S6-20_D16F1"]]
    a = [abs(x) for x in seq]
    steps = [a[i] - a[i + 1] for i in range(3)]
    tol = max(0.10 * max(a), SHIELD_ABS)
    if all(s > 0 for s in steps):
        v1 = "FALLS WITH n"
    elif all(s < 0 for s in steps):
        v1 = "GROWS WITH n"
    else:
        v1 = "NOT MONOTONIC"
    out["D16 READING 1"] = (v1, "steps %s, tol %.4f, thinnest %.4f"
                            % (", ".join("%+.4f" % s for s in steps), tol, min(abs(s) for s in steps)))
    # D16 READING 2
    d2 = abs(dsets["S6-20_D16F2"]) - abs(dsets["S6-10_D15"])
    out["D16 READING 2"] = ("CONSISTENT" if abs(d2) > TOL_VERDICT else "NOT DISTINGUISHED",
                            "difference %+.4f against tol %.4f, margin %.4f"
                            % (d2, TOL_VERDICT, abs(abs(d2) - TOL_VERDICT)))
    # D16 READING 3
    e20, e10, s20, s10, dsum10 = excesses(T6, dsets)
    d3 = e20 - e10
    if abs(d3) <= TOL_VERDICT:
        v3 = "FLAT"
    elif d3 < 0:
        v3 = "COOPERATIVITY FALLS"
    else:
        v3 = "COOPERATIVITY GROWS"
    out["D16 READING 3"] = (v3, "excess %+.6f and %+.6f, difference %+.6f against tol %.4f, "
                                "margin %.4f" % (e20, e10, d3, TOL_VERDICT, abs(abs(d3) - TOL_VERDICT)))
    # D15 READING 1 -- against Dsum = sum of magnitudes, not the signed sum
    dsum = dsum10
    d4 = abs(dsets["S6-10_D15"]) - dsum
    if abs(dsets["S6-10_D15"]) <= TOL_VERDICT and dsum <= TOL_VERDICT:
        v4 = "NO EFFECT TO BOUND"
    elif d4 > TOL_VERDICT:
        v4 = "NOT BOUNDED"
    else:
        v4 = "BOUNDED"
    out["D15 READING 1"] = (v4, "Dset %.4f - Dsum %.4f = %+.4f against tol %.4f, margin %.4f"
                            % (abs(dsets["S6-10_D15"]), dsum, d4, TOL_VERDICT,
                               abs(abs(d4) - TOL_VERDICT)))
    return out


def selftest(root):
    """THE NULL EXPERIMENT. Feed the arithmetic the committed Ne arms as if they were the new
    ones. Every difference must be exactly zero, f must be withheld everywhere, and all four
    verdicts must come back as the committed ones. A sign error, a swapped arm, a mixed control
    or a wrong condition mapping all fail this; none of them would be visible in a run that
    merely completed."""
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    TRX.controls(T6)                                   # D6's control table must parse
    NEG.check_correction(T6)                           # A1.2 must reproduce D14's numbers
    print("=" * 74)
    print("SELF TEST  the null experiment: the committed Ne arms read as if new")
    print("=" * 74)
    bad = 0
    dsets, ne_tab = {}, {}
    for nm, n, f, w, w0 in COND:
        wt, w0t = w.rsplit("/", 1), w0.rsplit("/", 1)
        E = {}
        for arm, src in (("wall", w), ("wall0", w0)):
            tree = src.split("/")[0]
            stem = src.split("/")[2]
            for st in ("R", "TS"):
                e = fsp(T6 / tree / "def2-SVPD" / ("%s_%s" % (stem, st)))
                if e is None:
                    stop("no committed energy for %s_%s" % (src, st))
                E[(arm, st)] = e
        bare = find_committed(T6, BARE[nm].split("/")[2])
        if bare is None:
            stop("no committed bare arm for %s" % nm)
        dz, draw, dd, fq = deltas(E, E, bare)
        z_q, d_q = NE_COMMITTED[nm]
        # the committed Ne values, re-derived, must match what READING 2 prints
        tree = w0.split("/")[0]
        nez = committed_ddE(T6, tree + "/def2-SVPD", w0.split("/")[2])
        new_wall = committed_ddE(T6, w.split("/")[0] + "/def2-SVPD", w.split("/")[2])
        ned = bare - new_wall + nez
        for lbl, got, want in (("ddE(wall0)", nez, z_q), ("d_shield", ned, d_q)):
            if abs(got - want) > 5e-5:
                print("   FAIL %-13s %s re-derives to %+.6f, READING 2 prints %+.4f"
                      % (nm, lbl, got, want)); bad += 1
        for lbl, got in (("dz", dz), ("d(wall,raw)", draw), ("dd", dd)):
            if got != 0.0:
                print("   FAIL %-13s %s is %+.3e in the null experiment, must be exactly 0"
                      % (nm, lbl, got)); bad += 1
        if fq is not None:
            print("   FAIL %-13s f was reported at dz = 0; it must be withheld" % nm); bad += 1
        dsets[nm] = ned
        ne_tab[nm] = (nez, ned)
    print("   %d conditions: every dz, d(wall,raw) and dd exactly zero, f withheld in all"
          % len(COND))
    v = recompute_verdicts(T6, dsets)
    for k in sorted(VERDICTS_COMMITTED):
        want, margin_want, _ = VERDICTS_COMMITTED[k]
        got, detail = v[k]
        if got != want:
            print("   FAIL %-14s recomputes to %r, committed is %r" % (k, got, want)); bad += 1
            continue
        mm = re.search(r"margin ([\d.]+)|thinnest ([\d.]+)", detail)
        margin_got = float(mm.group(1) or mm.group(2)) if mm else None
        if margin_got is None or abs(margin_got - margin_want) > 5e-5:
            print("   FAIL %-14s %s but its margin is %s, the committed margin is %.4f\n"
                  "                  %s" % (k, got, margin_got, margin_want, detail)); bad += 1
        else:
            print("   ok   %-14s %-22s margin %.4f   %s" % (k, got, margin_got, detail))
    print()
    if bad:
        stop("the null experiment failed in %d place(s). The arithmetic is wrong and no real\n"
             "   run may be read with it." % bad)
    print("SELF TEST PASSED - the arithmetic reproduces the committed record exactly and")
    print("returns zero for a null change.")
    return 0




def analyse(root, runs, outdir):
    """D23.3's gates and D23.4's readings, in the order the pre-registration fixes them."""
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    R = Path(runs)
    O = []

    def say(s=""):
        O.append(s)
        print(s)
        sys.stdout.flush()

    for rel in (PREREG, SICHK):
        if not (T6 / rel).is_file():
            stop("%s is not committed; nothing is read before it exists" % rel)
    C = TRX.controls(T6)                 # D6's control medians, for D6.3
    NEG.check_correction(T6)             # A1.2 must reproduce D14's numbers before use

    say("Stage 6, DIAGNOSTIC 23 RESULT.")
    say("Read against %s. Fifty single points at %s." % (PREREG, BASIS))
    say("")

    # ---------------------------------------------------------------- GATE 4
    dirs = sorted(d for d in R.iterdir() if d.is_dir())
    if len(dirs) != 50:
        stop("%d run directories at %s, expected the 50 D23.2 specifies" % (len(dirs), R))
    nonorm = [d.name for d in dirs
              if "****ORCA TERMINATED NORMALLY****" not in (d / "job.out").read_text(errors="replace")]
    if nonorm:
        stop("%d run(s) did not terminate normally: %s\n"
             "   A failed job is not committed; a note about it is. D23.3 GATE 4."
             % (len(nonorm), ", ".join(nonorm[:6])))
    say("GATE 4   all 50 runs terminated normally")

    # ---------------------------------------------------------------- GATE 0
    man = T6 / "d23" / "INPUT_CHECKSUMS.sha256"
    if not man.is_file():
        stop("%s is not committed; GATE 0 has nothing to check against" % man)
    want = {}
    for l in man.read_text().splitlines():
        if l.strip():
            h, rel = l.split(None, 1)
            want[rel.strip()] = h
    nck = 0
    for d in dirs:
        for f in sorted(d.iterdir()):
            key = "%s/%s" % (d.name, f.name)
            if key not in want:
                continue
            if D22.sha(f) != want[key]:
                stop("%s is not the committed input that was submitted.\n"
                     "   The outputs were produced by something other than what is in the\n"
                     "   repository. D23.3 GATE 0." % key)
            nck += 1
    if nck != len(want) - 1:                       # the manifest also lists run.pbs
        stop("GATE 0 checked %d of the manifest's %d input files; the run tree is incomplete.\n"
             "   D23.3 GATE 0." % (nck, len(want) - 1))
    say("GATE 0   all %d committed inputs byte-identical to what ran" % nck)

    # ---------------------------------------------------------------- GATE 2
    cE, ctl_nb, ctl_nel = {}, None, None
    for st in ("R", "TS"):
        new = fsp(R / ("ctl_%s" % st))
        old = fsp(T6 / "d16/def2-SVPD" / ("ctl_%s" % st))
        if None in (new, old):
            stop("the control energy is unreadable at %s" % st)
        if abs(new - old) > GATE_EH:
            stop("the re-run control at %s differs from the committed control by %.3e Eh,\n"
                 "   over the %.0e Eh D23.3 GATE 2 allows. Nothing new may be compared with\n"
                 "   anything committed. D23.3 GATE 2." % (st, abs(new - old), GATE_EH))
        cE[st] = new
        t = (R / ("ctl_%s" % st) / "job.out").read_text(errors="replace")
        nb = int(re.search(r"Number of basis functions\s+\.+\s+(\d+)", t).group(1))
        nel = int(re.search(r"Number of Electrons\s+NEL\s+\.+\s+(\d+)", t).group(1))
        ctl_nb = nb if ctl_nb is None else ctl_nb
        ctl_nel = nel if ctl_nel is None else ctl_nel
        if nb != ctl_nb or nel != ctl_nel:
            stop("the two control states disagree on nbasis or NEL")
    say("GATE 2   the control reproduces the committed control to %.1e Eh at both states; "
        "nbasis %d, NEL %d" % (GATE_EH, ctl_nb, ctl_nel))

    # ---------------------------------------------------------------- GATE 3
    nq = 0
    for d in dirs:
        if d.name.startswith("ctl_"):
            continue
        inp = (d / "job.inp").read_text()
        mem = mem_of(inp)
        t = (d / "job.out").read_text(errors="replace")
        csum = gate_centres_are_Q(d.name, t, len(mem), ctl_nb, ctl_nel)
        if csum is not None:
            exp = -2 + sum(q for _, q, _ in mem)
            if abs(csum - exp) > 1e-4:
                stop("%s prints a charge sum of %.6f; -2 plus its centre charges is %.6f.\n"
                     "   The centres are contributing something other than their point charge.\n"
                     "   D23.3 GATE 3." % (d.name, csum, exp))
        nq += len(mem)
    say("GATE 3   every walled output: nbasis %d, NEL %d, all %d centres labelled %r, charge "
        "sums exact" % (ctl_nb, ctl_nel, nq, AUD.ECP_SYMBOL))
    say("")

    # ------------------------------------------------------------- READING 1
    say("READING 1  CLASSIFICATION, BEFORE ANY ENERGY IS READ")
    say("")
    cls, notread = {}, {}
    for key in sorted(SETS):
        for src, arm, serves in arms():
            tree, stem = src.split("/")[0], src.split("/")[2]
            for st in Q13.GEOM:
                run = R / ("%s_%s_%s_%s" % (key, tree, stem, st))
                t = (run / "job.out").read_text(errors="replace")
                d = AUD.read(t)
                if d is None:
                    stop("%s is not a completed single point on this dianion" % run.name)
                mem = mem_of((run / "job.inp").read_text())
                charged = any(abs(q) > 0 for _, q, _ in mem)
                if charged:                       # A1.2 applies to charged arms only
                    d = dict(d, dip=NEG.corrected_dip(run, mem, NEG.com(T6, st)))
                v, why = TRX.verdict(d, C)
                cls[(key, src, st)] = (v, why, charged)
        for nm, n, f, w, w0 in COND:
            armcls = [cls[(key, s, st)][0] for s in (w, w0) for st in Q13.GEOM]
            if not all(x in ("CLEAN", "POLARISED") for x in armcls):
                notread.setdefault(key, []).append(nm)
        say("   %-3s %s" % (SETS[key]["el"],
            "every arm CLEAN or POLARISED, all 7 conditions readable" if key not in notread
            else "NOT READABLE: " + ", ".join(notread[key])))
        for (k2, s, st), (v, why, ch) in sorted(cls.items()):
            if k2 == key and v not in ("CLEAN", "POLARISED"):
                say("        %-34s %-3s %-10s %s" % (s, st, v, "; ".join(why)))
    say("")

    # --------------------------------------------------------- READING 2 / 3
    say("READING 2 and 3  THE TWO DIFFERENCES, AND THE PROPAGATION FACTOR")
    say("")
    say("   condition        set    dz        d(wall,raw)  dd        f         |dd| vs 0.05")
    res = {}
    for key in sorted(SETS):
        for nm, n, f_, w, w0 in COND:
            if nm in notread.get(key, []):
                say("   %-16s %-4s   NOT READABLE for this set" % (nm, SETS[key]["el"]))
                continue
            E, NE = {}, {}
            for arm, src in (("wall", w), ("wall0", w0)):
                tree, stem = src.split("/")[0], src.split("/")[2]
                for st in Q13.GEOM:
                    E[(arm, st)] = fsp(R / ("%s_%s_%s_%s" % (key, tree, stem, st)))
                    NE[(arm, st)] = fsp(T6 / tree / "def2-SVPD" / ("%s_%s" % (stem, st)))
                    if None in (E[(arm, st)], NE[(arm, st)]):
                        stop("missing energy for %s %s %s %s" % (key, nm, arm, st))
            bare = find_committed(T6, BARE[nm].split("/")[2])
            dz, draw, dd, fq = deltas(E, NE, bare)
            # the same two numbers formed the way READING 2 words them, each ddE against its
            # own control, so the control's contribution is visible rather than assumed
            newz = ddE_pair(E[("wall0", "R")], E[("wall0", "TS")], cE["R"], cE["TS"])
            neww = ddE_pair(E[("wall", "R")], E[("wall", "TS")], cE["R"], cE["TS"])
            nez, ned = NE_COMMITTED[nm]
            dz2, dd2 = newz - nez, (bare - neww + newz) - ned
            if abs(dz - dz2) > 1e-3 or abs(dd - dd2) > 1e-3:
                stop("for %s at %s the two routes to dz/dd differ by %.3e / %.3e kcal/mol,\n"
                     "   more than the control discrepancy GATE 2 bounds. READING 2."
                     % (nm, key, abs(dz - dz2), abs(dd - dd2)))
            res[(key, nm)] = (dz, draw, dd, fq, bare - neww + newz)
            say("   %-16s %-4s %+9.4f %+11.4f %+9.4f  %-9s %s"
                % (nm, SETS[key]["el"], dz, draw, dd,
                   "%+.3f" % fq if fq is not None else "withheld",
                   "over" if abs(dd) > SHIELD_ABS else "within"))
    say("")
    say("   f is withheld where |dz| <= %.3f kcal/mol, as READING 3 fixes." % DZ_FLOOR)
    say("   Both routes to each difference agreed to within the GATE 2 control bound.")
    say("")

    # ------------------------------------------------------------- READING 4
    say("READING 4  THE FOUR COMMITTED VERDICTS, RECOMPUTED")
    say("")
    changed = {}
    for key in sorted(SETS):
        if key in notread:
            say("   %-3s not evaluated: %s is NOT READABLE"
                % (SETS[key]["el"], ", ".join(notread[key])))
            continue
        dsets = {nm: res[(key, nm)][4] for nm, n, f_, w, w0 in COND}
        v = recompute_verdicts(T6, dsets)
        for k in sorted(VERDICTS_COMMITTED):
            want = VERDICTS_COMMITTED[k][0]
            got, detail = v[k]
            if got != want:
                changed.setdefault(key, []).append((k, want, got))
            say("   %-3s %-14s committed %-22s now %-22s %s"
                % (SETS[key]["el"], k, want, got, "CHANGED" if got != want else ""))
            say("       %s" % detail)
    say("")

    # ------------------------------------------------------------- READING 5
    say("READING 5  THE HE FACTORS AGAINST THE BASIS ANALOGUE (reported, no criterion)")
    say("")
    for arm, (dz_a, dd_a) in sorted(ANALOGUE.items()):
        say("   %-6s basis change: dz %+.4f, dd %+.4f, absorbed %.0f%%"
            % (arm, dz_a, dd_a, 100 * (1 - abs(dd_a / dz_a))))
    hef = [res[(k, nm)][3] for k, nm in res if k == "he" and res[(k, nm)][3] is not None]
    if hef:
        say("   He factors: %s" % ", ".join("%+.3f" % x for x in sorted(hef)))
        say("   mean %+.3f, which corresponds to %.0f%% absorbed by the charged arm"
            % (sum(hef) / len(hef), 100 * (1 - abs(sum(hef) / len(hef)))))
    else:
        say("   no He factor was reported; every |dz| sat at or below the floor")
    say("")

    # ------------------------------------------------------------- CRITERIA
    say("D23.5  CRITERIA")
    say("")
    def criterion(keys, label):
        dds = [abs(res[(k, nm)][2]) for k in keys for nm, n, f_, w, w0 in COND
               if (k, nm) in res]
        if not dds:
            return "INSUFFICIENT READABLE CONDITIONS"
        readable = {k: 7 - len(notread.get(k, [])) for k in keys}
        if any(readable[k] < 4 for k in keys):
            return ("INSUFFICIENT READABLE CONDITIONS (%s)"
                    % ", ".join("%s %d of 7" % (SETS[k]["el"], readable[k]) for k in keys))
        mx = max(dds)
        ch = [c for k in keys for c in changed.get(k, [])]
        if mx <= SHIELD_ABS:
            return "PARAMETER CHOICE DOES NOT REACH d_shield (max |dd| %.4f <= %.2f)" % (mx, SHIELD_ABS)
        if not ch:
            return "VALUES PROVISIONAL ON THE PARAMETER SET, VERDICTS STAND (max |dd| %.4f)" % mx
        return ("A PROPERTY OF THE PARAMETER SET: %s (max |dd| %.4f)"
                % ("; ".join("%s %s -> %s" % c for c in ch), mx))
    he_only = criterion(["he"], "He alone")
    both = criterion(sorted(SETS), "both sets")
    say("   on He alone  : %s" % he_only)
    say("   on both sets : %s" % both)
    say("")
    if he_only != both:
        say("   THE TWO DISAGREE. C23.5 quotes the He-only reading, because EXPERIMENT 0")
        say("   validated the He block against its own published fitting target to 1.381e-06 Eh")
        say("   and could not validate the Ar block better than 3.067e-03 Eh. The Ar result is")
        say("   a sensitivity probe, not a comparison of published operators.")
        say("   QUOTED: %s" % he_only)
    else:
        say("   The two readings agree, so no C23.5 preference arises.")
    say("")
    say("VERDICTS: %s" % he_only)
    say("")

    out = Path(outdir) / "DIAGNOSTIC23_RESULT.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(O) + "\n")
    print("written: %s" % out)
    return 0


def main(argv):
    if len(argv) == 5 and argv[1] == "analyse":
        return analyse(argv[2], argv[3], argv[4])
    if len(argv) == 5 and argv[1] == "prepare":
        return prepare(argv[2], argv[3], argv[4])
    if len(argv) == 3 and argv[1] == "selftest":
        return selftest(argv[2])
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
