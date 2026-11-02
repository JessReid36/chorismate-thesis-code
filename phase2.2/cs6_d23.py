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


def main(argv):
    if len(argv) == 5 and argv[1] == "prepare":
        return prepare(argv[2], argv[3], argv[4])
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
