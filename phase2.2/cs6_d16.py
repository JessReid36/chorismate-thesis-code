#!/usr/bin/env python3
"""
cs6_d16.py - STAGE6_DIAGNOSTIC16_20261029.txt.

Does the overpolarisation artefact grow with the NUMBER of charges, at a fixed total field?

  python3 cs6_d16.py prepare <results_repo> <OUT> <HPC_DIR>
  python3 cs6_d16.py analyse <results_repo> <RUNS> <OUT_DIR>

NOTHING ALREADY COMMITTED IS RETYPED OR RE-RUN.
  cs6_additivity  write_run, members, rname, one0, PBS, d7_table
  cs6_negatives   zeroed, com, qdip, corrected_dip, check_correction   (A1.2's correction)
  cs6_d15         sites()                                              (the sign split)
  cs6_transfer    design, head, controls, verdict, one, QCLEAN
  cs13_qscale     ECP_TAIL, FRAME_SHA, GEOM, read_xyz_raw, descriptors, SHIELD_REL, SHIELD_ABS
"""
import hashlib, re, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cs13_qscale as Q13
import cs6_transfer as TR
import cs6_additivity as AD
import cs6_negatives as NG
import cs6_d15 as D15

PREREG = "STAGE6_DIAGNOSTIC16_20261029.txt"
DESIGNS = (("S6-2", 2), ("S6-3", 3), ("S6-5", 5), ("S6-10", 10), ("S6-20", 20))
BIG, BIGN = "S6-20", 20
BASIS = "def2-SVPD"
F2 = 2.1900                      # D15's own condition; S6-10's committed Dset 0.4140 sits here
CTL_EH = 1e-6
D15_DSET10, D15_SIGNED10 = 0.4140, 0.2500      # DIAGNOSTIC 15's committed numbers at F2


def stop(m):
    print("STOP: " + m); sys.stdout.flush(); sys.exit(1)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def setup(T6):
    """The five designs, and the matched field D16.2 derives rather than chooses."""
    D = {nm: TR.design(T6, nm, n) for nm, n in DESIGNS}
    F1 = min(TR.QCLEAN * d["sabs"] / d["maxq"] for d in D.values())
    if abs(F1 - 1.30) > 5e-4:
        stop("the largest admissible matched field is %.4f e, not the 1.3000 D16.2 states" % F1)
    lam = {("F1", nm): F1 / D[nm]["sabs"] for nm, _ in DESIGNS}
    lam[("F2", BIG)] = F2 / D[BIG]["sabs"]
    for (f, nm), l in lam.items():
        if l > 1 + 1e-12 or l * D[nm]["maxq"] > TR.QCLEAN + 1e-9:
            stop("%s at %s needs lambda %.6f, largest |q| %.4f - outside what D16.2 admits"
                 % (nm, f, l, l * D[nm]["maxq"]))
    return D, F1, lam


def reusable_wall0(T6, D):
    """The wall-at-zero runs this diagnostic reuses. That arm does not depend on lambda, so a site
    that already has one keeps it - and how many there are decides whether the run count is the
    148 D16.2 states. Checked BEFORE the committer's first commit, because a rehearsal found that
    a missing one only surfaced as a wrong total after the code commit had landed."""
    have = {}
    for s in D[BIG]["site"]:
        for tree in ("additivity", "negatives"):
            if (T6 / tree / BASIS / ("s%d_wall0_R" % s)).is_dir():
                have[s] = tree; break
    want = sorted(set(D[BIG]["site"]) & set(D["S6-10"]["site"]))
    missing = [s for s in want if s not in have]
    if missing:
        stop("the committed wall-at-zero run for site%s %s is missing. D16.2 reuses %d of them "
             "and writes %d new; without these the run count is not the 148 it states."
             % ("" if len(missing) == 1 else "s", ", ".join(str(x) for x in missing),
                len(want), BIGN - len(want)))
    if len(have) != len(want):
        stop("%d of %s's sites have a committed wall-at-zero run; D16.2 expects %d"
             % (len(have), BIG, len(want)))
    return have


def check_existing(T6):
    D, F1, lam = setup(T6)
    have = reusable_wall0(T6, D)
    print("   %d of %s's sites reuse a committed wall-at-zero run (%s); %d are new"
          % (len(have), BIG, ", ".join(str(s) for s in sorted(have)), BIGN - len(have)))
    return D, F1, lam, have


def rn(nm, arm, field, st):
    return "%s_%s_%s" % (nm, arm, st) if arm == "wall0" else "%s_%s_%s_%s" % (nm, arm, field, st)


# ---------------------------------------------------------------- prepare

def prepare(root, out, hpc):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    out = Path(out)
    if out.exists():
        stop("%s already exists - this never overwrites" % out)
    if not (T6 / PREREG).is_file():
        stop("%s is not committed. The pre-registration goes in before any input." % PREREG)
    print("   %s is committed (sha256 %s)" % (PREREG, sha(T6 / PREREG)[:16]))
    geo, raw = {}, {}
    for st, fn in (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz")):
        p = T6 / "frame_42436" / fn
        if sha(p) != Q13.FRAME_SHA[fn]:
            stop("%s does not match the sha256 cs13_qscale records" % fn)
        geo[st], raw[st] = p.read_text(), Q13.read_xyz_raw(p)
    print("   both geometries match the sha256 the committed scripts record")
    NG.check_correction(T6)
    D, F1, lam = setup(T6)
    print("   matched field F1 = %.4f e, derived from the committed arrangements (S6-2 binds it)" % F1)
    for nm, n in DESIGNS:
        l = lam[("F1", nm)]
        print("      %-7s n %-3d lambda %.6f   largest |q| %.4f e" % (nm, n, l, l * D[nm]["maxq"]))
    print("   F2 = %.4f e: %s at lambda %.6f, against D15's committed %.4f for S6-10"
          % (F2, BIG, lam[("F2", BIG)], D15_DSET10))

    have = reusable_wall0(T6, D)
    new_w0 = [s for s in D[BIG]["site"] if s not in have]
    print("   %d of %s's sites reuse a committed wall-at-zero run; %d are new"
          % (len(have), BIG, len(new_w0)))

    out.mkdir(parents=True)
    plan, rs = [], []
    def emit(nm, arm, field, mem, walled):
        for st in Q13.GEOM:
            d = out / BASIS / rn(nm, arm, field, st)
            AD.write_run(d, BASIS, geo[st], raw[st], mem, arm)
            rs.append(rn(nm, arm, field, st))
            plan.append((nm, arm, field or "-", st, len(mem),
                         sum(abs(q) for _, q, _ in mem) if arm != "wall0" else 0.0))
    for nm, n in DESIGNS:                                   # five designs at F1
        m = AD.members(D[nm], D[nm]["site"], lam[("F1", nm)])
        emit(nm, "wall0", None, NG.zeroed(m), True)
        emit(nm, "bare", "F1", m, False)
        emit(nm, "wall", "F1", m, True)
    m2 = AD.members(D[BIG], D[BIG]["site"], lam[("F2", BIG)])   # S6-20 at F2; wall0 shared
    emit(BIG, "bare", "F2", m2, False)
    emit(BIG, "wall", "F2", m2, True)
    for s in D[BIG]["site"]:                                # the twenty singles at F2
        one = AD.members(D[BIG], [s], lam[("F2", BIG)])
        if s in new_w0:
            emit("s%d" % s, "wall0", None, NG.zeroed(one), True)
        emit("s%d" % s, "bare", "F2", one, False)
        emit("s%d" % s, "wall", "F2", one, True)
    for st in Q13.GEOM:                                     # controls
        AD.write_run(out / BASIS / ("ctl_%s" % st), BASIS, geo[st], raw[st], [], "none")
        rs.append("ctl_%s" % st); plan.append(("ctl", "none", "-", st, 0, 0.0))
    (out / BASIS / "run.pbs").write_text(
        AD.PBS % dict(tag="16", hpc=hpc, basis=BASIS, n=len(rs),
                      runs="\n".join("run %s" % r for r in rs)))
    with open(out / "PLAN.tsv", "w") as fh:
        fh.write("# %s D16.2; cs6_d16.py prepare\n" % PREREG)
        fh.write("# F1 %.4f e (derived), F2 %.4f e; %s\n" % (F1, F2, BASIS))
        for nm, _ in DESIGNS:
            fh.write("# %-7s lambda(F1) %.6f\n" % (nm, lam[("F1", nm)]))
        fh.write("# %-7s lambda(F2) %.6f\n" % (BIG, lam[("F2", BIG)]))
        fh.write("what\tarm\tfield\tstate\tn\tsum_abs_q\n")
        for r in plan:
            fh.write("%s\t%s\t%s\t%s\t%d\t%.4f\n" % r)
    files = sorted(p for p in out.rglob("*") if p.is_file())
    (out / "INPUT_CHECKSUMS.sha256").write_text(
        "".join("%s  %s\n" % (sha(p), p.relative_to(out)) for p in files))
    print("   %d run directories, %d files including the manifest" % (len(plan), len(files) + 1))
    if len(plan) != 148:
        stop("%d runs, but D16.2 states 148" % len(plan))
    print("   NOTHING WAS SUBMITTED BY THIS BUNDLE.")


# ---------------------------------------------------------------- analyse

def analyse(root, runs, outdir):
    T6, NEW, outdir = Path(root) / "phase2.2" / "cs_stage6", Path(runs), Path(outdir)
    if not (T6 / PREREG).is_file():
        stop("%s is not committed" % PREREG)
    D, F1, lam = setup(T6)
    C = TR.controls(T6)
    GEO = {st: TR.read_xyz(T6 / "frame_42436" / fn) for st, fn in
           (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz"))}
    O = {st: NG.com(T6, st) for st in Q13.GEOM}
    GLBL = Q13.read_xyz(T6 / "frame_42436" / "reactant_qm.xyz")

    cand = [m for m in NEW.rglob("INPUT_CHECKSUMS.sha256") if ".git" not in m.parts]
    if len(cand) != 1:
        stop("expected one INPUT_CHECKSUMS.sha256 under %s, found %d" % (NEW, len(cand)))
    base, bad = cand[0].parent, []
    for line in cand[0].read_text().splitlines():
        want, rel = line.split("  ", 1)
        f = base / rel
        if not f.is_file():
            bad.append("missing %s" % rel)
        elif sha(f) != want:
            bad.append("changed %s" % rel)
    if bad:
        stop("GATE 0 failed, %d input files do not match the manifest:\n   %s"
             % (len(bad), "\n   ".join(bad[:8])))

    def where(nm, arm):
        """Which committed tree holds this run. Only reused wall-at-zero runs live elsewhere."""
        if arm == "wall0" and nm.startswith("s"):
            s = int(nm[1:])
            for tree in ("additivity", "negatives"):
                if (T6 / tree / BASIS / ("%s_wall0_R" % nm)).is_dir():
                    return T6 / tree, tree
        return base, "new"

    E, CL, WHY, CTL = {}, {}, {}, {}
    for tag, tr in (("new", base), ("additivity", T6 / "additivity"), ("negatives", T6 / "negatives")):
        if (tr / BASIS / "ctl_R").is_dir():
            for st in Q13.GEOM:
                e, nb, dd = AD.one0(tr / BASIS / ("ctl_%s" % st), GEO[st], 0, False)
                CTL[(tag, st)] = e
                CL[("ctl_" + tag, "none", "-", st)] = TR.verdict(dd, C)[0]
                WHY[("ctl_" + tag, "none", "-", st)] = TR.verdict(dd, C)[1]

    def take(nm, arm, field, mem, walled):
        tr, tag = where(nm, arm)
        for st in Q13.GEOM:
            d = tr / BASIS / rn(nm, arm, field, st)
            if arm == "bare":
                e, nb, dd = TR.one(d, GEO[st], len(mem))
            else:
                e, nb, dd = AD.one0(d, GEO[st], len(mem), walled)
            if walled:
                dd = dict(dd); dd["dip"] = NG.corrected_dip(d, mem, O[st])
            E[(nm, arm, field or "-", st)] = (e, tag)
            CL[(nm, arm, field or "-", st)], WHY[(nm, arm, field or "-", st)] = TR.verdict(dd, C)

    for nm, _ in DESIGNS:
        m = AD.members(D[nm], D[nm]["site"], lam[("F1", nm)])
        take(nm, "wall0", None, NG.zeroed(m), True)
        take(nm, "bare", "F1", m, False)
        take(nm, "wall", "F1", m, True)
    m2 = AD.members(D[BIG], D[BIG]["site"], lam[("F2", BIG)])
    take(BIG, "bare", "F2", m2, False)
    take(BIG, "wall", "F2", m2, True)
    for s in D[BIG]["site"]:
        one = AD.members(D[BIG], [s], lam[("F2", BIG)])
        take("s%d" % s, "wall0", None, NG.zeroed(one), True)
        take("s%d" % s, "bare", "F2", one, False)
        take("s%d" % s, "wall", "F2", one, True)

    H = TR.H
    def ddE(nm, arm, field):
        d = {}
        for g in Q13.GEOM:
            e, tag = E[(nm, arm, field, g)]
            c = CTL[(tag, g)]
            z = 0.0
            if arm == "wall":
                e0, t0 = E[(nm, "wall0", "-", g)]
                z = e0 - CTL[(t0, g)]
            d[g] = e - c - z
        return (d["TS"] - d["R"]) * H
    def dsh(nm, field):
        return ddE(nm, "bare", field) - ddE(nm, "wall", field)
    def well(nm, field):
        ks = [(nm, "wall0", "-", g) for g in Q13.GEOM]
        for arm in ("bare", "wall"):
            ks += [(nm, arm, field, g) for g in Q13.GEOM]
        return all(CL[k] in ("CLEAN", "POLARISED") for k in ks)

    O_ = ["STAGE 6 - DIAGNOSTIC 16 RESULT: does the artefact grow with the number of charges?", "",
          "Measured by cs6_d16.py against %s, committed before any of this" % PREREG,
          "diagnostic's runs existed and applied unchanged. Classification is D6.3 with AMENDMENT 1",
          "A1.2's dipole correction. No geometry moved. In vacuo. %s throughout." % BASIS,
          "",
          "D16.4 states the largest limitation first and it is restated here: the five designs are",
          "the optimiser's best arrangement at each n, on different sites, not one arrangement with",
          "charges added. n and the arrangement vary together and this does not separate them.",
          "",
          "   matched field F1 = %.4f e, derived from the committed arrangements; F2 = %.4f e" % (F1, F2), ""]

    O_ += ["GATE 1: EVERY CONTROL MUST AGREE WITH DIAGNOSTIC 14'S, TO %g Eh" % CTL_EH, ""]
    g1 = True
    for tag in ("new", "additivity", "negatives"):
        for st in Q13.GEOM:
            if (tag, st) not in CTL:
                continue
            ok = abs(CTL[(tag, st)] - CTL[("additivity", st)]) <= CTL_EH
            g1 &= ok
            O_.append("   %-12s %-2s %18.10f   D14 %18.10f   %s"
                      % (tag, st, CTL[(tag, st)], CTL[("additivity", st)], "pass" if ok else "FAIL"))
    O_.append("")
    if not g1:
        return finish(O_, outdir, "GATE 1 failed: the trees disagree on the control. Nothing is read.")

    O_ += ["EVERY RUN CLASSIFIED BY D6.3 WITH A1.2'S CORRECTION, READ BEFORE ANY NUMBER", "",
           "   run                              class       note"]
    for k in sorted(CL, key=lambda k: (k[0], k[2], k[1], k[3])):
        O_.append("   %-32s %-10s  %s" % ("_".join(x for x in k if x != "-"), CL[k], "; ".join(WHY[k])))
    O_.append("")

    # READING 1
    O_ += ["READING 1: DOES Dset GROW WITH n, AT THE MATCHED FIELD %.4f e?" % F1, "",
           "   design   n    lambda     max|q|    Dset       class"]
    pts, unread = [], []
    for nm, n in DESIGNS:
        l = lam[("F1", nm)]
        if not well(nm, "F1"):
            O_.append("   %-7s %-4d %.6f  %.4f    NOT READABLE" % (nm, n, l, l * D[nm]["maxq"]))
            unread.append(nm); continue
        v = abs(dsh(nm, "F1"))
        pts.append((n, nm, v))
        O_.append("   %-7s %-4d %.6f  %.4f   %8.4f   %s"
                  % (nm, n, l, l * D[nm]["maxq"], v, CL[(nm, "wall", "F1", "R")]))
    O_.append("")
    if len(pts) < 3:
        v1 = "NOT READABLE (%d of 5 designs readable)" % len(pts)
        tol = None
    else:
        tol = max(Q13.SHIELD_REL * max(v for _, _, v in pts), Q13.SHIELD_ABS)
        vals = [v for _, _, v in sorted(pts)]
        spread = max(vals) - min(vals)
        O_ += ["   tol = max(0.10 x largest Dset, %.2f) = %.4f;  spread %.4f"
               % (Q13.SHIELD_ABS, tol, spread), ""]
        if max(vals) <= tol:
            v1 = "NO EFFECT TO COMPARE (every Dset <= tol)"
        elif spread <= tol:
            v1 = "FLAT"
        elif all(b >= a for a, b in zip(vals, vals[1:])):
            v1 = "GROWS WITH n"
        elif all(b <= a for a, b in zip(vals, vals[1:])):
            v1 = "FALLS WITH n"
        else:
            v1 = "NOT MONOTONIC"
    O_ += ["   READING 1: %s" % v1, ""]

    # READING 2
    O_ += ["READING 2: DOES IT HOLD AT THE STRONGER FIELD %.4f e?" % F2, ""]
    if not well(BIG, "F2"):
        v2 = "NOT READABLE"
        O_.append("   %s at F2 is unreadable" % BIG)
    else:
        d20 = abs(dsh(BIG, "F2"))
        t2 = tol if tol is not None else max(Q13.SHIELD_REL * max(d20, D15_DSET10), Q13.SHIELD_ABS)
        O_ += ["   %-7s n 20  Dset %8.4f      (lambda %.6f)" % (BIG, d20, lam[("F2", BIG)]),
               "   S6-10   n 10  Dset %8.4f      DIAGNOSTIC 15, committed" % D15_DSET10,
               "   difference %+.4f against tol %.4f" % (d20 - D15_DSET10, t2), ""]
        if abs(d20 - D15_DSET10) <= t2:
            v2 = "WITHIN TOLERANCE (the field does not decide it)"
        elif "GROWS" in v1:
            v2 = "CONSISTENT" if d20 > D15_DSET10 else "INCONSISTENT"
        elif "FALLS" in v1:
            v2 = "CONSISTENT" if d20 < D15_DSET10 else "INCONSISTENT"
        else:
            v2 = "REPORTED ONLY (Reading 1 gave no direction to be consistent with)"
    O_ += ["   READING 2: %s" % v2, ""]

    # READING 3
    O_ += ["READING 3: DOES THE COOPERATIVITY GROW WITH n?", ""]
    bad20 = [s for s in D[BIG]["site"] if not well("s%d" % s, "F2")]
    if bad20 or not well(BIG, "F2"):
        v3 = "NOT READABLE (%s)" % (("sites " + ", ".join(str(s) for s in bad20)) if bad20 else BIG)
        O_.append("   " + v3)
    else:
        S = {s: dsh("s%d" % s, "F2") for s in D[BIG]["site"]}
        O_.append("   site    |q| (e)   d_near    d_shield")
        for s in sorted(S, key=lambda s: -abs(S[s])):
            dn = Q13.descriptors(D[BIG]["X"][list(D[BIG]["site"]).index(s)], GLBL)["dn"]
            O_.append("   %-6d  %7.4f  %7.3f  %+9.4f"
                      % (s, abs(lam[("F2", BIG)] * D[BIG]["q"][list(D[BIG]["site"]).index(s)]), dn, S[s]))
        signed = abs(sum(S.values())); Dsum = sum(abs(v) for v in S.values())
        d20 = abs(dsh(BIG, "F2")); exc20 = d20 - signed
        exc10 = D15_DSET10 - D15_SIGNED10
        t3 = max(Q13.SHIELD_REL * max(Dsum, d20), Q13.SHIELD_ABS)
        O_ += ["",
               "   n = 20   Dsum %.4f   |signed sum| %.4f   Dset %.4f   excess %+.4f   ratio %.3f"
               % (Dsum, signed, d20, exc20, (d20 / signed if signed else float("nan"))),
               "   n = 10   (DIAGNOSTIC 15, committed)        Dset %.4f   excess %+.4f   ratio %.3f"
               % (D15_DSET10, exc10, D15_DSET10 / D15_SIGNED10),
               "   difference in excess %+.4f against tol %.4f" % (exc20 - exc10, t3), "",
               "   DIAGNOSTIC 14.3's rule applied to %s unchanged: %s"
               % (BIG, "BOUNDED" if d20 <= Dsum + t3 else "NOT BOUNDED"), ""]
        if exc20 > exc10 + t3:
            v3 = "COOPERATIVITY GROWS"
        elif exc20 < exc10 - t3:
            v3 = "COOPERATIVITY FALLS"
        else:
            v3 = "FLAT"
    O_ += ["   READING 3: %s" % v3, "",
           "d_near is reported and nothing is concluded from it (ERRATUM3 E3.1).", ""]
    return finish(O_, outdir, None)


def finish(O, outdir, tail):
    if tail:
        O += [tail, ""]
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "DIAGNOSTIC16_RESULT.txt").write_text("\n".join(O) + "\n")
    print("\n".join(O))
    print("\n   wrote %s" % (outdir / "DIAGNOSTIC16_RESULT.txt"))


if __name__ == "__main__":
    A_ = sys.argv[1:]
    if len(A_) == 4 and A_[0] == "prepare":
        prepare(A_[1], A_[2], A_[3])
    elif len(A_) == 4 and A_[0] == "analyse":
        analyse(A_[1], A_[2], A_[3])
    else:
        stop("usage: cs6_d16.py prepare <results_repo> <OUT> <HPC_DIR>\n"
             "       cs6_d16.py analyse <results_repo> <RUNS> <OUT_DIR>")
