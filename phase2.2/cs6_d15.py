#!/usr/bin/env python3
"""
cs6_d15.py - STAGE6_DIAGNOSTIC15_20261029.txt.

Does the sum of the ten per-site d_shield bound the d_shield of the whole mixed-sign design?

  python3 cs6_d15.py prepare <results_repo> <OUT> <HPC_DIR>
  python3 cs6_d15.py analyse <results_repo> <RUNS> <OUT_DIR>

NOTHING ALREADY COMMITTED IS RETYPED OR RE-RUN.
  cs6_additivity   write_run, members, rname, one0, the PBS template, ARM_LAM
  cs6_negatives    zeroed, com, qdip, corrected_dip, negatives  (A1.2's correction)
  cs6_transfer     design, lambdas, head, controls, verdict, one
  cs13_qscale      ECP_TAIL, FRAME_SHA, GEOM, read_xyz_raw, descriptors, SHIELD_REL, SHIELD_ABS

D15.2's conditions. The ten singles and the design come from three committed trees plus this one:
  additivity/     DIAGNOSTIC 14      the four positive singles and n10, both arms
  negatives/      AMENDMENT 1 A1.4   the six negative singles at arm B, and their wall-at-zero
                                     runs, which do not depend on lambda and so serve arm A too
  negatives_armA/ THIS DIAGNOSTIC    the six negatives bare and walled at arm A, def2-SVPD
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

PREREG = "STAGE6_DIAGNOSTIC15_20261029.txt"
DESIGN, NSPH = "S6-10", 10
DEST = "negatives_armA"
NEW_ARM, NEW_BASIS = "A", "def2-SVPD"
# D15.2: the two conditions that are read, and the one that is not
CONDS = (("def2-SVPD", "A"), ("def2-SVP", "B"))
UNREAD = ("def2-SVPD", "B", "n10_bare_B_R is STRAINED by D6.3 (gap 0.493 eV) in DIAGNOSTIC 14's "
                            "committed result; a gap is not A1.2's dipole artefact and no "
                            "correction reaches it")
CTL_EH, D7_TOL = 1e-6, 1e-4


def stop(m):
    print("STOP: " + m); sys.stdout.flush(); sys.exit(1)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def sites(D):
    """The ten, split by sign. Nothing is typed."""
    neg = NG.negatives(D)
    pos = [s for s in D["site"] if s not in neg]
    if len(pos) != 4:
        stop("expected four positive sites, parsed %d" % len(pos))
    return pos, neg


def check_existing(T6):
    """D15.2: everything this diagnostic does NOT re-run must already be committed. Separate from
    prepare() so the committer can run it BEFORE its first commit - a rehearsal found that a
    missing committed run stopped the generator only after the code commit had landed."""
    D = TR.design(T6, DESIGN, NSPH)
    pos, neg = sites(D)
    for rel, what in (("additivity/%s/n10_bare_%s_R" % (NEW_BASIS, NEW_ARM), "D14's n10 at arm A"),
                      ("additivity/%s/n10_wall0_R" % NEW_BASIS, "D14's n10 wall-at-zero")):
        if not (T6 / rel).is_dir():
            stop("%s is missing (%s)" % (rel, what))
    for s in pos:
        for a in ("bare", "wall"):
            if not (T6 / "additivity" / NEW_BASIS / ("s%d_%s_%s_R" % (s, a, NEW_ARM))).is_dir():
                stop("D14's positive single s%d (%s, arm A) is missing" % (s, a))
    for s in neg:
        if not (T6 / "negatives" / NEW_BASIS / ("s%d_wall0_R" % s)).is_dir():
            stop("A1.4's wall-at-zero run for s%d is missing; arm A reuses it" % s)
    print("   the committed runs D15 reads are all present; %d singles + n10 need no re-run"
          % (len(pos) + len(neg)))


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

    D = TR.design(T6, DESIGN, NSPH); REF = TR.design(T6, TR.REF_DESIGN, TR.REF_N)
    lam = TR.lambdas(D, REF)[AD.ARM_LAM[NEW_ARM]]
    pos, neg = sites(D)
    check_existing(T6)

    out.mkdir(parents=True)
    plan, rs = [], []
    i = {s: k for k, s in enumerate(D["site"])}
    for s in neg:
        nm, mem = "s%d" % s, AD.members(D, [s], lam)
        for st in Q13.GEOM:
            for arm in ("bare", "wall"):
                rn = "%s_%s" % (AD.rname(nm, arm, NEW_ARM), st)
                AD.write_run(out / NEW_BASIS / rn, NEW_BASIS, geo[st], raw[st], mem, arm)
                rs.append(rn); plan.append((NEW_BASIS, nm, arm, st, abs(mem[0][1])))
    for st in Q13.GEOM:
        rn = "ctl_%s" % st
        AD.write_run(out / NEW_BASIS / rn, NEW_BASIS, geo[st], raw[st], [], "none")
        rs.append(rn); plan.append((NEW_BASIS, "ctl", "none", st, 0.0))
    (out / NEW_BASIS / "run.pbs").write_text(
        AD.PBS % dict(tag="15", hpc=hpc, basis=NEW_BASIS, n=len(rs),
                      runs="\n".join("run %s" % r for r in rs)))
    with open(out / "PLAN.tsv", "w") as fh:
        fh.write("# %s D15.2; cs6_d15.py prepare\n" % PREREG)
        fh.write("# the six negatives of %s at arm %s, lambda %.6f, %s. wall-at-zero is reused\n"
                 % (DESIGN, NEW_ARM, lam, NEW_BASIS))
        fh.write("# from AMENDMENT 1 A1.4: it does not depend on lambda.\n")
        fh.write("basis\tsite\tarm\tstate\tabs_q\n")
        for r in plan:
            fh.write("%s\t%s\t%s\t%s\t%.4f\n" % r)
    files = sorted(p for p in out.rglob("*") if p.is_file())
    (out / "INPUT_CHECKSUMS.sha256").write_text(
        "".join("%s  %s\n" % (sha(p), p.relative_to(out)) for p in files))
    print("   arm %s, lambda %.6f, %s" % (NEW_ARM, lam, NEW_BASIS))
    for s in neg:
        print("      site %-5d q %+.4f -> %+.4f e" % (s, D["q"][i[s]], lam * D["q"][i[s]]))
    print("   %d run directories, %d files including the manifest" % (len(plan), len(files) + 1))
    print("   NOTHING WAS SUBMITTED BY THIS BUNDLE.")


# ---------------------------------------------------------------- analyse

def analyse(root, runs, outdir):
    T6, NEW, outdir = Path(root) / "phase2.2" / "cs_stage6", Path(runs), Path(outdir)
    if not (T6 / PREREG).is_file():
        stop("%s is not committed" % PREREG)
    ADD, NEGS = T6 / "additivity", T6 / "negatives"
    D = TR.design(T6, DESIGN, NSPH); REF = TR.design(T6, TR.REF_DESIGN, TR.REF_N)
    LAM = TR.lambdas(D, REF); C = TR.controls(T6)
    pos, neg = sites(D)
    GEO = {st: TR.read_xyz(T6 / "frame_42436" / fn) for st, fn in
           (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz"))}
    O = {st: NG.com(T6, st) for st in Q13.GEOM}
    idx = {s: k for k, s in enumerate(D["site"])}

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

    def tree_of(nm, arm, basis, a):
        """Which committed tree holds this run, per D15.2."""
        if nm == "n10" or nm.startswith("n") or int(nm[1:]) in pos:
            return ADD
        if arm == "wall0":
            return NEGS                       # lambda-independent, A1.4's
        return base if a == "A" else NEGS
    E, CL, WHY, NB = {}, {}, {}, {}
    def take(nm, arm, basis, a, st, mem, walled):
        tr = tree_of(nm, arm, basis, a)
        d = tr / basis / (("ctl_%s" % st) if nm == "ctl"
                          else "%s_%s" % (AD.rname(nm, arm, a if arm != "wall0" else None), st))
        if arm in ("bare", "none"):
            e, nb, dd = TR.one(d, GEO[st], len(mem)) if mem else AD.one0(d, GEO[st], 0, False)
        else:
            e, nb, dd = AD.one0(d, GEO[st], len(mem), walled)
        if walled:
            dd = dict(dd); dd["dip"] = NG.corrected_dip(d, mem, O[st])
        k = (basis, nm, arm, a, st)
        E[k], NB[k] = e, nb
        CL[k], WHY[k] = TR.verdict(dd, C)
        return tr

    CTL = {}
    for basis, a in CONDS:
        lam = LAM[AD.ARM_LAM[a]]
        for nm in ["n10"] + ["s%d" % s for s in pos + neg] + ["n2", "n3", "n4"]:
            mem = AD.members(D, (list(D["site"]) if nm == "n10" else
                                 [int(nm[1:])] if nm.startswith("s") else
                                 list(AD.POS[:int(nm[1:])])), lam)
            for st in Q13.GEOM:
                take(nm, "wall0", basis, a, st, NG.zeroed(mem), True)
                take(nm, "bare", basis, a, st, mem, False)
                take(nm, "wall", basis, a, st, mem, True)
        for tr, tag in ((ADD, "add"), (NEGS, "neg"), (base, "new")):
            if not (tr / basis / "ctl_R").is_dir():
                continue
            for st in Q13.GEOM:
                d = tr / basis / ("ctl_%s" % st)
                e, nb, dd = AD.one0(d, GEO[st], 0, False)
                CTL[(basis, tag, st)] = e
                CL[(basis, "ctl_" + tag, "none", a, st)], WHY[(basis, "ctl_" + tag, "none", a, st)] = \
                    TR.verdict(dd, C)

    H = TR.H
    def ctl_tag(nm, arm, a):
        tr = tree_of(nm, arm, "x", a)
        return "add" if tr == ADD else ("neg" if tr == NEGS else "new")
    def ddE(basis, nm, arm, a):
        c = {g: CTL[(basis, ctl_tag(nm, arm, a), g)] for g in Q13.GEOM}
        z = {g: (E[(basis, nm, "wall0", a, g)] - CTL[(basis, ctl_tag(nm, "wall0", a), g)])
             if arm == "wall" else 0.0 for g in Q13.GEOM}
        d = {g: E[(basis, nm, arm, a, g)] - c[g] - z[g] for g in Q13.GEOM}
        return (d["TS"] - d["R"]) * H
    def dsh(basis, nm, a):
        return ddE(basis, nm, "bare", a) - ddE(basis, nm, "wall", a)
    def well(basis, nm, a):
        return all(CL[(basis, nm, arm, a, g)] in ("CLEAN", "POLARISED")
                   for arm in ("wall0", "bare", "wall") for g in Q13.GEOM)

    O_ = ["STAGE 6 - DIAGNOSTIC 15 RESULT: does the per-site sum bound the whole mixed-sign design?",
          "",
          "Measured by cs6_d15.py against %s, committed before this" % PREREG,
          "diagnostic's own runs existed and applied unchanged. The criterion is DIAGNOSTIC 14.3's,",
          "carried over. Classification is D6.3 with AMENDMENT 1 A1.2's dipole correction. No",
          "geometry moved. In vacuo.",
          "",
          "D15.2 records that most of these runs were committed before the pre-registration was",
          "written: the four positive singles and n10 by DIAGNOSTIC 14, the six negative singles at",
          "arm B and every wall-at-zero run by AMENDMENT 1 A1.4. Only the six negatives bare and",
          "walled at arm A were made for this diagnostic. Every one is checked against its sha256",
          "manifest before it is read.", ""]

    O_ += ["GATE 1: EVERY CONTROL MUST AGREE WITH DIAGNOSTIC 14'S, TO %g Eh" % CTL_EH, ""]
    g1 = True
    for basis, a in CONDS:
        ref = CTL[(basis, "add", "R")]
        for tag in ("add", "neg", "new"):
            for st in Q13.GEOM:
                if (basis, tag, st) not in CTL:
                    continue
                ok = abs(CTL[(basis, tag, st)] - CTL[(basis, "add", st)]) <= CTL_EH
                g1 &= ok
                O_.append("   %-10s %-4s %-2s %18.10f   D14 %18.10f   %s"
                          % (basis, tag, st, CTL[(basis, tag, st)], CTL[(basis, "add", st)],
                             "pass" if ok else "FAIL"))
    O_.append("")
    if not g1:
        return finish(O_, outdir, "GATE 1 failed: the trees disagree on the control. Nothing is read.")

    W = AD.d7_table(T6)          # DIAGNOSTIC 7's own committed table, through cs6_additivity
    O_ += ["GATE 2: n10's BARE ddE MUST REPRODUCE DIAGNOSTIC 7's COMMITTED VALUE", ""]
    g2 = True
    for basis, a in CONDS:
        got, want = ddE(basis, "n10", "bare", a), W[(basis, AD.ARM_LAM[a])][1]
        ok = abs(got - want) <= D7_TOL; g2 &= ok
        O_.append("   %-10s arm %s   %+9.4f   D7 %+9.4f   %s" % (basis, a, got, want,
                                                                 "pass" if ok else "FAIL"))
    O_.append("")
    if not g2:
        return finish(O_, outdir, "GATE 2 failed: n10 does not reproduce DIAGNOSTIC 7. Nothing is read.")

    O_ += ["NOT READ, AND WHY (D15.2, stated before anything was computed)", "",
           "   %s arm %s: %s" % UNREAD, ""]

    O_ += ["EVERY RUN CLASSIFIED BY D6.3 WITH A1.2'S CORRECTION, READ BEFORE ANY NUMBER", "",
           "   run                                   class       note"]
    for k in sorted(CL, key=lambda k: (k[0], k[3], k[1], k[2], k[4])):
        O_.append("   %-37s %-10s  %s"
                  % ("%s/%s_%s_%s_%s" % (k[0], k[1], k[2], k[3], k[4]), CL[k], "; ".join(WHY[k])))
    O_.append("")

    for basis, a in CONDS:
        lam = LAM[AD.ARM_LAM[a]]
        O_ += ["CONDITION: %s, arm %s (lambda %.6f)" % (basis, a, lam), "",
               "   site    sign  |q| (e)   d_near    d_shield   class"]
        unread, S = [], {}
        for s in pos + neg:
            nm = "s%d" % s
            dn = Q13.descriptors(D["X"][idx[s]], Q13.read_xyz(T6 / "frame_42436" / "reactant_qm.xyz"))["dn"]
            if not well(basis, nm, a):
                O_.append("   %-6d  %-4s  %7.4f  %7.3f   NOT READABLE" %
                          (s, "+" if s in pos else "-", abs(lam * D["q"][idx[s]]), dn))
                unread.append(s); continue
            S[s] = dsh(basis, nm, a)
            O_.append("   %-6d  %-4s  %7.4f  %7.3f  %+9.4f   %s"
                      % (s, "+" if s in pos else "-", abs(lam * D["q"][idx[s]]), dn, S[s],
                         CL[(basis, nm, "wall", a, "R")]))
        O_.append("")
        if unread or not well(basis, "n10", a):
            O_ += ["   READING 1: NOT READABLE (%s)" % (("sites " + ", ".join(str(s) for s in unread))
                                                        if unread else "n10"), ""]
            continue
        Dsum = sum(abs(v) for v in S.values())
        Dset = abs(dsh(basis, "n10", a))
        tol = max(Q13.SHIELD_REL * Dsum, Q13.SHIELD_ABS)
        Dpos = sum(abs(S[s]) for s in pos); Dneg = sum(abs(S[s]) for s in neg)
        O_ += ["   Dsum %.4f   Dset %.4f   tol %.4f   Dset-Dsum %+.4f   Dset/Dsum %.3f"
               % (Dsum, Dset, tol, Dset - Dsum, (Dset / Dsum if Dsum else float("nan"))), ""]
        if Dsum <= tol and Dset <= tol:
            v1 = "NO EFFECT TO BOUND (Dsum %.4f and Dset %.4f both <= tol %.4f)" % (Dsum, Dset, tol)
        elif Dset <= Dsum + tol:
            v1 = "BOUNDED"
        else:
            v1 = "NOT BOUNDED"
        O_ += ["   READING 1: %s" % v1, "",
               "   READING 2 (reported, no criterion): Dpos %.4f over %d positives, Dneg %.4f over"
               % (Dpos, len(pos), Dneg),
               "              %d negatives, Dneg/Dsum %.3f" % (len(neg), (Dneg / Dsum if Dsum else float("nan"))),
               ""]
        rows = []
        for nm in ("n2", "n3", "n4"):
            mem = list(AD.POS[:int(nm[1:])])
            if not well(basis, nm, a) or any(s in unread for s in mem):
                rows.append((nm, None, None, None)); continue
            ds2 = sum(abs(S[s]) for s in mem); dt2 = abs(dsh(basis, nm, a))
            t2 = max(Q13.SHIELD_REL * ds2, Q13.SHIELD_ABS)
            rows.append((nm, ds2, dt2, "BOUNDED" if dt2 <= ds2 + t2 else "NOT BOUNDED"))
        O_ += ["   READING 3 (positive-only subsets; a NARROWER statement, not an answer to D15.1)",
               "      set   Dsum      Dset      verdict"]
        for nm, ds2, dt2, v in rows:
            O_.append("      %-4s  %s" % (nm, "NOT READABLE" if v is None
                                          else "%8.4f  %8.4f  %s" % (ds2, dt2, v)))
        O_.append("")
    return finish(O_, outdir, None)


def finish(O, outdir, tail):
    if tail:
        O += [tail, ""]
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "DIAGNOSTIC15_RESULT.txt").write_text("\n".join(O) + "\n")
    print("\n".join(O))
    print("\n   wrote %s" % (outdir / "DIAGNOSTIC15_RESULT.txt"))


if __name__ == "__main__":
    A_ = sys.argv[1:]
    if len(A_) == 4 and A_[0] == "prepare":
        prepare(A_[1], A_[2], A_[3])
    elif len(A_) == 4 and A_[0] == "analyse":
        analyse(A_[1], A_[2], A_[3])
    else:
        stop("usage: cs6_d15.py prepare <results_repo> <OUT> <HPC_DIR>\n"
             "       cs6_d15.py analyse <results_repo> <RUNS> <OUT_DIR>")
