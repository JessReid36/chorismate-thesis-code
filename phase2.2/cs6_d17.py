#!/usr/bin/env python3
"""
cs6_d17.py - STAGE6_DIAGNOSTIC17_20261030.txt.

Is the artefact linear in lambda on a fixed arrangement, and does n still matter once the
largest single charge is held equal?

  python3 cs6_d17.py prepare <results_repo> <OUT> <HPC_DIR>
  python3 cs6_d17.py analyse <results_repo> <RUNS> <OUT_DIR>

NOTHING ALREADY COMMITTED IS RETYPED OR RE-RUN.
  cs6_additivity  write_run, members, rname, one0, PBS
  cs6_negatives   zeroed, com, corrected_dip, check_correction      (A1.2's correction)
  cs6_transfer    design, head, controls, verdict, one, QCLEAN, H
  cs13_qscale     FRAME_SHA, GEOM, read_xyz_raw, SHIELD_REL, SHIELD_ABS
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

PREREG = "STAGE6_DIAGNOSTIC17_20261030.txt"
DESIGNS = (("S6-2", 2), ("S6-3", 3), ("S6-5", 5), ("S6-10", 10), ("S6-20", 20))
BASIS = "def2-SVPD"
CTL_EH = 1e-6
MATCHED = 0.538860          # D17.2: the matched largest single charge, derived not chosen
F1 = 1.3000                 # D16's matched total field, which fixes MATCHED through S6-3
REUSED = "S6-3"             # its D16 F1 runs ARE its D17 runs; lambda is identical

# D17.3 fixes these before any run exists. c is the line through the two COMMITTED points,
# so the new point is a prediction and not a fit.
PRED = {"S6-10": dict(c=0.8234, lam=0.979746, dset=0.8067),
        "S6-20": dict(c=0.7609, lam=0.979746, dset=0.7455)}
# D15's committed number, reproduced here from the committed runs rather than quoted (GATE 3)
D15_DSET10 = 0.4140


def stop(m):
    print("STOP: " + m); sys.stdout.flush(); sys.exit(1)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def setup(T6):
    """The five designs, the matched largest charge, and each design's lambda."""
    D = {nm: TR.design(T6, nm, n) for nm, n in DESIGNS}
    M = (F1 / D["S6-3"]["sabs"]) * D["S6-3"]["maxq"]
    if abs(M - MATCHED) > 5e-6:
        stop("the matched largest charge is %.6f e, not the %.6f D17.2 states" % (M, MATCHED))
    lam, s2 = {}, {}
    for nm, _ in DESIGNS:
        d = D[nm]
        q = np.asarray(d["q"], dtype=float)
        s2[nm] = float((q ** 2).sum())
        lam[nm] = M / d["maxq"]
        cap = min(1.0, TR.QCLEAN / d["maxq"], (2.0 / s2[nm]) ** 0.5)
        if lam[nm] > cap + 1e-9:
            stop("%s needs lambda %.6f, above its own cap %.6f - outside what D17.2 admits"
                 % (nm, lam[nm], cap))
    if abs(lam[REUSED] - F1 / D[REUSED]["sabs"]) > 1e-9:
        stop("%s's D17 lambda is not its committed D16 F1 lambda; its runs cannot be reused" % REUSED)
    return D, M, lam, s2


def pcfile(p):
    """charges and positions out of a committed pc.pc, sorted so two files can be compared."""
    rows = [l.split() for l in Path(p).read_text().split("\n")[1:] if l.strip()]
    a = np.array([[float(x) for x in r] for r in rows])
    o = np.lexsort(a[:, 1:].T)
    return a[o, 0], a[o, 1:]


def same_arrangement(pa, pb, expect):
    """Positions identical, charges in one constant ratio. Tolerances RELATIVE: pc.pc writes
    charges to six decimals, so an absolute test on the ratio is tighter than the file."""
    qa, Xa = pcfile(pa); qb, Xb = pcfile(pb)
    if len(qa) != len(qb):
        return False, "%d charges against %d" % (len(qa), len(qb))
    dX = float(np.abs(Xa - Xb).max())
    r = qa / qb
    rel = float((r.max() - r.min()) / abs(r.mean()))
    off = float(abs(r.mean() / expect - 1))
    ok = dX < 1e-6 and rel < 1e-4 and off < 1e-4
    return ok, ("positions %.1e A, charge ratio %.6f +/- %.1e relative, expected %.6f (off %.1e)"
                % (dX, r.mean(), rel, expect, off))


def rn(nm, arm):
    return "%s_%s_%s" % (nm, arm, "%s")       # filled with the state


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
    D, M, lam, s2 = setup(T6)
    print("   matched largest charge %.6f e, derived from %s at D16's %.4f e field" % (M, REUSED, F1))
    for nm, n in DESIGNS:
        print("      %-7s n %-3d lambda %.6f   sum |q| %7.4f e   max |q| %.4f e   sum q^2 %.4f%s"
              % (nm, n, lam[nm], lam[nm] * D[nm]["sabs"], lam[nm] * D[nm]["maxq"],
                 lam[nm] ** 2 * s2[nm], "   REUSED from D16" if nm == REUSED else ""))

    # every reused run must already be there, named, BEFORE anything is written
    missing = []
    for nm, _ in DESIGNS:
        for st in Q13.GEOM:
            if not (T6 / "d16" / BASIS / ("%s_wall0_%s" % (nm, st))).is_dir():
                missing.append("d16/%s/%s_wall0_%s" % (BASIS, nm, st))
        if nm == REUSED:
            for arm in ("bare", "wall"):
                for st in Q13.GEOM:
                    if not (T6 / "d16" / BASIS / ("%s_%s_F1_%s" % (nm, arm, st))).is_dir():
                        missing.append("d16/%s/%s_%s_F1_%s" % (BASIS, nm, arm, st))
    if missing:
        stop("D17 reuses runs that are not committed:\n   " + "\n   ".join(missing))
    print("   all five wall-at-zero arms and %s's bare and walled arms are committed in d16" % REUSED)

    out.mkdir(parents=True)
    plan, rs = [], []
    for nm, n in DESIGNS:
        if nm == REUSED:
            continue
        mem = AD.members(D[nm], D[nm]["site"], lam[nm])
        for arm, walled in (("bare", False), ("wall", True)):
            for st in Q13.GEOM:
                r = "%s_%s_M_%s" % (nm, arm, st)
                AD.write_run(out / BASIS / r, BASIS, geo[st], raw[st], mem, arm)
                rs.append(r)
                plan.append((nm, n, arm, st, len(mem), lam[nm],
                             sum(abs(q) for _, q, _ in mem), max(abs(q) for _, q, _ in mem)))
    for st in Q13.GEOM:
        AD.write_run(out / BASIS / ("ctl_%s" % st), BASIS, geo[st], raw[st], [], "none")
        rs.append("ctl_%s" % st)
        plan.append(("ctl", 0, "none", st, 0, 0.0, 0.0, 0.0))
    (out / BASIS / "run.pbs").write_text(
        AD.PBS % dict(tag="17", hpc=hpc, basis=BASIS, n=len(rs),
                      runs="\n".join("run %s" % r for r in rs)))
    with open(out / "PLAN.tsv", "w") as fh:
        fh.write("# %s D17.2; cs6_d17.py prepare\n" % PREREG)
        fh.write("# matched largest charge %.6f e; %s\n" % (M, BASIS))
        fh.write("# %s is reused from d16 unchanged and is not re-run\n" % REUSED)
        fh.write("what\tn\tarm\tstate\tnq\tlambda\tsum_abs_q\tmax_abs_q\n")
        for r in plan:
            fh.write("%s\t%d\t%s\t%s\t%d\t%.6f\t%.4f\t%.4f\n" % r)
    files = sorted(p for p in out.rglob("*") if p.is_file())
    (out / "INPUT_CHECKSUMS.sha256").write_text(
        "".join("%s  %s\n" % (sha(p), p.relative_to(out)) for p in files))
    print("   %d run directories, %d files including the manifest" % (len(plan), len(files) + 1))
    if len(plan) != 18:
        stop("%d runs, but D17.2 states 18" % len(plan))
    print("   NOTHING WAS SUBMITTED BY THIS BUNDLE.")


# ---------------------------------------------------------------- analyse

def analyse(root, runs, outdir):
    T6, NEW, outdir = Path(root) / "phase2.2" / "cs_stage6", Path(runs), Path(outdir)
    if not (T6 / PREREG).is_file():
        stop("%s is not committed" % PREREG)
    D, M, lam, s2 = setup(T6)
    C = TR.controls(T6)
    GEO = {st: TR.read_xyz(T6 / "frame_42436" / fn) for st, fn in
           (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz"))}
    O = {st: NG.com(T6, st) for st in Q13.GEOM}
    H = TR.H

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

    CL, WHY, CTL = {}, {}, {}
    TREE = {"new": base, "d16": T6 / "d16", "additivity": T6 / "additivity"}
    for tag, tr in TREE.items():
        if (tr / BASIS / "ctl_R").is_dir():
            for st in Q13.GEOM:
                e, nb, dd = AD.one0(tr / BASIS / ("ctl_%s" % st), GEO[st], 0, False)
                CTL[(tag, st)] = e
                kk = ("ctl_" + tag, "ctl_%s_%s" % (tag, st), st)
                CL[kk], WHY[kk] = TR.verdict(dd, C)

    def arm_energy(tag, run, mem, arm):
        """one arm at one state, classified; the walled arms get A1.2's dipole correction."""
        d = TREE[tag] / BASIS / run
        if not d.is_dir():
            stop("%s is not in the %s tree" % (run, tag))
        nsph = len(mem)
        if arm == "bare":
            e, nb, dd = TR.one(d, GEO[run[-2:].lstrip("_")], len(mem))
        else:
            e, nb, dd = AD.one0(d, GEO[run[-2:].lstrip("_")], nsph, True)
            dd = dict(dd); dd["dip"] = NG.corrected_dip(d, mem, O[run[-2:].lstrip("_")])
        return e, dd

    def point(label, nm, tag, bare_stem, wall_stem, w0_tag, w0_stem, lamv):
        """d_shield for one arrangement at one lambda, from whichever trees hold its runs."""
        mem = AD.members(D[nm], D[nm]["site"], lamv)
        zer = NG.zeroed(mem)
        e, cls = {}, []
        for st in Q13.GEOM:
            for tg, stem, mm, arm in ((tag, bare_stem, mem, "bare"),
                                      (tag, wall_stem, mem, "wall"),
                                      (w0_tag, w0_stem, zer, "wall")):
                run = "%s_%s" % (stem, st)
                d = TREE[tg] / BASIS / run
                if not d.is_dir():
                    stop("%s is not in the %s tree; D17 expects it there" % (run, tg))
                if arm == "bare":
                    ee, nb, dd = TR.one(d, GEO[st], len(mm))
                else:
                    ee, nb, dd = AD.one0(d, GEO[st], len(mm), True)
                    dd = dict(dd); dd["dip"] = NG.corrected_dip(d, mm, O[st])
                e[(stem, st)] = (ee, tg)
                v, w = TR.verdict(dd, C)
                CL[(label, run, st)] = v; WHY[(label, run, st)] = w
                cls.append(v)
        def dd_of(stem, zero_stem):
            out = {}
            for st in Q13.GEOM:
                ee, tg = e[(stem, st)]
                z = 0.0
                if zero_stem is not None:
                    e0, t0 = e[(zero_stem, st)]
                    z = e0 - CTL[(t0, st)]
                out[st] = ee - CTL[(tg, st)] - z
            return (out["TS"] - out["R"]) * H
        ddb = dd_of(bare_stem, None)
        ddw = dd_of(wall_stem, w0_stem)
        well = all(v in ("CLEAN", "POLARISED") for v in cls)
        return dict(lam=lamv, bare=ddb, wall=ddw, dset=abs(ddb - ddw), signed=ddb - ddw,
                    well=well, cls=cls, klass="/".join(sorted(set(cls))))

    O_ = ["STAGE 6 - DIAGNOSTIC 17 RESULT: is the artefact linear in lambda, and does n still",
          "matter at equal largest charge?", "",
          "Measured by cs6_d17.py against %s, committed before any of this" % PREREG,
          "diagnostic's runs existed and applied unchanged. Classification is D6.3 with AMENDMENT 1",
          "A1.2's dipole correction. No geometry moved. In vacuo. %s throughout." % BASIS, "",
          "D17.4 states the largest limitation first and it is restated here: Reading 2 holds the",
          "largest charge equal and lets everything else vary - the number of charges, where they",
          "sit, their signs, and the total field, 1.0777 e at n = 2 to 5.6433 e at n = 20. It removes",
          "one confound from DIAGNOSTIC 16 and puts another in its place.", "",
          "   matched largest charge %.6f e, derived from %s at D16's %.4f e field" % (M, REUSED, F1), ""]

    O_ += ["GATE 1: EVERY CONTROL MUST AGREE WITH DIAGNOSTIC 14'S, TO %g Eh" % CTL_EH, ""]
    g1 = True
    for tag in ("new", "d16", "additivity"):
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

    O_ += ["GATE 2: n10 IN additivity MUST BE THE SAME ARRANGEMENT AS S6-10", "",
           "   D17 takes S6-10's middle lambda from DIAGNOSTIC 14's n10 runs. The names differ and",
           "   n2 and n3 are NOT S6-2 and S6-3, so the arrangement is checked and not assumed.", ""]
    exp = 0.500000 / (F1 / D["S6-10"]["sabs"])
    ok2, msg2 = same_arrangement(T6 / "additivity" / BASIS / "n10_bare_A_R" / "pc.pc",
                                 T6 / "d16" / BASIS / "S6-10_bare_F1_R" / "pc.pc", exp)
    O_ += ["   %s" % msg2, "   GATE 2: %s" % ("pass" if ok2 else "FAIL"), ""]
    if not ok2:
        return finish(O_, outdir, "GATE 2 failed: n10 is not S6-10. Its lambda 0.50 point is not used.")

    SER = {"S6-10": [("d16", "S6-10_bare_F1", "S6-10_wall_F1", "d16", "S6-10_wall0", F1 / D["S6-10"]["sabs"]),
                     ("additivity", "n10_bare_A", "n10_wall_A", "additivity", "n10_wall0", 0.500000),
                     ("new", "S6-10_bare_M", "S6-10_wall_M", "d16", "S6-10_wall0", lam["S6-10"])],
           "S6-20": [("d16", "S6-20_bare_F1", "S6-20_wall_F1", "d16", "S6-20_wall0", F1 / D["S6-20"]["sabs"]),
                     ("d16", "S6-20_bare_F2", "S6-20_wall_F2", "d16", "S6-20_wall0", 2.1900 / D["S6-20"]["sabs"]),
                     ("new", "S6-20_bare_M", "S6-20_wall_M", "d16", "S6-20_wall0", lam["S6-20"])]}
    P = {nm: [point(nm, nm, tg, b, w, w0t, w0, lv) for tg, b, w, w0t, w0, lv in ser]
         for nm, ser in SER.items()}

    O_ += ["GATE 3: THE lambda 0.50 POINT MUST REPRODUCE DIAGNOSTIC 15'S COMMITTED %.4f" % D15_DSET10, ""]
    r = P["S6-10"][1]["dset"]
    ok3 = abs(r - D15_DSET10) <= 5e-4
    O_ += ["   recomputed from additivity's n10 arm A: %.4f   D15 %.4f   %s"
           % (r, D15_DSET10, "pass" if ok3 else "FAIL"), ""]
    if not ok3:
        return finish(O_, outdir, "GATE 3 failed: the committed runs no longer give D15's number.")

    O_ += ["GATE 4: THE TWO COMMITTED POINTS MUST REPRODUCE THE c D17.3 FIXED", "",
           "   D17.3's predictions were derived from these two points. If the committed runs no",
           "   longer give that c, the prediction is not the one the pre-registration states and",
           "   nothing is read against it.", ""]
    g4 = True
    for nm in ("S6-10", "S6-20"):
        two = P[nm][:2]
        if not all(q["well"] for q in two):
            O_.append("   %-6s a committed point is unreadable   FAIL" % nm); g4 = False; continue
        c = sum(q["dset"] / q["lam"] for q in two) / 2.0
        ok = abs(c - PRED[nm]["c"]) <= 1e-3
        g4 &= ok
        O_.append("   %-6s c recomputed %.4f   D17.3 %.4f   difference %+.4f   %s"
                  % (nm, c, PRED[nm]["c"], c - PRED[nm]["c"], "pass" if ok else "FAIL"))
    O_.append("")
    if not g4:
        return finish(O_, outdir, "GATE 4 failed: the committed points no longer give D17.3's c. "
                                  "Reading 1 is not read; the other readings are not attempted.")

    REST = {}
    for nm, _ in DESIGNS:
        if nm in P:
            REST[nm] = P[nm][-1]
        elif nm == REUSED:
            REST[nm] = point(nm, nm, "d16", "%s_bare_F1" % nm, "%s_wall_F1" % nm,
                             "d16", "%s_wall0" % nm, lam[nm])
        else:
            REST[nm] = point(nm, nm, "new", "%s_bare_M" % nm, "%s_wall_M" % nm,
                             "d16", "%s_wall0" % nm, lam[nm])

    O_ += ["EVERY RUN CLASSIFIED BY D6.3 WITH A1.2'S CORRECTION, READ BEFORE ANY NUMBER", "",
           "   run                                   class       note"]
    for k in sorted(CL, key=lambda k: (str(k[0]), str(k[1]))):
        nmk = k[1] if len(k) > 2 else k[0]
        O_.append("   %-37s %-10s  %s" % (nmk, CL[k], "; ".join(WHY[k])))
    O_.append("")

    # READING 1
    O_ += ["READING 1: IS Dset LINEAR IN lambda ON A FIXED ARRANGEMENT?", "",
           "   The line is fitted to the two COMMITTED points alone. The new point is a prediction",
           "   written into D17.3 before it existed, not a fit.", ""]
    dirs = []
    for nm in ("S6-10", "S6-20"):
        pts = P[nm]
        O_ += ["   %s" % nm,
               "      lambda      Dset      source",
               ]
        for p, src in zip(pts, ("DIAGNOSTIC 16, committed", "DIAGNOSTIC 14 n10 arm A, committed"
                                if nm == "S6-10" else "DIAGNOSTIC 16, committed", "new")):
            O_.append("      %.6f  %s   %s" % (p["lam"],
                      ("%8.4f" % p["dset"]) if p["well"] else "NOT READ", src))
        new = pts[-1]
        if not new["well"]:
            O_ += ["      the new point is unreadable: %s" % ", ".join(sorted(set(new["cls"]))), ""]
            dirs.append(None); continue
        pr = PRED[nm]
        tol = max(Q13.SHIELD_REL * max(p["dset"] for p in pts if p["well"]), Q13.SHIELD_ABS)
        diff = new["dset"] - pr["dset"]
        O_ += ["      c from the two committed points, fixed in D17.3: %.4f kcal/mol per unit lambda" % pr["c"],
               "      predicted %.4f   measured %.4f   difference %+.4f   tol %.4f"
               % (pr["dset"], new["dset"], diff, tol), ""]
        dirs.append(0 if abs(diff) <= tol else (1 if diff > 0 else -1))
    if any(d is None for d in dirs):
        v1 = "NOT READABLE (a design is unreadable at the new lambda)"
    elif all(d == 0 for d in dirs):
        v1 = "LINEAR"
    elif all(d > 0 for d in dirs):
        v1 = "SUPERLINEAR"
    elif all(d < 0 for d in dirs):
        v1 = "SUBLINEAR"
    else:
        v1 = "SPLIT (the two designs disagree in direction)"
    O_ += ["   READING 1: %s" % v1, ""]

    # READING 2
    O_ += ["READING 2: DOES THE ARTEFACT DEPEND ON n AT EQUAL LARGEST CHARGE (%.6f e)?" % M, "",
           "   design   n    lambda     sum|q|    Dset       class"]
    pts2 = []
    for nm, n in DESIGNS:
        p = REST[nm]
        if not p["well"]:
            O_.append("   %-7s %-4d %.6f  %7.4f    NOT READABLE   %s"
                      % (nm, n, p["lam"], p["lam"] * D[nm]["sabs"], ", ".join(sorted(set(p["cls"])))))
            continue
        pts2.append((n, nm, p["dset"]))
        O_.append("   %-7s %-4d %.6f  %7.4f   %8.4f   %s"
                  % (nm, n, p["lam"], p["lam"] * D[nm]["sabs"], p["dset"], p["klass"]))
    O_.append("")
    if len(pts2) < 3:
        v2 = "NOT READABLE (%d of 5 designs readable)" % len(pts2)
    else:
        tol2 = max(Q13.SHIELD_REL * max(v for _, _, v in pts2), Q13.SHIELD_ABS)
        vals = [v for _, _, v in sorted(pts2)]
        spread = max(vals) - min(vals)
        O_ += ["   tol = max(0.10 x largest Dset, %.2f) = %.4f;  spread %.4f"
               % (Q13.SHIELD_ABS, tol2, spread), ""]
        if max(vals) <= tol2:
            v2 = "NO EFFECT TO COMPARE (every Dset <= tol)"
        elif spread <= tol2:
            v2 = "FLAT"
        elif all(b >= a for a, b in zip(vals, vals[1:])):
            v2 = "GROWS WITH n"
        elif all(b <= a for a, b in zip(vals, vals[1:])):
            v2 = "FALLS WITH n"
        else:
            v2 = "NOT MONOTONIC"
    O_ += ["   READING 2: %s" % v2, ""]

    # READING 3
    O_ += ["READING 3: WHERE DOES A SINGLE CHARGE STOP BEING USABLE?", "",
           "   S6-2's bare arm was STRAINED at 0.6500 e in DIAGNOSTIC 16. Here it carries %.6f e." % M, ""]
    p2 = REST["S6-2"]
    if p2["well"]:
        v3 = "READABLE AT %.4f (the usable bare single-site charge lies between %.4f and 0.6500 e)" % (M, M)
    else:
        v3 = "STILL UNREADABLE at %.4f (it lies below that)" % M
    for st in Q13.GEOM:
        for stem in ("S6-2_bare_M", "S6-2_wall_M", "S6-2_wall0"):
            k = ("S6-2", "%s_%s" % (stem, st), st)
            if k in CL:
                O_.append("   %-28s %-10s  %s" % ("%s_%s" % (stem, st), CL[k], "; ".join(WHY[k])))
    O_ += ["", "   READING 3: %s" % v3, "",
           "   This bounds one number and nothing more. One site in one design at one geometry.", ""]
    return finish(O_, outdir, None)


def finish(O, outdir, tail):
    if tail:
        O += [tail, ""]
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "DIAGNOSTIC17_RESULT.txt").write_text("\n".join(O) + "\n")
    print("\n".join(O))
    print("\n   wrote %s" % (outdir / "DIAGNOSTIC17_RESULT.txt"))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        stop("usage: cs6_d17.py prepare|analyse ...")
    if sys.argv[1] == "prepare":
        prepare(*sys.argv[2:5])
    elif sys.argv[1] == "analyse":
        analyse(*sys.argv[2:5])
    else:
        stop("unknown command %s" % sys.argv[1])
