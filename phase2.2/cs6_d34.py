#!/usr/bin/env python3
"""DIAGNOSTIC 34: the production barriers at CAM-B3LYP/def2-SVPD, and the three readable
criterion legs on the committed frozen outputs.

Nothing here is reimplemented. The energies are read with the committed cs6_analyse.py's
read_run, which is where the bare / (MM) / (QM/MM) convention and the 118-electron, SCF,
echoed-coordinate and MM-point-charge checks live. The criterion legs are the committed
cs6_audit.py's read() at its committed thresholds, and the stability block and the level-
shift markers are the committed cs6_d33b.py's, imported rather than copied, so these
verdicts are comparable with D31's, D32's and D33's.

Stage A is POST HOC: the criterion postdates the frozen runs by two weeks. Its vocabulary
is fixed by the pre-registration as FLAGGED or NOT FLAGGED (3 LEGS) and can never be
RELIABLE, because no frozen output carries a stability analysis.
"""
import argparse, os, re, shutil, sys, tarfile, tempfile
from pathlib import Path

NQM = 24
STAGEA_TOL = 1e-4      # kcal/mol, on reproducing every committed ddE
CTL_TOL_EH = 1e-6      # G34.5, per leg
CTL_TOL_KCAL = 0.0100  # G34.5, on the barrier
D33A_BARRIER = 26.8248


def stop(m):
    sys.stderr.write("STOP  %s\n" % m)
    raise SystemExit(1)


def nspheres(prms):
    if not os.path.isfile(prms):
        return 0
    L = open(prms).read().splitlines()
    i = L.index("$atoms")
    return int(L[i + 1].split()[0]) - NQM


def sphere_q(prms):
    if not os.path.isfile(prms):
        return []
    L = open(prms).read().splitlines()
    i = L.index("$atoms")
    n = int(L[i + 1].split()[0])
    return [float(L[i + 2 + k].split()[2]) for k in range(n)][NQM:]


def stage_frozen(s6, work):
    """The committed frozen outputs, with their committed inputs beside them, so the
    committed readers can read them exactly as they read a live run directory."""
    tp = os.path.join(s6, "frozen", "sp_outputs.tar.xz")
    with tarfile.open(tp) as tf:
        for m in tf.getmembers():
            if re.match(r"\./[^/]+/(R|TS)_sp/", m.name):
                tf.extract(m, work, filter="data")
    got = []
    for a in sorted(os.listdir(work)):
        for leg in ("R", "TS"):
            d = os.path.join(work, a, "%s_sp" % leg)
            if not os.path.isdir(d):
                continue
            src = os.path.join(s6, "relax", a, "%s_sp" % leg)
            for fn in ("combined.xyz", "site.ORCAFF.prms"):
                p = os.path.join(src, fn)
                if os.path.isfile(p):
                    shutil.copyfile(p, os.path.join(d, fn))
            got.append((a, leg, d))
    if not got:
        stop("%s held no frozen output directories." % tp)
    return got


def legs(g, c, arm_label, rep, nm, leg, audit, flags):
    """The three readable legs, exactly as cs6_d33b.py applies them."""
    held = False
    if not g["norm"]:
        flags.append("%s abnormal termination" % leg)
    if g["notconv"]:
        flags.append("%s did not converge" % leg)
    rise = g["maxlw"] - c["maxlw"]
    if rise >= audit.D_LW:
        flags.append("%s relocation: rise %+.4f" % (leg, rise))
    if c["dip"] and g["dip"] and g["dip"] / c["dip"] >= audit.R_DIP:
        flags.append("%s dipole %.2fx control" % (leg, g["dip"] / c["dip"]))
    if g["shifted"] or g["gap"] is None or c["gap"] is None:
        held = True
    elif g["gap"] < audit.GAP_ABS and g["gap"] < audit.R_GAP * c["gap"]:
        flags.append("%s gap %.4f eV" % (leg, g["gap"]))
    return rise, held


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", required=True)
    ap.add_argument("--code", required=True)
    ap.add_argument("--prereg", required=True)
    ap.add_argument("--result", required=True)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--arm2note", required=True,
                    help="written only if arm 2 failed: ORCA's own words, for the committed note")
    a = ap.parse_args()
    S6 = "phase2.2/cs_stage6"
    s6 = os.path.join(a.res, S6)
    sys.path.insert(0, os.path.join(a.code, "phase2.2"))
    d33b = __import__("cs6_d33b")
    audit = d33b.load(os.path.join(a.code, "phase2.2", "cs6_audit.py"), "cs6_audit")
    anal = d33b.load(os.path.join(a.code, "phase2.2", "cs6_analyse.py"), "cs6_analyse")
    H = d33b.H2KCAL
    CONV = d33b.CONV
    gates, rep, rows = [], [], []

    say = rep.append
    say("STAGE 6, DIAGNOSTIC 34 - RESULT")
    say("The production barriers at CAM-B3LYP/def2-SVPD, and the three readable criterion")
    say("legs on the committed frozen outputs.")
    say("")
    say("Read by phase2.2/cs6_d34.py. Nothing is reimplemented:")
    say("   energies            cs6_analyse.py read_run      sha256 %s" % d33b.sha(os.path.join(a.code, "phase2.2", "cs6_analyse.py")))
    say("   criterion legs      cs6_audit.py read            sha256 %s" % d33b.sha(os.path.join(a.code, "phase2.2", "cs6_audit.py")))
    say("   stability, shift    cs6_d33b.py                  sha256 %s" % d33b.sha(os.path.join(a.code, "phase2.2", "cs6_d33b.py")))
    say("   pre-registration    %s" % os.path.basename(a.prereg))
    say("                                                    sha256 %s" % d33b.sha(a.prereg))
    say("   D_LW %.2f e   R_DIP %.1f   GAP_ABS %.1f eV   R_GAP %.1f   (committed thresholds)"
        % (audit.D_LW, audit.R_DIP, audit.GAP_ABS, audit.R_GAP))
    say("")

    # ---- the committed effects, and each arrangement's own charge budget
    want = {}
    for l in open(os.path.join(s6, "frozen", "frozen.tsv")):
        f = l.rstrip("\n").split("\t")
        if f and f[0] not in ("design", "") and not f[0].startswith("#"):
            want[f[0]] = (int(f[1]), float(f[2]))
    budget, nsph = {}, {}
    for arr in list(want) + ["S6-C"]:
        p = os.path.join(s6, "relax", arr, "R_sp", "site.ORCAFF.prms")
        q = sphere_q(p)
        budget[arr], nsph[arr] = sum(abs(x) for x in q), len(q)

    # =========================== STAGE A ===========================
    say("=" * 78)
    say("R34.1-R34.3  STAGE A. THE THREE READABLE LEGS ON THE COMMITTED FROZEN OUTPUTS.")
    say("             NO COMPUTE. POST HOC: the criterion postdates these runs.")
    say("=" * 78)
    work = tempfile.mkdtemp(prefix="d34a_")
    try:
        got = stage_frozen(s6, work)
        A, EA = {}, {}
        for arr, leg, d in got:
            c = d33b.run(d, audit)
            if c is None:
                gates.append("stage A: %s %s would not read" % (arr, leg))
                continue
            A[(arr, leg)] = c
            try:
                tot, mm, _ = anal.read_run(d, nsph[arr])
            except SystemExit as e:
                gates.append("stage A: %s %s failed the committed reader's own checks" % (arr, leg))
                continue
            EA[(arr, leg)] = (tot, mm)
        say("")
        say("   %d of 22 outputs read. Loewdin rows per output: the 24 QM atoms only, so the" % len(A))
        say("   relocation leg is the same quantity D31, D32 and D33 measure.")
        nstab = sum(1 for k in A if A[k]["stab"]["ran"])
        say("   stability analyses present: %d of %d. The fourth leg is therefore unreadable" % (nstab, len(A)))
        say("   and is NOT substituted. No verdict below is RELIABLE; none can be.")
        if nstab:
            gates.append("stage A: %d frozen outputs carry a stability analysis, "
                         "which the pre-registration states none does" % nstab)
        say("")
        # G34.8, a control test: the reader must reproduce every committed ddE
        if ("S6-C", "R") not in EA or ("S6-C", "TS") not in EA:
            stop("stage A: the control did not read; nothing downstream is computable.")
        cR, cT = EA[("S6-C", "R")][0], EA[("S6-C", "TS")][0]
        bar_old = (cT - cR) * H
        say("   G34.8  the control barrier recomputed from the committed outputs: %.4f kcal/mol" % bar_old)
        say("          FROZEN_REPORT.txt and A_REPORT.txt state 18.4610")
        if abs(bar_old - 18.4610) > 5e-4:
            gates.append("stage A: the control barrier reads %.4f, not 18.4610" % bar_old)
        worst, worst_nm = 0.0, "-"
        for arr, (n, dde_c) in want.items():
            if (arr, "R") not in EA or (arr, "TS") not in EA:
                gates.append("stage A: %s did not read on both legs" % arr)
                continue
            dde = ((EA[(arr, "TS")][0] - cT) - (EA[(arr, "R")][0] - cR)) * H
            if abs(dde - dde_c) > worst:
                worst, worst_nm = abs(dde - dde_c), arr
        say("   G34.8  every committed ddE reproduced to within %.1e kcal/mol (worst %s)"
            % (worst, worst_nm))
        if worst > STAGEA_TOL:
            gates.append("stage A: a committed ddE is not reproduced; worst %s by %.2e kcal/mol"
                         % (worst_nm, worst))
        say("")
        say("   %-22s %3s %8s  %10s %10s  %-5s %-30s" % ("arrangement", "n", "Sigma|q|", "max rise", "max dip", "shift", "verdict (ERRATUM 1)"))
        for arr in sorted(want, key=lambda x: (budget[x], x)):
            if (arr, "R") not in A or (arr, "TS") not in A:
                continue
            flags, held, rise_max, dip_max = [], False, None, None
            for leg in ("R", "TS"):
                r, h = legs(A[(arr, leg)], A[("S6-C", leg)], "frozen", rep, arr, leg, audit, flags)
                held = held or h
                rise_max = r if rise_max is None else max(rise_max, r)
                dr = (A[(arr, leg)]["dip"] / A[("S6-C", leg)]["dip"]
                      if A[("S6-C", leg)]["dip"] else None)
                if dr is not None:
                    dip_max = dr if dip_max is None else max(dip_max, dr)
            # ERRATUM 1: a run that converged under a level shift has a dark gap leg, so it
            # is not "not flagged on three legs". The term is narrower, never wider.
            if flags:
                v = "FLAGGED"
            elif held:
                v = "NOT FLAGGED (2 LEGS, GAP DARK)"
            else:
                v = "NOT FLAGGED (3 LEGS)"
            sh = "shift" if any(A[(arr, lg)]["shifted"] for lg in ("R", "TS")) else "-"
            say("   %-22s %3d %7.4f e  %+10.4f %9.2fx  %-5s %-30s"
                % (arr, want[arr][0], budget[arr], rise_max, dip_max or float("nan"), sh, v))
            for f in flags:
                say("        flag: %s" % f)
            rows.append(dict(stage="A", name=arr, n=want[arr][0], sumq=budget[arr],
                             verdict=v, rise=rise_max, dip=dip_max,
                             dde=((EA[(arr, "TS")][0] - cT) - (EA[(arr, "R")][0] - cR)) * H
                             if (arr, "TS") in EA and (arr, "R") in EA else None))
        nflag = sum(1 for r in rows if r["stage"] == "A" and r["verdict"] == "FLAGGED")
        say("")
        say("   R34.2  %d of %d designs FLAGGED on the three readable legs."
            % (nflag, sum(1 for r in rows if r["stage"] == "A")))
        ndark = sum(1 for r in rows if r["stage"] == "A" and "GAP DARK" in r["verdict"])
        say("          %d of them have a dark gap leg: the run converged under a level shift," % ndark)
        say("          so only relocation and dipole could be evaluated (ERRATUM 1).")
        say("   R34.3  a flag is a statement about the density. Nothing is withdrawn here;")
        say("          withdrawal of a committed number needs its own dated erratum.")
    finally:
        shutil.rmtree(work, ignore_errors=True)

    # =========================== STAGE B ===========================
    say("")
    say("=" * 78)
    say("R34.4-R34.10  STAGE B. 44 RUNS AT CAM-B3LYP/def2-SVPD, TWO ARMS.")
    say("=" * 78)
    d34 = os.path.join(s6, "d34")
    B, EB, missing = {}, {}, {"1": [], "2": []}
    for arm in ("1", "2"):
        for arr in list(want) + ["S6-C"]:
            for leg in ("R", "TS"):
                d = os.path.join(d34, "a%s_%s_%s" % (arm, arr, leg))
                if not os.path.isdir(d):
                    missing[arm].append("%s %s" % (arr, leg))
                    continue
                c = d33b.run(d, audit)
                if c is None:
                    missing[arm].append("%s %s" % (arr, leg))
                    continue
                B[(arm, arr, leg)] = c
                try:
                    tot, mm, _ = anal.read_run(d, nsph[arr])
                    EB[(arm, arr, leg)] = (tot, mm)
                except SystemExit:
                    missing[arm].append("%s %s (failed the committed reader)" % (arr, leg))
    say("")
    say("   arm 1, the lift            read %2d of 22" % sum(1 for k in EB if k[0] == "1"))
    say("   arm 2, with STABPerform    read %2d of 22" % sum(1 for k in EB if k[0] == "2"))

    if missing["1"]:
        for m in missing["1"]:
            gates.append("arm 1: %s did not read" % m)
        say("")
        say("   ARM 1 IS INCOMPLETE. Stage A above is complete and was read from committed")
        say("   outputs alone; stage B cannot be read without arm 1, so no result is written.")
        print("\n".join(rep))
        stop("arm 1 is incomplete (%d of 22 runs not read). Stage A is shown above and is "
             "complete. No result file was written and nothing is committable."
             % len(missing["1"]))

    arm2_ok = not missing["2"]
    arm2_reason = []
    if not arm2_ok:
        say("")
        say("   ARM 2 IS INCOMPLETE and is reported as a failure, as D34 section 4 fixed in")
        say("   advance. Not read: %s" % ", ".join(missing["2"][:8]) + (" ..." if len(missing["2"]) > 8 else ""))
        say("   Arm 1 stands. A stability-carrying frozen arm needs its own amendment.")
        arm2_reason.append("arm 2 runs that did not read: %s" % ", ".join(missing["2"]))
        say("")
        say("   ORCA's own words, from the first arm 2 output that carries them:")
        shown = 0
        for arr in list(want) + ["S6-C"]:
            for leg in ("R", "TS"):
                if shown:
                    break
                p2 = os.path.join(d34, "a2_%s_%s" % (arr, leg), "job.out")
                if not os.path.isfile(p2):
                    continue
                tx = open(p2, errors="replace").read()
                hit = [l for l in tx.splitlines()
                       if re.search(r"error|ERROR|aborting|not a valid|unknown|UNRECOGNIZED", l)]
                if hit:
                    for l in hit[:6]:
                        say("      %s" % l.strip()[:72])
                        arm2_reason.append(l.strip())
                    shown = 1
        if not shown:
            say("      none of the arm 2 outputs carries an error line.")
            arm2_reason.append("no arm 2 output carries an error line")
    else:
        nst = sum(1 for k in B if k[0] == "2" and B[k]["stab"]["ran"])
        say("   arm 2 stability analyses present: %d of 22" % nst)
        if nst != 22:
            arm2_ok = False
            say("")
            say("   ORCA ACCEPTED THE INPUT BUT RAN NO STABILITY ANALYSIS in %d of 22 runs, so" % (22 - nst))
            say("   the fourth leg is still unreadable and arm 2 carries no verdict. Arm 1 stands.")

    # G34.5, the control test
    say("")
    cR1, cT1 = EB[("1", "S6-C", "R")][0], EB[("1", "S6-C", "TS")][0]
    bar_new = (cT1 - cR1) * H
    dR, dT = None, None
    for leg, other in (("R", "sp_42436_R"), ("TS", "sp_42436_TS")):
        p = os.path.join(s6, "d33a", other, "job.out")
        if os.path.isfile(p):
            m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", open(p, errors="replace").read())
            if m:
                e = float(m[-1])
                d = abs(e - (cR1 if leg == "R" else cT1))
                if leg == "R":
                    dR = d
                else:
                    dT = d
    say("   G34.5  arm 1's control barrier %.4f kcal/mol against D33 stage A's %.4f"
        % (bar_new, D33A_BARRIER))
    if dR is not None and dT is not None:
        say("          its two energies differ from d33a's committed ones by %.2e and %.2e Eh" % (dR, dT))
        if max(dR, dT) > CTL_TOL_EH:
            gates.append("G34.5: arm 1's control energies differ from d33a's by up to %.2e Eh, "
                         "over the %.0e limit" % (max(dR, dT), CTL_TOL_EH))
    else:
        gates.append("G34.5: d33a's committed control outputs could not be read for comparison")
    if abs(bar_new - D33A_BARRIER) > CTL_TOL_KCAL:
        gates.append("G34.5: arm 1's control barrier is %.4f, not %.4f within %.4f"
                     % (bar_new, D33A_BARRIER, CTL_TOL_KCAL))
    if nsph["S6-C"] == 0 and EB[("1", "S6-C", "R")][1] != 0.0:
        gates.append("G34.6: the control reports an MM energy")

    # R34.4, R34.6, R34.8, R34.9
    say("")
    say("   R34.4 and R34.6. ddE at the corrected level, against the committed value.")
    say("   'left' is the barrier remaining against the reference at each level: 18.4610")
    say("   committed, %.4f here." % bar_new)
    say("")
    say("   %-22s %3s %8s  %10s %10s %9s   %8s %8s" %
        ("arrangement", "n", "Sigma|q|", "ddE new", "ddE old", "change", "left new", "left old"))
    for arr in sorted(want, key=lambda x: (budget[x], x)):
        dde = ((EB[("1", arr, "TS")][0] - cT1) - (EB[("1", arr, "R")][0] - cR1)) * H
        old = want[arr][1]
        say("   %-22s %3d %7.4f e  %+10.4f %+10.4f %+9.4f   %8.4f %8.4f"
            % (arr, want[arr][0], budget[arr], dde, old, dde - old,
               bar_new + dde, 18.4610 + old))
        rows.append(dict(stage="B1", name=arr, n=want[arr][0], sumq=budget[arr],
                         dde=dde, dde_old=old, left=bar_new + dde))
    say("")
    say("   R34.7. arm 2 minus arm 1: the size of the stability correction.")
    if arm2_ok:
        say("")
        say("   %-22s %4s  %16s %12s" % ("arrangement", "leg", "dE (Eh)", "dE (kcal/mol)"))
        worst = 0.0
        for arr in sorted(list(want) + ["S6-C"], key=lambda x: (budget[x], x)):
            for leg in ("R", "TS"):
                d = EB[("2", arr, leg)][0] - EB[("1", arr, leg)][0]
                worst = max(worst, abs(d))
                say("   %-22s %4s  %+16.9f %+12.6f" % (arr, leg, d, d * H))
        say("")
        say("   largest single correction %.3e Eh (%.4f kcal/mol)" % (worst, worst * H))
        say("")
        say("   R34.8. arm 2's four-leg verdicts, stability on the LAST analysis (AMENDMENT 5).")
        say("")
        say("   %-22s %3s %8s  %-12s %10s  %s" % ("arrangement", "n", "Sigma|q|", "verdict", "ddE", "quotable"))
        for arr in sorted(want, key=lambda x: (budget[x], x)):
            flags, held = [], False
            for leg in ("R", "TS"):
                r, h = legs(B[("2", arr, leg)], B[("2", "S6-C", leg)], "frozen", rep, arr, leg, audit, flags)
                held = held or h
                s = B[("2", arr, leg)]["stab"]
                if not s["ran"]:
                    gates.append("arm 2 %s %s ran no stability analysis" % (arr, leg))
                if s["last"] == "UNSTABLE":
                    flags.append("%s unstable at the last analysis, %s negative" % (leg, s["nneg"]))
                elif s["first"] == "UNSTABLE":
                    say("        note: %s %s first solution unstable (%s negative, lowest %+.6f Eh), corrected"
                        % (arr, leg, s["nneg"], s["lo"]))
            v = "UNRELIABLE" if flags else ("INCONCLUSIVE" if held else "RELIABLE")
            c2R, c2T = EB[("2", "S6-C", "R")][0], EB[("2", "S6-C", "TS")][0]
            dde2 = ((EB[("2", arr, "TS")][0] - c2T) - (EB[("2", arr, "R")][0] - c2R)) * H
            quot = (v == "RELIABLE") and (budget[arr] <= CONV)
            say("   %-22s %3d %7.4f e  %-12s %+10.4f  %s"
                % (arr, want[arr][0], budget[arr], v, dde2,
                   "yes" if quot else ("no, Sigma|q| > %.2f e" % CONV if v == "RELIABLE" else "no, %s" % v)))
            for f in flags:
                say("        flag: %s" % f)
            rows.append(dict(stage="B2", name=arr, n=want[arr][0], sumq=budget[arr],
                             verdict=v, dde=dde2, quotable=quot))
        say("")
        q = [r for r in rows if r["stage"] == "B2" and r.get("quotable")]
        say("   R34.9. quotable corrected barriers: %s"
            % (", ".join("%s %+.4f" % (r["name"], r["dde"]) for r in q) if q else "none"))
        say("")
        say("   R34.10. stage A's three-leg reading against arm 2's four-leg verdict:")
        for arr in sorted(want, key=lambda x: (budget[x], x)):
            sa = next((r["verdict"] for r in rows if r["stage"] == "A" and r["name"] == arr), "-")
            sb = next((r["verdict"] for r in rows if r["stage"] == "B2" and r["name"] == arr), "-")
            agree = ("agree" if (sa == "FLAGGED") == (sb == "UNRELIABLE") else "DISAGREE")
            say("   %-22s stage A %-22s arm 2 %-12s %s" % (arr, sa, sb, agree))
        say("")
        say("          They are different levels of theory, so a disagreement is a finding")
        say("          and not an error.")
    else:
        say("")
        say("   NOT MEASURED. Arm 2 did not produce a readable stability analysis, so R34.7,")
        say("   R34.8, R34.9 and R34.10 are not answered by this diagnostic. R34.9's quotable")
        say("   list is therefore EMPTY: a corrected ddE needs an arm 2 verdict, and there is")
        say("   none. Arm 1's numbers above stand as measurements and are not quotable.")

    # ---- scope, from the pre-registration
    say("")
    say("=" * 78)
    say("WHAT THIS DOES NOT SAY  (D34.1-D34.8, as pre-registered)")
    say("=" * 78)
    say("   Every sphere sits at 2.85-2.96 A, inside D11's r_safe of 3.50 A, and the route")
    say("   decision is open: \"No design measured inside r_safe is interpretable.\"")
    out = [arr for arr in want if budget[arr] > CONV]
    say("   Outside D30's confirmed basis window (Sigma|q| > %.2f e): %s"
        % (CONV, ", ".join("%s %.4f e" % (x, budget[x]) for x in sorted(out, key=lambda z: budget[z]))))
    say("   S6-20 at %.4f e is past the 5.6433 e vector where D30 found def2-SVPD not" % budget["S6-20"])
    say("   sufficient under this functional.")
    say("   These arrangements were placed by the OLD A matrix, which D33 stage A measured")
    say("   to be a different matrix. D34 re-reads the record; it designs nothing.")
    say("   D33 stage B measured that a verdict does not transfer between arrangements.")
    say("   Frozen geometries, one frame, in vacuo apart from the MM Lennard-Jones term.")
    zero = [arr for arr in want if abs(sum(sphere_q(os.path.join(s6, "relax", arr, "R_sp", "site.ORCAFF.prms")))) < 1e-12]
    say("   Sigma qi = 0 exactly: %s. The decision is untaken and D34 does not take it."
        % (", ".join(zero) if zero else "none"))
    say("   No budget is reinstated.")
    say("")
    say("GATES")
    if gates:
        for g in gates:
            say("   FAILED  %s" % g)
    else:
        say("   all gates passed")

    if not arm2_ok:
        Path(a.arm2note).write_text(
            "STAGE 6, DIAGNOSTIC 34 - ARM 2 DID NOT PRODUCE A READABLE RESULT\n"
            "\n"
            "Written by cs6_d34.py on reading the harvest. D34 section 4 fixed in advance that\n"
            "arm 2's outputs are not committed if it fails and that a note about the failure is.\n"
            "This is that note. Arm 1, the lift, is unaffected and its numbers are in\n"
            "STAGE6_DIAGNOSTIC34_RESULT_20261110.txt.\n"
            "\n"
            "A stability-carrying frozen arm needs its own dated amendment.\n"
            "\n" + "\n".join("   " + x for x in arm2_reason) + "\n")
    elif os.path.isfile(a.arm2note):
        os.remove(a.arm2note)
    Path(a.result).write_text("\n".join(rep) + "\n")
    with open(a.tsv, "w") as f:
        f.write("stage\tname\tn\tsumq\tverdict\tdde\tdde_old\tleft\trise\tdip\tquotable\n")
        for r in rows:
            f.write("%s\t%s\t%s\t%.4f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                r["stage"], r["name"], r.get("n", ""), r["sumq"], r.get("verdict", ""),
                "" if r.get("dde") is None else "%.6f" % r["dde"],
                "" if r.get("dde_old") is None else "%.6f" % r["dde_old"],
                "" if r.get("left") is None else "%.6f" % r["left"],
                "" if r.get("rise") is None else "%+.4f" % r["rise"],
                "" if r.get("dip") is None else "%.4f" % r["dip"],
                "" if "quotable" not in r else ("yes" if r["quotable"] else "no")))
    print("\n".join(rep))
    if gates:
        stop("%d gate(s) failed; the result records them." % len(gates))


if __name__ == "__main__":
    main()
