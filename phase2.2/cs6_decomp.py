#!/usr/bin/env python3
"""Stage 6 DIAGNOSTIC 5: decompose Diagnostic 4's sixteen committed single points.

NO NEW CALCULATION. Everything is read back out of the committed, sha256-manifested archives:
field/outputs/field_outputs.tar.xz for the sixteen runs and qscan/outputs/qscan_outputs.tar.xz for the bare
references, the same two files Diagnostic 4 itself measured from.

Diagnostic 4 reported ddE = (E_TS - E_TS,bare) - (E_R - E_R,bare) and nothing finer, so it could not say whether
def2-SVPD's turnover above lambda = 0.5 lives in the reactant, in the transition state, or in both. The sixteen
energies already on disk answer that, together with the orbital and dipole diagnostics ORCA printed in the same
files. The rule is STAGE6_DIAGNOSTIC5_20261026.txt D5.3, fixed before any of it was read.

usage: cs6_decomp.py <results-repo> <out-dir>
"""
import sys, re, tarfile, hashlib
from pathlib import Path
import numpy as np

H = 627.5094740631
EV = 27.211386245988
LAMBDAS = (0.25, 0.50, 0.75, 1.00)
BASES = (("def2-SVP", "C1"), ("def2-SVPD", "C5"))
NEL = 118
D4_DDE = {("def2-SVP", 0.25): -5.0969, ("def2-SVP", 0.50): -10.6573,
          ("def2-SVP", 0.75): -17.0581, ("def2-SVP", 1.00): -24.6898,
          ("def2-SVPD", 0.25): -5.5601, ("def2-SVPD", 0.50): -11.6329,
          ("def2-SVPD", 0.75): -10.2346, ("def2-SVPD", 1.00): -6.1422}
TIE_TOL = 1e-3          # kcal/mol: how closely this script must reproduce the committed D4 table
SMOOTH = 0.3            # kcal/mol: D5.3(b) leave-one-out error floor
OUT = 3.0               # D5.3(b): how much worse than the others a held-out point must be
AGREE = 1.0             # kcal/mol: D5.3(d) inter-basis agreement threshold
GAP_EV = 1.0            # eV:       D5.3(c) gap floor
POW = 2.5               # D5.3(c): lambda exponent above which the shift is steeper than polarisation
DIPR = 1.5              # D5.3(c): def2-SVPD dipole growth relative to def2-SVP
HOMO_EV = 0.5           # eV: D5.3(c) inter-basis HOMO departure
CANCEL = 3.0            # D5.3(a): state shift over Delta above which ddE is a cancellation


def stop(m):
    print("STOP: %s" % m, file=sys.stderr)
    raise SystemExit(1)


def manifest(p):
    M = {}
    for l in Path(p).read_text().splitlines():
        if l.strip():
            h, n = l.split(None, 1)
            M[n.strip()] = h
    return M


def grab(tf, name, man):
    """Pull one member out of a committed archive, checking it against that archive's own manifest."""
    key = name if name.startswith("./") else "./" + name.lstrip("./")
    try:
        raw = tf.extractfile(key).read()
    except KeyError:
        stop("%s is not in the archive" % key)
    if key not in man:
        stop("%s is not listed in the archive's manifest" % key)
    if hashlib.sha256(raw).hexdigest() != man[key]:
        stop("%s does not match its sha256 in the archive's manifest" % key)
    return raw.decode(errors="replace")


def parse(txt, who, npc):
    """Energy plus the diagnostics ORCA printed alongside it. Every field is required, never defaulted."""
    if "****ORCA TERMINATED NORMALLY****" not in txt:
        stop("%s: ORCA did not terminate normally" % who)
    if not re.search(r"SCF CONVERGED AFTER\s+\d+\s+CYCLES", txt) or "SCF NOT CONVERGED" in txt:
        stop("%s: the SCF did not converge" % who)
    m = re.search(r"Number of Electrons\s+NEL\s+\.+\s+(\d+)", txt)
    if not m or int(m.group(1)) != NEL:
        stop("%s: %s electrons, expected %d" % (who, m.group(1) if m else "no count", NEL))
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", txt)
    if len(e) != 1:
        stop("%s: %d final energies, expected 1" % (who, len(e)))
    if npc is None:
        if re.search(r"Reading point charge file", txt):
            stop("%s: a bare reference must carry no point-charge file" % who)
    elif not re.search(r"Reading point charge file\s*\.*\s*ok \(%d point charges\)" % npc, txt):
        stop("%s: ORCA did not confirm reading exactly %d point charges" % (who, npc))
    nb = re.search(r"Number of basis functions\s+\.+\s+(\d+)", txt)
    if not nb:
        stop("%s: no basis-function count" % who)
    # orbital energies: HOMO is the last occupied, LUMO the first empty
    k = txt.find("ORBITAL ENERGIES")
    if k < 0:
        stop("%s: no ORBITAL ENERGIES block" % who)
    occ, vir = [], []
    for l in txt[k:].splitlines()[4:]:
        f = l.split()
        if len(f) != 4:
            if occ and vir:
                break
            continue
        try:
            o, eh = float(f[1]), float(f[2])
        except ValueError:
            break
        (occ if o > 0.5 else vir).append(eh)
    if not occ or not vir:
        stop("%s: could not read occupied and virtual orbitals (%d/%d)" % (who, len(occ), len(vir)))
    if len(occ) != NEL // 2:
        stop("%s: %d occupied orbitals, expected %d" % (who, len(occ), NEL // 2))
    homo, lumo = occ[-1], vir[0]
    d = re.search(r"Magnitude \(Debye\)\s+:\s+(\d+\.\d+)", txt)
    if not d:
        stop("%s: no dipole magnitude" % who)
    lw = txt.find("LOEWDIN ATOMIC CHARGES")
    if lw < 0:
        stop("%s: no LOEWDIN ATOMIC CHARGES block" % who)
    Lq = []
    for l in txt[lw:].splitlines()[2:]:
        mm = re.match(r"\s*\d+\s+[A-Za-z]{1,2}\s*:\s*(-?\d+\.\d+)", l)
        if not mm:
            break
        Lq.append(float(mm.group(1)))
    if not Lq:
        stop("%s: could not read Loewdin charges" % who)
    return dict(E=float(e[0]), nb=int(nb.group(1)), homo=homo * EV, lumo=lumo * EV,
                gap=(lumo - homo) * EV, dip=float(d.group(1)), lq=np.array(Lq))


def main(root, out):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    fo, qo = T6 / "field" / "outputs", T6 / "qscan" / "outputs"
    for p in (fo / "field_outputs.tar.xz", fo / "field_outputs.MANIFEST.sha256",
              qo / "qscan_outputs.tar.xz", qo / "qscan_outputs.MANIFEST.sha256",
              T6 / "field" / "DIAGNOSTIC4_RESULT.txt"):
        if not p.is_file():
            stop("%s is missing - run the x97z bundle first" % p)
    FM, QM = manifest(fo / "field_outputs.MANIFEST.sha256"), manifest(qo / "qscan_outputs.MANIFEST.sha256")

    D, Bare = {}, {}
    with tarfile.open(str(qo / "qscan_outputs.tar.xz")) as t:
        for basis, cond in BASES:
            for st in ("R", "TS"):
                Bare[(basis, st)] = parse(grab(t, "./%s/bare_%s/job.out" % (cond, st), QM),
                                          "%s bare %s" % (basis, st), None)
    with tarfile.open(str(fo / "field_outputs.tar.xz")) as t:
        for basis, _ in BASES:
            for lam in LAMBDAS:
                for st in ("R", "TS"):
                    nm = "./%s/lam%03d_%s/job.out" % (basis, round(lam * 100), st)
                    D[(basis, lam, st)] = parse(grab(t, nm, FM), nm, 20)
    print("   16 field runs and 4 bare references read back out of the committed archives,")
    print("   every one checked against its archive's own sha256 manifest")

    # the tie to the published result: this script must reproduce D4's table or nothing below is readable
    dE, ddE = {}, {}
    for basis, _ in BASES:
        for lam in LAMBDAS:
            for st in ("R", "TS"):
                dE[(basis, lam, st)] = (D[(basis, lam, st)]["E"] - Bare[(basis, st)]["E"]) * H
            ddE[(basis, lam)] = dE[(basis, lam, "TS")] - dE[(basis, lam, "R")]
            w = D4_DDE[(basis, lam)]
            if abs(ddE[(basis, lam)] - w) > TIE_TOL:
                stop("%s lambda %.2f: this script gets ddE %+.4f, the committed DIAGNOSTIC4_RESULT.txt table "
                     "says %+.4f. %.4f kcal/mol apart, over the %.0e tolerance - the decomposition does not "
                     "belong to the published result and must be reconciled before anything is read"
                     % (basis, lam, ddE[(basis, lam)], w, abs(ddE[(basis, lam)] - w), TIE_TOL))
    print("   all 8 reconstructed ddE values reproduce the committed DIAGNOSTIC4_RESULT.txt table to < %.0e"
          % TIE_TOL)

    O = ["STAGE 6 - DIAGNOSTIC 5 RESULT: where Diagnostic 4's divergence lives", "",
         "Measured by cs6_decomp.py from the committed field and qscan output archives. NO NEW CALCULATION:",
         "the same sixteen single points Diagnostic 4 measured, decomposed into their two states, with the",
         "orbital and dipole diagnostics ORCA printed in the same files. The rule is",
         "STAGE6_DIAGNOSTIC5_20261026.txt D5.3, fixed before any of this was read.", "",
         "   dE_R  = E(R  + charges) - E(R  bare)        dE_TS = E(TS + charges) - E(TS bare)",
         "   ddE   = dE_TS - dE_R, which is what Diagnostic 4 reported", "",
         "      basis      lambda     dE_R       dE_TS        ddE      (kcal/mol)"]
    for basis, _ in BASES:
        for lam in LAMBDAS:
            O.append("      %-10s  %4.2f  %10.4f %10.4f %10.4f"
                     % (basis, lam, dE[(basis, lam, "R")], dE[(basis, lam, "TS")], ddE[(basis, lam)]))
    dR = {l: dE[("def2-SVPD", l, "R")] - dE[("def2-SVP", l, "R")] for l in LAMBDAS}
    dT = {l: dE[("def2-SVPD", l, "TS")] - dE[("def2-SVP", l, "TS")] for l in LAMBDAS}
    O += ["", "   the inter-basis shift of each state, def2-SVPD minus def2-SVP", "",
          "      lambda   delta dE_R   delta dE_TS   their difference (= Delta of Diagnostic 4)"]
    for lam in LAMBDAS:
        O.append("      %4.2f  %11.4f %13.4f %18.4f" % (lam, dR[lam], dT[lam], dT[lam] - dR[lam]))

    O += ["", "THE FOUR READINGS OF D5.3", "", "   (a) WHERE THE DIVERGENCE LIVES", "",
          "      delta dE_TS minus delta dE_R is identically Diagnostic 4's Delta, so these two are not",
          "      independent measurements: this reading asks how large the COMMON shift is beside Delta.", ""]
    aR, aT, aD = abs(dR[1.00]), abs(dT[1.00]), abs(dT[1.00] - dR[1.00])
    lo = min(aR, aT)
    O.append("      at lambda = 1.00 the reactant shifts %+.4f and the transition state %+.4f kcal/mol, against"
             % (dR[1.00], dT[1.00]))
    O.append("      a Delta of %+.4f" % (dT[1.00] - dR[1.00]))
    if lo >= CANCEL * max(aD, 1e-9):
        O += ["      CANCELLATION-DOMINATED: the smaller state shift is %.1f times Delta itself, so both states"
              % (lo / max(aD, 1e-9)),
              "      move far more than their difference does and ddE at def2-SVPD is a small difference of",
              "      large numbers. Neither basis can be quoted at this field, and what has to change is the",
              "      charge model rather than the basis."]
    else:
        O += ["      NOT A CANCELLATION: the smaller state shift is %.2f times Delta, under the factor %.0f this"
              % (lo / max(aD, 1e-9), CANCEL),
              "      reading uses, so the basis change does not mostly cancel in the difference and the",
              "      turnover is a real differential response of the two states.",
              "      %s carries the larger shift, which is where to look."
              % ("The reactant" if aR > aT else "The transition state")]
    O += ["", "   (b) HOW NON-QUADRATIC EACH BASIS IS", "",
          "      Four lambda values against a quadratic leave ONE degree of freedom. That is too little to say",
          "      WHICH point is an outlier: holding each lambda out in turn and predicting it from the other",
          "      three gives errors locked in a fixed ratio of about 3 to 1 between the endpoints and the",
          "      interior, set by the spacing alone, and a single bad point perturbs all four in that same",
          "      ratio. So this reading does NOT claim to localise an outlier, and it is not written to. A",
          "      fifth lambda - 0.875 is the natural one - would be needed for that, and that is recorded as",
          "      an open item rather than asserted here.",
          "      What four points DO settle is how far each basis departs from a quadratic at all, which is a",
          "      like-for-like comparison between the two bases on identical geometry and design.", ""]
    lamA = np.array(LAMBDAS)
    qr = {}
    for basis, _ in BASES:
        for st, lab in (("R", "reactant"), ("TS", "transition state")):
            y = np.array([dE[(basis, l, st)] for l in LAMBDAS])
            c = np.polyfit(lamA, y, 2)
            qr[(basis, st)] = float(np.abs(np.polyval(c, lamA) - y).max())
            O.append("      %-10s %-18s max quadratic residual %8.4f kcal/mol" % (basis, lab, qr[(basis, st)]))
    loo = []
    for st in ("R", "TS"):
        y = np.array([dE[("def2-SVPD", l, st)] for l in LAMBDAS])
        e = [abs(float(np.polyval(np.polyfit(lamA[[j for j in range(4) if j != i]],
                                             y[[j for j in range(4) if j != i]], 2), lamA[i])) - y[i])
             for i in range(4)]
        loo.append("      %-2s leave-one-out errors  %s  kcal/mol (reported as data, not as a test)"
                   % (st, "  ".join("%.3f" % v for v in e)))
    O += [""] + loo + [""]
    svpq = max(qr[("def2-SVP", s)] for s in ("R", "TS"))
    svpdq = max(qr[("def2-SVPD", s)] for s in ("R", "TS"))
    if svpq <= SMOOTH and svpdq > SMOOTH:
        O += ["      THE DIFFUSE BASIS INTRODUCES A RESPONSE THE SMALL BASIS DOES NOT HAVE AT ALL: def2-SVP is",
              "      quadratic in lambda to %.4f kcal/mol while def2-SVPD departs by %.4f, a factor %.0f. The"
              % (svpq, svpdq, svpdq / max(svpq, 1e-9)),
              "      higher-order response is therefore carried entirely by the diffuse functions, which is",
              "      what DIAGNOSTIC3_NOTE1 N3.2 measured at a single site and is reproduced here for twenty."]
    elif svpdq <= SMOOTH and svpq <= SMOOTH:
        O += ["      BOTH BASES ARE QUADRATIC to better than %.1f kcal/mol, so neither carries a higher-order"
              % SMOOTH,
              "      response and Diagnostic 4's turnover is not a departure from the quadratic model in either",
              "      state separately. It would then live wholly in the difference, which this cannot explain."]
    else:
        O += ["      BOTH BASES DEPART FROM A QUADRATIC (def2-SVP by %.4f, def2-SVPD by %.4f kcal/mol), so the"
              % (svpq, svpdq),
              "      higher-order response is not specific to the diffuse functions and the basis is not on its",
              "      own the explanation."]
    O += ["", "   (c) SPILL-OUT SIGNATURES", "",
          "      Reported first because it bears on the withdrawn Finding 8: the HOMO of this dianion is",
          "      POSITIVE in vacuo at every basis tested (def2-SVP +2.278, ma-def2-SVP +1.227, def2-SVPD +1.231",
          "      eV for the bare reactant), which is DFT self-interaction and not proof of an unbound electron -",
          "      the VDE of +0.5123 eV at two diffuse bases is the test that settled it.", "",
          "      basis      lambda  state     HOMO (eV)   gap (eV)   dipole (D)   max |Loewdin| shift"]
    for basis, _ in BASES:
        for lam in LAMBDAS:
            for st in ("R", "TS"):
                d = D[(basis, lam, st)]
                sh = float(np.abs(d["lq"] - D[("def2-SVP", lam, st)]["lq"]).max()) if basis == "def2-SVPD" else 0.0
                O.append("      %-10s  %4.2f   %-4s %11.4f %10.4f %12.4f %18.4f"
                         % (basis, lam, st, d["homo"], d["gap"], d["dip"], sh))
    # The decisive signature is the POWER with which the diffuse basis's extra stabilisation of each state grows
    # in lambda. Ordinary polarisation is second order in the field, so n ~ 2; anything steeper is response the
    # small basis cannot represent at all, which is what spill-out into diffuse tails looks like. A HOMO that
    # merely rises is NOT the test: these runs start from a POSITIVE HOMO in vacuo (DFT self-interaction, the
    # reason Finding 8 was withdrawn) and the design's mostly positive charges push orbital energies DOWN, so a
    # rise would not be expected even if density were leaking.
    sig = []
    lg = np.log(np.array(LAMBDAS))
    com = np.array([0.5 * (dR[l] + dT[l]) for l in LAMBDAS])
    if np.all(np.abs(com) > 1e-6):
        n = float(np.polyfit(lg, np.log(np.abs(com)), 1)[0])
        O.append("      the COMMON inter-basis shift, (delta dE_R + delta dE_TS)/2, grows as lambda^%.2f" % n)
        O.append("      (Delta's own shape is already known from Diagnostic 4, so the common part is the new")
        O.append("      information here. Ordinary polarisation is second order in the field, so n near 2.)")
        if n > POW:
            sig.append("the common inter-basis shift grows as lambda^%.2f, steeper than the lambda^2 of "
                       "ordinary polarisation" % n)
    for st in ("R", "TS"):
        gs = [D[("def2-SVPD", l, st)]["gap"] for l in LAMBDAS]
        if min(gs) < GAP_EV:
            sig.append("the %s gap falls to %.3f eV, below %.1f eV" % (st, min(gs), GAP_EV))
        ds = [D[("def2-SVPD", l, st)]["dip"] for l in LAMBDAS]
        d0 = [D[("def2-SVP", l, st)]["dip"] for l in LAMBDAS]
        g1 = (ds[-1] / ds[0]) / ((d0[-1] / d0[0]) or 1.0)
        O.append("      the %-2s dipole grows %.2f times as fast at def2-SVPD as at def2-SVP, gap floor %.3f eV"
                 % (st, g1, min(gs)))
        if g1 > DIPR:
            sig.append("the %s dipole grows %.2f times faster at def2-SVPD than at def2-SVP" % (st, g1))
        hb = [D[("def2-SVPD", l, st)]["homo"] - D[("def2-SVP", l, st)]["homo"] for l in LAMBDAS]
        O.append("      the %-2s HOMO sits %+.3f eV from def2-SVP at lambda 0.25 and %+.3f eV at lambda 1.00"
                 % (st, hb[0], hb[-1]))
        if abs(hb[-1]) > HOMO_EV and abs(hb[-1]) > 2 * abs(hb[0]):
            sig.append("the %s HOMO pulls %+.3f eV away from def2-SVP by lambda = 1, from %+.3f at lambda 0.25"
                       % (st, hb[-1], hb[0]))
    O.append("")
    if sig:
        O.append("      SPILL-OUT IS SUPPORTED. " + "; ".join(sig) + ".")
        O += ["      Taken with Finding 4 - our charge sites have a wall for nuclei and none for electrons - and",
              "      with Marefat Khah et al. JCTC 2020 doi:10.1021/acs.jctc.9b01162, which reports unphysical",
              "      leaking of the electron density in point-charge embedding 'in particular, when diffuse basis",
              "      functions are required', the def2-SVPD high-field points are NOT usable as physical",
              "      barriers. The bound Diagnostic 4 measured would then be a bound on the model, not on the",
              "      chemistry, and a charge model with a Pauli wall for electrons is needed before any",
              "      high-field number means anything."]
    else:
        O += ["      SPILL-OUT IS NOT SUPPORTED by these signatures: the common inter-basis shift is no steeper",
              "      than lambda^%.1f, no gap falls below %.1f eV, no dipole runs away at def2-SVPD relative to"
              % (POW, GAP_EV),
              "      def2-SVP, and no HOMO pulls away from it with field. Real differential",
              "      saturation then becomes the leading explanation of Diagnostic 4's turnover, and the",
              "      total-field bound it measured should be read as chemistry."]

    O += ["", "   (d) THE RANGE IN WHICH THE TWO BASES AGREE", ""]
    ok = [l for l in LAMBDAS if abs(ddE[("def2-SVPD", l)] - ddE[("def2-SVP", l)]) <= AGREE]
    SABS = 5.760
    if ok:
        O += ["      the two bases agree to within %.1f kcal/mol up to lambda = %.2f, i.e. sum |q| = %.3f e"
              % (AGREE, max(ok), max(ok) * SABS),
              "      (|Delta| = %.4f there). That is the honest domain of validity of the committed def2-SVP"
              % abs(ddE[("def2-SVPD", max(ok))] - ddE[("def2-SVP", max(ok))]),
              "      numbers, and it holds whichever of the two explanations above turns out to be right."]
    else:
        O += ["      the two bases do not agree to within %.1f kcal/mol at any lambda tested, the smallest"
              % AGREE,
              "      being %.4f at lambda = %.2f. There is then no field range in which the committed def2-SVP"
              % (min(abs(ddE[("def2-SVPD", l)] - ddE[("def2-SVP", l)]) for l in LAMBDAS), LAMBDAS[0]),
              "      numbers are basis-converged, and that must be said of all of them."]

    O += ["", "WHAT THIS DOES NOT TEST", "",
          "   Nothing new was computed, so this cannot test a hypothesis the sixteen existing points do not",
          "   already contain. It does not test a charge model with a Pauli wall for electrons, which is the",
          "   fix if spill-out is supported. It does not test def2-TZVPD, so it cannot say the def2-SVPD",
          "   numbers are themselves converged. One design, one frame. Nothing about Loss.", ""]
    Path(out).mkdir(parents=True, exist_ok=True)
    (Path(out) / "DIAGNOSTIC5_RESULT.txt").write_text("\n".join(O) + "\n")
    print()
    print("\n".join(O))
    print("   wrote %s" % (Path(out) / "DIAGNOSTIC5_RESULT.txt"))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        stop("usage: cs6_decomp.py <results-repo> <out-dir>")
    main(sys.argv[1], sys.argv[2])
