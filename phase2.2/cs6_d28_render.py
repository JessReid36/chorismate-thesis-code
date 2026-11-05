#!/usr/bin/env python3
"""Render DIAGNOSTIC 28's result document from the reader's values. Every number here comes from
d28_values.json and every one is independently re-derived from the outputs by verify.py."""
import argparse, json, sys
from pathlib import Path
FUNCS = ["B3LYP", "CAM-B3LYP", "wB97X-D3"]


def main():
    ap = argparse.ArgumentParser()
    for k in ("values", "out", "date", "prereg-commit", "reader-sha"): ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    V = json.loads(Path(a.values).read_text())
    H = V["H"]; pr, mech, sec = V["primary"], V["mechanism"], V["secondary"]
    L, w = [], None
    w = L.append
    def h(f, b="def2-SVPD", s="R"): return H.get("%s|%s|%s" % (f, b, s))

    w("STAGE 6 - DIAGNOSTIC 28, RESULT%s%s" % (" " * 18, a.date)); w("")
    if pr == "CONFIRMED":
        w("B3LYP WAS PLACING THE HIGHEST OCCUPIED ORBITAL TOO HIGH. WITH FULL LONG-RANGE EXACT")
        w("EXCHANGE IT IS BOUND, AS THE MEASURED DETACHMENT ENERGY SAYS IT SHOULD BE, SO THE FIELD")
        w("WINDOW DIAGNOSTIC 27 IMPLIES IS TOO PESSIMISTIC.")
    elif pr == "PARTIAL":
        w("B3LYP IS PLACING THE ORBITAL TOO HIGH, BUT CORRECTING THE LONG-RANGE EXCHANGE DOES NOT BY")
        w("ITSELF MAKE IT BOUND. THE ERROR IS REAL AND IS NOT THE WHOLE STORY.")
    elif pr == "NOT CONFIRMED":
        w("B3LYP'S PLACEMENT OF THE ORBITAL IS NOT THE EXPLANATION. DIAGNOSTIC 27'S GAP NARROWING")
        w("STANDS AS MEASURED, AND THE FUNCTIONAL COMES OFF THE LIST OF THINGS TO TRY.")
    else:
        w("NO PRIMARY READING IS AVAILABLE. THE RUNS THAT WOULD CARRY IT ENGAGED TRAH, SO THEIR")
        w("ORBITAL ENERGIES ARE LEVEL-SHIFTED AND NOT QUOTABLE.")
    w("")
    w("Pre-registered in STAGE6_DIAGNOSTIC28_20261105.txt, which is NOT edited and was committed")
    w("before any job was submitted (results %s). STAGE6_CRITERIA.txt is unchanged and the" % a.prereg_commit)
    w("production functional is NOT changed by this document. NOTHING IS DECIDED HERE.")
    w("")
    w("   reader   phase2.2/cs6_d28.py, sha256 %s" % a.reader_sha)
    w("   runs     22 ORCA jobs, phase2.2/cs_stage6/d28/, one PBS job")
    w("")
    w("")
    w("R28.1  THE PRIMARY READING")
    w("")
    w("Highest occupied orbital of the BARE substrate - no charges - at def2-SVPD, reactant. A")
    w("negative value means the orbital lies below the vacuum level, i.e. the electron is bound.")
    w("")
    w("   %-11s %-26s %s" % ("functional", "long-range exact exchange", "HOMO"))
    for f in FUNCS:
        v = h(f)
        w("   %-11s %-26.2f %s" % (f, V["LRX"][f], ("%+9.4f eV" % v) if v is not None else "NOT QUOTABLE (level-shifted)"))
    w("")
    w("   Probe 3's measured detachment energy is %+.4f eV, so for the exact functional this" % V["VDE"])
    w("   orbital should sit near %+.4f eV." % -V["VDE"])
    if V["drop"] is not None:
        w("   wB97X-D3 moves it by %+.4f eV from B3LYP's value." % -V["drop"])
    w("")
    w("VERDICT: %s" % pr)
    w("")
    w("")
    w("R28.2  THE MECHANISM CHECK")
    w("")
    w("D28.4 asked separately whether the orbital falls monotonically with the long-range exact")
    w("exchange fraction - 0.20, 0.65, 1.00. One functional moving proves much less than three in")
    w("order.")
    w("")
    w("   %-11s %-4s %11s %11s %11s   %s" % ("basis", "st", "B3LYP", "CAM-B3LYP", "wB97X-D3", "monotone?"))
    for r in V["mech_rows"]:
        cells = "".join(("%11.4f" % x) if x is not None else "         NA" for x in r["values"])
        w("   %-11s %-4s %s   %s" % (r["basis"], r["state"], cells,
          {True: "yes", False: "NO", None: "incomplete"}[r["monotone"]]))
    w("")
    w("VERDICT: %s" % mech)
    if mech == "NOT MONOTONE":
        w("")
        w("The three are not in order of long-range exact exchange, so that fraction is not the")
        w("mechanism, and the primary verdict above stands WITHOUT a mechanism behind it.")
    w("")
    w("")
    w("R28.3  THE SECONDARY READING")
    w("")
    w("The gap at sum|q| = 2.8800 e, def2-SVPD, reactant, against TOL = 1.0 eV:")
    w("")
    for f in FUNCS:
        v = V["gap_2880"].get(f)
        w("   %-11s %s" % (f, ("%7.4f eV" % v) if v is not None else "NOT QUOTABLE"))
    w("")
    w("VERDICT: %s" % sec)
    w("")
    w("")
    w("R28.4  WHAT ENGAGED TRAH, AND WHAT THAT WITHHELD")
    w("")
    w("This is recorded because DIAGNOSTIC 27 did not record it and a claim had to be retracted.")
    w("In a run where ORCA engages its TRAH converger, the printed orbital energies are")
    w("level-shifted and any negative-gap warning precedes convergence rather than describing the")
    w("converged state. No orbital energy from such a run is quoted anywhere above.")
    w("")
    w("   runs whose orbital energies were withheld : %d" % len(V["shifted"]))
    for s in V["shifted"]: w("      %s" % s.replace("|", " / "))
    w("   runs emitting a negative-gap warning before convergence : %d" % len(V["warn_before_conv"]))
    for s in V["warn_before_conv"]: w("      %s" % s.replace("|", " / "))
    w("")
    w("")
    w("R28.5  WHAT THIS DOES NOT SETTLE")
    w("")
    w("1. It does not establish which functional is correct, and does not change the production")
    w("   functional. That needs an amendment.")
    w("2. It compares orbital energies across functionals, which is weaker than within one. The")
    w("   detachment-energy comparison is the anchor because it is an energy difference.")
    w("3. It says nothing about solvation, which addresses the same headroom differently.")
    w("4. It did not run the strongest field at def2-TZVPD, so it cannot say whether the basis")
    w("   sequence converges under a different functional.")
    w("5. No barrier is quoted and no design is accepted or rejected.")
    w("6. It cannot recover a converged gap for any run that engaged TRAH. That needs a run with")
    w("   the level shift controlled, which is not this diagnostic.")
    if V["gates_failed"]:
        w(""); w("GATES THAT DID NOT PASS, committed as failures:")
        for x in V["gates_failed"]: w("   %s" % x)
    w("")
    Path(a.out).write_text("\n".join(L) + "\n")
    print("rendered %s  (%d lines)" % (a.out, len(L)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
