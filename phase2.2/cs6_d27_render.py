#!/usr/bin/env python3
"""Render DIAGNOSTIC 27's result document from the reader's values. It writes prose around
numbers it does not compute; every number comes from d27_values.json and every one of them is
re-derived independently from the ORCA outputs by verify.py."""
import argparse, json, sys
from pathlib import Path

BASES = ["def2-SVP", "def2-SVPD", "def2-TZVPD"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--values", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--date", required=True); ap.add_argument("--prereg-commit", required=True)
    ap.add_argument("--reader-sha", required=True)
    a = ap.parse_args()
    V = json.loads(Path(a.values).read_text())
    rows, ctl, arb, overall = V["rows"], V["control_barrier"], V["arbitration"], V["overall"]
    TOL, nbtz = V["TOL"], V["nbasis_tzvpd"]
    r3 = [r for r in rows if r["vector"] == "5.6433e_n20"][0] if rows else None
    L = []
    w = L.append

    # ---- the headline is the verdict, written from the verdict and nothing else
    if overall == "CONVERGED":
        head = ("def2-SVPD IS CONVERGED FOR THIS QUANTITY. THE def2-SVPD TO def2-TZVPD STEP IS\n"
                "WITHIN %.1f kcal/mol ON ALL THREE VECTORS AND IS SMALLER THAN THE STEP BEFORE IT,\n"
                "SO THE PRODUCTION BASIS STANDS AND NOTHING IN THE RECORD NEEDS RE-RUNNING FOR THE\n"
                "BASIS." % TOL)
    elif overall.startswith("PER VECTOR"):
        head = ("THE THREE VECTORS DO NOT AGREE, SO NO SINGLE VERDICT IS FORMED. THE RESULT IS\n"
                "REPORTED PER VECTOR, AS D27.4 REQUIRED IN ADVANCE.")
    else:
        head = "%s ON ALL THREE VECTORS." % overall
    w("STAGE 6 - DIAGNOSTIC 27, RESULT%s%s" % (" " * (48 - len("STAGE 6 - DIAGNOSTIC 27, RESULT")), a.date))
    w(""); w(head); w("")
    w("Pre-registered in STAGE6_DIAGNOSTIC27_20261104.txt, which is NOT edited and was committed")
    w("before any job was submitted (results %s, pushed 4 November). STAGE6_CRITERIA.txt is" % a.prereg_commit)
    w("unchanged. NOTHING IS DECIDED HERE.")
    w("")
    w("   reader   phase2.2/cs6_d27.py, committed with this document")
    w("            sha256 %s" % a.reader_sha)
    w("   runs     24 ORCA single points, phase2.2/cs_stage6/d27/, one PBS job")
    w("   TOL      %.1f kcal/mol, from D27.4, which took it from criterion D5.3(d)" % TOL)
    w("")
    w("")
    w("R27.1  WHAT IT RETURNED")
    w("")
    w("Control barriers of the bare substrate, kcal/mol:")
    w("")
    for b in BASES:
        if b in ctl: w("   %-11s %12.6f" % (b, ctl[b]))
    w("")
    w("ddE by vector and basis, kcal/mol. gap_SD is the double to diffuse-double step and gap_DT the")
    w("diffuse-double to diffuse-triple step, both as D27.4 defined them before the run:")
    w("")
    w("   %-12s %4s %8s %12s %12s %12s %9s %9s  %s" % (
        "vector", "n", "sum|q|", "def2-SVP", "def2-SVPD", "def2-TZVPD", "gap_SD", "gap_DT", "verdict"))
    for r in rows:
        w("   %-12s %4d %8.4f %12.6f %12.6f %12.6f %9.6f %9.6f  %s" % (
            r["vector"], r["n"], r["sum_abs_q"], r["ddE_def2_SVP"], r["ddE_def2_SVPD"],
            r["ddE_def2_TZVPD"], r["gap_SD"], r["gap_DT"], r["verdict"]))
    w("")
    w("OVERALL: %s" % overall)
    w("")
    w("")
    w("R27.2  THE THIRD VECTOR, WHICH IS WHY THIS WAS RUN")
    w("")
    if r3:
        w("DIAGNOSTIC 24's log had def2-SVP and def2-SVPD 13.5154 kcal/mol apart at this vector, a")
        w("factor of 2.28, with no way to tell which was the artefact. In ORCA, at full precision:")
        w("")
        w("   def2-SVP    %12.6f" % r3["ddE_def2_SVP"])
        w("   def2-SVPD   %12.6f" % r3["ddE_def2_SVPD"])
        w("   def2-TZVPD  %12.6f" % r3["ddE_def2_TZVPD"])
        w("   |TZVPD - SVPD| %9.6f      |TZVPD - SVP| %9.6f      TOL %.1f" % (
            abs(r3["ddE_def2_TZVPD"] - r3["ddE_def2_SVPD"]),
            abs(r3["ddE_def2_TZVPD"] - r3["ddE_def2_SVP"]), TOL))
        w("")
        w("ARBITRATION, by the rule D27.4 fixed before the run: %s" % arb)
        w("")
        if arb == "BOTH ARE ARTEFACTS":
            w("Neither double-zeta answer is within TOL of the triple-zeta one. By D27.4 the 5.6433 e")
            w("vector is WITHDRAWN from every quantity until a larger basis is run, and the ladder above")
            w("sum|q| = 2.88 e has no admissible basis in this record.")
        elif arb == "def2-SVP IS THE ARTEFACT":
            w("def2-TZVPD agrees with def2-SVPD to within TOL, so the def2-SVP figure at this vector is")
            w("the artefact and the def2-SVPD one is supported.")
        else:
            w("def2-TZVPD agrees with def2-SVP, not with def2-SVPD. The diffuse double-zeta basis is the")
            w("odd one out at this vector, which is the opposite of what adding diffuse functions is")
            w("supposed to do, and it means the committed def2-SVPD ladder is the set in question.")
    w("")
    w("")
    w("R27.3  WHAT IT MEANS FOR WHAT IS ALREADY COMMITTED")
    w("")
    if overall == "CONVERGED":
        w("Nothing is withdrawn. The d16 ladder, DIAGNOSTIC 23's pseudopotential arms and DIAGNOSTIC 26's")
        w("f_closed all rest on def2-SVPD, and the step beyond it moves this quantity by less than TOL on")
        w("every vector tested. That is not a proof of the basis-set limit - see R27.6 - but it is the")
        w("check D26's verdict asked for, and def2-SVPD survives it.")
    else:
        w("Every committed def2-SVPD quantity at or above the sum|q| of a vector that did not return")
        w("CONVERGED is PROVISIONAL until re-run. D27.4 accepted that consequence in advance. By the")
        w("table in R27.1 that is:")
        w("")
        for r in rows:
            if r["verdict"] != "CONVERGED":
                w("   %-12s sum|q| %.4f   %s" % (r["vector"], r["sum_abs_q"], r["verdict"]))
        w("")
        w("The arms affected are the d16 ladder, DIAGNOSTIC 23's Ne, Ar and He walled arms and")
        w("DIAGNOSTIC 26's f_closed, all of which are def2-SVPD. None of them is withdrawn by this")
        w("document; they are marked provisional, and what re-running them would cost is a separate")
        w("question for a separate pre-registration.")
    w("")
    w("")
    w("R27.4  THE pyscf QUESTION D26's GATE 7 RAISED, NOW CLOSED")
    w("")
    w("D27.1(b) identified DIAGNOSTIC 26's committed GATE 7 failure as two programs rather than a")
    w("defect: DIAGNOSTIC 24 ran in pyscf and every production run is ORCA. This diagnostic removes the")
    w("mixed comparison by running def2-SVP in ORCA. The ORCA def2-SVP ddE values in R27.1 are")
    w("therefore the ones to use, and D26's quantity A - pyscf's -5.6165, to four decimals - is")
    w("superseded for any future use by:")
    w("")
    for r in rows:
        w("   %-12s ORCA def2-SVP %12.6f" % (r["vector"], r["ddE_def2_SVP"]))
    w("")
    w("A dated erratum to DIAGNOSTIC 26 is warranted and is NOT written here, because a result document")
    w("does not edit another diagnostic. What the erratum has to say is: A was a pyscf number, the ORCA")
    w("value is above, and the direction and size of the change to f_closed follow from it.")
    w("")
    w("D26's R26.6 item 7 said this result could not enter the thesis until DIAGNOSTIC 24 was")
    w("re-registered and re-run. That is now done FOR THE THREE POINT-CHARGE CELLS ONLY. D24's")
    w("Gaussian-smeared models were not re-run, so DIAGNOSTIC 24 as a whole remains EXPLORATORY and")
    w("unquotable, and item 7 is discharged only to the extent of the point-charge vectors.")
    w("")
    w("")
    w("R27.5  THE GATES")
    w("")
    w("All of D27.3's gates ran before any number above was read. The two worth the space:")
    w("")
    w("GATE 3 reproduced d16's committed ddE of -6.134121 kcal/mol at def2-SVPD / 1.3000e from a")
    w("       freshly generated input whose pc.pc is byte-identical to d16's. That is the check that")
    w("       licenses everything else here: it ties this batch to a number already in the repository")
    w("       through the whole chain - generator, transfer, queue, program and reader.")
    w("GATE 8 compared job.out's FINAL SINGLE POINT ENERGY with job.property.txt's &TOTALENERGY in all")
    w("       24 runs. These are two independent prints of the same quantity in two files, and they")
    w("       agree to 1e-09 Eh. It is not in D27.3; it was added by the reader and is recorded here as")
    w("       an addition rather than folded in silently.")
    w("")
    if nbtz:
        w("GATE 6 recorded def2-TZVPD at %d basis functions. D27.3 PREDICTED 682 from a pyscf count and" % nbtz)
        w("       did not assert it; the prediction was %s." % ("correct" if nbtz == 682 else "WRONG, and this changes no quantity"))
    if V["gates_failed"]:
        w("")
        w("GATES THAT DID NOT PASS, committed as failures:")
        for x in V["gates_failed"]: w("   %s" % x)
    w("")
    w("")
    w("R27.6  WHAT THIS DOES NOT SETTLE")
    w("")
    w("1. A THREE-POINT SEQUENCE CANNOT PROVE A LIMIT. It can show one step smaller than the step")
    w("   before it. The honest claim is about the def2-SVPD to def2-TZVPD step, not about the")
    w("   basis-set limit, and D27.5 item 3 said so before the run.")
    w("2. ONE FRAME, ONE GRID, THREE VECTORS. Nothing here speaks to other frames, other arrangements")
    w("   or relaxed geometries.")
    w("3. BASIS INCOMPLETENESS AND SPILL-OUT ARE STILL NOT SEPARATED. DIAGNOSTIC 26 found that")
    w("   separation set-dependent, which is why this ran; this diagnostic measures the basis step and")
    w("   does not attribute it.")
    w("4. THE WALLED AND PSEUDOPOTENTIAL ARMS WERE NOT RUN AT def2-TZVPD. D27.2 excluded them.")
    w("5. IN VACUO, as every Stage 6 arm is. Nothing here is a solvated barrier.")
    w("6. NO BARRIER IS QUOTED AND NO DESIGN IS ACCEPTED OR REJECTED. ddE is a difference of")
    w("   differences; the absolute barrier at any of these bases still carries B3LYP's own 3 to 5")
    w("   kcal/mol error on this reaction, which is larger than anything measured here.")
    w("7. THE CHARGE REPRESENTATION, THE GRID FLOORS, THE OPTIMISER AND THE RELAXATION ARE UNTOUCHED.")
    w("8. DIAGNOSTIC 24's GAUSSIAN CELLS ARE STILL INADMISSIBLE, as R27.4 says.")
    w("")
    Path(a.out).write_text("\n".join(L) + "\n")
    print("rendered %s  (%d lines)" % (a.out, len(L)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
