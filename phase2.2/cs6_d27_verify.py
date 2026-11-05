#!/usr/bin/env python3
"""D27.3 GATE 7, run as its own program.

Every figure in the rendered result document is re-derived here FROM THE OUTPUTS, not from the
reader's values file, and by a deliberately different route: the reader took each energy from
job.out's FINAL SINGLE POINT ENERGY, and this takes it from job.property.txt's &TOTALENERGY.
Two files, two regexes, no shared code. Then every numeric token in the document must be
accounted for - matched to a derived quantity or to the small list of figures the document
legitimately quotes from elsewhere. An unaccounted number is a failure.
"""
import argparse, re, sys
from pathlib import Path

H = 627.5094740631
TOL = 1.0
BASES = ["def2-SVP", "def2-SVPD", "def2-TZVPD"]
VECTORS = ["1.3000e_n10", "2.8800e_n10", "5.6433e_n20"]
SUMQ = {"1.3000e_n10": 1.3000, "2.8800e_n10": 2.8800, "5.6433e_n20": 5.6433}
NCH = {"1.3000e_n10": 10, "2.8800e_n10": 10, "5.6433e_n20": 20}

# figures the document quotes from OUTSIDE this run, each with where it comes from
QUOTED = {
    13.5154: "D24's inter-basis gap at the third vector, from d24_smear.log",
    2.28:    "the factor between the two bases there, from the same log",
    -5.6165: "D26's quantity A, pyscf's four-decimal print",
    -6.134121: "d16's committed def2-SVPD ddE, which GATE 3 reproduces",
    18.332256: "d16's committed def2-SVPD control barrier",
    2.88:    "the sum|q| above which the ladder is in question, from D20/D5.3(d)",
    1.0:     "TOL, from criterion D5.3(d)",
    682:     "the predicted def2-TZVPD basis count",
    402:     "def2-SVPD's basis count, ORCA's own",
    24:      "the number of single points",
    27: "the diagnostic number", 26: "the diagnostic it follows", 23: "ditto", 5: "ditto",
    3: "B3LYP's error range, lower", 6: "a section number", 7: "a section number",
    8: "a section number", 4: "a section number", 2: "a section number", 1: "a section number",
}
ALLOWED_INT = set(range(0, 29)) | {118, 264, 402, 682, 2026}


def energies(tree):
    """primary source: job.property.txt's &TOTALENERGY. Deliberately NOT job.out."""
    E = {}
    for b in BASES:
        for n in ["ctl_R", "ctl_TS"] + ["%s_%s" % (v, s) for v in VECTORS for s in ("R", "TS")]:
            p = tree / b / n / "job.property.txt"
            if not p.exists(): return None, "%s/%s has no job.property.txt" % (b, n)
            m = re.findall(r"&TOTALENERGY\s*\[&Type\s+\"Double\"\]\s*(-?[\d.]+[eE][+-]\d+)", p.read_text(errors="replace"))
            if not m: return None, "%s/%s job.property.txt has no &TOTALENERGY" % (b, n)
            E[(b, n)] = float(m[-1])
    return E, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True); ap.add_argument("--doc", required=True)
    a = ap.parse_args()
    tree, doc = Path(a.tree), Path(a.doc)
    D = doc.read_text()
    rows, fails = [], []

    def chk(ok, msg):
        rows.append("   %s  GATE 7  %s" % ("pass" if ok else "FAIL", msg))
        if not ok: fails.append(msg)

    E, err = energies(tree)
    if E is None:
        print("   FAIL  GATE 7  %s" % err); return 1

    derived = {}
    ctl, ddE = {}, {}
    for b in BASES:
        ctl[b] = (E[(b, "ctl_TS")] - E[(b, "ctl_R")]) * H
        derived[round(ctl[b], 6)] = "control barrier, %s" % b
        for v in VECTORS:
            ddE[(b, v)] = (E[(b, v + "_TS")] - E[(b, v + "_R")]) * H - ctl[b]
            derived[round(ddE[(b, v)], 6)] = "ddE %s %s" % (b, v)
    verdicts = {}
    for v in VECTORS:
        gsd, gdt = abs(ddE[("def2-SVPD", v)] - ddE[("def2-SVP", v)]), abs(ddE[("def2-TZVPD", v)] - ddE[("def2-SVPD", v)])
        derived[round(gsd, 6)] = "gap_SD %s" % v
        derived[round(gdt, 6)] = "gap_DT %s" % v
        derived[round(SUMQ[v], 4)] = "sum|q| %s" % v
        verdicts[v] = ("NOT CONVERGING" if gdt >= gsd else ("NOT CONVERGED" if gdt > TOL else "CONVERGED"))
    v3 = VECTORS[2]
    d_svpd = abs(ddE[("def2-TZVPD", v3)] - ddE[("def2-SVPD", v3)])
    d_svp = abs(ddE[("def2-TZVPD", v3)] - ddE[("def2-SVP", v3)])
    derived[round(d_svpd, 6)] = "|TZVPD-SVPD| at the third vector"
    derived[round(d_svp, 6)] = "|TZVPD-SVP| at the third vector"
    arb = ("def2-SVP IS THE ARTEFACT" if d_svpd <= TOL else
           ("def2-SVPD IS THE ARTEFACT" if d_svp <= TOL else "BOTH ARE ARTEFACTS"))
    overall = (list(verdicts.values())[0] if len(set(verdicts.values())) == 1
               else "PER VECTOR; NO SINGLE VERDICT IS FORMED")

    # ---- the document must carry each derived quantity, at full precision
    for val, what in sorted(derived.items()):
        fmt6, fmt4 = "%.6f" % val, "%.4f" % val
        chk(fmt6 in D or fmt4 in D, "%-34s %s present in the document" % (what, fmt6))
    # ---- and the verdicts and the arbitration, as words
    for v in VECTORS:
        chk(re.search(re.escape(v) + r".*" + re.escape(verdicts[v]), D),
            "%-12s is reported as %s" % (v, verdicts[v]))
    chk("ARBITRATION, by the rule D27.4 fixed before the run: %s" % arb in D,
        "the arbitration sentence reads %s" % arb)
    chk("OVERALL: %s" % overall in D, "OVERALL reads %s" % overall)
    # ---- the headline must agree with the verdict, not merely sit above it
    # the headline is everything between the title line and the provenance paragraph; a fixed
    # line slice ran into "which is NOT edited" and made the test fire on the wrong text
    head = D.split("\n\n")[1] if "\n\n" in D else ""
    chk("Pre-registered" not in head and len(head) > 40, "the headline block was located cleanly")
    if overall == "CONVERGED":
        chk("IS CONVERGED" in head and "NOT" not in head.replace("NOTHING", ""),
            "the headline says converged, and the verdict is CONVERGED")
    else:
        chk("DO NOT AGREE" in head or "NOT CONVERG" in head,
            "the headline does not claim convergence, and the verdict is %s" % overall)
    # ---- no number in the document may be unaccounted for
    toks = re.findall(r"(?<![\w.])-?\d+(?:\.\d+)?(?![\w.])", D)
    unaccounted = []
    for t in set(toks):
        f = float(t)
        if "." not in t:
            if int(f) in ALLOWED_INT or int(f) in QUOTED: continue
            unaccounted.append(t); continue
        if any(abs(f - k) < 5e-7 for k in derived) or any(abs(f - k) < 5e-7 for k in QUOTED): continue
        if any(abs(f - round(k, 4)) < 5e-5 for k in derived): continue
        unaccounted.append(t)
    chk(not unaccounted, "every numeric token in the document is accounted for%s"
        % ("" if not unaccounted else "; unaccounted: " + ", ".join(sorted(unaccounted)[:8])))

    print("\n".join(rows)); print()
    if fails:
        print("GATE 7 FAILED: %d" % len(fails))
        for f in fails: print("   " + f)
        return 1
    print("GATE 7 PASSES: %d independent checks, every figure re-derived from job.property.txt." % len(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
