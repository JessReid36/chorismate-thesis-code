#!/usr/bin/env python3
"""
cs4b_explore.py - EXPLORATORY analysis of Stage 4b (not pre-registered; written 8 October 2026 after Stage 4b's
results were seen). It does not change any Stage 4b verdict; it separates the two effects the pre-registered b
criterion mixed, to inform Stage 4c (see STAGE4B_EXPLORATORY.txt).
  curvature from the negative side only (q = -1, -0.5, 0): a negative sphere cannot trap electrons, so the
      def2-SVPD / def2-SVP ratio here measures basis-set polarisation;
  curvature from the positive side only (q = 0, +0.5, +1); the excess of positive over negative side measures the
      trap; HOMO-LUMO gap at q = +1.
Reads each job.out, or job.out.xz if that is what is present (committed outputs). Python 3.6, standard library.
USAGE  python3 cs4b_explore.py <phase2.2 folder>   -> writes cs_stage4b/exploratory/STAGE4B_EXPLORATORY_REPORT.txt
"""
import lzma, os, re, sys
P2 = sys.argv[1]
B4, S4 = os.path.join(P2, "cs_stage4b"), os.path.join(P2, "cs_stage4")
H2KCAL = 627.5094740631
def text(d):
    p = os.path.join(d, "job.out")
    if os.path.exists(p):
        return open(p, errors="replace").read()
    if os.path.exists(p + ".xz"):
        return lzma.open(p + ".xz", "rt", errors="replace").read()
    raise SystemExit("missing output in " + d)
def parse(d, qmmm):
    t = text(d)
    assert "ORCA TERMINATED NORMALLY" in t, d
    pat = r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)" if qmmm else r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)"
    e = float(re.findall(pat, t)[-1])
    i = t.rfind("ORBITAL ENERGIES"); orbs = []
    for l in t[i:i + 40000].splitlines()[4:]:
        m = re.match(r"\s*(\d+)\s+(\d\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", l)
        if m: orbs.append((float(m.group(2)), float(m.group(4))))
        elif orbs: break
    gap = min(o[1] for o in orbs if o[0] < 0.5) - max(o[1] for o in orbs if o[0] > 0.5)
    return e, gap
def measures(E):
    e = lambda q: E[q][0]
    D1 = (e(1.0) - e(-1.0)) * H2KCAL; D5 = (e(0.5) - e(-0.5)) * H2KCAL
    S1 = (e(1.0) + e(-1.0) - 2 * e(0.0)) * H2KCAL; S5 = (e(0.5) + e(-0.5) - 2 * e(0.0)) * H2KCAL
    a = (8 * D5 - D1) / 6; b = (16 * S5 - S1) / 6
    neg = 2 * ((e(-1.0) - e(0.0)) * H2KCAL - 2 * (e(-0.5) - e(0.0)) * H2KCAL)
    pos = 2 * ((e(1.0) - e(0.0)) * H2KCAL - 2 * (e(0.5) - e(0.0)) * H2KCAL)
    return a, b, neg, pos, E[1.0][1]
QS = ((-1.0, "qm1"), (-0.5, "qm05"), (0.0, "q0"), (0.5, "qp05"), (1.0, "qp1"))
rows = []; R = {}
for dt, d in (("d2674", 2.674), ("d290", 2.90), ("d320", 3.20), ("d350", 3.50)):
    for bt in ("svp", "svpd"):
        E = {}
        for q, qt in QS:
            run = os.path.join(S4, "sp_%s_%s" % (bt, qt)) if dt == "d2674" else os.path.join(B4, "A_%s_%s_%s" % (bt, dt, qt))
            E[q] = parse(run, True)
        R["bare", d, bt] = measures(E)
for bt in ("svp", "svpd"):
    E = {}
    for q, qt in QS:
        run = os.path.join(B4, "B2_inline_%s_q0" % bt) if qt == "q0" else os.path.join(B4, "B_ecp_%s_%s" % (bt, qt))
        E[q] = parse(run, False)
    R["ecp", 2.674, bt] = measures(E)
out = ["Stage 4b EXPLORATORY report (cs4b_explore.py) - not a pre-registered test; Stage 4b's verdicts stand as recorded.",
       "Curvatures in kcal/mol (second-order coefficient of E(q)); gap at q = +1 in eV.",
       "", "%-26s %-5s %9s %9s %9s %9s %9s %8s" % ("sphere, distance from O2", "basis", "a", "b", "neg-side", "pos-side", "excess", "gap(+1)")]
for key in [("bare", d, bt) for d in (2.674, 2.90, 3.20, 3.50) for bt in ("svp", "svpd")] + [("ecp", 2.674, "svp"), ("ecp", 2.674, "svpd")]:
    a, b, neg, pos, gap = R[key]
    lab = "%s %.3f A" % ("bare sphere" if key[0] == "bare" else "Ne-type ECP sphere", key[1])
    out.append("%-26s %-5s %9.3f %9.3f %9.3f %9.3f %+9.3f %8.2f" % (lab, key[2], a, b, neg, pos, pos - neg, gap))
out.append("")
out.append("negative-side curvature, def2-SVPD / def2-SVP (basis-set polarisation; no trap possible):")
for kind, d in [("bare", 2.674), ("bare", 2.90), ("bare", 3.20), ("bare", 3.50), ("ecp", 2.674)]:
    out.append("   %-26s %.3f" % ("%s %.3f A" % ("bare sphere" if kind == "bare" else "Ne-type ECP sphere", d), R[kind, d, "svpd"][2] / R[kind, d, "svp"][2]))
out.append("pseudopotential effect on def2-SVP (ECP vs bare, 2.674 A): a %+.2f%%, b %+.2f%%" % (
    100 * (R["ecp", 2.674, "svp"][0] - R["bare", 2.674, "svp"][0]) / abs(R["bare", 2.674, "svp"][0]),
    100 * (R["ecp", 2.674, "svp"][1] - R["bare", 2.674, "svp"][1]) / abs(R["bare", 2.674, "svp"][1])))
os.makedirs(os.path.join(B4, "exploratory"), exist_ok=True)
open(os.path.join(B4, "exploratory", "STAGE4B_EXPLORATORY_REPORT.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
