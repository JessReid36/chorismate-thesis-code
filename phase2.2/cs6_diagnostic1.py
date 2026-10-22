#!/usr/bin/env python3
"""
cs6_diagnostic1.py - measure DIAGNOSTIC 1 and apply the verdict rule that was fixed before the runs went in.

  python3 cs6_diagnostic1.py <results_repo> <outputs_dir> <OUT>

<outputs_dir> holds the eight downloaded ma-def2-SVP runs as DIFF-<design>/{R,TS}_sp/.

STAGE6_DIAGNOSTIC1_20261022.txt D1.3 fixes, in advance, what counts as which answer:

  STAND      every |ddE| moves by less than 15 per cent (decision O6's method-and-basis figure, used as a
             yardstick and not as a criterion), AND the sign of the measured-minus-predicted gap is unchanged for
             every design, AND the order of |gap| by n is unchanged.
  QUALIFIED  a |ddE| moves by 15 per cent or more, but every sign and the order survive.
  NOT USABLE a gap's sign reverses, or the order by n changes.

This script does not choose. It computes the three tests and reports which of the three the result is, and it
prints every number it used so the verdict can be checked by hand.

The energy reading is cs6_analyse.py's, deliberately unchanged and for its stated reason: a QM/MM job prints
three FINAL SINGLE POINT ENERGY lines and a plain one prints one, so the three are read apart rather than by
allowing an optional suffix, and the QM/MM TOTAL is the quantity because the score it is measured against
contains L.

Checks, each a stop:
  (a) all eight runs present, ORCA TERMINATED NORMALLY, the SCF converged, 118 electrons in the QM region, the
      energy lines in the structure the job kind requires, the echoed coordinates equal to the committed
      combined.xyz, and the MM point-charge count equal to the arrangement's sphere count;
  (b) the basis ACTUALLY used is larger than the committed def2-SVP run's for the same folder - which is what
      establishes that the diagnostic changed the basis and not merely the input text;
  (c) the committed frozen.tsv holds a row for each design measured here, and its control energies are the ones
      the def2-SVP comparison is made against.
"""
import hashlib, lzma, re, sys, tarfile
from pathlib import Path
import numpy as np

H = 627.5094740631
NQM, NEL = 24, 118
NSPH = {"S6-C": 0, "S6-P": 1, "S6-2": 2, "S6-20": 20}
ORDER = ["S6-C", "S6-P", "S6-2", "S6-20"]
ROOT, OUTS, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
T6 = Path(ROOT) / "phase2.2" / "cs_stage6"


def stop(m):
    sys.exit("STOP: " + m)


def read_xyz(p):
    L = Path(p).read_text().splitlines()
    return np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + int(L[0].split()[0])]])


def homo(t):
    """The highest orbital with an occupation of 2, in eV, from ORCA's ORBITAL ENERGIES table."""
    k = t.find("ORBITAL ENERGIES")
    if k < 0:
        stop("no ORBITAL ENERGIES block")
    h = None
    for l in t[k:].splitlines()[1:]:
        if "MOLECULAR ORBITALS" in l or "DENSITY" in l:
            break
        f = l.split()
        if len(f) == 4:
            try:
                no, occ, eh, ev = int(f[0]), float(f[1]), float(f[2]), float(f[3])
            except ValueError:
                continue
            if occ >= 1.9:
                h = ev
    if h is None:
        stop("no occupied orbital found in the ORBITAL ENERGIES block")
    return h


def nbas(t):
    m = re.search(r"Number of basis functions\s+\.\.\.\s+(\d+)", t)
    if not m:
        stop("no basis-function count in the output")
    return int(m.group(1))


def read_run(d, nsph, xyz):
    """cs6_analyse.py's reader, unchanged. Returns (total energy in Eh, HOMO in eV, basis functions).

    xyz is the COMMITTED combined.xyz of the same folder, not a copy that came down with the output: the check is
    that the geometry ORCA echoed is the geometry the repository holds.
    """
    d = Path(d)
    if not Path(xyz).is_file():
        stop("%s is missing - the committed diagnostic input is not where it should be" % xyz)
    if (d / "job.out").is_file():
        t = (d / "job.out").read_text(errors="replace")
    elif (d / "job.out.xz").is_file():
        t = lzma.open(str(d / "job.out.xz"), "rt", errors="replace").read()
    else:
        stop("%s: no job.out or job.out.xz - has the job run?" % d)
    if "****ORCA TERMINATED NORMALLY****" not in t:
        stop("%s: ORCA did not terminate normally" % d)
    if not re.search(r"SCF CONVERGED AFTER\s+\d+\s+CYCLES", t) or "SCF NOT CONVERGED" in t:
        stop("%s: the SCF did not converge" % d)
    m = re.search(r"Number of Electrons\s+NEL\s+\.\.\.\.\s+(\d+)", t)
    if not m or int(m.group(1)) != NEL:
        stop("%s: %s electrons in the QM region, expected %d" % (d, m.group(1) if m else "no count", NEL))
    qm = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    mm = re.findall(r"FINAL SINGLE POINT ENERGY \(MM\)\s+(-?\d+\.\d+)", t)
    tot = re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", t)
    if nsph == 0:
        if mm or tot:
            stop("%s: the control reports MM or QM/MM energies; it is meant to be a plain QM job" % d)
        if len(qm) != 1:
            stop("%s: %d final energies, expected 1" % (d, len(qm)))
        e_tot = float(qm[0])
    else:
        if not (len(qm) == len(mm) == len(tot) == 1):
            stop("%s: %d bare, %d (MM) and %d (QM/MM) final energies, expected one of each"
                 % (d, len(qm), len(mm), len(tot)))
        e_tot = float(tot[0])
        if abs(float(qm[0]) + float(mm[0]) - e_tot) > 1e-9:
            stop("%s: the QM energy plus the MM energy is not the QM/MM total" % d)
    k = t.find("CARTESIAN COORDINATES (ANGSTROEM)")
    if k < 0:
        stop("%s: no echoed coordinates" % d)
    X = np.array([[float(v) for v in l.split()[1:4]] for l in t[k:].splitlines()[2:2 + NQM]])
    C = read_xyz(xyz)[:NQM]
    if X.shape != (NQM, 3) or np.abs(X - C).max() > 1e-5:
        stop("%s: the echoed coordinates differ from the committed combined.xyz by %.1e A"
             % (d, np.abs(X - C).max()))
    pc = re.search(r"Point charges in QM calc\.*\s*from MM atoms\.*\s*(\d+)", t)
    if nsph == 0:
        if "Multiscale model" in t or pc:
            stop("%s: the control ran as QM/MM; it is meant to be a plain QM job" % d)
    elif not pc:
        if "Multiscale model" not in t:
            stop("%s: ORCA reports no MM region, but this arrangement has %d sphere(s)" % (d, nsph))
    elif int(pc.group(1)) != nsph:
        stop("%s: ORCA read %s MM point charges, the arrangement has %d" % (d, pc.group(1), nsph))
    return e_tot, homo(t), nbas(t)


# ---- (c) the committed def2-SVP record this is compared against
ftsv = T6 / "frozen" / "frozen.tsv"
if not ftsv.is_file():
    stop("%s is missing" % ftsv)
rows = [l.split("\t") for l in ftsv.read_text().splitlines() if l.strip()]
hd, rows = rows[0], rows[1:]
F = {r[0]: dict(zip(hd, r)) for r in rows}
for k in ORDER:
    if k != "S6-C" and k not in F:
        stop("the committed frozen.tsv holds no row for %s" % k)

# ---- (b) the basis the committed def2-SVP runs actually used, per folder
svp_nbas = {}
with tarfile.open(str(T6 / "frozen" / "sp_outputs.tar.xz")) as t:
    for n in t.getnames():
        m = re.match(r"^\.?/?(S6-[A-Z0-9]+)/(R|TS)_sp/job\.out$", n)
        if m and m.group(1) in NSPH:
            svp_nbas[(m.group(1), m.group(2))] = nbas(t.extractfile(n).read().decode(errors="replace"))
if len(svp_nbas) != 8:
    stop("the committed frozen archive holds %d of the 8 def2-SVP outputs for these four designs" % len(svp_nbas))

# ---- (a) read the eight diffuse runs
E, HO, NB = {}, {}, {}
for des in ORDER:
    for g in ("R", "TS"):
        d = Path(OUTS) / ("DIFF-%s" % des) / ("%s_sp" % g)
        if not d.is_dir():
            stop("%s is missing - were all eight downloaded?" % d)
        xyz = T6 / "diffuse" / ("DIFF-%s" % des) / ("%s_sp" % g) / "combined.xyz"
        e, h, nb = read_run(d, NSPH[des], xyz)
        if nb <= svp_nbas[(des, g)]:
            stop("%s used %d basis functions, the committed def2-SVP run used %d: the basis did not grow, so this "
                 "is not a diffuse-basis calculation" % (d, nb, svp_nbas[(des, g)]))
        E[(des, g)], HO[(des, g)], NB[(des, g)] = e, h, nb
        print("   %-13s %-2s  E %.9f Eh  HOMO %+7.4f eV  %d basis functions (def2-SVP: %d)"
              % ("DIFF-" + des, g, e, h, nb, svp_nbas[(des, g)]))

# ---- the measurement
bar = {d: (E[(d, "TS")] - E[(d, "R")]) * H for d in ORDER}
rep = []
for des in ("S6-P", "S6-2", "S6-20"):
    ma = bar[des] - bar["S6-C"]
    svp = float(F[des]["ddE_frozen_kcal"])
    score = float(F[des]["score_kcal"])
    rep.append({"design": des, "n": int(F[des]["n"]), "ddE_ma": ma, "ddE_svp": svp, "score": score,
                "gap_ma": ma - score, "gap_svp": svp - score,
                "move_pct": 100.0 * abs(ma - svp) / abs(svp)})

# ---- the three tests of D1.3
worst = max(r["move_pct"] for r in rep)
signs_ok = all(np.sign(r["gap_ma"]) == np.sign(r["gap_svp"]) for r in rep)
by_n = sorted(rep, key=lambda r: r["n"])
ord_ma = [abs(r["gap_ma"]) for r in by_n]
ord_svp = [abs(r["gap_svp"]) for r in by_n]
order_ok = (all(ord_ma[i] < ord_ma[i + 1] for i in range(len(ord_ma) - 1))
            == all(ord_svp[i] < ord_svp[i + 1] for i in range(len(ord_svp) - 1))
            and all((ord_ma[i] < ord_ma[i + 1]) == (ord_svp[i] < ord_svp[i + 1]) for i in range(len(ord_ma) - 1)))

if not (signs_ok and order_ok):
    verdict, why = "NOT USABLE AS THEY STAND", []
    if not signs_ok:
        why.append("the sign of the measured-minus-predicted gap reverses for at least one design")
    if not order_ok:
        why.append("the order of the gap by n changes")
elif worst >= 15.0:
    verdict = "QUALIFIED"
    why = ["every sign and the order survive, but the largest ddE moves by %.1f per cent, at or above decision "
           "O6's 15 per cent method-and-basis figure" % worst]
else:
    verdict = "STAND"
    why = ["the largest ddE moves by %.1f per cent, below decision O6's 15 per cent, and every sign and the "
           "order are unchanged" % worst]

# ---- the report
L = []
L.append("STAGE 6 - DIAGNOSTIC 1 RESULT: ma-def2-SVP against the committed def2-SVP frozen single points")
L.append("")
L.append("Measured by cs6_diagnostic1.py from the eight committed outputs. The verdict rule is")
L.append("STAGE6_DIAGNOSTIC1_20261022.txt D1.3, fixed before the runs went in and applied here unchanged.")
L.append("STAGE6_CRITERIA.txt's B3LYP-D3BJ/def2-SVP is unchanged and no number below replaces a reported one.")
L.append("")
L.append("BARE BARRIERS AND THE HOMO, ma-def2-SVP")
L.append("")
L.append("  design    n  E(R) Eh            E(TS) Eh           TS-R kcal  HOMO(R) eV  HOMO(TS) eV  basis")
for des in ORDER:
    L.append("  %-8s %2d  %.9f  %.9f  %9.4f  %+10.4f  %+11.4f  %4d/%d"
             % (des, NSPH[des], E[(des, "R")], E[(des, "TS")], bar[des], HO[(des, "R")], HO[(des, "TS")],
                NB[(des, "R")], svp_nbas[(des, "R")]))
L.append("")
L.append("  The control's HOMO at def2-SVP was +2.278 eV (Finding 8). A larger basis does not make a gas-phase")
L.append("  dianion bound; it makes the calculation honest about it, so a more positive value here is the")
L.append("  expected direction and not a failure.")
L.append("")
L.append("THE QUANTITY D1.3 NAMES: ddE(frozen) AGAINST THE CONTROL OF THE SAME BASIS")
L.append("")
L.append("  design    n  ddE ma-def2-SVP  ddE def2-SVP  move      score (def2-SVP)  gap ma    gap def2-SVP")
for r in rep:
    L.append("  %-8s %2d  %14.4f  %12.4f  %6.1f%%  %16.4f  %8.4f  %12.4f"
             % (r["design"], r["n"], r["ddE_ma"], r["ddE_svp"], r["move_pct"], r["score"], r["gap_ma"], r["gap_svp"]))
L.append("")
L.append("  gap = measured ddE minus the linear score the optimiser minimised. The score is NOT recomputed: it")
L.append("  comes from the committed A and L matrices at def2-SVP, which is the point - the question is whether")
L.append("  the gap Finding 1 read as polarisation survives a basis that lets loosely-held density go.")
L.append("")
L.append("THE THREE TESTS")
L.append("")
L.append("  largest |ddE| move          %.1f%%  (D1.3 yardstick: 15%%)" % worst)
L.append("  gap signs unchanged         %s" % ("yes" if signs_ok else "NO"))
L.append("  order of |gap| by n         %s" % ("unchanged" if order_ok else "CHANGED"))
L.append("")
L.append("VERDICT PER D1.3: THE FINDINGS %s" % verdict)
for w in why:
    L.append("  - %s" % w)
L.append("")
L.append("WHAT THIS DOES NOT TOUCH")
L.append("")
L.append("  Finding 7. The control's 33.8 kcal/mol descent from frame 42436's reactant geometry is a statement")
L.append("  about the geometry, not the basis, and is unchanged by anything here. The relaxation question and the")
L.append("  OptTS/Cartesian protocol conflict are both still open.")
L.append("")
L.append("FILES")
L.append("  phase2.2/cs_stage6/diffuse/outputs/diff_outputs.tar.xz   the eight outputs, committed before measuring")
L.append("  phase2.2/cs_stage6/diffuse/outputs/diff_outputs.MANIFEST.sha256")
L.append("  phase2.2/cs_stage6/diffuse/DIAGNOSTIC1_RESULT.txt        this file")
L.append("  phase2.2/cs_stage6/diffuse/diagnostic1.tsv               the same numbers machine-readable")

Path(OUT).mkdir(parents=True, exist_ok=True)
(Path(OUT) / "DIAGNOSTIC1_RESULT.txt").write_text("\n".join(L) + "\n")
with open(Path(OUT) / "diagnostic1.tsv", "w") as fh:
    fh.write("design\tn\tddE_ma_kcal\tddE_svp_kcal\tmove_pct\tscore_kcal\tgap_ma_kcal\tgap_svp_kcal\n")
    for r in rep:
        fh.write("%s\t%d\t%.6f\t%.6f\t%.4f\t%.6f\t%.6f\t%.6f\n"
                 % (r["design"], r["n"], r["ddE_ma"], r["ddE_svp"], r["move_pct"], r["score"], r["gap_ma"], r["gap_svp"]))
    fh.write("#\tHOMO_eV\n")
    for des in ORDER:
        fh.write("# %s\tR %+.4f\tTS %+.4f\n" % (des, HO[(des, "R")], HO[(des, "TS")]))
    fh.write("# verdict\t%s\n" % verdict)
print()
print("\n".join(L[L.index("THE THREE TESTS"):]))
print()
print("   wrote %s and diagnostic1.tsv" % (Path(OUT) / "DIAGNOSTIC1_RESULT.txt"))
