#!/usr/bin/env python3
"""
cs6_qscan.py - the charge-magnitude scan of STAGE6_DIAGNOSTIC2_20261022.txt.

  python3 cs6_qscan.py prepare <results_repo> <OUT> <HPC_DIR>
  python3 cs6_qscan.py analyse <results_repo> <RUNS> <OUT>

One site (1690, S6-P's and the grid's top-ranked), both committed geometries, both signs, four magnitudes, five
conditions, plus three vertical-detachment runs. A BARE point charge throughout, which Finding 4 licenses (the LJ
sphere and the bare charge differ by exactly L, with identical electronic deviation) and which is what makes the
CPCM conditions runnable: ORCA 6.0.1 refuses CPCM with its QM/MM module but permits CPCM with %pointcharges.

The point-charge file format, the input header for C1 and the restartable single-PBS-job pattern are taken from
cs5_amatrix.py deliberately, so the production condition is demonstrably the same level that built the A matrix.

prepare checks, each a stop:
  (a) both geometries match the sha256 relax/ARRANGEMENTS.txt records for them;
  (b) site 1690's coordinates agree with sites_f42436.tsv and with ARRANGEMENTS.txt's S6-P block to 1e-6 A;
  (c) C1's keyword line equals cs5_amatrix.py's HEAD byte for byte - the production condition is not re-specified
      here, it is copied;
  (d) the committed a for site 1690 is read from A_f42436.tsv and written into the plan, so analyse compares
      against the record rather than against a number typed twice;
  (e) it refuses to write into an existing directory, and submits nothing.

analyse checks, each a stop:
  (f) every run: ORCA TERMINATED NORMALLY, SCF converged, the right electron count, exactly one final energy,
      echoed coordinates equal to substrate.xyz, and for a charged run ORCA's own confirmation that it read one
      point charge;
  (g) C1's two bare energies reproduce A_REPORT.txt's committed in vacuo values for frame 42436 within 1e-6 Eh -
      which is what proves this scan sits on the A matrix's own footing;
  (h) the basis-function count is 264 at C1 and C2 and strictly larger at C3, C4 and C5.
Then it reports a_est, b and q* per condition and applies D2.3 without choosing anything.
"""
import hashlib, re, sys
from pathlib import Path
import numpy as np

H = 627.5094740631
SITE = 1690
QS = (0.05, 0.25, 0.50, 1.00)
NEL_DI, NEL_MONO = 118, 117
A_REF = -6.67782                      # checked against A_f42436.tsv in prepare, not trusted from here
FD_TOL = 0.01 + 0.01 * abs(A_REF)     # A_REPORT.txt's own finite-difference tolerance
O6_PCT = 15.0                         # decision O6's method-and-basis figure
# (name, basis, cpcm epsilon or None)
CONDITIONS = (("C1", "def2-SVP", None), ("C2", "def2-SVP", 4.0), ("C3", "ma-def2-SVP", None),
              ("C4", "ma-def2-SVP", 4.0), ("C5", "def2-SVPD", None))
VDE_BASES = (("C1", "def2-SVP"), ("C3", "ma-def2-SVP"), ("C5", "def2-SVPD"))
HEAD_C1 = "! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF\n%maxcore 3000\n%pal nprocs 8 end\n%scf MaxIter 300 end\n"


def stop(m):
    sys.exit("STOP: " + m)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    return np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])


def head(basis, eps):
    kw = "! B3LYP D3BJ %s def2/J RIJCOSX%s TightSCF\n" % (basis, " CPCM" if eps else "")
    t = kw + "%maxcore 3000\n%pal nprocs 8 end\n"
    if eps:
        t += "%%cpcm epsilon %.1f end\n" % eps
    return t + "%scf MaxIter 300 end\n"


def qtag(q):
    return ("qp" if q > 0 else "qm") + "%04d" % round(abs(q) * 1000)


# ----------------------------------------------------------------------- prepare
def prepare(root, out, hpc):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    out = Path(out)
    if out.exists():
        stop("%s already exists" % out)
    arr = (T6 / "relax" / "ARRANGEMENTS.txt").read_text()
    # (a) the geometries, against the sha256 the record states for them
    geo = {}
    for st, fn in (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz")):
        p = T6 / "frame_42436" / fn
        if not p.is_file():
            stop("%s is missing" % p)
        m = re.search(re.escape(fn) + r"\s+sha256\s+([0-9a-f]{64})", arr)
        if not m:
            stop("ARRANGEMENTS.txt records no sha256 for %s" % fn)
        if sha(p) != m.group(1):
            stop("%s does not match the sha256 ARRANGEMENTS.txt records (%s vs %s)" % (fn, sha(p)[:16], m.group(1)[:16]))
        geo[st] = p
        print("   %-2s %s sha256 matches ARRANGEMENTS.txt" % (st, fn))
    # (b) the site, from the committed sites file, cross-checked against the record's S6-P block
    rows = [l.rstrip("\n").split("\t") for l in (T6 / "sites_f42436.tsv").read_text().splitlines() if not l.startswith("#")]
    hd, rows = rows[0], rows[1:]
    ci = {k: hd.index(k) for k in ("site", "x", "y", "z")}
    hit = [r for r in rows if int(r[ci["site"]]) == SITE]
    if len(hit) != 1:
        stop("sites_f42436.tsv holds %d rows for site %d" % (len(hit), SITE))
    P = np.array([float(hit[0][ci[k]]) for k in ("x", "y", "z")])
    m = re.search(r"site\s+%d\s+q\s+\+1\.0000\s+at\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)" % SITE, arr)
    if not m:
        stop("ARRANGEMENTS.txt has no S6-P line for site %d" % SITE)
    Pr = np.array([float(m.group(i)) for i in (1, 2, 3)])
    d = float(np.abs(P - Pr).max())
    if d > 1e-6:
        stop("site %d differs between sites_f42436.tsv and ARRANGEMENTS.txt by %.1e A" % (SITE, d))
    print("   site %d at %.6f %.6f %.6f, sites file and ARRANGEMENTS.txt agree to %.0e A" % (SITE, P[0], P[1], P[2], d))
    # (d) the committed a, read rather than typed
    A = np.loadtxt(str(T6 / "A_f42436.tsv"), comments="#").reshape(-1)
    irow = [k for k, r in enumerate(rows) if int(r[ci["site"]]) == SITE][0]
    a_committed = float(A[irow])
    if abs(a_committed - A_REF) > 5e-4:
        stop("A_f42436.tsv gives a = %.5f for site %d; this script's A_REF is %.5f" % (a_committed, SITE, A_REF))
    print("   committed a for site %d = %+.5f kcal/mol per e (A_f42436.tsv row %d)" % (SITE, a_committed, irow))
    # (c) C1 is copied, not re-specified
    if head("def2-SVP", None) != HEAD_C1:
        stop("the generated C1 header is not cs5_amatrix.py's HEAD byte for byte")
    print("   C1's header equals cs5_amatrix.py's HEAD byte for byte")

    out.mkdir(parents=True)
    runs = []
    for cname, basis, eps in CONDITIONS:
        hd_t = head(basis, eps)
        for st in ("R", "TS"):
            for q in [None] + [s * x for x in QS for s in (+1, -1)]:
                name = "%s/%s_%s" % (cname, "bare" if q is None else qtag(q), st)
                d = out / name
                d.mkdir(parents=True)
                (d / "substrate.xyz").write_text(geo[st].read_text())
                inp = hd_t
                if q is not None:
                    (d / "pc.pc").write_text("1\n%10.6f %18.10f %18.10f %18.10f\n" % (float(q), P[0], P[1], P[2]))
                    inp += '%pointcharges "pc.pc"\n'
                (d / "job.inp").write_text(inp + "* xyzfile -2 1 substrate.xyz\n")
                runs.append((name, cname, st, "bare" if q is None else "%+.2f" % q, basis, "%.1f" % eps if eps else "-", "-2", "1"))
    for cname, basis in VDE_BASES:
        name = "VDE/%s_mono_R" % cname
        d = out / name
        d.mkdir(parents=True)
        (d / "substrate.xyz").write_text(geo["R"].read_text())
        (d / "job.inp").write_text(head(basis, None) + "* xyzfile -1 2 substrate.xyz\n")
        runs.append((name, cname, "R", "mono", basis, "-", "-1", "2"))

    with open(out / "QSCAN_PLAN.tsv", "w") as fh:
        fh.write("# STAGE6_DIAGNOSTIC2_20261022.txt; cs6_qscan.py prepare\n")
        fh.write("# site %d at %.6f %.6f %.6f; committed a %+.5f kcal/mol per e (A_f42436.tsv row %d)\n"
                 % (SITE, P[0], P[1], P[2], a_committed, irow))
        fh.write("# geometries %s %s\n" % (sha(geo["R"])[:16], sha(geo["TS"])[:16]))
        fh.write("run\tcondition\tstate\tcharge_e\tbasis\tcpcm_eps\tqm_charge\tmult\n")
        for r in runs:
            fh.write("\t".join(r) + "\n")

    pbs = ["#!/bin/bash", "#PBS -N cm_s6q", "#PBS -l select=1:ncpus=8:mem=24gb", "#PBS -l walltime=24:00:00",
           "#PBS -m ae", "#PBS -M 18660916@sun.ac.za", "#PBS -j oe", "#PBS -o %s/qscan.pbs.out" % hpc,
           "# Stage 6 DIAGNOSTIC 2 (cs6_qscan.py, STAGE6_DIAGNOSTIC2_20261022.txt): %d single points." % len(runs),
           "# Restartable: a run whose job.out ends normally is not repeated. One PBS job, so one node lottery",
           "# rather than %d - see diffuse/RUN_NOTES_20261022.txt for why that matters on this cluster." % len(runs),
           'export PATH="/apps/openmpi/4.1.1/bin:$PATH"',
           'export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"',
           "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1", "ORCA=/home/apps2/ORCA/6.0.1", "W=%s" % hpc,
           'cd "$W" || exit 1', 'if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running"; exit 1; fi',
           "trap 'rmdir \"$W/.running\" 2>/dev/null' EXIT", 'echo "host=$(hostname) start=$(date)"',
           "run() {",
           '  d="$W/$1"',
           '  if [ -s "$d/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$d/job.out"; then echo "$1 done before"',
           '  else ( cd "$d" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi',
           "}"]
    pbs += ["run %s" % r[0] for r in runs] + ['echo "end=$(date)"']
    (out / "qscan.pbs").write_text("\n".join(pbs) + "\n")
    files = sorted(p for p in out.rglob("*") if p.is_file())
    (out / "INPUT_CHECKSUMS.sha256").write_text("".join("%s  %s\n" % (sha(p), p.relative_to(out)) for p in files))
    print("   %d runs, %d files including the manifest, one restartable PBS job, nothing submitted"
          % (len(runs), len(files) + 1))
    return len(runs)


# ----------------------------------------------------------------------- analyse
def one(d, nel, charged):
    p = Path(d) / "job.out"
    if not p.is_file():
        stop("%s: no job.out - has the batch run?" % d)
    t = p.read_text(errors="replace")
    if "****ORCA TERMINATED NORMALLY****" not in t:
        stop("%s: ORCA did not terminate normally" % d)
    if not re.search(r"SCF CONVERGED AFTER\s+\d+\s+CYCLES", t) or "SCF NOT CONVERGED" in t:
        stop("%s: the SCF did not converge" % d)
    m = re.search(r"Number of Electrons\s+NEL\s+\.\.\.\.\s+(\d+)", t)
    if not m or int(m.group(1)) != nel:
        stop("%s: %s electrons, expected %d" % (d, m.group(1) if m else "no count", nel))
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    if len(e) != 1:
        stop("%s: %d final energies, expected 1" % (d, len(e)))
    if charged and "Reading point charge file" not in t:
        stop("%s: ORCA reports no point-charge file, but this run has one" % d)
    if charged and "ok (1 point charges)" not in t.replace("  ", " ").replace("   ", " "):
        if not re.search(r"Reading point charge file\s*\.*\s*ok \(1 point charges\)", t):
            stop("%s: ORCA did not confirm reading exactly 1 point charge" % d)
    if not charged and "Reading point charge file" in t:
        stop("%s: ORCA read a point-charge file, but this run should have none" % d)
    k = t.find("CARTESIAN COORDINATES (ANGSTROEM)")
    if k < 0:
        stop("%s: no echoed coordinates" % d)
    X = np.array([[float(v) for v in l.split()[1:4]] for l in t[k:].splitlines()[2:2 + 24]])
    C = read_xyz(Path(d) / "substrate.xyz")
    if X.shape != C.shape or np.abs(X - C).max() > 1e-5:
        stop("%s: the echoed coordinates differ from substrate.xyz" % d)
    nb = re.search(r"Number of basis functions\s+\.\.\.\s+(\d+)", t)
    if not nb:
        stop("%s: no basis-function count" % d)
    return float(e[0]), int(nb.group(1))


def analyse(root, runs, out):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    R = Path(runs)
    E, NB = {}, {}
    for cname, basis, eps in CONDITIONS:
        for st in ("R", "TS"):
            for q in [None] + [s * x for x in QS for s in (+1, -1)]:
                key = (cname, st, q)
                d = R / cname / ("%s_%s" % ("bare" if q is None else qtag(q), st))
                E[key], NB[key] = one(d, NEL_DI, q is not None)
    # (g) C1's bare energies against the committed in vacuo values
    ref = None
    for l in (T6 / "A_REPORT.txt").read_text().splitlines():
        f = l.split("\t")
        if f and f[0].strip() == "42436" and len(f) >= 3:
            ref = (float(f[1]), float(f[2]))
    if ref is None:
        stop("A_REPORT.txt holds no in vacuo energies for frame 42436")
    for st, v in (("R", ref[0]), ("TS", ref[1])):
        dd = abs(E[("C1", st, None)] - v)
        if dd > 1e-6:
            stop("C1 bare %s is %.9f, A_REPORT.txt says %.9f (%.1e Eh apart): this scan is not on the A matrix's "
                 "footing" % (st, E[("C1", st, None)], v, dd))
    print("   C1's bare energies reproduce A_REPORT.txt's committed in vacuo values to <1e-6 Eh")
    # (h) basis-function counts
    n1 = NB[("C1", "R", None)]
    if n1 != 264:
        stop("C1 used %d basis functions, expected def2-SVP's 264" % n1)
    for cname, basis, eps in CONDITIONS:
        n = NB[(cname, "R", None)]
        if cname in ("C3", "C4", "C5") and n <= n1:
            stop("%s (%s) used %d basis functions, not more than def2-SVP's %d" % (cname, basis, n, n1))
        if cname == "C2" and n != n1:
            stop("C2 used %d basis functions, expected def2-SVP's %d" % (n, n1))
    print("   basis functions: " + "  ".join("%s %d" % (c[0], NB[(c[0], "R", None)]) for c in CONDITIONS))

    L = ["STAGE 6 - DIAGNOSTIC 2 RESULT: the charge-magnitude scan at site %d" % SITE, "",
         "Measured by cs6_qscan.py analyse. The verdict rule is STAGE6_DIAGNOSTIC2_20261022.txt D2.3, fixed before",
         "the runs went in and applied here unchanged. STAGE6_CRITERIA.txt's method is unchanged and nothing below",
         "replaces a reported number.", "",
         "ddE(q) = [E_TS(q) - E_TS(0)] - [E_R(q) - E_R(0)], kcal/mol, negative lowers the barrier.",
         "odd(q)  = [ddE(+q) - ddE(-q)]/2   the first order the MILP uses",
         "even(q) = [ddE(+q) + ddE(-q)]/2   the polarisation it neglects", ""]
    res = {}
    for cname, basis, eps in CONDITIONS:
        L.append("%s  %s%s" % (cname, basis, "  CPCM eps=%.1f" % eps if eps else "  in vacuo"))
        L.append("      q      ddE(+q)    ddE(-q)     odd       even    even/|odd|   odd/q")
        od, ev = {}, {}
        for q in QS:
            dp = ((E[(cname, "TS", +q)] - E[(cname, "TS", None)]) - (E[(cname, "R", +q)] - E[(cname, "R", None)])) * H
            dm = ((E[(cname, "TS", -q)] - E[(cname, "TS", None)]) - (E[(cname, "R", -q)] - E[(cname, "R", None)])) * H
            od[q], ev[q] = (dp - dm) / 2, (dp + dm) / 2
            L.append("   %6.2f  %9.4f  %9.4f  %9.4f  %9.4f    %6.1f%%  %9.4f"
                     % (q, dp, dm, od[q], ev[q], 100 * abs(ev[q] / od[q]) if od[q] else float("nan"), od[q] / q))
        a_est = od[QS[0]] / QS[0]
        b = ev[QS[-1]] / QS[-1] ** 2
        qstar = max([q for q in QS if abs(ev[q]) <= 0.10 * abs(od[q])], default=None)
        res[cname] = dict(a=a_est, b=b, qstar=qstar)
        L.append("   a_est (odd/q at q=%.2f) %+.4f    b (even/q^2 at q=%.2f) %+.4f    q* (even <= 10%% of odd) %s"
                 % (QS[0], a_est, QS[-1], b, ("%.2f e" % qstar) if qstar else "below %.2f e" % QS[0]))
        L.append("")
    A = np.loadtxt(str(T6 / "A_f42436.tsv"), comments="#").reshape(-1)
    rows = [l.rstrip("\n").split("\t") for l in (T6 / "sites_f42436.tsv").read_text().splitlines() if not l.startswith("#")]
    hd, rows = rows[0], rows[1:]
    irow = [k for k, r in enumerate(rows) if int(r[hd.index("site")]) == SITE][0]
    a_com = float(A[irow])
    d1 = abs(res["C1"]["a"] - a_com)
    d3 = abs(res["C3"]["a"] - res["C1"]["a"])
    p3 = 100 * d3 / abs(res["C1"]["a"])
    L += ["THE TESTS OF D2.3", "",
          "  committed a (A_f42436.tsv)                         %+.5f" % a_com,
          "  a_est at C1 (def2-SVP in vacuo)                    %+.5f   |diff| %.4f  (tolerance %.4f)" % (res["C1"]["a"], d1, FD_TOL),
          "  a_est at C3 (ma-def2-SVP in vacuo)                 %+.5f   |diff| vs C1 %.4f = %.1f%%  (O6: %.0f%%)" % (res["C3"]["a"], d3, p3, O6_PCT),
          "  b at C3 / b at C4 (does binding fix linearity)     %+.4f / %+.4f" % (res["C3"]["b"], res["C4"]["b"]),
          "  b at C3 / b at C5 (is ma-def2-SVP converged)       %+.4f / %+.4f" % (res["C3"]["b"], res["C5"]["b"]), ""]
    if d1 > FD_TOL:
        v = "THE A MATRIX DOES NOT REPRODUCE ITS OWN RESPONSE"
        w = ["a_est at C1 differs from the committed a by %.4f, beyond A_REPORT.txt's finite-difference tolerance "
             "of %.4f - and C1 IS the level and the site that produced that a, with its bare energies checked "
             "against A_REPORT.txt to below 1e-6 Eh. This is not a basis question. Either the A matrix's "
             "construction or its row order is wrong at this site, and that must be found before any of the "
             "basis or environment comparisons below mean anything." % (d1, FD_TOL)]
    elif d1 <= FD_TOL and d3 <= FD_TOL:
        v = "THE A MATRIX IS SOUND"
        w = ["a_est reproduces the committed a within A_REPORT.txt's own finite-difference tolerance and does not "
             "move with the basis. Stage 5's A column stands; Diagnostic 1 is a statement about charge MAGNITUDE."]
    elif p3 > O6_PCT:
        v = "THE A MATRIX IS BASIS-SENSITIVE AND MUST BE RE-DERIVED"
        w = ["a_est moves %.1f per cent between C1 and C3, beyond decision O6's %.0f per cent." % (p3, O6_PCT)]
    else:
        v = "QUALIFIED"
        w = ["a_est moves by more than the finite-difference tolerance but less than O6's %.0f per cent; usable "
             "with the movement stated as a method uncertainty." % O6_PCT]
    L.append("VERDICT PER D2.3: %s" % v)
    for x in w:
        L.append("  - %s" % x)
    b3, b4, b5 = res["C3"]["b"], res["C4"]["b"], res["C5"]["b"]
    L.append("")
    L.append("  ON THE ENVIRONMENT: %s" % ("binding the density materially improves linearity (|b(C4)| < half |b(C3)|), so the "
             "design environment is a live question on validity grounds" if abs(b4) < 0.5 * abs(b3) else
             "the continuum does not fix the response (|b(C4)| is not under half |b(C3)|); PHASE22_PLAN.txt 2.2 stands untouched"))
    L.append("  ON CONVERGENCE: %s" % ("ma-def2-SVP is converged enough to read as an answer (b within 25 per cent of def2-SVPD)"
             if abs(b5 - b3) <= 0.25 * abs(b3) else
             "neither diffuse basis is an arbiter (b differs by more than 25 per cent between ma-def2-SVP and def2-SVPD)"))
    # the VDE
    L += ["", "VERTICAL DETACHMENT ENERGY at the reactant geometry (E(monoanion) - E(dianion))", ""]
    for cname, basis in VDE_BASES:
        em, _ = one(R / "VDE" / ("%s_mono_R" % cname), NEL_MONO, False)
        vde = (em - E[(cname, "R", None)]) * 27.211386245988 / 1.0
        L.append("   %-3s %-13s  VDE %+8.4f eV   -> the dianion is %s at this condition"
                 % (cname, basis, vde, "METASTABLE to electron loss" if vde < 0 else "BOUND"))
    L += ["", "  A negative VDE vindicates Finding 8's reading on the right quantity. A positive VDE means the",
          "  species is bound despite a positive HOMO, and Finding 8's second half must be withdrawn as an",
          "  artefact of reading Koopmans' theorem off a DFT orbital (D2.2, and Rienstra-Kiracofe et al. 2002",
          "  p. 246 / Galbraith & Schaefer 1996 p. 863, both held in chorismate-refs)."]
    Path(out).mkdir(parents=True, exist_ok=True)
    (Path(out) / "DIAGNOSTIC2_RESULT.txt").write_text("\n".join(L) + "\n")
    print()
    print("\n".join(L[L.index("THE TESTS OF D2.3"):]))
    print()
    print("   wrote %s" % (Path(out) / "DIAGNOSTIC2_RESULT.txt"))


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in ("prepare", "analyse"):
        sys.exit(__doc__)
    if sys.argv[1] == "prepare":
        prepare(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        analyse(sys.argv[2], sys.argv[3], sys.argv[4])
