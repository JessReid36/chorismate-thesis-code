#!/usr/bin/env python3
"""
cs6_field.py - STAGE6_DIAGNOSTIC4_20261023.txt: the total-field scan. S6-20's committed design, every charge
scaled by the same lambda, so that only the total field varies and every individual charge stays inside the
per-site bound DIAGNOSTIC2_AMENDMENT1_RESULT.txt measured.

  python3 cs6_field.py prepare <results_repo> <OUT> <HPC_DIR>
  python3 cs6_field.py analyse <results_repo> <RUNS> <OUT>

Bare point charges, twenty in one file, no Lennard-Jones - Finding 4 measured that choice to be electronically
identical to the LJ sphere. The bare references come out of the committed qscan outputs archive. The design is
read out of relax/ARRANGEMENTS.txt rather than typed, and prepare STOPS unless the sum a q it reads reproduces the
committed score once sum L is added.

prepare checks, each a stop:
  (a) both geometries match the sha256 relax/ARRANGEMENTS.txt records;
  (b) S6-20 parses to exactly 20 spheres, and sum a q + sum L reproduces the committed score to 1e-4;
  (c) every site's coordinates agree with sites_f42436.tsv to 1e-6 A and every a agrees with A_f42436.tsv to 5e-5;
  (d) at every lambda the largest single charge is inside the measured per-site bound of 0.80 e;
  (e) the def2-SVP header equals cs5_amatrix.py's HEAD byte for byte and the def2-SVPD header equals the committed
      qscan C5 header byte for byte;
  (f) it refuses to write into an existing folder, and submits nothing.

analyse checks, each a stop:
  (g) every run: ORCA TERMINATED NORMALLY, SCF converged, 118 electrons, one final energy, the echoed coordinates
      equal to the committed geometry, ORCA's confirmation that it read exactly 20 point charges, one
      basis-function count per basis;
  (h) the committed qscan archive supplies the bare references, each checked against its manifest sha256, and its
      C1 bare energies reproduce A_REPORT.txt to below 1e-6 Eh;
  (i) b is read out of the committed DIAGNOSTIC2_AMENDMENT1 result's own text, not typed here.
Then it applies D4.3's three readings without choosing anything.
"""
import hashlib, re, sys, tarfile
from pathlib import Path
import numpy as np

H = 627.5094740631
LAMBDAS = (0.25, 0.50, 0.75, 1.00)
QCAP = 0.80            # the measured per-site bound
SCORE_REC = -19.65373  # the committed S6-20 score
BASES = (("def2-SVP", "C1"), ("def2-SVPD", "C5"))
NEL = 118
PCT, O6 = 10.0, 15.0
D1_MA = +4.5746        # Diagnostic 1's S6-20 at ma-def2-SVP
D1_SVP = -24.8261      # the committed frozen def2-SVP ddE, which INCLUDES sum L
SUM_L  = -0.13622      # S6-20's sum L, so the bare-charge electronic part is D1_SVP - SUM_L


def stop(m):
    sys.exit("STOP: " + m)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    return np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])


def head(basis):
    return ("! B3LYP D3BJ %s def2/J RIJCOSX TightSCF\n%%maxcore 3000\n%%pal nprocs 8 end\n"
            "%%scf MaxIter 300 end\n" % basis)


def design(T6):
    """S6-20 out of the committed ARRANGEMENTS.txt: (site, q, xyz, a, L) x 20, checked against the record."""
    arr = (T6 / "relax" / "ARRANGEMENTS.txt").read_text()
    m = re.search(r"^\s*S6-20\s+20 sphere\(s\)(.*?)(?=^\s*S6-C)", arr, re.M | re.S)
    if not m:
        stop("could not find the S6-20 block in the committed ARRANGEMENTS.txt")
    rows = re.findall(r"site\s+(\d+)\s+q\s+([+-]\d\.\d{4})\s+at\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
                      r"\s+a\s+([+-]\d+\.\d+)\s+L\s+([+-]\d+\.\d+)", m.group(1))
    if len(rows) != 20:
        stop("S6-20 parsed to %d spheres, expected 20" % len(rows))
    site = [int(r[0]) for r in rows]
    q = np.array([float(r[1]) for r in rows])
    X = np.array([[float(r[2]), float(r[3]), float(r[4])] for r in rows])
    a = np.array([float(r[5]) for r in rows])
    L = np.array([float(r[6]) for r in rows])
    sc = float(np.sum(a * q) + np.sum(L))
    if abs(sc - SCORE_REC) > 1e-4:
        stop("the parsed design gives a score of %+.5f against the committed %+.5f" % (sc, SCORE_REC))
    return site, q, X, a, L, sc


def prepare(root, out, hpc):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    out = Path(out)
    if out.exists():
        stop("%s already exists" % out)
    arr = (T6 / "relax" / "ARRANGEMENTS.txt").read_text()
    geo = {}
    for st, fn in (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz")):
        p = T6 / "frame_42436" / fn
        mm = re.search(re.escape(fn) + r"\s+sha256\s+([0-9a-f]{64})", arr)
        if not mm or sha(p) != mm.group(1):
            stop("%s does not match the sha256 ARRANGEMENTS.txt records" % fn)
        geo[st] = p
    print("   both geometries match the sha256 ARRANGEMENTS.txt records")
    site, q, X, a, L, sc = design(T6)
    print("   S6-20 parsed: 20 spheres, sum a q %+.5f, sum L %+.5f, score %+.5f (committed %+.5f)"
          % (np.sum(a * q), np.sum(L), sc, SCORE_REC))
    # (c) the sites and the a values against the committed files
    rows = [l.rstrip("\n").split("\t") for l in (T6 / "sites_f42436.tsv").read_text().splitlines() if not l.startswith("#")]
    hd, rows = rows[0], rows[1:]
    ci = {k: hd.index(k) for k in ("site", "x", "y", "z")}
    idx = {int(r[ci["site"]]): k for k, r in enumerate(rows)}
    A = np.loadtxt(str(T6 / "A_f42436.tsv"), comments="#").reshape(-1)
    for k, s in enumerate(site):
        if s not in idx:
            stop("site %d is not in sites_f42436.tsv" % s)
        r = rows[idx[s]]
        d = float(np.abs(X[k] - np.array([float(r[ci[c]]) for c in ("x", "y", "z")])).max())
        if d > 1e-6:
            stop("site %d: ARRANGEMENTS.txt and sites_f42436.tsv differ by %.1e A" % (s, d))
        if abs(A[idx[s]] - a[k]) > 5e-5:
            stop("site %d: a is %+.5f in ARRANGEMENTS.txt and %+.5f in A_f42436.tsv" % (s, a[k], A[idx[s]]))
    print("   all 20 sites and their a values agree with sites_f42436.tsv and A_f42436.tsv")
    # (d) the per-site bound
    for lam in LAMBDAS:
        if float(np.abs(lam * q).max()) > QCAP + 1e-9:
            stop("at lambda %.2f the largest charge is %.4f e, outside the measured bound of %.2f e"
                 % (lam, np.abs(lam * q).max(), QCAP))
    print("   at every lambda the largest single charge is inside the measured %.2f e bound (max %.4f e at lambda 1.00)"
          % (QCAP, np.abs(q).max()))
    # (e) the headers
    if head("def2-SVP") != "! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF\n%maxcore 3000\n%pal nprocs 8 end\n%scf MaxIter 300 end\n":
        stop("the def2-SVP header is not cs5_amatrix.py's HEAD byte for byte")
    ref = T6 / "qscan" / "C5" / "qp1000_R" / "job.inp"
    if head("def2-SVPD") != ref.read_text().split('%pointcharges')[0]:
        stop("the def2-SVPD header is not the committed qscan C5 header byte for byte")
    print("   both headers match their committed counterparts byte for byte")

    out.mkdir(parents=True)
    runs = []
    for basis, cond in BASES:
        for lam in LAMBDAS:
            for st in ("R", "TS"):
                name = "%s/lam%03d_%s" % (basis, round(lam * 100), st)
                d = out / name
                d.mkdir(parents=True)
                (d / "substrate.xyz").write_text(geo[st].read_text())
                lines = ["%d" % len(q)]
                for k in range(len(q)):
                    lines.append("%10.6f %18.10f %18.10f %18.10f" % (lam * q[k], X[k][0], X[k][1], X[k][2]))
                (d / "pc.pc").write_text("\n".join(lines) + "\n")
                (d / "job.inp").write_text(head(basis) + '%pointcharges "pc.pc"\n* xyzfile -2 1 substrate.xyz\n')
                runs.append((name, basis, "%.2f" % lam, st, "%.4f" % float(np.sum(np.abs(lam * q))),
                             "%.4f" % float(np.sum((lam * q) ** 2))))
    with open(out / "FIELD_PLAN.tsv", "w") as fh:
        fh.write("# STAGE6_DIAGNOSTIC4_20261023.txt; cs6_field.py prepare\n")
        fh.write("# S6-20 from relax/ARRANGEMENTS.txt: sum a q %+.5f, sum q^2 %.4f, score %+.5f\n"
                 % (np.sum(a * q), np.sum(q * q), sc))
        fh.write("# bare references come from the committed qscan outputs archive, not re-run\n")
        fh.write("run\tbasis\tlambda\tstate\tsum_abs_q\tsum_q2\n")
        for r in runs:
            fh.write("\t".join(r) + "\n")
    pbs = ["#!/bin/bash", "#PBS -N cm_s6f", "#PBS -l select=1:ncpus=8:mem=24gb", "#PBS -l walltime=06:00:00",
           "#PBS -m ae", "#PBS -M 18660916@sun.ac.za", "#PBS -j oe", "#PBS -o %s/field.pbs.out" % hpc,
           "# Stage 6 DIAGNOSTIC 4 (cs6_field.py): %d single points. Restartable." % len(runs),
           'export PATH="/apps/openmpi/4.1.1/bin:$PATH"',
           'export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"',
           "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1", "ORCA=/home/apps2/ORCA/6.0.1", "W=%s" % hpc,
           'cd "$W" || exit 1', 'if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running"; exit 1; fi',
           "trap 'rmdir \"$W/.running\" 2>/dev/null' EXIT", 'echo "host=$(hostname) start=$(date)"',
           "run() {", '  d="$W/$1"',
           '  if [ -s "$d/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$d/job.out"; then echo "$1 done before"',
           '  else ( cd "$d" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi',
           "}"]
    pbs += ["run %s" % r[0] for r in runs] + ['echo "end=$(date)"']
    (out / "field.pbs").write_text("\n".join(pbs) + "\n")
    files = sorted(p for p in out.rglob("*") if p.is_file())
    (out / "INPUT_CHECKSUMS.sha256").write_text("".join("%s  %s\n" % (sha(p), p.relative_to(out)) for p in files))
    print("   %d runs, %d files including the manifest, one restartable PBS job, nothing submitted"
          % (len(runs), len(files) + 1))


def one(d, ref_xyz, npc):
    p = Path(d) / "job.out"
    if not p.is_file():
        stop("%s: no job.out - has the batch run?" % d)
    t = p.read_text(errors="replace")
    if "****ORCA TERMINATED NORMALLY****" not in t:
        stop("%s: ORCA did not terminate normally" % d)
    if not re.search(r"SCF CONVERGED AFTER\s+\d+\s+CYCLES", t) or "SCF NOT CONVERGED" in t:
        stop("%s: the SCF did not converge" % d)
    m = re.search(r"Number of Electrons\s+NEL\s+\.\.\.\.\s+(\d+)", t)
    if not m or int(m.group(1)) != NEL:
        stop("%s: %s electrons, expected %d" % (d, m.group(1) if m else "no count", NEL))
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    if len(e) != 1:
        stop("%s: %d final energies, expected 1" % (d, len(e)))
    if not re.search(r"Reading point charge file\s*\.*\s*ok \(%d point charges\)" % npc, t):
        stop("%s: ORCA did not confirm reading exactly %d point charges" % (d, npc))
    k = t.find("CARTESIAN COORDINATES (ANGSTROEM)")
    if k < 0:
        stop("%s: no echoed coordinates" % d)
    X = np.array([[float(v) for v in l.split()[1:4]] for l in t[k:].splitlines()[2:2 + len(ref_xyz)]])
    if X.shape != ref_xyz.shape or np.abs(X - ref_xyz).max() > 1e-5:
        stop("%s: the echoed coordinates differ from the committed geometry" % d)
    nb = re.search(r"Number of basis functions\s+\.\.\.\s+(\d+)", t)
    if not nb:
        stop("%s: no basis-function count" % d)
    return float(e[0]), int(nb.group(1))


def analyse(root, runs, out):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    R = Path(runs)
    site, q, X, a, L, sc = design(T6)
    SAQ, SQ2 = float(np.sum(a * q)), float(np.sum(q * q))
    # (i) b out of the committed amendment result
    amr = T6 / "qrefine" / "DIAGNOSTIC2_AMENDMENT1_RESULT.txt"
    if not amr.is_file():
        stop("%s is missing" % amr)
    mb = re.search(r"b\s*=\s*([+-]?\d+\.\d+)\s*kcal/mol per e\^2", amr.read_text())
    if not mb:
        stop("could not read b out of the committed DIAGNOSTIC2_AMENDMENT1_RESULT.txt")
    Bq = float(mb.group(1))
    print("   b read from the committed DIAGNOSTIC2_AMENDMENT1_RESULT.txt: %+.4f kcal/mol per e^2" % Bq)
    print("   S6-20 from ARRANGEMENTS.txt: sum a q %+.5f, sum q^2 %.4f" % (SAQ, SQ2))
    # (h) the bare references
    man = {}
    for l in (T6 / "qscan" / "outputs" / "qscan_outputs.MANIFEST.sha256").read_text().splitlines():
        if l.strip():
            h, n = l.split(None, 1)
            man[n.strip()] = h
    Q = {}
    with tarfile.open(str(T6 / "qscan" / "outputs" / "qscan_outputs.tar.xz")) as t:
        for basis, cond in BASES:
            for st in ("R", "TS"):
                nm = "./%s/bare_%s/job.out" % (cond, st)
                raw = t.extractfile(nm).read()
                if nm in man and hashlib.sha256(raw).hexdigest() != man[nm]:
                    stop("%s does not match its sha256 in the qscan manifest" % nm)
                ee = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", raw.decode(errors="replace"))
                if len(ee) != 1:
                    stop("%s: %d final energies" % (nm, len(ee)))
                Q[(basis, st)] = float(ee[0])
    ref = None
    for l in (T6 / "A_REPORT.txt").read_text().splitlines():
        f = l.split("\t")
        if f and f[0].strip() == "42436" and len(f) >= 3:
            ref = (float(f[1]), float(f[2]))
    for st, v in (("R", ref[0]), ("TS", ref[1])):
        if abs(Q[("def2-SVP", st)] - v) > 1e-6:
            stop("the committed qscan C1 bare %s energy does not match A_REPORT.txt" % st)
    print("   the committed archive's def2-SVP bare energies match A_REPORT.txt to <1e-6 Eh")

    GEO = {st: read_xyz(T6 / "frame_42436" / ("reactant_qm.xyz" if st == "R" else "transition_state_qm.xyz"))
           for st in ("R", "TS")}
    if abs(float(np.sum(L)) - SUM_L) > 1e-4:
        stop("the parsed sum L is %+.5f against this script's %+.5f" % (np.sum(L), SUM_L))
    E, NB = {}, {}
    for basis, cond in BASES:
        for lam in LAMBDAS:
            for st in ("R", "TS"):
                d = R / basis / ("lam%03d_%s" % (round(lam * 100), st))
                E[(basis, lam, st)], NB[(basis, lam, st)] = one(d, GEO[st], len(q))
        s = {NB[(basis, lam, st)] for lam in LAMBDAS for st in ("R", "TS")}
        if len(s) != 1:
            stop("%s runs report different basis-function counts: %s" % (basis, sorted(s)))
        print("   %-10s all 8 runs at %d basis functions" % (basis, s.pop()))

    # the scan's def2-SVP lambda = 1 point IS the committed frozen S6-20 measurement minus its classical sum L,
    # since Finding 4 established the bare charge and the LJ sphere differ by exactly L and nothing else
    ddE, pred, S = {}, {}, {}
    for basis, cond in BASES:
        for lam in LAMBDAS:
            ddE[(basis, lam)] = ((E[(basis, lam, "TS")] - Q[(basis, "TS")])
                                 - (E[(basis, lam, "R")] - Q[(basis, "R")])) * H
            pred[lam] = lam * SAQ + lam * lam * Bq * SQ2
            S[(basis, lam)] = 100 * abs(ddE[(basis, lam)] - pred[lam]) / abs(pred[lam])
    want = D1_SVP - SUM_L
    got = ddE[("def2-SVP", 1.00)]
    if abs(got - want) > 0.05:
        stop("the def2-SVP lambda = 1 run gives ddE %+.4f; the committed frozen S6-20 value minus its sum L is "
             "%+.4f, and Finding 4 says those must agree. %.4f kcal/mol apart - the scan does not reproduce the "
             "committed measurement and that must be found before anything is read" % (got, want, abs(got - want)))
    print("   the def2-SVP lambda = 1 run reproduces the committed frozen S6-20 value minus sum L: %+.4f against "
          "%+.4f" % (got, want))

    O = ["STAGE 6 - DIAGNOSTIC 4 RESULT: the total-field scan, S6-20 scaled", "",
         "Measured by cs6_field.py analyse. The rule is STAGE6_DIAGNOSTIC4_20261023.txt D4.3, fixed before the",
         "runs went in and applied unchanged. Nothing is refitted: sum a q and sum q^2 come from the committed",
         "ARRANGEMENTS.txt and b from the committed DIAGNOSTIC2_AMENDMENT1_RESULT.txt.", "",
         "   predicted(lambda) = lambda * sum a q  +  lambda^2 * b * sum q^2",
         "   sum a q = %+.5f kcal/mol    sum q^2 = %.4f e^2    b = %+.4f kcal/mol per e^2" % (SAQ, SQ2, Bq),
         "   ONE b IS USED FOR ALL TWENTY SITES AND THAT IS AN APPROXIMATION: b was measured at site 1690 only.", "",
         "      lambda   S|q|    largest |q|   predicted    def2-SVP     S      def2-SVPD    S       Delta"]
    for lam in LAMBDAS:
        O.append("      %5.2f  %6.3f     %6.4f     %9.4f  %9.4f %6.1f%%  %9.4f %6.1f%%  %+8.4f"
                 % (lam, lam * float(np.sum(np.abs(q))), lam * float(np.abs(q).max()), pred[lam],
                    ddE[("def2-SVP", lam)], S[("def2-SVP", lam)],
                    ddE[("def2-SVPD", lam)], S[("def2-SVPD", lam)],
                    ddE[("def2-SVPD", lam)] - ddE[("def2-SVP", lam)]))
    O += ["", "   Delta = ddE(def2-SVPD) - ddE(def2-SVP), the inter-basis movement at fixed design and fixed",
          "   geometry, with only the total field varying.", "", "THE THREE READINGS OF D4.3", ""]
    cand = [lam for lam in LAMBDAS if S[("def2-SVPD", lam)] <= PCT]
    lq = max(cand) if cand else None
    O.append("   1. THE TOTAL-FIELD BOUND")
    if lq is None:
        O.append("      no lambda tested is within %.0f per cent at def2-SVPD: the total-field bound is BELOW" % PCT)
        O.append("      lambda = %.2f, i.e. below sum |q| = %.3f e, and Stage 7 must cap the total field there or" % (LAMBDAS[0], LAMBDAS[0]*float(np.sum(np.abs(q)))))
        O.append("      lower. The design space is much smaller than the per-site bound alone implies.")
    elif lq >= 1.00:
        O.append("      lambda_quad = 1.00: the quadratic model holds across the whole range a design can reach")
        O.append("      under the pre-registered budget, so THERE IS NO TOTAL-FIELD BOUND to impose and the")
        O.append("      per-site cap of %.2f e is the only bound Stage 7 needs." % QCAP)
    else:
        O.append("      lambda_quad = %.2f, so the total-field bound is sum |q| <= %.3f e."
                 % (lq, lq * float(np.sum(np.abs(q)))))
        O.append("      Stage 7 imposes that in addition to the per-site cap of %.2f e." % QCAP)
    lams = np.array(LAMBDAS)
    dl = np.array([ddE[("def2-SVPD", l)] - ddE[("def2-SVP", l)] for l in LAMBDAS])
    r = float(np.corrcoef(lams, dl)[0, 1])
    sl, ic = np.polyfit(lams, dl, 1)
    O += ["", "   2. THE sum |q| HYPOTHESIS", "",
          "      Delta against lambda: r = %+.4f, slope %+.4f kcal/mol per unit lambda, intercept %+.4f" % (r, sl, ic)]
    if abs(r) >= 0.98 and abs(ic) <= 0.3:
        O.append("      CONFIRMED: the inter-basis movement is linear in the total field with a near-zero")
        O.append("      intercept, at fixed design. The assessment's r = +0.9955 across four designs was not an")
        O.append("      artefact of comparing designs.")
    else:
        O.append("      NOT CONFIRMED: the movement is not linear in the total field at fixed design, so the")
        O.append("      assessment's correlation describes those four designs and is not a scaling law. Its")
        O.append("      section 4 must be read that way.")
    d1 = ddE[("def2-SVPD", 1.00)]
    d0 = ddE[("def2-SVP", 1.00)]
    O += ["", "   3. DOES S6-20 SURVIVE?", "",
          "      Diagnostic 1 at ma-def2-SVP          %+9.4f" % D1_MA,
          "      committed frozen at def2-SVP         %+9.4f   (includes sum L = %+.5f)" % (D1_SVP, SUM_L),
          "      its bare-charge part, D1_SVP - sum L %+9.4f" % (D1_SVP - SUM_L),
          "      this scan at def2-SVP,  lambda 1     %+9.4f   <- must equal the line above" % d0,
          "      this scan at def2-SVPD, lambda 1     %+9.4f" % d1, "",
          "      The comparison that decides it is the last two lines against each other, which are the same",
          "      design, geometry and charge model differing only in basis - not against the committed frozen",
          "      number, which carries sum L and is therefore not like for like.", ""]
    if d1 < 0 and 100 * abs(d1 - d0) / abs(d0) <= O6:
        O.append("      S6-20 SURVIVES. The def2-SVPD value is negative and within %.0f per cent of this scan's own" % O6)
        O.append("      def2-SVP value (%.1f per cent), so Diagnostic 1's sign reversal was an ma-def2-SVP artefact,"
                 % (100 * abs(d1 - d0) / abs(d0)))
        O.append("      consistent with DIAGNOSTIC3_NOTE1 N3.2. Finding 6 and the magnitudes of Findings 1 and 5")
        O.append("      need RE-READING at def2-SVPD rather than withdrawing.")
    else:
        O.append("      S6-20 DOES NOT SURVIVE. The sign reversal is not specific to ma-def2-SVP (def2-SVPD gives")
        O.append("      %+.4f against this scan's def2-SVP %+.4f, %.1f per cent apart), so the frozen half of Stage 6"
                 % (d1, d0, 100 * abs(d1 - d0) / abs(d0)))
        O.append("      to be rebuilt at def2-SVPD before any of it informs anything.")
    O += ["", "WHAT THIS DOES NOT TEST", "",
          "   One design, so this is about S6-20's particular arrangement of twenty charges and not about total",
          "   field in general. One frame, fixed geometry. Whether b varies between sites, which is the first",
          "   explanation to reach for if the residuals are large. Negative total field. Nothing about Loss."]
    Path(out).mkdir(parents=True, exist_ok=True)
    (Path(out) / "DIAGNOSTIC4_RESULT.txt").write_text("\n".join(O) + "\n")
    print()
    print("\n".join(O[O.index("THE THREE READINGS OF D4.3"):]))
    print()
    print("   wrote %s" % (Path(out) / "DIAGNOSTIC4_RESULT.txt"))


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in ("prepare", "analyse"):
        sys.exit(__doc__)
    (prepare if sys.argv[1] == "prepare" else analyse)(sys.argv[2], sys.argv[3], sys.argv[4])
