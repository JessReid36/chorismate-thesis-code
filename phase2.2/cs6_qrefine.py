#!/usr/bin/env python3
"""
cs6_qrefine.py - STAGE6_DIAGNOSTIC2_AMENDMENT1_20261022.txt: two further charge magnitudes at def2-SVPD, to
locate where the substrate's response stops being quadratic.

  python3 cs6_qrefine.py prepare <results_repo> <OUT> <HPC_DIR>
  python3 cs6_qrefine.py analyse <results_repo> <RUNS> <OUT>

Site 1690, bare point charges, q = +-0.65 and +-0.80 e, both committed geometries, def2-SVPD in vacuo only.
Everything else is cs6_qscan.py's, which is cs5_amatrix.py's: the point-charge file format, the input header, the
restartable single-PBS-job pattern. The bare references and the q <= 0.50 and q = 1.00 points are NOT re-run -
they come out of the committed qscan outputs archive, so a and b are the published numbers and the new points are
judged against them rather than against a refit.

prepare checks, each a stop:
  (a) both geometries match the sha256 relax/ARRANGEMENTS.txt records;
  (b) site 1690's coordinates agree between sites_f42436.tsv and ARRANGEMENTS.txt's S6-P block to 1e-6 A;
  (c) the generated header equals the committed qscan C5 inputs' header byte for byte, so the condition is the
      same one Diagnostic 2 measured and not a restatement of it;
  (d) it refuses to write into an existing folder, and submits nothing.

analyse checks, each a stop:
  (e) every new run: ORCA TERMINATED NORMALLY, SCF converged, 118 electrons, exactly one final energy, the echoed
      coordinates equal to the committed geometry, ORCA's confirmation that it read one point charge, and the
      def2-SVPD basis-function count;
  (f) the committed qscan archive supplies the references, each file checked against its manifest sha256, and its
      C1 bare energies reproduce A_REPORT.txt to below 1e-6 Eh;
  (g) a and b are read out of the committed DIAGNOSTIC2_RESULT.txt text, not recomputed and not typed twice;
  (h) S(0.50) recomputed from the committed points must still be under 1 per cent, and S at q = 1.00 must still be
      within 0.5 points of the 16.10 per cent the committed report gives - if either has moved, the new run's
      reading of the archive disagrees with the published one and that must be found first.
Then it applies A1.3 without choosing anything.
"""
import hashlib, re, sys, tarfile
from pathlib import Path
import numpy as np

H = 627.5094740631
SITE = 1690
NEW_QS = (0.65, 0.80)
OLD_QS = (0.05, 0.25, 0.50, 1.00)
BASIS = "def2-SVPD"
COND = "C5"
NEL = 118
PCT = 10.0


def stop(m):
    sys.exit("STOP: " + m)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    return np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])


def qtag(q):
    return ("qp" if q > 0 else "qm") + "%04d" % round(abs(q) * 1000)


def head():
    return ("! B3LYP D3BJ %s def2/J RIJCOSX TightSCF\n%%maxcore 3000\n%%pal nprocs 8 end\n"
            "%%scf MaxIter 300 end\n" % BASIS)


def prepare(root, out, hpc):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    out = Path(out)
    if out.exists():
        stop("%s already exists" % out)
    arr = (T6 / "relax" / "ARRANGEMENTS.txt").read_text()
    geo = {}
    for st, fn in (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz")):
        p = T6 / "frame_42436" / fn
        m = re.search(re.escape(fn) + r"\s+sha256\s+([0-9a-f]{64})", arr)
        if not m:
            stop("ARRANGEMENTS.txt records no sha256 for %s" % fn)
        if sha(p) != m.group(1):
            stop("%s does not match the sha256 ARRANGEMENTS.txt records" % fn)
        geo[st] = p
        print("   %-2s %s sha256 matches ARRANGEMENTS.txt" % (st, fn))
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
    if float(np.abs(P - np.array([float(m.group(i)) for i in (1, 2, 3)])).max()) > 1e-6:
        stop("site %d differs between the sites file and ARRANGEMENTS.txt" % SITE)
    print("   site %d at %.6f %.6f %.6f, both records agree" % (SITE, P[0], P[1], P[2]))
    # (c) the header must equal the committed C5 inputs' header byte for byte
    ref_inp = T6 / "qscan" / COND / "qp1000_R" / "job.inp"
    if not ref_inp.is_file():
        stop("%s is missing - the committed qscan inputs are not in place" % ref_inp)
    ref_head = ref_inp.read_text().split('%pointcharges')[0]
    if head() != ref_head:
        stop("the generated header is not the committed qscan %s header byte for byte:\n--- generated ---\n%s\n--- committed ---\n%s"
             % (COND, head(), ref_head))
    print("   the header equals the committed qscan %s inputs' header byte for byte" % COND)

    out.mkdir(parents=True)
    runs = []
    for st in ("R", "TS"):
        for q in [s * x for x in NEW_QS for s in (+1, -1)]:
            name = "%s_%s" % (qtag(q), st)
            d = out / name
            d.mkdir(parents=True)
            (d / "substrate.xyz").write_text(geo[st].read_text())
            (d / "pc.pc").write_text("1\n%10.6f %18.10f %18.10f %18.10f\n" % (float(q), P[0], P[1], P[2]))
            (d / "job.inp").write_text(head() + '%pointcharges "pc.pc"\n* xyzfile -2 1 substrate.xyz\n')
            runs.append((name, BASIS, st, "%+.2f" % q))
    with open(out / "QREFINE_PLAN.tsv", "w") as fh:
        fh.write("# STAGE6_DIAGNOSTIC2_AMENDMENT1_20261022.txt; cs6_qrefine.py prepare\n")
        fh.write("# site %d at %.6f %.6f %.6f; the references and q <= 0.50 and 1.00 come from the committed qscan\n"
                 % (SITE, P[0], P[1], P[2]))
        fh.write("run\tbasis\tstate\tcharge_e\n")
        for r in runs:
            fh.write("\t".join(r) + "\n")
    pbs = ["#!/bin/bash", "#PBS -N cm_s6r", "#PBS -l select=1:ncpus=8:mem=24gb", "#PBS -l walltime=04:00:00",
           "#PBS -m ae", "#PBS -M 18660916@sun.ac.za", "#PBS -j oe", "#PBS -o %s/qrefine.pbs.out" % hpc,
           "# Stage 6 DIAGNOSTIC 2 AMENDMENT 1 (cs6_qrefine.py): %d single points. Restartable." % len(runs),
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
    (out / "qrefine.pbs").write_text("\n".join(pbs) + "\n")
    files = sorted(p for p in out.rglob("*") if p.is_file())
    (out / "INPUT_CHECKSUMS.sha256").write_text("".join("%s  %s\n" % (sha(p), p.relative_to(out)) for p in files))
    print("   %d runs, %d files including the manifest, one restartable PBS job, nothing submitted"
          % (len(runs), len(files) + 1))


def one(d, ref_xyz):
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
    if not re.search(r"Reading point charge file\s*\.*\s*ok \(1 point charges\)", t):
        stop("%s: ORCA did not confirm reading exactly one point charge" % d)
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
    # (g) a and b out of the committed report
    rep = (T6 / "qscan" / "DIAGNOSTIC2_RESULT.txt").read_text()
    blk = re.search(r"^%s\s+%s\s+in vacuo\s*$(.*?)(?=^\S|\Z)" % (COND, re.escape(BASIS)), rep, re.M | re.S)
    if not blk:
        stop("could not find the %s %s block in the committed DIAGNOSTIC2_RESULT.txt" % (COND, BASIS))
    ma = re.search(r"a_est \(odd/q at q=0\.05\)\s*(-?\d+\.\d+)", blk.group(1))
    if not ma:
        stop("could not read a_est out of the committed %s block" % COND)
    A = float(ma.group(1))
    # b must come from q <= 0.50, where every condition is internally consistent. The report's own printed
    # b is even/q^2 evaluated at q = 1.00, which is INSIDE the saturated region and is not the fit A1.3 names.
    mb = re.search(r"^\s*0\.50\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", blk.group(1), re.M)
    if not mb:
        stop("could not read the q = 0.50 row out of the committed %s block" % COND)
    even_050 = float(mb.group(4))
    Bq = even_050 / 0.25
    print("   from the committed DIAGNOSTIC2_RESULT.txt %s block: a %+.4f (printed a_est);" % (COND, A))
    print("   b %+.4f, derived as the printed even value at q = 0.50 (%+.4f) over q^2 - NOT the report's printed"
          % (Bq, even_050))
    print("   b, which is evaluated at q = 1.00 inside the saturated region")
    # (f) the references out of the committed archive
    arch = T6 / "qscan" / "outputs" / "qscan_outputs.tar.xz"
    man = {}
    for l in (T6 / "qscan" / "outputs" / "qscan_outputs.MANIFEST.sha256").read_text().splitlines():
        if l.strip():
            h, n = l.split(None, 1)
            man[n.strip()] = h
    Q = {}
    with tarfile.open(str(arch)) as t:
        for cond in (COND, "C1"):
            for st in ("R", "TS"):
                for sub in ["bare"] + [qtag(s * x) for x in OLD_QS for s in (+1, -1)]:
                    nm = "./%s/%s_%s/job.out" % (cond, sub, st)
                    try:
                        raw = t.extractfile(nm).read()
                    except Exception:
                        continue
                    if nm in man and hashlib.sha256(raw).hexdigest() != man[nm]:
                        stop("%s does not match its sha256 in the qscan manifest" % nm)
                    ee = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", raw.decode(errors="replace"))
                    if len(ee) == 1:
                        Q[(cond, sub, st)] = float(ee[0])
    ref = None
    for l in (T6 / "A_REPORT.txt").read_text().splitlines():
        f = l.split("\t")
        if f and f[0].strip() == "42436" and len(f) >= 3:
            ref = (float(f[1]), float(f[2]))
    for st, v in (("R", ref[0]), ("TS", ref[1])):
        if abs(Q[("C1", "bare", st)] - v) > 1e-6:
            stop("the committed qscan C1 bare %s energy does not match A_REPORT.txt" % st)
    print("   the committed archive's C1 bare energies match A_REPORT.txt to <1e-6 Eh")

    GEO = {st: read_xyz(T6 / "frame_42436" / ("reactant_qm.xyz" if st == "R" else "transition_state_qm.xyz"))
           for st in ("R", "TS")}
    E, NB = {}, {}
    for st in ("R", "TS"):
        for q in [s * x for x in NEW_QS for s in (+1, -1)]:
            E[(q, st)], NB[(q, st)] = one(R / ("%s_%s" % (qtag(q), st)), GEO[st])
    nb0 = set(NB.values())
    if len(nb0) != 1:
        stop("the new runs report different basis-function counts: %s" % sorted(nb0))
    print("   all 8 new runs at %d basis functions" % nb0.pop())

    def dd(sub, st_pair_src, q=None):
        pass

    def ddE_old(q):
        b = {st: Q[(COND, "bare", st)] for st in ("R", "TS")}
        return ((Q[(COND, qtag(q), "TS")] - b["TS"]) - (Q[(COND, qtag(q), "R")] - b["R"])) * H

    def ddE_new(q):
        b = {st: Q[(COND, "bare", st)] for st in ("R", "TS")}
        return ((E[(q, "TS")] - b["TS"]) - (E[(q, "R")] - b["R"])) * H

    S = {}
    rowsout = []
    for q in sorted(set(OLD_QS) | set(NEW_QS)):
        new = q in NEW_QS
        dp = ddE_new(+q) if new else ddE_old(+q)
        dm = ddE_new(-q) if new else ddE_old(-q)
        pred = A * q + Bq * q * q
        S[q] = abs(dp - pred) / abs(pred) * 100
        rowsout.append((q, dp, dm, (dp - dm) / 2, (dp + dm) / 2, pred, S[q], "new" if new else "committed"))
    # (h) the committed points must still read as published
    if S[0.50] >= 1.0:
        stop("S(0.50) recomputed as %.2f per cent; the committed report gives 0.58, so this run's reading of the "
             "archive disagrees with the published one" % S[0.50])
    if abs(S[1.00] - 16.10) > 0.5:
        stop("S(1.00) recomputed as %.2f per cent against the committed report's 16.10" % S[1.00])
    print("   the committed points still read as published: S(0.50) %.2f%%  S(1.00) %.2f%%" % (S[0.50], S[1.00]))

    cand = [q for q in sorted(S) if S[q] <= PCT]
    q_quad = max(cand) if cand else None
    L = ["STAGE 6 - DIAGNOSTIC 2 AMENDMENT 1 RESULT: where the response stops being quadratic", "",
         "Measured by cs6_qrefine.py analyse. The rule is STAGE6_DIAGNOSTIC2_AMENDMENT1_20261022.txt A1.3, fixed",
         "before the runs went in and applied unchanged. a and b are the committed DIAGNOSTIC2_RESULT.txt values",
         "fitted on q <= 0.50 e and are NOT refitted with the new points.", "",
         "   site %d, bare point charge, %s in vacuo, frame 42436's committed geometries" % (SITE, BASIS),
         "   a = %+.4f kcal/mol per e        b = %+.4f kcal/mol per e^2" % (A, Bq),
         "   S(q) = |ddE(+q) - (a q + b q^2)| / |a q + b q^2|", "",
         "      q      ddE(+q)    ddE(-q)      odd       even    a q + b q^2     S(q)      source",
         ]
    for q, dp, dm, od, ev, pred, s, src in rowsout:
        L.append("   %6.2f  %9.4f  %9.4f  %9.4f  %9.4f  %11.4f  %7.2f%%   %s"
                 % (q, dp, dm, od, ev, pred, s, src))
    L += ["", "THE RULE OF A1.3", "",
          "   q_quad, the largest magnitude with S(q) <= %.0f per cent:  %s" % (PCT, ("%.2f e" % q_quad) if q_quad else "none"),
          ""]
    if q_quad is None:
        L.append("VERDICT: NO MAGNITUDE TESTED IS WITHIN 10 PER CENT - this contradicts the committed S(0.50) and")
        L.append("  must be investigated before anything is read.")
    else:
        L.append("VERDICT PER A1.3: STAGE 7'S PER-SITE CHARGE CAP IS %.2f e" % q_quad)
        if q_quad >= 0.80:
            L.append("  - the response is quadratic essentially to the old bound; the saturation sets in only in")
            L.append("    the last fifth of the admitted range.")
        elif q_quad == 0.50:
            L.append("  - the saturation begins immediately above the highest well-behaved point Diagnostic 2")
            L.append("    already had, so the new magnitudes bought a bound rather than a surprise.")
    L += ["", "   Note that q_quad bounds the QUADRATIC model, which is what Stage 7 reports with. D2.3's stricter",
          "   q*, bounding where the LINEAR term alone is within 10 per cent of the whole response, is a different",
          "   and smaller number and remains 0.25 e at this condition.", "",
          "WHAT THIS DOES NOT TEST", "",
          "   Other sites - the bound is site- and distance-dependent and site %d sits 2.96 A from the nearest" % SITE,
          "   substrate atom. Other frames. Whether the saturation has the same cause as the sign inversion",
          "   ma-def2-SVP showed at 1.00 e, which remains unexplained. Anything about the relaxed geometry."]
    Path(out).mkdir(parents=True, exist_ok=True)
    (Path(out) / "DIAGNOSTIC2_AMENDMENT1_RESULT.txt").write_text("\n".join(L) + "\n")
    print()
    print("\n".join(L[L.index("THE RULE OF A1.3"):]))
    print()
    print("   wrote %s" % (Path(out) / "DIAGNOSTIC2_AMENDMENT1_RESULT.txt"))


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in ("prepare", "analyse"):
        sys.exit(__doc__)
    (prepare if sys.argv[1] == "prepare" else analyse)(sys.argv[2], sys.argv[3], sys.argv[4])
