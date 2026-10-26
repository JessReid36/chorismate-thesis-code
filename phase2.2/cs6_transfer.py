#!/usr/bin/env python3
"""
cs6_transfer.py - STAGE6_DIAGNOSTIC7_20261026.txt: does the defensible domain transfer? A SECOND design,
scaled, against the total fields S6-20 was measured at.

  python3 cs6_transfer.py prepare <results_repo> <DESIGN> <N> <OUT> <HPC_DIR>
  python3 cs6_transfer.py analyse <results_repo> <DESIGN> <N> <RUNS> <OUT>

cs6_field.py IS NOT MODIFIED AND IS NOT IMPORTED. HANDOVER_20261026.md section 10 names it the committed
generator of record for the D4 inputs, checked byte for byte by the D4 result bundle. This script is a separate
generator for Diagnostic 7 and carries cs6_field.py's stops across in substance, with these differences, each
required by the pre-registration:

  - the design is an ARGUMENT, not S6-20. Nothing design-specific is typed: the committed score, sum L, sum q^2
    and the sites all come out of relax/ARRANGEMENTS.txt, and the frozen ddE for the consistency stop comes out
    of frozen/frozen.tsv at full precision.
  - the lambda grid is D4's four plus TWO MATCHED values, computed as (S6-20's sum |q| at lambda 0.25 and 0.50)
    divided by this design's sum |q|. S6-20's own design is parsed for that, so the targets 1.440 and 2.880 e
    are derived and not typed either.
  - the per-site check does NOT use the 0.80 e cap. STAGE6_ERRATUM3_20261026.txt E3.1 withdrew it. The check is
    that the design's largest charge at lambda = 1 is at or below 0.65 e, which E3.1 records as the largest
    single magnitude CLEAN at def2-SVPD at site 1690 at 2.96 A. That stop REFUSES any of the four 1.0000 e
    designs, which is intended.
  - every run is classified by DIAGNOSTIC 6's rule, which D4 could not do because D6 did not exist. D6's parser
    cs6_audit.read() and all five of its thresholds are IMPORTED from cs6_audit, never retyped. Its verdict()
    is a closure inside its main() and cannot be imported, so the five comparisons are transcribed here
    structurally unchanged, against the control medians published in the committed DIAGNOSTIC6_RESULT.txt.

prepare checks, each a stop:
  (a) both geometries match the sha256 relax/ARRANGEMENTS.txt records;
  (b) the design parses to exactly N spheres and sum a q + sum L reproduces its own committed score to 1e-4;
  (c) every site's coordinates agree with sites_f42436.tsv to 1e-6 A and every a with A_f42436.tsv to 5e-5;
  (d) the design's largest charge is at or below E3.1's 0.65 e, and at every lambda it is at or below the
      design's own lambda = 1 maximum;
  (e) the def2-SVP header equals cs5_amatrix.py's HEAD byte for byte and the def2-SVPD header equals the
      committed qscan C5 header byte for byte;
  (f) it refuses to write into an existing folder, and submits nothing.

analyse checks, each a stop:
  (g) every run: ORCA TERMINATED NORMALLY, SCF converged, 118 electrons, one final energy, the echoed
      coordinates equal to the committed geometry, ORCA's confirmation that it read exactly N point charges,
      one basis-function count per basis;
  (h) the committed qscan archive supplies the bare references, each checked against its manifest sha256, and
      its C1 bare energies reproduce A_REPORT.txt to below 1e-6 Eh;
  (i) b is read out of the committed DIAGNOSTIC2_AMENDMENT1 result's own text, not typed here;
  (j) S6-20's measured Delta at the two matched total fields is read out of the committed DIAGNOSTIC4_RESULT.txt
      table, not typed here;
  (k) the consistency stop of D7.3: the def2-SVP lambda = 1 point must equal the committed frozen ddE minus
      sum L within 0.05 kcal/mol, or nothing else is read.
Then it applies D7.3's four readings without choosing anything.
"""
import hashlib, re, statistics as stats, sys, tarfile
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cs6_audit as AUD          # read() and every D6 threshold: imported, never retyped

H = 627.5094740631
NEL = 118
BASES = (("def2-SVP", "C1"), ("def2-SVPD", "C5"))
BASE_LAMS = (0.25, 0.50, 0.75, 1.00)      # D4's grid, so the two scans compare point for point
REF_DESIGN, REF_N = "S6-20", 20           # the design whose total fields the matched lambdas target
REF_LAMS = (0.25, 0.50)                   # its lambda 0.25 and 0.50, i.e. sum |q| 1.440 and 2.880 e
PCT = 10.0            # Diagnostic 2's convention, unchanged (D7.3 reading 2)
AGREE_KCAL = 1.0      # E3.3 criterion (i): the two bases agree on ddE to within this
DELTA_TOL = 0.3       # D4.3's intercept tolerance, carried over BY CHOICE for reading 3; D7.3 discloses that
STOP_TOL = 0.05       # D4.3's own consistency-stop tolerance
QCLEAN = 0.65         # E3.1: largest single magnitude CLEAN at def2-SVPD, site 1690, 2.96 A


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


def design(T6, name, n):
    """One design out of the committed ARRANGEMENTS.txt. Its score and sum q^2 are READ from its own block."""
    arr = (T6 / "relax" / "ARRANGEMENTS.txt").read_text()
    m = re.search(r"^\s*%s\s+%d sphere\(s\)(.*?)(?=^\s*S6-\S+\s+\d+ sphere\(s\)|\Z)" % (re.escape(name), n),
                  arr, re.M | re.S)
    if not m:
        stop("could not find the %s block (%d spheres) in the committed ARRANGEMENTS.txt" % (name, n))
    body = m.group(1)
    ms = re.search(r"score\s+([+-]?\d+\.\d+)\s+kcal/mol", body)
    mq = re.search(r"sum q\^2\s+(\d+\.\d+)", body)
    if not ms or not mq:
        stop("%s: could not read its committed score and sum q^2 out of ARRANGEMENTS.txt" % name)
    score_rec, sq2_rec = float(ms.group(1)), float(mq.group(1))
    rows = re.findall(r"site\s+(\d+)\s+q\s+([+-]\d\.\d{4})\s+at\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
                      r"\s+a\s+([+-]\d+\.\d+)\s+L\s+([+-]\d+\.\d+)", body)
    if len(rows) != n:
        stop("%s parsed to %d spheres, expected %d" % (name, len(rows), n))
    site = [int(r[0]) for r in rows]
    q = np.array([float(r[1]) for r in rows])
    X = np.array([[float(r[2]), float(r[3]), float(r[4])] for r in rows])
    a = np.array([float(r[5]) for r in rows])
    L = np.array([float(r[6]) for r in rows])
    sc = float(np.sum(a * q) + np.sum(L))
    if abs(sc - score_rec) > 1e-4:
        stop("%s: the parsed design gives a score of %+.5f against its committed %+.5f"
             % (name, sc, score_rec))
    if abs(float(np.sum(q * q)) - sq2_rec) > 1e-4:
        stop("%s: the parsed sum q^2 is %.5f against its committed %.5f"
             % (name, np.sum(q * q), sq2_rec))
    return dict(name=name, site=site, q=q, X=X, a=a, L=L, score=score_rec,
                saq=float(np.sum(a * q)), sq2=float(np.sum(q * q)),
                sabs=float(np.sum(np.abs(q))), sumL=float(np.sum(L)), maxq=float(np.abs(q).max()))


def lambdas(D, REF):
    """D4's four, plus the two matched to the reference design's total field at REF_LAMS. Nothing typed."""
    out = {}
    for lam in BASE_LAMS:
        out["lam%03d" % round(lam * 100)] = lam
    for k, rl in enumerate(REF_LAMS, start=1):
        target = rl * REF["sabs"]
        lam = target / D["sabs"]
        if lam > 1.0 + 1e-12:
            stop("matching %s's sum |q| of %.4f e would need lambda %.6f on %s, above 1.00; D7.2 rejects that"
                 % (REF["name"], target, lam, D["name"]))
        out["lamM%d" % k] = lam
    return dict(sorted(out.items(), key=lambda kv: kv[1]))


def controls(T6):
    """D6's own published control medians, per basis-set size, out of the committed DIAGNOSTIC6_RESULT.txt."""
    p = T6 / "DIAGNOSTIC6_RESULT.txt"
    if not p.is_file():
        stop("%s is missing: reading 4 cannot be anchored" % p)
    t = p.read_text()
    k = t.find("THE CONTROLS EVERY THRESHOLD IS ANCHORED TO")
    if k < 0:
        stop("DIAGNOSTIC6_RESULT.txt does not carry its control table")
    C = {}
    for l in t[k:].splitlines():
        m = re.match(r"\s*(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$", l)
        if m:
            C[int(m.group(1))] = dict(basis=m.group(2), n=int(m.group(3)), cyc=int(m.group(4)),
                                      gap=float(m.group(5)), maxlw=float(m.group(6)), dip=float(m.group(7)))
        elif C and not l.strip():
            break
    if not C:
        stop("could not parse D6's control table")
    return C


def verdict(d, C):
    """A structural transcription of cs6_audit.py main()'s verdict(), which is a closure and cannot be imported.
    Every threshold is AUD's, imported above. Do not re-tune anything here."""
    c = C.get(d["nb"])
    if not c:
        return "NO CONTROL", []
    moved, sick = [], []
    if d["maxlw"] is not None and d["maxlw"] - c["maxlw"] >= AUD.D_LW:
        moved.append("Loewdin +%.3f e over control" % (d["maxlw"] - c["maxlw"]))
    if d["dip"] is not None and c["dip"] and d["dip"] / c["dip"] >= AUD.R_DIP:
        moved.append("dipole %.1fx control" % (d["dip"] / c["dip"]))
    if d["gap"] is not None and d["gap"] < AUD.GAP_ABS and d["gap"] < AUD.R_GAP * c["gap"]:
        sick.append("gap %.3f eV" % d["gap"])
    if d["cyc"] is not None and c["cyc"] and d["cyc"] >= AUD.R_CYC * c["cyc"]:
        sick.append("%d SCF cycles, %.1fx control" % (d["cyc"], d["cyc"] / c["cyc"]))
    if moved and sick:
        return "COLLAPSED", moved + sick
    if sick:
        return "STRAINED", sick
    if moved:
        return "POLARISED", moved
    return "CLEAN", []


def prepare(root, name, n, out, hpc):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    out = Path(out)
    if out.exists():
        stop("%s already exists" % out)
    arr = (T6 / "relax" / "ARRANGEMENTS.txt").read_text()
    # (a) the geometries
    geo = {}
    for st, fn in (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz")):
        p = T6 / "frame_42436" / fn
        mm = re.search(re.escape(fn) + r"\s+sha256\s+([0-9a-f]{64})", arr)
        if not mm or sha(p) != mm.group(1):
            stop("%s does not match the sha256 ARRANGEMENTS.txt records" % fn)
        geo[st] = p
    print("   both geometries match the sha256 ARRANGEMENTS.txt records")
    # (b) this design and the reference design
    D = design(T6, name, n)
    REF = design(T6, REF_DESIGN, REF_N)
    print("   %s parsed: %d spheres, sum a q %+.5f, sum L %+.5f, sum q^2 %.5f, sum |q| %.4f e, score %+.5f"
          % (D["name"], n, D["saq"], D["sumL"], D["sq2"], D["sabs"], D["score"]))
    print("   %s parsed as the reference: sum |q| %.4f e, so the matched total fields are %.4f and %.4f e"
          % (REF["name"], REF["sabs"], REF_LAMS[0] * REF["sabs"], REF_LAMS[1] * REF["sabs"]))
    # (c) sites and a
    rows = [l.rstrip("\n").split("\t") for l in (T6 / "sites_f42436.tsv").read_text().splitlines()
            if not l.startswith("#")]
    hd, rows = rows[0], rows[1:]
    ci = {k: hd.index(k) for k in ("site", "x", "y", "z")}
    idx = {int(r[ci["site"]]): k for k, r in enumerate(rows)}
    A = np.loadtxt(str(T6 / "A_f42436.tsv"), comments="#").reshape(-1)
    for k, s in enumerate(D["site"]):
        if s not in idx:
            stop("site %d is not in sites_f42436.tsv" % s)
        r = rows[idx[s]]
        dd = float(np.abs(D["X"][k] - np.array([float(r[ci[c]]) for c in ("x", "y", "z")])).max())
        if dd > 1e-6:
            stop("site %d: ARRANGEMENTS.txt and sites_f42436.tsv differ by %.1e A" % (s, dd))
        if abs(A[idx[s]] - D["a"][k]) > 5e-5:
            stop("site %d: a is %+.5f in ARRANGEMENTS.txt and %+.5f in A_f42436.tsv"
                 % (s, D["a"][k], A[idx[s]]))
    print("   all %d sites and their a values agree with sites_f42436.tsv and A_f42436.tsv" % n)
    # (d) the per-site magnitude. NOT the withdrawn 0.80 e cap: E3.1's 0.65 e CLEAN observation.
    if D["maxq"] > QCLEAN + 1e-9:
        stop("%s's largest charge is %.4f e, above the %.2f e that E3.1 records as the largest single magnitude "
             "CLEAN at def2-SVPD. D7.2 does not scale such a design; the 0.80 e cap is WITHDRAWN and is not "
             "used here" % (D["name"], D["maxq"], QCLEAN))
    LAM = lambdas(D, REF)
    for nm, lam in LAM.items():
        if float(np.abs(lam * D["q"]).max()) > D["maxq"] + 1e-9:
            stop("at %s (lambda %.6f) the largest charge exceeds the design's own lambda = 1 maximum" % (nm, lam))
    print("   largest charge %.4f e at lambda 1.00, inside E3.1's %.2f e CLEAN observation; every lambda is below it"
          % (D["maxq"], QCLEAN))
    # (e) the headers
    if head("def2-SVP") != ("! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF\n%maxcore 3000\n"
                            "%pal nprocs 8 end\n%scf MaxIter 300 end\n"):
        stop("the def2-SVP header is not cs5_amatrix.py's HEAD byte for byte")
    ref = T6 / "qscan" / "C5" / "qp1000_R" / "job.inp"
    if head("def2-SVPD") != ref.read_text().split('%pointcharges')[0]:
        stop("the def2-SVPD header is not the committed qscan C5 header byte for byte")
    print("   both headers match their committed counterparts byte for byte")
    # the D6 anchors must exist now, not at analyse time
    C = controls(T6)
    for nb in (264, 402):
        if nb not in C:
            stop("D6's control table has no nb %d row; reading 4 could not be anchored" % nb)
    print("   D6's control medians are present for nb 264 and nb 402 (reading 4 can be anchored)")

    out.mkdir(parents=True)
    runs = []
    for basis, cond in BASES:
        for nm, lam in LAM.items():
            for st in ("R", "TS"):
                rn = "%s/%s_%s" % (basis, nm, st)
                d = out / rn
                d.mkdir(parents=True)
                (d / "substrate.xyz").write_text(geo[st].read_text())
                lines = ["%d" % n]
                for k in range(n):
                    lines.append("%10.6f %18.10f %18.10f %18.10f"
                                 % (lam * D["q"][k], D["X"][k][0], D["X"][k][1], D["X"][k][2]))
                (d / "pc.pc").write_text("\n".join(lines) + "\n")
                (d / "job.inp").write_text(head(basis) + '%pointcharges "pc.pc"\n'
                                           '* xyzfile -2 1 substrate.xyz\n')
                runs.append((rn, basis, "%.6f" % lam, st,
                             "%.4f" % (lam * D["sabs"]), "%.5f" % (lam * lam * D["sq2"]),
                             "%.4f" % (lam * D["maxq"])))
    with open(out / "TRANSFER_PLAN.tsv", "w") as fh:
        fh.write("# STAGE6_DIAGNOSTIC7_20261026.txt; cs6_transfer.py prepare\n")
        fh.write("# %s from relax/ARRANGEMENTS.txt: sum a q %+.5f, sum q^2 %.5f, sum L %+.5f, score %+.5f\n"
                 % (D["name"], D["saq"], D["sq2"], D["sumL"], D["score"]))
        fh.write("# matched lambdas target %s's sum |q| at lambda %s: %.4f and %.4f e\n"
                 % (REF["name"], " and ".join("%.2f" % x for x in REF_LAMS),
                    REF_LAMS[0] * REF["sabs"], REF_LAMS[1] * REF["sabs"]))
        fh.write("# bare references come from the committed qscan outputs archive, not re-run\n")
        fh.write("run\tbasis\tlambda\tstate\tsum_abs_q\tsum_q2\tmax_abs_q\n")
        for r in runs:
            fh.write("\t".join(r) + "\n")
    pbs = ["#!/bin/bash", "#PBS -N cm_s6t", "#PBS -l select=1:ncpus=8:mem=24gb", "#PBS -l walltime=12:00:00",
           "#PBS -m ae", "#PBS -M 18660916@sun.ac.za", "#PBS -j oe", "#PBS -o %s/transfer.pbs.out" % hpc,
           "# Stage 6 DIAGNOSTIC 7 (cs6_transfer.py): %d single points. Restartable." % len(runs),
           'export PATH="/apps/openmpi/4.1.1/bin:$PATH"',
           'export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:'
           '/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"',
           "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1", "ORCA=/home/apps2/ORCA/6.0.1", "W=%s" % hpc,
           'cd "$W" || exit 1',
           'if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running"; exit 1; fi',
           "trap 'rmdir \"$W/.running\" 2>/dev/null' EXIT", 'echo "host=$(hostname) start=$(date)"',
           "run() {", '  d="$W/$1"',
           '  if [ -s "$d/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$d/job.out"; then echo "$1 done before"',
           '  else ( cd "$d" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); '
           'echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi',
           "}"]
    pbs += ["run %s" % r[0] for r in runs] + ['echo "end=$(date)"']
    (out / "transfer.pbs").write_text("\n".join(pbs) + "\n")
    files = sorted(p for p in out.rglob("*") if p.is_file())
    (out / "INPUT_CHECKSUMS.sha256").write_text(
        "".join("%s  %s\n" % (sha(p), p.relative_to(out)) for p in files))
    print("   %d runs over %d lambda, %d files including the manifest, one restartable PBS job, nothing submitted"
          % (len(runs), len(LAM), len(files) + 1))
    for nm, lam in LAM.items():
        print("      %-7s lambda %.6f   sum |q| %.4f e   sum q^2 %.5f   largest |q| %.4f e"
              % (nm, lam, lam * D["sabs"], lam * lam * D["sq2"], lam * D["maxq"]))


def one(d, ref_xyz, npc):
    p = Path(d) / "job.out"
    if not p.is_file():
        stop("%s: no job.out - has the batch run?" % d)
    t = p.read_text(errors="replace")
    if "****ORCA TERMINATED NORMALLY****" not in t:
        stop("%s: ORCA did not terminate normally" % d)
    if not re.search(r"SCF CONVERGED AFTER\s+\d+\s+CYCLES", t) or "SCF NOT CONVERGED" in t:
        stop("%s: the SCF did not converge" % d)
    m = re.search(r"Number of Electrons\s+NEL\s+\.+\s+(\d+)", t)
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
    nb = re.search(r"Number of basis functions\s+\.+\s+(\d+)", t)
    if not nb:
        stop("%s: no basis-function count" % d)
    dd = AUD.read(t)        # D6's own parser, imported
    if dd is None:
        stop("%s: cs6_audit.read() does not accept this output as a completed single point on this dianion" % d)
    return float(e[0]), int(nb.group(1)), dd


def frozen_ddE(T6, name):
    """The committed frozen ddE and sum L for this design, at full precision, from frozen/frozen.tsv."""
    p = T6 / "frozen" / "frozen.tsv"
    if not p.is_file():
        stop("%s is missing: the consistency stop of D7.3 cannot be applied" % p)
    L = [l.rstrip("\n").split("\t") for l in p.read_text().splitlines() if l.strip()]
    hd = L[0]
    for f in L[1:]:
        if f[0].strip() == name:
            return float(f[hd.index("ddE_frozen_kcal")]), float(f[hd.index("sum_L_kcal")])
    stop("frozen.tsv has no row for %s" % name)


def ref_deltas(T6):
    """S6-20's measured Delta at each of its lambda, read out of the committed DIAGNOSTIC4_RESULT.txt table."""
    p = T6 / "field" / "DIAGNOSTIC4_RESULT.txt"
    if not p.is_file():
        stop("%s is missing: reading 3 has nothing to compare against" % p)
    out = {}
    for l in p.read_text().splitlines():
        m = re.match(r"\s*(\d\.\d{2})\s+(\d+\.\d+)\s+\d+\.\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+\d+\.\d+%"
                     r"\s+[+-]?\d+\.\d+\s+\d+\.\d+%\s+([+-]\d+\.\d+)\s*$", l)
        if m:
            out[round(float(m.group(1)), 2)] = (float(m.group(2)), float(m.group(3)))
    for rl in REF_LAMS:
        if round(rl, 2) not in out:
            stop("could not read %s's Delta at lambda %.2f out of DIAGNOSTIC4_RESULT.txt" % (REF_DESIGN, rl))
    return out


def analyse(root, name, n, runs, out):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    R = Path(runs)
    D = design(T6, name, n)
    REF = design(T6, REF_DESIGN, REF_N)
    LAM = lambdas(D, REF)
    C = controls(T6)
    # (i) b
    amr = T6 / "qrefine" / "DIAGNOSTIC2_AMENDMENT1_RESULT.txt"
    if not amr.is_file():
        stop("%s is missing" % amr)
    mb = re.search(r"b\s*=\s*([+-]?\d+\.\d+)\s*kcal/mol per e\^2", amr.read_text())
    if not mb:
        stop("could not read b out of the committed DIAGNOSTIC2_AMENDMENT1_RESULT.txt")
    Bq = float(mb.group(1))
    print("   b read from the committed DIAGNOSTIC2_AMENDMENT1_RESULT.txt: %+.4f kcal/mol per e^2" % Bq)
    print("   %s from ARRANGEMENTS.txt: sum a q %+.5f, sum q^2 %.5f, sum |q| %.4f e"
          % (D["name"], D["saq"], D["sq2"], D["sabs"]))
    print("   ONE b IS USED FOR ALL %d SITES AND THAT IS AN APPROXIMATION: b was measured at site 1690 only,"
          " at 2.96 A" % n)
    # (h) bare references
    man = {}
    for l in (T6 / "qscan" / "outputs" / "qscan_outputs.MANIFEST.sha256").read_text().splitlines():
        if l.strip():
            h, nm = l.split(None, 1)
            man[nm.strip()] = h
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
    if ref is None:
        stop("A_REPORT.txt has no frame 42436 row to check the bare energies against")
    for st, v in (("R", ref[0]), ("TS", ref[1])):
        if abs(Q[("def2-SVP", st)] - v) > 1e-6:
            stop("the committed qscan C1 bare %s energy does not match A_REPORT.txt" % st)
    print("   the committed archive's def2-SVP bare energies match A_REPORT.txt to <1e-6 Eh")
    # (g) the runs
    GEO = {st: read_xyz(T6 / "frame_42436" / ("reactant_qm.xyz" if st == "R" else "transition_state_qm.xyz"))
           for st in ("R", "TS")}
    E, NB, AU = {}, {}, {}
    for basis, cond in BASES:
        for nm, lam in LAM.items():
            for st in ("R", "TS"):
                d = R / basis / ("%s_%s" % (nm, st))
                E[(basis, nm, st)], NB[(basis, nm, st)], AU[(basis, nm, st)] = one(d, GEO[st], n)
        s = {NB[(basis, nm, st)] for nm in LAM for st in ("R", "TS")}
        if len(s) != 1:
            stop("%s runs report different basis-function counts: %s" % (basis, sorted(s)))
        print("   %-10s all %d runs at %d basis functions" % (basis, 2 * len(LAM), s.pop()))
    # ddE, the prediction, S
    ddE, pred, S = {}, {}, {}
    for basis, cond in BASES:
        for nm, lam in LAM.items():
            ddE[(basis, nm)] = ((E[(basis, nm, "TS")] - Q[(basis, "TS")])
                                - (E[(basis, nm, "R")] - Q[(basis, "R")])) * H
            pred[nm] = lam * D["saq"] + lam * lam * Bq * D["sq2"]
            S[(basis, nm)] = 100 * abs(ddE[(basis, nm)] - pred[nm]) / abs(pred[nm])
    # (k) THE CONSISTENCY STOP
    fz, fzL = frozen_ddE(T6, D["name"])
    if abs(fzL - D["sumL"]) > 1e-4:
        stop("frozen.tsv's sum L for %s is %+.6f, ARRANGEMENTS.txt gives %+.5f" % (D["name"], fzL, D["sumL"]))
    want = fz - fzL
    got = ddE[("def2-SVP", "lam100")]
    if abs(got - want) > STOP_TOL:
        stop("the def2-SVP lambda = 1 run gives ddE %+.4f; the committed frozen %s value minus its sum L is "
             "%+.6f, and Finding 4 says those must agree. %.4f kcal/mol apart - the scan does not reproduce the "
             "committed measurement and that must be found before anything is read"
             % (got, D["name"], want, abs(got - want)))
    print("   the def2-SVP lambda = 1 run reproduces the committed frozen %s value minus sum L: %+.4f against "
          "%+.6f" % (D["name"], got, want))
    # (j) S6-20's Delta at the matched fields
    RD = ref_deltas(T6)
    # reading 4 first: the verdicts, because no energy may be quoted without them
    V = {k: verdict(AU[k], C) for k in AU}
    # reading 1
    dom_target = REF_LAMS[1] * REF["sabs"]
    inside = [nm for nm, lam in LAM.items() if lam * D["sabs"] <= dom_target + 1e-9]
    ok_i = all(abs(ddE[("def2-SVPD", nm)] - ddE[("def2-SVP", nm)]) <= AGREE_KCAL for nm in inside)
    ok_ii = all(V[(b, nm, st)][0] == "CLEAN" for b, _ in BASES for nm in inside for st in ("R", "TS"))
    r1 = "TRANSFERS" if (ok_i and ok_ii) else ("DOES NOT TRANSFER" if not (ok_i or ok_ii) else "PARTIAL")
    # reading 2
    names = list(LAM)
    flags = [S[("def2-SVPD", nm)] <= PCT for nm in names]
    crossings = sum(1 for k in range(1, len(flags)) if flags[k] != flags[k - 1])
    if crossings > 1:
        r2, r2int = "NOT MONOTONE", None
    else:
        idx = [k for k, f in enumerate(flags) if f]
        if not idx:
            r2, r2int = "TIGHTER", (0.0, LAM[names[0]] * D["sabs"])
        else:
            kq = max(idx)
            lo = LAM[names[kq]] * D["sabs"]
            hi = LAM[names[kq + 1]] * D["sabs"] if kq + 1 < len(names) else float("inf")
            r2int = (lo, hi)
            rlo, rhi = REF_LAMS[0] * REF["sabs"], REF_LAMS[1] * REF["sabs"]
            r2 = "AGREES" if (lo < rhi and hi > rlo) else ("TIGHTER" if hi <= rlo else "LOOSER")
    # reading 3
    r3pairs = []
    for k, rl in enumerate(REF_LAMS, start=1):
        nm = "lamM%d" % k
        dD = ddE[("def2-SVPD", nm)] - ddE[("def2-SVP", nm)]
        dR = RD[round(rl, 2)][1]
        r3pairs.append((rl * REF["sabs"], dD, dR, abs(dD - dR)))
    r3 = "MATCHED" if all(p[3] <= DELTA_TOL for p in r3pairs) else "NOT MATCHED"

    O = ["STAGE 6 - DIAGNOSTIC 7 RESULT: does the defensible domain transfer? %s scaled" % D["name"], "",
         "Measured by cs6_transfer.py analyse. The rule is STAGE6_DIAGNOSTIC7_20261026.txt D7.3, fixed before",
         "the runs went in and applied unchanged. Nothing is refitted: sum a q, sum q^2 and sum L come from the",
         "committed ARRANGEMENTS.txt, b from the committed DIAGNOSTIC2_AMENDMENT1_RESULT.txt, %s's Delta from"
         % REF["name"],
         "the committed DIAGNOSTIC4_RESULT.txt, and D6's control medians from the committed",
         "DIAGNOSTIC6_RESULT.txt.", "",
         "   predicted(lambda) = lambda * sum a q  +  lambda^2 * b * sum q^2",
         "   sum a q = %+.5f kcal/mol   sum q^2 = %.5f e^2   b = %+.4f kcal/mol per e^2" % (D["saq"], D["sq2"], Bq),
         "   ONE b IS USED FOR ALL %d SITES AND THAT IS AN APPROXIMATION: b was measured at site 1690 only." % n,
         "   %s's closest sphere-to-substrate distance is recorded in ARRANGEMENTS.txt; no cap as a function of"
         % D["name"],
         "   distance is measured here (D7.4).", "",
         "   lambda      S|q|    largest|q|  predicted   def2-SVP      S     def2-SVPD      S       Delta"]
    for nm in names:
        lam = LAM[nm]
        O.append("   %-7s %.6f  %6.4f  %6.4f  %9.4f  %9.4f %6.1f%%  %9.4f %6.1f%%  %+8.4f"
                 % (nm, lam, lam * D["sabs"], lam * D["maxq"], pred[nm],
                    ddE[("def2-SVP", nm)], S[("def2-SVP", nm)],
                    ddE[("def2-SVPD", nm)], S[("def2-SVPD", nm)],
                    ddE[("def2-SVPD", nm)] - ddE[("def2-SVP", nm)]))
    O += ["",
          "   lamM1 and lamM2 are the matched points: %s carries sum |q| %.4f and %.4f e there, which are"
          % (D["name"], REF_LAMS[0] * REF["sabs"], REF_LAMS[1] * REF["sabs"]),
          "   %s's own total fields at its lambda %s." % (REF["name"], " and ".join("%.2f" % x for x in REF_LAMS)),
          "",
          "READING 4 FIRST: THE D6 VERDICTS. No energy above may be quoted without its verdict in the same",
          "sentence (E3.4's rule, applied here in advance of the runs by D7.3).",
          "   D6's thresholds and parser are cs6_audit.py's, imported unchanged; the control medians are the",
          "   committed DIAGNOSTIC6_RESULT.txt ones, per basis-set size.", "",
          "   run                        nb   cyc   gap(eV)   dip(D)  max|Lw|   verdict and why"]
    for basis, cond in BASES:
        for nm in names:
            for st in ("R", "TS"):
                d, (v, w) = AU[(basis, nm, st)], V[(basis, nm, st)]
                O.append("   %-24s %4s %5s %9s %8s %8s   %s%s"
                         % ("%s/%s_%s" % (basis, nm, st), d["nb"], d["cyc"],
                            "%.4f" % d["gap"] if d["gap"] is not None else "-",
                            "%.4f" % d["dip"] if d["dip"] is not None else "-",
                            "%.3f" % d["maxlw"] if d["maxlw"] is not None else "-",
                            v, (": " + "; ".join(w)) if w else ""))
    O += ["", "   control medians used: " + ", ".join(
        "nb %d (%s) cyc %d gap %.3f maxLw %.3f dip %.2f"
        % (nb, C[nb]["basis"], C[nb]["cyc"], C[nb]["gap"], C[nb]["maxlw"], C[nb]["dip"])
        for nb in sorted({NB[k] for k in NB}) if nb in C)]
    O += ["", "THE FOUR READINGS OF D7.3", "",
          "   1. DOES THE DEFENSIBLE DOMAIN TRANSFER?  %s" % r1,
          "      evaluated at the lambda with sum |q| <= %.4f e: %s" % (dom_target, ", ".join(inside)),
          "      (i)  inter-basis agreement, |Delta| <= %.1f kcal/mol at each:  %s"
          % (AGREE_KCAL, "HOLDS" if ok_i else "FAILS"),
          "      (ii) every run CLEAN at BOTH bases by D6.3's rule:             %s"
          % ("HOLDS" if ok_ii else "FAILS")]
    if r1 == "PARTIAL":
        O += ["      PARTIAL: the domain rests on ONE criterion for %s where it rested on two for %s. E3.3"
              % (D["name"], REF["name"]),
              "      rested the claim on the coincidence of two unrelated criteria, and half of it is a weaker",
              "      claim, not the same one. That must be stated wherever this is quoted."]
    O += ["", "   2. DOES THE QUADRATIC-MODEL BOUND TRANSFER?  %s" % r2]
    if r2int:
        O += ["      lambda_quad(%s) puts the bound in [%.4f, %s) e; %s's own interval is [%.4f, %.4f) e."
              % (D["name"], r2int[0], "inf" if r2int[1] == float("inf") else "%.4f" % r2int[1],
                 REF["name"], REF_LAMS[0] * REF["sabs"], REF_LAMS[1] * REF["sabs"]),
              "      The bound is an INTERVAL and never a value: E3.3 requires that it never be quoted as",
              "      sharply determined."]
    else:
        O += ["      S crosses %g per cent %d times across the %d lambda, so NO bound is stated for %s at all."
              % (PCT, crossings, len(names), D["name"]),
              "      The crossings are in the table above. This reading returns no verdict (D7.3)."]
    O += ["", "   3. IS Delta A FUNCTION OF TOTAL FIELD ALONE?  %s" % r3,
          "      at matched total field, %s against %s's committed Delta:" % (D["name"], REF["name"]),
          "         sum |q| (e)    %s      %s     |difference|   tolerance %.1f" % (D["name"], REF["name"], DELTA_TOL)]
    for sq, dD, dR, df in r3pairs:
        O.append("         %8.4f    %+8.4f   %+8.4f     %8.4f      %s"
                 % (sq, dD, dR, df, "within" if df <= DELTA_TOL else "OUTSIDE"))
    O += ["      The %.1f kcal/mol is carried over from D4.3, where it bounds the INTERCEPT of a linearity fit"
          % DELTA_TOL,
          "      at fixed design. Here it bounds the DIFFERENCE between two designs at matched total field,",
          "      which is a different quantity. D7.3 discloses that it is carried over by choice.",
          "", "   4. THE D6 TABLE  - above, before any energy was read.", "",
          "WHAT THIS DOES NOT TEST", "",
          "   Two designs, so this is about %s and %s and not about arrangements in general. One frame, one"
          % (D["name"], REF["name"]),
          "   geometry pair, fixed geometry, so nothing about Loss or a relaxed stationary point. Whether b",
          "   varies between sites. NO CAP AS A FUNCTION OF DISTANCE, so Stage 7's requirement that a cap be a",
          "   magnitude AND a distance is not met here. Negative total field. Nothing about the 1.0000 e",
          "   designs, ma-def2-SVP, or a charge model with a Pauli wall for electrons."]
    Path(out).write_text("\n".join(O) + "\n")
    print("   wrote %s" % out)
    print("   reading 1: %s    reading 2: %s    reading 3: %s" % (r1, r2, r3))


if __name__ == "__main__":
    A = sys.argv[1:]
    if len(A) == 6 and A[0] == "prepare":
        prepare(A[1], A[2], int(A[3]), A[4], A[5])
    elif len(A) == 6 and A[0] == "analyse":
        analyse(A[1], A[2], int(A[3]), A[4], A[5])
    else:
        sys.exit(__doc__.strip().splitlines()[2].strip() + "\n"
                 + "  python3 cs6_transfer.py prepare <results_repo> <DESIGN> <N> <OUT> <HPC_DIR>\n"
                 + "  python3 cs6_transfer.py analyse <results_repo> <DESIGN> <N> <RUNS> <OUT>")
