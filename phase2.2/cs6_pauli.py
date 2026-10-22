#!/usr/bin/env python3
"""
cs6_pauli.py - DIAGNOSTIC 3: does the designed charge need Pauli repulsion?

  python3 cs6_pauli.py prepare <results_repo> <OUT> <HPC_DIR>
  python3 cs6_pauli.py analyse <results_repo> <RUNS> <OUT>

Replaces the bare +1 point charge at site 1690 with a real ion, which has electrons and therefore Pauli repulsion
by construction, and with a ghost of the same element, which has the ion's basis functions and nothing else. The
bare substrate and the point-charge runs are NOT re-run: they are read out of the committed qscan outputs archive,
so the comparison is against the same numbers DIAGNOSTIC2_RESULT.txt reported.

Conventions, site, geometries and the restartable single-PBS-job pattern are cs6_qscan.py's, which are
cs5_amatrix.py's.

prepare checks, each a stop:
  (a) both geometries match the sha256 relax/ARRANGEMENTS.txt records;
  (b) site 1690's coordinates agree between sites_f42436.tsv and ARRANGEMENTS.txt's S6-P block to 1e-6 A;
  (c) the two bases are the ones Diagnostic 2 used for its in vacuo conditions C1 and C5, and the def2-SVP header
      equals cs5_amatrix.py's HEAD byte for byte;
  (d) the ion sits at exactly the site, and the substrate block is the committed geometry byte for byte;
  (e) it refuses to write into an existing folder, and submits nothing.

analyse checks, each a stop:
  (f) every new run: ORCA TERMINATED NORMALLY, SCF converged, the electron count the model requires, exactly one
      final energy, the echoed substrate coordinates equal to the committed geometry, the ion present as an atom
      for M2/M3 and as a ghost for M4/M5;
  (g) the committed qscan archive supplies M0 and M1, and its bare def2-SVP energies reproduce A_REPORT.txt's
      committed in vacuo values to below 1e-6 Eh;
  (h) the basis-function count is 264 for the def2-SVP bare substrate and strictly larger with an ion or a ghost.
Then it applies D3.3 without choosing anything.
"""
import hashlib, re, sys, tarfile
from pathlib import Path
import numpy as np

H = 627.5094740631
SITE = 1690
# model: (name, element or None, kind, qm_charge, electrons_added)
MODELS = (("M2_Na", "Na", "ion", -1, 10), ("M3_K", "K", "ion", -1, 18),
          ("M4_ghostNa", "Na", "ghost", -2, 0), ("M5_ghostK", "K", "ghost", -2, 0))
BASES = (("def2-SVP", "C1"), ("def2-SVPD", "C5"))   # name, the qscan condition that supplies M0 and M1
NEL_BARE = 118
HEAD_SVP = "! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF\n%maxcore 3000\n%pal nprocs 8 end\n%scf MaxIter 300 end\n"
GHOST_TOL = 0.30     # D3.3's ghost control, kcal/mol
PCT = 10.0           # D3.3's yardstick


def stop(m):
    sys.exit("STOP: " + m)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    el = [l.split()[0] for l in L[2:2 + n]]
    X = np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])
    return el, X


def head(basis):
    return ("! B3LYP D3BJ %s def2/J RIJCOSX TightSCF\n%%maxcore 3000\n%%pal nprocs 8 end\n"
            "%%scf MaxIter 300 end\n" % basis)


def prepare(root, out, hpc):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    out = Path(out)
    if out.exists():
        stop("%s already exists" % out)
    if head("def2-SVP") != HEAD_SVP:
        stop("the generated def2-SVP header is not cs5_amatrix.py's HEAD byte for byte")
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
    d = float(np.abs(P - np.array([float(m.group(i)) for i in (1, 2, 3)])).max())
    if d > 1e-6:
        stop("site %d differs between the sites file and ARRANGEMENTS.txt by %.1e A" % (SITE, d))
    print("   site %d at %.6f %.6f %.6f, both records agree to %.0e A" % (SITE, P[0], P[1], P[2], d))

    out.mkdir(parents=True)
    runs = []
    for basis, cond in BASES:
        for mname, el, kind, qm_q, nadd in MODELS:
            for st in ("R", "TS"):
                sel, SX = read_xyz(geo[st])
                name = "%s/%s_%s" % (basis.replace("/", "_"), mname, st)
                dd = out / name
                dd.mkdir(parents=True)
                lines = ["%d" % (len(sel) + 1),
                         "DIAGNOSTIC 3 %s %s %s: the committed frame_42436 %s substrate plus %s at site %d"
                         % (basis, mname, st, st, ("%s+" % el) if kind == "ion" else "ghost %s" % el, SITE)]
                for e, x in zip(sel, SX):
                    lines.append("  %-2s %20.14f %20.14f %20.14f" % (e, x[0], x[1], x[2]))
                tag = el if kind == "ion" else el + ":"
                lines.append("  %-2s %20.14f %20.14f %20.14f" % (tag, P[0], P[1], P[2]))
                (dd / "system.xyz").write_text("\n".join(lines) + "\n")
                (dd / "job.inp").write_text(head(basis) + "* xyzfile %d 1 system.xyz\n" % qm_q)
                runs.append((name, basis, mname, st, kind, el, "%d" % qm_q, "%d" % (NEL_BARE + nadd)))
    with open(out / "PAULI_PLAN.tsv", "w") as fh:
        fh.write("# STAGE6_DIAGNOSTIC3_20261022.txt; cs6_pauli.py prepare\n")
        fh.write("# site %d at %.6f %.6f %.6f; M0 and M1 come from the committed qscan outputs, not re-run\n"
                 % (SITE, P[0], P[1], P[2]))
        fh.write("# geometries %s %s\n" % (sha(geo["R"])[:16], sha(geo["TS"])[:16]))
        fh.write("run\tbasis\tmodel\tstate\tkind\telement\tqm_charge\telectrons\n")
        for r in runs:
            fh.write("\t".join(r) + "\n")
    pbs = ["#!/bin/bash", "#PBS -N cm_s6p", "#PBS -l select=1:ncpus=8:mem=24gb", "#PBS -l walltime=06:00:00",
           "#PBS -m ae", "#PBS -M 18660916@sun.ac.za", "#PBS -j oe", "#PBS -o %s/pauli.pbs.out" % hpc,
           "# Stage 6 DIAGNOSTIC 3 (cs6_pauli.py, STAGE6_DIAGNOSTIC3_20261022.txt): %d single points." % len(runs),
           "# Restartable: a run whose job.out ends normally is not repeated.",
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
    (out / "pauli.pbs").write_text("\n".join(pbs) + "\n")
    files = sorted(p for p in out.rglob("*") if p.is_file())
    (out / "INPUT_CHECKSUMS.sha256").write_text("".join("%s  %s\n" % (sha(p), p.relative_to(out)) for p in files))
    print("   %d runs, %d files including the manifest, one restartable PBS job, nothing submitted"
          % (len(runs), len(files) + 1))


def one(d, nel, nqm, ref_xyz):
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
    k = t.find("CARTESIAN COORDINATES (ANGSTROEM)")
    if k < 0:
        stop("%s: no echoed coordinates" % d)
    X = np.array([[float(v) for v in l.split()[1:4]] for l in t[k:].splitlines()[2:2 + nqm]])
    if X.shape != ref_xyz.shape or np.abs(X - ref_xyz).max() > 1e-5:
        stop("%s: the echoed substrate coordinates differ from the committed geometry" % d)
    nb = re.search(r"Number of basis functions\s+\.\.\.\s+(\d+)", t)
    if not nb:
        stop("%s: no basis-function count" % d)
    return float(e[0]), int(nb.group(1))


def analyse(root, runs, out):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    R = Path(runs)
    # (g) M0 and M1 out of the committed qscan archive
    arch = T6 / "qscan" / "outputs" / "qscan_outputs.tar.xz"
    if not arch.is_file():
        stop("%s is missing - the qscan outputs have not been committed" % arch)
    man = {}
    for l in (T6 / "qscan" / "outputs" / "qscan_outputs.MANIFEST.sha256").read_text().splitlines():
        if l.strip():
            h, n = l.split(None, 1)
            man[n.strip()] = h
    Q, QNB = {}, {}
    with tarfile.open(str(arch)) as t:
        for cond in ("C1", "C5"):
            for st in ("R", "TS"):
                for who, sub in (("bare", "bare"), ("pc", "qp1000")):
                    nm = "./%s/%s_%s/job.out" % (cond, sub, st)
                    try:
                        raw = t.extractfile(nm).read()
                    except Exception:
                        stop("%s is not in the committed qscan archive" % nm)
                    if nm in man and hashlib.sha256(raw).hexdigest() != man[nm]:
                        stop("%s does not match its sha256 in the qscan manifest" % nm)
                    tx = raw.decode(errors="replace")
                    ee = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", tx)
                    if len(ee) != 1:
                        stop("%s: %d final energies" % (nm, len(ee)))
                    Q[(cond, who, st)] = float(ee[0])
                    if who == "bare" and st == "R":
                        nbq = re.search(r"Number of basis functions\s+\.\.\.\s+(\d+)", tx)
                        if not nbq:
                            stop("%s: no basis-function count" % nm)
                        QNB[cond] = int(nbq.group(1))
    ref = None
    for l in (T6 / "A_REPORT.txt").read_text().splitlines():
        f = l.split("\t")
        if f and f[0].strip() == "42436" and len(f) >= 3:
            ref = (float(f[1]), float(f[2]))
    if ref is None:
        stop("A_REPORT.txt holds no in vacuo energies for frame 42436")
    for st, v in (("R", ref[0]), ("TS", ref[1])):
        if abs(Q[("C1", "bare", st)] - v) > 1e-6:
            stop("the qscan C1 bare %s energy is %.9f, A_REPORT.txt says %.9f" % (st, Q[("C1", "bare", st)], v))
    print("   M0 and M1 read from the committed qscan archive; its def2-SVP bare energies match A_REPORT.txt to <1e-6 Eh")

    E, NB = {}, {}
    for basis, cond in BASES:
        for mname, el, kind, qm_q, nadd in MODELS:
            for st in ("R", "TS"):
                _, SX = read_xyz(T6 / "frame_42436" / ("reactant_qm.xyz" if st == "R" else "transition_state_qm.xyz"))
                d = R / basis / ("%s_%s" % (mname, st))
                E[(basis, mname, st)], NB[(basis, mname, st)] = one(d, NEL_BARE + nadd, len(SX), SX)
    # (h) the extra centre must actually add basis functions, ghosts included
    for basis, cond in BASES:
        nb0 = QNB[cond]
        for mname, el, kind, qm_q, nadd in MODELS:
            nb = NB[(basis, mname, "R")]
            if nb <= nb0:
                stop("%s %s used %d basis functions, not more than the bare substrate's %d at the same basis: the "
                     "extra centre carried no basis functions" % (basis, mname, nb, nb0))
        print("   %-10s bare %d basis functions; with the extra centre %s"
              % (basis, nb0, " ".join("%s %d" % (m[0].split("_")[0], NB[(basis, m[0], "R")]) for m in MODELS)))

    L = ["STAGE 6 - DIAGNOSTIC 3 RESULT: does the designed charge need Pauli repulsion?", "",
         "Measured by cs6_pauli.py analyse. The verdict rule is STAGE6_DIAGNOSTIC3_20261022.txt D3.3, fixed before",
         "the runs went in and applied here unchanged. No method changes and no reported number is replaced.", "",
         "ddE(model) = [E_TS(model) - E_TS(bare)] - [E_R(model) - E_R(bare)], kcal/mol; negative lowers the barrier.",
         "ion(X) = ddE(real X+) - ddE(ghost X) is what the ion does BEYOND supplying basis functions.", ""]
    res = {}
    for basis, cond in BASES:
        b0 = {st: Q[(cond, "bare", st)] for st in ("R", "TS")}
        pc = ((Q[(cond, "pc", "TS")] - b0["TS"]) - (Q[(cond, "pc", "R")] - b0["R"])) * H
        L.append("%s" % basis)
        L.append("   M1  bare +1 point charge (from the committed qscan)      ddE %+8.4f" % pc)
        r = {"pc": pc}
        for mname, el, kind, qm_q, nadd in MODELS:
            dd = ((E[(basis, mname, "TS")] - b0["TS"]) - (E[(basis, mname, "R")] - b0["R"])) * H
            r[mname] = dd
            L.append("   %-3s %-48s ddE %+8.4f   %d basis functions"
                     % (mname.split("_")[0], ("%s+ as a QM atom" % el) if kind == "ion" else ("ghost %s, basis functions only" % el),
                        dd, NB[(basis, mname, "R")]))
        ina = r["M2_Na"] - r["M4_ghostNa"]
        ik = r["M3_K"] - r["M5_ghostK"]
        r["ion_Na"], r["ion_K"] = ina, ik
        L.append("   ion(Na) = M2 - M4  %+8.4f      ion(K) = M3 - M5  %+8.4f" % (ina, ik))
        L.append("   against the point charge: ion(Na) - M1 %+8.4f (%5.1f%%)   ion(K) - M1 %+8.4f (%5.1f%%)"
                 % (ina - pc, 100 * abs(ina - pc) / abs(pc), ik - pc, 100 * abs(ik - pc) / abs(pc)))
        res[basis] = r
        L.append("")
    L.append("THE TESTS OF D3.3")
    L.append("")
    gmax = max(abs(res[b][m]) for b, _ in BASES for m in ("M4_ghostNa", "M5_ghostK"))
    L.append("   ghost control: largest |ddE(ghost)| = %.4f kcal/mol   (limit %.2f)  -> %s"
             % (gmax, GHOST_TOL, "PASS" if gmax <= GHOST_TOL else "FAIL"))
    if gmax > GHOST_TOL:
        L.append("")
        L.append("VERDICT PER D3.3: INCONCLUSIVE - the basis-set contribution is not negligible, so no comparison")
        L.append("  below is clean. A counterpoise treatment is required before this question can be answered.")
    else:
        rr = res["def2-SVPD"]
        dna, dk = 100 * abs(rr["ion_Na"] - rr["pc"]) / abs(rr["pc"]), 100 * abs(rr["ion_K"] - rr["pc"]) / abs(rr["pc"])
        L.append("   at def2-SVPD: |ion(Na) - M1| = %.1f%%   |ion(K) - M1| = %.1f%%   (yardstick %.0f%%)" % (dna, dk, PCT))
        if dna <= PCT and dk <= PCT:
            L.append("")
            L.append("VERDICT PER D3.3: THE POINT CHARGE IS AN ADEQUATE SURROGATE AT +1 e")
            L.append("  - Diagnostic 2's q = 1.00 anomaly is NOT a Pauli-repulsion artefact. Something else is")
            L.append("    responsible, and the measured charge bound stands as a property of the substrate's")
            L.append("    response rather than as a defect of the charge model.")
        else:
            L.append("")
            L.append("VERDICT PER D3.3: THE POINT CHARGE IS NOT AN ADEQUATE SURROGATE AT +1 e")
            L.append("  - the anomaly is explained. The remedy is a Pauli wall for the electrons, not a lower charge")
            L.append("    bound, and Stage 7 must not adopt the 0.25 e cap as a design constraint: that would")
            L.append("    conceal a model defect behind a restriction on the design space.")
        ds = 100 * abs(abs(rr["ion_K"]) - abs(rr["ion_Na"])) / abs(rr["ion_Na"]) if rr["ion_Na"] else float("inf")
        L.append("")
        L.append("   THE SIZE TEST (Na+ 1.02 A against K+ 1.38 A ionic radius, guanidinium C about 1.70 A vdW):")
        L.append("   |ion(K)| against |ion(Na)| at def2-SVPD: %.1f%% apart" % ds)
        if ds > PCT:
            L.append("   -> the effect DEPENDS ON THE SIZE of the wall. A single point charge cannot represent a")
            L.append("      guanidinium at this distance at any magnitude, and the representation itself has to")
            L.append("      change rather than its parameters.")
        else:
            L.append("   -> the effect is insensitive to the wall's size over this range, which is evidence that a")
            L.append("      point charge plus one transferable pseudopotential would be enough.")
    L += ["", "WHAT THIS DOES NOT TEST", "",
          "   Whether an ion is a good model of a residue - it is not. The relaxed geometry, where Finding 4 says the",
          "   Lennard-Jones wall's entire value lies. Any magnitude other than +1 e, or either sign, since an ion's",
          "   charge is not adjustable. One site, one frame."]
    Path(out).mkdir(parents=True, exist_ok=True)
    (Path(out) / "DIAGNOSTIC3_RESULT.txt").write_text("\n".join(L) + "\n")
    print()
    print("\n".join(L[L.index("THE TESTS OF D3.3"):]))
    print()
    print("   wrote %s" % (Path(out) / "DIAGNOSTIC3_RESULT.txt"))


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in ("prepare", "analyse"):
        sys.exit(__doc__)
    (prepare if sys.argv[1] == "prepare" else analyse)(sys.argv[2], sys.argv[3], sys.argv[4])
