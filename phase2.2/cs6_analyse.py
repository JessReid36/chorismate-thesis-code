#!/usr/bin/env python3
"""
cs6_analyse.py - Stage 6's measurements (phase2.2, cs_stage6/STAGE6_CRITERIA.txt with STAGE6_AMENDMENT1.txt,
STAGE6_ERRATUM1.txt, STAGE6_ERRATUM2.txt and relax/ADDITIVITY_RUNS.txt).

  python3 cs6_analyse.py selftest
  python3 cs6_analyse.py frozen <relax_dir> <results_repo> <OUT>

FROZEN is the half that needs only the single points: the design's effect at the fixed geometries and the added
gradient it puts on the substrate. The relaxed half - Loss, the two survival criteria, G0 and kappa - waits for the
relaxations and is a separate mode, because nothing about it can be computed until they have run.

WHAT IT COMPUTES, all of it from STAGE6_CRITERIA.txt
  ddE(frozen)   [E(TS,D) - E(R,D)] - [E(TS,0) - E(R,0)], with 0 the no-sphere control S6-C. Negative means the
                design lowers the barrier at the fixed geometries. This is the quantity the MILP's score estimates.
  G_add(X)      sqrt( (1/3N) sum_i |g_i(D) - g_i(0)|^2 ) over the N = 24 substrate atoms, in kcal/mol/A - an RMS
                over all 3N Cartesian components, as the criteria defines it, of the ADDED gradient.
  the score     Sum_i a_i q_i + Sum_i L_i for the same design, read from relax/ARRANGEMENTS.txt, so the measured
                effect can be set against the linear estimate the optimiser worked from. The difference is the
                polarisation the score leaves out (decision O1) plus whatever the first-order model misses.
  additivity    for S6-2 and S6-3: the sum of their one-sphere measured effects against their measured effect, from
                the ADD-* runs. The criteria asks for this at fixed geometry, single points only.

CHECKS, every one a stop
  (a) every run: ORCA TERMINATED NORMALLY, an SCF that converged, 118 electrons in the QM region, exactly one final
      energy, and the echoed coordinates equal to combined.xyz within 1e-5 A;
  (b) the number of MM point charges ORCA reports equals the number of spheres the arrangement has - so a design
      whose force field was misread cannot pass unnoticed. The control must report no QM/MM section at all;
  (c) every job.engrad holds exactly the 24 QM atoms PLUS that arrangement's spheres - an EnGrad job writes the
      whole system, unlike the Opt job WP4 left behind, which writes only the atoms it optimised - and which rows
      are the substrate is read from the engrad's own atom block rather than assumed: its atomic numbers must be
      the substrate's in order, the rows after them must be nitrogen, its coordinates must equal combined.xyz, and
      its energy must equal job.out's QM/MM total;
  (d) every input file is unchanged since it was generated (relax/INPUT_CHECKSUMS.sha256);
  (e) the control's two energies equal the in vacuo energies the A-matrix job already computed for this frame at
      these geometries (cs_stage6/A_REPORT.txt), within 1e-6 Eh. That is an independent check: the A column's
      single points and the control are the same method, basis and geometry by two different routes, so if they
      disagree either the control or the A column is wrong. It is also what settles that the control being a plain
      QM job rather than QM/MM with an empty MM region costs nothing.

Units: ORCA writes the gradient in Eh/bohr; 1 Eh/bohr = 627.5094740631 / 0.5291772083 kcal/mol/A, the project's own
two constants (cs5_amatrix.py).

Python 3.8+, numpy.
"""
import hashlib
import lzma
import re
import sys
from pathlib import Path

import numpy as np

H = 627.5094740631                      # kcal/mol per Eh
BOHR = 0.5291772083                     # A per bohr
GCONV = H / BOHR                        # kcal/mol/A per Eh/bohr
NQM = 24
NEL = 118
NAMES = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
ELS = [n.rstrip("0123456789") for n in NAMES]
ZS = [{"C": 6, "H": 1, "O": 8}[e] for e in ELS]
STATES = ("R", "TS")


def stop(m):
    sys.exit("STOP: " + m)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    return np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])


def read_engrad(p):
    """(natoms, energy in Eh, gradient (N,3) in Eh/bohr, atomic numbers, coordinates in Bohr (N,3)).

    ORCA's .engrad ends with "# The atomic numbers and current coordinates in Bohr", one row per atom. That block
    is what makes the file self-describing: it says which row is which atom, so the substrate's rows can be
    identified from the file itself instead of assumed from the order combined.xyz happens to use."""
    v = [l.strip() for l in Path(p).read_text().splitlines() if l.strip() and not l.strip().startswith("#")]
    if len(v) < 2:
        stop("%s is too short to be an engrad" % p)
    n = int(v[0])
    e = float(v[1])
    if len(v) < 2 + 3 * n:
        stop("%s has fewer gradient values than %d atoms need" % (p, n))
    g = np.array([float(x) for x in v[2:2 + 3 * n]], dtype=float).reshape(n, 3)
    rest = [r.split() for r in v[2 + 3 * n:]]
    rows = [r for r in rest if len(r) == 4]
    if len(rows) != n:
        stop("%s: its atom block has %d rows, not the %d atoms it declares - this reader needs that block to tell"
             " which rows are the substrate" % (p, len(rows), n))
    Z = [int(float(r[0])) for r in rows]
    X = np.array([[float(x) for x in r[1:4]] for r in rows], dtype=float)
    return n, e, g, Z, X


def read_run(d, nsph):
    """the checks of (a), (b) and (c); returns (total energy Eh, MM part Eh, gradient (24,3) Eh/bohr)"""
    d = Path(d)
    eg = d / "job.engrad"
    # job.out is committed compressed, as WP4's outputs are, so a committed tree can be re-analysed without hpc1.
    if (d / "job.out").is_file():
        t = (d / "job.out").read_text(errors="replace")
    elif (d / "job.out.xz").is_file():
        t = lzma.open(str(d / "job.out.xz"), "rt", errors="replace").read()
    else:
        stop("%s: no job.out or job.out.xz - has the job run?" % d)
    if not eg.is_file():
        stop("%s: no job.engrad - did it ask for EnGrad?" % d)
    if "****ORCA TERMINATED NORMALLY****" not in t:
        stop("%s: ORCA did not terminate normally" % d)
    if not re.search(r"SCF CONVERGED AFTER\s+\d+\s+CYCLES", t) or "SCF NOT CONVERGED" in t:
        stop("%s: the SCF did not converge" % d)
    m = re.search(r"Number of Electrons\s+NEL\s+\.\.\.\.\s+(\d+)", t)
    if not m or int(m.group(1)) != NEL:
        stop("%s: %s electrons in the QM region, expected %d" % (d, m.group(1) if m else "no count", NEL))
    # A QM/MM job prints THREE energy lines per evaluation and a plain job one, so the three are read separately:
    #   FINAL SINGLE POINT ENERGY (MM)      the molecular-mechanics part - here the sphere's Lennard-Jones energy
    #   FINAL SINGLE POINT ENERGY           the QM energy, with the MM point charges already embedded in it
    #   FINAL SINGLE POINT ENERGY (QM/MM)   their sum, which is the total and what job.engrad holds
    # Taking the QM/MM total is forced, not a preference: the score this is measured against contains L, the
    # sphere's Lennard-Jones energy at TS minus at R, so the measurement must contain it too. Checked against WP4's
    # committed armb_freecart/job.out, where bare + (MM) equals (QM/MM) exactly and the engrad energy is the total.
    # A pattern with an OPTIONAL "(QM/MM)" suffix matches the bare line and the total both, which is why this reads
    # the three apart rather than allowing the suffix.
    qm = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    mm = re.findall(r"FINAL SINGLE POINT ENERGY \(MM\)\s+(-?\d+\.\d+)", t)
    tot = re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", t)
    if nsph == 0:
        if mm or tot:
            stop("%s: the control reports MM or QM/MM energies; it is meant to be a plain QM job" % d)
        if len(qm) != 1:
            stop("%s: %d final energies, expected 1" % (d, len(qm)))
        e_mm, e_tot = 0.0, float(qm[0])
    else:
        if not (len(qm) == len(mm) == len(tot) == 1):
            stop("%s: %d bare, %d (MM) and %d (QM/MM) final energies, expected one of each for a single point"
                 % (d, len(qm), len(mm), len(tot)))
        e_mm, e_tot = float(mm[0]), float(tot[0])
        if abs(float(qm[0]) + e_mm - e_tot) > 1e-9:
            stop("%s: the QM energy plus the MM energy is not the QM/MM total (%.12f + %.12f != %.12f)"
                 % (d, float(qm[0]), e_mm, e_tot))
    e = [e_tot]
    k = t.find("CARTESIAN COORDINATES (ANGSTROEM)")
    if k < 0:
        stop("%s: no echoed coordinates" % d)
    X = np.array([[float(v) for v in l.split()[1:4]] for l in t[k:].splitlines()[2:2 + NQM]])
    C = read_xyz(d / "combined.xyz")[:NQM]
    if X.shape != (NQM, 3) or np.abs(X - C).max() > 1e-5:
        stop("%s: the echoed coordinates differ from combined.xyz by %.1e A" % (d, np.abs(X - C).max()))
    # (b) the MM point charges ORCA actually read
    mm = re.search(r"Point charges in QM calc\. from MM atoms\.*\s*(\d+)", t)
    if nsph == 0:
        if "Multiscale model" in t or mm:
            stop("%s: the control ran as QM/MM; it is meant to be a plain QM job" % d)
    else:
        if not mm:
            stop("%s: ORCA reports no MM point charges, but this arrangement has %d sphere(s)" % (d, nsph))
        if int(mm.group(1)) != nsph:
            stop("%s: ORCA read %s MM point charges, the arrangement has %d sphere(s)" % (d, mm.group(1), nsph))
    n, ee, g, Z, Xb = read_engrad(eg)
    # An EnGrad job writes the gradient of EVERY atom in the system - the 24 QM atoms and the MM spheres - so the
    # count is 24 + the sphere count. WP4's committed engrad holds 24 because it is an Opt job, where ORCA writes
    # only the atoms it was optimising; that file verified the format but not this, and the difference is the job
    # type, which is why this is checked against the arrangement rather than against WP4.
    if n != NQM + nsph:
        stop("%s: job.engrad holds %d atoms, expected %d (the %d QM atoms and %d sphere(s))"
             % (d, n, NQM + nsph, NQM, nsph))
    # WHICH rows are the substrate is read from the engrad's own atom block, not assumed: its atomic numbers must be
    # the substrate's, in order, and its coordinates must be combined.xyz's to within rounding. That identifies the
    # rows and re-checks the geometry at the same time, at the engrad's own precision.
    if Z[:NQM] != ZS:
        stop("%s: the first %d rows of job.engrad are not the substrate's atoms.\n   got      %s\n   expected %s"
             % (d, NQM, " ".join(map(str, Z[:NQM])), " ".join(map(str, ZS))))
    if nsph and any(z != 7 for z in Z[NQM:]):
        stop("%s: the rows after the substrate are not the spheres (expected nitrogen, 7): %s"
             % (d, " ".join(map(str, Z[NQM:]))))
    Xall = read_xyz(d / "combined.xyz")
    dx = float(np.abs(Xb * BOHR - Xall).max())
    if dx > 1e-6:
        stop("%s: job.engrad's coordinates differ from combined.xyz by %.1e A" % (d, dx))
    g = g[:NQM]
    if abs(ee - e_tot) > 1e-9:
        stop("%s: job.engrad's energy %.9f differs from job.out's total %.9f" % (d, ee, e_tot))
    return e_tot, e_mm, g


def read_arrangements(p):
    """the designs as relax/ARRANGEMENTS.txt records them: sites, charges, a, L, the score, the sphere count."""
    txt = Path(p).read_text()
    out = {}
    for block in re.split(r"\n  (?=S6-\S|ADD-\S)", txt):
        m = re.match(r"(\S+)\s+(\d+) sphere\(s\)(.*)", block)
        if not m:
            continue
        name, n = m.group(1), int(m.group(2))
        bare = "BARE point charge" in block
        sph = [(int(a), float(b), float(c), float(d_)) for a, b, c, d_ in re.findall(
            r"site\s+(\d+)\s+q\s+([+-][\d.]+)\s+at\s+\S+\s+\S+\s+\S+\s+a\s+([+-][\d.]+)\s+L\s+([+-][\d.]+)", block)]
        if len(sph) != n:
            stop("ARRANGEMENTS.txt: %s says %d spheres and lists %d" % (name, n, len(sph)))
        s = re.search(r"score ([+-][\d.]+)", block)
        lj = 0.0 if bare else 1.0
        out[name] = {"n": n, "bare": bare, "sph": sph,
                     "score": float(s.group(1)) if s else None,
                     "aq": sum(q * a for _s, q, a, _l in sph),
                     "sl": lj * sum(l for *_x, l in sph)}
    return out


def g_add(gD, g0):
    """the criteria's RMS over all 3N Cartesian components of the ADDED gradient, in kcal/mol/A."""
    d = (gD - g0) * GCONV
    return float(np.sqrt((d ** 2).sum() / d.size))


def frozen(relax, root, out_p):
    relax, root, out = Path(relax), Path(root), Path(out_p)
    if out.exists():
        stop("%s exists" % out)
    arr = read_arrangements(relax / "ARRANGEMENTS.txt")
    # The ADD-* runs are NOT in ARRANGEMENTS.txt: it was deliberately not regenerated when they were added, so the
    # committed record would keep its proved certified gap (STAGE6_ERRATUM2.txt E6.6, relax/ADDITIVITY_RUNS.txt).
    # They are discovered from the folders instead, their site and charge read back from the folder name and then
    # CHECKED against the spheres of S6-2 and S6-3, so a folder whose name does not match a real sphere of those two
    # designs is a stop rather than something quietly measured.
    want = {}
    for nm in ("S6-2", "S6-3"):
        for site, q, a_, l_ in arr.get(nm, {"sph": []})["sph"]:
            want["ADD-s%d_q%s%05d" % (site, "p" if q >= 0 else "m", round(abs(q) * 10000))] = (site, q, a_, l_)
    for d_ in sorted(relax.glob("ADD-*")):
        tag = d_.name
        if tag not in want:
            stop("%s is not a one-sphere run of any S6-2 or S6-3 sphere" % tag)
        site, q, a_, l_ = want[tag]
        # the folder's own input must really hold that one sphere at that charge
        prms = (d_ / "R_sp" / "site.ORCAFF.prms").read_text().splitlines()
        rows = [r.split() for r in prms if len(r.split()) == 7]
        if len(rows) != NQM + 1 or abs(float(rows[-1][2]) - q) > 1e-6:
            stop("%s does not hold exactly one sphere at q %+.4f" % (tag, q))
        arr[tag] = {"n": 1, "bare": False, "sph": [(site, q, a_, l_)], "score": None,
                    "aq": q * a_, "sl": l_}
    miss = [t for t in want if t not in arr]
    if miss:
        print("   note: %d of the %d additivity runs are not present yet: %s"
              % (len(miss), len(want), " ".join(sorted(miss))))
    rep = []

    def say(s=""):
        print(s)
        rep.append(s)

    say("Stage 6, the frozen half: the design's effect and its added gradient at the fixed geometries.")
    say("cs6_analyse.py sha256 %s" % sha(__file__))
    say("relax/ARRANGEMENTS.txt sha256 %s" % sha(relax / "ARRANGEMENTS.txt"))
    say("")

    # (d) every input unchanged since it was generated
    bad = []
    for line in (relax / "INPUT_CHECKSUMS.sha256").read_text().splitlines():
        h, rel = line.split("  ", 1)
        f = relax / rel
        if not f.is_file():
            bad.append("missing: " + rel)
        elif sha(f) != h:
            bad.append("changed: " + rel)
    if bad:
        for b in bad[:10]:
            say("   " + b)
        stop("%d input files have changed or gone missing since they were generated" % len(bad))
    say("(d) every input file is unchanged since it was generated (%d entries in INPUT_CHECKSUMS.sha256)"
        % len((relax / "INPUT_CHECKSUMS.sha256").read_text().splitlines()))

    # read every frozen single point
    E, MM, G, pending = {}, {}, {}, []
    for name, d in sorted(arr.items()):
        folders = [relax / name / ("%s_sp" % st) for st in STATES]
        for f in folders:
            if not f.is_dir():
                stop("%s has no %s folder" % (name, f.name))
        if not all((f / "job.out").is_file() or (f / "job.out.xz").is_file() for f in folders):
            # An S6-* design is required; an ADD-* run may still be in the queue, and the additivity test below
            # then reports itself as not computable rather than this failing.
            if not name.startswith("ADD-"):
                stop("%s has no output yet - its frozen single points have not finished" % name)
            pending.append(name)
            continue
        for st, f in zip(STATES, folders):
            E[(name, st)], MM[(name, st)], G[(name, st)] = read_run(f, d["n"])
    done = [k for k in arr if k not in pending]
    say("(a,b,c) all %d frozen single points of %d arrangements read: terminated normally, SCF converged, %d"
        % (2 * len(done), len(done), NEL))
    say("        electrons, coordinates echoed as generated, the MM point-charge count equal to the sphere count,")
    say("        and each job.engrad holding the %d QM atoms plus its own spheres, identified by its own atom" % NQM)
    say("        block, with coordinates matching combined.xyz and an energy matching job.out's QM/MM total")
    if pending:
        say("        STILL RUNNING, so not read: %s" % " ".join(sorted(pending)))

    # (e) the control against the A-matrix job's own in vacuo energies
    ar = (root / "phase2.2" / "cs_stage6" / "A_REPORT.txt").read_text().splitlines()
    row = [l for l in ar if l.startswith("42436\t")]
    if not row:
        stop("A_REPORT.txt has no row for frame 42436")
    f = row[0].split("\t")
    want = {"R": float(f[1]), "TS": float(f[2])}
    for st in STATES:
        d_ = abs(E[("S6-C", st)] - want[st])
        if d_ > 1e-6:
            stop("the control's %s energy %.9f differs from the A-matrix job's in vacuo %.9f by %.1e Eh"
                 % (st, E[("S6-C", st)], want[st], d_))
    say("(e) the control's energies equal the A-matrix job's own in vacuo values for this frame to within"
        " %.1e Eh (R) and %.1e Eh (TS)"
        % (abs(E[("S6-C", "R")] - want["R"]), abs(E[("S6-C", "TS")] - want["TS"])))
    bare_barrier = (E[("S6-C", "TS")] - E[("S6-C", "R")]) * H
    say("    the bare substrate's barrier at these frozen geometries: %.4f kcal/mol (A_REPORT.txt: %.4f)"
        % (bare_barrier, (want["TS"] - want["R"]) * H))
    say("")

    # the measurements
    say("WHAT EACH ARRANGEMENT GIVES, at the fixed geometries")
    say("  ddE(frozen) negative lowers the barrier. 'score' is the MILP's linear estimate, sum a*q + sum L.")
    say("  The LJ columns are ORCA's own MM energy difference against the score's L term: cs5_grid.py defines")
    say("  L[i] as the sphere's Lennard-Jones energy at TS minus at R, and a QM/MM run reports exactly that as its")
    say("  (MM) energy, so the two are the same quantity measured and estimated. The control has no MM region.")
    say("  %-18s %5s %12s %12s %10s   %9s %9s   %9s %9s" % ("", "n", "ddE(frozen)", "score", "difference",
                                                            "LJ meas.", "sum L", "G_add,R", "G_add,TS"))
    res = {}
    for name in sorted(arr, key=lambda k: (k.startswith("ADD-"), arr[k]["n"], k)):
        d = arr[name]
        if name == "S6-C" or name in pending:
            continue
        dd = ((E[(name, "TS")] - E[(name, "R")]) - (E[("S6-C", "TS")] - E[("S6-C", "R")])) * H
        gr = g_add(G[(name, "R")], G[("S6-C", "R")])
        gt = g_add(G[(name, "TS")], G[("S6-C", "TS")])
        sc = d["aq"] + d["sl"]
        lj = (MM[(name, "TS")] - MM[(name, "R")]) * H
        res[name] = {"ddE": dd, "score": sc, "GR": gr, "GTS": gt, "lj": lj, "sl": d["sl"]}
        say("  %-18s %5d %+12.4f %+12.4f %+10.4f   %+9.4f %+9.4f   %9.4f %9.4f"
            % (name, d["n"], dd, sc, dd - sc, lj, d["sl"], gr, gt))
    say("")
    wl = [(n, r["lj"] - r["sl"]) for n, r in res.items() if arr[n]["n"] and not arr[n]["bare"]]
    if wl:
        w = max(wl, key=lambda t: abs(t[1]))
        say("")
        say("  THE L MATRIX, MEASURED: across the %d sphere arrangements the measured Lennard-Jones contribution"
            % len(wl))
        say("  differs from the score's sum L by at most %+.4f kcal/mol (%s). cs5_grid.py computes L classically"
            % (w[1], w[0]))
        say("  from the same Amber N parameters ORCA uses here, so agreement is a check of the grid's own")
        say("  arithmetic and of the force field it was handed, not of any physics.")
        say("")
    say("  G_add is an RMS over all %d Cartesian components of the added gradient, in kcal/mol/A. For scale only,"
        % (3 * NQM))
    say("  and not as a bound: a relaxed GOCAT design reaches about 1 kcal/mol/A and Dittner's frozen-path designs")
    say("  sat near 11.3 on the NORM, which for 24 atoms is 11.3/sqrt(72) = %.2f as an RMS." % (11.3 / np.sqrt(72)))
    say("")

    # additivity, the criteria's own test
    say("FIRST-ORDER ADDITIVITY AT FIXED GEOMETRY (STAGE6_CRITERIA.txt; the ADD-* runs)")
    for name in ("S6-2", "S6-3"):
        if name not in arr:
            continue
        parts, miss = [], []
        for site, q, _a, _l in arr[name]["sph"]:
            tag = "ADD-s%d_q%s%05d" % (site, "p" if q >= 0 else "m", round(abs(q) * 10000))
            if tag not in res:
                miss.append(tag)
            else:
                parts.append((tag, res[tag]["ddE"]))
        if miss:
            say("  %s: cannot be computed, these one-sphere runs are absent: %s" % (name, " ".join(miss)))
            continue
        s = sum(v for _t, v in parts)
        m = res[name]["ddE"]
        say("  %s measured %+.4f, sum of its one-sphere effects %+.4f, additivity error %+.4f kcal/mol (%.1f %% of"
            " the measured effect)" % (name, m, s, m - s, 100.0 * abs(m - s) / abs(m) if m else float("nan")))
        for t, v in parts:
            say("      %-18s %+.4f" % (t, v))
    say("")
    say("NOTHING HERE IS A BARRIER. These are effects at fixed geometries, the quantity the MILP's score estimates;")
    say("the catalytic verdict remains A2.5's relaxed chain, and Loss, the two survival criteria, G0 and kappa all")
    say("wait on the relaxations.")

    out.mkdir(parents=True)
    (out / "FROZEN_REPORT.txt").write_text("\n".join(rep) + "\n")
    with open(out / "frozen.tsv", "w") as fh:
        fh.write("design\tn\tddE_frozen_kcal\tscore_kcal\tdifference_kcal\tLJ_measured_kcal\tsum_L_kcal\tG_add_R\tG_add_TS\n")
        for name in sorted(res, key=lambda k: (k.startswith("ADD-"), arr[k]["n"], k)):
            r = res[name]
            fh.write("%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n"
                     % (name, arr[name]["n"], r["ddE"], r["score"], r["ddE"] - r["score"],
                        r["lj"], r["sl"], r["GR"], r["GTS"]))
    print("\nwrote %s/FROZEN_REPORT.txt and frozen.tsv" % out)
    return 0


def _fake(d, nsph, energy, grad, qmmm=True, nel=NEL, normal=True, coords=None, mmq=None, egrad_e=None,
          nat=NQM):
    """a minimal job.out/job.engrad pair in ORCA's real layout, for the selftest"""
    d = Path(d)
    d.mkdir(parents=True, exist_ok=True)
    X = np.zeros((NQM + nsph, 3)) if coords is None else coords
    (d / "combined.xyz").write_text("%d\nx\n" % len(X) + "".join(
        "C %.8f %.8f %.8f\n" % tuple(r) for r in X))
    t = []
    if qmmm:
        t.append("Multiscale model                       ... QM/MM")
        t.append("Point charges in QM calc. from MM atoms... %d" % (nsph if mmq is None else mmq))
    t.append("Number of Electrons    NEL             ....    %d" % nel)
    t.append("SCF CONVERGED AFTER   9 CYCLES")
    t.append("CARTESIAN COORDINATES (ANGSTROEM)")
    t.append("---")
    for r in X[:NQM]:
        t.append("  C   %14.8f %14.8f %14.8f" % tuple(r))
    t.append("FINAL SINGLE POINT ENERGY%s %.9f" % (" (QM/MM)" if qmmm else "", energy))
    if normal:
        t.append("****ORCA TERMINATED NORMALLY****")
    (d / "job.out").write_text("\n".join(t) + "\n")
    g = np.asarray(grad, dtype=float)
    (d / "job.engrad").write_text(
        "#\n# Number of atoms\n#\n %d\n#\n# The current total energy in Eh\n#\n  %.9f\n"
        "#\n# The current gradient in Eh/bohr\n#\n" % (nat, energy if egrad_e is None else egrad_e)
        + "".join("  %18.12f\n" % v for v in g.reshape(-1)))


def selftest():
    import tempfile
    ok = True

    def say(good, text):
        nonlocal ok
        ok &= bool(good)
        print("   [%s] %s" % ("PASS" if good else "FAIL", text))

    print("1. the arithmetic")
    g0 = np.zeros((NQM, 3))
    gD = np.zeros((NQM, 3)); gD[0, 0] = 1.0                      # one component differs by 1 Eh/bohr
    want = GCONV / np.sqrt(3 * NQM)
    say(abs(g_add(gD, g0) - want) < 1e-9,
        "G_add for a single unit component: %.4f kcal/mol/A, closed form %.4f" % (g_add(gD, g0), want))
    say(abs(g_add(g0, g0)) < 1e-12, "G_add of a design against itself is zero")
    gs = np.full((NQM, 3), 2.0)
    say(abs(g_add(gs, g0) - 2.0 * GCONV) < 1e-9,
        "a uniform difference of 2 Eh/bohr gives 2 x %.3f = %.2f kcal/mol/A" % (GCONV, g_add(gs, g0)))
    say(abs(GCONV - 627.5094740631 / 0.5291772083) < 1e-9,
        "the unit factor is the project's own two constants: %.4f kcal/mol/A per Eh/bohr" % GCONV)

    print("\n2. the readers, on fabricated ORCA output in its real layout")
    with tempfile.TemporaryDirectory() as T:
        T = Path(T)
        _fake(T / "ok", 2, -836.5, np.zeros((NQM, 3)))
        e, g = read_run(T / "ok", 2)
        say(abs(e + 836.5) < 1e-12 and g.shape == (NQM, 3),
            "a QM/MM run reads back: energy %.3f, gradient %s" % (e, g.shape))
        _fake(T / "plain", 0, -836.1, np.zeros((NQM, 3)), qmmm=False)
        e2, _ = read_run(T / "plain", 0)
        say(abs(e2 + 836.1) < 1e-12,
            "a PLAIN run reads back too, so the (QM/MM) suffix being optional is exercised both ways")

        def fires(make, nsph, what):
            d = T / ("bad_" + what.replace(" ", "_")[:22])
            make(d)
            try:
                read_run(d, nsph)
            except SystemExit:
                return True
            return False

        say(fires(lambda d: _fake(d, 2, -836.5, np.zeros((NQM, 3)), normal=False), 2, "no normal termination"),
            "stops when ORCA did not terminate normally")
        say(fires(lambda d: _fake(d, 2, -836.5, np.zeros((NQM, 3)), nel=117), 2, "wrong electrons"),
            "stops on the wrong electron count")
        say(fires(lambda d: _fake(d, 2, -836.5, np.zeros((NQM, 3)), mmq=1), 2, "wrong mm count"),
            "stops when ORCA read 1 MM point charge for a 2-sphere arrangement")
        say(fires(lambda d: _fake(d, 0, -836.1, np.zeros((NQM, 3)), qmmm=True), 0, "control as qmmm"),
            "stops if the no-sphere control ran as QM/MM")
        say(fires(lambda d: _fake(d, 2, -836.5, np.zeros((NQM, 3)), egrad_e=-836.4), 2, "energy mismatch"),
            "stops when job.engrad's energy disagrees with job.out's")
        say(fires(lambda d: _fake(d, 2, -836.5, np.zeros((25, 3)), nat=25), 2, "atom count"),
            "stops when job.engrad holds 25 atoms instead of the 24 QM atoms")

        def shifted(d):
            # job.out echoes the geometry it ran; combined.xyz is then moved under it by 0.001 A, which is what a
            # mismatch between the input and what ORCA actually used would look like.
            _fake(d, 2, -836.5, np.zeros((NQM, 3)))
            X = np.zeros((NQM + 2, 3)); X[0, 0] = 0.001
            (d / "combined.xyz").write_text("%d\nx\n" % len(X) + "".join(
                "C %.8f %.8f %.8f\n" % tuple(r) for r in X))
        say(fires(shifted, 2, "coords"), "stops when the echoed coordinates differ from combined.xyz")

    print("\n%s" % ("ALL PASS" if ok else "SOMETHING FAILED"))
    return 0 if ok else 1


def main():
    if len(sys.argv) == 2 and sys.argv[1] == "selftest":
        return selftest()
    if len(sys.argv) != 5 or sys.argv[1] != "frozen":
        sys.exit(__doc__)
    return frozen(sys.argv[2], sys.argv[3], sys.argv[4])


if __name__ == "__main__":
    sys.exit(main())
