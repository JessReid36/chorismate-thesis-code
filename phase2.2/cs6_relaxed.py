#!/usr/bin/env python3
"""
cs6_relaxed.py - Stage 6's RELAXED half (phase2.2, cs_stage6/STAGE6_CRITERIA.txt with STAGE6_AMENDMENT1.txt,
STAGE6_ERRATUM1.txt, STAGE6_ERRATUM2.txt, STAGE6_AMENDMENT3_20261103.txt and STAGE6_AMENDMENT4_20261103.txt).

  python3 cs6_relaxed.py selftest <results_repo>
  python3 cs6_relaxed.py relaxed <relax_dir> <results_repo> <OUT>

cs6_analyse.py's docstring already said this half would be a separate mode "because nothing about it can be
computed until they have run". This is that mode, in its own file rather than inside cs6_analyse.py, because
frozen/FROZEN_REPORT.txt records cs6_analyse.py's own sha256 and editing it would make that committed line false.
cs6_analyse.py and wp4_analyse.py are therefore IMPORTED and their definitions reused unaltered: the energy and
gradient readers and G_add from the first, the bond list and the superposition from the second. Their sha256 are
printed in the report so a reader knows which versions produced it.

WHAT IT COMPUTES, all of it from STAGE6_CRITERIA.txt
  ddE(frozen)    [E(TS,D) - E(R,D)] - [E(TS,0) - E(R,0)] at the FIXED geometries, from the four _sp runs.
  ddE(relaxed)   the same four energies after relaxation: R by Opt, TS by OptTS, the spheres frozen.
  Loss           ddE(frozen) - ddE(relaxed). The criteria add "with the control's own change subtracted so that
                 Loss is the design's share"; because ddE is a difference against the control at BOTH states, that
                 subtraction is already in the expression, and the report proves the identity numerically by also
                 computing Loss as (the design's own barrier change) - (the control's own barrier change) and
                 showing the two agree. Both routes are printed so neither reading is hidden.
  the tolerance  Loss below 15 % of |ddE(frozen)|; and if |ddE(frozen)| < 3.2 kcal/mol the arrangement is NOT
                 ASSESSABLE and does not enter the fit for G0, exactly as the floor clause says.
  G_add,R/TS     at the fixed geometries, from the _sp engrads, by cs6_analyse.g_add - the quantity decision O2
                 asks whether Loss can be predicted from.
  the geometry   centroid displacement, rotation and internal RMSD, the forming C1-C6 distance before and after,
                 the largest change in any substrate bond, and the closest sphere-to-substrate-atom distance, by
                 wp4_analyse.py's own definitions, which is what "as WP4 reported them" means.
  the TS         the imaginary mode count after relaxation, from the frequency block that follows the stationary
                 point, not the one computed on the way in.
  the verdicts   EFFECT and GEOMETRY, printed as two lines and never combined into one.

CHECKS, every one a stop
  (a) every run: ORCA TERMINATED NORMALLY, and for a relaxation THE OPTIMIZATION HAS CONVERGED and exactly one
      "FINAL ENERGY EVALUATION AT THE STATIONARY POINT". The scheduler's and ORCA's exit codes are READ AND
      PRINTED BUT ARE NOT CRITERIA: on 3 November 2026 four jobs that failed all carried Exit_status 0, and one
      of them reported orca_exit=0 in the same output that says "ORCA finished by error termination in LEANSCF".
      An exit code cannot license a number here; the termination banner can.
  (b) the input ORCA echoed at the top of every output equals that folder's job.inp line for line, so each number
      is attributed to the input it actually came from. check_input_echo.py's method; this reader reproduces that
      script's own answer on the committed WP4 outputs in selftest.
  (c) AMENDMENT 4's specification, read off the echoed input rather than trusted: the TS relaxations carry no
      coordsys cartesian and no ProjectTR, carry Calc_Hess true with Hybrid_Hess over the 24 substrate atoms and
      NumFreq; the reactant relaxations carry coordsys cartesian and ProjectTR false, as ERRATUM 1 set them.
  (d) the relaxation started where the frozen single point sits: each _opt run's echoed starting coordinates equal
      the matching _sp run's combined.xyz within 1e-5 A, so Loss is a difference between two states of ONE system.
  (e) the design's runs really carried their spheres and the control's really carried none: the MM point-charge
      count equals the sphere count, and the control reports no multiscale model.
  (f) the sphere count is read from the committed site.ORCAFF.prms, and for n10_A each sphere's charge is checked
      to be exactly half of S6-10's charge at the same site - AMENDMENT 3's construction, re-derived here from
      the committed files rather than taken on trust.
  (g) every committed input unchanged (relax/INPUT_CHECKSUMS.sha256), and the TS_opt2 inputs unchanged against
      the sha256 list in STAGE6_RELAX_TS_NOTE1_20261103.txt, which is where they are recorded.

NOTHING HERE IS THE CATALYTIC VERDICT. The criteria and AMENDMENT 3 both say so: the verdict is
SP_AMENDMENT2.txt A2.5's chain. These are Loss and the two survival criteria, and no barrier is quoted.

Python 3.8+, numpy.
"""
import hashlib
import importlib.util
import lzma
import re
import shutil
import sys
import tarfile
import tempfile
from pathlib import Path

import numpy as np

LOSS_FRACTION = 0.15            # STAGE6_CRITERIA.txt, THE TOLERANCE
FLOOR_KCAL = 3.2                # STAGE6_CRITERIA.txt, THE FLOOR CLAUSE
COLLAPSE_A = 1.2                # ibid., FAILURE, taken from WP4
INTEGRITY_A = 0.15              # ibid.
DESIGN = "n10_A"                # AMENDMENT 3's arrangement
CONTROL = "S6-C"
PARENT = "S6-10"                # the arrangement n10_A is half of
HALF = 0.5                      # AMENDMENT 3 / cs6_relaxA.py's SCALE


def stop(m):
    sys.exit("STOP: " + m)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _load(name):
    """Load a sibling module by explicit path. Not "import name": run.sh runs this under -I, which also implies
    -P, so the script's own directory is NOT on sys.path and a plain import would fail. An explicit path is also
    the safer form, because nothing else on the machine can shadow it."""
    p = Path(__file__).resolve().parent / (name + ".py")
    if not p.is_file():
        stop("%s is not beside this script. cs6_relaxed.py reuses cs6_analyse.py and wp4_analyse.py and must sit"
             " in the same folder as both." % p)
    spec = importlib.util.spec_from_file_location("_cs6_" + name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod, sha(p)


A, A_SHA = _load("cs6_analyse")
W4, W4_SHA = _load("wp4_analyse")

H = A.H
NQM = A.NQM
NEL = A.NEL
NAMES = A.NAMES
ELS = A.ELS
STATES = A.STATES


def text_of(d, stem="job.out"):
    d = Path(d)
    if (d / stem).is_file():
        return (d / stem).read_text(errors="replace"), str(d / stem)
    if (d / (stem + ".xz")).is_file():
        return lzma.open(str(d / (stem + ".xz")), "rt", errors="replace").read(), str(d / (stem + ".xz"))
    return None, None


def echoed_input(t):
    """the input ORCA echoed, as check_input_echo.py reads it: the lines between INPUT FILE and END OF INPUT."""
    i = t.find("INPUT FILE")
    j = t.find("****END OF INPUT****", i)
    if i < 0 or j < 0:
        return None
    e = [m.group(1).rstrip() for m in re.finditer(r"^\|\s*\d+> ?(.*)$", t[i:j], re.M)]
    while e and not e[-1]:
        e.pop()
    return e


def check_echo(d, t):
    """(b): the echoed input is this folder's job.inp, line for line."""
    d = Path(d)
    e = echoed_input(t)
    if e is None:
        stop("%s: job.out has no echoed INPUT FILE block, so its numbers cannot be attributed to an input" % d)
    want = [l.rstrip() for l in (d / "job.inp").read_text().splitlines()]
    while want and not want[-1]:
        want.pop()
    if e != want:
        for k in range(max(len(e), len(want))):
            a = e[k] if k < len(e) else "<absent>"
            b = want[k] if k < len(want) else "<absent>"
            if a != b:
                stop("%s: the input ORCA echoed differs from the committed job.inp at line %d.\n"
                     "   echoed    %r\n   committed %r" % (d, k + 1, a, b))
    return e


def freq_modes(t, after):
    """the modes of the frequency block that follows `after`. ORCA prints one block for the Hessian it computes on
    the way IN and one at the stationary point; 18c_reduced_region/ts_reduced.out has both, and they disagree -
    two imaginary modes before, one after. Taking the wrong one would be a wrong verdict, so the block is located
    by the stationary point rather than by position."""
    marks = [m.start() for m in re.finditer(r"-+\nVIBRATIONAL FREQUENCIES\n-+", t)]
    late = [i for i in marks if i > after]
    return marks, late


def parse_modes(seg):
    out = []
    for line in seg.splitlines():
        m = re.match(r"\s*(\d+):\s*(-?\d+\.\d+)\s*cm\*\*-1(.*)$", line)
        if not m:
            if out and not line.strip():
                break
            continue
        out.append((int(m.group(1)), float(m.group(2)), "imaginary mode" in m.group(3)))
    return out


def read_opt(d, nsph, is_ts):
    """every check of (a) to (e) for one relaxation; returns a dict."""
    d = Path(d)
    t, src = text_of(d)
    if t is None:
        stop("%s: no job.out or job.out.xz - has the relaxation run?" % d)

    # (a) the exit codes are read and reported, and are not used to license anything.
    exit_note = "no job.pbs.out"
    pt, _ = text_of(d, "job.pbs.out")
    if pt is not None:
        m = re.search(r"orca_exit=(\d+)", pt)
        exit_note = ("orca_exit=%s, which this reader reports and does not use" % m.group(1)) if m \
            else "job.pbs.out carries no orca_exit line"

    if "****ORCA TERMINATED NORMALLY****" not in t:
        stop("%s: ORCA did not terminate normally (%s). %s" % (d, src, exit_note))
    if "THE OPTIMIZATION HAS CONVERGED" not in t:
        stop("%s: the optimisation did not converge. %s" % (d, exit_note))
    anchors = [m.start() for m in re.finditer(r"FINAL ENERGY EVALUATION AT THE STATIONARY POINT", t)]
    if len(anchors) != 1:
        stop("%s: %d stationary-point evaluations, expected exactly 1 - an output holding more than one is not a"
             " single relaxation and may be two jobs that shared the folder" % (d, len(anchors)))
    anc = anchors[0]
    tail = t[anc:]
    mc = re.search(r"\(AFTER\s+(\d+)\s+CYCLES\)", tail)
    if not mc:
        stop("%s: the stationary-point banner does not give a cycle count" % d)
    cycles = int(mc.group(1))
    if "SCF NOT CONVERGED" in t:
        stop("%s: an SCF did not converge" % d)
    m = re.search(r"Number of Electrons\s+NEL\s+\.\.\.\.\s+(\d+)", t)
    if not m or int(m.group(1)) != NEL:
        stop("%s: %s electrons in the QM region, expected %d" % (d, m.group(1) if m else "no count", NEL))

    echo = check_echo(d, t)                                     # (b)

    # (e) the spheres, as ORCA read them
    mm = re.search(r"Point charges in QM calc\.*\s*from MM atoms\.*\s*(\d+)", t)
    if nsph == 0:
        if mm or "Multiscale model" in t:
            stop("%s: the control ran as QM/MM; it is meant to be a plain QM job" % d)
    else:
        if not mm:
            stop("%s: ORCA reports no MM point charges, but this arrangement has %d sphere(s)" % (d, nsph))
        if int(mm.group(1)) != nsph:
            stop("%s: ORCA read %s MM point charges, the arrangement has %d sphere(s)"
                 % (d, mm.group(1), nsph))

    # the energies at the stationary point, read the three apart exactly as cs6_analyse.read_run does
    qm = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", tail)
    emm = re.findall(r"FINAL SINGLE POINT ENERGY \(MM\)\s+(-?\d+\.\d+)", tail)
    tot = re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", tail)
    if nsph == 0:
        if emm or tot:
            stop("%s: the control reports MM or QM/MM energies at the stationary point" % d)
        if len(qm) != 1:
            stop("%s: %d final energies after the stationary point, expected 1" % (d, len(qm)))
        e_mm, e_tot = 0.0, float(qm[0])
    else:
        if not (len(qm) == len(emm) == len(tot) == 1):
            stop("%s: %d bare, %d (MM) and %d (QM/MM) energies after the stationary point, expected one of each"
                 % (d, len(qm), len(emm), len(tot)))
        e_mm, e_tot = float(emm[0]), float(tot[0])
        if abs(float(qm[0]) + e_mm - e_tot) > 1e-9:
            stop("%s: the QM energy plus the MM energy is not the QM/MM total (%.12f + %.12f != %.12f)"
                 % (d, float(qm[0]), e_mm, e_tot))

    # the relaxed geometry: the coordinate block at the stationary point. An Opt with ActiveAtoms writes only the
    # atoms it optimised, which here is the substrate - WP4's committed armb_freecart does exactly that - so the
    # block is required to be the 24 substrate atoms and its elements are checked against their names.
    k = tail.find("CARTESIAN COORDINATES (ANGSTROEM)")
    if k < 0:
        stop("%s: no coordinates printed at the stationary point" % d)
    rows = []
    for line in tail[k:].splitlines()[2:]:
        f = line.split()
        if len(f) == 4:
            rows.append(f)
        elif rows:
            break
    if len(rows) != NQM:
        stop("%s: the relaxed coordinate block holds %d atoms, expected the %d active substrate atoms"
             % (d, len(rows), NQM))
    if [r[0] for r in rows] != ELS:
        stop("%s: the relaxed coordinate block is not the substrate in order.\n   got      %s\n   expected %s"
             % (d, " ".join(r[0] for r in rows), " ".join(ELS)))
    X = np.array([[float(v) for v in r[1:4]] for r in rows], dtype=float)

    # the starting geometry ORCA echoed, for (d)
    k0 = t.find("CARTESIAN COORDINATES (ANGSTROEM)")
    if k0 < 0 or k0 >= anc:
        stop("%s: no starting coordinates printed before the stationary point" % d)
    X0 = np.array([[float(v) for v in l.split()[1:4]] for l in t[k0:].splitlines()[2:2 + NQM]], dtype=float)
    if X0.shape != (NQM, 3):
        stop("%s: the starting coordinate block is not %d atoms" % (d, NQM))

    modes = None
    marks, late = freq_modes(t, anc)
    if is_ts:
        if not late:
            stop("%s: no frequency block after the stationary point - a TS without frequencies cannot be tested"
                 " for exactly one imaginary mode" % d)
        if len(late) != 1:
            stop("%s: %d frequency blocks after the stationary point, expected 1" % (d, len(late)))
        modes = parse_modes(t[late[0]:late[0] + 200000])
        if not modes:
            stop("%s: the frequency block after the stationary point holds no modes this reader can parse" % d)
    else:
        if marks:
            stop("%s: a reactant relaxation printed %d frequency block(s); an Opt should print none, so this is"
                 " not the job the protocol specifies" % (d, len(marks)))
    return {"src": src, "cycles": cycles, "e_tot": e_tot, "e_mm": e_mm, "X": X, "X0": X0,
            "modes": modes, "exit_note": exit_note, "echo": echo, "nfreq_blocks": len(marks)}


def spec_check(d, echo, is_ts):
    """(c) AMENDMENT 4's specification, off the echoed input."""
    low = [l.strip().lower() for l in echo]
    kw = [l for l in echo if l.strip().startswith("!")]
    if len(kw) != 1:
        stop("%s: %d keyword lines in the echoed input, expected 1" % (d, len(kw)))
    words = kw[0].lstrip("! ").split()
    bad = []
    if is_ts:
        if "OptTS" not in words:
            bad.append("no OptTS on the keyword line")
        if "NumFreq" not in words:
            bad.append("no NumFreq on the keyword line")
        if "Freq" in words:
            bad.append("Freq is on the keyword line, which is the analytical path AMENDMENT 4 removes")
        if any(l == "coordsys cartesian" for l in low):
            bad.append("coordsys cartesian is present; AMENDMENT 4 removes it from the TS side")
        if any(l.startswith("projecttr") for l in low):
            bad.append("ProjectTR is present; AMENDMENT 4 removes it from the TS side")
        if not any(l == "calc_hess true" for l in low):
            bad.append("no Calc_Hess true")
        hh = [l for l in echo if l.strip().lower().startswith("hybrid_hess")]
        if len(hh) != 1:
            bad.append("%d Hybrid_Hess lines, expected 1" % len(hh))
        else:
            got = sorted(int(x) for x in re.findall(r"\d+", hh[0].split("{", 1)[1].split("}", 1)[0]))
            if got != list(range(NQM)):
                bad.append("Hybrid_Hess is not the %d substrate atoms: %s" % (NQM, hh[0].strip()))
    else:
        if "Opt" not in words:
            bad.append("no Opt on the keyword line")
        if "OptTS" in words:
            bad.append("OptTS on a reactant relaxation")
        if not any(l == "coordsys cartesian" for l in low):
            bad.append("no coordsys cartesian; ERRATUM 1 E6.1 sets the reactant side in Cartesians")
        if not any(l == "projecttr false" for l in low):
            bad.append("no ProjectTR false; ERRATUM 1 E6.1 sets it")
    if bad:
        stop("%s does not carry the specification:\n   %s" % (d, "\n   ".join(bad)))


def spheres_of(prms_path):
    """the sphere charges of an arrangement, from its committed force field: the rows past the 24 substrate ones."""
    rows = [r.split() for r in Path(prms_path).read_text().splitlines() if len(r.split()) == 7]
    if len(rows) <= NQM:
        stop("%s holds %d parameter rows, which is not %d substrate atoms plus at least one sphere"
             % (prms_path, len(rows), NQM))
    return [float(r[2]) for r in rows[NQM:]], np.array(
        [[float(x) for x in r[4:7]] for r in rows[NQM:]], dtype=float)


def loss_forms(bf, br):
    """every form of Loss the criteria's two lines can mean, from the four barriers. Kept as its own function so
    the arithmetic can be tested on constructed numbers whose answer is known, rather than only on the
    relaxations it is meant to judge."""
    dd_f = bf[DESIGN] - bf[CONTROL]
    dd_r = br[DESIGN] - br[CONTROL]
    own = {a: br[a] - bf[a] for a in (CONTROL, DESIGN)}
    return {"dd_f": dd_f, "dd_r": dd_r,
            # line 109's expression, exactly as written
            "loss": dd_f - dd_r,
            # line 109's own clause: the design's own barrier change with the control's subtracted. This is
            # identically dd_r - dd_f, i.e. MINUS the expression above, which is the first half of the sign
            # question STAGE6_LOSS_SIGN_FINDING_20261103.txt records.
            "loss_clause": own[DESIGN] - own[CONTROL],
            # the magnitude of effect that did not survive. Equal to loss_clause whenever both ddE are negative,
            # which criterion (i) requires anyway.
            "lost": abs(dd_f) - abs(dd_r), "own": own}


def effect_verdict(dd_f, dd_r, value):
    """criterion (i) applied to ONE reading of Loss: ddE(relaxed) still negative, and that reading below the
    tolerance; NOT ASSESSABLE under the floor clause."""
    if abs(dd_f) < FLOOR_KCAL:
        return "NOT ASSESSABLE"
    return "PASS" if (dd_r < 0 and value < LOSS_FRACTION * abs(dd_f)) else "FAIL"


def geometry(X0, X1, sph):
    """wp4_analyse.py's own metrics, by its own functions."""
    L0 = [list(map(float, p)) for p in X0]
    L1 = [list(map(float, p)) for p in X1]
    c0, c1 = W4.centroid(L0), W4.centroid(L1)
    shift = float(np.linalg.norm(np.array(c1) - np.array(c0)))
    ang, rms = W4.superpose(L0, L1)
    c16 = (W4.d(L0[0], L0[12]), W4.d(L1[0], L1[12]))          # C1 and C6, wp4_analyse.py line 129
    worst, dmax = None, 0.0
    for i, j in W4.BONDS:
        c = abs(W4.d(L1[i], L1[j]) - W4.d(L0[i], L0[j]))
        if c > dmax:
            dmax, worst = c, (i, j)
    signed = (W4.d(L1[worst[0]], L1[worst[1]]) - W4.d(L0[worst[0]], L0[worst[1]])) if worst else 0.0
    near_name, near_d = None, float("inf")
    if len(sph):
        dd = np.linalg.norm(X1[:, None, :] - sph[None, :, :], axis=2)
        a, s = np.unravel_index(int(np.argmin(dd)), dd.shape)
        near_name, near_d = "%s..sphere%d" % (NAMES[a], s + 1), float(dd[a, s])
    return {"shift": shift, "rot": ang, "rmsd": rms, "c16": c16,
            "bond": (NAMES[worst[0]] + "-" + NAMES[worst[1]]) if worst else "none",
            "dbond": dmax, "sbond": signed, "near_name": near_name, "near_d": near_d}


def stage_sp(relax, root, arm, state, tmp):
    """put one frozen single point where cs6_analyse.read_run can read it: its committed inputs, plus its outputs
    from wherever they are - loose on disk if they have just come off hpc1, or out of the committed
    frozen/sp_outputs.tar.xz if this is the control, whose 22 outputs are already in the repository."""
    src = Path(relax) / arm / ("%s_sp" % state)
    dst = Path(tmp) / arm / ("%s_sp" % state)
    dst.mkdir(parents=True)
    for f in ("combined.xyz", "job.inp", "job.pbs", "site.ORCAFF.prms"):
        if (src / f).is_file():
            shutil.copy2(src / f, dst / f)
    where = None
    if (src / "job.out").is_file() or (src / "job.out.xz").is_file():
        for f in ("job.out", "job.out.xz", "job.engrad", "job.pbs.out"):
            if (src / f).is_file():
                shutil.copy2(src / f, dst / f)
        where = "loose files in %s" % src
    else:
        tar = Path(root) / "phase2.2" / "cs_stage6" / "frozen" / "sp_outputs.tar.xz"
        if not tar.is_file():
            stop("%s/%s_sp has no output and %s does not exist" % (arm, state, tar))
        pre = "./%s/%s_sp/" % (arm, state)
        got = 0
        with tarfile.open(tar, "r:xz") as th:
            for mem in th.getmembers():
                if mem.name.startswith(pre) and mem.isfile():
                    data = th.extractfile(mem).read()
                    (dst / Path(mem.name).name).write_bytes(data)
                    got += 1
        if not got:
            stop("%s holds nothing under %s" % (tar, pre))
        where = "%d file(s) from the committed frozen/sp_outputs.tar.xz" % got
    return dst, where


def report(say, out, nsph, E, MM, G, R, JOBS, xyz_d, rep):
    """Everything after the reading: the quantities, the tolerance, the geometry, the two verdicts and
    the two files. Separated from relaxed() so the selftest can drive it with CONSTRUCTED results and
    check the figures and the verdicts without fabricating ORCA outputs - the glue gets tested, not just
    the parts."""
    # ---- the quantities
    say("")
    say("WHAT THE TWO STATES GIVE")
    bf = {a: (E[(a, "TS")] - E[(a, "R")]) * H for a in (CONTROL, DESIGN)}
    br = {a: (R[(a, "TS")]["e_tot"] - R[(a, "R")]["e_tot"]) * H for a in (CONTROL, DESIGN)}
    F = loss_forms(bf, br)
    dd_f, dd_r, loss, loss2, own = F["dd_f"], F["dd_r"], F["loss"], F["loss_clause"], F["own"]
    if abs(loss + loss2) > 1e-9:
        stop("line 109's expression and its own clause are not each other's negative (%+.9f and %+.9f), which"
             " cannot happen if both are arithmetic on the same eight energies" % (loss, loss2))
    say("  barrier at the frozen geometries   %-8s %+10.4f   %-8s %+10.4f  kcal/mol"
        % (CONTROL, bf[CONTROL], DESIGN, bf[DESIGN]))
    say("  barrier after relaxation           %-8s %+10.4f   %-8s %+10.4f"
        % (CONTROL, br[CONTROL], DESIGN, br[DESIGN]))
    say("  each arm's own change              %-8s %+10.4f   %-8s %+10.4f"
        % (CONTROL, own[CONTROL], DESIGN, own[DESIGN]))
    say("")
    lost = F["lost"]
    say("  ddE(frozen)   %+10.4f kcal/mol" % dd_f)
    say("  ddE(relaxed)  %+10.4f" % dd_r)
    say("  THREE READINGS OF Loss, because STAGE6_CRITERIA.txt supports three and they do not agree in sign.")
    say("  (1) as line 109 writes it      %+10.4f   = ddE(frozen) - ddE(relaxed)" % loss)
    say("  (2) as line 109's clause says  %+10.4f   = the design's own barrier change with the control's"
        % loss2)
    say("                                              subtracted; identically -(1), to %.1e"
        % abs(loss + loss2))
    say("  (3) the magnitude not kept     %+10.4f   = |ddE(frozen)| - |ddE(relaxed)|" % lost)
    agree23 = abs(loss2 - lost) <= 1e-9
    say("  (2) and (3) %s here%s, and both are positive exactly when the design loses effect, which is what"
        % ("agree" if agree23 else "DISAGREE", " (they coincide while both ddE are negative)" if agree23
           else ""))
    say("  line 139's tolerance has to mean. (1) is their negative. SEE STAGE6_LOSS_SIGN_FINDING_20261103.txt:")
    say("  two of the three readings agree with each other and with the stated intent, and the bare expression")
    say("  of line 109 is the one out of step - but correcting a pre-registered formula is an amendment, so")
    say("  this reader reports all three and decides none.")
    if abs(dd_f) > 0:
        say("  retained      %9.2f %%   = |ddE(relaxed)| / |ddE(frozen)|" % (100.0 * abs(dd_r) / abs(dd_f)))
    say("")
    gr = A.g_add(G[(DESIGN, "R")], G[(CONTROL, "R")])
    gt = A.g_add(G[(DESIGN, "TS")], G[(CONTROL, "TS")])
    say("  G_add,R %9.4f    G_add,TS %9.4f   kcal/mol/A, at the FIXED geometries, by cs6_analyse.g_add"
        % (gr, gt))
    say("")

    # ---- the tolerance and its floor
    say("THE TOLERANCE (STAGE6_CRITERIA.txt: Loss below %.0f %% of |ddE(frozen)|, and the floor clause)"
        % (100 * LOSS_FRACTION))
    assessable = abs(dd_f) >= FLOOR_KCAL
    allowed = LOSS_FRACTION * abs(dd_f)
    neg = dd_r < 0
    if not assessable:
        say("  |ddE(frozen)| = %.4f kcal/mol is below the %.1f floor, so this arrangement's survival is NOT"
            % (abs(dd_f), FLOOR_KCAL))
        say("  ASSESSABLE and it does not enter the fit for G0. The figures above are still reported.")
        eff_written = eff_lost = "NOT ASSESSABLE"
    else:
        say("  |ddE(frozen)| = %.4f, so the tolerance allows at most %.4f kcal/mol." % (abs(dd_f), allowed))
        say("  ddE(relaxed) is %s, which criterion (i) requires to be negative."
            % ("negative" if neg else "NOT negative"))
        eff_written = effect_verdict(dd_f, dd_r, loss)
        eff_lost = effect_verdict(dd_f, dd_r, lost)
        say("  under the expression as written      Loss %+.4f against %.4f   ->  %s"
            % (loss, allowed, eff_written))
        say("  under the magnitude of effect lost   %+.4f against %.4f   ->  %s"
            % (lost, allowed, eff_lost))
        if eff_written != eff_lost:
            say("  THE TWO READINGS DISAGREE. STAGE6_LOSS_SIGN_FINDING_20261103.txt records why and does not")
            say("  resolve it; resolving it is an amendment to the criteria and is not this reader's to make.")
    say("")

    # ---- the geometry
    say("THE RELAXED GEOMETRY, by wp4_analyse.py's own definitions")
    geo, gfail = {}, []
    for arm in (CONTROL, DESIGN):
        sph = xyz_d if arm == DESIGN else np.zeros((0, 3))
        for st in STATES:
            r = R[(arm, st)]
            g = geometry(r["start"], r["X"], sph)
            geo[(arm, st)] = g
            say("  %-8s %-3s  centroid %6.3f A, rotated %6.2f deg, internal RMSD %6.3f A; C1..C6 %.3f -> %.3f A"
                % (arm, st, g["shift"], g["rot"], g["rmsd"], g["c16"][0], g["c16"][1]))
            say("                largest bond change %-8s %+.3f A%s"
                % (g["bond"], g["sbond"], "" if g["dbond"] <= INTEGRITY_A else "   BROKEN INTEGRITY"))
            if g["near_name"]:
                say("                closest sphere contact %-16s %.3f A%s"
                    % (g["near_name"], g["near_d"], "" if g["near_d"] >= COLLAPSE_A else "   COLLAPSE"))
            if g["dbond"] > INTEGRITY_A:
                gfail.append("%s %s: a bond moved %.3f A, over the %.2f limit" % (arm, st, g["dbond"],
                                                                                  INTEGRITY_A))
            if g["near_name"] and g["near_d"] < COLLAPSE_A:
                gfail.append("%s %s: a sphere sits %.3f A from the substrate, under the %.1f limit"
                             % (arm, st, g["near_d"], COLLAPSE_A))
    say("  the design's C1..C6 change against the control's: %+.3f A against %+.3f A"
        % (geo[(DESIGN, "TS")]["c16"][1] - geo[(DESIGN, "TS")]["c16"][0],
           geo[(CONTROL, "TS")]["c16"][1] - geo[(CONTROL, "TS")]["c16"][0]))
    say("")

    say("THE TRANSITION STATES AFTER RELAXATION")
    for arm in (CONTROL, DESIGN):
        md = R[(arm, "TS")]["modes"]
        im = [(i, v) for i, v, f in md if f]
        say("  %-8s %d mode(s) parsed at the stationary point, %d imaginary%s"
            % (arm, len(md), len(im), (": " + ", ".join("%.2f cm-1" % v for _i, v in im)) if im else ""))
        if len(im) != 1:
            gfail.append("%s TS: %d imaginary modes, not the one the criteria require" % (arm, len(im)))
    say("  The block read is the one AFTER the stationary point. 18c_reduced_region/ts_reduced.out prints two")
    say("  blocks that disagree - two imaginary modes on the way in, one at the end - so which block is read")
    say("  decides the verdict, and this reader locates it by the stationary point rather than by position.")
    say("")

    say("THE TWO SURVIVAL CRITERIA, REPORTED SEPARATELY")
    if eff_written == eff_lost:
        say("  (i)  EFFECT   : %s   (both readings of Loss agree)" % eff_written)
    else:
        say("  (i)  EFFECT   : %s under the expression as written, %s under the magnitude of effect lost."
            % (eff_written, eff_lost))
        say("                  NOT DECIDED HERE. One of the two is the criteria's intent and the other is a")
        say("                  defect in it; which, is an amendment Jessica makes.")
    say("  (ii) GEOMETRY : %s" % ("PASS" if not gfail else "FAIL"))
    for f in gfail:
        say("                  %s" % f)
    say("")
    say("NEITHER IS THE CATALYTIC VERDICT. STAGE6_CRITERIA.txt and AMENDMENT 3 both place that in")
    say("SP_AMENDMENT2.txt A2.5's chain - an optimised transition state with one counted imaginary mode, an IRC")
    say("in both directions, the ends minimised, and B against the reference by the same chain in the same")
    say("medium. No barrier is quoted here, and one relaxed design cannot fix G0 or kappa: the criteria's fit")
    say("needs the arrangements that meet both criteria, and this adds one candidate to that set.")

    out.mkdir(parents=True)
    (out / "RELAXED_REPORT.txt").write_text("\n".join(rep) + "\n")
    with open(out / "relaxed.tsv", "w") as fh:
        fh.write("design\tn\tddE_frozen_kcal\tddE_relaxed_kcal\tloss_as_written_kcal\teffect_lost_kcal\t"
                 "loss_allowed_kcal\tretained_pct\tassessable\tG_add_R\tG_add_TS\teffect_as_written\t"
                 "effect_lost_reading\tgeometry\tcycles_R\tcycles_TS\timag_TS\timag_control_TS\n")
        fh.write("%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f\t%s\t%.6f\t%.6f\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n"
                 % (DESIGN, nsph, dd_f, dd_r, loss, lost, allowed,
                    100.0 * abs(dd_r) / abs(dd_f) if dd_f else float("nan"),
                    "yes" if assessable else "no", gr, gt, eff_written, eff_lost,
                    "PASS" if not gfail else "FAIL", R[(DESIGN, "R")]["cycles"], R[(DESIGN, "TS")]["cycles"],
                    sum(1 for _i, _v, f in R[(DESIGN, "TS")]["modes"] if f),
                    sum(1 for _i, _v, f in R[(CONTROL, "TS")]["modes"] if f)))
    print("\nwrote %s/RELAXED_REPORT.txt and relaxed.tsv" % out)
    return 0


def relaxed(relax_p, root_p, out_p):
    relax, root, out = Path(relax_p), Path(root_p), Path(out_p)
    if out.exists():
        stop("%s exists" % out)
    cs6 = root / "phase2.2" / "cs_stage6"
    rep = []

    def say(s=""):
        print(s)
        rep.append(s)

    say("Stage 6, the relaxed half: Loss, the two survival criteria, and the geometry of the relaxed points.")
    say("cs6_relaxed.py sha256 %s" % sha(__file__))
    say("cs6_analyse.py  sha256 %s   (imported: the energy and gradient readers, and G_add)" % A_SHA)
    say("wp4_analyse.py  sha256 %s   (imported: the bond list and the superposition)" % W4_SHA)
    say("")

    # ---- (g) the inputs are what they were
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
        stop("%d committed input files have changed or gone missing" % len(bad))
    n_ck = len((relax / "INPUT_CHECKSUMS.sha256").read_text().splitlines())
    note = cs6 / "STAGE6_RELAX_TS_NOTE1_20261103.txt"
    if not note.is_file():
        stop("%s is absent; the TS_opt2 inputs' checksums live in it" % note)
    ts2 = [(m.group(1), m.group(2)) for m in
           re.finditer(r"^(relax/\S*TS_opt2/\S+)\s+([0-9a-f]{64})\s*$", note.read_text(), re.M)]
    if len(ts2) != 7:
        stop("%s lists %d TS_opt2 files, expected 7" % (note.name, len(ts2)))
    for rel, h in ts2:
        f = relax.parent / rel
        if not f.is_file():
            stop("%s is absent" % rel)
        if sha(f) != h:
            stop("%s differs from the sha256 NOTE 1 records for it" % rel)
    say("(g) %d committed inputs unchanged, and the 7 TS_opt2 files match the sha256 NOTE 1 records for them"
        % n_ck)

    # ---- (f) the design's spheres, from its own force field, re-derived against its parent
    q_d, xyz_d = spheres_of(relax / DESIGN / "R_sp" / "site.ORCAFF.prms")
    q_p, xyz_p = spheres_of(relax / PARENT / "R_sp" / "site.ORCAFF.prms")
    nsph = len(q_d)
    if len(q_p) != nsph:
        stop("%s has %d spheres and %s has %d; AMENDMENT 3 builds the first from the second"
             % (DESIGN, nsph, PARENT, len(q_p)))
    if np.abs(xyz_d - xyz_p).max() > 1e-9:
        stop("%s's spheres are not at %s's positions (largest difference %.1e A)"
             % (DESIGN, PARENT, np.abs(xyz_d - xyz_p).max()))
    worst = max(abs(a - HALF * b) for a, b in zip(q_d, q_p))
    if worst > 1e-9:
        stop("%s's charges are not exactly %.2f of %s's: largest difference %.3e e"
             % (DESIGN, HALF, PARENT, worst))
    for st in STATES:
        q2, _x2 = spheres_of(relax / DESIGN / ("%s_sp" % st) / "site.ORCAFF.prms")
        if q2 != q_d:
            stop("%s's %s force field does not carry the same charges as its R one" % (DESIGN, st))
    say("(f) %s carries %d spheres, each at %s's site and at exactly %.2f of its charge (largest departure"
        " %.1e e), and R and TS carry the same field" % (DESIGN, nsph, PARENT, HALF, worst))

    # ---- the inventory, before anything is read: all eight outputs, and where each one is
    JOBS = {(CONTROL, "R"): ("R_opt", False), (CONTROL, "TS"): ("TS_opt2", True),
            (DESIGN, "R"): ("R_opt", False), (DESIGN, "TS"): ("TS_opt2", True)}
    tarp = cs6 / "frozen" / "sp_outputs.tar.xz"
    intar = set()
    if tarp.is_file():
        with tarfile.open(tarp, "r:xz") as th:
            intar = {m.name.split("/")[1] + "/" + m.name.split("/")[2] for m in th.getmembers()
                     if m.isfile() and m.name.count("/") >= 3}
    inv, absent = [], []
    for arm in (CONTROL, DESIGN):
        for st in STATES:
            d = relax / arm / ("%s_sp" % st)
            key = "%s/%s_sp" % (arm, st)
            if text_of(d)[0] is not None:
                inv.append((key, "downloaded, loose in the repository"))
            elif key in intar:
                inv.append((key, "committed, inside frozen/sp_outputs.tar.xz"))
            else:
                inv.append((key, "ABSENT - still on hpc1, not downloaded"))
                absent.append(key)
    for (arm, st), (folder, _i) in sorted(JOBS.items()):
        key = "%s/%s" % (arm, folder)
        if text_of(relax / arm / folder)[0] is not None:
            inv.append((key, "downloaded, loose in the repository"))
        else:
            inv.append((key, "ABSENT - running, queued or not downloaded"))
            absent.append(key)
    say("")
    say("THE EIGHT OUTPUTS Loss NEEDS, and where each one is")
    for k, w in inv:
        say("      %-18s %s" % (k, w))
    if absent:
        say("")
        say("%d of the 8 are absent, so there is nothing to compute. Loss is a difference between two states of"
            % len(absent))
        say("two arms; from seven corners it would be a guess, and this reader does not report a partial Loss.")
        say("Download the absent outputs - every output file of a finished job, by the rule of 18 October 2026 -")
        say("and run this again.")
        out.mkdir(parents=True)
        (out / "RELAXED_REPORT.txt").write_text("\n".join(rep) + "\n")
        print("\nwrote %s/RELAXED_REPORT.txt (an inventory only: %d of 8 outputs present)"
              % (out, 8 - len(absent)))
        return 0
    say("      all eight are present")

    # ---- the frozen half, by cs6_analyse's own reader
    tmp = Path(tempfile.mkdtemp(prefix="cs6_relaxed_"))
    E, MM, G, prov = {}, {}, {}, []
    try:
        for arm, ns in ((CONTROL, 0), (DESIGN, nsph)):
            for st in STATES:
                d, where = stage_sp(relax, root, arm, st, tmp)
                prov.append("      %-8s %-3s  %s" % (arm, st, where))
                E[(arm, st)], MM[(arm, st)], G[(arm, st)] = A.read_run(d, ns)
                t, _ = text_of(d)
                check_echo(d, t)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    say("(a,b,c,e) the four frozen single points read by cs6_analyse.read_run - terminated normally, SCF")
    say("      converged, %d electrons, the MM point-charge count equal to the sphere count, each engrad" % NEL)
    say("      self-describing and matching combined.xyz - and each output's echoed input equal to its job.inp")
    for p in prov:
        say(p)

    # ---- the relaxed half
    R = {}
    for (arm, st), (folder, is_ts) in sorted(JOBS.items()):
        ns = 0 if arm == CONTROL else nsph
        d = relax / arm / folder
        r = read_opt(d, ns, is_ts)
        spec_check(d, r["echo"], is_ts)
        # (d) it started where the frozen point sits
        C = A.read_xyz(relax / arm / ("%s_sp" % st) / "combined.xyz")[:NQM]
        dx = float(np.abs(r["X0"] - C).max())
        if dx > 1e-5:
            stop("%s started %.1e A from %s/%s_sp's geometry; Loss must compare two states of one system"
                 % (d, dx, arm, st))
        r["start"] = C
        r["d_start"] = dx
        R[(arm, st)] = r
    say("")
    say("(a,b,c,d) the four relaxations read: terminated normally, converged, exactly one stationary-point")
    say("      evaluation each, the echoed input equal to the committed job.inp, AMENDMENT 4's specification")
    say("      present on the TS side and ERRATUM 1's Cartesians on the reactant side, and each one started")
    say("      at its own frozen geometry.")
    say("      %-8s %-9s %7s  %-11s  %s" % ("", "folder", "cycles", "start diff", "exit code, reported not used"))
    for (arm, st), (folder, _i) in sorted(JOBS.items()):
        r = R[(arm, st)]
        say("      %-8s %-9s %7d  %9.1e A  %s" % (arm, folder, r["cycles"], r["d_start"], r["exit_note"]))

    return report(say, out, nsph, E, MM, G, R, JOBS, xyz_d, rep)
    return 0


# --------------------------------------------------------------------------------------------------------------
# selftest: every check is made to fire, on CONTROLS - committed outputs of jobs that are not the ones this reader
# licenses - never on the relaxations themselves.
# --------------------------------------------------------------------------------------------------------------
def _fires(fn, label, expect):
    """A negative must stop FOR ITS OWN REASON. `expect` is a fragment of the message this doctoring should
    produce; a stop that fires on something else is a failure, not a pass, because it proves nothing about the
    check it was meant to exercise."""
    try:
        fn()
    except SystemExit as e:
        m = str(e)
        if not m.startswith("STOP:"):
            return False, "   DID NOT STOP (exited %r): %s" % (m, label)
        first = m[6:].splitlines()[0]
        if expect not in m:
            return False, "   WRONG REASON: %-44s wanted %r, got %r" % (label, expect, first[:70])
        return True, "   fires: %-56s %s" % (label, first[:80])
    return False, "   DID NOT STOP: %s" % label


def selftest(root_p):
    root = Path(root_p)
    ok = True
    out = []

    def say(good, line):
        nonlocal ok
        ok = ok and good
        print(line)
        out.append(line)

    wp4 = root / "phase2.2" / "wp4_s13_rerun"
    red = root / "05_qmmm" / "18c_reduced_region" / "ts_reduced.out"
    for p in (wp4 / "armb_freecart" / "job.out.xz", red):
        if not p.is_file():
            sys.exit("STOP: the selftest needs the committed control output %s" % p)

    print("cs6_relaxed.py selftest - every check made to fire on committed outputs of OTHER jobs")
    print("cs6_analyse.py sha256 %s" % A_SHA)
    print("wp4_analyse.py sha256 %s" % W4_SHA)
    print()

    t_opt = lzma.open(str(wp4 / "armb_freecart" / "job.out.xz"), "rt", errors="replace").read()
    t_ts = red.read_text(errors="replace")

    # 1 the echo reader reproduces check_input_echo.py's own answer on a committed pair
    e = echoed_input(t_opt)
    want = [l.rstrip() for l in (wp4 / "armb_freecart" / "job.inp").read_text().splitlines()]
    while want and not want[-1]:
        want.pop()
    say(e == want, "   %s  the echoed-input reader returns armb_freecart's committed job.inp exactly"
        % ("ok  " if e == want else "FAIL"))

    # 2 the stationary-point anchor and its cycle count, against WP4's committed table
    row = [l.split("\t") for l in (wp4 / "wp4_results.tsv").read_text().splitlines()
           if l.startswith("armb_freecart\t")]
    anc = t_opt.find("FINAL ENERGY EVALUATION AT THE STATIONARY POINT")
    cyc = int(re.search(r"\(AFTER\s+(\d+)\s+CYCLES\)", t_opt[anc:]).group(1))
    say(cyc == int(row[0][1]), "   %s  the cycle count at the stationary point is %d, which is wp4_results.tsv's"
        % ("ok  " if cyc == int(row[0][1]) else "FAIL", cyc))

    # 3 the frequency block is located by the stationary point, not by position
    marks, late = freq_modes(t_ts, t_ts.find("FINAL ENERGY EVALUATION AT THE STATIONARY POINT"))
    first = [(i, v) for i, v, f in parse_modes(t_ts[marks[0]:marks[0] + 200000]) if f]
    last = [(i, v) for i, v, f in parse_modes(t_ts[late[0]:late[0] + 200000]) if f]
    good = len(marks) == 2 and len(late) == 1 and len(first) == 2 and len(last) == 1
    say(good, "   %s  ts_reduced.out: %d blocks, %d after the stationary point; %d imaginary before and %d after"
        % ("ok  " if good else "FAIL", len(marks), len(late), len(first), len(last)))
    g2 = abs(last[0][1] + 313.07) < 5e-3 and abs(first[0][1] + 313.30) < 5e-3
    say(g2, "   %s  and the modes read are %.2f before and %.2f after, the committed pair"
        % ("ok  " if g2 else "FAIL", first[0][1], last[0][1]))

    # 4 the geometry metrics reproduce WP4's committed row for armb_freecart
    X0 = A.read_xyz(wp4 / "armb_freecart" / "combined.xyz")[:NQM]
    k = t_opt.find("CARTESIAN COORDINATES (ANGSTROEM)", anc)
    rows = []
    for line in t_opt[k:].splitlines()[2:]:
        f = line.split()
        if len(f) == 4:
            rows.append(f)
        elif rows:
            break
    X1 = np.array([[float(v) for v in r[1:4]] for r in rows], dtype=float)
    g = geometry(X0, X1, np.zeros((0, 3)))
    w = {k2: float(v) for k2, v in zip((wp4 / "wp4_results.tsv").read_text().splitlines()[0].split("\t"),
                                       row[0]) if _isnum(v)}
    checks = [("centroid_shift", g["shift"], 5e-4), ("rotation", g["rot"], 5e-3),
              ("internal_rmsd", g["rmsd"], 5e-4), ("C1C6_start", g["c16"][0], 5e-4),
              ("C1C6_final", g["c16"][1], 5e-4), ("max_bond_change", g["dbond"], 5e-4)]
    for name, got, tol in checks:
        good = name in w and abs(got - w[name]) <= tol
        say(good, "   %s  %-16s %10.4f against wp4_results.tsv's %10.4f"
            % ("ok  " if good else "FAIL", name, got, w.get(name, float("nan"))))

    # 5 the Loss arithmetic, on constructed barriers whose answer is known by hand. This is where the criteria's
    #   two lines are shown to disagree: the expression of line 109 and the tolerance of line 139 give opposite
    #   verdicts on the same relaxation, which STAGE6_LOSS_SIGN_FINDING_20261103.txt records.
    print()
    print("   the Loss arithmetic, on constructed barriers (control frozen 0, control relaxed 0):")
    print("   %-26s %8s %8s %9s %9s %8s | %-14s %-14s"
          % ("case", "ddE_f", "ddE_r", "as writ.", "lost", "allowed", "verdict written", "verdict lost"))
    CASES = [("keeps a tenth of it", -20.0, -2.0, "PASS", "FAIL"),
             ("keeps nine tenths", -20.0, -18.0, "PASS", "PASS"),
             ("effect grows by half", -20.0, -30.0, "FAIL", "PASS"),
             ("effect reverses sign", -20.0, +5.0, "FAIL", "FAIL"),
             ("under the 3.2 floor", -2.0, -0.2, "NOT ASSESSABLE", "NOT ASSESSABLE")]
    for name, f_, r_, w_want, l_want in CASES:
        # barriers chosen so that the design's frozen barrier is f_ and its relaxed barrier is r_, with the
        # control flat at zero: then ddE(frozen) = f_ and ddE(relaxed) = r_ exactly.
        F = loss_forms({CONTROL: 0.0, DESIGN: f_}, {CONTROL: 0.0, DESIGN: r_})
        vw = effect_verdict(F["dd_f"], F["dd_r"], F["loss"])
        vl = effect_verdict(F["dd_f"], F["dd_r"], F["lost"])
        good = (abs(F["dd_f"] - f_) < 1e-12 and abs(F["dd_r"] - r_) < 1e-12
                and abs(F["loss"] + F["loss_clause"]) < 1e-12          # (1) is minus (2), always
                and (r_ >= 0 or abs(F["loss_clause"] - F["lost"]) < 1e-12)   # (2) is (3) while both are negative
                and vw == w_want and vl == l_want)
        say(good, "   %s %-26s %+8.2f %+8.2f %+9.2f %+9.2f %8.2f | %-14s %-14s"
            % ("ok " if good else "FAIL", name, F["dd_f"], F["dd_r"], F["loss"], F["lost"],
               LOSS_FRACTION * abs(F["dd_f"]), vw, vl))
    dis = sum(1 for n, f_, r_, w, l in CASES if w != l)
    say(dis == 2, "   %s  the two readings disagree on %d of the %d constructed cases, which is the defect the"
                  " finding records" % ("ok " if dis == 2 else "FAIL", dis, len(CASES)))

    # 6 report() end to end, on CONSTRUCTED results whose answer is known by hand. This tests the glue between
    #   the readers and the report - the dictionaries, the formatting and the two files - without fabricating
    #   ORCA outputs, and it is run on numbers that are not any relaxation's.
    print()
    print("   report() on constructed results:")
    X0 = A.read_xyz(root / "phase2.2" / "cs_stage6" / "relax" / CONTROL / "R_sp" / "combined.xyz")[:NQM]
    gc = 10.0 / (A.H / A.BOHR)                      # an added gradient of exactly 10 kcal/mol/A per component
    JOBS = {(CONTROL, "R"): ("R_opt", False), (CONTROL, "TS"): ("TS_opt2", True),
            (DESIGN, "R"): ("R_opt", False), (DESIGN, "TS"): ("TS_opt2", True)}

    def mk(bf_c, bf_d, br_c, br_d, modes_d, sph_at, X1=None):
        E = {(CONTROL, "R"): 0.0, (CONTROL, "TS"): bf_c / A.H,
             (DESIGN, "R"): 0.0, (DESIGN, "TS"): bf_d / A.H}
        MM = {k: 0.0 for k in E}
        G = {(CONTROL, s): np.zeros((NQM, 3)) for s in STATES}
        G.update({(DESIGN, s): np.full((NQM, 3), gc) for s in STATES})
        R = {}
        for arm, st, e in ((CONTROL, "R", 0.0), (CONTROL, "TS", br_c / A.H),
                           (DESIGN, "R", 0.0), (DESIGN, "TS", br_d / A.H)):
            R[(arm, st)] = {"cycles": 42, "e_tot": e, "e_mm": 0.0,
                            "X": X0.copy() if X1 is None else X1, "X0": X0.copy(), "start": X0.copy(),
                            "d_start": 0.0, "exit_note": "constructed", "nfreq_blocks": 1,
                            "modes": ([(6, -300.0, True)] + [(7, 10.0, False)]) if st == "TS" else None}
        R[(DESIGN, "TS")]["modes"] = modes_d
        return E, MM, G, R, np.full((10, 3), sph_at)

    def run(tag, *a, **k):
        d = tmp0 / tag
        lines = []
        E, MM, G, R, sph = mk(*a, **k)
        report(lambda s="": lines.append(s), d, 10, E, MM, G, R, JOBS, sph, lines)
        row = dict(zip((d / "relaxed.tsv").read_text().splitlines()[0].split("\t"),
                       (d / "relaxed.tsv").read_text().splitlines()[1].split("\t")))
        return row

    tmp0 = Path(tempfile.mkdtemp(prefix="cs6_relaxed_rep_"))
    try:
        # case A: frozen effect -20, relaxed effect -2, so a tenth survives; spheres 100 A away; one imaginary
        rowA = run("A", 18.0, -2.0, 16.0, 14.0, [(6, -300.0, True), (7, 9.0, False)], 100.0)
        want = {"ddE_frozen_kcal": -20.0, "ddE_relaxed_kcal": -2.0, "loss_as_written_kcal": -18.0,
                "effect_lost_kcal": 18.0, "loss_allowed_kcal": 3.0, "retained_pct": 10.0,
                "G_add_R": 10.0, "G_add_TS": 10.0}
        for k, v in want.items():
            good = abs(float(rowA[k]) - v) < 5e-4
            say(good, "   %s  %-22s %+12.4f, constructed to be %+12.4f"
                % ("ok " if good else "FAIL", k, float(rowA[k]), v))
        for k, v in (("assessable", "yes"), ("effect_as_written", "PASS"),
                     ("effect_lost_reading", "FAIL"), ("geometry", "PASS"),
                     ("imag_TS", "1"), ("imag_control_TS", "1")):
            good = rowA[k] == v
            say(good, "   %s  %-22s %-6s, expected %s" % ("ok " if good else "FAIL", k, rowA[k], v))

        # case B: a sphere 0.9 A from the substrate and two imaginary modes - GEOMETRY must fail, twice
        rowB = run("B", 18.0, -2.0, 16.0, 14.0,
                   [(6, -300.0, True), (7, -18.0, True)], float(X0[0, 0]) + 0.9)
        good = rowB["geometry"] == "FAIL" and rowB["imag_TS"] == "2"
        say(good, "   %s  a collapsed sphere and two imaginary modes give GEOMETRY %s with imag_TS %s"
            % ("ok " if good else "FAIL", rowB["geometry"], rowB["imag_TS"]))
    finally:
        shutil.rmtree(tmp0, ignore_errors=True)

    # 7 every stop fires, each on a doctored copy of a committed output
    print()
    print("   the stops, each on a doctored copy of a committed control output:")
    tmp = Path(tempfile.mkdtemp(prefix="cs6_relaxed_st_"))
    try:
        def make(text, inp=None, folder="w"):
            d = tmp / folder
            if d.exists():
                shutil.rmtree(d)
            d.mkdir(parents=True)
            (d / "job.out").write_text(text)
            src = inp if inp is not None else (wp4 / "armb_freecart" / "job.inp").read_text()
            (d / "job.inp").write_text(src)
            return d

        inp_ts = "\n".join(echoed_input(t_ts)) + "\n"            # ts_reduced's own input, from its own echo
        # the final coordinate block with one atom row removed, and nothing else touched
        kfin = t_opt.find("CARTESIAN COORDINATES (ANGSTROEM)", anc)
        head, blk = t_opt[:kfin], t_opt[kfin:]
        bl = blk.splitlines(True)
        short = "".join(bl[:2 + NQM - 1] + bl[2 + NQM:])
        say(*_fires(lambda: read_opt(make(t_opt.replace("****ORCA TERMINATED NORMALLY****", "")), 1, False),
                    "no termination banner", "did not terminate normally"))
        say(*_fires(lambda: read_opt(make(t_opt.replace("THE OPTIMIZATION HAS CONVERGED", "x")), 1, False),
                    "the optimisation did not converge", "did not converge"))
        say(*_fires(lambda: read_opt(make(t_opt + t_opt), 1, False),
                    "two stationary-point evaluations in one output", "stationary-point evaluations"))
        say(*_fires(lambda: read_opt(make(t_opt.replace("NEL             ....  118",
                                                        "NEL             ....  117")), 1, False),
                    "the wrong electron count", "electrons in the QM region"))
        say(*_fires(lambda: read_opt(make(t_opt, inp="! wrong\n"), 1, False),
                    "the echoed input is not job.inp", "differs from the committed job.inp"))
        say(*_fires(lambda: read_opt(make(t_opt), 0, False),
                    "a control that ran as QM/MM", "ran as QM/MM"))
        say(*_fires(lambda: read_opt(make(t_opt), 99, False),
                    "the MM point-charge count is not the sphere count", "MM point charges, the arrangement has"))
        say(*_fires(lambda: read_opt(make(t_opt), 1, True),
                    "a TS with no frequency block after the stationary point", "no frequency block after"))
        # a valid reactant output with ts_reduced.out's own frequency block spliced on: everything else passes, so
        # the only thing left to stop it is the check that an Opt prints no frequencies.
        fb = t_ts[marks[-1]:marks[-1] + 1200]
        say(*_fires(lambda: read_opt(make(t_opt + "\n" + fb), 1, False),
                    "a reactant relaxation that printed frequencies", "should print none"))
        say(*_fires(lambda: read_opt(make(head + short), 1, False),
                    "the relaxed coordinate block is one atom short", "coordinate block holds 23 atoms"))
        e_ok = echoed_input(t_opt)
        say(*_fires(lambda: spec_check("w", [l for l in e_ok if "cartesian" not in l.lower()], False),
                    "a reactant input with no coordsys cartesian", "no coordsys cartesian"))
        say(*_fires(lambda: spec_check("w", e_ok, True),
                    "a TS input that is not AMENDMENT 4's", "no OptTS on the keyword line"))
        say(*_fires(lambda: spec_check("w", ["! B3LYP OptTS NumFreq QMMM", "%geom", "  Calc_Hess true",
                                             "  Hybrid_Hess {0 1 2} end", "end"], True),
                    "Hybrid_Hess that is not the 24 substrate atoms", "Hybrid_Hess is not the 24"))
        say(*_fires(lambda: spec_check("w", ["! B3LYP OptTS Freq QMMM", "%geom", "  Calc_Hess true",
                                             "  Hybrid_Hess {%s} end" % " ".join(map(str, range(NQM))),
                                             "  coordsys cartesian", "end"], True),
                    "a TS input still carrying Freq and Cartesians", "coordsys cartesian is present"))
        prms = tmp / "short.prms"
        prms.write_text("\n".join("a b 0.0 d 0.0 0.0 0.0" for _ in range(NQM)) + "\n")
        say(*_fires(lambda: spheres_of(prms), "a force field with no sphere rows", "parameter rows"))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print("\n%s" % ("ALL PASS" if ok else "SOMETHING FAILED"))
    return 0 if ok else 1


def _isnum(v):
    try:
        float(v)
        return True
    except (TypeError, ValueError):
        return False


def main():
    if len(sys.argv) == 3 and sys.argv[1] == "selftest":
        return selftest(sys.argv[2])
    if len(sys.argv) != 5 or sys.argv[1] != "relaxed":
        sys.exit(__doc__)
    return relaxed(sys.argv[2], sys.argv[3], sys.argv[4])


if __name__ == "__main__":
    sys.exit(main())
