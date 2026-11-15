#!/usr/bin/env python3
"""cs6_d39.py  DIAGNOSTIC 39: read the detachment energy from six outputs.

THE CONVENTION IS THE COMMITTED ONE AND IS NOT REDEFINED HERE.
STAGE6_PROBE2_20261027.txt states it and both probe results repeat it:

    VDE = E(monoanion, charge -1, multiplicity 2) - E(dianion, charge -2, multiplicity 1),
    both at the same geometry.  Positive means bound.

D39 uses it unchanged and in the same direction. The constant is the one
vde_probe/cs11_vde.py already uses.

THE THREE DIANION ENERGIES ARE NOT IN THIS FILE. They are read from the committed
d38/N, d38/B1 and d38/A1 outputs, whose sha256 are checked first. A diagnostic that
retyped them as constants could not notice if the file they came from had changed.

THERE IS NO BARRIER HERE AND NO MEANS TO COMPUTE ONE. Six single points at one
geometry. The only energy conversion in this file is the one to electron volts, which
the committed convention requires; G39.12 fails the bundle if any other appears.

D39 PRE-REGISTERED READINGS AND NO VERDICT. This reports the readings.

Usage: cs6_d39.py <results-repo> <out-dir>
"""
import hashlib, io, os, re, sys

S6 = "phase2.2/cs_stage6"

# vde_probe/cs11_vde.py
EV = 27.211386245988

# The three committed dianion outputs D39 reuses instead of re-running, with the sha256
# STAGE6_DIAGNOSTIC39_20261115.txt section 6 records for each.
DIANION_SHA = {
    "N":  "18a58ac48f7890b6bbd7622a507cd560a75ce8f13cb3742fa6cec9ce6ef28a9c",
    "B1": "307eab716f3f719ebff3afc9d07b58dbfc891abb7ace04599b6fcb04e6835929",
    "A1": "64013339b8e8f5ee8442638f785cebd32434512cec7a4551335ce796d8228422",
}

# STAGE6_DIAGNOSTIC38_RESULT_20261115.txt, R38.5 and R38.7. The Hartree-Fock highest
# occupied orbital at def2-SVPD, which is the proxy G38.10 decided interpretability on.
D38_HOMO = {"N": -2.0444, "B1": -4.8267, "A1": +0.3539}

# The committed DFT comparators, each with the document it comes from and its basis.
# None means no committed run measures that arrangement.
CMP = {
    "N":  (0.5013, "STAGE6_PROBE1_ERRATUM1_20261027.txt, C5", "def2-SVPD",
           "same geometry AND same basis; the method is the only difference"),
    "B1": (3.4366, "vde_designs/PROBE3_RESULT.txt, S6-P (= vde_probe dref_p1)",
           "ma-def2-SVP",
           "same geometry, same site, same charge; the basis differs too"),
    "A1": (None, None, None,
           "no committed run measures -1 e at site 678 alone"),
}

# Context for the bare case only, from the same Probe 1 erratum: the same quantity at
# two other bases, which is the sign flip this project already measured.
P1_CONTEXT = ((-0.3167, "C1", "def2-SVP"), (0.5123, "C3", "ma-def2-SVP"))

# Probe 3's S6-2 carries BOTH spheres, so it is not A1's comparator. Carried so the
# result can say what it is instead of leaving the reader to assume.
S6_2 = (1.1058, "vde_designs/PROBE3_RESULT.txt, S6-2", "ma-def2-SVP")

# STAGE6_DIAGNOSTIC39_AMENDMENT1_20261115.txt. A convention adopted before the outputs
# were read, not a literature value, and the result says so.
SINGLES_MIN = 90.0

MROUTE = ("M_N", "M_B1", "M_A1")
XROUTE = ("X_N", "X_B1", "X_A1")
SYS = (("N", "no embedding charge"),
       ("B1", "+1.000000 e at site 1690"),
       ("A1", "-1.000000 e at site 678"))
NEL_M, NEL_X = 117, 118
MULT_M, MULT_X = 2, 1

gates = []
OUT = []


def say(s=""):
    OUT.append(s)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def rd(p):
    return io.open(p, encoding="utf-8", errors="replace").read()


def num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def last(pat, txt, grp=1):
    # re.M, because every energy pattern here is anchored to the start of a line and
    # without it ^ only matches the start of the whole file. Every energy read as
    # absent until this was added.
    m = re.findall(pat, txt, re.M)
    if not m:
        return None
    v = m[-1]
    if isinstance(v, tuple):
        v = v[grp - 1]
    return num(v)


def read_scf(txt):
    """The charge, multiplicity, electron count and SCF cycle count ORCA reports."""
    d = {}
    d["charge"] = last(r"Total Charge\s+Charge\s+\.\.\.\.\s+(-?\d+)", txt)
    d["mult"] = last(r"Multiplicity\s+Mult\s+\.\.\.\.\s+(\d+)", txt)
    d["nel"] = last(r"Number of Electrons\s+NEL\s+\.\.\.\.\s+(\d+)", txt)
    d["cycles"] = last(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES", txt)
    d["converged"] = d["cycles"] is not None
    d["normal"] = "ORCA TERMINATED NORMALLY" in txt
    m = re.search(r"TOTAL RUN TIME:\s*(.+)", txt)
    d["walltime"] = m.group(1).strip() if m else None
    return d


def read_mdci(txt):
    """The correlated energies, as ORCA labels them."""
    return {
        "e0": last(r"^E\(0\)\s+\.\.\.\s+(-?\d+\.\d+)", txt),
        "ecorr": last(r"^E\(CORR\)\(corrected\)\s+\.\.\.\s+(-?\d+\.\d+)", txt),
        "eccsd": last(r"^E\(CCSD\)\s+\.\.\.\s+(-?\d+\.\d+)", txt),
        "eccsdt": last(r"^E\(CCSD\(T\)\)\s+\.\.\.\s+(-?\d+\.\d+)", txt),
        "t1": last(r"^T1 diagnostic\s+\.\.\.\s+(-?\d+\.\d+)", txt),
    }


def read_s2(txt):
    """The expectation value of S squared, the ideal value and ORCA's own deviation.

    ORCA writes the ideal line's spin as a fraction in some builds and a decimal in
    others ("for S=1/2" against "for S=0.5"), so nothing is matched after S* (S+1)
    except the colon. A closed-shell output prints none of these, so None here is not a
    failure for route X.
    """
    return {
        "s2": last(r"Expectation value of\s*<S\*\*2>\s*:\s*(-?\d+\.\d+)", txt),
        "ideal": last(r"Ideal value\s*S\*\(S\+1\)[^:\n]*:\s*(-?\d+\.\d+)", txt),
        "orca_dev": last(r"^\s*Deviation\s*:\s*(-?\d+\.\d+)", txt),
    }


def read_eom(txt):
    """The EOM roots, in the form the ORCA 6.0 manual documents:

        IROOT=  1:  0.322826 au     8.785 eV   70852.1 cm**-1

    followed, for each root, by an amplitude table ending in

        Percentage singles character=     96.89

    Returns a list of (iroot, au, ev, singles) in the order printed. The singles value
    is matched to the root it follows, so a missing one stays None rather than sliding
    onto the wrong root.
    """
    roots = []
    pat = re.compile(r"IROOT=\s*(\d+):\s*(-?\d+\.\d+)\s*au\s+(-?\d+\.\d+)\s*eV")
    sing = re.compile(r"Percentage singles character=\s*(-?\d+\.\d+)")
    hits = list(pat.finditer(txt))
    for i, m in enumerate(hits):
        end = hits[i + 1].start() if i + 1 < len(hits) else len(txt)
        s = sing.search(txt, m.end(), end)
        roots.append((int(m.group(1)), num(m.group(2)), num(m.group(3)),
                      num(s.group(1)) if s else None))
    return roots


def fmt(v, w=16, p=9, plus=False):
    if v is None:
        return " " * (w - 1) + "-"
    return ("%+*.*f" if plus else "%*.*f") % (w, p, v)


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("usage: cs6_d39.py <results-repo> <out-dir>\n")
        return 2
    res, out = sys.argv[1:3]
    s6 = os.path.join(res, S6)
    d38 = os.path.join(s6, "d38")
    d39 = os.path.join(s6, "d39")

    # ---- the dianion halves, read from the committed files after a hash check -------
    di = {}
    for k, want in sorted(DIANION_SHA.items()):
        p = os.path.join(d38, k, "job.out")
        if not os.path.isfile(p):
            sys.stderr.write("REFUSED: %s is missing.\n" % p)
            return 1
        got = sha(p)
        if got != want:
            sys.stderr.write(
                "REFUSED: d38/%s/job.out has sha256 %s, the pre-registration records "
                "%s.\n   The dianion half of a subtraction has changed. Nothing is "
                "read.\n" % (k, got, want))
            return 1
        t = rd(p)
        di[k] = read_mdci(t)
        di[k].update(read_scf(t))
        di[k]["sha"] = got

    # ---- the six D39 outputs -------------------------------------------------------
    M, X = {}, {}
    missing = []
    for r in MROUTE:
        p = os.path.join(d39, r, "job.out")
        if not os.path.isfile(p):
            missing.append(r)
            continue
        t = rd(p)
        d = read_scf(t)
        d.update(read_mdci(t))
        d.update(read_s2(t))
        d["sha"] = sha(p)
        M[r] = d
    for r in XROUTE:
        p = os.path.join(d39, r, "job.out")
        if not os.path.isfile(p):
            missing.append(r)
            continue
        t = rd(p)
        d = read_scf(t)
        d["roots"] = read_eom(t)
        d["sha"] = sha(p)
        X[r] = d

    say("=" * 78)
    say("STAGE 6, DIAGNOSTIC 39.  RESULT")
    say("THE DETACHMENT ENERGY AT THE BENCHMARK LEVEL")
    say("=" * 78)
    say()
    say("Pre-registered in STAGE6_DIAGNOSTIC39_20261115.txt and amended once, in")
    say("STAGE6_DIAGNOSTIC39_AMENDMENT1_20261115.txt, before any of these outputs was")
    say("read. Generated by cs6_d39.py from the six outputs and from Diagnostic 38's")
    say("three committed dianion outputs.")
    say()
    say("D39 PRE-REGISTERED READINGS AND NO VERDICT. This document reports the readings.")
    say("Section 10 of the pre-registration lists what it does not do, and that list is")
    say("restated at the end.")
    say()
    say("THE CONVENTION, which is the committed one and is not redefined here:")
    say("   VDE = E(monoanion, charge -1, multiplicity 2)")
    say("       - E(dianion,   charge -2, multiplicity 1), same geometry.")
    say("   Positive means bound. STAGE6_PROBE2_20261027.txt states it; both probe")
    say("   results repeat it. Constant EV = %.12f, from vde_probe/cs11_vde.py." % EV)
    say()

    if missing:
        say("=" * 78)
        say("OUTPUTS NOT PRESENT")
        say("=" * 78)
        say()
        say("   %s" % ", ".join(missing))
        say()
        for r in missing:
            if r in XROUTE and "X_N" in X and not X["X_N"].get("normal"):
                say("   %s  C39.2: X_N did not terminate normally, so the remaining route"
                    % r)
                say("         X runs were not attempted. This is a recorded outcome.")
            else:
                gates.append("%s produced no output and C39.2 does not explain it" % r)
        say()

    # ---- R39.1  route M ------------------------------------------------------------
    say("=" * 78)
    say("R39.1  ROUTE M, THE THREE MONOANIONS")
    say("=" * 78)
    say()
    say("   run     chg mult  NEL  cyc              E(0)       E(CORR)corr"
        "           E(CCSD)        E(CCSD(T))           T1")
    for r in MROUTE:
        d = M.get(r)
        if not d:
            say("   %-7s %s" % (r, "no output"))
            continue
        say("   %-7s %3s %4s %4s %4s %s %s %s %s %s"
            % (r,
               "%d" % d["charge"] if d["charge"] is not None else "-",
               "%d" % d["mult"] if d["mult"] is not None else "-",
               "%d" % d["nel"] if d["nel"] is not None else "-",
               "%d" % d["cycles"] if d["cycles"] is not None else "-",
               fmt(d["e0"]), fmt(d["ecorr"], 17), fmt(d["eccsd"], 17),
               fmt(d["eccsdt"], 17), fmt(d["t1"], 12, 9)))
        if d["charge"] != -1:
            gates.append("R39.1 %s reports charge %s, expected -1" % (r, d["charge"]))
        if d["mult"] != MULT_M:
            gates.append("R39.1 %s reports multiplicity %s, expected %d"
                         % (r, d["mult"], MULT_M))
        if d["nel"] != NEL_M:
            gates.append("R39.1 %s reports %s electrons, expected %d"
                         % (r, d["nel"], NEL_M))
        if not d["converged"]:
            say("           C39.3: SCF did not converge. No energy is read from this run.")
        if not d["normal"]:
            say("           ORCA did not terminate normally.")
    say()
    say("   The charge, multiplicity and electron count above are ORCA's own report of")
    say("   what it ran, not a restatement of the input. 117 electrons against the")
    say("   dianion's 118 is the detachment this diagnostic is about.")
    say()
    for r in MROUTE:
        d = M.get(r)
        if d and d.get("walltime"):
            say("   %-7s %s" % (r, d["walltime"]))
    say()

    # ---- C39.1 ---------------------------------------------------------------------
    no_t = [r for r in MROUTE if M.get(r) and M[r]["eccsdt"] is None
            and M[r]["eccsd"] is not None]
    any_ccsd = [r for r in MROUTE if M.get(r) and M[r]["eccsd"] is not None]
    say("=" * 78)
    say("C39.1  WAS THE OPEN-SHELL TRIPLES CORRECTION AVAILABLE")
    say("=" * 78)
    say()
    if not M:
        say("   No route M output to read.")
    elif not any_ccsd:
        say("   NOT ANSWERABLE. No route M run reported E(CCSD) either, so the absence")
        say("   of E(CCSD(T)) says nothing about whether open-shell triples exist.")
        gates.append("C39.1 no route M run reported E(CCSD); the outputs were not read")
    elif no_t:
        say("   NO. These runs reported E(CCSD) but no E(CCSD(T)):")
        say("      %s" % ", ".join(no_t))
        say()
        say("   C39.1 fired. The detachment energy is reported at the CCSD level from")
        say("   R39.4 and R39.3 is not produced for those systems. The pre-registration's")
        say("   section 4 took plain open-shell (T) in 6.0.1 as available by inference")
        say("   from the manual's 'Since ORCA 4.2' statement about the iterative variant.")
        say("   THAT INFERENCE WAS WRONG and this is it going on the record as wrong.")
    else:
        say("   YES. Every route M run that converged reported E(CCSD(T)).")
        say("   The pre-registration's section 4 inference holds: plain open-shell (T)")
        say("   is available in 6.0.1. That was an inference and it is now a measurement.")
    say()

    # ---- R39.5  the open-shell reference -------------------------------------------
    say("=" * 78)
    say("R39.5  THE OPEN-SHELL REFERENCE, MEASURED RATHER THAN ASSUMED")
    say("=" * 78)
    say()
    say("   The manual's surviving warning is that 'an UHF reference wavefunction")
    say("   suffers from spin-contamination which can spoil the results'. It is general")
    say("   and not specific to DLPNO. The pre-registration records that I had wrongly")
    say("   attributed a second, stronger warning to this method; that one belongs to")
    say("   the older UHF-LPNO scheme. This reading is why the concern does not have to")
    say("   be argued about.")
    say()
    say("   run        <S**2>      ideal   deviation  ORCA's own")
    worst = 0.0
    for r in MROUTE:
        d = M.get(r)
        if not d:
            continue
        s2, idl, od = d.get("s2"), d.get("ideal"), d.get("orca_dev")
        if s2 is None:
            say("   %-7s %s" % (r, "not printed"))
            gates.append("R39.5 %s is an open-shell run with no <S**2> printed" % r)
            continue
        assumed = False
        if idl is None:
            idl, assumed = 0.75, True
        dev = s2 - idl
        worst = max(worst, abs(dev))
        say("   %-7s %10.6f %10.6f  %+.6f  %s%s"
            % (r, s2, idl, dev,
               "%+.6f" % od if od is not None else "        -",
               "   ideal not printed; 0.75 assumed for a doublet" if assumed else ""))
        # ORCA prints the deviation itself. If both are present they must agree, or one
        # of the two numbers is not what this script thinks it is.
        if od is not None and abs(abs(od) - abs(dev)) > 5e-6:
            gates.append("R39.5 %s: ORCA's deviation %+.6f and the difference computed "
                         "here %+.6f do not agree" % (r, od, dev))
    say()
    if M and worst:
        say("   Worst deviation from the ideal doublet value: %+.6f" % worst)
        if worst < 0.01:
            say("   That is negligible. The UHF reference is essentially uncontaminated")
            say("   at this geometry, so route X's cross-check below is a confirmation")
            say("   rather than a correction. This is a measurement, not a reassurance.")
        elif worst < 0.1:
            say("   Small but not zero. Route X's cross-check below is the thing that")
            say("   says whether it mattered for the detachment energy.")
        else:
            say("   NOT NEGLIGIBLE. The UHF reference is contaminated and route M's")
            say("   numbers carry that. R39.6 is the measurement of what it cost.")
    say()

    # ---- R39.2 and R39.10  route X -------------------------------------------------
    say("=" * 78)
    say("R39.2  ROUTE X, THE EOM ROOTS")
    say("R39.10 THE SINGLES CHARACTER OF EACH ROOT  (AMENDMENT 1)")
    say("=" * 78)
    say()
    say("   The lowest root of an IP-EOM calculation is the detachment energy directly:")
    say("   E(N-1) - E(N) from the closed-shell dianion, with no open-shell reference")
    say("   anywhere. Same definition as the convention above, so the numbers are")
    say("   directly comparable with no sign change.")
    say()
    xlow = {}
    for r in XROUTE:
        d = X.get(r)
        if not d:
            say("   %-7s no output" % r)
            continue
        if d["charge"] != -2:
            gates.append("R39.2 %s reports charge %s, expected -2" % (r, d["charge"]))
        if d["nel"] != NEL_X:
            gates.append("R39.2 %s reports %s electrons, expected %d"
                         % (r, d["nel"], NEL_X))
        if not d["roots"]:
            say("   %-7s terminated %s, but no IROOT line was found."
                % (r, "normally" if d["normal"] else "abnormally"))
            say("           Route X produced no readable root for this system. Recorded.")
            # C39.2 requires the refusal to be recorded with ORCA's own message, so it
            # is quoted here rather than described.
            t = rd(os.path.join(d39, r, "job.out"))
            msgs = [l.strip() for l in t.split("\n")
                    if re.search(r"(?i)\b(error|aborting|abort|not available|"
                                 r"no auxiliary|unknown|error termination)\b", l)]
            if msgs:
                say("           ORCA's own words:")
                for l in msgs[-6:]:
                    say("              %s" % l[:72])
            else:
                say("           ORCA printed no error line this script recognises.")
            continue
        say("   %-7s charge %d  NEL %d  SCF %s cycles"
            % (r, d["charge"], d["nel"],
               "%d" % d["cycles"] if d["cycles"] is not None else "-"))
        for (ir, au, ev, sg) in d["roots"]:
            say("           IROOT=%2d  %+.6f au  %+8.3f eV   singles %s"
                % (ir, au, ev, "%6.2f %%" % sg if sg is not None else "     -"))
        lo = min(d["roots"], key=lambda x: x[2])
        xlow[r] = lo
        say("           lowest root: IROOT=%d at %+.4f eV" % (lo[0], lo[2]))
        if lo[3] is None:
            say("           singles character not printed for it. Recorded.")
            gates.append("G39.18 %s has a lowest root with no singles character" % r)
        elif lo[3] < SINGLES_MIN:
            say("           NOT A CLEAN ONE-ELECTRON DETACHMENT: singles %.2f %% is"
                % lo[3])
            say("           below the %.0f %% bar Amendment 1 set before these outputs"
                % SINGLES_MIN)
            say("           were read. The number stands but is not read as evidence")
            say("           about route M's reference.")
        else:
            say("           clean one-electron detachment: singles %.2f %% is at or"
                % lo[3])
            say("           above the %.0f %% bar." % SINGLES_MIN)
        say()
    say("   The %.0f %% bar is a convention adopted in Amendment 1 before these outputs"
        % SINGLES_MIN)
    say("   were read. It is NOT a literature value. No committed document in this")
    say("   project sets a threshold for this quantity and I could not find one. The")
    say("   manual's worked example reports 96.89 % for a normal ionisation, which is")
    say("   context and not a source.")
    say()

    # ---- the reused halves ---------------------------------------------------------
    say("=" * 78)
    say("THE DIANION HALVES, REUSED AND NOT RE-RUN")
    say("=" * 78)
    say()
    say("   Read from Diagnostic 38's committed outputs after checking each sha256")
    say("   against the pre-registration's section 6. Not retyped as constants.")
    say()
    say("   d38 run             E(CCSD)        E(CCSD(T))           T1  sha256")
    for k, _ in SYS:
        d = di[k]
        say("   %-7s %s %s %s  %s"
            % (k, fmt(d["eccsd"], 17), fmt(d["eccsdt"], 17), fmt(d["t1"], 12, 9),
               d["sha"][:16]))
    say()

    # ---- R39.3 and R39.4  the detachment energies ----------------------------------
    rows = []
    VT, VC = {}, {}       # system -> VDE at CCSD(T) and at CCSD, in eV
    say("=" * 78)
    say("R39.3  THE DETACHMENT ENERGY AT CCSD(T)")
    say("R39.4  THE DETACHMENT ENERGY AT CCSD")
    say("=" * 78)
    say()
    say("   system  embedding charge              VDE CCSD(T)       VDE CCSD")
    for k, desc in SYS:
        r = "M_" + k
        d = M.get(r)
        dd = di[k]
        vt = vc = None
        if d and d["converged"]:
            if d["eccsdt"] is not None and dd["eccsdt"] is not None:
                vt = (d["eccsdt"] - dd["eccsdt"]) * EV
            if d["eccsd"] is not None and dd["eccsd"] is not None:
                vc = (d["eccsd"] - dd["eccsd"]) * EV
        say("   %-7s %-28s %s %s"
            % (k, desc,
               ("%+12.4f eV" % vt) if vt is not None else "           -",
               ("%+12.4f eV" % vc) if vc is not None else "           -"))
        rows.append((k, desc, vt, vc))
        VT[k], VC[k] = vt, vc
    say()
    say("   Each is one subtraction: this diagnostic's monoanion minus Diagnostic 38's")
    say("   dianion, same geometry, same basis, same auxiliary treatment, same SCF")
    say("   settings. Both halves measured; nothing here is estimated.")
    say()

    # ---- R39.6  the cross-check ----------------------------------------------------
    say("=" * 78)
    say("R39.6  THE TWO CCSD NUMBERS AGAINST EACH OTHER")
    say("=" * 78)
    say()
    say("   Route M's CCSD detachment energy goes through a UHF reference. Route X's")
    say("   lowest root does not. Same physical quantity, two formalisms.")
    say()
    say("   system     route M CCSD      route X root 1     difference   read as")
    any_x = False
    for k, _ in SYS:
        vc = VC.get(k)
        lo = xlow.get("X_" + k)
        if vc is None or lo is None:
            say("   %-7s %s" % (k, "not comparable: one side is missing"))
            continue
        any_x = True
        diff = vc - lo[2]
        if lo[3] is None or lo[3] < SINGLES_MIN:
            how = "NOT A CLEAN ONE-ELECTRON DETACHMENT"
        else:
            how = "evidence about the UHF reference"
        say("   %-7s %+12.4f eV  %+12.4f eV  %+11.4f eV   %s"
            % (k, vc, lo[2], diff, how))
    say()
    if any_x:
        say("   Where the label is 'evidence about the UHF reference': a small")
        say("   difference means spin contamination is not affecting this quantity at")
        say("   this geometry, and route M's CCSD(T) number carries that support. A")
        say("   large one is the measurement of the problem. No preference is asserted")
        say("   between the two and route X is not substituted for route M (D39.5).")
    else:
        say("   Route X produced nothing comparable, so route M's open-shell reference")
        say("   is uncontrolled and R39.5's deviation is all there is on it.")
    say()

    # ---- R39.7  against the committed DFT comparators ------------------------------
    say("=" * 78)
    say("R39.7  AGAINST THE COMMITTED DFT COMPARATORS")
    say("=" * 78)
    say()
    say("   Every comparator below is B3LYP-D3BJ, from a committed document, at the")
    say("   identical geometry. A method change is expected to move a detachment")
    say("   energy; the size of the move is the reading. None of this is a gate.")
    say()
    for k, _ in SYS:
        val, srcf, basis, note = CMP[k]
        vt = VT.get(k)
        say("   %s" % k)
        if val is None:
            say("      no committed comparator. %s." % note)
            say("      Probe 3's S6-2 reports %+.4f eV at %s but its charge file holds"
                % (S6_2[0], S6_2[2]))
            say("      two charges, the positive sphere as well, so it is the pair and")
            say("      not this arrangement. It is not used as a comparator here.")
            if vt is not None:
                say("      D39 measures %+.4f eV. It is the first of its kind in this"
                    % vt)
                say("      record and has nothing to reproduce.")
        else:
            say("      committed  %+.4f eV   %s" % (val, basis))
            say("                 %s" % srcf)
            say("                 %s" % note)
            if vt is not None:
                say("      D39        %+.4f eV   def2-SVPD, CCSD(T)" % vt)
                say("      difference %+.4f eV" % (vt - val))
            else:
                say("      D39        no CCSD(T) value for this system")
        say()
    say("   For the bare case only, the same quantity at two other bases, from the same")
    say("   Probe 1 erratum. This is the sign flip this project has already measured on")
    say("   this very quantity, and it is why a single-basis result is a real limit:")
    for v, nm, b in P1_CONTEXT:
        say("      %-3s %-14s %+.4f eV" % (nm, b, v))
    say()

    # ---- R39.8  bound or metastable, and the proxy ---------------------------------
    say("=" * 78)
    say("R39.8  BOUND OR METASTABLE, AND WHAT THE PROXY SAID")
    say("=" * 78)
    say()
    say("   The classification uses the committed sign rule: positive means bound.")
    say("   The proxy column is -HOMO from Diagnostic 38's R38.5 and R38.7, which is")
    say("   the Koopmans-style quantity G38.10 decided interpretability on. The")
    say("   comparison is what R39.8 pre-registered.")
    say()
    say("   system   VDE CCSD(T)   class        -HOMO (proxy)   proxy class   agree")
    for k, _ in SYS:
        vt = VT.get(k)
        pr = -D38_HOMO[k]
        pc = "BOUND" if pr > 0 else "METASTABLE"
        if vt is None:
            say("   %-7s %s" % (k, "no CCSD(T) value"))
            continue
        cl = "BOUND" if vt > 0 else "METASTABLE"
        say("   %-7s %+11.4f eV  %-11s %+11.4f eV   %-11s   %s"
            % (k, vt, cl, pr, pc, "yes" if cl == pc else "NO"))
    say()
    offs = []
    for k, _ in SYS:
        vt = VT.get(k)
        if vt is not None:
            offs.append((k, vt - (-D38_HOMO[k])))
    if offs:
        say("   proxy error, VDE minus -HOMO:")
        for k, o in offs:
            say("      %-3s %+.4f eV" % (k, o))
        say()
        say("   This is the quantitative form of the agreement R39.8 asked for. It says")
        say("   how good the proxy G38.10 relied on actually is on this quantity, at")
        say("   this geometry and this basis, and nothing more than that.")
    say()
    a1 = VT.get("A1")
    if a1 is not None:
        say("   C39.4  A1")
        if a1 <= 0:
            say("      A1's detachment energy is negative, which is what the")
            say("      pre-registration expected. D38 measured its dianion's highest")
            say("      occupied orbital at %+.4f eV and G38.10 called it NOT"
                % D38_HOMO["A1"])
            say("      INTERPRETABLE on that ground. The observable agrees with the")
            say("      proxy.")
            say()
            say("      This number is a bound-state method's description of what the")
            say("      record's citation of Dreuw and Cederbaum calls a resonance, one")
            say("      that 'can in principle not be correctly computed with bound state")
            say("      methods'. Basis-set dependence of it is expected rather than")
            say("      anomalous. IT IS NOT PLACED BESIDE THE BOUND VALUES WITHOUT THIS")
            say("      LABEL, and it does not retract G38.10.")
        else:
            say("      A1's detachment energy is POSITIVE, which the pre-registration")
            say("      did not expect. The observable DISAGREES with the proxy: D38")
            say("      measured the orbital at %+.4f eV and G38.10 excluded the arm on"
                % D38_HOMO["A1"])
            say("      that basis.")
            say()
            say("      C39.4 states what follows and what does not. This is a finding")
            say("      about G38.10's rule: the orbital energy was the wrong test for")
            say("      this arrangement. It is NOT a licence to reinterpret Diagnostic")
            say("      38's numbers here, and nothing in this document does so.")
        say()

    # ---- R39.9  scratch ------------------------------------------------------------
    say("=" * 78)
    say("R39.9  SCRATCH, WITH THE PEAK AND THE RESIDUE KEPT APART")
    say("=" * 78)
    say()
    sp = os.path.join(d39, "run.pbs.out")
    if os.path.isfile(sp):
        t = rd(sp)
        shown = False
        for pat in (r"R39\.9[^\n]*", r"host=[^\n]*"):
            for m in re.findall(pat, t):
                say("   %s" % m.strip())
                shown = True
        if not shown:
            say("   the scheduler output is present but carries no R39.9 line")
    else:
        say("   no scheduler output placed, so nothing to report")
    say()
    say("   Diagnostic 38's result reported its footprint as 11 MB and said that figure")
    say("   was what a later job should use. It was the residue, measured after ORCA had")
    say("   deleted its temporaries. STAGE6_DIAGNOSTIC38_RESULT_ERRATUM1_20261115.txt")
    say("   withdraws that and states the corrected bound, 31 GB to 43 GB per")
    say("   closed-shell run. D39's sampler runs while ORCA runs and the two quantities")
    say("   are reported under separate names above so they cannot be conflated again.")
    say()

    # ---- limits --------------------------------------------------------------------
    say("=" * 78)
    say("WHAT THIS DOES NOT DO")
    say("=" * 78)
    say()
    say("   NO BARRIER (D39.1). Six single points at one geometry. This script carries")
    say("   no conversion out of hartree but the one to electron volts that the")
    say("   committed convention requires. G39.12 fails the bundle if another appears.")
    say()
    say("   IT REVISITS NO D34 VERDICT (D39.2). ADD-s1690_qp10000 and ADD-s678_qm10000")
    say("   keep whatever reliability classification D34 gave them.")
    say()
    say("   NO DESIGN CONCLUSION (D39.3). Nothing here says either arrangement is a good")
    say("   design or one worth pursuing. They were chosen because D38 already ran their")
    say("   dianions.")
    say()
    say("   IT DOES NOT LICENSE THE ddE TEST (D39.4). D38 licensed that. D39 neither")
    say("   extends nor withdraws it.")
    say()
    say("   ROUTE X CANNOT PROMOTE ITSELF (D39.5). It is a control on route M's")
    say("   open-shell reference. Route X's number is not substituted for route M's and")
    say("   route M's is not adjusted toward route X's.")
    say()
    say("   ONE GEOMETRY, AND IT IS THE REACTANT (D39.6). Nothing here says whether")
    say("   boundedness changes along the reaction coordinate.")
    say()
    say("   ONE BASIS (D39.7). def2-SVPD only, on a quantity Probe 1 already showed")
    say("   flipping sign between def2-SVP and def2-SVPD.")
    say()
    say("   IT DOES NOT COMPARE TO EXPERIMENT (D39.8). E7.11 forbids it and nothing here")
    say("   is a free energy or a barrier in any case.")
    say()
    say("   THE PARSER WAS TESTED BEFORE IT READ THESE OUTPUTS. Neither an open-shell")
    say("   nor an EOM output had ever been read by this project, so this script was")
    say("   exercised against nine synthetic ORCA outputs covering every branch above:")
    say("   the clean case, open-shell triples absent, route X refusing its auxiliary")
    say("   basis, an SCF that does not converge, A1 coming out bound, a lowest root")
    say("   below the singles bar, no roots at all, no spin-contamination block, and a")
    say("   badly contaminated reference. Three real bugs were found that way, the worst")
    say("   of which made every energy read as absent while the document still reported")
    say("   a conclusion. The fixtures are NOT committed: nothing fabricated belongs in")
    say("   this record, and they are in the bundle instead.")
    say()
    say("   NOTHING WAS SHAPED (D39.9). Method, basis, auxiliary treatment, SCF settings")
    say("   and geometry are D38's. The two arrangements are D38's. The convention and")
    say("   the constant are the committed probes'. The comparators are quoted from")
    say("   committed documents with their bases named. The 90 per cent singles bar was")
    say("   set in Amendment 1 before these outputs were read.")
    say()

    say("=" * 78)
    if gates:
        say("GATE FAILURES: %d" % len(gates))
        for g in gates:
            say("   %s" % g)
        say()
        say("Read this document before committing anything.")
    else:
        say("ALL GATES PASS")
    say("=" * 78)

    if not os.path.isdir(out):
        os.makedirs(out)
    doc = "\n".join(OUT) + "\n"
    io.open(os.path.join(out, "STAGE6_DIAGNOSTIC39_RESULT_20261115.txt"), "w",
            encoding="utf-8").write(doc)

    with io.open(os.path.join(out, "d39.tsv"), "w", encoding="utf-8") as f:
        f.write("system\tarrangement\tVDE_CCSDT_eV\tVDE_CCSD_eV\tVDE_EOM_eV\t"
                "EOM_singles_pct\tproxy_minus_HOMO_eV\tDFT_comparator_eV\t"
                "DFT_basis\tclass\n")
        for k, desc in SYS:
            vt = VT.get(k)
            vc = VC.get(k)
            lo = xlow.get("X_" + k)
            val, _srcf, basis, _n = CMP[k]
            f.write("%s\t%s\t%s\t%s\t%s\t%s\t%+.4f\t%s\t%s\t%s\n"
                    % (k, desc,
                       "%+.4f" % vt if vt is not None else "",
                       "%+.4f" % vc if vc is not None else "",
                       "%+.4f" % lo[2] if lo else "",
                       "%.2f" % lo[3] if lo and lo[3] is not None else "",
                       -D38_HOMO[k],
                       "%+.4f" % val if val is not None else "",
                       basis or "",
                       ("BOUND" if vt > 0 else "METASTABLE") if vt is not None else ""))

    sys.stdout.write(doc)
    return 1 if gates else 0


if __name__ == "__main__":
    sys.exit(main())
