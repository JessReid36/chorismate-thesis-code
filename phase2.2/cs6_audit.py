#!/usr/bin/env python3
"""Stage 6 DIAGNOSTIC 6: a wavefunction-health audit of every committed ORCA output on the chorismate dianion.

NO NEW CALCULATION AND NO NETWORK. It opens every committed *.tar.xz in the results repo, finds every ORCA
output inside, keeps the ones that are this dianion (118 electrons, 59 occupied orbitals, a final single point
energy present), and reads four fields that no earlier diagnostic read: the SCF cycle count, the HOMO-LUMO gap,
the dipole magnitude and the largest Loewdin charge on a light atom. Diagnostics 1 to 5 read energies only.

Why: in the committed frozen runs, at def2-SVP with no diffuse functions anywhere, the SCF took 12 to 16 cycles
for every arrangement up to three charges, 21 for five, 36 to 40 for ten and 42 to 47 for twenty, while the gap
fell 3.40, 3.19, 2.29, 2.13, 1.76, 0.85, 0.61 eV. The breakdown was recorded in the original outputs from the
first design run onward and was never read.

The rule is STAGE6_DIAGNOSTIC6_20261026.txt D6.3. Every threshold is anchored to the charge-free control runs
of the SAME basis-set size in the same corpus, never to an absolute value chosen freely, and the control
distributions are printed so the separation is checkable by eye.

usage: cs6_audit.py <results-repo> <out-dir>
"""
import sys, re, tarfile, statistics as st
from pathlib import Path

EV = 27.211386245988
NEL, NOCC = 118, 59
D_LW = 0.25      # e      Loewdin rise over the same-basis control that counts as the density having MOVED
R_DIP = 2.0      #        dipole over that control that counts as the same
GAP_ABS = 1.0    # eV     near-degenerate in absolute terms
R_GAP = 0.5      #        and also less than half its own control's gap, so a basis whose controls already sit
                 #        low is not condemned for sitting low
R_CYC = 2.0      #        SCF cycles over control that count as convergence difficulty
MM_NPC = 200     #        above this the point charges are an MM environment in a QM/MM run, not a design
METALS = {"Na", "K", "Mg", "Ca", "Li", "Zn", "Fe", "Cl", "Br", "I"}


def stop(m):
    print("STOP: %s" % m, file=sys.stderr)
    raise SystemExit(1)


def read(s):
    """Only ORCA's own printed fields. Returns None unless this is a completed single point on THIS dianion."""
    if "FINAL SINGLE POINT ENERGY" not in s:
        return None
    m = re.search(r"Number of Electrons\s+NEL\s+\.+\s+(\d+)", s)
    if not m or int(m.group(1)) != NEL:
        return None
    nb = re.search(r"Number of basis functions\s+\.+\s+(\d+)", s)
    pc = re.search(r"Reading point charge file\s*\.*\s*ok \((\d+) point charges\)", s)
    cy = re.findall(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES", s)
    dp = re.findall(r"Magnitude \(Debye\)\s+:\s+(\d+\.\d+)", s)
    en = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", s)
    k = s.rfind("ORBITAL ENERGIES")
    occ, vir = [], []
    if k >= 0:
        for l in s[k:].splitlines()[4:]:
            f = l.split()
            if len(f) != 4:
                if occ and vir:
                    break
                continue
            try:
                o, eh = float(f[1]), float(f[2])
            except ValueError:
                break
            (occ if o > 0.5 else vir).append(eh)
    if len(occ) != NOCC:
        return None
    k = s.rfind("LOEWDIN ATOMIC CHARGES")
    q = []
    if k >= 0:
        for l in s[k:].splitlines()[2:]:
            mm = re.match(r"\s*\d+\s+([A-Za-z]{1,2})\s*:\s*(-?\d+\.\d+)", l)
            if not mm:
                break
            if mm.group(1) not in METALS:
                q.append(abs(float(mm.group(2))))
    bas = re.search(r"^\|\s*\d+>\s*!.*?(def2-SVPD|ma-def2-SVP|def2-SVP|def2-TZVPD|def2-TZVP)", s, re.M | re.I)
    if not bas:   # a composite method names no basis of its own; use the method keyword instead
        bas = re.search(r"^\|\s*\d+>\s*!\s*(?:QMMM\s+)?(\S*-3c|\S*3c)\b", s, re.M | re.I)
    return dict(nb=int(nb.group(1)) if nb else None, npc=int(pc.group(1)) if pc else 0,
                cyc=int(cy[-1]) if cy else None, dip=float(dp[-1]) if dp else None,
                gap=(vir[0] - occ[-1]) * EV if vir else None, homo=occ[-1] * EV,
                maxlw=max(q) if q else None, E=float(en[-1]),
                nopt=len(re.findall(r"GEOMETRY OPTIMIZATION CYCLE", s)),
                bas=bas.group(1) if bas else "?",
                norm="****ORCA TERMINATED NORMALLY****" in s, notconv="SCF NOT CONVERGED" in s)


def main(root, out):
    root = Path(root)
    if not (root / ".git").is_dir():
        stop("%s is not a git repo" % root)
    R = []
    for a in sorted(p for p in root.rglob("*.tar.xz") if ".git" not in p.parts):
        try:
            t = tarfile.open(str(a))
        except Exception:
            continue
        for m in t.getmembers():
            if not m.isfile() or not m.name.endswith(".out") or m.name.endswith(".pbs.out"):
                continue
            try:
                s = t.extractfile(m).read().decode(errors="replace")
            except Exception:
                continue
            d = read(s)
            if d:
                d["arc"] = str(a.relative_to(root)); d["run"] = m.name
                R.append(d)
    if not R:
        stop("found no completed ORCA single points on this dianion inside the committed archives")
    good = [d for d in R if d["norm"] and not d["notconv"]]

    CTL = {}
    for nb in {d["nb"] for d in good}:
        g = [d for d in good if d["nb"] == nb and d["npc"] == 0]
        if g:
            CTL[nb] = {k: st.median([x[k] for x in g if x[k] is not None])
                       for k in ("cyc", "gap", "maxlw", "dip") if any(x[k] is not None for x in g)}
            CTL[nb]["n"] = len(g)

    def verdict(d):
        c = CTL.get(d["nb"])
        if not c or "gap" not in c:
            return "NO CONTROL", []
        moved, sick = [], []
        if d["maxlw"] is not None and "maxlw" in c and d["maxlw"] - c["maxlw"] >= D_LW:
            moved.append("Loewdin +%.3f e over control" % (d["maxlw"] - c["maxlw"]))
        if d["dip"] is not None and c.get("dip") and d["dip"] / c["dip"] >= R_DIP:
            moved.append("dipole %.1fx control" % (d["dip"] / c["dip"]))
        if d["gap"] is not None and d["gap"] < GAP_ABS and d["gap"] < R_GAP * c["gap"]:
            sick.append("gap %.3f eV" % d["gap"])
        if d["cyc"] is not None and c.get("cyc") and d["cyc"] >= R_CYC * c["cyc"]:
            sick.append("%d SCF cycles, %.1fx control" % (d["cyc"], d["cyc"] / c["cyc"]))
        if moved and sick:
            return "COLLAPSED", moved + sick
        if sick:
            return "STRAINED", sick
        if moved:
            return "POLARISED", moved
        return "CLEAN", []

    for d in good:
        d["v"], d["w"] = verdict(d)

    cf = [d for d in good if d["npc"] == 0]
    # npc in the thousands is an MM environment in a QM/MM run, not a designed charge arrangement
    sp = [d for d in good if 0 < d["npc"] <= MM_NPC and d["nopt"] == 0]
    mm = [d for d in good if d["npc"] > MM_NPC]
    go = [d for d in good if d["npc"] > 0 and d["nopt"] > 0]

    O = ["STAGE 6 - DIAGNOSTIC 6 RESULT: a wavefunction-health audit of the whole committed charge-design record",
         "",
         "Measured by cs6_audit.py. NO NEW CALCULATION. Four fields are read out of every committed ORCA output",
         "on this dianion that Diagnostics 1 to 5 never read: the SCF cycle count, the HOMO-LUMO gap, the dipole",
         "magnitude and the largest Loewdin charge on a light atom. The rule is STAGE6_DIAGNOSTIC6_20261026.txt",
         "D6.3, applied unchanged.", "",
         "   %d completed single points and optimisations on this dianion (118 electrons, 59 occupied)" % len(good),
         "   %d carry no point charges and are the controls; %d are point-charge single points; %d are"
         % (len(cf), len(sp), len(go)),
         "   point-charge geometry optimisations, reported separately because a charge-free control is NOT",
         "   like-for-like once the geometry has moved. A further %d runs carry more than %d point charges and" % (len(mm), MM_NPC),
         "   are QM/MM calculations whose charges are an MM environment rather than a designed arrangement;",
         "   they are excluded from every table below.",
         "   %d outputs were skipped as not completed or not this system." % (len(R) - len(good)), "",
         "FOUR WORDS, AND WHAT EACH MEANS", "",
         "   CLEAN      no indicator departs from the same-basis charge-free control.",
         "   POLARISED  the density has moved - Loewdin up at least %.2f e over control, or the dipole at least"
         % D_LW,
         "              %.1f times it - but the SCF is healthy. THIS IS NOT A FAULT. It is what a strong"
         % R_DIP,
         "              external field is supposed to do, and it is reported so that it is not mistaken for one.",
         "   STRAINED   the SCF is unwell - gap below %.1f eV AND below %.0f per cent of its own control's gap,"
         % (GAP_ABS, 100 * R_GAP),
         "              or at least %.0f times the control's cycle count - without demonstrable relocation."
         % R_CYC,
         "              A near-degeneracy: single-reference closed-shell DFT is not reliable there.",
         "   COLLAPSED  both at once. The density has moved AND the SCF is unwell. This is the signature of",
         "              density being captured by a positive point charge that has no Pauli wall for electrons",
         "              (Finding 4), the effect Marefat Khah et al. JCTC 2020 doi:10.1021/acs.jctc.9b01162",
         "              describe as unphysical leaking of the electron density.", "",
         "THE CONTROLS EVERY THRESHOLD IS ANCHORED TO", "",
         "   %5s %6s %5s %10s %9s %9s   %s" % ("nb", "basis", "n", "med cyc", "med gap", "med maxLw", "med dip")]
    for nb in sorted(CTL):
        c = CTL[nb]
        bs = [d["bas"] for d in cf if d["nb"] == nb]
        bs = max(set(bs), key=bs.count) if bs else "?"
        O.append("   %5d %-6s %5d %10.0f %9.3f %9.3f   %6.2f"
                 % (nb, bs, c["n"], c.get("cyc", 0), c.get("gap", 0), c.get("maxlw", 0), c.get("dip", 0)))
    nbad = sum(1 for d in cf if d["v"] != "CLEAN")
    O += ["",
          "   SPECIFICITY CHECK: of the %d charge-free controls, %d are classified anything but CLEAN."
          % (len(cf), nbad),
          "   %s" % ("A control that trips the classifier would mean the thresholds are too tight; none does."
                     if nbad == 0 else "THAT IS NOT ZERO AND THE THRESHOLDS ARE TOO TIGHT - say so in any use of this table."),
          "",
          "POINT-CHARGE SINGLE POINTS, BY BASIS SET", "",
          "   %5s %-12s %6s %8s %7s %9s %8s %9s" % ("nb", "basis", "runs", "npc med", "clean", "polarised", "strained", "collapsed")]
    for nb in sorted({d["nb"] for d in sp}):
        g = [d for d in sp if d["nb"] == nb]
        bs = [d["bas"] for d in g]
        n = lambda v: sum(1 for d in g if d["v"] == v)
        O.append("   %5d %-12s %6d %8d %7d %9d %8d %9d"
                 % (nb, max(set(bs), key=bs.count), len(g), st.median([d["npc"] for d in g]),
                    n("CLEAN"), n("POLARISED"), n("STRAINED"), n("COLLAPSED")))

    O += ["", "AT HOW MANY CHARGES DOES EACH BASIS COLLAPSE?", "",
          "   This is the reading that answers why the representation passed its own tests and then failed.", "",
          "   %5s %-12s %10s %28s" % ("nb", "basis", "collapsed", "number of point charges in them")]
    for nb in sorted({d["nb"] for d in sp if d["v"] == "COLLAPSED"}):
        g = [d for d in sp if d["nb"] == nb and d["v"] == "COLLAPSED"]
        v = sorted(d["npc"] for d in g)
        bs = [d["bas"] for d in g]
        O.append("   %5d %-12s %10d %28s" % (nb, max(set(bs), key=bs.count), len(g),
                 "%d to %d, median %d" % (v[0], v[-1], st.median(v))))
    for nb in sorted({d["nb"] for d in sp}):
        g = [d for d in sp if d["nb"] == nb]
        if not any(d["v"] == "COLLAPSED" for d in g):
            v = sorted(d["npc"] for d in g)
            bs = [d["bas"] for d in g]
            O.append("   %5d %-12s %10d %28s" % (nb, max(set(bs), key=bs.count), 0,
                     "none, over npc %d to %d" % (v[0], v[-1])))

    O += ["", "INDICATORS AGAINST THE NUMBER OF POINT CHARGES, def2-SVP ONLY (nb 264)", "",
          "   def2-SVP has no diffuse functions, so it is the basis every committed Stage 6 design was",
          "   evaluated at, and the one whose quiet failure matters most.", "",
          "   %8s %6s %10s %10s %11s %12s %10s"
          % ("npc", "runs", "med cyc", "med gap", "worst gap", "med maxLw", "collapsed")]
    for lo, hi, lab in ((1, 1, "1"), (2, 4, "2-4"), (5, 9, "5-9"), (10, 13, "10-13"),
                        (14, 29, "14-29"), (30, 10 ** 6, ">=30")):
        g = [d for d in sp if d["nb"] == 264 and lo <= d["npc"] <= hi]
        if not g:
            continue
        md = lambda k: (lambda v: "%.3f" % st.median(v) if v else "-")([x[k] for x in g if x[k] is not None])
        O.append("   %8s %6d %10s %10s %11s %12s %10d"
                 % (lab, len(g), md("cyc"), md("gap"),
                    "%.3f" % min([x["gap"] for x in g if x["gap"] is not None], default=0),
                    md("maxlw"), sum(1 for d in g if d["v"] == "COLLAPSED")))

    O += ["", "EVERY COLLAPSED OR STRAINED SINGLE POINT", "",
          "   %-38s %-28s %4s %4s %4s  %s" % ("archive", "run", "npc", "nb", "cyc", "verdict and why")]
    for d in sorted([x for x in sp if x["v"] in ("COLLAPSED", "STRAINED")],
                    key=lambda x: (x["v"], x["arc"], x["run"])):
        O.append("   %-38s %-28s %4d %4d %4s  %s: %s"
                 % (d["arc"][-38:], d["run"].split("/")[-1][:28], d["npc"], d["nb"],
                    d["cyc"] if d["cyc"] is not None else "-", d["v"], "; ".join(d["w"])))

    O += ["", "THE POINT-CHARGE GEOMETRY OPTIMISATIONS, REPORTED APART", "",
          "   A charge-free control is not like-for-like once the geometry has moved, so these verdicts carry",
          "   much less weight than the single points above and must not be quoted as collapses on their own.", ""]
    for v in ("CLEAN", "POLARISED", "STRAINED", "COLLAPSED", "NO CONTROL"):
        n = sum(1 for d in go if d["v"] == v)
        if n:
            O.append("   %-11s %4d" % (v, n))

    O += ["", "WHAT THIS DOES NOT TEST", "",
          "   Nothing was computed, so this cannot test any hypothesis the committed outputs do not already",
          "   contain. It cannot say what the correct barrier is for a run it calls COLLAPSED or STRAINED - only",
          "   that the number should not be quoted. It does not establish a usable charge model; a Pauli wall",
          "   for electrons is the fix and it is untested here. The four thresholds were chosen AFTER the Stage 6",
          "   outputs had been inspected, which is why they are anchored to in-corpus controls and why the",
          "   control table and the specificity check are printed: the classification is reproducible and",
          "   checkable, but it is not a blind test. A Loewdin charge is basis-dependent and is used here only",
          "   as a difference against a same-basis control, never as an absolute quantity.", ""]

    Path(out).mkdir(parents=True, exist_ok=True)
    (Path(out) / "DIAGNOSTIC6_RESULT.txt").write_text("\n".join(O) + "\n")
    print("\n".join(O))
    print()
    print("   wrote %s" % (Path(out) / "DIAGNOSTIC6_RESULT.txt"))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        stop("usage: cs6_audit.py <results-repo> <out-dir>")
    main(sys.argv[1], sys.argv[2])
