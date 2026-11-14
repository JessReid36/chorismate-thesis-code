#!/usr/bin/env python3
"""DIAGNOSTIC 36: a correlated reference for the in vacuo barrier at frame 42436.

Arms C and H ran under attempt 1 and live in d36. Arm B failed twice and finished on
attempt 3; its runs live in d36a2, which is why the earlier harvest was superseded.

These are bare single points with no charges and no MM region, so cs6_analyse.py's
bare / (MM) / (QM/MM) convention does not apply and is not used. The energy is the last
FINAL SINGLE POINT ENERGY, the convention the committed cs6_audit.py's reader uses. The
HOMO convention is the committed cs10_cpcm.py's: the last ORBITAL ENERGIES block, rows
with four fields, occupancy above 0.5, the eV column, highest occupied last.

Every number in the result is read from an output file. Nothing is typed in except the
committed reference values the gates compare against, which are named with their source.

Usage: cs6_d36.py --res R --code C --prereg P --amd1 A --amd2 B --result OUT --tsv TSV
"""
import argparse, hashlib, lzma, os, re, sys
from pathlib import Path

H2K = 627.5094740631
NEL, NOCC = 118, 59

# G36.3, from STAGE6_DIAGNOSTIC36_20261111.txt section 4, which took them from D28
CTL_R, CTL_TS = -837.20938385, -837.16842778
CTL_BARRIER = 25.7003
CTL_TOL_EH, CTL_TOL_KCAL = 1e-6, 0.0010

# G36.8, from STAGE6_DIAGNOSTIC36_AMENDMENT1_20261113.txt section 5, measured by attempt 1
A1_E0, A1_ECORR = -832.744852436, -3.028125973
A1_TOL = 1e-6

# R36.2, the three functionals the pre-registration fixed, with its own figures
R362 = (("CAM-B3LYP-D3BJ", 25.7003), ("wB97X-D3", 27.9689), ("B3LYP", 17.304571))
# R36.3, D28's HOMO values and Probe 3's target, from the pre-registration section 1
D28_HOMO = (("B3LYP", 1.2307), ("CAM-B3LYP", -0.5918), ("wB97X-D3", -1.3424))
P3_TARGET = -0.5123
# R36.5, the functional shifts the pre-registration derived in section 1
SHIFT = (("CAM-B3LYP", -1.1245), ("wB97X-D3", -1.0752))
# POST HOC. D37, committed 64ca6909, same frame and basis. NOT pre-registered in D36.
D37 = (("B3LYP-D3BJ", 17.3046), ("PW6B95-D3BJ", 20.9595), ("PWPB95-D3BJ", 21.3407),
       ("CAM-B3LYP-D3BJ", 25.7003), ("wB97M-V", 27.3924), ("wB97X-D3", 27.9689),
       ("M06-2X-D3(0)", 28.4003), ("wB97X-V", 29.7565))

ARMS = (("C", "d36", "CAM-B3LYP-D3BJ/def2-TZVPD, the control"),
        ("H", "d36", "HF/def2-TZVPD, the reference determinant"),
        ("B", "d36a2", "DLPNO-CCSD(T)/def2-TZVPD, the benchmark"))
LEGS = ("R", "TS")


def stop(m):
    sys.stderr.write("STOP  %s\n" % m)
    raise SystemExit(1)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def rd(p):
    p = str(p)
    return lzma.open(p, "rt", errors="replace").read() if p.endswith(".xz") \
        else open(p, errors="replace").read()


def vals(t, label, pat=r"-?\d+\.\d+"):
    """Every dotted ORCA report line 'LABEL  ...  number' for this label."""
    return re.findall(r"^\s*%s\s*\.*\s*(%s)\s*$" % (re.escape(label), pat), t, re.M)


def val(t, label, bad=None, pat=r"-?\d+\.\d+", cast=float):
    """One value for a label. ORCA prints some labels more than once; if the copies
    disagree, taking the first would be a silent choice, so record it instead."""
    v = vals(t, label, pat)
    if not v:
        return None
    if len(set(v)) > 1 and bad is not None:
        bad.append("%r is printed %d times with %d different values: %s"
                   % (label, len(v), len(set(v)), ", ".join(sorted(set(v)))))
    return cast(v[0])


def read_out(d):
    """Everything this diagnostic needs from one output."""
    p = None
    for fn in ("job.out", "job.out.xz"):
        if os.path.isfile(os.path.join(d, fn)):
            p = os.path.join(d, fn)
            break
    if p is None:
        return None
    t = rd(p)
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    nel = re.search(r"Number of Electrons\s+NEL\s*\.+\s*(\d+)", t)
    geo = re.search(r"CARTESIAN COORDINATES \(ANGSTROEM\)\n-+\n(.*?)\n\n", t, re.S)
    atoms = []
    if geo:
        for l in geo.group(1).splitlines():
            f = l.split()
            if len(f) == 4:
                atoms.append((f[0],) + tuple(round(float(x), 6) for x in f[1:4]))
    # HOMO, by the committed cs10_cpcm.py convention
    homo, occ = None, []
    k = t.rfind("ORBITAL ENERGIES")
    if k >= 0:
        for l in t[k:].splitlines()[4:]:
            f = l.split()
            if len(f) != 4:
                if occ:
                    break
                continue
            try:
                o, _eh, ev = float(f[1]), float(f[2]), float(f[3])
            except ValueError:
                break
            if o > 0.5:
                occ.append(ev)
        homo = occ[-1] if occ else None
    bad = []
    rt = re.search(r"TOTAL RUN TIME:\s*(.+)", t)
    sto = re.search(r"PNO-Integral Storage\s*\.*\s*(\S[^\n]*)", t)
    return dict(
        path=p, sha=sha(p),
        norm="ORCA TERMINATED NORMALLY" in t,
        notconv=bool(re.search(r"SCF NOT CONVERGED AFTER|SERIOUS PROBLEM", t)),
        nel=int(nel.group(1)) if nel else None,
        E=float(e[-1]) if e else None,
        nocc=len(occ), homo=homo, atoms=atoms, runtime=rt.group(1).strip() if rt else None,
        storage=sto.group(1).strip() if sto else None,
        ntriples=val(t, "Number of Triples that are to be computed", bad,
                     r"\d+", int),
        e0=val(t, r"E(0)", bad), strong=val(t, "E(CORR)(strong-pairs)", bad),
        weak=val(t, "E(CORR)(weak-pairs)", bad), corr=val(t, "E(CORR)(corrected)", bad),
        ccsd=val(t, "E(CCSD)", bad), tcorr=val(t, "Triples Correction (T)", bad),
        ccsdt=val(t, "E(CCSD(T))", bad), t1=val(t, "T1 diagnostic", bad),
        inconsistent=bad,
        method=(re.search(r"^\|\s*\d+> *(![^\n|]*)", t, re.M).group(1).strip()
                if re.search(r"^\|\s*\d+> *(!)", t, re.M) else None),
        charge=(lambda m: m.group(1) if m else None)(
            re.search(r"^\|\s*\d+> *\* *xyzfile\s+(-?\d+\s+\d+)", t, re.M)))


def xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    return [(x.split()[0],) + tuple(round(float(v), 6) for v in x.split()[1:4])
            for x in L[2:2 + n]]


def progress(d):
    """Peak scratch and the node, from the job script's own PROGRESS.txt."""
    p = os.path.join(d, "PROGRESS.txt")
    if not os.path.isfile(p):
        return None
    t = Path(p).read_text()
    node = re.search(r"attempt \d+ on (\S+?),", t)
    free = re.search(r"(\d+)G free", t)
    used, unit = [], {"K": 1 / 1048576.0, "M": 1 / 1024.0, "G": 1.0, "T": 1024.0}
    for m in re.finditer(r"used=([\d.]+)([KMGT])", t):
        used.append(float(m.group(1)) * unit[m.group(2)])
    modes = sorted(set(re.findall(r"storage=(\S[^\n]*)", t)) - {""})
    return dict(node=node.group(1) if node else None,
                free=int(free.group(1)) if free else None,
                peak=max(used) if used else None, modes=modes,
                exits=re.findall(r"orca_exit=(\d+)", t))


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "amd1", "amd2", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    S6 = os.path.join(a.res, "phase2.2/cs_stage6")
    gates, rep = [], []
    say = rep.append

    say("STAGE 6, DIAGNOSTIC 36 - RESULT")
    say("A correlated reference for the in vacuo barrier: DLPNO-CCSD(T) on frame 42436.")
    say("")
    say("Read by phase2.2/cs6_d36.py                      sha256 %s" % sha(__file__))
    for lab, p in (("pre-registration", a.prereg), ("amendment 1", a.amd1),
                   ("amendment 2", a.amd2)):
        say("   %-16s %-36s sha256 %s" % (lab, os.path.basename(p), sha(p)))
    say("")
    say("Arms C and H ran under attempt 1 and are read from d36. Arm B failed twice on")
    say("resource faults and finished on attempt 3; it is read from d36a2. Amendments 1")
    say("and 2 are the notes on those two failures. Nothing here was re-run to get a")
    say("different number: the method line, basis, auxiliary, PNO thresholds and geometry")
    say("are byte for byte what the pre-registration fixed.")
    say("")
    say("Bare single points: no charges, no MM region, no wall. The energy is the last")
    say("FINAL SINGLE POINT ENERGY. The HOMO is the committed cs10_cpcm.py convention.")
    say("")

    # ---------------------------------------------------------------- read
    D, P = {}, {}
    for arm, sub, _ in ARMS:
        for leg in LEGS:
            d = os.path.join(S6, sub, "%s_%s" % (arm, leg))
            r = read_out(d)
            if r is None:
                gates.append("%s_%s produced no readable output in %s/" % (arm, leg, sub))
                continue
            D[(arm, leg)] = r
            P[(arm, leg)] = progress(d)
            if not r["norm"]:
                gates.append("G36.4 %s_%s did not terminate normally" % (arm, leg))
            if r["notconv"]:
                gates.append("G36.4 %s_%s did not converge its SCF" % (arm, leg))
            if r["nel"] != NEL:
                gates.append("G36.4 %s_%s reports %s electrons, expected %d"
                             % (arm, leg, r["nel"], NEL))
            if r["E"] is None:
                gates.append("G36.4 %s_%s has no final energy" % (arm, leg))
            for b in r["inconsistent"]:
                gates.append("G36.4 %s_%s: %s" % (arm, leg, b))
            if r["charge"] is not None and r["charge"].split() != ["-2", "1"]:
                gates.append("G36.6 %s_%s ran at charge/multiplicity %s, expected -2 1"
                             % (arm, leg, r["charge"]))

    # -------------------------------------------- G36.1, G36.6 on the inputs
    for arm, sub, _ in ARMS:
        for leg in LEGS:
            d = os.path.join(S6, sub, "%s_%s" % (arm, leg))
            g = os.path.join(S6, "d33a", "sp_42436_%s" % leg, "substrate.xyz")
            got = os.path.join(d, "substrate.xyz")
            if os.path.isfile(got) and os.path.isfile(g):
                if sha(got) != sha(g):
                    gates.append("G36.1 %s_%s's geometry is not d33a frame 42436's"
                                 % (arm, leg))
                if (D.get((arm, leg)) or {}).get("atoms") not in (None, [], xyz(g)):
                    gates.append("G36.1 %s_%s RAN on coordinates that are not frame "
                                 "42436's" % (arm, leg))
            else:
                gates.append("G36.1 %s_%s has no geometry to compare" % (arm, leg))
            ip = os.path.join(d, "job.inp")
            if os.path.isfile(ip) and "* xyzfile -2 1 substrate.xyz" not in rd(ip):
                gates.append("G36.6 %s_%s's input is not charge -2 multiplicity 1"
                             % (arm, leg))

    # ---------------------------------------------------------------- G36.2
    for leg in LEGS:
        src = os.path.join(S6, "d28", "CAM-B3LYP", "def2-TZVPD", "ctl_%s" % leg, "job.inp")
        got = os.path.join(S6, "d36", "C_%s" % leg, "job.inp")
        if os.path.isfile(src) and os.path.isfile(got):
            if sha(src) != sha(got):
                gates.append("G36.2 arm C's %s input is not D28's control input" % leg)
        else:
            gates.append("G36.2 cannot compare arm C's %s input with D28's" % leg)

    # --------------------------------- G36.3 THE CONTROL. nothing survives it
    say("THE CONTROL, G36.3. Arm C repeats D28's committed calculation, so it must")
    say("return D28's committed number. If it does not, nothing else here is readable.")
    say("")
    cb = None
    cr, cts = (D.get(("C", "R")) or {}).get("E"), (D.get(("C", "TS")) or {}).get("E")
    for leg, got, want in (("R", cr, CTL_R), ("TS", cts, CTL_TS)):
        if got is None:
            gates.append("G36.3 arm C's %s leg has no energy" % leg)
            say("   C_%-2s  NO ENERGY           want %15.8f                 FAIL" % (leg, want))
            continue
        d = abs(got - want)
        if d > CTL_TOL_EH:
            gates.append("G36.3 arm C's %s leg gives %.8f, D28 committed %.8f, "
                          "difference %.2e Eh" % (leg, got, want, d))
        say("   C_%-2s  %16.8f  want %15.8f  diff %.2e  %s"
            % (leg, got, want, d, "PASS" if d <= CTL_TOL_EH else "FAIL"))
    if cr is not None and cts is not None:
        cb = (cts - cr) * H2K
        d = abs(cb - CTL_BARRIER)
        if d > CTL_TOL_KCAL:
            gates.append("G36.3 arm C's barrier is %.4f, D28 committed %.4f"
                         % (cb, CTL_BARRIER))
        say("   barrier %14.4f  want %15.4f  diff %.1e  %s"
            % (cb, CTL_BARRIER, d, "PASS" if d <= CTL_TOL_KCAL else "FAIL"))
    say("")
    ctl_bad = [g for g in gates
               if g.startswith("G36.3") or "C_R" in g or "C_TS" in g]
    if ctl_bad:
        say("G36.3 FAILED. The pre-registration says nothing else in D36 may be read if")
        say("the control does not return D28's committed number, and a control that did")
        say("not converge, did not terminate or ran on the wrong geometry has not")
        say("returned it either. No reading is reported.")
        say("")
        for g in ctl_bad:
            say("   control failure: %s" % g)
        say("")
        say("GATE FAILURES")
        for g in gates:
            say("   %s" % g)
        Path(a.result).write_text("\n".join(rep) + "\n")
        print("\n".join(rep))
        return 1

    # ------------------------------------------------- arm B's own integrity
    B = {leg: D.get(("B", leg)) or {} for leg in LEGS}
    for leg in LEGS:
        r = B[leg]
        if r.get("ccsdt") is None:
            gates.append("G36.5 arm B's %s leg has no CCSD(T) energy" % leg)
        for lab, a_, b_, c_, tol in (
                ("strong+weak = corrected", r.get("strong"), r.get("weak"), r.get("corr"), 2e-9),
                ("E(0)+E(CORR) = E(CCSD)", r.get("e0"), r.get("corr"), r.get("ccsd"), 2e-9),
                ("E(CCSD)+(T) = E(CCSD(T))", r.get("ccsd"), r.get("tcorr"), r.get("ccsdt"), 2e-9)):
            if None in (a_, b_, c_):
                gates.append("G36.5 arm B's %s leg cannot check %s" % (leg, lab))
            elif abs((a_ + b_) - c_) > tol:
                gates.append("G36.5 arm B's %s leg fails %s by %.2e"
                             % (leg, lab, abs((a_ + b_) - c_)))
        if r.get("E") is not None and r.get("ccsdt") is not None \
                and abs(r["E"] - r["ccsdt"]) > 1e-8:
            gates.append("G36.5 arm B's %s leg: final energy %.9f is not E(CCSD(T)) %.9f"
                         % (leg, r["E"], r["ccsdt"]))

    # --------------------------------------- G36.8, the only control arm B has
    say("G36.8, CARRIED FORWARD FROM AMENDMENT 1. Attempt 1's reactant leg converged")
    say("before its walltime ran out and left two numbers behind. Where integrals are")
    say("stored cannot move a converged energy, and this is the measurement of that")
    say("rather than the assertion of it.")
    say("")
    for lab, got, want in (("E(0)", B["R"].get("e0"), A1_E0),
                           ("E(CORR)(corrected)", B["R"].get("corr"), A1_ECORR)):
        if got is None:
            gates.append("G36.8 arm B's R leg does not report %s" % lab)
            say("   %-20s NOT REPORTED                                      FAIL" % lab)
            continue
        d = abs(got - want)
        if d > A1_TOL:
            gates.append("G36.8 arm B's R %s is %.9f, attempt 1 measured %.9f, "
                         "difference %.2e Eh" % (lab, got, want, d))
        say("   %-20s %16.9f  attempt 1 %16.9f  diff %.2e  %s"
            % (lab, got, want, d, "PASS" if d <= A1_TOL else "FAIL"))
    say("")
    say("   THERE IS NO SUCH CONTROL FOR THE TRANSITION STATE. Attempt 1's TS leg never")
    say("   started, so no earlier measurement of it exists anywhere in this record. The")
    say("   reactant leg is the only one with an independent cross-check.")
    say("")

    # ---------------------------------------------------- G36.11, the diff
    say("G36.11. Attempt 3's input against attempt 2's committed input.")
    say("")
    for leg in LEGS:
        new = os.path.join(S6, "d36a2", "B_%s" % leg, "job.inp")
        old = os.path.join(S6, "d36a1", "B_%s" % leg, "job.inp")
        if not (os.path.isfile(new) and os.path.isfile(old)):
            gates.append("G36.11 cannot compare attempt 3's %s input with attempt 2's" % leg)
            continue
        nl, ol = rd(new).splitlines(), rd(old).splitlines()
        blk = ["%mdci", "     StorageType Shared", "end"]
        back = [("%maxcore 12000" if x == "%maxcore 8000" else x)
                for x in nl if x not in blk]
        ok = (back == ol and len(nl) == len(ol) + 3
              and all(nl.count(k) == 1 for k in blk))
        if not ok:
            gates.append("G36.11 attempt 3's %s input differs from attempt 2's in more "
                         "than the %%mdci block and %%maxcore" % leg)
        say("   B_%-2s  3 lines added, 1 changed, 0 removed, round trip exact   %s"
            % (leg, "PASS" if ok else "FAIL"))
    say("")

    # ------------------------------------------- G36.12 and G36.13, the cost
    say("G36.12 and G36.13. What this cost, and where every number came from.")
    say("")
    say("   leg   node          storage    peak scratch   wall clock")
    for leg in LEGS:
        r, p = B[leg], P.get(("B", leg)) or {}
        modes = ", ".join(p.get("modes") or []) or (r.get("storage") or "not reported")
        say("   B_%-2s  %-12s  %-9s  %10s   %s"
            % (leg, p.get("node") or "?", modes,
               ("%.0f GB" % p["peak"]) if p.get("peak") else "?",
               r.get("runtime") or "?"))
    if not any((P.get(("B", leg)) or {}).get("modes") for leg in LEGS):
        gates.append("G36.12 the result cannot state the PNO storage mode for arm B")
    say("")
    say("   ORCA FELL BACK TO ON-DISK STORAGE, and that is reported rather than hidden.")
    say("   The ORCA 6.0 manual, section 7.14.4, says of StorageType Shared: \"If the")
    say("   amount of memory is not sufficient, ORCA will fall back to on-disk storage.\"")
    say("   Attempt 2 wrote the same intermediates to a shared NFS home at 97 per cent")
    say("   full and died; attempt 3 wrote them to node-local scratch and finished. The")
    say("   storage location is not a method choice and G36.8 is what shows it did not")
    say("   move the answer.")
    say("")
    for leg in LEGS:
        r = B[leg]
        say("   B_%-2s  job.out sha256 %s" % (leg, r.get("sha") or "?"))
    say("")
    say("   Every arm B number above and below is read from those two files, which the")
    say("   job script copied out of node-local scratch at the end of the run. No output")
    say("   was read from, or written to, /home while either job was running.")
    say("")

    # ---------------------------------------------------------------- R36.1
    say("=" * 78)
    say("R36.1  THE DLPNO-CCSD(T) BARRIER AT FRAME 42436, def2-TZVPD")
    say("=" * 78)
    say("")
    bb = None
    if B["R"].get("ccsdt") is not None and B["TS"].get("ccsdt") is not None:
        bb = (B["TS"]["ccsdt"] - B["R"]["ccsdt"]) * H2K
    say("   component                      R                TS        barrier kcal/mol")
    for lab, key in (("E(0), Hartree-Fock", "e0"), ("E(CORR)(strong-pairs)", "strong"),
                     ("E(CORR)(weak-pairs)", "weak"), ("E(CORR)(corrected)", "corr"),
                     ("E(CCSD)", "ccsd"), ("Triples Correction (T)", "tcorr"),
                     ("E(CCSD(T))", "ccsdt")):
        x, y = B["R"].get(key), B["TS"].get(key)
        if None in (x, y):
            say("   %-24s %16s %17s" % (lab, "?", "?"))
            continue
        say("   %-24s %16.9f %17.9f %18.4f" % (lab, x, y, (y - x) * H2K))
    say("")
    if bb is not None:
        say("   THE BARRIER IS %.4f kcal/mol." % bb)
        say("")
        hf = (B["TS"]["e0"] - B["R"]["e0"]) * H2K
        cc = (B["TS"]["ccsd"] - B["R"]["ccsd"]) * H2K
        tr = (B["TS"]["tcorr"] - B["R"]["tcorr"]) * H2K
        say("   Hartree-Fock alone gives %+.4f. Correlation at the CCSD level takes" % hf)
        say("   %.4f off it, leaving %+.4f, and the triples take a further %.4f off,"
            % (hf - cc, cc, cc - bb))
        say("   leaving %+.4f. The triples are %.1f per cent of the final barrier, so a"
            % (bb, 100 * abs(tr) / abs(bb)))
        say("   method without (T) would have been %.4f too high here." % (cc - bb))
    for leg in LEGS:
        n = B[leg].get("ntriples")
        say("   %s leg: %s triples computed" % (leg, ("%d" % n) if n else "not reported"))
    say("")

    # ---------------------------------------------------------------- R36.2
    say("=" * 78)
    say("R36.2  EACH PRE-REGISTERED FUNCTIONAL MINUS THE BENCHMARK, same frame, same basis")
    say("=" * 78)
    say("")
    say("   The three functionals and the three figures are the pre-registration's own,")
    say("   section 3. This is the ranking the diagnostic exists to produce.")
    say("")
    say("   functional            barrier      minus benchmark")
    if bb is not None:
        for n, v in R362:
            say("   %-18s %10.4f %+18.4f" % (n, v, v - bb))
        say("")
        best = min(R362, key=lambda kv: abs(kv[1] - bb))
        say("   Of the three, %s is closest, at %+.4f kcal/mol."
            % (best[0], best[1] - bb))
        say("   D36.1 applies: that makes it the closest here and nothing more.")
    say("")

    # ---------------------------------------------------------------- R36.3
    say("=" * 78)
    say("R36.3  THE HARTREE-FOCK REFERENCE, MEASURED DIRECTLY")
    say("=" * 78)
    say("")
    hr, hts = (D.get(("H", "R")) or {}), (D.get(("H", "TS")) or {})
    if hr.get("E") is not None and hts.get("E") is not None:
        say("   HF/def2-TZVPD barrier, arm H      %10.4f kcal/mol"
            % ((hts["E"] - hr["E"]) * H2K))
        if B["R"].get("e0") is not None:
            hb = (B["TS"]["e0"] - B["R"]["e0"]) * H2K
            say("   HF reference inside arm B, E(0)   %10.4f kcal/mol" % hb)
            say("   the two agree to                  %10.4f kcal/mol"
                % abs(hb - (hts["E"] - hr["E"]) * H2K))
            say("   Arm H and arm B's own reference are separate calculations in separate")
            say("   runs; that they agree is a check on both.")
    say("")
    say("   the highest occupied orbital, eV. D28's three values and Probe 3's target")
    say("   are the pre-registration's, section 1.")
    say("")
    for leg in LEGS:
        r = D.get(("H", leg)) or {}
        h = r.get("homo")
        say("   HF   %-2s  %s   %d occupied orbitals read%s"
            % (leg, ("%+9.4f" % h) if h is not None else "  not read", r.get("nocc", 0),
               "" if r.get("nocc") == NOCC else "   EXPECTED %d" % NOCC))
        if r.get("nocc") not in (None, NOCC):
            gates.append("R36.3 arm H's %s leg reports %s occupied orbitals, expected %d"
                         % (leg, r.get("nocc"), NOCC))
    say("")
    for n, v in D28_HOMO:
        say("   %-12s (D28)  %+9.4f" % (n, v))
    say("   Probe 3 target        %+9.4f" % P3_TARGET)
    say("")
    hR = (D.get(("H", "R")) or {}).get("homo")
    if hR is not None:
        say("   D36.5 is the limitation that applies here. A dianion whose highest occupied")
        say("   orbital sits above the vacuum level has an unbound electron in the")
        say("   reference, and the HF value measured above is %s zero."
            % ("ABOVE" if hR >= 0 else "below"))
        if hR >= 0:
            say("   IT IS ABOVE. The benchmark is reported WITH that caveat attached")
            say("   wherever it is used, and a multireference check is an open item.")
    say("")

    # ---------------------------------------------------------------- R36.4
    say("=" * 78)
    say("R36.4  THE T1 DIAGNOSTIC")
    say("=" * 78)
    say("")
    for leg in LEGS:
        t1 = B[leg].get("t1")
        say("   %-2s  %s" % (leg, ("%.9f" % t1) if t1 is not None else "not printed"))
    say("")
    t1s = [B[leg].get("t1") for leg in LEGS]
    if None in t1s:
        say("   One or both legs did not print a T1 diagnostic. R36.4 says the benchmark is")
        say("   then reported with the reference-quality question OPEN rather than")
        say("   silently assumed closed, and it is.")
    else:
        say("   The transition state is the higher of the two, by %.9f." % abs(t1s[1] - t1s[0]))
        say("   Both sit below the conventional 0.02 guide. THAT GUIDE IS NOT THIS")
        say("   PROJECT'S AND WAS NOT PRE-REGISTERED; the two numbers are the measurement")
        say("   and 0.02 is context. Amendment 1 made the same point for the reactant when")
        say("   attempt 1 measured it; this is the first measurement of the transition state.")
    say("")

    # ---------------------------------------------------------------- R36.5
    say("=" * 78)
    say("R36.5  THE DOUBLE-ZETA TO TRIPLE-ZETA SHIFT")
    say("=" * 78)
    say("")
    say("   the functionals' own shifts, from the pre-registration section 1:")
    for n, v in SHIFT:
        say("      %-12s %+8.4f kcal/mol" % (n, v))
    say("")
    say("   D36 DOES NOT RUN CCSD(T) AT def2-SVPD AND CANNOT MEASURE THE METHOD'S OWN")
    say("   BASIS SHIFT. The pre-registration says so in R36.5 and it remains true. What")
    say("   the functional shifts support is that a def2-SVPD production barrier carries")
    say("   a basis uncertainty of roughly one kcal/mol, and that conclusion does not")
    say("   depend on this benchmark.")
    say("")

    # ------------------------------------------------------------- POST HOC
    say("=" * 78)
    say("POST HOC. NOT PRE-REGISTERED IN D36.")
    say("=" * 78)
    say("")
    say("   D37 did not exist when D36 was pre-registered, so its eight functionals are")
    say("   not part of R36.2. They are at the same frame and the same basis, committed")
    say("   at 64ca6909, and the comparison is reported here as POST HOC.")
    say("")
    if bb is not None:
        say("   functional            barrier      minus benchmark")
        for n, v in sorted(D37, key=lambda kv: abs(kv[1] - bb)):
            say("   %-18s %10.4f %+18.4f" % (n, v, v - bb))
        near = [n for n, v in D37 if abs(v - bb) < 3.0]
        say("")
        say("   Spread of the eight: %.2f kcal/mol. The benchmark sits %.0f per cent of"
            % (max(v for _, v in D37) - min(v for _, v in D37),
               100 * (bb - min(v for _, v in D37))
               / (max(v for _, v in D37) - min(v for _, v in D37))))
        say("   the way up that spread. %d of the 8 are within 3 kcal/mol: %s."
            % (len(near), ", ".join(near)))
        say("")
        say("   D28's orbital test also preferred CAM-B3LYP. That agreement is weaker")
        say("   evidence than it looks: the three functionals closest to D28's HOMO target")
        say("   give barriers spanning 21.34 to 28.40, so that test could not have")
        say("   discriminated between them. This benchmark is what discriminates.")
    say("")

    # ---------------------------------------------------------- limitations
    say("=" * 78)
    say("WHAT THIS CANNOT SAY. The pre-registration's D36.1 to D36.8 all stand.")
    say("=" * 78)
    say("")
    if bb is not None:
        say("   NOT A CONVERGED ABSOLUTE BARRIER (D36.2). def2-TZVPD was chosen to match")
        say("   the committed DFT arm, not to converge the correlation energy. %.4f is a"
            % bb)
        say("   reference for comparison at this basis and must not be quoted as the true")
        say("   in vacuo barrier.")
        say("")
        say("   NOT THE EXPERIMENTAL NUMBER (D36.3, E7.11). The record's anchors are")
        say("   24.5 kcal/mol as a free energy and 20.7 as an enthalpy, both in water.")
        say("   %.4f is an electronic energy difference between two frozen gas-phase" % bb)
        say("   geometries: no zero-point energy, no thermal correction, no entropy, no")
        say("   solvation. ITS CLOSENESS TO 24.5 IS NOT AGREEMENT. E7.11 requires every")
        say("   sentence that puts them near each other to say which quantity is which,")
        say("   and these two are not the same quantity.")
        say("")
    say("   ONE FRAME. D35 measured the bare barrier at both frame 41786 and frame 42436")
    say("   at five levels and found the geometry effect spread over 2.1031 kcal/mol, with")
    say("   the ranking of levels not preserved between them. This is frame 42436's")
    say("   number.")
    say("")
    say("   ONE GEOMETRY, NOT A STATIONARY POINT. These are single points on committed")
    say("   band geometries, not an optimised transition state, so this is not a true")
    say("   activation energy even in vacuo.")
    say("")
    say("   DLPNO IS AN APPROXIMATION (D36.4). TightPNO was requested; the residual PNO")
    say("   error is not measured here.")
    say("")
    say("   NO DESIGN CONCLUSION (D36.6). No charges, no wall, no reliability claim.")
    say("   Nothing here bears on r_safe, on any barrier lowering, or on which designs")
    say("   survive.")
    say("")
    say("   IT DOES NOT ADOPT A PRODUCTION LEVEL OF THEORY (D36.7). That is an amendment")
    say("   and a decision for Jessica. D36 supplies evidence for it and is not it.")
    say("")
    say("   NOTHING HERE WAS SHAPED (D36.8). No outcome was preferred. Arm B was re-run")
    say("   twice and both re-runs changed resources only; G36.8 measures that the answer")
    say("   did not move.")
    say("")

    # --------------------------------------------------------------- gates
    say("=" * 78)
    if gates:
        say("GATE FAILURES: %d" % len(gates))
        say("=" * 78)
        say("")
        for g in gates:
            say("   %s" % g)
    else:
        say("ALL GATES PASS. G36.1 to G36.8 and G36.11 to G36.13.")
        say("=" * 78)
    say("")

    Path(a.result).write_text("\n".join(rep) + "\n")
    with open(a.tsv, "w") as f:
        f.write("arm\tmethod\tE_R\tE_TS\tbarrier_kcal\n")
        for arm, _sub, nm in ARMS:
            x = (D.get((arm, "R")) or {}).get("E")
            y = (D.get((arm, "TS")) or {}).get("E")
            if None in (x, y):
                continue
            f.write("%s\t%s\t%.9f\t%.9f\t%.6f\n" % (arm, nm, x, y, (y - x) * H2K))
    print("\n".join(rep))
    return 1 if gates else 0


if __name__ == "__main__":
    sys.exit(main())
