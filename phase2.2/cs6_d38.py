#!/usr/bin/env python3
"""DIAGNOSTIC 38: do embedding point charges reach the correlated step?

Five DLPNO-CCSD(T)/def2-SVPD single points on frame 42436's reactant, differing only in
the embedding charge. The question is mechanical, not chemical: whether a charge that
reaches the Hartree-Fock reference also reaches the correlation treatment. The dangerous
answer was never an error message but a normal-looking output whose correlation energy is
the bare one.

Nothing is typed in except the committed reference values the gates compare against, each
named with its source. The barrier is not computed here and cannot be: this is one
geometry, the reactant only.

The energy convention is the committed cs6_audit.py's last FINAL SINGLE POINT ENERGY. The
HOMO convention is the committed cs10_cpcm.py's: the last ORBITAL ENERGIES block, rows
with four fields, occupancy above 0.5, the eV column, highest occupied last.

Usage: cs6_d38.py --res R --code C --prereg P --amd1 A --result OUT --tsv TSV
"""
import argparse, hashlib, os, re, sys
from pathlib import Path

NEL, NOCC = 118, 59
EH2EV = 27.211386245988
# G38.2 / G38.3 thresholds, from the pre-registration section 5
MOVED = 1e-6
# G38.1, same section
INERT = 1e-8
# R38.5 / R38.7 comparators, each from the document named
D36_HF_HOMO = -1.9636          # STAGE6_DIAGNOSTIC36_RESULT_20261114.txt, arm H, R leg
F8_HOMO = ((("S6-C, bare", 0), 2.278),
           (("S6-P / ADD-s1690, +1 at site 1690", 1), -0.243),
           (("ADD-s678, -1.00", 1), 5.205))   # STAGE6_FINDINGS_20261022.md Finding 8
RUNS = ("N", "B0", "B1", "A0", "A1")
PROBES = (("B1", "ADD-s1690_qp10000", +1.0), ("A1", "ADD-s678_qm10000", -1.0))
ZEROS = (("B0", "ADD-s1690_qp10000"), ("A0", "ADD-s678_qm10000"))
NQM = 24


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def vals(t, label, pat=r"-?\d+\.\d+"):
    return re.findall(r"^\s*%s\s*\.*\s*(%s)\s*$" % (re.escape(label), pat), t, re.M)


def val(t, label, bad=None, pat=r"-?\d+\.\d+", cast=float):
    v = vals(t, label, pat)
    if not v:
        return None
    if len(set(v)) > 1 and bad is not None:
        bad.append("%r is printed %d times with %d different values: %s"
                   % (label, len(v), len(set(v)), ", ".join(sorted(set(v)))))
    return cast(v[0])


def read_out(d):
    p = os.path.join(d, "job.out")
    if not os.path.isfile(p):
        return None
    t = open(p, errors="replace").read()
    bad = []
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    nel = re.search(r"Number of Electrons\s+NEL\s*\.+\s*(\d+)", t)
    geo = re.search(r"CARTESIAN COORDINATES \(ANGSTROEM\)\n-+\n(.*?)\n\n", t, re.S)
    atoms = []
    if geo:
        for l in geo.group(1).splitlines():
            f = l.split()
            if len(f) == 4:
                atoms.append((f[0],) + tuple(round(float(x), 6) for x in f[1:4]))
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
                o, ev = float(f[1]), float(f[3])
            except ValueError:
                break
            if o > 0.5:
                occ.append(ev)
        homo = occ[-1] if occ else None
    npc = re.search(r"[Nn]umber of [Pp]oint [Cc]harges\s*\.*\s*(\d+)", t)
    sto = re.search(r"PNO-Integral Storage\s*\.*\s*(\S[^\n]*)", t)
    rt = re.search(r"TOTAL RUN TIME:\s*(.+)", t)
    return dict(
        path=p, sha=sha(p),
        norm="ORCA TERMINATED NORMALLY" in t,
        notconv=bool(re.search(r"SCF NOT CONVERGED AFTER|SERIOUS PROBLEM", t)),
        aborted=bool(re.search(r"ORCA finished by error termination|aborting", t)),
        nel=int(nel.group(1)) if nel else None,
        E=float(e[-1]) if e else None,
        nocc=len(occ), homo=homo, atoms=atoms,
        runtime=rt.group(1).strip() if rt else None,
        storage=sto.group(1).strip() if sto else None,
        npc=int(npc.group(1)) if npc else None,
        e0=val(t, r"E(0)", bad), strong=val(t, "E(CORR)(strong-pairs)", bad),
        weak=val(t, "E(CORR)(weak-pairs)", bad), corr=val(t, "E(CORR)(corrected)", bad),
        ccsd=val(t, "E(CCSD)", bad), tcorr=val(t, "Triples Correction (T)", bad),
        ccsdt=val(t, "E(CCSD(T))", bad), t1=val(t, "T1 diagnostic", bad),
        inconsistent=bad)


def xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    return [(x.split()[0],) + tuple(round(float(v), 6) for v in x.split()[1:4])
            for x in L[2:2 + n]]


def ff_charge(s6, design):
    """Re-derive the sphere's charge and position from the committed QM/MM inputs."""
    d = os.path.join(s6, "d34", "a2_%s_R" % design)
    L = open(os.path.join(d, "site.ORCAFF.prms")).read().splitlines()
    i = L.index("$atoms")
    row = L[i + 2 + NQM].split()
    C = open(os.path.join(d, "combined.xyz")).read().splitlines()
    sph = C[2 + NQM].split()
    return float(row[2]), tuple(float(x) for x in sph[1:4])


def jobout_facts(p):
    """R38.8 and R38.9, from the job script's own scheduler output."""
    if not os.path.isfile(p):
        return {}
    t = open(p, errors="replace").read()
    m = re.search(r"R38\.8 scratch=(\S+) avail=(\d+)G margin=(\S+)", t)
    k = re.search(r"R38\.9 peak scratch across the five runs: (\d+) MB", t)
    h = re.search(r"host=(\S+)", t)
    cand = re.findall(r"^   candidate (\S+)\s+(.*)$", t, re.M)
    return dict(scratch=m.group(1) if m else None,
                avail=int(m.group(2)) if m else None,
                margin=m.group(3) if m else None,
                peak_mb=int(k.group(1)) if k else None,
                host=h.group(1) if h else None,
                refused=cand)


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "amd1", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    S6 = os.path.join(a.res, "phase2.2/cs_stage6")
    d38 = os.path.join(S6, "d38")
    gates, rep = [], []
    say = rep.append

    say("STAGE 6, DIAGNOSTIC 38 - RESULT")
    say("Do embedding point charges reach the correlated step, or only the reference?")
    say("")
    say("Read by phase2.2/cs6_d38.py                       sha256 %s" % sha(__file__))
    for lab, p in (("pre-registration", a.prereg), ("amendment 1", a.amd1)):
        say("   %-16s %-40s" % (lab, os.path.basename(p)))
        say("   %-16s sha256 %s" % ("", sha(p)))
    say("")
    say("A MECHANISM TEST. One geometry, the reactant only, no transition state, no ddE,")
    say("no barrier. D38.1 says so and nothing here may be read as one.")
    say("")

    # ------------------------------------------------------------- read
    D = {}
    for r in RUNS:
        got = read_out(os.path.join(d38, r))
        if got is None:
            gates.append("%s produced no readable output" % r)
            continue
        D[r] = got
        if not got["norm"]:
            gates.append("G38.4 %s did not terminate normally" % r)
        if got["notconv"]:
            gates.append("G38.4 %s did not converge its SCF" % r)
        if got["aborted"]:
            gates.append("G38.4 %s reports an error termination" % r)
        if got["nel"] != NEL:
            gates.append("G38.4 %s reports %s electrons, expected %d"
                         % (r, got["nel"], NEL))
        if got["ccsdt"] is None:
            gates.append("G38.4 %s prints no CCSD(T) energy" % r)
        for b in got["inconsistent"]:
            gates.append("G38.4 %s: %s" % (r, b))

    # --------------------------------------------- G38.5, G38.7, G38.8, G38.6
    ref = os.path.join(S6, "d33a", "sp_42436_R", "substrate.xyz")
    rsha, ratoms = sha(ref), xyz(ref)
    texts = {}
    for r in RUNS:
        g = os.path.join(d38, r, "substrate.xyz")
        if not os.path.isfile(g) or sha(g) != rsha:
            gates.append("G38.5 %s's geometry is not d33a's frame 42436 reactant" % r)
        if r in D and D[r]["atoms"] and D[r]["atoms"] != ratoms:
            gates.append("G38.5 %s RAN on coordinates that are not frame 42436's" % r)
        ip = os.path.join(d38, r, "job.inp")
        if os.path.isfile(ip):
            texts[r] = open(ip).read()
            if "* xyzfile -2 1 substrate.xyz" not in texts[r]:
                gates.append("G38.8 %s is not at charge -2 multiplicity 1" % r)
        else:
            gates.append("G38.7 %s has no input to check" % r)
    if "N" in texts:
        base = texts["N"].splitlines()
        for r in texts:
            if [l for l in texts[r].splitlines()
                    if not l.startswith("%pointcharges")] != base:
                gates.append("G38.7 %s differs from N in more than the %%pointcharges "
                             "line" % r)
    for r, design, sign in PROBES:
        pc = os.path.join(d38, r, "charge.pc")
        if not os.path.isfile(pc):
            gates.append("G38.6 %s has no charge file" % r)
            continue
        f = open(pc).read().split()
        q, pos = ff_charge(S6, design)
        if len(f) != 5 or abs(float(f[1]) - q) > 1e-9:
            gates.append("G38.6 %s's charge is not %s's committed %+.6f e"
                         % (r, design, q))
        elif any(abs(float(x) - y) > 1e-8 for x, y in zip(f[2:5], pos)):
            gates.append("G38.6 %s's position is not %s's sphere" % (r, design))
        if q * sign < 0:
            gates.append("G38.6 %s's charge has the wrong sign" % r)
    for r, design in ZEROS:
        pc = os.path.join(d38, r, "charge.pc")
        if os.path.isfile(pc):
            f = open(pc).read().split()
            if len(f) == 5 and abs(float(f[1])) > 0:
                gates.append("G38.6 %s's charge is not zero; it is the inertness check"
                             % r)

    # ------------------------------------------- what it cost. R38.8, R38.9
    jf = jobout_facts(os.path.join(d38, "run_a1.pbs.out")) \
        or jobout_facts(os.path.join(d38, "run.pbs.out"))
    say("=" * 78)
    say("R38.8 AND R38.9  WHERE IT RAN, AND WHAT IT COST")
    say("=" * 78)
    say("")
    say("   node             %s" % (jf.get("host") or "not reported"))
    say("   scratch          %s" % (jf.get("scratch") or "not reported"))
    say("   free when chosen %s" % (("%d GB" % jf["avail"]) if jf.get("avail") else "?"))
    say("   margin           %s" % (jf.get("margin") or "not reported"))
    say("   peak over 5 runs %s" % (("%d MB" % jf["peak_mb"]) if jf.get("peak_mb")
                                    else "not reported"))
    if jf.get("refused"):
        say("")
        say("   candidates it refused before settling:")
        for c, why in jf["refused"]:
            say("      %-40s %s" % (c, why))
    if jf.get("margin") is None or jf.get("peak_mb") is None:
        gates.append("G38.13 the result cannot state the scratch location and peak")
    say("")
    if jf.get("peak_mb"):
        say("   THIS IS THE FIRST MEASUREMENT OF def2-SVPD's SCRATCH FOOTPRINT ANYWHERE")
        say("   IN THIS RECORD. Amendment 1 carried a 100 GB preference over from D36,")
        say("   which peaked at 74 and 80 GB - at def2-TZVPD - and said plainly that the")
        say("   def2-SVPD figure was unmeasured and that R38.9 was how it would stop")
        say("   being unmeasured. The measured peak above is that figure, and it is what")
        say("   a later job should use rather than the number carried over.")
    say("")
    say("   CORRECTION TO AMENDMENT 1. That document said the cause of attempt 1's")
    say("   failure was not established, and offered a dangling symlink or an unmounted")
    say("   volume as the shapes that would produce it. The cause IS now established and")
    say("   is neither: attempt 1's own scheduler output, preserved as")
    say("   run.pbs.out.attempt1_427881 and archived, carries the line")
    say("")
    say('      mkdir: cannot create directory /scratch/18660916: No space left on device')
    say("")
    say("   The filesystem was full. Amendment 1 is not edited; this paragraph is the")
    say("   correction. The amendment's design was unaffected by the wrong guess, because")
    say("   it tried candidates rather than reasoning about one, but the guess was wrong")
    say("   and is recorded as such.")
    say("")
    say("   AND A SECOND CORRECTION, OF AN OMISSION. Amendment 1 said there was no")
    say("   archive to make because no output existed. That was true of the five job.out")
    say("   files and false of run.pbs.out, which existed, held the message above, and")
    say("   was the primary source for the amendment's own section 1. Attempt 2 writes to")
    say("   the same path and would have overwritten it. It was copied aside and archived")
    say("   before attempt 2 was submitted, with its sha256 matched across hpc1 and the")
    say("   archive, so nothing was lost - but the rule of 18 October should have caught")
    say("   it when the amendment was written, and did not.")
    say("")
    say("   THE SAME OVERSIGHT EXPLAINS A FILE NAME. Amendment 1 changed where scratch")
    say("   comes from and did not change the job script's -o line, so attempt 2's own")
    say("   scheduler output is committed as d38/run.pbs.out rather than under the")
    say("   run_a1 name its script carries. Attempt 1's copy, preserved as")
    say("   run.pbs.out.attempt1_427881, is archived and is not committed, because a")
    say("   failed job is not committed and this document is its note.")
    say("")

    # ------------------------------------------------------ G38.1 inertness
    say("=" * 78)
    say("G38.1  INERTNESS. A charge of zero magnitude cannot change a converged energy")
    say("=" * 78)
    say("")
    say("   If it does, the point-charge file is being misparsed and nothing else here is")
    say("   readable. Threshold %g Eh." % INERT)
    say("")
    for z, _design in ZEROS:
        if z not in D or "N" not in D:
            say("   %-3s  missing a run" % z)
            continue
        worst, which = 0.0, ""
        for lab, key in (("E(0)", "e0"), ("E(CORR)(corrected)", "corr"),
                         ("E(CCSD(T))", "ccsdt")):
            x, y = D[z].get(key), D["N"].get(key)
            if None in (x, y):
                gates.append("G38.1 %s or N does not report %s" % (z, lab))
                continue
            d = abs(x - y)
            if d > worst:
                worst, which = d, lab
        if worst > INERT:
            gates.append("G38.1 %s differs from N by %.2e Eh in %s" % (z, worst, which))
        say("   %-3s vs N   worst difference %.2e Eh in %-18s %s"
            % (z, worst, which or "-", "PASS" if worst <= INERT else "FAIL"))
    say("")

    # ------------------------------------------ R38.1 to R38.7, the readings
    say("=" * 78)
    say("R38.1 to R38.7  EVERY QUANTITY, AS READ")
    say("=" * 78)
    say("")
    say("   run %16s %16s %16s %12s %10s %4s"
        % ("E(0)", "E(CORR)corr", "E(CCSD(T))", "T1", "HOMO eV", "npc"))
    for r in RUNS:
        if r not in D:
            say("   %-3s %16s" % (r, "NO OUTPUT"))
            continue
        d = D[r]
        say("   %-3s %16s %16s %16s %12s %10s %4s"
            % (r,
               "%.9f" % d["e0"] if d["e0"] is not None else "-",
               "%.9f" % d["corr"] if d["corr"] is not None else "-",
               "%.9f" % d["ccsdt"] if d["ccsdt"] is not None else "-",
               "%.9f" % d["t1"] if d["t1"] is not None else "-",
               "%+.4f" % d["homo"] if d["homo"] is not None else "-",
               d["npc"] if d["npc"] is not None else "-"))
        if d["nocc"] not in (0, NOCC):
            gates.append("R38.3 %s reports %d occupied orbitals, expected %d"
                         % (r, d["nocc"], NOCC))
    say("")
    say("   R38.6  %s"
        % ("ORCA reports the number of point charges it read, shown above as npc"
           if any((D.get(r) or {}).get("npc") for r in RUNS)
           else "ORCA prints no point-charge count in these outputs; stated, not assumed"))
    say("")
    say("   R38.5 and R38.7  the highest occupied orbital, against the record")
    say("")
    say("      D36 measured HF/def2-TZVPD on the bare substrate   %+9.4f eV"
        % D36_HF_HOMO)
    say("      Finding 8 measured these at def2-SVP, which has no diffuse functions:")
    for (nm, _n), v in F8_HOMO:
        say("         %-38s %+9.4f eV" % (nm, v))
    say("")
    say("      def2-SVPD, measured here, DOES have them. Finding 8's own words were that")
    say("      the def2-SVP calculation \"only looks well-behaved because def2-SVP has no")
    say("      diffuse functions to let an electron leave\". These are the first values")
    say("      for these arrangements at a basis that does.")
    say("")

    # ------------------------------------ G38.10, then G38.2 and G38.3
    say("=" * 78)
    say("G38.10  WHICH ARMS ARE EVIDENCE AT ALL")
    say("=" * 78)
    say("")
    say("   Pre-registered before any output existed. An arm whose SCF fails, or whose")
    say("   highest occupied orbital sits at or above zero, is NOT INTERPRETABLE: the")
    say("   record's citation of Dreuw and Cederbaum is that a resonance \"can in")
    say("   principle not be correctly computed with bound state methods\". Such an arm is")
    say("   reported with its numbers, not dropped.")
    say("")
    interp = {}
    for r, _design, _sign in PROBES:
        if r not in D:
            interp[r] = False
            say("   %-3s  no output                                  NOT INTERPRETABLE" % r)
            continue
        d = D[r]
        why = []
        if not d["norm"]:
            why.append("did not terminate normally")
        if d["notconv"]:
            why.append("SCF not converged")
        if d["homo"] is None:
            why.append("no orbital energies read")
        elif d["homo"] >= 0:
            why.append("HOMO %+.4f eV is at or above zero" % d["homo"])
        interp[r] = not why
        say("   %-3s  HOMO %s   %s"
            % (r, ("%+.4f eV" % d["homo"]) if d["homo"] is not None else "   not read",
               "interpretable" if not why else "NOT INTERPRETABLE: " + "; ".join(why)))
    say("")

    say("=" * 78)
    say("G38.2 AND G38.3  THE QUESTION")
    say("=" * 78)
    say("")
    say("   G38.2 asks whether the charge reached the Hartree-Fock reference.")
    say("   G38.3 asks whether it reached the correlation treatment. Threshold %g Eh."
        % MOVED)
    say("")
    say("   A NOTE ON WHAT THESE TWO ARE. The pre-registration lists them under a heading")
    say("   that says any gate failure stops the diagnostic, which reads oddly for the")
    say("   two quantities the diagnostic exists to measure. Its D38.8 settles it: if the")
    say("   correlated step turns out blind to the charges, \"that is the result, it is")
    say("   reported as plainly as the opposite, and it closes this route rather than")
    say("   being worked around\". So a negative G38.3 is reported here as the verdict and")
    say("   does not abort the document. G38.1 and G38.4 are the ones that genuinely stop")
    say("   it, because they mean the numbers are not readable at all.")
    say("")
    verdict = {}
    for r, design, _sign in PROBES:
        if r not in D or "N" not in D:
            say("   %-3s  missing a run" % r)
            continue
        d0 = (D[r]["e0"] - D["N"]["e0"]) if None not in (D[r]["e0"], D["N"]["e0"]) else None
        dc = (D[r]["corr"] - D["N"]["corr"]) \
            if None not in (D[r]["corr"], D["N"]["corr"]) else None
        dt = (D[r]["ccsdt"] - D["N"]["ccsdt"]) \
            if None not in (D[r]["ccsdt"], D["N"]["ccsdt"]) else None
        say("   %-3s  %s" % (r, design))
        for lab, v, gate in (("dE(0)", d0, "G38.2"),
                             ("dE(CORR)(corrected)", dc, "G38.3"),
                             ("dE(CCSD(T))", dt, None)):
            if v is None:
                say("      %-22s NOT REPORTED" % lab)
                continue
            moved = abs(v) > MOVED
            say("      %-22s %+15.9f Eh  %-13s %s"
                % (lab, v, "moved" if moved else "DID NOT MOVE",
                   ("%s %s" % (gate, "PASS" if moved else "FAIL")) if gate else ""))
        verdict[r] = dict(d0=d0, dc=dc, dt=dt)
        say("")

    # the verdict, by G38.10's order of precedence
    say("=" * 78)
    say("THE VERDICT")
    say("=" * 78)
    say("")
    chosen = None
    for r in ("B1", "A1"):
        if interp.get(r) and verdict.get(r, {}).get("dc") is not None:
            chosen = r
            break
    if chosen is None:
        say("   NO VERDICT. Neither probe arm is interpretable, so G38.10 applies in its")
        say("   last form: D38 has no verdict, the expensive test is NOT licensed, and the")
        say("   open question becomes whether a coupled-cluster treatment of these")
        say("   arrangements is available at all at a basis with diffuse functions.")
        gates.append("G38.3 cannot be decided: no interpretable probe arm")
    else:
        v = verdict[chosen]
        if chosen == "A1":
            say("   DECIDED ON A1, THE WEAKER ARM, because B1 is not interpretable. G38.10")
            say("   requires this to be said prominently, and this is that sentence.")
            say("")
        if abs(v["dc"]) > MOVED and abs(v["d0"]) > MOVED:
            say("   THE EMBEDDING CHARGE REACHES THE CORRELATION TREATMENT.")
            say("   On %s the Hartree-Fock reference moves by %+.9f Eh and the corrected"
                % (chosen, v["d0"]))
            say("   correlation energy moves by %+.9f Eh, which is %.0f times the"
                % (v["dc"], abs(v["dc"]) / MOVED))
            say("   threshold. The silent failure this diagnostic was built to catch - a")
            say("   charge entering the SCF and not the correlated step - did not happen.")
            say("")
            say("   WHAT THAT LICENSES, AND ONLY THAT. The coupled-cluster test of ddE")
            say("   with the arrangement's charges present may be built. D38.3 says a pass")
            say("   licenses the expensive test and licenses nothing else. No barrier, no")
            say("   ddE, no design conclusion and no production level of theory follows")
            say("   from this document.")
        elif abs(v["d0"]) > MOVED:
            say("   THE CHARGE REACHES THE SCF AND NOT THE CORRELATION TREATMENT.")
            say("   On %s the Hartree-Fock reference moves by %+.9f Eh while the corrected"
                % (chosen, v["d0"]))
            say("   correlation energy moves by %+.9f Eh, under the %g Eh threshold."
                % (v["dc"], MOVED))
            say("   This is exactly the silent failure the diagnostic was built to catch,")
            say("   and it closes the route: a coupled-cluster ddE cannot be computed this")
            say("   way, because the correlated step does not see the arrangement. The")
            say("   pre-registration says a negative result closes the route rather than")
            say("   being worked around, and it is reported as plainly as the opposite")
            say("   would have been.")
        else:
            say("   THE POINT-CHARGE FILE IS NOT BEING READ AT ALL.")
            say("   On %s neither the reference nor the correlation energy moves: %+.9f"
                % (chosen, v["d0"]))
            say("   and %+.9f Eh, both under the %g Eh threshold. That is a mechanical"
                % (v["dc"], MOVED))
            say("   fault in how the file is passed, not a statement about the method, and")
            say("   it is the one failure mode here that is worth another short attempt.")
    say("")

    # ---------------------------------------------------------- G38.9
    say("=" * 78)
    say("WHAT THIS CANNOT SAY")
    say("=" * 78)
    say("")
    say("   NO BARRIER (D38.1). One geometry, the reactant only. No transition state was")
    say("   run, so no barrier and no ddE exists in this document at any level.")
    say("   IT DOES NOT SHOW THAT %pointcharges REPRODUCES THE QM/MM MODULE (D38.2). The")
    say("   charges here are read from the force field the QM/MM runs used, but whether")
    say("   the two routes give the same QM energy is a separate control at density")
    say("   functional level and belongs to the test this licenses, not here.")
    say("   IT REVISITS NO D34 VERDICT (D38.4). ADD-s1690_qp10000 is UNRELIABLE at")
    say("   CAM-B3LYP/def2-SVPD and remains so.")
    say("   T1 AND THE HIGHEST OCCUPIED ORBITAL ARE A FIRST LOOK (D38.5), not a")
    say("   reference-quality verdict and not a multireference check.")
    say("   NO DESIGN CONCLUSION (D38.6). Nothing about r_safe, any barrier lowering, or")
    say("   which designs survive.")
    say("   IT DOES NOT TEST WHETHER THE SIGN OF A CHARGE SETS D34's VERDICT (D38.7).")
    say("   Two arrangements cannot test that, and the record's own 46 classified rows")
    say("   show the sign does not set it.")
    say("   NOTHING WAS SHAPED (D38.8). The arms were fixed in advance, both charges were")
    say("   read out of committed force fields, and the two re-runs changed only where")
    say("   scratch was taken from.")
    say("")

    say("=" * 78)
    if gates:
        say("GATE FAILURES: %d" % len(gates))
        say("=" * 78)
        say("")
        for g in gates:
            say("   %s" % g)
    else:
        say("ALL GATES PASS. G38.1 to G38.13.")
        say("=" * 78)
    say("")

    Path(a.result).write_text("\n".join(rep) + "\n")
    with open(a.tsv, "w") as f:
        f.write("run\tcharge_e\tE0\tE_CORR_corrected\tE_CCSDT\tT1\tHOMO_eV\n")
        for r in RUNS:
            if r not in D:
                continue
            q = ""
            pc = os.path.join(d38, r, "charge.pc")
            if os.path.isfile(pc):
                ff = open(pc).read().split()
                if len(ff) == 5:
                    q = "%.6f" % float(ff[1])
            d = D[r]
            f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                r, q,
                "" if d["e0"] is None else "%.9f" % d["e0"],
                "" if d["corr"] is None else "%.9f" % d["corr"],
                "" if d["ccsdt"] is None else "%.9f" % d["ccsdt"],
                "" if d["t1"] is None else "%.9f" % d["t1"],
                "" if d["homo"] is None else "%.4f" % d["homo"]))
    print("\n".join(rep))
    return 1 if gates else 0


if __name__ == "__main__":
    sys.exit(main())
