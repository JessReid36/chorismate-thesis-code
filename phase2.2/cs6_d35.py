#!/usr/bin/env python3
"""DIAGNOSTIC 35: Stage 4e's five levels on frame 42436's stationary points.

These are bare single points with no charges and no MM region, so cs6_analyse.py's
bare / (MM) / (QM/MM) convention does not apply and is not used. The energy is the last
FINAL SINGLE POINT ENERGY, which is the convention the committed cs6_audit.py's own
reader uses, and the sha256 of every script is printed with the result.
"""
import argparse, hashlib, lzma, os, re, sys
from pathlib import Path

H2K = 627.5094740631
NEL = 118
CTL_R, CTL_TS = -836.1397955187, -836.1103760199
CTL_BARRIER = 18.461014
CTL_TOL_EH, CTL_TOL_KCAL = 1e-6, 0.0010
LEVELS = ("L1", "L2", "L3", "L4", "L5")
NAMES = {"L1": "B3LYP-D3BJ/def2-SVP", "L2": "B3LYP/6-31G(d)", "L3": "HF/6-31G(d)",
         "L4": "MP2/6-31G(d)", "L5": "wB97X-D3/def2-SVP"}


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


def read_out(d):
    """Everything this diagnostic needs from one output, with its own checks."""
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
    return dict(
        norm="ORCA TERMINATED NORMALLY" in t,
        notconv=bool(re.search(r"SCF NOT CONVERGED AFTER|SERIOUS PROBLEM", t)),
        nel=int(nel.group(1)) if nel else None,
        E=float(e[-1]) if e else None,
        atoms=atoms,
        method=(re.search(r"^\|\s*\d+> *(![^\n|]*)", t, re.M).group(1).strip()
                if re.search(r"^\|\s*\d+> *(!)", t, re.M) else None))


def xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    return [(x.split()[0],) + tuple(round(float(v), 6) for v in x.split()[1:4]) for x in L[2:2 + n]]


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    S6 = os.path.join(a.res, "phase2.2/cs_stage6")
    S4E = os.path.join(a.res, "phase2.2/cs_stage4e")
    d35 = os.path.join(S6, "d35")
    gates, rep = [], []
    say = rep.append

    say("STAGE 6, DIAGNOSTIC 35 - RESULT")
    say("Stage 4e's five levels of theory, on frame 42436's stationary points.")
    say("")
    say("Read by phase2.2/cs6_d35.py            sha256 %s" % sha(__file__))
    say("   pre-registration %-38s sha256 %s"
        % (os.path.basename(a.prereg), sha(a.prereg)))
    say("")
    say("These are bare single points: no charges, no MM region, no wall. The energy is the")
    say("last FINAL SINGLE POINT ENERGY, the convention the committed cs6_audit.py uses.")
    say("")

    # ---- read, with the gates that belong to reading
    E, A = {}, {}
    for k in LEVELS:
        for leg in ("R", "TS"):
            d = os.path.join(d35, "%s_%s" % (k, leg))
            r = read_out(d)
            if r is None:
                gates.append("%s_%s produced no output" % (k, leg))
                continue
            if not r["norm"]:
                gates.append("%s_%s did not terminate normally" % (k, leg))
            if r["notconv"]:
                gates.append("%s_%s did not converge its SCF" % (k, leg))
            if r["nel"] != NEL:
                gates.append("%s_%s reports %s electrons, expected %d" % (k, leg, r["nel"], NEL))
            if r["E"] is None:
                gates.append("%s_%s has no final energy" % (k, leg))
            E[(k, leg)] = r["E"]
            A[(k, leg)] = r

    # ---- G35.1, G35.2, G35.5: the inputs and the geometry, after the run
    for k in LEVELS:
        for leg in ("R", "TS"):
            d = os.path.join(d35, "%s_%s" % (k, leg))
            src = os.path.join(S4E, "%s_bare_%s" % (k, leg), "job.inp")
            if os.path.isfile(src) and os.path.isfile(os.path.join(d, "job.inp")):
                if sha(src) != sha(os.path.join(d, "job.inp")):
                    gates.append("G35.1 %s_%s's input is no longer Stage 4e's" % (k, leg))
            g = os.path.join(S6, "d33a", "sp_42436_%s" % leg, "substrate.xyz")
            if os.path.isfile(g) and os.path.isfile(os.path.join(d, "substrate.xyz")):
                if sha(g) != sha(os.path.join(d, "substrate.xyz")):
                    gates.append("G35.2 %s_%s's geometry is not frame 42436's" % (k, leg))
            want = xyz(g)
            got = A.get((k, leg), {}).get("atoms") or []
            if got and got != want:
                gates.append("G35.5 %s_%s ran on coordinates that are not frame 42436's" % (k, leg))
    s4e_geo = xyz(os.path.join(S4E, "L1_bare_R", "substrate.xyz"))
    for k in LEVELS:
        if A.get((k, "R"), {}).get("atoms") == s4e_geo:
            gates.append("G35.5 %s_R ran on frame 41786's geometry" % k)

    if gates:
        say("GATES")
        for g in gates:
            say("   FAILED  %s" % g)
        print("\n".join(rep))
        stop("%d gate(s) failed before any reading. Nothing is written." % len(gates))

    # ---- G35.3, the control
    say("=" * 78)
    say("G35.3  THE CONTROL. L1 is a level that already exists at this frame.")
    say("=" * 78)
    say("")
    dR = abs(E[("L1", "R")] - CTL_R)
    dT = abs(E[("L1", "TS")] - CTL_TS)
    bar = {k: (E[(k, "TS")] - E[(k, "R")]) * H2K for k in LEVELS}
    say("   reactant     %.10f Eh   committed %.10f   differ %.2e" % (E[("L1", "R")], CTL_R, dR))
    say("   transition   %.10f Eh   committed %.10f   differ %.2e" % (E[("L1", "TS")], CTL_TS, dT))
    say("   barrier      %.6f kcal/mol    committed %.6f    differ %.2e"
        % (bar["L1"], CTL_BARRIER, abs(bar["L1"] - CTL_BARRIER)))
    if max(dR, dT) > CTL_TOL_EH:
        gates.append("G35.3: L1's energies differ from the committed ones by up to %.2e Eh, "
                     "over the %.0e limit" % (max(dR, dT), CTL_TOL_EH))
    if abs(bar["L1"] - CTL_BARRIER) > CTL_TOL_KCAL:
        gates.append("G35.3: L1's barrier is %.6f, not %.6f within %.4f"
                     % (bar["L1"], CTL_BARRIER, CTL_TOL_KCAL))
    say("")
    if gates:
        say("   THE CONTROL FAILED. Nothing below is readable.")
        for g in gates:
            say("   FAILED  %s" % g)
        print("\n".join(rep))
        stop("the control failed.")
    say("   The control holds, so the pipeline and the geometry are what D35 says they are.")

    # ---- R35.1 and R35.2
    b41 = {}
    for k in LEVELS:
        r = read_out(os.path.join(S4E, "%s_bare_R" % k))
        t = read_out(os.path.join(S4E, "%s_bare_TS" % k))
        if r and t and r["E"] and t["E"]:
            b41[k] = (t["E"] - r["E"]) * H2K
        else:
            gates.append("Stage 4e's %s bare pair could not be read for the comparison" % k)
    say("")
    say("=" * 78)
    say("R35.1 and R35.2. THE FIVE LEVELS, AT BOTH FRAMES.")
    say("=" * 78)
    say("")
    say("   frame 41786 is Stage 4e's own geometry, which is NEB band image 0;")
    say("   frame 42436 is the SP package's optimised stationary points.")
    say("")
    say("   %-4s %-22s %12s %12s %12s" % ("", "level", "41786", "42436", "difference"))
    for k in LEVELS:
        say("   %-4s %-22s %12.4f %12.4f %+12.4f" % (k, NAMES[k], b41[k], bar[k], bar[k] - b41[k]))
    say("")
    sp = [bar[k] - b41[k] for k in LEVELS]
    say("   the geometry moves the barrier by %+.4f to %+.4f kcal/mol across the five levels"
        % (min(sp), max(sp)))

    # ---- R35.3
    o41 = sorted(LEVELS, key=lambda k: b41[k])
    o42 = sorted(LEVELS, key=lambda k: bar[k])
    say("")
    say("=" * 78)
    say("R35.3. DOES THE ORDER OF THE LEVELS SURVIVE THE GEOMETRY?")
    say("=" * 78)
    say("")
    say("   at 41786, lowest barrier first:  %s" % "  ".join(o41))
    say("   at 42436, lowest barrier first:  %s" % "  ".join(o42))
    say("")
    if o41 == o42:
        say("   PRESERVED. The five levels rank the same way at both geometries, so Stage 4e's")
        say("   method comparison is ordinally robust to the change of frame. That is a")
        say("   statement about the ORDER only; every absolute barrier moved.")
    else:
        say("   NOT PRESERVED. The ranking changes with the geometry, which is a finding about")
        say("   Stage 4e: a method comparison made at band image 0 does not carry over to the")
        say("   stationary points unchanged.")
        say("   moved: %s" % ", ".join("%s %d->%d" % (k, o41.index(k) + 1, o42.index(k) + 1)
                                       for k in LEVELS if o41.index(k) != o42.index(k)))

    # ---- R35.4
    say("")
    say("=" * 78)
    say("R35.4. EACH LEVEL AGAINST L4 (MP2/6-31G(d)) AT FRAME 42436.")
    say("=" * 78)
    say("")
    say("   THIS IS NOT STAGE 4e's k. Stage 4e's k = 0.82 to 0.89 is a ratio of the BARRIER")
    say("   CHANGES produced by eleven residues. Nothing here recomputes it.")
    say("")
    say("   %-4s %-22s %12s %12s" % ("", "level", "barrier", "minus L4"))
    for k in LEVELS:
        say("   %-4s %-22s %12.4f %+12.4f" % (k, NAMES[k], bar[k], bar[k] - bar["L4"]))
    say("")
    say("   MP2/6-31G(d) is Stage 4e's chosen reference at a small basis. It is not a")
    say("   benchmark: 6-31G(d) carries no diffuse functions, and this is a dianion.")

    # ---- R35.5
    say("")
    say("=" * 78)
    say("R35.5. WHAT MAY NOW BE READ AGAINST WHAT.")
    say("=" * 78)
    say("")
    say("   All five levels now exist at frame 42436 beside the ten already there, so the")
    say("   reference-barrier table can be read as one frame throughout. The frame 41786")
    say("   column stays in the record as Stage 4e's own measurement and is NOT deleted; it")
    say("   is simply not comparable with anything at 42436, and this result is the measure")
    say("   of how far apart they are.")
    say("")
    say("   What is still NOT settled by D35:")
    say("   - Stage 4e's conclusion itself. That rests on eleven residues at five levels,")
    say("     170 runs, whose coordinates are frame 41786's. D35 does not transfer it.")
    say("   - Which functional is right. D35 compares five levels against MP2 at a small")
    say("     basis on one geometry; it does not benchmark them.")
    say("   - Anything about charges, designs, reliability or r_safe.")

    say("")
    say("GATES")
    if gates:
        for g in gates:
            say("   FAILED  %s" % g)
    else:
        say("   all gates passed")

    Path(a.result).write_text("\n".join(rep) + "\n")
    with open(a.tsv, "w") as f:
        f.write("level\tmethod\tbarrier_41786\tbarrier_42436\tdifference\tminus_L4_42436\n")
        for k in LEVELS:
            f.write("%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\n"
                    % (k, NAMES[k], b41[k], bar[k], bar[k] - b41[k], bar[k] - bar["L4"]))
    print("\n".join(rep))
    if gates:
        stop("%d gate(s) failed; the result records them." % len(gates))


if __name__ == "__main__":
    main()
