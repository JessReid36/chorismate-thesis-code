#!/usr/bin/env python3
"""Stage 6 PROBE 3: are the committed designs electronically bound?

Probe 2 measured the vertical detachment energy with ONE designed charge present and found that a
-1 charge at 2.96 A drives the chorismate dianion to VDE -2.10 eV - metastable, Dreuw & Cederbaum's
resonance, which "can in principle not be correctly computed with bound state methods". Probe 2 was
pre-registered as NOT a clearance for the many-charge designs, and this is the measurement it
deferred.

Where the committed designs actually put their negative charge, from relax/ARRANGEMENTS.txt:

    S6-P    none            one +1 at 2.960 A
    S6-2    -1.00 e         at 3.485 A
    S6-20   -2.70 e         across eight sites between 3.388 and 3.485 A

If S6-2 or S6-20 is metastable, its ddE is not a barrier effect, and Diagnostic 10's sign reversal
and ranking inversion are artefacts of computing barriers for species that are not bound.

FOUR RUNS. The dianion single points already exist and are committed - Diagnostic 10's
cpcm/S6-*_R_vac - so only the monoanion partners are missing. No geometry moves. ma-def2-SVP, in
vacuo, no continuum.

  generate   build the four monoanion inputs from the committed D10 inputs
  measure    read them against D10's committed dianion energies, write PROBE3_RESULT.txt

usage: cs11_design_vde.py generate <results-repo> <out-dir>
       cs11_design_vde.py measure  <run-dir> <results-repo> <out-dir>
"""
import sys, os, re, hashlib
from pathlib import Path

EV = 27.211386245988
ARR = ("S6-C", "S6-P", "S6-2", "S6-20")
CPCM = "phase2.2/cs_stage6/cpcm"
GATE_EH = 1e-6
GATE_VDE_EV = 0.01
D1_CTL_R = -836.227913439
D2_VDE_C3 = 0.5123
NEAR = 3.5      # the stand-off Diagnostic 11 returned as r_safe


def qneg(res):
    """Negative charge each design places within r_safe of the substrate, DERIVED from the
    committed ARRANGEMENTS.txt and geometry rather than asserted here."""
    import math
    rd = Path(res) / "phase2.2/cs_stage6"
    L = (rd / "frame_42436/reactant_qm.xyz").read_text().splitlines()
    R = [(f[0], float(f[1]), float(f[2]), float(f[3]))
         for f in (l.split() for l in L[2:2 + int(L[0].split()[0])])]
    out, cur = {a: [0, 0.0] for a in ARR}, None
    for l in (rd / "relax/ARRANGEMENTS.txt").read_text().splitlines():
        m = re.match(r"\s{2}(S6-\S+)\s+(\d+) sphere", l)
        if m:
            cur = m.group(1)
            continue
        m = re.match(r"\s+site\s+\d+\s+q\s+([-+]\d+\.\d+)\s+at\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", l)
        if m and cur in out:
            q = float(m.group(1))
            pos = tuple(float(m.group(i)) for i in (2, 3, 4))
            if q < 0 and min(math.dist(pos, a[1:4]) for a in R) < NEAR:
                out[cur][0] += 1
                out[cur][1] += q
    return {a: tuple(v) for a, v in out.items()}


def stop(m):
    print("STOP: %s" % m, file=sys.stderr)
    raise SystemExit(1)


def src_dir(res, arr):
    d = Path(res) / CPCM / ("%s_R_vac" % arr)
    if not d.is_dir():
        stop("committed D10 input not found: %s" % d)
    return d


def generate(res, out):
    res, out = Path(res), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    man = []
    for arr in ARR:
        s = src_dir(res, arr)
        inp = (s / "job.inp").read_text()
        if "ma-def2-SVP" not in inp or "CPCM" in inp:
            stop("%s/job.inp is not the ma-def2-SVP vacuum input this probe expects" % s)
        d = out / ("%s_mono" % arr)
        d.mkdir(parents=True, exist_ok=True)
        # the geometry and the charges are copied byte for byte from the committed D10 input, so
        # the monoanion differs from its dianion partner in the charge line and nothing else
        (d / "substrate.xyz").write_bytes((s / "substrate.xyz").read_bytes())
        body = ["%maxcore 3000", "%pal nprocs 8 end", "%scf MaxIter 300 end"]
        if (s / "pc.pc").is_file():
            (d / "pc.pc").write_bytes((s / "pc.pc").read_bytes())
            body.append('%pointcharges "pc.pc"')
        (d / "job.inp").write_text(
            "! B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF\n"
            + "\n".join(body) + "\n* xyzfile -1 2 substrate.xyz\n")
        for f in sorted(os.listdir(d)):
            man.append("%s  %s/%s" % (hashlib.sha256((d / f).read_bytes()).hexdigest(), d.name, f))
    (out / "MANIFEST.sha256").write_text("\n".join(man) + "\n")
    QN = qneg(res)
    print("   generated %d monoanion inputs in %s" % (len(ARR), out))
    for arr in ARR:
        n, q = QN[arr]
        print("     %-6s_mono   negative charge within 3.5 A: %d site(s), %+.2f e" % (arr, n, q))


RE_E = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")


def energy(p, what):
    s = Path(p).read_text(errors="replace")
    if "****ORCA TERMINATED NORMALLY****" not in s:
        stop("%s (%s) did not terminate normally" % (p, what))
    e = RE_E.findall(s)
    if not e:
        stop("%s (%s) has no final energy" % (p, what))
    return float(e[-1])


def measure(run, res, out):
    run, res, out = Path(run), Path(res), Path(out)
    if str(out.resolve()).startswith(str(run.resolve())):
        stop("refusing to write inside the run directory being read: %s" % out)
    # GATE 0: the runs must be the committed inputs, byte for byte.
    cand = [m for m in Path(res).rglob("MANIFEST.sha256")
            if ".git" not in m.parts and (m.parent / "S6-C_mono").is_dir()]
    if len(cand) != 1:
        stop("expected exactly one MANIFEST.sha256 beside the Probe 3 runs, found %d" % len(cand))
    bad = []
    for l in cand[0].read_text().splitlines():
        if not l.strip():
            continue
        want, rel = l.split(None, 1)
        f = run / rel.strip()
        if not f.is_file():
            bad.append("%s is missing" % rel.strip())
        elif hashlib.sha256(f.read_bytes()).hexdigest() != want:
            bad.append("%s does not match the committed input" % rel.strip())
    if bad:
        stop("the runs are not the committed inputs:\n   " + "\n   ".join(bad))

    di = {a: energy(src_dir(res, a) / "job.out", "%s dianion, committed D10" % a) for a in ARR}
    mo = {a: energy(run / ("%s_mono" % a) / "job.out", "%s monoanion" % a) for a in ARR}
    vde = {a: (mo[a] - di[a]) * EV for a in ARR}

    O = ["STAGE 6 - PROBE 3 RESULT: are the committed designs electronically bound?",
         "",
         "VDE = E(monoanion) - E(dianion) at the committed reactant geometry, the design's charges",
         "present in both. Positive means bound. ma-def2-SVP, in vacuo, no continuum.",
         "",
         "The dianion energies are Diagnostic 10's committed cpcm/S6-*_R_vac outputs, not re-run.",
         "Only the monoanion partners are new, and they differ from their dianion partners in the",
         "charge line and nothing else - same substrate.xyz, same pc.pc, byte for byte.",
         "",
         "GATE 1: THE CONTROL DIANION MUST REPRODUCE DIAGNOSTIC 1'S ENERGY", ""]
    g1 = abs(di["S6-C"] - D1_CTL_R) <= GATE_EH
    O += ["   measured %.9f   D1 %.9f   diff %+.2e Eh   %s"
          % (di["S6-C"], D1_CTL_R, di["S6-C"] - D1_CTL_R, "OK" if g1 else "FAIL"),
          "", "   GATE 1 %s" % ("PASSED" if g1 else "FAILED"), ""]
    if not g1:
        return finish(O, out, "   The control is not D1's calculation; nothing below is interpretable.")

    O += ["GATE 2: THE CONTROL VDE MUST REPRODUCE DIAGNOSTIC 2'S COMMITTED VALUE", "",
          "   D2 committed +%.4f eV for the bare dianion at ma-def2-SVP (C3)." % D2_VDE_C3, ""]
    g2 = abs(vde["S6-C"] - D2_VDE_C3) <= GATE_VDE_EV
    O += ["   measured %+.4f eV   D2 %+.4f eV   diff %+.4f eV   %s"
          % (vde["S6-C"], D2_VDE_C3, vde["S6-C"] - D2_VDE_C3, "OK" if g2 else "FAIL"),
          "", "   GATE 2 %s" % ("PASSED" if g2 else "FAILED"), ""]
    if not g2:
        return finish(O, out, "   This probe does not reproduce Diagnostic 2 on the bare substrate,\n"
                              "   so its designs cannot be compared with it. Nothing further.")

    QN = qneg(res)
    O += ["THE READING", "",
          "   design   neg. charge within %.1f A    VDE / eV    bound?" % NEAR,
          "            sites      sum q-"]
    for a in ARR:
        n, q = QN[a]
        O.append("   %-7s  %4d    %8.2f e       %+9.4f    %s"
                 % (a, n, q, vde[a], "bound" if vde[a] > 0 else "METASTABLE"))
    meta = [a for a in ARR if a != "S6-C" and vde[a] <= 0]
    O += ["", "VERDICT", ""]
    if not meta:
        O += ["   Every committed design is bound at the reactant geometry. Dreuw & Cederbaum's",
              "   resonance is not realised by any of them, and their ddE values are not",
              "   disqualified on this ground. Diagnostic 11's separate finding stands unaffected:",
              "   all three sit inside r_safe and their ddE values remain uninterpretable as",
              "   design evaluations for THAT reason."]
    else:
        O += ["   METASTABLE: %s." % ", ".join(meta),
              "",
              "   These designs are not bound states at the reactant geometry. A barrier computed",
              "   for them is not a barrier: it is a bound-state calculation of a resonance, which",
              "   Dreuw & Cederbaum say \"can in principle not be correctly computed\" that way.",
              "   Their Diagnostic 1 and Diagnostic 10 ddE values are WITHDRAWN as barrier effects,",
              "   and the sign reversal and ranking inversion built on them are withdrawn with",
              "   them - not as wrong measurements, but as measurements of something else.",
              "",
              "   This implies a constraint the optimiser does not have. Its budget is sum q_i^2 <= B,",
              "   which is BLIND TO SIGN. A bound on negative charge near the dianion is a physical",
              "   constraint, and the grid's admitted negative ranges need one before any further",
              "   design work."]
    O += ["", "WHAT THIS DOES NOT TEST", "",
          "   The transition-state geometry - this and D2's VDE are both at the reactant, so",
          "   neither says whether boundedness changes along the reaction coordinate. The basis",
          "   dependence Dreuw & Cederbaum warn of, since this runs one basis. And it does not",
          "   locate the threshold: it says which of four arrangements are bound, not how much",
          "   negative charge is too much.", ""]
    return finish(O, out, None)


def finish(O, out, tail):
    if tail:
        O += [tail, ""]
    out.mkdir(parents=True, exist_ok=True)
    (out / "PROBE3_RESULT.txt").write_text("\n".join(O) + "\n")
    print("\n".join(O))


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "generate":
        generate(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 5 and sys.argv[1] == "measure":
        measure(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        stop("usage: cs11_design_vde.py generate <results-repo> <out-dir>\n"
             "       cs11_design_vde.py measure  <run-dir> <results-repo> <out-dir>")
