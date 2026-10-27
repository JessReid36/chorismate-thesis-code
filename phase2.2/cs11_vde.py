#!/usr/bin/env python3
"""Stage 6 PROBE 2: does a designed charge push the dianion from bound to metastable?

Diagnostic 2 measured the vertical detachment energy of the BARE chorismate dianion at the reactant
geometry and committed it on 22 October:

    C1  def2-SVP      VDE  -0.3167 eV   metastable
    C3  ma-def2-SVP   VDE  +0.5123 eV   BOUND
    C5  def2-SVPD     VDE  +0.5013 eV   BOUND

So at the basis Stage 6 actually runs in, the bare substrate is bound. What is NOT known is whether
it stays bound once a designed charge is present. Dreuw & Cederbaum (cr0104227, held) give the
reason to ask: "at short distances the monoanion-point charge system represents a resonance state
that can in principle not be correctly computed with bound state methods". A negative charge near a
dianion is the dangerous case, and S6-2 and S6-20 carry them.

VDE = E(monoanion, charge -1, doublet) - E(dianion, charge -2, singlet) at the SAME geometry with
the SAME external charge present. Positive means bound. The method copies the committed, already-run
phase2.2/cs_stage6/qscan/VDE/C3_mono_R/job.inp exactly.

THE HOMO IS NOT THE TEST. Probe 1 used it and was wrong to; see STAGE6_PROBE1_ERRATUM1_20261027.
Diagnostic 2's own words: "A positive VDE means the species is bound despite a positive HOMO, and
Finding 8's second half must be withdrawn as an artefact of reading Koopmans' theorem off a DFT
orbital." HOMO values are printed below as data, never as a verdict.

Eight single points at the reactant geometry only, in vacuo, no continuum. Self-contained: each
condition carries its own dianion and monoanion run, so this does not wait on Diagnostic 11.

  generate   build the eight input directories from committed files only
  measure    read them and write PROBE2_RESULT.txt

usage: cs11_vde.py generate <results-repo> <out-dir>
       cs11_vde.py measure  <run-dir> <results-repo> <out-dir>
"""
import sys, os, re, math, hashlib
from pathlib import Path

EV = 27.211386245988
GATE_EH = 1e-6                 # the control dianion must reproduce D1's own energy
GATE_VDE_EV = 0.01             # the control VDE must reproduce D2's committed +0.5123 eV
D1_CTL_R = -836.227913439
D2_VDE_C3 = 0.5123
FRAME_DIR = "phase2.2/cs_stage6/frame_42436"
FRAME_SHA = {"reactant_qm.xyz":
             "6e8f41c7da3bcf4ed9f5c650dc9fcf12c2a36e5af7a900e164892c276497aa88"}

# condition -> (label on the A1 axis or None, charge)
COND = {"ctl":      (None,   0.0),
        "dref_p1":  ("dref", +1.0),
        "dref_m1":  ("dref", -1.0),
        "d250_m1":  ("d250", -1.0)}
SPECIES = {"di": (-2, 1), "mono": (-1, 2)}


def stop(m):
    print("STOP: %s" % m, file=sys.stderr)
    raise SystemExit(1)


def one(root, name):
    h = [p for p in Path(root).rglob(name) if ".git" not in p.parts]
    if len(h) != 1:
        stop("expected exactly one %s under %s, found %d" % (name, root, len(h)))
    return h[0]


def frame(res):
    d = Path(res) / FRAME_DIR
    for f, want in FRAME_SHA.items():
        p = d / f
        if not p.is_file():
            stop("not found: %s" % p)
        got = hashlib.sha256(p.read_bytes()).hexdigest()
        if got != want:
            stop("%s sha256 is %s, ARRANGEMENTS.txt records %s" % (p, got, want))
    return d


def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    return [(f[0], float(f[1]), float(f[2]), float(f[3]))
            for f in (l.split() for l in L[2:2 + n])]


def read_xyz_raw(p):
    L = Path(p).read_text().splitlines()
    return L[2:2 + int(L[0].split()[0])]


def site1690(res):
    for l in Path(one(res, "ARRANGEMENTS.txt")).read_text().splitlines():
        m = re.match(r"\s+site\s+1690\s+q\s+[-+]\d+\.\d+\s+at\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", l)
        if m:
            return tuple(float(m.group(i)) for i in (1, 2, 3))
    stop("site 1690 is not in ARRANGEMENTS.txt")


def position(res, label):
    """The same axis Diagnostic 11 uses, so the two are directly comparable."""
    A = site1690(res)
    if label == "dref":
        return A
    R = read_xyz(frame(res) / "reactant_qm.xyz")
    k, d0 = min(((i, math.dist(A, a[1:4])) for i, a in enumerate(R)), key=lambda t: t[1])
    o = R[k][1:4]
    u = [(A[j] - o[j]) / d0 for j in range(3)]
    d = float(label[1:]) / 100.0
    return tuple(o[j] + u[j] * d for j in range(3))


def generate(res, out):
    res, out = Path(res), Path(out)
    raw = read_xyz_raw(frame(res) / "reactant_qm.xyz")
    if len(raw) != 24:
        stop("reactant_qm.xyz has %d atoms, expected 24" % len(raw))
    out.mkdir(parents=True, exist_ok=True)
    man = []
    for cname, (label, q) in sorted(COND.items()):
        for sp, (chg, mult) in SPECIES.items():
            d = out / ("%s_%s" % (cname, sp))
            d.mkdir(parents=True, exist_ok=True)
            (d / "substrate.xyz").write_text("24\n%s %s, R, from %s\n%s\n"
                                             % (cname, sp, frame(res), "\n".join(raw)))
            body = ["%maxcore 3000", "%pal nprocs 8 end", "%scf MaxIter 300 end"]
            if label is not None:
                p = position(res, label)
                (d / "pc.pc").write_text("1\n%12.6f %18.8f %18.8f %18.8f\n" % (q, *p))
                body.append('%pointcharges "pc.pc"')
            (d / "job.inp").write_text(
                "! B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF\n"
                + "\n".join(body) + "\n* xyzfile %d %d substrate.xyz\n" % (chg, mult))
            for f in sorted(os.listdir(d)):
                man.append("%s  %s/%s" % (hashlib.sha256((d / f).read_bytes()).hexdigest(),
                                          d.name, f))
    (out / "MANIFEST.sha256").write_text("\n".join(man) + "\n")
    print("   generated %d input directories in %s" % (len(COND) * 2, out))
    for cname, (label, q) in sorted(COND.items()):
        if label is None:
            print("     %-9s no charge" % cname)
        else:
            p = position(res, label)
            R = read_xyz(frame(res) / "reactant_qm.xyz")
            print("     %-9s q %+.1f at %s, %.3f A from the nearest substrate atom"
                  % (cname, q, label, min(math.dist(p, a[1:4]) for a in R)))


RE_E = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")


def read_out(p):
    s = Path(p).read_text(errors="replace")
    e = RE_E.findall(s)
    k, occ = s.rfind("ORBITAL ENERGIES"), []
    homo = None
    if k >= 0:
        for l in s[k:].splitlines()[4:]:
            f = l.split()
            if len(f) != 4:
                if occ:
                    break
                continue
            try:
                o_, _, ev = float(f[1]), float(f[2]), float(f[3])
            except ValueError:
                break
            if o_ > 0.5:
                occ.append(ev)
        homo = occ[-1] if occ else None
    cy = re.findall(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES", s)
    return dict(E=float(e[-1]) if e else None, homo=homo,
                cyc=int(cy[-1]) if cy else None,
                norm="****ORCA TERMINATED NORMALLY****" in s)


def measure(run, res, out):
    run, res, out = Path(run), Path(res), Path(out)
    if str(out.resolve()) == str(run.resolve()) or \
       str(out.resolve()).startswith(str(run.resolve()) + os.sep):
        stop("refusing to write inside the run directory being read: %s" % out)
    cand = [m for m in Path(res).rglob("MANIFEST.sha256")
            if ".git" not in m.parts and (m.parent / "ctl_di").is_dir()]
    if len(cand) != 1:
        stop("expected exactly one MANIFEST.sha256 beside the Probe 2 runs, found %d" % len(cand))
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
    nman = len([l for l in cand[0].read_text().splitlines() if l.strip()])

    D = {}
    for cname in COND:
        for sp in SPECIES:
            p = run / ("%s_%s" % (cname, sp)) / "job.out"
            if not p.is_file():
                stop("missing output: %s" % p)
            d = read_out(p)
            if not d["norm"] or d["E"] is None:
                stop("%s did not terminate normally or has no final energy" % p)
            D[(cname, sp)] = d

    def vde(c):
        return (D[(c, "mono")]["E"] - D[(c, "di")]["E"]) * EV

    O = ["STAGE 6 - PROBE 2 RESULT: does a designed charge unbind the dianion?",
         "",
         "VDE = E(monoanion) - E(dianion) at the committed reactant geometry, same external charge",
         "in both. Positive means bound. ma-def2-SVP, in vacuo, no continuum.",
         "",
         "THE HOMO IS NOT THE TEST - Diagnostic 2 withdrew a finding for using it, and Probe 1",
         "repeated the error (STAGE6_PROBE1_ERRATUM1_20261027). HOMO is printed as data only.",
         "",
         "GATE 0: THE RUNS MUST BE THE COMMITTED INPUTS",
         "", "   %d files checked against MANIFEST.sha256, all match." % nman,
         "", "   GATE 0 PASSED", "",
         "GATE 1: THE CONTROL DIANION MUST REPRODUCE DIAGNOSTIC 1'S ENERGY", ""]
    got = D[("ctl", "di")]["E"]
    g1 = abs(got - D1_CTL_R) <= GATE_EH
    O += ["   measured %.9f   D1 %.9f   diff %+.2e Eh   %s"
          % (got, D1_CTL_R, got - D1_CTL_R, "OK" if g1 else "FAIL"),
          "", "   GATE 1 %s" % ("PASSED" if g1 else "FAILED"), ""]
    if not g1:
        return finish(O, out, "   The control is not D1's calculation; nothing below is interpretable.")

    O += ["GATE 2: THE CONTROL VDE MUST REPRODUCE DIAGNOSTIC 2'S COMMITTED VALUE", "",
          "   D2 committed +%.4f eV for the bare dianion at ma-def2-SVP (C3)." % D2_VDE_C3, ""]
    v0 = vde("ctl")
    g2 = abs(v0 - D2_VDE_C3) <= GATE_VDE_EV
    O += ["   measured %+.4f eV   D2 %+.4f eV   diff %+.4f eV   %s"
          % (v0, D2_VDE_C3, v0 - D2_VDE_C3, "OK" if g2 else "FAIL"),
          "", "   GATE 2 %s" % ("PASSED" if g2 else "FAILED"), ""]
    if not g2:
        return finish(O, out,
                      "   This probe does not reproduce Diagnostic 2 on the bare substrate, so its\n"
                      "   charged conditions cannot be compared with it. Nothing further is reported.")

    O += ["THE READING", "",
          "   condition   q      d / A     VDE / eV    bound?      HOMO(di) / eV   SCF cyc"]
    R = read_xyz(frame(res) / "reactant_qm.xyz")
    verdicts = {}
    for cname, (label, q) in sorted(COND.items()):
        v = vde(cname)
        dd = ("%8.3f" % min(math.dist(position(res, label), a[1:4]) for a in R)) if label else "       -"
        h = D[(cname, "di")]["homo"]
        verdicts[cname] = v
        O.append("   %-10s %+4.1f %s   %+9.4f    %-10s %13s   %7d"
                 % (cname, q, dd, v, "bound" if v > 0 else "METASTABLE",
                    ("%+.4f" % h) if h is not None else "-", D[(cname, "di")]["cyc"]))
    worst = min(verdicts, key=lambda c: verdicts[c])
    O += ["", "VERDICT", ""]
    if all(v > 0 for v in verdicts.values()):
        O += ["   No condition tested drives the dianion metastable. The lowest is %s at %+.4f eV."
              % (worst, verdicts[worst]),
              "   Dreuw & Cederbaum's resonance hazard is not realised at these charges and",
              "   separations, so a bound-state treatment remains defensible for them.",
              "   This says nothing about the many-charge designs: S6-20 carries ten negative",
              "   charges totalling -3.00 e and is not tested here."]
    else:
        bad = [c for c in verdicts if verdicts[c] <= 0]
        O += ["   METASTABLE at: %s." % ", ".join(sorted(bad)),
              "   A designed charge drives the dianion to or past the detachment threshold, which is",
              "   Dreuw & Cederbaum's case: the system is then a resonance that \"can in principle not",
              "   be correctly computed with bound state methods\". Any ddE measured at such a",
              "   condition is not a barrier effect and must not be reported as one."]
    O += ["", "WHAT THIS DOES NOT TEST", "",
          "   The transition-state geometry - D2's VDE and this one are both at the reactant.",
          "   Multi-charge arrangements. The relaxed geometries. And it does not revisit any",
          "   Diagnostic 1, 2 or 10 number.", ""]
    return finish(O, out, None)


def finish(O, out, tail):
    if tail:
        O += [tail, ""]
    out.mkdir(parents=True, exist_ok=True)
    (out / "PROBE2_RESULT.txt").write_text("\n".join(O) + "\n")
    print("\n".join(O))


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "generate":
        generate(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 5 and sys.argv[1] == "measure":
        measure(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        stop("usage: cs11_vde.py generate <results-repo> <out-dir>\n"
             "       cs11_vde.py measure  <run-dir> <results-repo> <out-dir>")
