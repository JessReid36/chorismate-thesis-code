#!/usr/bin/env python3
"""Stage 6 DIAGNOSTIC 11: how far must a designed charge stand off before a bare point charge is a
faithful model of it?

Senn & Thiel (anie.200802019, held in chorismate-refs) state the hazard in terms this project's
designs meet on both axes at once: the overpolarisation artefact "becomes more severe the closer the
point charge approaches the QM density and the more spatially flexible the density is", and "Larger
basis sets, which include polarization and possibly diffuse functions, are more flexible and hence
more prone to overpolarization." Stage 6's charges sit 2.85-2.96 A from the substrate and are run in
ma-def2-SVP.

Everything downstream of that - the a coefficients, the grid's admitted sites, the search itself -
assumes a point charge at ~2.9 A is a valid stand-in for a catalytic group. This diagnostic measures
whether it is, by walking ONE charge outward along a fixed axis and comparing a bare point charge
against the same charge behind a Pauli wall.

NO GEOMETRY MOVES. Single points on Diagnostic 1's committed R and TS.

  generate   build the input directories from committed files only
  measure    read the outputs and write DIAGNOSTIC11_RESULT.txt
  stage      copy the outputs and the result into the results repo

usage: cs11_standoff.py generate <results-repo> <out-dir>
       cs11_standoff.py measure  <run-dir> <results-repo> <out-dir>
       cs11_standoff.py stage    <run-dir> <results-repo> <out-dir>
"""
import sys, os, re, math, hashlib
from pathlib import Path

HARTREE = 627.5094740631
GEOM = ("R", "TS")
LADDER = ("d250", "dref", "d350", "d450", "d600")
GATE1_EH = 1e-6          # the ma control must reproduce D1's own energies
GATE2_KCAL = 0.01        # the committed-position bare arm must reproduce D10's ddE
WALL0_KCAL = 0.30        # C-A1.1: a wall carrying no charge must do nothing
SHIELD_REL = 0.10        # r_safe (i): |dshield| <= 10% of |ddE(bare)| ...
SHIELD_ABS = 0.05        #            ... or this absolute floor, whichever is larger
PLATEAU_REL = 0.10       # r_safe (iii): measured/classical changes <= 10% between successive d
D10_BARE_296 = -3.857852 # D10's S6-P vacuum ddE, which d=2.96 bare/ma must reproduce
BOHR = 0.529177210903

# The Ne coreless-ECP block, copied byte for byte from the committed, already-run input
# phase2.2/cs_stage4c/L0_E5_R/job.inp. ORCA calls such a centre a "coreless ECP center with
# (optional) point charge". STAGE6_ERRATUM7 records that these are Marefat Khah's published Ne
# parameters used outside their octet assumption; that caveat is carried, not resolved, here.
ECP_TAIL = """  N_core 10
  lmax f
  s 1
    1  2.0475   54.5100  2
  p 1
    1  0.44815   1.46500  2
  d 1
    1  0.49205  -0.83800  2
  f 1
    1  1.0000    0.0000  2
end"""

BASIS = {"ma": "ma-def2-SVP", "svp": "def2-SVP"}


def stop(m):
    print("STOP: %s" % m, file=sys.stderr)
    raise SystemExit(1)


def one(root, name):
    h = [p for p in Path(root).rglob(name) if ".git" not in p.parts]
    if len(h) != 1:
        stop("expected exactly one %s under %s, found %d" % (name, root, len(h)))
    return h[0]


FRAME_DIR = "phase2.2/cs_stage6/frame_42436"
FRAME_SHA = {"reactant_qm.xyz":
             "6e8f41c7da3bcf4ed9f5c650dc9fcf12c2a36e5af7a900e164892c276497aa88",
             "transition_state_qm.xyz":
             "2babeb9000ca5220a908a70eaec80ec4a64358c8f5677e9a16cb104e0dffc5bf"}


def frame(res):
    """Stage 6's geometries, by their committed path AND their committed hashes.

    Five directories named frame_42436 exist in the results repo, so the name alone will not do.
    STAGE6_CRITERIA.txt names this one, and relax/ARRANGEMENTS.txt records both sha256 values in
    its INPUTS block; both are checked here rather than trusted.
    """
    d = Path(res) / FRAME_DIR
    for f, want in FRAME_SHA.items():
        p = d / f
        if not p.is_file():
            stop("not found: %s" % p)
        got = hashlib.sha256(p.read_bytes()).hexdigest()
        if got != want:
            stop("%s sha256 is %s, ARRANGEMENTS.txt records %s" % (p, got, want))
    return d


def read_xyz_raw(p):
    """The committed coordinate lines, verbatim.

    Diagnostic 10 wrote its substrate.xyz by copying these lines through unchanged, so carrying
    them rather than reformatting makes the dref rung byte-identical to S6-P and gate 2 an exact
    reproduction instead of a near one.
    """
    L = Path(p).read_text().splitlines()
    return L[2:2 + int(L[0].split()[0])]


def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    return [(f[0], float(f[1]), float(f[2]), float(f[3]))
            for f in (l.split() for l in L[2:2 + n])]


def axis(res):
    """The A1 axis, derived from committed files alone.

    The charge walks outward from the substrate atom nearest to site 1690 at R, along the unit
    vector pointing from that atom to site 1690. At d = 2.9596 A the position IS site 1690, which
    is S6-P's committed position, so the ladder passes exactly through the arrangement already
    measured by Diagnostic 10. The same positions are used at R and at TS, because a designed
    charge does not move with the substrate.
    """
    fr = frame(res)
    R = read_xyz(fr / "reactant_qm.xyz")
    A = read_arrangement_site(one(res, "ARRANGEMENTS.txt"))
    k, d0 = min(((i, math.dist(A, a[1:4])) for i, a in enumerate(R)), key=lambda t: t[1])
    o = R[k][1:4]
    u = tuple((A[j] - o[j]) / d0 for j in range(3))
    return o, u, k, d0, R[k][0]


def read_arrangement_site(p):
    """Site 1690's committed position, from ARRANGEMENTS.txt - the authority."""
    for l in Path(p).read_text().splitlines():
        m = re.match(r"\s+site\s+1690\s+q\s+[-+]\d+\.\d+\s+at\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", l)
        if m:
            return tuple(float(m.group(i)) for i in (1, 2, 3))
    stop("site 1690 is not in %s" % p)


# name -> (d_centi_angstrom, arm, q, basis).  arm: none | bare | wall | wall0
def runset():
    R = {}
    for b in ("ma", "svp"):
        R["ctl_%s" % b] = (None, "none", 0.0, b)
    for d in LADDER:                                          # main ladder, ma, q = +1
        for arm in ("bare", "wall"):
            R["%s_%s_p1_ma" % (d, arm)] = (d, arm, +1.0, "ma")
        R["%s_wall0_ma" % d] = (d, "wall0", 0.0, "ma")        # C-A1.1 wall inertness
    for d in ("dref", "d450"):                                # C-A1.4 sign control
        for arm in ("bare", "wall"):
            R["%s_%s_m1_ma" % (d, arm)] = (d, arm, -1.0, "ma")
    for d in ("dref", "d450"):                                # C-A1.5 basis control
        for arm in ("bare", "wall"):
            R["%s_%s_p1_svp" % (d, arm)] = (d, arm, +1.0, "svp")
    R["d1500_bare_p1_ma"] = ("d1500", "bare", +1.0, "ma")     # A3.4 far-field zero
    return R


def position(res, label):
    """Where the charge sits for a ladder rung.

    'dref' is placed at site 1690's COMMITTED coordinates, exactly, so that the rung is byte-for-byte
    the arrangement Diagnostic 10 measured as S6-P and gate 2 is a true reproduction rather than a
    near miss. Every other rung is o + u*d along the same axis.
    """
    A = read_arrangement_site(one(res, "ARRANGEMENTS.txt"))
    if label == "dref":
        return A
    o, u, _, _, _ = axis(res)
    d = float(label[1:]) / 100.0
    return tuple(o[j] + u[j] * d for j in range(3))


def generate(res, out):
    res, out = Path(res), Path(out)
    fr = frame(res)
    G = {"R": read_xyz(fr / "reactant_qm.xyz"),
         "TS": read_xyz(fr / "transition_state_qm.xyz")}
    RAW = {"R": read_xyz_raw(fr / "reactant_qm.xyz"),
           "TS": read_xyz_raw(fr / "transition_state_qm.xyz")}
    for g, a in G.items():
        if len(a) != 24:
            stop("%s geometry has %d atoms, expected 24" % (g, len(a)))
    o, u, k, d0, el = axis(res)
    if abs(d0 - 2.9596) > 1e-3:
        stop("the nearest substrate atom to site 1690 is at %.4f A, expected 2.9596" % d0)
    out.mkdir(parents=True, exist_ok=True)
    man, rows = [], []
    for name, (dc, arm, q, b) in sorted(runset().items()):
        for g in GEOM:
            d = out / ("%s_%s" % (name, g))
            d.mkdir(parents=True, exist_ok=True)
            atoms = G[g]
            (d / "substrate.xyz").write_text(
                "24\n%s %s, from %s\n%s\n" % (name, g, fr, "\n".join(RAW[g])))
            head = "! B3LYP D3BJ %s def2/J RIJCOSX TightSCF" % BASIS[b]
            body = ["%maxcore 3000", "%pal nprocs 8 end", "%scf MaxIter 250 end"]
            if arm == "none":
                coord = "* xyzfile -2 1 substrate.xyz\n"
            else:
                p = position(res, dc)
                mind = min(math.dist(p, a[1:4]) for a in atoms)
                rows.append((name, g, math.dist(o, p), mind))
                if arm == "bare":
                    (d / "pc.pc").write_text("1\n%12.6f %18.8f %18.8f %18.8f\n" % (q, *p))
                    body.append('%pointcharges "pc.pc"')
                    coord = "* xyzfile -2 1 substrate.xyz\n"
                else:                                   # wall or wall0: coordinate-block ECP centre
                    qq = 0.0 if arm == "wall0" else q
                    lines = ["* xyz -2 1"]
                    lines += list(RAW[g])
                    lines.append("  Ne>  %.6f  %18.8f %18.8f %18.8f  NewECP" % (qq, *p))
                    lines.append(ECP_TAIL)
                    lines.append("*")
                    coord = "\n".join(lines) + "\n"
            (d / "job.inp").write_text(head + "\n" + "\n".join(body) + "\n" + coord)
            for f in sorted(os.listdir(d)):
                man.append("%s  %s/%s" % (hashlib.sha256((d / f).read_bytes()).hexdigest(),
                                          d.name, f))
    (out / "MANIFEST.sha256").write_text("\n".join(man) + "\n")
    n = len(runset()) * 2
    print("   generated %d input directories in %s" % (n, out))
    print("   axis: from substrate atom %d (%s) at R, |o->site1690| = %.4f A" % (k, el, d0))
    print("   d / A   actual min distance to the substrate")
    seen = set()
    for name, g, d, mind in rows:
        if (d, g) in seen:
            continue
        seen.add((d, g))
        print("     %6.2f  %-2s  %.3f" % (d, g, mind))


RE_E = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
RE_NORM = "****ORCA TERMINATED NORMALLY****"
RE_LW = re.compile(r"^\s*(\d+)\s+([A-Za-z]{1,2})\s*:\s*(-?\d+\.\d+)")


def read_out(p):
    s = Path(p).read_text(errors="replace")
    e = RE_E.findall(s)
    k = s.rfind("ORBITAL ENERGIES")
    homo, occ = None, []
    if k >= 0:
        for l in s[k:].splitlines()[4:]:
            f = l.split()
            if len(f) != 4:
                if occ:
                    break
                continue
            try:
                o, _, ev = float(f[1]), float(f[2]), float(f[3])
            except ValueError:
                break
            if o > 0.5:
                occ.append(ev)
        homo = occ[-1] if occ else None
    cy = re.findall(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES", s)
    lw = []
    j = s.rfind("LOEWDIN ATOMIC CHARGES")
    if j >= 0:
        for l in s[j:].splitlines()[2:]:
            m = RE_LW.match(l)
            if not m:
                if lw:
                    break
                continue
            lw.append(float(m.group(3)))
    return dict(E=float(e[-1]) if e else None, homo=homo, lw=lw[:24],
                cyc=int(cy[-1]) if cy else None, norm=RE_NORM in s)


def d1_control(res):
    r = one(res, "DIAGNOSTIC1_RESULT.txt")
    for l in Path(r).read_text(errors="replace").splitlines():
        m = re.match(r"\s+S6-C\s+0\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", l)
        if m:
            return {"R": float(m.group(1)), "TS": float(m.group(2))}
    stop("could not read D1's committed control energies from %s" % r)


def measure(run, res, out):
    run, res, out = Path(run), Path(res), Path(out)
    if str(out.resolve()) == str(run.resolve()) or \
       str(out.resolve()).startswith(str(run.resolve()) + os.sep):
        stop("refusing to write inside the run directory being read: %s" % out)
    RS = runset()
    o, u, _, _, _ = axis(res)
    DVAL = {lab: math.dist(o, position(res, lab)) for lab in LADDER}
    fr = frame(res)
    G = {"R": read_xyz(fr / "reactant_qm.xyz"),
         "TS": read_xyz(fr / "transition_state_qm.xyz")}

    # GATE 0: the runs must be the committed inputs, byte for byte.
    cand = [m for m in Path(res).rglob("MANIFEST.sha256")
            if ".git" not in m.parts and (m.parent / "ctl_ma_R").is_dir()]
    if len(cand) != 1:
        stop("expected exactly one MANIFEST.sha256 beside the D11 runs, found %d" % len(cand))
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
    for name in RS:
        for g in GEOM:
            p = run / ("%s_%s" % (name, g)) / "job.out"
            if not p.is_file():
                stop("missing output: %s" % p)
            d = read_out(p)
            if not d["norm"] or d["E"] is None:
                stop("%s did not terminate normally or has no final energy" % p)
            D[(name, g)] = d

    def ddE(name):
        b = RS[name][3]
        return ((D[(name, "TS")]["E"] - D[(name, "R")]["E"])
                - (D[("ctl_%s" % b, "TS")]["E"] - D[("ctl_%s" % b, "R")]["E"])) * HARTREE

    def classical(name):
        dc, _, q, b = RS[name]
        p = position(res, dc)
        tot = {}
        for g in GEOM:
            Q = D[("ctl_%s" % b, g)]["lw"]
            if len(Q) != 24:
                return None
            tot[g] = sum(q * Qa / (math.dist(p, a[1:4]) / BOHR) for Qa, a in zip(Q, G[g]))
        return (tot["TS"] - tot["R"]) * HARTREE

    O = ["STAGE 6 - DIAGNOSTIC 11 RESULT: how far must a designed charge stand off?",
         "",
         "Measured by cs11_standoff.py. NO GEOMETRY MOVED; single points on Diagnostic 1's",
         "committed R and TS. No continuum anywhere (STAGE6_AMENDMENT2_20261027). The rule is",
         "STAGE6_DIAGNOSTIC11_20261027.txt D11.3, fixed before the runs and applied unchanged.",
         "STAGE6_CRITERIA.txt is untouched.",
         "",
         "GATE 0: THE RUNS MUST BE THE COMMITTED INPUTS",
         "", "   %d files checked against MANIFEST.sha256, all match." % nman,
         "", "   GATE 0 PASSED", ""]

    O += ["GATE 1: THE ma-def2-SVP CONTROL MUST REPRODUCE DIAGNOSTIC 1'S OWN ENERGIES", ""]
    c1, g1 = d1_control(res), True
    for g in GEOM:
        got, want = D[("ctl_ma", g)]["E"], c1[g]
        ok = abs(got - want) <= GATE1_EH
        g1 &= ok
        O.append("   %-3s  measured %.9f   D1 %.9f   diff %+.2e Eh   %s"
                 % (g, got, want, got - want, "OK" if ok else "FAIL"))
    O += ["", "   GATE 1 %s" % ("PASSED" if g1 else "FAILED"), ""]
    if not g1:
        return finish(O, out, "   The control is not D1's calculation; nothing below is interpretable.")

    O += ["GATE 2: THE COMMITTED POSITION MUST REPRODUCE DIAGNOSTIC 10", "",
          "   dref is site 1690 exactly - S6-P - which D10 measured at %.6f kcal/mol."
          % D10_BARE_296, ""]
    got = ddE("dref_bare_p1_ma")
    g2 = abs(got - D10_BARE_296) <= GATE2_KCAL
    O += ["   measured %.6f   D10 %.6f   diff %+.4f kcal/mol   %s"
          % (got, D10_BARE_296, got - D10_BARE_296, "OK" if g2 else "FAIL"),
          "", "   GATE 2 %s" % ("PASSED" if g2 else "FAILED"), ""]
    if not g2:
        return finish(O, out, "   The harness has drifted from D10; nothing below is interpretable.")

    O += ["C-A1.1  A WALL CARRYING NO CHARGE MUST DO NOTHING", "",
          "   Diagnostic 9 found the wall inert in general. This re-establishes it on this axis.",
          "   Tolerance %.2f kcal/mol; a distance that fails has its walled result voided." % WALL0_KCAL,
          "", "   d / A    ddE(wall, q=0)   verdict"]
    inert = {}
    for dc in LADDER:
        v = ddE("%s_wall0_ma" % dc)
        inert[dc] = abs(v) <= WALL0_KCAL
        O.append("   %6s   %+13.4f    %s" % (dc, v, "inert" if inert[dc] else "NOT INERT"))
    O.append("")

    O += ["THE LADDER: BARE POINT CHARGE AGAINST THE SAME CHARGE BEHIND A PAULI WALL", "",
          "   q = +1, ma-def2-SVP, in vacuo. 'classical' is the same charge acting on the",
          "   control's frozen Loewdin density - it omits real polarisation, so the ratio is not",
          "   expected to reach 1; what matters is whether it plateaus.", "",
          "   d / A   ddE bare   ddE wall   d_shield   bar      classical  ratio   wall"]
    L = {}
    for dc in LADDER:
        b = ddE("%s_bare_p1_ma" % dc)
        w = ddE("%s_wall_p1_ma" % dc)
        cl = classical("%s_bare_p1_ma" % dc)
        bar = max(SHIELD_REL * abs(b), SHIELD_ABS)
        rat = (b / cl) if cl else float("nan")
        L[dc] = dict(b=b, w=w, s=b - w, bar=bar, cl=cl, rat=rat)
        O.append("   %6s  %9.4f  %9.4f  %+9.4f  %7.4f  %9.4f  %6.3f  %s"
                 % (dc, b, w, b - w, bar, cl if cl else float("nan"), rat,
                    "ok" if inert[dc] else "VOID"))
    O.append("")

    ds = list(LADDER)
    r_safe, r_at = None, None
    for i in range(len(ds) - 1):
        d, nxt = ds[i], ds[i + 1]
        if not (inert[d] and inert[nxt]):
            continue
        c_i = abs(L[d]["s"]) <= L[d]["bar"] and abs(L[nxt]["s"]) <= L[nxt]["bar"]
        c_ii = abs(L[nxt]["s"]) <= abs(L[d]["s"]) + 1e-9
        if c_i and c_ii:
            r_safe = DVAL[d]
            r_at = (d, nxt)
            break

    O += ["D11.3  THE PRE-REGISTERED READING", "",
          "   r_safe is the smallest d at which BOTH hold, AT THAT d AND THE NEXT d OUT:",
          "     (i)   |d_shield| <= max(10%% of |ddE(bare)|, %.2f kcal/mol)" % SHIELD_ABS,
          "     (ii)  |d_shield| is non-increasing outward",
          "",
          "   The classical-ratio plateau is a QUALIFIER, reported below, not a third gate. It",
          "   measures whether the physics has settled, not whether the charge model is faithful,",
          "   and AND-ing it would let a poor classical reference at short range void a genuine",
          "   r_safe. Shielding is what r_safe is named for and what it is decided by.",
          ""]
    if r_safe is None:
        O += ["   r_safe NOT REACHED within the ladder.",
              "",
              "   The bare point-charge model does not converge onto the shielded one anywhere",
              "   between 2.50 and 6.00 A. On D11.3 this means bare point charges are not a usable",
              "   model for designed charges at these separations, and the charge model - not the",
              "   grid spacing - is what has to change."]
    else:
        O += ["   r_safe = %.2f A" % r_safe, ""]
        if r_safe <= 2.96:
            O += ["   The committed arrangements sit at or beyond r_safe. Overpolarisation is NOT",
                  "   the explanation for Stage 6's behaviour, and the surrogate's failure is",
                  "   physics rather than artefact. The grid stands."]
        else:
            O += ["   r_safe is OUTSIDE the committed arrangements, whose closest sphere sits at",
                  "   2.85-2.96 A. On D11.3 the grid must be rebuilt with this stand-off and every",
                  "   a coefficient recomputed. No design measured inside r_safe is interpretable."]
    if r_safe is not None:
        d, nxt = r_at
        pl = (L[d]["cl"] and L[nxt]["cl"]
              and abs(L[nxt]["rat"] - L[d]["rat"]) <= PLATEAU_REL * abs(L[d]["rat"]))
        O += ["", "   QUALIFIER: the classical ratio %s between %s (%.3f) and %s (%.3f)."
              % ("HAS plateaued" if pl else "has NOT plateaued", d, L[d]["rat"], nxt, L[nxt]["rat"])]
        if not pl:
            O += ["   r_safe stands on the shielding test, but the response is still changing with",
                  "   distance there, so the magnitude at r_safe should not be treated as converged."]
    O.append("")

    O += ["C-A1.4  SIGN CONTROL: THE SAME LADDER WITH q = -1", "",
          "   Dreuw & Cederbaum (cr0104227) warn that a negative charge near an anion can drive the",
          "   system into a resonance that bound-state methods cannot describe. If the pathology is",
          "   that, it appears here and not on the q = +1 arm.", "",
          "   d / A   ddE bare   ddE wall   d_shield   HOMO(R) bare / eV"]
    for dc in ("dref", "d450"):
        b = ddE("%s_bare_m1_ma" % dc)
        w = ddE("%s_wall_m1_ma" % dc)
        h = D[("%s_bare_m1_ma" % dc, "R")]["homo"]
        O.append("   %6s  %9.4f  %9.4f  %+9.4f   %s"
                 % (dc, b, w, b - w, ("%+.4f" % h) if h is not None else "-"))
    O += ["", "C-A1.5  BASIS CONTROL: THE SAME POINTS AT def2-SVP", "",
          "   Senn & Thiel predict the bare-minus-walled gap SHRINKS with a compact basis. If it",
          "   does, def2-SVP conceals the artefact rather than being free of it.", "",
          "   d / A   d_shield ma   d_shield svp   ratio svp/ma"]
    for dc in ("dref", "d450"):
        sma = L[dc]["s"]
        ssv = ddE("%s_bare_p1_svp" % dc) - ddE("%s_wall_p1_svp" % dc)
        O.append("   %6s  %12.4f  %13.4f  %12s"
                 % (dc, sma, ssv, ("%.3f" % (ssv / sma)) if abs(sma) > 1e-9 else "-"))

    far = ddE("d1500_bare_p1_ma")
    fcl = classical("d1500_bare_p1_ma")
    O += ["", "A3.4  FAR-FIELD ZERO CONTROL", "",
          "   A charge at 15 A must do almost nothing, and must agree with the classical estimate.",
          "", "   ddE measured %+.4f   classical %+.4f   diff %+.4f kcal/mol   %s"
          % (far, fcl, far - fcl, "OK" if abs(far - fcl) <= 0.1 else "FAIL"), ""]

    O += ["WHAT THIS DOES NOT TEST", "",
          "   The environment. These are in vacuo and the control's density is unbound (D10 gave",
          "   HOMO +1.2269 eV); that is the separate question STAGE6_PROBE1 begins on.",
          "   It does not test additivity, which needs more than one charge. It does not revisit",
          "   any Diagnostic 1 or Diagnostic 10 number. And the Ne ECP parameters are used outside",
          "   the octet assumption they were published under, as STAGE6_ERRATUM7 records.", ""]
    return finish(O, out, None)


def finish(O, out, tail):
    if tail:
        O += [tail, ""]
    out.mkdir(parents=True, exist_ok=True)
    (out / "DIAGNOSTIC11_RESULT.txt").write_text("\n".join(O) + "\n")
    print("\n".join(O))


def stage(run, res, out):
    run, res, out = Path(run), Path(res), Path(out)
    dest = res / "phase2.2" / "cs_stage6" / "standoff"
    if not dest.is_dir():
        stop("not found: %s" % dest)
    r = out / "DIAGNOSTIC11_RESULT.txt"
    if not r.is_file():
        stop("no DIAGNOSTIC11_RESULT.txt in %s - run measure first" % out)
    lines, n = [], 0
    for name in runset():
        for g in GEOM:
            s, d = run / ("%s_%s" % (name, g)), dest / ("%s_%s" % (name, g))
            if not d.is_dir():
                stop("%s is not a committed input directory" % d)
            if not (s / "job.out").is_file():
                stop("missing output: %s" % (s / "job.out"))
            d.joinpath("job.out").write_bytes((s / "job.out").read_bytes())
            lines.append("%s  %s/job.out"
                         % (hashlib.sha256((d / "job.out").read_bytes()).hexdigest(), d.name))
            n += 1
    (dest / "DIAGNOSTIC11_RESULT.txt").write_bytes(r.read_bytes())
    lines.append("%s  DIAGNOSTIC11_RESULT.txt"
                 % hashlib.sha256((dest / "DIAGNOSTIC11_RESULT.txt").read_bytes()).hexdigest())
    (dest / "OUTPUTS.sha256").write_text("\n".join(lines) + "\n")
    print("   staged %d job.out files and the result into %s" % (n, dest))


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "generate":
        generate(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 5 and sys.argv[1] == "measure":
        measure(sys.argv[2], sys.argv[3], sys.argv[4])
    elif len(sys.argv) == 5 and sys.argv[1] == "stage":
        stage(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        stop(__doc__.strip().splitlines()[-3].strip())
