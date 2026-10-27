#!/usr/bin/env python3
"""
STAGE 6 - DIAGNOSTIC 13.  Does the overpolarisation artefact depend on the charge magnitude?

Pre-registered in STAGE6_DIAGNOSTIC13_20261027.txt, which is not edited by this script.

Diagnostics 11 and 12 both run at |q| = 1 and say so.  Diagnostic 2's committed C3 scan at site
1690 shows the linear coefficient flat to 0.28 per cent from q = 0.05 to 0.50 and then collapsing
33 per cent at q = 1.00, at 2.96 A - inside D11's r_safe.  This asks whether that holds across
distance, on eight of D12's sites spanning its full range.

The sites are PROBES, not candidate designs: seven of the eight have q_min_pos = 1.00, so the grid
would refuse them a 0.25 charge.  That rule governs where a resting sphere may sit and nothing here
moves, which is the same reasoning under which D2 scanned q at site 1690.

No geometry moves.  No continuum anywhere (STAGE6_AMENDMENT2_20261027).

  cs13_qscale.py generate <results-repo> <out-dir>
  cs13_qscale.py measure  <run-dir> <results-repo> <out-dir>
  cs13_qscale.py stage    <run-dir> <results-repo> <out-dir>
"""
import sys, os, re, math, hashlib
from pathlib import Path

HARTREE = 627.5094740631
BOHR = 0.529177210903
GEOM = ("R", "TS")
SITES_SHA = "d2a750077d79e34dfecf448c30c2899ae0cea98f1edc5de3d10361b7c357cfc7"
LBL = ("C1", "H1", "H2", "C2", "C3", "O1", "O2", "O3", "C4", "H3", "C5", "H4",
       "C6", "C7", "H5", "C8", "H6", "C9", "H7", "O4", "H8", "C10", "O5", "O6")
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
FRAME_DIR = "phase2.2/cs_stage6/frame_42436"
FRAME_SHA = {"reactant_qm.xyz":
             "6e8f41c7da3bcf4ed9f5c650dc9fcf12c2a36e5af7a900e164892c276497aa88",
             "transition_state_qm.xyz":
             "2babeb9000ca5220a908a70eaec80ec4a64358c8f5677e9a16cb104e0dffc5bf"}

GATE1_EH = 1e-6           # the control must reproduce D1's own energies
GATE2_KCAL = 0.01         # site 1690 must reproduce D2, D10 and D11
WALL0_KCAL = 0.30         # C-A13.1: a wall carrying no charge must do nothing
SHIELD_REL = 0.10         # PASS: |d_shield| <= 10% of |ddE(bare)| ...
SHIELD_ABS = 0.05         #       ... or this absolute floor, whichever is larger
NOISE_KCAL = 0.02         # below this at two or more magnitudes -> no exponent is fitted

# Diagnostic 2's committed C3 (ma-def2-SVP, in vacuo) ddE at site 1690, and D11's dref wall.
D2_C3_BARE = {0.25: -1.8935, 0.50: -3.9781, 1.00: -3.8579}
D11_DREF_WALL = -7.8221
GATE_SITE = "1690"

# THE FROZEN SAMPLE.  Eight of Diagnostic 12's thirty-six, spanning its full nearest-atom range at
# R (2.588 to 3.732 A), site 1690 forced in.  Selected once and fixed here; not re-derived.
SITES = ("2009", "1981", "1464", "1690", "1710", "1679", "2330", "4")
QLEVELS = (0.25, 0.50, 1.00)

def stop(m):
    print("STOP: %s" % m, file=sys.stderr)
    raise SystemExit(1)

def one(root, name):
    h = [p for p in Path(root).rglob(name) if ".git" not in p.parts]
    if len(h) != 1:
        stop("expected exactly one %s under %s, found %d" % (name, root, len(h)))
    return h[0]

def frame(res):
    """Stage 6's geometries, by committed path AND committed hashes (five dirs share the name)."""
    d = Path(res) / FRAME_DIR
    for f, want in FRAME_SHA.items():
        p = d / f
        if not p.is_file():
            stop("not found: %s" % p)
        got = hashlib.sha256(p.read_bytes()).hexdigest()
        if got != want:
            stop("%s sha256 is %s, ARRANGEMENTS.txt records %s" % (p, got, want))
    return d

def sites_file(res):
    """GATE 3.  The site coordinates are this diagnostic's input, so they are checked against an
    independent committed record - ARRANGEMENTS.txt's INPUTS block - and not merely trusted."""
    p = Path(res) / "phase2.2" / "cs_stage6" / "sites_f42436.tsv"
    if not p.is_file():
        stop("not found: %s" % p)
    got = hashlib.sha256(p.read_bytes()).hexdigest()
    arr = Path(one(res, "ARRANGEMENTS.txt")).read_text(errors="replace")
    m = re.search(r"sites_f42436\.tsv\s+sha256\s+([0-9a-f]{64})", arr)
    if not m:
        stop("ARRANGEMENTS.txt does not record a sha256 for sites_f42436.tsv")
    if m.group(1) != got:
        stop("sites_f42436.tsv sha256 is %s, ARRANGEMENTS.txt records %s" % (got, m.group(1)))
    if got != SITES_SHA:
        stop("sites_f42436.tsv sha256 %s does not match the pre-registered %s" % (got, SITES_SHA))
    return p

def read_sites(res):
    """The eight probe positions, from the committed grid after GATE 3.

    q_min_pos is READ AND REPORTED but NOT enforced. Seven of the eight sites admit only
    |q| >= 1 by the grid's rule, which governs where a resting sphere may sit - a statement about
    forces. Nothing in this diagnostic moves, and Diagnostic 2 scanned q = 0.05 to 1.00 at site
    1690 under the same reasoning. These are probes; D13.2 says so and none may be read as a design.
    """
    rows = [l.split("\t") for l in sites_file(res).read_text().splitlines()
            if not l.startswith("#")]
    h = rows[0]
    ix, iy, iz = h.index("x"), h.index("y"), h.index("z")
    il, ip, ineg = h.index("layer"), h.index("q_min_pos"), h.index("q_min_neg")
    out = {}
    for r in rows[1:]:
        if not r or not r[0].strip():
            continue
        out[r[0]] = dict(p=(float(r[ix]), float(r[iy]), float(r[iz])),
                         layer=r[il], qpos=r[ip], qneg=r[ineg])
    missing = [s for s in SITES if s not in out]
    if missing:
        stop("sites not present in the committed grid: %s" % ", ".join(missing))
    return out

def read_xyz_raw(p):
    """The committed coordinate lines, verbatim.

    Diagnostic 10 wrote its substrate.xyz by copying these lines through unchanged, so carrying
    them rather than reformatting makes site 1690's COORDINATE BLOCK byte-identical to S6-P's and
    gate 2 an exact reproduction rather than a near one. The files are not wholly identical: the
    xyz comment line carries the run name and the generating path, which ORCA ignores. D10, D11
    and this diagnostic all share coordinate-block sha256 5740f38f924d4330...
    """
    L = Path(p).read_text().splitlines()
    return L[2:2 + int(L[0].split()[0])]

def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    return [(f[0], float(f[1]), float(f[2]), float(f[3]))
            for f in (l.split() for l in L[2:2 + n])]

def descriptors(p, g):
    """The five pre-registered descriptors' geometric part, at one geometry.

    Computed from the coordinates of THIS geometry.  The d_ columns in sites_f42436.tsv are
    per-atom MINIMA over R, TS and P - correct for a clearance criterion, wrong for a
    configuration - so they are not used here.
    """
    d = {LBL[i]: math.dist(p, a[1:4]) for i, a in enumerate(g)}
    return dict(dn=min(d.values()),
                near=min(d, key=d.get),
                inv2=sum(1.0 / (x * x) for x in d.values()),
                invO=sum(1.0 / (d[a] ** 2) for a in LBL if a[0] == "O"),
                dO=min(d[a] for a in LBL if a[0] == "O"))

RE_E = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
RE_NORM = "****ORCA TERMINATED NORMALLY****"
RE_LW = re.compile(r"^\s*(\d+)\s+([A-Za-z]{1,2})\s*:\s*(-?\d+\.\d+)")


def read_out(p):
    s = Path(p).read_text(errors="replace")
    e = RE_E.findall(s)
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
    return dict(E=float(e[-1]) if e else None, lw=lw[:24],
                cyc=int(cy[-1]) if cy else None, norm=RE_NORM in s)

def d1_control(res):
    r = one(res, "DIAGNOSTIC1_RESULT.txt")
    for l in Path(r).read_text(errors="replace").splitlines():
        m = re.match(r"\s+S6-C\s+0\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", l)
        if m:
            return {"R": float(m.group(1)), "TS": float(m.group(2))}
    stop("could not read D1's committed control energies from %s" % r)

def finish(O, out, tail):
    if tail:
        O += [tail, ""]
    out.mkdir(parents=True, exist_ok=True)
    (out / "DIAGNOSTIC13_RESULT.txt").write_text("\n".join(O) + "\n")
    print("\n".join(O))

def qtag(q):
    """q = 0.25 -> 'q025'. Used in run names, so the name says the magnitude."""
    return "q%03d" % round(q * 100)


def runset():
    """name -> (site or None, arm, q)."""
    R = {"ctl_ma": (None, "none", 0.0)}
    for s in SITES:
        for q in QLEVELS:
            R["s%s_%s_bare" % (s, qtag(q))] = (s, "bare", q)
            R["s%s_%s_wall" % (s, qtag(q))] = (s, "wall", q)
        R["s%s_wall0" % s] = (s, "wall0", 0.0)
    return R


def generate(res, out):
    res, out = Path(res), Path(out)
    fr = frame(res)
    S = read_sites(res)
    G = {"R": read_xyz(fr / "reactant_qm.xyz"),
         "TS": read_xyz(fr / "transition_state_qm.xyz")}
    RAW = {"R": read_xyz_raw(fr / "reactant_qm.xyz"),
           "TS": read_xyz_raw(fr / "transition_state_qm.xyz")}
    for g, a in G.items():
        if len(a) != 24:
            stop("%s geometry has %d atoms, expected 24" % (g, len(a)))
    out.mkdir(parents=True, exist_ok=True)
    man = []
    for name, (site, arm, q) in sorted(runset().items()):
        for g in GEOM:
            d = out / ("%s_%s" % (name, g))
            d.mkdir(parents=True, exist_ok=True)
            (d / "substrate.xyz").write_text(
                "24\n%s %s, from %s\n%s\n" % (name, g, fr, "\n".join(RAW[g])))
            head = "! B3LYP D3BJ %s def2/J RIJCOSX TightSCF" % BASIS["ma"]
            body = ["%maxcore 3000", "%pal nprocs 8 end", "%scf MaxIter 250 end"]
            if arm == "none":
                coord = "* xyzfile -2 1 substrate.xyz\n"
            else:
                p = S[site]["p"]
                if arm == "bare":
                    (d / "pc.pc").write_text("1\n%12.6f %18.8f %18.8f %18.8f\n" % (q, *p))
                    body.append('%pointcharges "pc.pc"')
                    coord = "* xyzfile -2 1 substrate.xyz\n"
                else:
                    qq = 0.0 if arm == "wall0" else q
                    lines = ["* xyz -2 1"] + list(RAW[g])
                    lines.append("  Ne>  %.6f  %18.8f %18.8f %18.8f  NewECP" % (qq, *p))
                    lines.append(ECP_TAIL)
                    lines.append("*")
                    coord = "\n".join(lines) + "\n"
            (d / "job.inp").write_text(head + "\n" + "\n".join(body) + "\n" + coord)
            for f in sorted(os.listdir(d)):
                man.append("%s  %s/%s"
                           % (hashlib.sha256((d / f).read_bytes()).hexdigest(), d.name, f))
    (out / "MANIFEST.sha256").write_text("\n".join(man) + "\n")
    n = len(runset()) * 2
    print("   generated %d input directories in %s" % (n, out))
    print("   %d sites x %d charge magnitudes, bare and walled, plus wall0 and the control"
          % (len(SITES), len(QLEVELS)))
    print("   site   layer   q_min_pos   d_nearest(R)   admits 0.25 by the grid's rule?")
    for s in SITES:
        dn = descriptors(S[s]["p"], G["R"])["dn"]
        adm = S[s]["qpos"] != "none" and float(S[s]["qpos"]) <= 0.25
        print("   %-5s  %-7s %-9s   %6.3f         %s" % (s, S[s]["layer"], S[s]["qpos"], dn,
                                                         "yes" if adm else "no - probe only"))


def fit_exponent(qs, vals):
    """Least-squares slope of log|d_shield| against log q. Returns (p, note)."""
    pts = [(math.log(q), math.log(abs(v))) for q, v in zip(qs, vals) if abs(v) > 0.0]
    if len(pts) < 2:
        return (None, "too few non-zero points")
    n = len(pts)
    mx = sum(a for a, _ in pts) / n
    my = sum(b for _, b in pts) / n
    den = sum((a - mx) ** 2 for a, _ in pts)
    if den <= 0:
        return (None, "degenerate")
    return (sum((a - mx) * (b - my) for a, b in pts) / den, "")


def measure(run, res, out):
    run, res, out = Path(run), Path(res), Path(out)
    if str(out.resolve()) == str(run.resolve()) or \
       str(out.resolve()).startswith(str(run.resolve()) + os.sep):
        stop("refusing to write inside the run directory being read: %s" % out)
    RS = runset()
    fr = frame(res)
    S = read_sites(res)
    G = {"R": read_xyz(fr / "reactant_qm.xyz"),
         "TS": read_xyz(fr / "transition_state_qm.xyz")}

    cand = [m for m in Path(run).rglob("MANIFEST.sha256") if ".git" not in m.parts]
    if len(cand) != 1:
        stop("expected exactly one MANIFEST.sha256 in the run directory, found %d" % len(cand))
    base = cand[0].parent
    bad = []
    for line in cand[0].read_text().splitlines():
        want, rel = line.split("  ", 1)
        f = base / rel
        if not f.is_file():
            bad.append("missing %s" % rel)
        elif hashlib.sha256(f.read_bytes()).hexdigest() != want:
            bad.append("changed %s" % rel)
    if bad:
        stop("GATE 0 failed, %d input files do not match the manifest:\n   %s"
             % (len(bad), "\n   ".join(bad[:10])))

    D = {}
    for name in RS:
        for g in GEOM:
            p = base / ("%s_%s" % (name, g)) / "job.out"
            if not p.is_file():
                stop("missing output: %s" % p)
            r = read_out(p)
            if not r["norm"] or r["E"] is None:
                stop("%s_%s did not terminate normally" % (name, g))
            D[(name, g)] = r

    def ddE(name):
        return ((D[(name, "TS")]["E"] - D[(name, "R")]["E"])
                - (D[("ctl_ma", "TS")]["E"] - D[("ctl_ma", "R")]["E"])) * HARTREE

    O = ["STAGE 6 - DIAGNOSTIC 13 RESULT: does the artefact depend on the charge magnitude?",
         "",
         "Measured by cs13_qscale.py against STAGE6_DIAGNOSTIC13_20261027.txt, committed before",
         "the runs. NO GEOMETRY MOVED; single points on Diagnostic 1's committed R and TS. No",
         "continuum (STAGE6_AMENDMENT2_20261027). Eight probe positions, three charge magnitudes,",
         "bare point charge against the same charge behind a Pauli wall.",
         "",
         "These are PROBES, not designs: seven of the eight sites would be refused a 0.25 charge",
         "by the grid's own admissibility rule, which governs resting spheres. Nothing moved here.",
         ""]

    d1 = d1_control(res)
    O += ["GATE 1: THE CONTROL MUST REPRODUCE DIAGNOSTIC 1'S COMMITTED ENERGIES", ""]
    g1 = True
    for g in GEOM:
        got, want = D[("ctl_ma", g)]["E"], d1[g]
        ok = abs(got - want) <= GATE1_EH
        g1 &= ok
        O.append("   %-2s  measured %.9f   D1 %.9f   diff %+.2e Eh   %s"
                 % (g, got, want, got - want, "OK" if ok else "FAIL"))
    O += ["", "   GATE 1 %s" % ("PASSED" if g1 else "FAILED"), ""]
    if not g1:
        return finish(O, out, "GATE 1 failed. Nothing below is read.")

    O += ["GATE 2: SITE %s MUST REPRODUCE D2, D10 AND D11" % GATE_SITE, "",
          "   Four committed numbers from three diagnostics.", ""]
    g2 = True
    for q in QLEVELS:
        got = ddE("s%s_%s_bare" % (GATE_SITE, qtag(q)))
        want = D2_C3_BARE[q]
        ok = abs(got - want) <= GATE2_KCAL
        g2 &= ok
        tag = "D2 C3" + ("  = D10's S6-P" if q == 1.00 else "")
        O.append("   bare q=%.2f  measured %+.4f   %s %+.4f   diff %+.4f   %s"
                 % (q, got, tag, want, got - want, "OK" if ok else "FAIL"))
    gotw = ddE("s%s_q100_wall" % GATE_SITE)
    okw = abs(gotw - D11_DREF_WALL) <= GATE2_KCAL
    g2 &= okw
    O.append("   wall q=1.00  measured %+.4f   D11 dref wall %+.4f   diff %+.4f   %s"
             % (gotw, D11_DREF_WALL, gotw - D11_DREF_WALL, "OK" if okw else "FAIL"))
    O += ["", "   GATE 2 %s" % ("PASSED" if g2 else "FAILED"), ""]
    if not g2:
        return finish(O, out, "GATE 2 failed: the harness has drifted. Nothing below is read.")

    void = set()
    O += ["C-A13.1  THE WALL MUST BE INERT AT EACH SITE", "",
          "   |ddE| <= %.2f kcal/mol for an ECP centre carrying no charge. A site that fails has"
          % WALL0_KCAL,
          "   its walled results voided at EVERY magnitude.", ""]
    for s in SITES:
        v = ddE("s%s_wall0" % s)
        if abs(v) > WALL0_KCAL:
            void.add(s)
            O.append("   site %-5s wall0 ddE %+.4f   VOID" % (s, v))
    O += ["   %d of %d sites voided%s" % (len(void), len(SITES),
                                          "" if void else " - the wall is inert everywhere"), ""]

    rows = []
    for s in SITES:
        dn = descriptors(S[s]["p"], G["R"])["dn"]
        rec = dict(site=s, dn=dn, voided=s in void, q={})
        for q in QLEVELS:
            b = ddE("s%s_%s_bare" % (s, qtag(q)))
            w = ddE("s%s_%s_wall" % (s, qtag(q)))
            ds = b - w
            tol = max(SHIELD_REL * abs(b), SHIELD_ABS)
            rec["q"][q] = dict(bare=b, wall=w, ds=ds, tol=tol, ok=abs(ds) <= tol)
        rows.append(rec)
    rows.sort(key=lambda r: r["dn"])

    O += ["THE READING: d_shield AT EACH MAGNITUDE", "",
          "   PASS means |d_shield| <= max(0.10 x |ddE bare|, %.2f) - Diagnostic 11's tolerance."
          % SHIELD_ABS, "",
          "   site   d_near |   q=0.25                  | q=0.50                  | q=1.00",
          "                  |   ddE bare  d_shield  v  |  ddE bare  d_shield  v  |  ddE bare"
          "  d_shield  v"]
    for r in rows:
        cells = []
        for q in QLEVELS:
            e = r["q"][q]
            v = "V" if r["voided"] else ("p" if e["ok"] else "F")
            cells.append("%+9.4f %+9.4f  %s" % (e["bare"], e["ds"], v))
        O.append("   %-5s  %6.3f | %s | %s | %s" % (r["site"], r["dn"], *cells))
    O.append("")
    for q in QLEVELS:
        live = [r for r in rows if not r["voided"]]
        n = sum(1 for r in live if r["q"][q]["ok"])
        O.append("   q = %.2f : %d of %d live sites pass" % (q, n, len(live)))
    O.append("")

    O += ["HOW THE ARTEFACT SCALES WITH CHARGE", "",
          "   p is the least-squares slope of log|d_shield| against log q over the three",
          "   magnitudes. The two-point slopes are printed beside it: a pure power law has them",
          "   equal, a switch-on does not. A site below %.2f kcal/mol at two or more magnitudes"
          % NOISE_KCAL,
          "   is BELOW NOISE and is not fitted - a log fit through numerical noise is not a",
          "   measurement.", "",
          "   site   d_near    p      0.25->0.50   0.50->1.00   note"]
    exps = []
    for r in rows:
        if r["voided"]:
            O.append("   %-5s  %6.3f    -          -            -          VOIDED" % (r["site"], r["dn"]))
            continue
        vals = [r["q"][q]["ds"] for q in QLEVELS]
        if sum(1 for v in vals if abs(v) < NOISE_KCAL) >= 2:
            O.append("   %-5s  %6.3f    -          -            -          BELOW NOISE"
                     % (r["site"], r["dn"]))
            continue
        p, note = fit_exponent(list(QLEVELS), vals)
        def two(a, b):
            if abs(vals[a]) <= 0 or abs(vals[b]) <= 0:
                return float("nan")
            return math.log(abs(vals[b]) / abs(vals[a])) / math.log(QLEVELS[b] / QLEVELS[a])
        p1, p2 = two(0, 1), two(1, 2)
        if p is not None:
            exps.append(p)
        O.append("   %-5s  %6.3f  %6.2f     %6.2f       %6.2f       %s"
                 % (r["site"], r["dn"], p if p is not None else float("nan"), p1, p2, note))
    O.append("")

    O += ["VERDICT", ""]
    live = [r for r in rows if not r["voided"]]
    if not live:
        O += ["   Every site was voided by C-A13.1. Nothing is read.", ""]
    else:
        f1 = [r for r in live if not r["q"][1.00]["ok"]]
        rescued = [r for r in f1 if r["q"][0.25]["ok"]]
        if f1 and len(rescued) == len(f1):
            O += ["   EVERY SITE THAT FAILS AT |q| = 1 PASSES AT |q| = 0.25.",
                  "",
                  "   r_safe is a distance AT A CHARGE, not a distance. Diagnostic 11's 3.50 A and",
                  "   whatever threshold Diagnostic 12 returns are |q| = 1 results and must not be",
                  "   applied to weaker charges. The admissible region is a (distance, |q|) region",
                  "   and the grid rebuild must carry both.", ""]
            passing = [r for r in live if r["q"][0.25]["ok"]]
            if passing:
                O += ["   Smallest nearest-atom distance passing at 0.25 e in this sample: %.3f A"
                      % min(r["dn"] for r in passing), ""]
        elif f1 and not rescued:
            O += ["   NO SITE IS RESCUED BY REDUCING THE CHARGE.",
                  "",
                  "   The artefact is not a large-|q| phenomenon. Diagnostic 2's observation at",
                  "   site 1690 does not generalise, D11's floor stands at every magnitude, and",
                  "   the many-weak-charges route buys nothing on this axis.", ""]
        elif f1:
            O += ["   PARTIAL: %d of %d sites failing at |q| = 1 pass at 0.25."
                  % (len(rescued), len(f1)),
                  "",
                  "   The magnitude dependence is real but site-specific. A single cap cannot be",
                  "   stated from this sample; the screen would have to be per-site, which the",
                  "   grid cannot carry.", ""]
        else:
            O += ["   No live site fails at |q| = 1, so there is nothing for a lower charge to",
                  "   rescue. That conflicts with Diagnostic 11 at site 1690, which GATE 2 just",
                  "   reproduced, so the two readings must be reconciled before either is used.", ""]
        if exps:
            O += ["   Fitted exponents span %.2f to %.2f (mean %.2f) over %d sites."
                  % (min(exps), max(exps), sum(exps) / len(exps), len(exps)), ""]

    O += ["WHAT THIS DOES NOT TEST", "",
          "   ADDITIVITY, and it is the limitation that matters most. Every run carries ONE",
          "   charge. A design of twenty-six charges at 0.25 e is not twenty-six independent",
          "   probes. Finding 2 measured ddE additivity at def2-SVP (3.1 % at n = 2, 6.1 % at",
          "   n = 3); the ARTEFACT's additivity has never been measured at any n, and this does",
          "   not measure it. Nothing here supports a claim that weak charges are safe in a",
          "   many-charge design.",
          "",
          "   Eight sites, a sample of Diagnostic 12's sample. The Ne ECP parameters remain used",
          "   outside their published octet assumption (ERRATUM7). In vacuo. The product geometry.",
          "   And whether any of these positions is catalytic: they are probes, and seven of the",
          "   eight are positions the grid would refuse to a 0.25 charge."]
    return finish(O, out, None)


def stage(run, res, out):
    run, res, out = Path(run), Path(res), Path(out)
    dest = res / "phase2.2" / "cs_stage6" / "qscale"
    if not dest.is_dir():
        stop("not found: %s" % dest)
    r = out / "DIAGNOSTIC13_RESULT.txt"
    if not r.is_file():
        stop("no DIAGNOSTIC13_RESULT.txt in %s - run measure first" % out)
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
    (dest / "DIAGNOSTIC13_RESULT.txt").write_bytes(r.read_bytes())
    lines.append("%s  DIAGNOSTIC13_RESULT.txt"
                 % hashlib.sha256((dest / "DIAGNOSTIC13_RESULT.txt").read_bytes()).hexdigest())
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
        stop("usage: cs13_qscale.py generate|measure|stage ...")
