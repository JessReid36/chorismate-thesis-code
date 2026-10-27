#!/usr/bin/env python3
"""
STAGE 6 - DIAGNOSTIC 12.  What predicts the overpolarisation artefact, and can the grid be
screened on it?

Pre-registered in STAGE6_DIAGNOSTIC12_20261027.txt, which is not edited by this script.

Diagnostic 11 measured r_safe = 3.50 A along one axis through one site.  This measures the same
artefact at 36 REAL grid sites and asks which geometric property of a site predicts it.  The
per-atom floor table first proposed was abandoned before it was written, for the three reasons
recorded in D12.1; the sample here replaces it.

No geometry moves.  No continuum anywhere (STAGE6_AMENDMENT2_20261027).

  cs12_sites.py generate <results-repo> <out-dir>
  cs12_sites.py measure  <run-dir> <results-repo> <out-dir>
  cs12_sites.py stage    <run-dir> <results-repo> <out-dir>
"""
import sys, os, re, math, hashlib
from pathlib import Path

HARTREE = 627.5094740631
BOHR = 0.529177210903
GEOM = ("R", "TS")

GATE1_EH = 1e-6           # the control must reproduce D1's own energies
GATE2_KCAL = 0.01         # site 1690 must reproduce D10's S6-P and D11's dref
WALL0_KCAL = 0.30         # C-A12.1: a wall carrying no charge must do nothing
SHIELD_REL = 0.10         # PASS: |d_shield| <= 10% of |ddE(bare)| ...
SHIELD_ABS = 0.05         #       ... or this absolute floor, whichever is larger
FF_KCAL = 0.05            # far-field: |d_shield| must vanish
RHO_BAR = 0.80            # adequacy (i)

D10_S6P = -3.857852       # D10's S6-P vacuum ddE  = D11's dref bare
D11_DREF_WALL = -7.8221   # D11's dref wall
D11_FF1_BARE = -0.2860    # D11's d1500 bare

# sha256 that phase2.2/cs_stage6/relax/ARRANGEMENTS.txt records for the sites file.  Carried as a
# tripwire only; the check below reads the value out of ARRANGEMENTS.txt and compares to the file
# itself, so the authority is the committed record, not this constant.
SITES_SHA = "d2a750077d79e34dfecf448c30c2899ae0cea98f1edc5de3d10361b7c357cfc7"

# THE FROZEN SAMPLE.  Selected once by the procedure in D12.2 and fixed here.  The procedure is
# NOT re-run: a pre-registration that re-derives its own sample is not a pre-registration.
SITES = ("2009", "2051", "1988", "1975", "1959", "1991", "1981", "1967", "1463", "2271",
         "1766", "1469", "1464", "1494", "1690", "2363", "2308", "1732", "1710", "1715",
         "1698", "2090", "2373", "1679", "55", "96", "77", "2328", "0", "2330",
         "158", "326", "396", "1043", "397", "4")
SIGN_SITES = ("55", "96", "77", "0", "158", "326", "396", "1043", "397", "4")
BASIS_SITES = ("2009", "1464", "55", "4")
GATE_SITE = "1690"

LBL = ("C1", "H1", "H2", "C2", "C3", "O1", "O2", "O3", "C4", "H3", "C5", "H4",
       "C6", "C7", "H5", "C8", "H6", "C9", "H7", "O4", "H8", "C10", "O5", "O6")

# The Ne coreless-ECP block, copied byte for byte from the committed, already-run input
# phase2.2/cs_stage4c/L0_E5_R/job.inp.  STAGE6_ERRATUM7 records that these are Marefat Khah's
# published Ne parameters used outside their octet assumption; that caveat is carried, not
# resolved, here.
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
    for s in SITES:
        if out[s]["qpos"] == "none":
            stop("site %s does not admit a positive charge; the sample is wrong" % s)
    for s in SIGN_SITES:
        if out[s]["qneg"] != "-1.00":
            stop("site %s does not admit q = -1; the sign arm is wrong" % s)
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


def far_field(res, S, GR):
    """The two 15 A positions, both from the committed R geometry.

    ff1 is Diagnostic 11's d1500 exactly: 15.00 A from the substrate atom nearest site 1690
    along the vector to it.  Its bare arm therefore reproduces D11's committed -0.2860.
    ff2 is 15.00 A from the substrate centroid along the centroid -> site 0 vector, a second
    direction so the control is not a single line.
    """
    A = S[GATE_SITE]["p"]
    k, d0 = min(((i, math.dist(A, a[1:4])) for i, a in enumerate(GR)), key=lambda t: t[1])
    if abs(d0 - 2.9596) > 1e-3:
        stop("nearest substrate atom to site %s is %.4f A, D11 recorded 2.9596" % (GATE_SITE, d0))
    o = GR[k][1:4]
    u = [(A[j] - o[j]) / d0 for j in range(3)]
    ff1 = tuple(o[j] + u[j] * 15.0 for j in range(3))
    c = [sum(a[1 + j] for a in GR) / len(GR) for j in range(3)]
    B = S["0"]["p"]
    n = math.dist(c, B)
    v = [(B[j] - c[j]) / n for j in range(3)]
    ff2 = tuple(c[j] + v[j] * 15.0 for j in range(3))
    return {"ff1": ff1, "ff2": ff2}


def runset():
    """name -> (where, arm, q, basis).  'where' is a site id or 'ff1'/'ff2'/None."""
    R = {"ctl_ma": (None, "none", 0.0, "ma")}
    for s in SITES:
        R["s%s_bare_p1_ma" % s] = (s, "bare", +1.0, "ma")
        R["s%s_wall_p1_ma" % s] = (s, "wall", +1.0, "ma")
        R["s%s_wall0_ma" % s] = (s, "wall0", 0.0, "ma")
    for s in SIGN_SITES:
        R["s%s_bare_m1_ma" % s] = (s, "bare", -1.0, "ma")
        R["s%s_wall_m1_ma" % s] = (s, "wall", -1.0, "ma")
    for s in BASIS_SITES:
        R["s%s_bare_p1_svp" % s] = (s, "bare", +1.0, "svp")
        R["s%s_wall_p1_svp" % s] = (s, "wall", +1.0, "svp")
    for f in ("ff1", "ff2"):
        R["%s_bare_p1_ma" % f] = (f, "bare", +1.0, "ma")
        R["%s_wall_p1_ma" % f] = (f, "wall", +1.0, "ma")
    return R


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
    FF = far_field(res, S, G["R"])
    POS = {s: S[s]["p"] for s in SITES}
    POS.update(FF)

    out.mkdir(parents=True, exist_ok=True)
    man = []
    for name, (where, arm, q, b) in sorted(runset().items()):
        for g in GEOM:
            d = out / ("%s_%s" % (name, g))
            d.mkdir(parents=True, exist_ok=True)
            (d / "substrate.xyz").write_text(
                "24\n%s %s, from %s\n%s\n" % (name, g, fr, "\n".join(RAW[g])))
            head = "! B3LYP D3BJ %s def2/J RIJCOSX TightSCF" % BASIS[b]
            body = ["%maxcore 3000", "%pal nprocs 8 end", "%scf MaxIter 250 end"]
            if arm == "none":
                coord = "* xyzfile -2 1 substrate.xyz\n"
            else:
                p = POS[where]
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
    print("   %d sites, %d sign, %d basis, 2 far field, 1 control" %
          (len(SITES), len(SIGN_SITES), len(BASIS_SITES)))
    dmin = min(min(math.dist(POS[s], a[1:4]) for a in G["R"]) for s in SITES)
    dmax = max(min(math.dist(POS[s], a[1:4]) for a in G["R"]) for s in SITES)
    print("   sample d_nearest at R: %.3f to %.3f A" % (dmin, dmax))
    for f in ("ff1", "ff2"):
        print("   %s: %.3f A from the nearest substrate atom at R"
              % (f, min(math.dist(FF[f], a[1:4]) for a in G["R"])))


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


def spearman(x, y):
    """Rank correlation.  Ties take their average rank."""
    def rank(v):
        o = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[o[k]] = avg
            i = j + 1
        return r
    a, b = rank(x), rank(y)
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((p - ma) * (q - mb) for p, q in zip(a, b))
    den = math.sqrt(sum((p - ma) ** 2 for p in a) * sum((q - mb) ** 2 for q in b))
    return num / den if den else 0.0


def separates(x, ok):
    """Adequacy (ii): is there a threshold on x that splits PASS from FAIL with no mistakes?

    Both orientations are tried, because some descriptors grow with safety (distance) and others
    shrink with it (proximity sums).  Every midpoint between consecutive distinct values is a
    candidate; the best one is returned whether or not it is perfect, so the number of sites a
    descriptor misorders is reported rather than just whether it failed.

    Returns (threshold, orientation, misclassified).  misclassified == 0 means adequate on (ii).
    """
    pts = sorted(set(x))
    if len(pts) < 2:
        return (None, None, len(x))
    best = None
    for i in range(len(pts) - 1):
        t = (pts[i] + pts[i + 1]) / 2.0
        for orient in (">=", "<="):
            if orient == ">=":
                miss = sum(1 for xi, oki in zip(x, ok) if (xi >= t) != oki)
            else:
                miss = sum(1 for xi, oki in zip(x, ok) if (xi <= t) != oki)
            if best is None or miss < best[2]:
                best = (t, orient, miss)
    return best


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
    FF = far_field(res, S, G["R"])
    POS = {s: S[s]["p"] for s in SITES}
    POS.update(FF)

    # GATE 0: the runs must be the committed inputs, byte for byte.
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
        b = RS[name][3]
        c = "ctl_ma" if b == "ma" else None
        if c is None:
            return None
        return ((D[(name, "TS")]["E"] - D[(name, "R")]["E"])
                - (D[(c, "TS")]["E"] - D[(c, "R")]["E"])) * HARTREE

    def classical(name):
        where, _, q, b = RS[name]
        if b != "ma" or where is None:
            return None
        p = POS[where]
        tot = {}
        for g in GEOM:
            Q = D[("ctl_ma", g)]["lw"]
            if len(Q) != 24:
                return None
            tot[g] = sum(q * Qa / (math.dist(p, a[1:4]) / BOHR) for Qa, a in zip(Q, G[g]))
        return (tot["TS"] - tot["R"]) * HARTREE

    O = ["STAGE 6 - DIAGNOSTIC 12 RESULT: what predicts the overpolarisation artefact?",
         "",
         "Measured by cs12_sites.py against STAGE6_DIAGNOSTIC12_20261027.txt, which was committed",
         "before the runs. NO GEOMETRY MOVED; single points on Diagnostic 1's committed R and TS.",
         "No continuum anywhere (STAGE6_AMENDMENT2_20261027). 36 real grid sites, one +1 charge",
         "each, bare point charge against the same charge behind a Pauli wall.",
         ""]

    # ---- GATE 1 -------------------------------------------------------------------------
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

    # ---- GATE 2 -------------------------------------------------------------------------
    gb = ddE("s%s_bare_p1_ma" % GATE_SITE)
    gw = ddE("s%s_wall_p1_ma" % GATE_SITE)
    O += ["GATE 2: SITE %s MUST REPRODUCE WHAT IS ALREADY COMMITTED FOR THAT POSITION" % GATE_SITE,
          "",
          "   Site %s is Diagnostic 10's S6-P and Diagnostic 11's dref rung." % GATE_SITE, ""]
    g2 = True
    for lab, got, want in (("bare", gb, D10_S6P), ("wall", gw, D11_DREF_WALL)):
        ok = abs(got - want) <= GATE2_KCAL
        g2 &= ok
        O.append("   %-4s  measured %+.4f   committed %+.4f   diff %+.4f kcal/mol   %s"
                 % (lab, got, want, got - want, "OK" if ok else "FAIL"))
    O += ["", "   GATE 2 %s" % ("PASSED" if g2 else "FAILED"), ""]
    if not g2:
        return finish(O, out, "GATE 2 failed: the harness has drifted. Nothing below is read.")

    # ---- C-A12.1 wall inertness, per site -----------------------------------------------
    void = set()
    O += ["C-A12.1  THE WALL MUST BE INERT AT EACH SITE", "",
          "   An ECP centre carrying no charge must give |ddE| <= %.2f kcal/mol. Diagnostic 11"
          % WALL0_KCAL,
          "   found the wall NOT inert at its closest rung, so this is checked at every site.",
          "   A site that fails has its WALLED result voided; it enters no fit.", ""]
    for s in SITES:
        v = ddE("s%s_wall0_ma" % s)
        if abs(v) > WALL0_KCAL:
            void.add(s)
            O.append("   site %-5s wall0 ddE %+.4f   VOID" % (s, v))
    O += ["   %d of %d sites voided%s" % (len(void), len(SITES),
                                          "" if void else " - the wall is inert everywhere"), ""]

    # ---- the reading --------------------------------------------------------------------
    rows = []
    for s in SITES:
        b, w = ddE("s%s_bare_p1_ma" % s), ddE("s%s_wall_p1_ma" % s)
        ds = b - w
        tol = max(SHIELD_REL * abs(b), SHIELD_ABS)
        d = descriptors(POS[s], G["R"])
        rows.append(dict(site=s, bare=b, wall=w, ds=ds, tol=tol, ok=abs(ds) <= tol,
                         voided=s in void, cls=classical("s%s_bare_p1_ma" % s), **d))
    rows.sort(key=lambda r: r["dn"])

    O += ["THE 36 SITES", "",
          "   PASS means |d_shield| <= max(0.10 x |ddE bare|, %.2f kcal/mol) - Diagnostic 11's"
          % SHIELD_ABS,
          "   tolerance, unchanged. It is relative because an absolute one passes automatically",
          "   wherever the quantities are small.", "",
          "   site   near  d_near   d_O    S1/d2  S_O1/d2  ddE bare  ddE wall  d_shield   tol"
          "    verdict"]
    for r in rows:
        O.append("   %-5s  %-4s %6.3f %6.3f  %6.3f  %6.3f  %+8.4f  %+8.4f  %+8.4f %6.3f   %s"
                 % (r["site"], r["near"], r["dn"], r["dO"], r["inv2"], r["invO"],
                    r["bare"], r["wall"], r["ds"], r["tol"],
                    "VOID" if r["voided"] else ("pass" if r["ok"] else "FAIL")))
    live = [r for r in rows if not r["voided"]]
    npass = sum(1 for r in live if r["ok"])
    O += ["", "   %d live, %d pass, %d fail, %d voided"
          % (len(live), npass, len(live) - npass, len(rows) - len(live)), ""]

    # ---- adequacy of the five descriptors ------------------------------------------------
    O += ["WHICH DESCRIPTOR, IF ANY, PREDICTS THE ARTEFACT", "",
          "   Five descriptors, named in the pre-registration before the runs. ADEQUATE means",
          "   BOTH |Spearman rho| >= %.2f against |d_shield| AND a threshold that separates pass" % RHO_BAR,
          "   from fail with ZERO misclassifications. The separation is the operative test.", ""]
    DESCS = (("X1", "d_nearest", lambda r: r["dn"]),
             ("X2", "sum 1/d^2", lambda r: r["inv2"]),
             ("X3", "sum_O 1/d^2", lambda r: r["invO"]),
             ("X4", "d_O", lambda r: r["dO"]),
             ("X5", "|ddE classical|", lambda r: abs(r["cls"]) if r["cls"] is not None else None))
    adequate = []
    if len(live) < 3 or npass == 0 or npass == len(live):
        O += ["   NOT EVALUABLE: the live sites are %s, so no threshold can be tested."
              % ("all pass" if npass == len(live) else
                 ("all fail" if npass == 0 else "too few")), ""]
    else:
        O.append("   id  descriptor         rho     threshold      misclassified   adequate")
        y = [abs(r["ds"]) for r in live]
        for tag, nm, f in DESCS:
            x = [f(r) for r in live]
            if any(v is None for v in x):
                O.append("   %-3s %-17s  (not available)" % (tag, nm))
                continue
            rho = spearman(x, y)
            ok = [r["ok"] for r in live]
            thr, orient, miss = separates(x, ok)
            sep = (miss == 0)
            good = abs(rho) >= RHO_BAR and sep
            if good:
                adequate.append((tag, nm, thr, orient))
            O.append("   %-3s %-17s %+.3f   %-22s   %2d              %s"
                     % (tag, nm, rho,
                        ("pass if X %s %.4f" % (orient, thr)) if sep else "none separates cleanly",
                        miss, "YES" if good else "no"))
        O.append("")

    # ---- far field ----------------------------------------------------------------------
    O += ["FAR-FIELD CONTROL", "",
          "   At 15 A the bare and walled arms must agree. This is a check on d_shield only.",
          "   STAGE6_DIAGNOSTIC11_ERRATUM1 E11.1 records that D11's A3.4 was mis-specified",
          "   because it required measured ~ classical, and the classical reference omits a",
          "   polarisation enhancement present at every distance. That cannot recur here: the",
          "   classical reference does not enter this control.", ""]
    ffok = True
    for f in ("ff1", "ff2"):
        b, w = ddE("%s_bare_p1_ma" % f), ddE("%s_wall_p1_ma" % f)
        ok = abs(b - w) <= FF_KCAL
        ffok &= ok
        O.append("   %s  bare %+.4f   wall %+.4f   d_shield %+.4f   %s"
                 % (f, b, w, b - w, "OK" if ok else "FAIL"))
    b1 = ddE("ff1_bare_p1_ma")
    ok1 = abs(b1 - D11_FF1_BARE) <= GATE2_KCAL
    ffok &= ok1
    O += ["   ff1 is D11's d1500: bare %+.4f against D11's committed %+.4f, diff %+.4f   %s"
          % (b1, D11_FF1_BARE, b1 - D11_FF1_BARE, "OK" if ok1 else "FAIL"), ""]

    # ---- sign arm -----------------------------------------------------------------------
    O += ["C-A12.2  THE SIGN ARM (reported, not scored)", "",
          "   The same ten anion-layer positions at q = -1. D11's C-A1.4 found the artefact",
          "   specific to positive charges; ten sites in one layer cannot separate sign from",
          "   position, so this is reported only.", "",
          "   site   d_near   d_shield(+1)   d_shield(-1)"]
    for s in SIGN_SITES:
        dp = ddE("s%s_bare_p1_ma" % s) - ddE("s%s_wall_p1_ma" % s)
        dm = ddE("s%s_bare_m1_ma" % s) - ddE("s%s_wall_m1_ma" % s)
        dn = descriptors(POS[s], G["R"])["dn"]
        O.append("   %-5s  %6.3f     %+8.4f      %+8.4f" % (s, dn, dp, dm))
    O.append("")

    # ---- basis arm ----------------------------------------------------------------------
    O += ["BASIS ARM (reported, not scored)", "",
          "   D11 found at its dref rung that def2-SVP reports d_shield -0.067 where ma-def2-SVP",
          "   reports +3.964, a ratio of 0.017, and concluded the compact basis conceals the",
          "   artefact. This asks whether that is the basis or that one axis.", "",
          "   site   d_shield(ma)   d_shield(svp)   ratio"]
    for s in BASIS_SITES:
        dm = ddE("s%s_bare_p1_ma" % s) - ddE("s%s_wall_p1_ma" % s)
        es = ((D[("s%s_bare_p1_svp" % s, "TS")]["E"] - D[("s%s_bare_p1_svp" % s, "R")]["E"])
              - (D[("s%s_wall_p1_svp" % s, "TS")]["E"] - D[("s%s_wall_p1_svp" % s, "R")]["E"])) \
             * HARTREE
        O.append("   %-5s    %+8.4f       %+8.4f      %s"
                 % (s, dm, es, ("%.3f" % (es / dm)) if abs(dm) > 1e-9 else "n/a"))
    O.append("")

    # ---- verdict ------------------------------------------------------------------------
    O += ["VERDICT", ""]
    if not ffok:
        O += ["   The far-field control FAILED. The oracle reports shielding where none can",
              "   exist, so no reading below it is trustworthy and the screen is not issued.", ""]
    elif npass == 0:
        O += ["   EVERY LIVE SITE FAILS.",
              "",
              "   No position the grid admits is a faithful one for a bare +1 charge, so there is",
              "   nothing for a screen to select and the descriptors were not testable. D12.3 gives",
              "   this the same reading as no descriptor being adequate, reached sooner: the charge",
              "   model must change, not the grid spacing.", ""]
    elif npass == len(live):
        O += ["   EVERY LIVE SITE PASSES.",
              "",
              "   This contradicts Diagnostic 11's measurement at site 1690, which is in the sample",
              "   and which D11 measured as a FAIL. GATE 2 passed, so the harness reproduces D11",
              "   at that position - the two cannot both be right, and nothing here is read until",
              "   that is resolved.", ""]
    elif not adequate:
        O += ["   NO DESCRIPTOR IS ADEQUATE.",
              "",
              "   Geometry does not predict the artefact. No screen can be built from the",
              "   quantities a grid has. Every candidate site would need its own QM check, which",
              "   the grid-and-budget formulation cannot absorb. The charge model must change -",
              "   not the spacing, and not the surface. Moving the grid outward does not help,",
              "   because what fails is the prediction, not the distance.", ""]
    else:
        tag, nm, t, orient = next(((a, b, c, d) for a, b, c, d in adequate if a == "X1"),
                                  adequate[0])
        O += ["   ADEQUATE: %s" % ", ".join("%s (%s)" % (a, b) for a, b, _, _ in adequate),
              "   THE SCREEN USES %s (%s): admit a site iff %s %s %.4f" % (tag, nm, nm, orient, t),
              ""]
        if tag == "X1":
            O += ["   Nearest-atom distance is sufficient. A single distance floor is a valid",
                  "   screen, Diagnostic 11's r_safe generalises beyond its one axis, and the",
                  "   grid rebuild applies this threshold to all 2626 sites.", ""]
        else:
            O += ["   Nearest-atom distance is NOT sufficient and is abandoned as the admission",
                  "   criterion. Every a coefficient computed at a site the new screen rejects is",
                  "   withdrawn.", ""]
        keep = [r for r in rows if not r["voided"]]
        O += ["   Applied to the sample: %d of %d live sites admitted."
              % (sum(1 for r in keep if r["ok"]), len(keep)), ""]

    O += ["WHAT THIS DOES NOT TEST", "",
          "   Additivity. Every run here carries ONE charge; a design carries several, and",
          "   nothing here says the artefact adds. The Ne ECP parameters remain used outside",
          "   their published octet assumption (ERRATUM7), so the walled arm is a contrast, not",
          "   a truth. The environment - these are in vacuo. The product geometry. Charge",
          "   magnitudes other than +-1, so any screen issued here is a screen AT |q| = 1, the",
          "   worst case the grid admits, and applies to smaller charges only as an upper bound.",
          "   And whether any site that passes is catalytic: a faithful charge model is a",
          "   precondition for reading a design, not a design."]
    return finish(O, out, None)


def finish(O, out, tail):
    if tail:
        O += [tail, ""]
    out.mkdir(parents=True, exist_ok=True)
    (out / "DIAGNOSTIC12_RESULT.txt").write_text("\n".join(O) + "\n")
    print("\n".join(O))


def stage(run, res, out):
    run, res, out = Path(run), Path(res), Path(out)
    dest = res / "phase2.2" / "cs_stage6" / "sitescan"
    if not dest.is_dir():
        stop("not found: %s" % dest)
    r = out / "DIAGNOSTIC12_RESULT.txt"
    if not r.is_file():
        stop("no DIAGNOSTIC12_RESULT.txt in %s - run measure first" % out)
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
    (dest / "DIAGNOSTIC12_RESULT.txt").write_bytes(r.read_bytes())
    lines.append("%s  DIAGNOSTIC12_RESULT.txt"
                 % hashlib.sha256((dest / "DIAGNOSTIC12_RESULT.txt").read_bytes()).hexdigest())
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
        stop("usage: cs12_sites.py generate|measure|stage ...")
