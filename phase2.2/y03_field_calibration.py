#!/usr/bin/env python3
"""
y03_field_calibration.py -- READ-ONLY re-derivation of the Stage 6 charge
distributions, of the electrostatic field the committed designs apply, and of
the same field for the Stage 4c arrangements, so that both sit on one axis.

Nothing is recalculated with ORCA. Every number below is derived from files
already committed in the results repo:

    phase2.2/cs_stage6/relax/ARRANGEMENTS.txt        (designs: site, q, xyz, a, L)
    phase2.2/cs_stage6/frame_42436/reactant_qm.xyz   (24-atom substrate at R)
    phase2.2/cs_stage6/frame_42436/transition_state_qm.xyz   (at TS, if present)
    phase2.2/cs_stage4c/L1_<arr>_R/combined.xyz      (substrate + sphere positions)
    phase2.2/cs_stage4c/L0_<arr>_R/job.out.xz        (the Ne> charges for <arr>)

The Stage 4c arrangements are read as a pair: the geometry from the L1 run and
the charges from its L0 counterpart, which places the same arrangement on Ne>
ECP centres. The sphere coordinates are checked to agree between the two before
anything is computed, and the run refuses if they do not.

The script writes ONE file into the output directory given on the command line
and touches nothing inside the repo.

    python3 y03_field_calibration.py <results-repo> <outdir> [arr ...]

<arr> defaults to E5 E7. Any Stage 4c arrangement name may be given.

Pure standard library: no numpy, so it is safe on the hpc1 login node.
"""

import hashlib
import lzma
import math
import os
import re
import sys

# e^2 / (4 pi eps0), in eV*Angstrom. Potential of +1 e at 1 A is 14.3996 V, so
# the field of q (in e) at r (in A) is 14.3996*q/r^2 in V/A.
KE = 14.39964548

# --------------------------------------------------------------------------
# reference values taken from outside this repo. Printed for comparison only;
# nothing in this script is derived from them.
# --------------------------------------------------------------------------
EXTERNAL = [
    ("WT BsCM field along the breaking C-O, no solvent",      "-0.20 V/A",
     "Siddiqui & Dubey, PCCP 2022, 24, 1974. doi:10.1039/d1cp03978d"),
    ("WT BsCM field along the breaking C-O, with solvent",    "-0.53 V/A",
     "same"),
    ("R90A variant",                                          "-0.133 V/A", "same"),
    ("R90G variant",                                          "-0.095 V/A", "same"),
    ("KSI, field on the bound steroid C=O (vibrational Stark)", "1.44 V/A",
     "Fried, Bagchi & Boxer, Science 2014, 346, 1510. doi:10.1126/science.1259802"),
    ("B3LYP critical field, heterolysis of Me3C-Cl",          "1.54 V/A",
     "Shaik et al., Chem Soc Rev 2018, 47, 5125, ESI. doi:10.1039/C8CS00354H"),
    ("B3LYP critical field, heterolysis of I2",               "2.57 V/A", "same"),
    ("glutarate dianion, experimental electron binding energy", "0.550 +/- 0.090 eV",
     "Ding, Wang & Wang, JPCA 1998, 102, 8633, via NIST WebBook"),
]
# empirical fit for linear dicarboxylate dianions, EB and repulsive Coulomb
# barrier against the charge-centroid separation R (Angstrom).
# SECONDARY SOURCE -- read from a figure caption via a summary, NOT from the
# paper. Verify against Verlet, Horke & Chatterley, PCCP 2014, 16, 15043
# (doi:10.1039/c4cp01667j) before quoting.
DICARB_EB = lambda R: 3.21 - 16.7 / R
DICARB_RCB = lambda R: 16.8 / R


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def read_xyz(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh]
    n = int(lines[0].split()[0])
    atoms = []
    for l in lines[2:2 + n]:
        p = l.split()
        atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
    if len(atoms) != n:
        raise SystemExit("xyz: declared %d atoms, found %d in %s" % (n, len(atoms), path))
    return atoms


def parse_arrangements(path):
    """Return (designs, ranking, recorded_sha) from ARRANGEMENTS.txt."""
    designs, ranking, recorded = {}, [], {}
    cur = None
    re_hdr = re.compile(r"^\s{2}(S6-\S+)\s+(\d+)\s+sphere")
    re_site = re.compile(
        r"^\s+site\s+(\d+)\s+q\s+([+-]?[\d.]+)\s+at\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
        r"\s+a\s+([+-]?[\d.]+)\s+L\s+([+-]?[\d.]+)")
    re_rank = re.compile(
        r"^\s+site\s+(\d+)\s+([-\d.]+)\s+at q\s+([+-][\d.]+)\s+admits\s+\[([-\d.]+),\s*([-\d.]+)\]")
    re_sq2 = re.compile(r"sum q\^2\s+([\d.]+)\s+\(true\)")
    re_gap = re.compile(r"strategy\s+(\S+)\s+certified gap\s+([\d.]+)")
    re_sha = re.compile(r"^\s+(\S+\.xyz)\s+sha256\s+([0-9a-f]{64})")
    for line in open(path):
        m = re_sha.match(line)
        if m:
            recorded[m.group(1)] = m.group(2)
            continue
        m = re_rank.match(line)
        if m:
            ranking.append((int(m.group(1)), float(m.group(2)),
                            float(m.group(4)), float(m.group(5))))
            continue
        m = re_hdr.match(line)
        if m:
            cur = m.group(1).strip()
            designs[cur] = {"n_declared": int(m.group(2)), "sites": [],
                            "sq2_recorded": None, "strategy": None, "gap": None}
            continue
        if cur:
            m = re_sq2.search(line)
            if m:
                designs[cur]["sq2_recorded"] = float(m.group(1))
            m = re_gap.search(line)
            if m:
                designs[cur]["strategy"] = m.group(1)
                designs[cur]["gap"] = float(m.group(2))
            m = re_site.match(line)
            if m:
                designs[cur]["sites"].append(dict(
                    site=int(m.group(1)), q=float(m.group(2)),
                    xyz=(float(m.group(3)), float(m.group(4)), float(m.group(5))),
                    a=float(m.group(6)), L=float(m.group(7))))
    return designs, ranking, recorded


def efield(charges, pt):
    """Field vector at pt from a list of (q, x, y, z). V/A."""
    E = [0.0, 0.0, 0.0]
    for q, x, y, z in charges:
        d = (pt[0] - x, pt[1] - y, pt[2] - z)
        r = math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
        if r < 1e-9:
            raise SystemExit("a charge coincides with an atom; refusing")
        f = KE * q / r ** 3
        E[0] += f * d[0]; E[1] += f * d[1]; E[2] += f * d[2]
    return E


def norm(v):
    return math.sqrt(sum(c * c for c in v))


def unit(a, b):
    d = [b[i] - a[i] for i in range(3)]
    n = norm(d)
    return tuple(c / n for c in d), n


def midpoint(a, b):
    return tuple((a[i] + b[i]) / 2.0 for i in range(3))


def identify_bonds(R, TS):
    """Breaking C-O and forming C-C, from the R->TS change where TS is
    available; otherwise from connectivity at R alone."""
    pos = lambda A, i: (A[i][1], A[i][2], A[i][3])
    n = len(R)
    if TS is not None:
        if len(TS) != n:
            raise SystemExit("R and TS have different atom counts")
        best_break, best_form = None, None
        for i in range(n):
            for j in range(i + 1, n):
                e = {R[i][0], R[j][0]}
                dR = math.dist(pos(R, i), pos(R, j))
                dT = math.dist(pos(TS, i), pos(TS, j))
                if e == {"C", "O"} and dR < 1.8:
                    if best_break is None or dT - dR > best_break[2]:
                        best_break = (i, j, dT - dR, dR, dT)
                if e == {"C"} and 2.0 < dR < 4.0:
                    if best_form is None or dT - dR < best_form[2]:
                        best_form = (i, j, dT - dR, dR, dT)
        return best_break, best_form, "R->TS change"
    # fallback: the ether oxygen is the O bonded to two carbons; the breaking
    # C-O is the longer of its two C-O bonds.
    ether = None
    for i in range(n):
        if R[i][0] != "O":
            continue
        cs = [j for j in range(n) if R[j][0] == "C"
              and math.dist(pos(R, i), pos(R, j)) < 1.8]
        if len(cs) == 2:
            ether = (i, cs)
            break
    if ether is None:
        raise SystemExit("no ether oxygen found and no TS geometry supplied")
    i, cs = ether
    j = max(cs, key=lambda c: math.dist(pos(R, i), pos(R, c)))
    brk = (i, j, float("nan"), math.dist(pos(R, i), pos(R, j)), float("nan"))
    # forming C-C: the closest non-bonded C...C pair in 3.0-3.5 A
    cand = []
    for a in range(n):
        for b in range(a + 1, n):
            if R[a][0] == R[b][0] == "C":
                d = math.dist(pos(R, a), pos(R, b))
                if 3.0 < d < 3.5:
                    cand.append((a, b, float("nan"), d, float("nan")))
    frm = min(cand, key=lambda t: abs(t[3] - 3.225)) if cand else None
    return brk, frm, "connectivity at R only (no TS geometry)"


def carboxylate_centroids(A):
    """Return the two carboxylate O-pair centroids and their separation."""
    pos = lambda i: (A[i][1], A[i][2], A[i][3])
    n = len(A)
    groups = []
    for c in range(n):
        if A[c][0] != "C":
            continue
        os_ = [o for o in range(n) if A[o][0] == "O" and math.dist(pos(c), pos(o)) < 1.6]
        if len(os_) == 2:
            cen = tuple(sum(pos(o)[k] for o in os_) / 2.0 for k in range(3))
            groups.append((c, os_, cen))
    return groups


def read_text_maybe_xz(path):
    if path.endswith(".xz"):
        with lzma.open(path, "rt", errors="replace") as fh:
            return fh.read()
    with open(path, errors="replace") as fh:
        return fh.read()


RE_NE = re.compile(
    r"Ne>\s+([+-]?\d+\.\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)")


def read_stage4c(repo, arr, state="R"):
    """Geometry from the L1 run, charges from its L0 counterpart.

    Returns (substrate_atoms, charges, note) or raises ValueError with a reason.
    The sphere coordinates must agree between the two levels.
    """
    base = os.path.join(repo, "phase2.2", "cs_stage4c")
    l1 = os.path.join(base, "L1_%s_%s" % (arr, state), "combined.xyz")
    if not os.path.isfile(l1):
        raise ValueError("no L1 geometry at %s" % l1)
    atoms = read_xyz(l1)
    # the substrate is the leading 24 atoms; the spheres are whatever follows
    n_sub = 24
    if len(atoms) <= n_sub:
        raise ValueError("%s holds %d atoms, expected more than %d"
                         % (l1, len(atoms), n_sub))
    sub, sph = atoms[:n_sub], atoms[n_sub:]
    l0 = None
    for cand in ("job.out.xz", "job.out"):
        p = os.path.join(base, "L0_%s_%s" % (arr, state), cand)
        if os.path.isfile(p):
            l0 = p
            break
    if l0 is None:
        raise ValueError("no L0 counterpart for %s_%s, so the charges are unknown"
                         % (arr, state))
    hits = RE_NE.findall(read_text_maybe_xz(l0))
    if len(hits) != len(sph):
        raise ValueError("L1 has %d spheres but L0 prints %d Ne> centres"
                         % (len(sph), len(hits)))
    charges = []
    for (q, x, y, z), s in zip(hits, sph):
        q, x, y, z = float(q), float(x), float(y), float(z)
        if math.dist((x, y, z), (s[1], s[2], s[3])) > 1e-4:
            raise ValueError("sphere coordinates disagree between L1 and L0 for %s" % arr)
        charges.append((q, x, y, z))
    return sub, charges, "geometry %s, charges %s" % (
        os.path.relpath(l1, repo), os.path.relpath(l0, repo))


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    repo, outdir = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    ARRS = sys.argv[3:] or ["E5", "E7"]
    if outdir.startswith(repo + os.sep) or outdir == repo:
        raise SystemExit("REFUSING: the output directory is inside the results repo. "
                         "This script is read-only with respect to the repo.")
    os.makedirs(outdir, exist_ok=True)
    s6 = os.path.join(repo, "phase2.2", "cs_stage6")
    arr = os.path.join(s6, "relax", "ARRANGEMENTS.txt")
    rxyz = os.path.join(s6, "frame_42436", "reactant_qm.xyz")
    txyz = os.path.join(s6, "frame_42436", "transition_state_qm.xyz")
    for p in (arr, rxyz):
        if not os.path.isfile(p):
            raise SystemExit("missing required committed file: %s" % p)

    out = []
    w = out.append
    w("STAGE 6 - y03: FIELD CALIBRATION OF THE COMMITTED ARRANGEMENTS")
    w("")
    w("Re-derived by y03_field_calibration.py from committed files only. NO NEW")
    w("CALCULATION. Nothing here is a barrier and nothing here is quoted as a result.")
    w("")

    # ---- 1. inputs and their checksums -----------------------------------
    w("1. INPUTS AND CHECKSUMS")
    w("")
    designs, ranking, recorded = parse_arrangements(arr)
    w("   ARRANGEMENTS.txt      sha256 %s" % sha256(arr))
    ok = True
    for name, path in (("reactant_qm.xyz", rxyz), ("transition_state_qm.xyz", txyz)):
        if not os.path.isfile(path):
            w("   %-24s ABSENT" % name)
            continue
        got = sha256(path)
        exp = recorded.get(name)
        if exp is None:
            w("   %-24s sha256 %s   (no value recorded in ARRANGEMENTS)" % (name, got))
        elif got == exp:
            w("   %-24s sha256 %s   MATCHES ARRANGEMENTS" % (name, got))
        else:
            ok = False
            w("   %-24s sha256 %s   *** DOES NOT MATCH %s ***" % (name, got, exp))
    if not ok:
        w("")
        w("   A GEOMETRY DOES NOT MATCH ITS RECORDED CHECKSUM. That is the finding;")
        w("   everything below is void until it is resolved.")
        open(os.path.join(outdir, "Y03_FIELD_CALIBRATION_RESULT.txt"), "w").write("\n".join(out) + "\n")
        raise SystemExit(1)

    R = read_xyz(rxyz)
    TS = read_xyz(txyz) if os.path.isfile(txyz) else None
    w("")
    w("   substrate: %d atoms at R%s" % (len(R), ", %d at TS" % len(TS) if TS else ""))

    # ---- 2. the charge distributions -------------------------------------
    w("")
    w("2. ARE THE CHARGES EVEN, SATURATED, OR PROPORTIONAL TO -a ?")
    w("")
    w("   decision 3 of OPTIMISER_DECISIONS_20261020.md predicts that under an L2")
    w("   budget with a linear objective and no bound active, q_i = -a_i*sqrt(B/sum a^2)")
    w("   and sum a q = -sqrt(B * sum a^2). Tested here against each committed design.")
    w("")
    w("   %-7s %3s %8s %8s %9s %9s %7s %8s" %
      ("design", "n", "sum|q|", "max|q|", "sum q^2", "recorded", "#|q|", "gap"))
    order = sorted(designs, key=lambda d: designs[d]["n_declared"])
    for d in order:
        D = designs[d]
        qs = [s["q"] for s in D["sites"]]
        if not qs:
            w("   %-7s %3d %8s %8s %9s %9s %7s %8s" %
              (d, D["n_declared"], "-", "-", "-", "-", "-",
               "-" if D["gap"] is None else "%.5f" % D["gap"]))
            continue
        s1 = sum(abs(q) for q in qs); s2 = sum(q * q for q in qs)
        rec = D["sq2_recorded"]
        flag = "" if rec is None or abs(s2 - rec) < 5e-5 else "  <-- MISMATCH"
        w("   %-7s %3d %8.4f %8.4f %9.5f %9s %7d %8s%s" %
          (d, D["n_declared"], s1, max(abs(q) for q in qs), s2,
           "-" if rec is None else "%.5f" % rec,
           len(set(round(abs(q), 4) for q in qs)),
           "-" if D["gap"] is None else "%.5f" % D["gap"], flag))
    w("")
    w("   '#|q|' is the number of DISTINCT charge magnitudes in the design. A value")
    w("   of 1 would mean the charge is split evenly; a value equal to n with every")
    w("   |q| at a bound would mean saturation.")
    w("")
    w("   proportionality test, q_i / (-a_i) per design:")
    w("")
    w("   %-7s %10s %10s %10s %9s %12s %12s" %
      ("design", "mean", "min", "max", "spread%", "pred sum a q", "actual"))
    for d in order:
        D = designs[d]
        if not D["sites"]:
            continue
        rat = [s["q"] / -s["a"] for s in D["sites"] if abs(s["a"]) > 1e-9]
        if not rat:
            continue
        mu = sum(rat) / len(rat)
        sd = math.sqrt(sum((r - mu) ** 2 for r in rat) / len(rat))
        B = sum(s["q"] ** 2 for s in D["sites"])
        sa2 = sum(s["a"] ** 2 for s in D["sites"])
        pred = -math.sqrt(B * sa2)
        act = sum(s["a"] * s["q"] for s in D["sites"])
        w("   %-7s %10.4f %10.4f %10.4f %9.1f %12.4f %12.4f" %
          (d, mu, min(rat), max(rat), sd / abs(mu) * 100 if mu else 0, pred, act))
    w("")
    w("   A near-constant ratio means the charges ARE the analytic interior optimum")
    w("   of the L2 budget, i.e. proportional to each site's own -a. Outliers are")
    w("   sites sitting at an admitted bound.")

    # ---- 3. what the top-ranked sites admit ------------------------------
    if ranking:
        w("")
        w("3. WHAT THE TOP-RANKED SINGLE SITES ADMIT")
        w("")
        degenerate = [r for r in ranking if abs(r[2] - r[3]) < 1e-9]
        w("   %d sites in the ranking block; %d of them admit a SINGLE charge value"
          % (len(ranking), len(degenerate)))
        w("   (lower bound equal to upper bound), so the charge there is not a")
        w("   decision variable at all.")
        w("")
        seen = {}
        for s, sc, lo, hi in ranking:
            seen["[%.2f, %.2f]" % (lo, hi)] = seen.get("[%.2f, %.2f]" % (lo, hi), 0) + 1
        for k in sorted(seen, key=lambda k: -seen[k]):
            w("      admits %-14s %d site(s)" % (k, seen[k]))
        used_top = {}
        topset = set(r[0] for r in ranking)
        for d in order:
            ov = [s["site"] for s in designs[d]["sites"] if s["site"] in topset]
            used_top[d] = ov
        w("")
        w("   which designs use a top-ranked site:")
        for d in order:
            if designs[d]["sites"]:
                w("      %-7s %s" % (d, used_top[d] if used_top[d] else "none"))
        w("")
        w("   CONSEQUENCE: a per-site charge cap below the degenerate value makes")
        w("   those sites infeasible, i.e. it deletes them from the grid rather than")
        w("   restraining them. And the Stage 5 floors are monotone in q, so changing")
        w("   the cap also moves the layer. A per-site cap is a STAGE 5 input.")

    # ---- 4. the reacting bonds -------------------------------------------
    w("")
    w("4. THE REACTING BONDS")
    w("")
    brk, frm, how = identify_bonds(R, TS)
    w("   identified by: %s" % how)
    lab = lambda A, i: "%s%d" % (A[i][0], i + 1)
    pos = lambda A, i: (A[i][1], A[i][2], A[i][3])
    if brk:
        i, j, dd, dR, dT = brk
        w("   breaking C-O : %s-%s  %.3f A at R%s" %
          (lab(R, i), lab(R, j), dR, "" if TS is None else "  %.3f A at TS" % dT))
    if frm:
        i, j, dd, dR, dT = frm
        w("   forming  C-C : %s-%s  %.3f A at R%s" %
          (lab(R, i), lab(R, j), dR, "" if TS is None else "  %.3f A at TS" % dT))
    w("")
    w("   Cross-check these against STAGE6_FINDINGS_20261022.md Finding 7, which")
    w("   records the forming C1-C6 at 3.225 A and the breaking O3-C4 at 1.495 A at R.")

    # ---- 5. the field --------------------------------------------------
    w("")
    w("5. THE FIELD THE DESIGNS APPLY")
    w("")
    w("   E = %.5f * q / r^2  V/A, in vacuo, no screening. Computed at every" % KE)
    w("   substrate atom, and projected onto each reacting bond at its midpoint.")
    w("   1 V/A = 100 MV/cm.")
    w("")
    hdr = "   %-7s %8s %9s %9s" % ("design", "sum|q|", "mean|E|", "max|E|")
    if brk: hdr += " %11s" % "E//(C-O)"
    if frm: hdr += " %11s" % "E//(C-C)"
    w(hdr)
    FIELDS = []   # (label, sum|q|, max|q|, sum q^2, mean|E|, max|E|, n)
    for d in order:
        D = designs[d]
        if not D["sites"]:
            continue
        ch = [(s["q"],) + s["xyz"] for s in D["sites"]]
        mags = [norm(efield(ch, pos(R, k))) for k in range(len(R))]
        qs = [s["q"] for s in D["sites"]]
        FIELDS.append((d, sum(abs(q) for q in qs), max(abs(q) for q in qs),
                       sum(q * q for q in qs), sum(mags) / len(mags), max(mags),
                       len(qs)))
        row = "   %-7s %8.4f %9.3f %9.3f" % (
            d, sum(abs(s["q"]) for s in D["sites"]),
            sum(mags) / len(mags), max(mags))
        for bond in (brk, frm):
            if not bond:
                continue
            i, j = bond[0], bond[1]
            u, _ = unit(pos(R, i), pos(R, j))
            E = efield(ch, midpoint(pos(R, i), pos(R, j)))
            row += " %11.3f" % sum(E[k] * u[k] for k in range(3))
        w(row)
    w("")
    w("   single charge at the top-ranked site, for the qrefine series:")
    if ranking:
        top = ranking[0][0]
        site = None
        for d in order:
            for s in designs[d]["sites"]:
                if s["site"] == top:
                    site = s
                    break
            if site:
                break
        if site:
            dmin = min(math.dist(site["xyz"], pos(R, k)) for k in range(len(R)))
            w("      site %d, closest substrate atom at %.4f A" % (top, dmin))
            for q in (0.25, 0.50, 0.65, 0.80, 1.00):
                w("         q = %.2f e  ->  %6.3f V/A  = %4.0f MV/cm"
                  % (q, KE * q / dmin ** 2, KE * q / dmin ** 2 * 100))

    # ---- 6. carboxylate separation and the dianion -----------------------
    w("")
    w("6. CHARGE SEPARATION AND THE DIANION'S ELECTRONIC STABILITY")
    w("")
    groups = carboxylate_centroids(R)
    if len(groups) == 2:
        (c1, o1, cen1), (c2, o2, cen2) = groups
        sep = math.dist(cen1, cen2)
        w("   two carboxylate groups found: C%d(O%d,O%d) and C%d(O%d,O%d)"
          % (c1 + 1, o1[0] + 1, o1[1] + 1, c2 + 1, o2[0] + 1, o2[1] + 1))
        w("   O-pair centroid separation R = %.3f A" % sep)
        w("")
        w("   empirical relation for linear dicarboxylate dianions:")
        w("      EB  = 3.21 - 16.7/R  ->  %+.3f eV" % DICARB_EB(sep))
        w("      RCB =        16.8/R  ->  %.2f eV" % DICARB_RCB(sep))
        w("")
        w("   *** THE FIT COEFFICIENTS ARE A SECONDARY SOURCE. They were read from a")
        w("   figure caption via a summary, not from the paper. Verify against Verlet,")
        w("   Horke & Chatterley, PCCP 2014, 16, 15043, doi:10.1039/c4cp01667j before")
        w("   this estimate is quoted anywhere. The relation is also fitted to LINEAR")
        w("   alkyl dicarboxylates, and chorismate is neither linear nor alkyl. ***")
    else:
        w("   expected 2 carboxylate groups, found %d; skipping." % len(groups))

    # ---- 7. external comparators ----------------------------------------
    w("")
    w("7. EXTERNAL COMPARATORS (printed for comparison; nothing above derives from them)")
    w("")
    for what, val, src in EXTERNAL:
        w("   %-52s %-20s" % (what, val))
        w("      %s" % src)

    # ---- 8. the Stage 4c arrangements on the same axis -------------------
    w("")
    w("8. THE STAGE 4c ARRANGEMENTS ON THE SAME AXIS")
    w("")
    w("   Erratum 4 E4.5 records that cs6_audit.py opens only tarballs, so 415 completed")
    w("   runs in loose .xz files were never classified, and that four Stage 4c runs at")
    w("   three point charges are STRAINED at def2-SVP - putting the def2-SVP onset at")
    w("   THREE charges, not the ten the published table shows. Those arrangements are")
    w("   read here and placed beside the Stage 6 designs.")
    w("")
    for arr in ARRS:
        try:
            sub4, ch4, note = read_stage4c(repo, arr)
        except ValueError as e:
            w("   %-7s could not be read: %s" % (arr, e))
            continue
        mags = [norm(efield(ch4, (a[1], a[2], a[3]))) for a in sub4]
        qs = [c[0] for c in ch4]
        FIELDS.append((arr, sum(abs(q) for q in qs), max(abs(q) for q in qs),
                       sum(q * q for q in qs), sum(mags) / len(mags), max(mags),
                       len(qs)))
        w("   %-7s %s" % (arr, note))
        w("           %d spheres, charges %s" %
          (len(ch4), " ".join("%+.2f" % q for q in qs)))
        near = min(math.dist((c[1], c[2], c[3]), (a[1], a[2], a[3]))
                   for c in ch4 for a in sub4)
        w("           closest sphere-to-substrate %.3f A" % near)
        w("")

    if len(FIELDS) > 1:
        w("   %-7s %4s %8s %8s %9s %9s %9s" %
          ("label", "n", "sum|q|", "max|q|", "sum q^2", "mean|E|", "max|E|"))
        for f in FIELDS:
            w("   %-7s %4d %8.3f %8.2f %9.3f %9.3f %9.3f" %
              (f[0], f[6], f[1], f[2], f[3], f[4], f[5]))
        w("")
        w("   V/A. NOTE the Stage 6 designs are frame 42436 and the Stage 4c")
        w("   arrangements are a different frame, so frame is confounded with")
        w("   arrangement in any comparison across the two blocks.")
        w("")
        w("   WHICH SCALAR, IF ANY, ORDERS THE VERDICTS is left to the reader and to a")
        w("   corpus-wide test once cs6_audit.py is fixed. This file supplies the")
        w("   quantities and states the verdicts it was given; it fits nothing and")
        w("   asserts no threshold. Verdicts, from DIAGNOSTIC6_RESULT.txt and from")
        w("   Erratum 4 E4.5, at def2-SVP:")
        w("      CLEAN     S6-C S6-P S6-2 S6-3 S6-5")
        w("      STRAINED  S6-10 S6-20, and the four Stage 4c npc=3 runs E5/E7 at R and TS")

    w("")
    w("WHAT THIS DOES NOT TEST")
    w("")
    w("   Nothing was computed with ORCA, so this cannot test any hypothesis the")
    w("   committed files do not already contain. The field is computed in vacuo with")
    w("   no screening, from point charges at the committed sphere centres, at the")
    w("   REACTANT geometry; ddE is a TS-minus-R differential, so the two are not")
    w("   like for like. Field at a bond midpoint is not the whole of the differential")
    w("   stabilisation of a pericyclic reaction, in which charge redistributes over")
    w("   the entire substrate. Two frames, one basis, a handful of arrangements, and")
    w("   the Stage 4c block is the only one that breaks the Stage 6 designs' shared")
    w("   budget, so it carries the whole weight of any cross-block comparison.")
    w("   NOTHING HERE IS A PREDICTOR. Two electrostatic predictors have already been")
    w("   proposed and tested and failed in this project (ERRATUM3 E3.5, and sum q.V_R")
    w("   against the D6 table). A third must be tested against the whole corpus, after")
    w("   cs6_audit.py's walker and charge detection are fixed together, before it is")
    w("   believed.")

    dest = os.path.join(outdir, "Y03_FIELD_CALIBRATION_RESULT.txt")
    with open(dest, "w") as fh:
        fh.write("\n".join(out) + "\n")
    print("\n".join(out))
    print()
    print("wrote %s" % dest)


if __name__ == "__main__":
    main()
