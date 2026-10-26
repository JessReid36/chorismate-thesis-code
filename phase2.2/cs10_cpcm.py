#!/usr/bin/env python3
"""Stage 6 DIAGNOSTIC 10: does ddE come back when the dianion is bound?

Diagnostic 1 measured ddE(frozen) at ma-def2-SVP against the committed def2-SVP values and returned
NOT USABLE AS THEY STAND: the moves are 50 to 118 per cent against a 15 per cent yardstick and
S6-20's sign reverses. SOLVATION_DECISION_NOTE.md predicted that failure on 29 September, before any
of it was run: the bare dianion is unbound in vacuum, so "designed charges would then be doing two
jobs at once - binding the excess electron density and catalysing the reaction".

This diagnostic asks the one question that separates those two jobs. BIND THE DENSITY WITH A
CONTINUUM, CHANGE NOTHING ELSE, AND SEE WHETHER ddE RETURNS.

HOW CPCM BECOMES AVAILABLE AT ALL. ORCA will not combine CPCM with its QMMM module, and D1's design
runs are QMMM jobs because the spheres carry Lennard-Jones parameters. Finding 4 is the way through
and it is exact: "the Lennard-Jones term never enters the QM Hamiltonian", and at a fixed geometry
the difference between an LJ sphere and a bare point charge "is exactly L, to zero". So at FIXED
GEOMETRY a design is re-expressible as %pointcharges, which CPCM does support - the qscan's C4
condition already runs ma-def2-SVP with CPCM and %pointcharges, 18 committed runs - with an offset
of sum L that ARRANGEMENTS.txt states per site. That offset is not assumed here. It is a GATE.

NO GEOMETRY MOVES. NO OPTIMISATION. Sixteen single points on the committed D1 geometries.

  generate   build the 16 input directories from committed files only
  measure    read the 16 outputs and write DIAGNOSTIC10_RESULT.txt

usage: cs10_cpcm.py generate <results-repo> <out-dir>
       cs10_cpcm.py measure  <run-dir> <results-repo> <out-dir>
"""
import sys, os, re, math, hashlib
from pathlib import Path

HARTREE = 627.5094740631          # kcal/mol per Hartree
ARR = ("S6-C", "S6-P", "S6-2", "S6-20")
GEOM = ("R", "TS")
COND = ("vac", "cpcm")
EPS = 4.0
GATE1_EH = 1e-6                   # control must reproduce D1's own energies to this
GATE2_KCAL = 0.05                 # ddE(pointcharge) must equal ddE(QMMM) - sum L to this


def stop(m):
    print("STOP: %s" % m, file=sys.stderr)
    raise SystemExit(1)


def one(root, name):
    h = [p for p in Path(root).rglob(name) if ".git" not in p.parts]
    if len(h) != 1:
        stop("expected exactly one %s under %s, found %d" % (name, root, len(h)))
    return h[0]


def read_arrangements(p):
    """Per design: list of (site, q, x, y, z, a, L). Straight from the committed file."""
    out, cur = {}, None
    for l in Path(p).read_text(errors="replace").splitlines():
        m = re.match(r"\s{2}(S6-\S+)\s+(\d+) sphere", l)
        if m:
            cur = m.group(1); out[cur] = []
            continue
        m = re.match(r"\s+site\s+(\d+)\s+q\s+([-+]\d+\.\d+)\s+at\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
                     r"\s+a\s+([-+]\d+\.\d+)\s+L\s+([-+]\d+\.\d+)", l)
        if m and cur:
            out[cur].append((int(m.group(1)), float(m.group(2)), float(m.group(3)),
                             float(m.group(4)), float(m.group(5)), float(m.group(6)), float(m.group(7))))
    return out


def read_d1_inputs(diff, arr, geom):
    """The committed D1 input pair: substrate atoms, and the sphere charges and positions."""
    d = Path(diff) / ("DIFF-" + arr) / (geom + "_sp")
    xyz = (d / "combined.xyz").read_text().splitlines()
    n = int(xyz[0].split()[0])
    body = [l for l in xyz[2:2 + n]]
    atoms = [l for l in body[:24]]
    sph = []
    for l in body[24:]:
        f = l.split()
        sph.append((float(f[1]), float(f[2]), float(f[3])))
    q = []
    pr = d / "site.ORCAFF.prms"
    if sph and not pr.is_file():
        stop("%s has %d spheres but no site.ORCAFF.prms" % (d, len(sph)))
    if pr.is_file():                      # the control has no spheres and no force field
        prms = pr.read_text().splitlines()
        k = prms.index("$atoms")
        for l in prms[k + 2: k + 2 + n]:
            f = l.split()
            q.append(float(f[2]))
        q = q[24:]
    return atoms, sph, q, d


def generate(res, out):
    res, out = Path(res), Path(out)
    diff = one(res, "DIAGNOSTIC1_RESULT.txt").parent
    arrp = one(res, "ARRANGEMENTS.txt")
    A = read_arrangements(arrp)
    out.mkdir(parents=True, exist_ok=True)
    man, made = [], 0
    for arr in ARR:
        for geom in GEOM:
            atoms, sph, q, src = read_d1_inputs(diff, arr, geom)
            if len(atoms) != 24:
                stop("%s %s: expected 24 QM atoms, found %d" % (arr, geom, len(atoms)))
            if len(sph) != len(q):
                stop("%s %s: %d sphere positions but %d charges" % (arr, geom, len(sph), len(q)))
            # cross-check against ARRANGEMENTS.txt, which is the authority for q and for L
            if arr != "S6-C":
                ref = A.get(arr)
                if ref is None:
                    stop("%s is not in %s" % (arr, arrp))
                if len(ref) != len(sph):
                    stop("%s: ARRANGEMENTS lists %d spheres, the D1 input has %d"
                         % (arr, len(ref), len(sph)))
                for (s_, qa, x, y, z, a_, L_), (px, py, pz), qp in zip(ref, sph, q):
                    if abs(qa - qp) > 1e-6:
                        stop("%s: charge %.6f in ARRANGEMENTS, %.6f in the D1 input" % (arr, qa, qp))
                    if max(abs(x - px), abs(y - py), abs(z - pz)) > 1e-4:
                        stop("%s: sphere position disagrees between ARRANGEMENTS and the D1 input" % arr)
            for cond in COND:
                d = out / ("%s_%s_%s" % (arr, geom, cond))
                d.mkdir(parents=True, exist_ok=True)
                (d / "substrate.xyz").write_text("24\n%s %s, from %s\n%s\n"
                                                 % (arr, geom, src, "\n".join(atoms)))
                head = "! B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF"
                body = ["%maxcore 3000", "%pal nprocs 8 end", "%scf MaxIter 250 end"]
                if cond == "cpcm":
                    head += " CPCM"
                    body.insert(0, "%%cpcm epsilon %.1f end" % EPS)
                if sph:
                    (d / "pc.pc").write_text("%d\n" % len(sph) + "".join(
                        "%12.6f %18.8f %18.8f %18.8f\n" % (qq, p[0], p[1], p[2])
                        for qq, p in zip(q, sph)))
                    body.append('%pointcharges "pc.pc"')
                (d / "job.inp").write_text(head + "\n" + "\n".join(body)
                                           + "\n* xyzfile -2 1 substrate.xyz\n")
                for f in sorted(os.listdir(d)):
                    man.append("%s  %s/%s" % (hashlib.sha256((d / f).read_bytes()).hexdigest(),
                                              d.name, f))
                made += 1
    (out / "MANIFEST.sha256").write_text("\n".join(man) + "\n")
    print("   generated %d input directories in %s" % (made, out))
    for arr in ARR:
        if arr == "S6-C":
            print("     %-6s control, no spheres" % arr)
        else:
            L = sum(r[6] for r in A[arr])
            print("     %-6s %2d sphere(s)   sum L = %+.5f kcal/mol" % (arr, len(A[arr]), L))


RE_E = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
RE_NORM = "****ORCA TERMINATED NORMALLY****"


def read_out(p):
    s = Path(p).read_text(errors="replace")
    e = RE_E.findall(s)
    homo = None
    k = s.rfind("ORBITAL ENERGIES")
    if k >= 0:
        occ = []
        for l in s[k:].splitlines()[4:]:
            f = l.split()
            if len(f) != 4:
                if occ:
                    break
                continue
            try:
                o, eh, ev = float(f[1]), float(f[2]), float(f[3])
            except ValueError:
                break
            if o > 0.5:
                occ.append(ev)
        homo = occ[-1] if occ else None
    cy = re.findall(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES", s)
    return dict(E=float(e[-1]) if e else None, homo=homo,
                cyc=int(cy[-1]) if cy else None, norm=RE_NORM in s)


def d1_values(res):
    """D1's committed numbers: the control's energies and every design's two ddE values."""
    r = one(res, "DIAGNOSTIC1_RESULT.txt")
    ctl = {}
    for l in Path(r).read_text(errors="replace").splitlines():
        m = re.match(r"\s+S6-C\s+0\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", l)
        if m:
            ctl = {"R": float(m.group(1)), "TS": float(m.group(2))}
    t = r.parent / "diagnostic1.tsv"
    dd = {}
    for l in t.read_text(errors="replace").splitlines():
        f = l.split("\t")
        if len(f) >= 4 and f[0].startswith("S6-"):
            dd[f[0]] = (float(f[2]), float(f[3]))        # ddE_ma, ddE_svp
    if not ctl or not dd:
        stop("could not read D1's committed values from %s" % r)
    return ctl, dd


def measure(run, res, out):
    run, res, out = Path(run), Path(res), Path(out)
    if str(out.resolve()).startswith(str(res.resolve())):
        stop("refusing to write inside the results repo: %s" % out)
    A = read_arrangements(one(res, "ARRANGEMENTS.txt"))
    ctl1, dd1 = d1_values(res)

    R = {}
    for arr in ARR:
        for geom in GEOM:
            for cond in COND:
                p = run / ("%s_%s_%s" % (arr, geom, cond)) / "job.out"
                if not p.is_file():
                    stop("missing output: %s" % p)
                d = read_out(p)
                if not d["norm"] or d["E"] is None:
                    stop("%s did not terminate normally or has no final energy" % p)
                R[(arr, geom, cond)] = d

    def ddE(arr, cond):
        b = (R[(arr, "TS", cond)]["E"] - R[(arr, "R", cond)]["E"])
        c = (R[("S6-C", "TS", cond)]["E"] - R[("S6-C", "R", cond)]["E"])
        return (b - c) * HARTREE

    O = ["STAGE 6 - DIAGNOSTIC 10 RESULT: does ddE come back when the dianion is bound?",
         "",
         "Measured by cs10_cpcm.py from sixteen single points on the committed Diagnostic 1 geometries.",
         "NO GEOMETRY MOVED. The verdict rule is STAGE6_DIAGNOSTIC10_20261026.txt D10.3, fixed before",
         "the runs went in and applied here unchanged. STAGE6_CRITERIA.txt is untouched.",
         "",
         "GATE 1: THE CONTROL MUST REPRODUCE DIAGNOSTIC 1'S OWN ENERGIES",
         "",
         "   S6-C carries no spheres, so its %pointcharges job is the same calculation D1 ran. Any",
         "   difference means the two are not the same calculation and nothing below is interpretable.",
         ""]
    g1 = True
    for geom in GEOM:
        got, want = R[("S6-C", geom, "vac")]["E"], ctl1[geom]
        ok = abs(got - want) <= GATE1_EH
        g1 &= ok
        O.append("   %-3s  measured %.9f   D1 %.9f   diff %+.2e Eh   %s"
                 % (geom, got, want, got - want, "OK" if ok else "FAIL"))
    O += ["", "   GATE 1 %s" % ("PASSED" if g1 else "FAILED"), ""]
    if not g1:
        O += ["   Nothing further is reported. The control does not reproduce D1, so the",
              "   QMMM-to-pointcharges conversion is not the identity Finding 4 says it is, and no",
              "   comparison below would mean anything.", ""]
        (out).mkdir(parents=True, exist_ok=True)
        (out / "DIAGNOSTIC10_RESULT.txt").write_text("\n".join(O) + "\n")
        print("\n".join(O))
        return

    O += ["GATE 2: THE VACUUM ARM MUST DIFFER FROM D1 BY EXACTLY sum L", "",
          "   Finding 4: at fixed geometry a bare point charge and an LJ sphere differ by exactly L.",
          "   So ddE(pointcharge) must equal D1's ddE(QMMM) minus sum L, to %.2f kcal/mol." % GATE2_KCAL,
          "",
          "   design   ddE vac    D1 ddE(ma)   sum L      expected   diff      "]
    g2 = True
    for arr in ARR:
        if arr == "S6-C":
            continue
        L = sum(r[6] for r in A[arr])
        exp = dd1[arr][0] - L
        got = ddE(arr, "vac")
        ok = abs(got - exp) <= GATE2_KCAL
        g2 &= ok
        O.append("   %-7s %9.4f  %10.4f  %+8.5f  %9.4f  %+8.4f  %s"
                 % (arr, got, dd1[arr][0], L, exp, got - exp, "OK" if ok else "FAIL"))
    O += ["", "   GATE 2 %s" % ("PASSED" if g2 else "FAILED"), ""]
    if not g2:
        O += ["   Nothing further is reported. Finding 4 says a bare point charge and an LJ sphere",
              "   differ by exactly L at a fixed geometry. They do not here, so the re-expression of",
              "   these designs as point charges is not the identity it has to be for the continuum",
              "   arm to be comparable with anything, and no verdict is available.", ""]
        out.mkdir(parents=True, exist_ok=True)
        (out / "DIAGNOSTIC10_RESULT.txt").write_text("\n".join(O) + "\n")
        print("\n".join(O))
        return

    O += ["THE READING: ddE WITH THE DENSITY BOUND, AGAINST THE COMMITTED def2-SVP VALUES", "",
          "   design   ddE cpcm   ddE vac    D1 ddE(svp)   move vs svp   sign vs svp"]
    moves, signs = {}, {}
    for arr in ARR:
        if arr == "S6-C":
            continue
        c, v, s = ddE(arr, "cpcm"), ddE(arr, "vac"), dd1[arr][1]
        mv = abs(c - s) / abs(s) * 100 if s else float("inf")
        sg = (c < 0) == (s < 0)
        moves[arr], signs[arr] = mv, sg
        O.append("   %-7s %9.4f  %9.4f  %12.4f  %10.1f%%   %s"
                 % (arr, c, v, s, mv, "same" if sg else "REVERSED"))
    O += ["", "   HOMO, which says whether the continuum did its job", "",
          "   design   HOMO(R) vac eV   HOMO(R) cpcm eV   bound at cpcm?"]
    for arr in ARR:
        hv, hc = R[(arr, "R", "vac")]["homo"], R[(arr, "R", "cpcm")]["homo"]
        O.append("   %-7s %14s  %17s   %s"
                 % (arr, "%+.4f" % hv if hv is not None else "-",
                    "%+.4f" % hc if hc is not None else "-",
                    "yes" if (hc is not None and hc < 0) else "NO"))
    unbound = [a for a in ARR
               if R[(a, "R", "cpcm")]["homo"] is None or R[(a, "R", "cpcm")]["homo"] >= 0]
    worst = max(moves.values()) if moves else 0.0
    allsign = all(signs.values()) if signs else False
    if unbound:
        verdict = "NO VERDICT - THE TEST WAS NOT PERFORMED"
        why = ["The continuum did not bind the species: %s %s a non-negative HOMO at CPCM."
               % (", ".join(unbound), "has" if len(unbound) == 1 else "have"),
               "D10.3 names this as a fourth outcome and not as a verdict. The instrument did not do",
               "the one thing it was introduced to do, so nothing can be concluded about whether ddE",
               "returns when the density is bound. A larger dielectric or an explicit environment",
               "would be the next thing to try, and that is a decision, not a re-run."]
    elif not allsign:
        verdict = "NOT RECOVERED"
        why = ["A sign still reverses against the committed def2-SVP value. Binding the density does",
               "not restore the design numbers, so the def2-SVP frozen barriers are not recoverable by",
               "a continuum and the frozen half must be rebuilt, as D1.3 already required."]
    elif worst <= 15.0:
        verdict = "RECOVERED"
        why = ["Every sign holds and the largest move is within the 15 per cent yardstick decision O6",
               "already allows. On this evidence the vacuum unboundedness was the artefact Diagnostic 1",
               "detected, and the frozen numbers are usable WITH the continuum, not without it."]
    else:
        verdict = "PARTIALLY RECOVERED"
        why = ["Every sign holds but at least one move exceeds 15 per cent. The ordering survives and",
               "the magnitudes are continuum-dependent; they would have to be reported with that said."]
    O += ["", "VERDICT PER D10.3: %s" % verdict, ""] + ["   " + w for w in why]
    O += ["", "WHAT THIS DOES NOT TEST", "",
          "   The geometries. These are single points on D1's committed R and TS, which are stationary",
          "   points of a different Hamiltonian; nothing here relaxes. It does not test whether a",
          "   continuum is the right environment - SOLVATION_DECISION_NOTE.md section 3 argues a",
          "   continuum cannot represent the solvent screening that matters, and that argument is",
          "   untouched by this result. It does not revisit any wavefunction-health verdict. And a",
          "   RECOVERED verdict would license the frozen numbers only in the condition tested here,",
          "   not in vacuo, which is the condition Stage 6 actually used.", ""]
    out.mkdir(parents=True, exist_ok=True)
    (out / "DIAGNOSTIC10_RESULT.txt").write_text("\n".join(O) + "\n")
    print("\n".join(O))


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "generate" and len(sys.argv) == 4:
        generate(sys.argv[2], sys.argv[3])
    elif len(sys.argv) >= 2 and sys.argv[1] == "measure" and len(sys.argv) == 5:
        measure(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        stop("usage: cs10_cpcm.py generate <results-repo> <out-dir>\n"
             "       cs10_cpcm.py measure  <run-dir> <results-repo> <out-dir>")
