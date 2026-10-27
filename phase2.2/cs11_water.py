#!/usr/bin/env python3
"""Stage 6 PROBE 1: what does an explicit water shell cost, and does it bind the dianion?

This is a COST AND BINDING probe, not a diagnostic. It runs the bare substrate only - no designed
charges, no arrangements, no ddE - so there is nothing here that could be read as a result about a
design. It answers two planning questions:

  1. How long does a ma-def2-SVP single point take at 24, 42, 60 and 84 QM atoms?
  2. Does an explicit first shell bind the chorismate dianion, whose HOMO is +1.2269 eV in vacuo?

Question 2 is the one that matters. Diagnostic 10 removed CPCM as the way of binding the density
(STAGE6_AMENDMENT2_20261027). If a small explicit shell binds it, Phase A's environment question has
a cheap answer; if it does not, the evaluation has to move to Phase 1's full QM/MM, which is a
decision for the supervisor and not for this code.

WHAT THIS IS NOT. The real first shell around the substrate in the enzyme is mostly PROTEIN, not
water: within 3.5 A of the substrate at frame 42436 there are 63 non-water atoms and 17 water atoms.
Substrate-plus-waters in vacuo is therefore a device for binding the electron density, NOT a model
of the active site. Treat any energy it produces accordingly.

  generate   build the input directories from committed files only
  measure    report wall time, SCF cycles and HOMO. It reports no ddE by design.

usage: cs11_water.py generate <results-repo> <out-dir>
       cs11_water.py measure  <run-dir> <out-dir>
"""
import sys, os, re, math, hashlib, tarfile
from pathlib import Path

SHELLS = (0, 6, 12, 20)          # waters; 0 is the in-vacuo control
ARCHIVE = "05_qmmm/24_sp_barriers/f42436_r8_a1/hpc_text_outputs_2_r_irc.tar.xz"
MEMBER = "05_qmmm/24_sp_barriers/f42436_r8_a1/r_irc/final.xyz"
OH_MAX = 1.2                     # an O with exactly two H inside this is a water
CLASH = 1.5                      # refuse if an R-frame water sits this close to a TS atom


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


def read_xyz_text(t):
    L = t.splitlines()
    n = int(L[0].split()[0])
    return [(f[0], float(f[1]), float(f[2]), float(f[3]))
            for f in (l.split() for l in L[2:2 + n])]


def waters(res):
    """The ranked water shell, derived from committed files alone.

    The solvated frame is read from the committed archive. The substrate is located inside it by
    exact coordinate match against frame_42436/reactant_qm.xyz - BUILD_RECORD.txt records that the
    two agree to 0.0e+00 A, and this re-checks it rather than trusting the record. A water is an
    oxygen outside the substrate carrying exactly two hydrogens within %.1f A. Waters are ranked by
    the shortest distance from their oxygen to any substrate atom, and the first N are taken.

    The SAME water positions are used at R and at TS. A designed charge does not move with the
    substrate, and neither does this shell; that keeps the environment a fixed, preorganised one
    rather than a reorganising one.
    """ % OH_MAX
    res = Path(res)
    arc = res / ARCHIVE
    if not arc.is_file():
        stop("committed archive not found: %s" % arc)
    with tarfile.open(arc, "r:xz") as t:
        m = t.getmember(MEMBER)
        full = read_xyz_text(t.extractfile(m).read().decode("utf-8", "replace"))
    fr = frame(res)
    sub = read_xyz_text((fr / "reactant_qm.xyz").read_text())
    if len(sub) != 24:
        stop("reactant_qm.xyz has %d atoms, expected 24" % len(sub))
    idx = []
    for s in sub:
        hit = [i for i, a in enumerate(full)
               if a[0] == s[0] and math.dist(a[1:4], s[1:4]) < 1e-4]
        if len(hit) != 1:
            stop("substrate atom %s not uniquely located in %s (%d matches)" % (s[0], MEMBER, len(hit)))
        idx.append(hit[0])
    if idx != list(range(idx[0], idx[0] + 24)):
        stop("the substrate is not contiguous inside %s" % MEMBER)
    S = set(idx)
    subxyz = [a[1:4] for a in sub]
    W = []
    for i, a in enumerate(full):
        if a[0] != "O" or i in S:
            continue
        hs = [j for j in (i + 1, i + 2) if j < len(full) and full[j][0] == "H"
              and math.dist(a[1:4], full[j][1:4]) < OH_MAX]
        if len(hs) == 2:
            W.append((min(math.dist(a[1:4], s) for s in subxyz),
                      [full[i], full[hs[0]], full[hs[1]]]))
    W.sort(key=lambda t: t[0])
    return sub, W


def generate(res, out):
    res, out = Path(res), Path(out)
    sub, W = waters(res)
    fr = frame(res)
    TS = read_xyz_text((fr / "transition_state_qm.xyz").read_text())
    out.mkdir(parents=True, exist_ok=True)
    man, made = [], 0
    jobs = [("R", n) for n in SHELLS] + [("TS", 6)]
    for geom, n in jobs:
        base = sub if geom == "R" else TS
        shell = [a for _, trio in W[:n] for a in trio]
        if shell:
            c = min(math.dist(w[1:4], a[1:4]) for w in shell for a in base)
            if c < CLASH:
                stop("%s with %d waters: closest water-substrate contact %.2f A, below %.2f"
                     % (geom, n, c, CLASH))
        d = out / ("w%02d_%s" % (n, geom))
        d.mkdir(parents=True, exist_ok=True)
        atoms = list(base) + shell
        (d / "system.xyz").write_text(
            "%d\n%s, substrate + %d waters from %s\n%s\n"
            % (len(atoms), geom, n, MEMBER,
               "\n".join("%-3s %18.8f %18.8f %18.8f" % a for a in atoms)))
        (d / "job.inp").write_text(
            "! B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF\n"
            "%maxcore 3000\n%pal nprocs 8 end\n%scf MaxIter 250 end\n"
            "* xyzfile -2 1 system.xyz\n")
        for f in sorted(os.listdir(d)):
            man.append("%s  %s/%s" % (hashlib.sha256((d / f).read_bytes()).hexdigest(), d.name, f))
        made += 1
    (out / "MANIFEST.sha256").write_text("\n".join(man) + "\n")
    rep = ["STAGE 6 - PROBE 1: the water shell as built, from committed files only", "",
           "Archive  %s" % ARCHIVE, "Member   %s" % MEMBER,
           "Rule     an O outside the substrate with exactly two H within %.1f A; ranked by the" % OH_MAX,
           "         shortest O-to-substrate distance; first N taken; same positions at R and TS.",
           "", "   rank   O-to-substrate / A"]
    for k, (dd, _) in enumerate(W[:24], 1):
        mark = "   <- shell boundary" if k in SHELLS else ""
        rep.append("   %4d   %8.3f%s" % (k, dd, mark))
    rep += ["", "The gap between rank 6 (%.3f A) and rank 7 (%.3f A) is %.3f A, which is why 6 is"
            % (W[5][0], W[6][0], W[6][0] - W[5][0]),
            "taken as the hydrogen-bonded first shell; 12 and 20 test convergence beyond it.", ""]
    (out / "SHELL_REPORT.txt").write_text("\n".join(rep) + "\n")
    print("   generated %d input directories in %s" % (made, out))
    for geom, n in jobs:
        print("     w%02d_%-2s  %3d QM atoms" % (n, geom, 24 + 3 * n))
    print("   shell boundary gap: rank 6 at %.3f A, rank 7 at %.3f A" % (W[5][0], W[6][0]))


RE_E = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
RE_T = re.compile(r"TOTAL RUN TIME:\s*(\d+) days (\d+) hours (\d+) minutes (\d+) seconds")


def measure(run, out):
    run, out = Path(run), Path(out)
    O = ["STAGE 6 - PROBE 1 RESULT: cost of an explicit water shell, and whether it binds",
         "",
         "Bare substrate only. NO designed charges, NO arrangements, NO ddE - by design, so that",
         "nothing here can be read as a result about a design. See STAGE6_PROBE1_20261027.txt.",
         "",
         "The real first shell in the enzyme is mostly PROTEIN (63 non-water against 17 water atoms",
         "within 3.5 A at frame 42436). This is a device for binding the density, not a model of",
         "the active site.",
         "",
         "   run        QM atoms   wall time   SCF cycles   HOMO / eV   bound?"]
    rows = 0
    for geom, n in [("R", s) for s in SHELLS] + [("TS", 6)]:
        p = run / ("w%02d_%s" % (n, geom)) / "job.out"
        if not p.is_file():
            O.append("   w%02d_%-2s     %8d   MISSING" % (n, geom, 24 + 3 * n))
            continue
        s = p.read_text(errors="replace")
        if "****ORCA TERMINATED NORMALLY****" not in s:
            O.append("   w%02d_%-2s     %8d   DID NOT TERMINATE NORMALLY" % (n, geom, 24 + 3 * n))
            continue
        t = RE_T.search(s)
        wall = ("%dd %02dh %02dm %02ds" % tuple(int(x) for x in t.groups())) if t else "?"
        cy = re.findall(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES", s)
        k, homo, occ = s.rfind("ORBITAL ENERGIES"), None, []
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
        O.append("   w%02d_%-2s     %8d   %-11s %10s   %9s   %s"
                 % (n, geom, 24 + 3 * n, wall, cy[-1] if cy else "?",
                    ("%+.4f" % homo) if homo is not None else "-",
                    ("yes" if homo < 0 else "NO") if homo is not None else "?"))
        rows += 1
    O += ["", "READING", "",
          "   Binding is the question that matters. In vacuo the bare dianion's HOMO is +1.2269 eV",
          "   (Diagnostic 10's control). A shell that brings it below zero binds the density by a",
          "   physical mechanism rather than a continuum one. If the shells disagree with each",
          "   other, the shell is not converged and the environment question goes to Phase 1's",
          "   full QM/MM - a decision for the supervisor, not for this probe.",
          "",
          "   Cost is the secondary output: the wall times above size Phase A's environment step.",
          ""]
    out.mkdir(parents=True, exist_ok=True)
    (out / "PROBE1_RESULT.txt").write_text("\n".join(O) + "\n")
    print("\n".join(O))


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "generate":
        generate(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 4 and sys.argv[1] == "measure":
        measure(sys.argv[2], sys.argv[3])
    else:
        stop("usage: cs11_water.py generate <results-repo> <out-dir>\n"
             "       cs11_water.py measure  <run-dir> <out-dir>")
