#!/usr/bin/env python3
"""
STAGE 6 - CROSS-CHECK between Diagnostic 12 and Diagnostic 13.

D13's q = 1.00 arm and D12's main arm measure the SAME quantity - ddE for a +1 point charge, bare
and behind a Pauli wall, at R and TS - at the eight sites D13 shares with D12.  The two input sets
were generated independently by two harnesses written a day apart.  Nothing in either diagnostic
compares them, so this does.

It is a harness check, not physics.  Agreement proves only that the two generators placed the same
charge at the same place and that ORCA returned the same number twice; disagreement proves one of
them is wrong, and says which quantity.

  cs_crosscheck.py <d12-run-dir> <d13-run-dir> <results-repo> <out-dir>
"""
import sys, os, re, math, hashlib
from pathlib import Path

TOL_EH = 1e-6          # the two controls are the same calculation; they must agree to this
TOL_KCAL = 0.01        # ddE agreement, the same tolerance D12 and D13 use for their own gates

# Constants and helpers below are D13's own code, extracted verbatim from cs13_qscale.py
# so that the cross-check cannot disagree with the diagnostics about what it is reading.
HARTREE = 627.5094740631
GEOM = ("R", "TS")
SITES_SHA = "d2a750077d79e34dfecf448c30c2899ae0cea98f1edc5de3d10361b7c357cfc7"
LBL = ("C1", "H1", "H2", "C2", "C3", "O1", "O2", "O3", "C4", "H3", "C5", "H4",
       "C6", "C7", "H5", "C8", "H6", "C9", "H7", "O4", "H8", "C10", "O5", "O6")
FRAME_DIR = "phase2.2/cs_stage6/frame_42436"
FRAME_SHA = {"reactant_qm.xyz":
             "6e8f41c7da3bcf4ed9f5c650dc9fcf12c2a36e5af7a900e164892c276497aa88",
             "transition_state_qm.xyz":
             "2babeb9000ca5220a908a70eaec80ec4a64358c8f5677e9a16cb104e0dffc5bf"}
RE_E = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
RE_NORM = "****ORCA TERMINATED NORMALLY****"
RE_LW = re.compile(r"^\s*(\d+)\s+([A-Za-z]{1,2})\s*:\s*(-?\d+\.\d+)")
SITES = ("2009", "1981", "1464", "1690", "1710", "1679", "2330", "4")
GATE_SITE = "1690"
D2_C3_BARE = {0.25: -1.8935, 0.50: -3.9781, 1.00: -3.8579}
D11_DREF_WALL = -7.8221



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

def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    return [(f[0], float(f[1]), float(f[2]), float(f[3]))
            for f in (l.split() for l in L[2:2 + n])]

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


def manifest_check(run, label):
    """GATE 0 on each run directory, separately."""
    cand = [m for m in Path(run).rglob("MANIFEST.sha256") if ".git" not in m.parts]
    if len(cand) != 1:
        stop("%s: expected exactly one MANIFEST.sha256, found %d" % (label, len(cand)))
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
        stop("%s: GATE 0 failed, %d files do not match the manifest:\n   %s"
             % (label, len(bad), "\n   ".join(bad[:8])))
    return base


def energies(base, names, label):
    D = {}
    for n in names:
        for g in GEOM:
            p = base / ("%s_%s" % (n, g)) / "job.out"
            if not p.is_file():
                stop("%s: missing output %s" % (label, p))
            r = read_out(p)
            if not r["norm"] or r["E"] is None:
                stop("%s: %s_%s did not terminate normally" % (label, n, g))
            D[(n, g)] = r["E"]
    return D


def main(d12, d13, res, out):
    d12, d13, res, out = Path(d12), Path(d13), Path(res), Path(out)
    for p in (out,):
        for r in (d12, d13):
            if str(p.resolve()) == str(r.resolve()) or \
               str(p.resolve()).startswith(str(r.resolve()) + os.sep):
                stop("refusing to write inside a run directory being read: %s" % p)
    b12 = manifest_check(d12, "D12")
    b13 = manifest_check(d13, "D13")

    # the eight sites D13 uses are a subset of D12's thirty-six; SITES here is D13's list
    n12 = ["ctl_ma"] + ["s%s_%s_p1_ma" % (s, a) for s in SITES for a in ("bare", "wall")]
    n13 = ["ctl_ma"] + ["s%s_q100_%s" % (s, a) for s in SITES for a in ("bare", "wall")]
    E12 = energies(b12, n12, "D12")
    E13 = energies(b13, n13, "D13")

    O = ["STAGE 6 - CROSS-CHECK: DIAGNOSTIC 12 AGAINST DIAGNOSTIC 13", "",
         "D13's q = 1.00 arm repeats D12's main arm at the eight sites they share, from inputs",
         "generated independently by two harnesses. Neither diagnostic compares them. This does.",
         "",
         "A harness check, not physics: agreement shows the two generators placed the same charge",
         "in the same place and ORCA returned the same number twice.", ""]

    O += ["THE CONTROL", ""]
    ok_all = True
    for g in GEOM:
        a, b = E12[("ctl_ma", g)], E13[("ctl_ma", g)]
        ok = abs(a - b) <= TOL_EH
        ok_all &= ok
        O.append("   %-2s  D12 %.9f   D13 %.9f   diff %+.2e Eh   %s"
                 % (g, a, b, a - b, "OK" if ok else "FAIL"))
    O.append("")

    def ddE(E, name):
        return ((E[(name, "TS")] - E[(name, "R")])
                - (E[("ctl_ma", "TS")] - E[("ctl_ma", "R")])) * HARTREE

    O += ["ddE AT THE EIGHT SHARED SITES, q = +1", "",
          "   site    arm     D12 ddE     D13 ddE      diff       verdict"]
    worst = 0.0
    for s in SITES:
        for arm in ("bare", "wall"):
            a = ddE(E12, "s%s_%s_p1_ma" % (s, arm))
            b = ddE(E13, "s%s_q100_%s" % (s, arm))
            d = a - b
            worst = max(worst, abs(d))
            ok = abs(d) <= TOL_KCAL
            ok_all &= ok
            O.append("   %-6s  %-5s  %+9.4f   %+9.4f   %+9.6f   %s"
                     % (s, arm, a, b, d, "OK" if ok else "FAIL"))
    O += ["", "   largest disagreement %.6f kcal/mol against a tolerance of %.2f"
          % (worst, TOL_KCAL), ""]

    O += ["VERDICT", ""]
    if ok_all:
        O += ["   THE TWO HARNESSES AGREE.",
              "",
              "   D12 and D13 return the same ddE at every shared site and arm, and the same",
              "   control energies. Two independently generated input sets, 114 and 282 runs,",
              "   reproduce each other. Neither diagnostic's reading rests on this, but any",
              "   conclusion drawn from both of them together now has it.", ""]
    else:
        O += ["   THE TWO HARNESSES DISAGREE.",
              "",
              "   The same charge at the same position returned different energies, or the two",
              "   controls differ. One of the two generators is wrong and nothing that compares",
              "   D12 with D13 may be read until it is found. Each diagnostic's own gates are",
              "   unaffected and stand on their own.", ""]
    out.mkdir(parents=True, exist_ok=True)
    (out / "CROSSCHECK_D12_D13.txt").write_text("\n".join(O) + "\n")
    print("\n".join(O))


if __name__ == "__main__":
    if len(sys.argv) == 5:
        main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        stop("usage: cs_crosscheck.py <d12-run-dir> <d13-run-dir> <results-repo> <out-dir>")
