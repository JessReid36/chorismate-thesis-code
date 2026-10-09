#!/usr/bin/env python3
"""
cs4c_explore.py - EXPLORATORY analysis of Stage 4c (not pre-registered; written 9 October 2026 after Stage 4c's
results were seen). It changes no Stage 4c verdict; it locates the sources of the L1-L3 differences to inform
Amendment 1 and the validation protocol (see STAGE4C_EXPLORATORY.txt).
Reads each run's job.out, or job.out.xz (committed outputs). Python 3.6, standard library.
USAGE  python3 cs4c_explore.py <cs_stage4c folder>  -> writes <folder>/exploratory/STAGE4C_EXPLORATORY_REPORT.txt
"""
import lzma, os, re, sys
D = sys.argv[1]
H = 627.5094740631
ENVS = ["E%d" % k for k in range(1, 9)]
LEV = {"L1": "svp", "L0": "svp", "L2": "svpd", "L3": "tzvpd"}
N = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
def text(run):
    p = os.path.join(D, run, "job.out")
    if os.path.exists(p):
        return open(p, errors="replace").read()
    if os.path.exists(p + ".xz"):
        return lzma.open(p + ".xz", "rt", errors="replace").read()
    raise SystemExit("missing output: " + run)
def parse(run):
    t = text(run)
    assert "ORCA TERMINATED NORMALLY" in t, run
    e = float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)[-1])
    i = t.rfind("ORBITAL ENERGIES"); orbs = []
    for l in t[i:i + 80000].splitlines()[4:]:
        m = re.match(r"\s*\d+\s+(\d\.\d+)\s+-?\d+\.\d+\s+(-?\d+\.\d+)\s*$", l)
        if m: orbs.append((float(m.group(1)), float(m.group(2))))
        elif orbs: break
    homo = max(o[1] for o in orbs if o[0] > 0.5); lumo = min(o[1] for o in orbs if o[0] < 0.5)
    cyc = re.findall(r"SCF CONVERGED AFTER\s+(\d+)", t)
    i = t.rfind("LOEWDIN ATOMIC CHARGES"); q = {}
    for l in t[i:i + 4000].splitlines()[2:]:
        m = re.match(r"\s*(\d+)\s+\S+\s*:\s*(-?\d+\.\d+)", l)
        if not m: break
        q[int(m.group(1))] = float(m.group(2))
    return {"E": e, "homo": homo, "lumo": lumo, "cyc": cyc[-1] if cyc else "?", "low": [q[k] for k in range(24)]}
FIRST = {}
for l in open(os.path.join(D, "arrangements.tsv")):
    if l.startswith("# E") and "first-order" in l:
        FIRST[l.split()[1]] = float(re.findall(r"=\s*([-+]?\d+\.\d+)", l)[0])
R = {}
for b in ("svp", "svpd", "tzvpd"):
    for g in ("R", "TS"): R["bare", b, g] = parse("bare_%s_%s" % (b, g))
for L in LEV:
    for e in ENVS:
        for g in ("R", "TS"): R[L, e, g] = parse("%s_%s_%s" % (L, e, g))
def dde(L, e):
    b = LEV[L]
    return ((R[L, e, "TS"]["E"] - R["bare", b, "TS"]["E"]) - (R[L, e, "R"]["E"] - R["bare", b, "R"]["E"])) * H
X = {(L, e): dde(L, e) for L in LEV for e in ENVS}
def ranks(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]: j += 1
        for t in range(i, j + 1): r[o[t]] = (i + j) / 2.0 + 1
        i = j + 1
    return r
def spearman(x, y):
    rx, ry = ranks(x), ranks(y); n = len(x); mx = sum(rx) / n; my = sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    return num / (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
def comp(test, ref, envs):
    x = [X[test, e] for e in envs]; y = [X[ref, e] for e in envs]
    k = sum(a * b for a, b in zip(x, y)) / sum(a * a for a in x)
    res = (sum((b - k * a) ** 2 for a, b in zip(x, y)) / len(x)) ** 0.5
    md = max(abs(b - a) for a, b in zip(x, y)); rms = (sum((b - a) ** 2 for a, b in zip(x, y)) / len(x)) ** 0.5
    return k, res, spearman(x, y), md, rms
out = []; say = out.append
say("Stage 4c EXPLORATORY report (cs4c_explore.py) - not pre-registered; Stage 4c's verdict stands as recorded.")
say("ddE in kcal/mol (negative lowers the barrier); orbital energies in eV; Loewdin charges in e.")
say("\n1. ddE recomputed from the outputs, and the bare barrier at the fixed geometries")
for e in ENVS:
    say("   %s  first-order %+8.3f  L1 %+8.3f  L0 %+8.3f  L2 %+8.3f  L3 %+8.3f" % (e, FIRST[e], X["L1", e], X["L0", e], X["L2", e], X["L3", e]))
say("   bare barrier: " + ", ".join("%s %.3f" % (b, (R["bare", b, "TS"]["E"] - R["bare", b, "R"]["E"]) * H) for b in ("svp", "svpd", "tzvpd")))
say("\n2. comparisons between levels (k: fit ref = k x test through the origin; residual RMS after scaling; rho Spearman)")
for lab, envs in (("all eight", ENVS), ("without E5", [e for e in ENVS if e != "E5"])):
    say("   " + lab)
    for t, r in (("L1", "L3"), ("L1", "L2"), ("L2", "L3"), ("L0", "L2"), ("L1", "L0")):
        k, res, rho, md, rms = comp(t, r, envs)
        say("      %s vs %s  k %.3f  residual %.3f  rho %.3f  max|diff| %.2f  RMS diff %.2f" % (t, r, k, res, rho, md, rms))
say("\n3. which state moves from def2-SVPD to def2-TZVPD: shift of [E(spheres) - E(bare)] for R and TS")
for e in ENVS:
    sh = [((R["L3", e, g]["E"] - R["bare", "tzvpd", g]["E"]) - (R["L2", e, g]["E"] - R["bare", "svpd", g]["E"])) * H for g in ("R", "TS")]
    say("   %s  R %+6.2f  TS %+6.2f  (R minus TS %+6.2f)" % (e, sh[0], sh[1], sh[0] - sh[1]))
say("\n4. diffuse functions at fixed representation (L2 - L0, pseudopotential spheres in both)")
say("   " + "  ".join("%s %+.2f" % (e, X["L2", e] - X["L0", e]) for e in ENVS))
cat = [e for e in ENVS if X["L0", e] < 0]
say("   catalytic arrangements (L0 < 0) made more catalytic by diffuse functions: %d of %d" % (
    sum(1 for e in cat if X["L2", e] < X["L0", e]), len(cat)))
say("\n5. polarisation at production level: L1 minus the first-order prediction")
say("   " + "  ".join("%s %+.2f (%+.0f%%)" % (e, X["L1", e] - FIRST[e], 100 * (X["L1", e] - FIRST[e]) / abs(FIRST[e])) for e in ENVS))
a = 4 * X["L1", "E2"] - X["L1", "E1"]; b = 2 * (X["L1", "E1"] - 2 * X["L1", "E2"])
say("   site 117 from E1 (q = +1) and E2 (q = +0.5) at L1: a %.3f (A matrix %.3f), b %.3f" % (a, FIRST["E1"], b))
say("\n6. E5 and E7 (same three sites, signs reversed): odd part (E5 - E7)/2, even part (E5 + E7)/2")
for L in ("L1", "L0", "L2", "L3"):
    say("   %s  odd %+8.3f  even %+8.3f" % (L, (X[L, "E5"] - X[L, "E7"]) / 2, (X[L, "E5"] + X[L, "E7"]) / 2))
say("   first-order odd part %+8.3f" % FIRST["E5"])
say("\n7. electronic diagnostics: SCF cycles, HOMO, LUMO, gap (eV)")
runs = [("bare", b, g) for b in ("svp", "svpd", "tzvpd") for g in ("R", "TS")] + [(L, e, g) for L in ("L1", "L0", "L2", "L3") for e in ENVS for g in ("R", "TS")]
for k in runs:
    r = R[k]
    say("   %-14s cycles %-3s HOMO %7.2f  LUMO %7.2f  gap %5.2f" % ("_".join(k), r["cyc"], r["homo"], r["lumo"], r["lumo"] - r["homo"]))
say("\n8. Loewdin charge shifts against the bare substrate at the same basis (largest five, and O3)")
for L, e, g in (("L1", "E5", "R"), ("L2", "E5", "R"), ("L3", "E5", "R"), ("L1", "E5", "TS"), ("L2", "E5", "TS"), ("L3", "E5", "TS"),
                ("L2", "E7", "R"), ("L3", "E7", "R")):
    a = R[L, e, g]["low"]; b = R["bare", LEV[L], g]["low"]; d = [x - y for x, y in zip(a, b)]
    top = sorted(range(24), key=lambda k: -abs(d[k]))[:5]
    say("   %s_%s_%-2s sum %+.3f; O3 %+.3f; largest: %s" % (L, e, g, sum(a), d[7], ", ".join("%s %+.3f" % (N[k], d[k]) for k in top)))
os.makedirs(os.path.join(D, "exploratory"), exist_ok=True)
open(os.path.join(D, "exploratory", "STAGE4C_EXPLORATORY_REPORT.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
