#!/usr/bin/env python3
"""
cs4c_amend1_explore.py - EXPLORATORY analysis of Stage 4c Amendment 1 (not pre-registered; written 10 October 2026
after STAGE4C_AMENDMENT1_REPORT.txt was seen). It changes no verdict. It looks at where, after the zero-charge
correction, diffuse functions change the barrier change, and checks the reading carried forward from
STAGE4C_EXPLORATORY.txt item 3 (see STAGE4C_AMENDMENT1_EXPLORATORY.txt).
Reads each run's job.out, or job.out.xz (committed outputs). Python 3.6, standard library.
USAGE  python3 cs4c_amend1_explore.py <cs_stage4c folder>
       -> writes <folder>/exploratory/STAGE4C_AMENDMENT1_EXPLORATORY_REPORT.txt
"""
import lzma, os, re, sys
D = sys.argv[1]
H = 627.5094740631
ENVS = ["E%d" % k for k in range(1, 9)]
LEV = {"L1": "svp", "L0": "svp", "L2": "svpd", "L3": "tzvpd"}
SETOF = {"E1": "S1", "E2": "S1", "E3": "S3", "E4": "S4", "E5": "S5", "E6": "S6", "E7": "S5", "E8": "S8"}
SRC = {"S1": "E1", "S3": "E3", "S4": "E4", "S5": "E5", "S6": "E6", "S8": "E8"}
def text(run):
    p = os.path.join(D, run, "job.out")
    if os.path.exists(p):
        return open(p, errors="replace").read()
    if os.path.exists(p + ".xz"):
        return lzma.open(p + ".xz", "rt", errors="replace").read()
    raise SystemExit("missing output: " + run)
def parse(run):
    t = text(run)
    if "ORCA TERMINATED NORMALLY" not in t:
        raise SystemExit("not terminated normally: " + run)
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    mm = re.findall(r"FINAL SINGLE POINT ENERGY \(MM\)\s+(-?\d+\.\d+)", t)
    return {"E": float(e[-1]), "mm": float(mm[-1]) if mm else None}
E = {}
for b in ("svp", "svpd", "tzvpd"):
    for g in ("R", "TS"): E["bare", b, g] = parse("bare_%s_%s" % (b, g))["E"]
MM = {}
for L in LEV:
    for e in ENVS:
        for g in ("R", "TS"):
            p = parse("%s_%s_%s" % (L, e, g)); E[L, e, g] = p["E"]
            if L == "L1": MM[e, g] = p["mm"]
for L in ("L0", "L2", "L3"):
    for s in SRC:
        for g in ("R", "TS"): E["Z" + L, s, g] = parse("Z_%s_%s_%s" % (L, s, g))["E"]
def dde(L, e):
    b = LEV[L]
    return ((E[L, e, "TS"] - E["bare", b, "TS"]) - (E[L, e, "R"] - E["bare", b, "R"])) * H
def ddz(L, s):
    b = LEV[L]
    return ((E["Z" + L, s, "TS"] - E["bare", b, "TS"]) - (E["Z" + L, s, "R"] - E["bare", b, "R"])) * H
RAW = {(L, e): dde(L, e) for L in LEV for e in ENVS}
X = {("L1", e): RAW["L1", e] for e in ENVS}
for L in ("L0", "L2", "L3"):
    for e in ENVS: X[L, e] = RAW[L, e] - ddz(L, SETOF[e])
LJ = {e: (MM[e, "TS"] - MM[e, "R"]) * H for e in ENVS}
NEAR = {}
for l in open(os.path.join(D, "arrangements.tsv")):
    if l.startswith("#") or not l.strip(): continue
    f = l.rstrip("\n").split("\t")
    NEAR.setdefault(f[0], []).append((f[1], f[2], f[4]))
out = []; say = out.append
def rms(v): return (sum(x * x for x in v) / len(v)) ** 0.5
say("Stage 4c Amendment 1 - EXPLORATORY report (cs4c_amend1_explore.py; not pre-registered; verdicts unchanged)")
say("Corrected ddE as in cs4c_amend1_analyse.py: ddE_corr(L, E) = ddE(L, E) - ddE(zero-charge spheres of E's site set).")
say("\n1. Diffuse functions at fixed representation (pseudopotential spheres; L2 against L0), kcal/mol. 'change' > 0 means")
say("   the barrier change moves up (a catalytic arrangement less catalytic, an anti-catalytic one more anti-catalytic).")
say("   %-4s %9s %9s %8s %7s   %9s %11s" % ("", "L0 corr", "L2 corr", "change", "% |L0|", "uncorr.", "z(L2)-z(L0)"))
nless = 0; ncat = 0
for e in ENVS:
    c = X["L2", e] - X["L0", e]; u = RAW["L2", e] - RAW["L0", e]; dz = ddz("L2", SETOF[e]) - ddz("L0", SETOF[e])
    say("   %-4s %+9.3f %+9.3f %+8.3f %+6.1f%%   %+9.3f %+11.3f" % (e, X["L0", e], X["L2", e], c, 100 * c / abs(X["L0", e]), u, dz))
    if X["L0", e] < 0:
        ncat += 1; nless += c > 0
say("   catalytic arrangements (L0 corr < 0) made less catalytic by diffuse functions: corrected %d of %d; uncorrected %d of %d"
    % (nless, ncat, sum(RAW["L2", e] - RAW["L0", e] > 0 for e in ENVS if X["L0", e] < 0), ncat))
k = sum(X["L0", e] * X["L2", e] for e in ENVS) / sum(X["L0", e] ** 2 for e in ENVS)
say("   common scale k (L2 corr = k x L0 corr): %.3f" % k)
say("\n2. Production against the check level (L2 corr - L1), kcal/mol, and as a percentage of |L1|:")
for e in ENVS:
    c = X["L2", e] - X["L1", e]
    say("   %-4s L1 %+9.3f  L2 corr %+9.3f  diff %+7.3f  (%+5.1f%%)" % (e, X["L1", e], X["L2", e], c, 100 * c / abs(X["L1", e])))
say("\n3. Grouping by the committed nearest-atom column of arrangements.tsv (a rule chosen after the numbers were seen):")
say("   group C = some sphere's nearest substrate atom is a C10 carboxylate oxygen (O5 or O6); group N = none is.")
G = {"C": [], "N": []}
for e in ENVS:
    g = "C" if any(a.split()[0] in ("O5", "O6") for _, _, a in NEAR[e]) else "N"
    G[g].append(e)
    say("   %-4s group %s   spheres: %s" % (e, g, "; ".join("site %s q %s nearest %s A" % s for s in NEAR[e])))
for g in ("C", "N"):
    say("   group %s (%s): RMS (L2 corr - L0 corr) %.3f; RMS (L2 corr - L1) %.3f; max |L2 corr - L0 corr| %.3f" % (
        g, " ".join(G[g]), rms([X["L2", e] - X["L0", e] for e in G[g]]), rms([X["L2", e] - X["L1", e] for e in G[g]]),
        max(abs(X["L2", e] - X["L0", e]) for e in G[g])))
say("\n4. Order from most to least barrier-lowering at each level, and the gaps between neighbours (kcal/mol):")
for L in ("L1", "L0", "L2", "L3"):
    o = sorted(ENVS, key=lambda e: X[L, e])
    say("   %-8s %s" % (L if L == "L1" else L + " corr", "  ".join(
        "%s%s" % (o[i], "" if i == len(o) - 1 else " (%.2f)" % (X[L, o[i + 1]] - X[L, o[i]])) for i in range(len(o)))))
say("\n5. The pseudopotential's zero-charge term against the sphere LJ term (dE_LJ from L1's MM energy), kcal/mol:")
same = 0; tot = 0
for s in SRC:
    zs = [ddz(L, s) for L in ("L0", "L2", "L3")]; lj = LJ[SRC[s]]
    same += sum((z > 0) == (lj > 0) for z in zs); tot += 3
    say("   %s  z: L0 %+7.3f  L2 %+7.3f  L3 %+7.3f   dE_LJ %+7.3f   z/dE_LJ %s" % (
        s, zs[0], zs[1], zs[2], lj, " ".join("%.2f" % (z / lj) for z in zs)))
say("   same sign as dE_LJ: %d of %d" % (same, tot))
def fit(t, r, envs, tab):
    x = [tab[t, e] for e in envs]; y = [tab[r, e] for e in envs]
    kk = sum(a * b for a, b in zip(x, y)) / sum(a * a for a in x)
    return kk, rms([b - kk * a for a, b in zip(x, y)]), max(abs(b - a) for a, b in zip(x, y))
say("\n6. Effect of the correction on the representation check (L1 against L0; fit L0 = k x L1):")
say("   uncorrected  k %.3f  residual %.3f  max |diff| %.3f" % fit("L1", "L0", ENVS, RAW))
say("   corrected    k %.3f  residual %.3f  max |diff| %.3f" % fit("L1", "L0", ENVS, X))
os.makedirs(os.path.join(D, "exploratory"), exist_ok=True)
open(os.path.join(D, "exploratory", "STAGE4C_AMENDMENT1_EXPLORATORY_REPORT.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
