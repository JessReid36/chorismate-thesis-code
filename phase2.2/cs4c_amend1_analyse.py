#!/usr/bin/env python3
"""
cs4c_amend1_analyse.py - Stage 4c, Amendment 1 (STAGE4C_AMENDMENT1.txt). Removes the pseudopotential's own
charge-independent reactant/TS energy from the pseudopotential levels' ddE, using the zero-charge runs Z_*, and applies
STAGE4C_CRITERIA.txt's comparisons and thresholds unchanged. The original verdict (STAGE4C_REPORT.txt) is not touched.
  ddE_corr(L, E) = ddE(L, E) - ddE(Z_L_set(E)),  L = L0, L2, L3;  L1 (bare spheres, LJ excluded) is unchanged.
Reads job.out, or job.out.xz once committed. Python 3.6, standard library. Writes STAGE4C_AMENDMENT1_REPORT.txt.
USAGE  python3 cs4c_amend1_analyse.py [cs_stage4c folder]
"""
import lzma, os, re, sys
D = sys.argv[1] if len(sys.argv) > 1 else "/home/18660916/system_development/phase2.2/cs_stage4c"
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
    return ""
def parse(run):
    t = text(run)
    if "ORCA TERMINATED NORMALLY" not in t:
        return None
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    if not e:
        return None
    mm = re.findall(r"FINAL SINGLE POINT ENERGY \(MM\)\s+(-?\d+\.\d+)", t)
    nel = re.findall(r"Number of Electrons\s+NEL\s+\.+\s+(\d+)", t)
    nbf = re.findall(r"Number of basis functions\s+\.+\s+(\d+)", t)
    i = t.rfind("ORBITAL ENERGIES"); orbs = []
    for l in t[i:i + 80000].splitlines()[4:]:
        m = re.match(r"\s*\d+\s+(\d\.\d+)\s+-?\d+\.\d+\s+(-?\d+\.\d+)\s*$", l)
        if m: orbs.append((float(m.group(1)), float(m.group(2))))
        elif orbs: break
    gap = (min(o[1] for o in orbs if o[0] < 0.5) - max(o[1] for o in orbs if o[0] > 0.5)) if orbs else None
    return {"E": float(e[-1]), "mm": float(mm[-1]) if mm else None, "nel": int(nel[-1]) if nel else None,
            "nbf": int(nbf[0]) if nbf else None, "necp": len(re.findall(r"Atom\s+\d+Ne>\s+ECP group", t)), "gap": gap}
NSPH = {}
for l in open(os.path.join(D, "arrangements.tsv")):
    if not l.startswith("#") and l.strip():
        NSPH[l.split("\t")[0]] = NSPH.get(l.split("\t")[0], 0) + 1
FIRST = {}
for l in open(os.path.join(D, "arrangements.tsv")):
    if l.startswith("# E") and "first-order" in l:
        FIRST[l.split()[1]] = float(re.findall(r"=\s*([-+]?\d+\.\d+)", l)[0])
rep = []; say = rep.append
say("Stage 4c Amendment 1 report (cs4c_amend1_analyse.py) - STAGE4C_AMENDMENT1.txt; thresholds from STAGE4C_CRITERIA.txt")
R = {}
for b in ("svp", "svpd", "tzvpd"):
    for g in ("R", "TS"): R["bare", b, g] = parse("bare_%s_%s" % (b, g))
for L in LEV:
    for e in ENVS:
        for g in ("R", "TS"): R[L, e, g] = parse("%s_%s_%s" % (L, e, g))
Z = {}
for L in ("L0", "L2", "L3"):
    for s in SRC:
        for g in ("R", "TS"): Z[L, s, g] = parse("Z_%s_%s_%s" % (L, s, g))
missing = ["_".join(k) for k, v in R.items() if v is None] + ["Z_" + "_".join(k) for k, v in Z.items() if v is None]
def finish():
    open(os.path.join(D, "STAGE4C_AMENDMENT1_REPORT.txt"), "w").write("\n".join(rep) + "\n"); print("\n".join(rep)); sys.exit(0)
if missing:
    say("  %d run(s) missing or not terminated normally: %s -> not assessable" % (len(missing), " ".join(missing[:12])))
    finish()
nbf_ref = {b: {R[L, e, g]["nbf"] for L in LEV if LEV[L] == b for e in ENVS for g in ("R", "TS")} for b in ("svp", "svpd", "tzvpd")}
ok = all(v["nel"] == 118 for v in Z.values()) and all(len(s) == 1 for s in nbf_ref.values()) and \
    all({Z[L, s, g]["nbf"]} == nbf_ref[LEV[L]] and Z[L, s, g]["necp"] == NSPH[SRC[s]] for (L, s, g) in Z)
say("check 0 for the 36 zero-charge runs (all terminated; 118 electrons; basis functions as in the original runs of that"
    " basis; one Ne-type ECP centre per sphere): %s" % ("PASS" if ok else "FAIL"))
if not ok:
    finish()
def dde(L, e):
    b = LEV[L]
    return ((R[L, e, "TS"]["E"] - R["bare", b, "TS"]["E"]) - (R[L, e, "R"]["E"] - R["bare", b, "R"]["E"])) * H
def ddz(L, s):
    b = LEV[L]
    return ((Z[L, s, "TS"]["E"] - R["bare", b, "TS"]["E"]) - (Z[L, s, "R"]["E"] - R["bare", b, "R"]["E"])) * H
X = {}
for e in ENVS:
    X["L1", e] = dde("L1", e)
    for L in ("L0", "L2", "L3"):
        X[L, e] = dde(L, e) - ddz(L, SETOF[e])
X.update({("first", e): FIRST[e] for e in ENVS})
say("\nthe pseudopotential's own reactant/TS term, ddE at zero charge (kcal/mol), beside the sphere LJ term at L1:")
for s in SRC:
    say("  %s (%d sphere%s)  L0 %+7.3f  L2 %+7.3f  L3 %+7.3f   dE_LJ %+7.3f" % (s, NSPH[SRC[s]], "" if NSPH[SRC[s]] == 1 else "s",
        ddz("L0", s), ddz("L2", s), ddz("L3", s), (R["L1", SRC[s], "TS"]["mm"] - R["L1", SRC[s], "R"]["mm"]) * H))
say("\ncorrected ddE by arrangement (kcal/mol; negative lowers the barrier; L1 unchanged):")
say("  %-4s %10s %10s %10s %10s %10s" % ("", "first-ord", "L1", "L0 corr", "L2 corr", "L3 corr"))
for e in ENVS:
    say("  %-4s %+10.3f %+10.3f %+10.3f %+10.3f %+10.3f" % (e, X["first", e], X["L1", e], X["L0", e], X["L2", e], X["L3", e]))
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
def compare(test, ref, envs):
    x = [X[test, e] for e in envs]; y = [X[ref, e] for e in envs]
    k = sum(a * b for a, b in zip(x, y)) / sum(a * a for a in x)
    res = (sum((b - k * a) ** 2 for a, b in zip(x, y)) / len(x)) ** 0.5
    rho = spearman(x, y)
    signs = all((a > 0) == (b > 0) for a, b in zip(x, y) if abs(b) >= 0.5)
    return k, res, rho, signs, max(abs(b - a) for a, b in zip(x, y))
def verdict(k, res, rho, signs):
    if not (rho >= 0.95 and signs):
        return "FAIL"
    return "PASS" if (res <= 0.5 and 0.8 <= k <= 1.25) else "RANKING ONLY"
def line(test, ref, envs, role):
    k, res, rho, signs, md = compare(test, ref, envs); v = verdict(k, res, rho, signs)
    say("  %-5s vs %-3s k %.3f  residual %.3f  rho %.3f  signs %s  max |diff| %.3f -> %s  (%s)" % (
        test, ref, k, res, rho, "agree" if signs else "DIFFER", md, v, role))
    return v
say("\ncomparisons on the corrected values (fit ref = k x test through the origin; thresholds unchanged):")
vA = line("L1", "L3", ENVS, "PRIMARY, corrected")
vF = line("L2", "L3", ENVS, "fallback level, corrected")
vR = line("L1", "L0", ENVS, "REPRESENTATION CHECK: pseudopotential spheres against production at the same basis")
line("L1", "L2", ENVS, "informative: production against the proposed check level")
line("L0", "L2", ENVS, "informative: diffuse functions at fixed representation")
line("first", "L1", ENVS, "informative: Tier 2, first order against full L1")
noE5 = [e for e in ENVS if e != "E5"]
say("informative, chosen after E5's def2-TZVPD reactant was seen to change electronic state (STAGE4C_EXPLORATORY.txt):")
line("L1", "L3", noE5, "without E5")
line("L1", "L2", noE5, "without E5")
say("\nzero-charge runs, HOMO-LUMO gap (eV): " + ", ".join("%s_%s_%s %.2f" % (k[0], k[1], k[2], v["gap"]) for k, v in sorted(Z.items())))
say("\nS4c verdict after Amendment 1 (L1 against L3, corrected): %s" % vA)
say("representation check (L1 against L0, corrected): %s -> %s" % (vR,
    "the corrected pseudopotential-sphere single point is accepted as the diffuse-basis check for final designs"
    if vR == "PASS" else "the corrected pseudopotential-sphere single point is NOT accepted as a check tool; the basis "
    "uncertainty of final designs is stated from Stage 4c instead"))
say("decision rule applied (STAGE4C_AMENDMENT1.txt): production design and relaxed validation stay at def2-SVP with bare"
    " LJ spheres; see the note for the final-design check and the enzyme reference.")
finish()
