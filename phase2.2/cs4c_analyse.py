#!/usr/bin/env python3
"""
cs4c_analyse.py - Stage 4c: applies STAGE4C_CRITERIA.txt to phase2.2/cs_stage4c/. Python 3.6, standard library.
ddE(arrangement, level) = [E_TS(spheres) - E_TS(bare)] - [E_R(spheres) - E_R(bare)], kcal/mol, from the QM energy
(ORCA's plain FINAL SINGLE POINT ENERGY: in the QM/MM runs this is the QM part, without the LJ term). Writes
STAGE4C_REPORT.txt.  USAGE  python3 cs4c_analyse.py [cs_stage4c folder]
"""
import os, re, sys
D = sys.argv[1] if len(sys.argv) > 1 else "/home/18660916/system_development/phase2.2/cs_stage4c"
H = 627.5094740631
ENVS = ["E%d" % k for k in range(1, 9)]
def parse(name):
    p = os.path.join(D, name, "job.out")
    t = open(p, errors="replace").read() if os.path.exists(p) else ""
    if "ORCA TERMINATED NORMALLY" not in t:
        return None
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    mm = re.findall(r"FINAL SINGLE POINT ENERGY \(MM\)\s+(-?\d+\.\d+)", t)
    nel = re.findall(r"Number of Electrons\s+NEL\s+\.+\s+(\d+)", t)
    nbf = re.findall(r"Number of basis functions\s+\.+\s+(\d+)", t)
    if not e:
        return None
    return {"E": float(e[-1]), "mm": float(mm[-1]) if mm else None, "nel": int(nel[-1]) if nel else None,
            "nbf": int(nbf[0]) if nbf else None, "necp": len(re.findall(r"Atom\s+\d+Ne>\s+ECP group", t))}
NSPH, FIRST = {}, {}
for l in open(os.path.join(D, "arrangements.tsv")):
    if l.startswith("# E") and "first-order" in l:
        FIRST[l.split()[1]] = float(re.findall(r"=\s*([-+]?\d+\.\d+)", l)[0])
    elif not l.startswith("#") and l.strip():
        NSPH[l.split("\t")[0]] = NSPH.get(l.split("\t")[0], 0) + 1
rep = []; say = rep.append
say("Stage 4c report (cs4c_analyse.py) - criteria from STAGE4C_CRITERIA.txt")
LEVELS = {"L1": ("svp", "QM/MM, def2-SVP, bare spheres (production)"), "L0": ("svp", "def2-SVP, Ne-type pseudopotential spheres"),
          "L2": ("svpd", "def2-SVPD, pseudopotential spheres"), "L3": ("tzvpd", "def2-TZVPD, pseudopotential spheres (reference)")}
R = {}
for b in ("svp", "svpd", "tzvpd"):
    for g in ("R", "TS"): R["bare", b, g] = parse("bare_%s_%s" % (b, g))
for L in LEVELS:
    for e in ENVS:
        for g in ("R", "TS"): R[L, e, g] = parse("%s_%s_%s" % (L, e, g))
missing = [k for k, v in R.items() if v is None]
if missing:
    say("  %d run(s) missing or not terminated normally: %s -> not assessable" % (len(missing), " ".join("_".join(k) for k in missing[:12])))
    open(os.path.join(D, "STAGE4C_REPORT.txt"), "w").write("\n".join(rep) + "\n"); print("\n".join(rep)); sys.exit(0)
# check 0
nel = {v["nel"] for v in R.values()}
nbf = {}
for k, v in R.items():
    b = k[1] if k[0] == "bare" else LEVELS[k[0]][0]
    nbf.setdefault(b, set()).add(v["nbf"])
ecp_ok = all(R[L, e, g]["necp"] == NSPH[e] for L in ("L0", "L2", "L3") for e in ENVS for g in ("R", "TS")) and \
    all(R[L, e, g]["necp"] == 0 for L in ("L1",) for e in ENVS for g in ("R", "TS")) and all(R["bare", b, g]["necp"] == 0 for b in ("svp", "svpd", "tzvpd") for g in ("R", "TS"))
ck0 = nel == {118} and all(len(s) == 1 for s in nbf.values()) and ecp_ok
say("check 0 (70 runs; 118 electrons; one basis-function count per basis, so spheres add none; one Ne-type ECP per "
    "sphere in L0/L2/L3, none elsewhere): %s  [NEL %s; basis functions %s; ECP centres %s]" % (
        "PASS" if ck0 else "FAIL", sorted(nel), {b: sorted(s) for b, s in nbf.items()}, "as expected" if ecp_ok else "MISMATCH"))
def dde(L, e):
    b = LEVELS[L][0]
    return ((R[L, e, "TS"]["E"] - R["bare", b, "TS"]["E"]) - (R[L, e, "R"]["E"] - R["bare", b, "R"]["E"])) * H
X = {L: [dde(L, e) for e in ENVS] for L in LEVELS}
X["first"] = [FIRST[e] for e in ENVS]
lj = [(R["L1", e, "TS"]["mm"] - R["L1", e, "R"]["mm"]) * H for e in ENVS]
say("\nbare-substrate barrier at the fixed geometries (E_TS - E_R): " + ", ".join(
    "%s %.3f" % (b, (R["bare", b, "TS"]["E"] - R["bare", b, "R"]["E"]) * H) for b in ("svp", "svpd", "tzvpd")) + " kcal/mol")
say("\nddE by arrangement (kcal/mol; negative lowers the barrier):")
say("  %-4s %3s %10s %10s %10s %10s %10s %9s" % ("", "n", "first-ord", "L1", "L0", "L2", "L3", "dE_LJ(L1)"))
for k, e in enumerate(ENVS):
    say("  %-4s %3d %+10.3f %+10.3f %+10.3f %+10.3f %+10.3f %+9.3f" % (e, NSPH[e], X["first"][k], X["L1"][k], X["L0"][k], X["L2"][k], X["L3"][k], lj[k]))
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
def compare(test, ref):
    x, y = X[test], X[ref]
    k = sum(a * b for a, b in zip(x, y)) / sum(a * a for a in x)
    res = (sum((b - k * a) ** 2 for a, b in zip(x, y)) / len(x)) ** 0.5
    rho = spearman(x, y)
    signs = all((a > 0) == (b > 0) for a, b in zip(x, y) if abs(b) >= 0.5)
    maxd = max(abs(b - a) for a, b in zip(x, y))
    return k, res, rho, signs, maxd
def verdict(k, res, rho, signs):
    if not (rho >= 0.95 and signs):                 # ranking or sign not kept
        return "FAIL"
    return "PASS" if (res <= 0.5 and 0.8 <= k <= 1.25) else "RANKING ONLY"
say("\ncomparisons (fit ref = k x test through the origin; residual = RMS of ref - k x test; Spearman rho over the 8;"
    " signs compared where |ref| >= 0.5):")
for test, ref, role in (("L1", "L3", "PRIMARY"), ("L2", "L3", "fallback level"), ("L0", "L3", "informative"),
                        ("L1", "L0", "informative: pseudopotential at def2-SVP"), ("first", "L1", "informative: Tier 2, first order against full L1")):
    k, res, rho, signs, maxd = compare(test, ref)
    v = verdict(k, res, rho, signs)
    say("  %-5s vs %-3s k %.3f  residual %.3f  rho %.3f  signs %s  max |diff| %.3f -> %s  (%s)" % (
        test, ref, k, res, rho, "agree" if signs else "DIFFER", maxd, v, role))
k, res, rho, signs, _ = compare("L1", "L3"); v = verdict(k, res, rho, signs)
say("\nS4c verdict (L1 against L3): %s" % v)
say({"PASS": "  -> design and relaxed validation stay at def2-SVP with bare spheres; every final design gets an L3 single-point check.",
     "RANKING ONLY": "  -> def2-SVP ranks designs; magnitudes are judged at L3 (D5 reference 2 then recomputed at L3).",
     "FAIL": "  -> validation of top designs moves to the diffuse level: L2 if its comparison against L3 is PASS, otherwise L3."}[v])
open(os.path.join(D, "STAGE4C_REPORT.txt"), "w").write("\n".join(rep) + "\n")
print("\n".join(rep))
