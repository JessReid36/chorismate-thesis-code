#!/usr/bin/env python3
"""
cs4c_amend1_direction.py - EXPLORATORY (not pre-registered; written 10 October 2026, after
STAGE4C_AMENDMENT1_REPORT.txt was seen). In which direction does Amendment 1's zero-charge correction move each
arrangement's ddE? It answers one question from the audit for anything that could favour a lower barrier: does the
correction systematically make arrangements look more barrier-lowering? Changes no verdict. The energy reading and
the correction are those of cs4c_amend1_explore.py (copied, unchanged).
Reads each run's job.out, or job.out.xz (committed outputs). Python 3.6, standard library.
USAGE  python3 cs4c_amend1_direction.py <cs_stage4c folder>
       -> writes <folder>/exploratory/STAGE4C_AMENDMENT1_DIRECTION_REPORT.txt
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
out = []; say = out.append
say("Stage 4c Amendment 1 - direction of the zero-charge correction (cs4c_amend1_direction.py; EXPLORATORY; verdicts unchanged)")
say("corrected - uncorrected ddE, kcal/mol; < 0 means the correction makes the arrangement look more barrier-lowering.")
say("L1 (bare spheres, def2-SVP, the production level) is not corrected, so production values are untouched.")
for L in ("L0", "L2", "L3"):
    d = {e: X[L, e] - RAW[L, e] for e in ENVS}; cat = [e for e in ENVS if X["L1", e] < 0]
    say("  %s  %s" % (L, " ".join("%s %+.3f" % (e, d[e]) for e in ENVS)))
    say("      catalytic six (L1 < 0): %d made to look more lowering, %d less; sum %+.3f" % (
        sum(d[e] < 0 for e in cat), sum(d[e] > 0 for e in cat), sum(d[e] for e in cat)))
os.makedirs(os.path.join(D, "exploratory"), exist_ok=True)
open(os.path.join(D, "exploratory", "STAGE4C_AMENDMENT1_DIRECTION_REPORT.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
