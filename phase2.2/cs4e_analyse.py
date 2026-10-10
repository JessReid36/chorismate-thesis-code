#!/usr/bin/env python3
"""
cs4e_analyse.py - Stage 4e analysis (method dependence of the sign beside the carboxylates; STAGE4E_CRITERIA.txt).
Reads runs.tsv, probes.tsv, each run's job.out (or job.out.xz once committed) and, for check R0, Stage 4d's outputs in
the sibling folder cs_stage4d. Python 3.6, standard library. Writes STAGE4E_REPORT.txt.
USAGE  python3 cs4e_analyse.py [cs_stage4e folder]
"""
import lzma, os, re, sys
D = sys.argv[1] if len(sys.argv) > 1 else "/home/18660916/system_development/phase2.2/cs_stage4e"
D4 = os.path.join(os.path.dirname(os.path.abspath(D)), "cs_stage4d")
H = 627.5094740631
LEV = ["L1", "L2", "L3", "L4", "L5"]
NAME = {"L1": "B3LYP-D3BJ/def2-SVP", "L2": "B3LYP/6-31G(d)", "L3": "HF/6-31G(d)", "L4": "MP2/6-31G(d)", "L5": "wB97X-D3/def2-SVP"}
RESIDUES = ["R90", "R7", "E78", "R116", "R63p", "K60p"]
RNAME = {"R90": "Arg90", "R7": "Arg7", "E78": "Glu78", "R116": "Arg116", "R63p": "Arg63'", "K60p": "Lys60'"}
SITES = [37, 52, 87, 209, 211]
C = ["ALONE_R7", "ALONE_R63p", "a52", "a87", "a209", "a211"]
K = ["ALONE_R90", "a37"]
ALL11 = ["ALONE_" + r for r in RESIDUES] + ["a%d" % s for s in SITES]
SZ_SCF = {"R90": -11.67, "R7": -6.66, "E78": -6.12, "R116": -3.70, "R63p": -0.79, "K60p": 1.32}   # Table 1, dSCF
GAP_FLAG = 1.0
def text(d, run):
    p = os.path.join(d, run, "job.out")
    if os.path.exists(p): return open(p, errors="replace").read()
    if os.path.exists(p + ".xz"): return lzma.open(p + ".xz", "rt", errors="replace").read()
    return ""
def parse(d, run):
    t = text(d, run)
    if "ORCA TERMINATED NORMALLY" not in t: return None
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    if not e: return None
    nel = re.findall(r"Number of Electrons\s+NEL\s+\.+\s+(\d+)", t)
    nbf = re.findall(r"Number of basis functions\s+\.+\s+(\d+)", t)
    npc = re.findall(r"Reading point charge file\s+\.+\s+ok \((\d+) point charges\)", t)
    i = t.rfind("ORBITAL ENERGIES"); orbs = []
    for l in t[i:i + 80000].splitlines()[4:]:
        m = re.match(r"\s*\d+\s+(\d\.\d+)\s+-?\d+\.\d+\s+(-?\d+\.\d+)\s*$", l)
        if m: orbs.append((float(m.group(1)), float(m.group(2))))
        elif orbs: break
    gap = (min(o[1] for o in orbs if o[0] < 0.5) - max(o[1] for o in orbs if o[0] > 0.5)) if orbs else None
    return {"E": float(e[-1]), "nel": int(nel[-1]) if nel else None, "nbf": int(nbf[0]) if nbf else None,
            "npc": int(npc[-1]) if npc else 0, "gap": gap}
RUNS = [l.rstrip("\n").split("\t") for l in open(os.path.join(D, "runs.tsv")) if not l.startswith("#") and l.strip()]
PROBE = {}
for l in open(os.path.join(D, "probes.tsv")):
    if not l.startswith("#") and l.strip():
        f = l.rstrip("\n").split("\t"); PROBE[int(f[0])] = (f[3], float(f[4]))
rep = []; say = rep.append
def finish():
    open(os.path.join(D, "STAGE4E_REPORT.txt"), "w").write("\n".join(rep) + "\n"); print("\n".join(rep)); sys.exit(0)
say("Stage 4e report (cs4e_analyse.py) - criteria from STAGE4E_CRITERIA.txt")
P = {r[0]: parse(D, r[0]) for r in RUNS}
miss = [k for k, v in P.items() if v is None]
if miss:
    say("%d run(s) missing or not terminated normally: %s -> not assessable" % (len(miss), " ".join(miss[:12]))); finish()
nbf = {L: sorted({P[r[0]]["nbf"] for r in RUNS if r[1] == L}) for L in LEV}
ok0 = all(P[r[0]]["nel"] == 118 and P[r[0]]["npc"] == int(r[4]) for r in RUNS) and all(len(v) == 1 for v in nbf.values())
say("check 0 (%d runs terminated; 118 electrons; one basis-function count per level %s; point charges read = runs.tsv): %s"
    % (len(RUNS), " ".join("%s %s" % (L, nbf[L][0] if len(nbf[L]) == 1 else nbf[L]) for L in LEV), "PASS" if ok0 else "FAIL"))
r0 = []
for e in ["bare"] + ["ALONE_" + r for r in RESIDUES]:
    for g in ("R", "TS"):
        ref = parse(D4, "%s_%s" % (e, g))
        r0.append(abs(P["L1_%s_%s" % (e, g)]["E"] - ref["E"]) if ref else float("inf"))
okR = max(r0) <= 1e-6
say("R0 L1 bare and residue-alone energies against Stage 4d's (same inputs): max |diff| %.2e Eh (tolerance 1e-6) -> %s"
    % (max(r0), "PASS" if okR else "FAIL"))
if not (ok0 and okR):
    say("\nmachinery not verified -> not assessable"); finish()
def dd(L, e):
    return ((P["%s_%s_TS" % (L, e)]["E"] - P["%s_bare_TS" % L]["E"]) - (P["%s_%s_R" % (L, e)]["E"] - P["%s_bare_R" % L]["E"])) * H
V = {}
for L in LEV:
    for r in RESIDUES: V[L, "ALONE_" + r] = dd(L, "ALONE_" + r)
    for s in SITES:
        p, m = dd(L, "P%dp" % s), dd(L, "P%dm" % s)
        V[L, "a%d" % s] = (p - m) / 2; V[L, "pol%d" % s] = (p + m) / 2
def label(x):
    if x.startswith("ALONE_"): return RNAME[x[6:]] + " alone"
    s = int(x[1:]); return "a(%d) [%s]" % (s, PROBE[s][0].split(",")[0])
say("\nbarrier change by level (kcal/mol, QM part; negative lowers the barrier; a = first-order part of a unit charge):")
say("  %-26s %s" % ("", " ".join("%10s" % L for L in LEV)))
for x in ALL11:
    tag = " C" if x in C else (" K" if x in K else "  ")
    say("  %-24s%s %s" % (label(x), tag, " ".join("%+10.3f" % V[L, x] for L in LEV)))
say("  polarisation part [ddE(+1) + ddE(-1)]/2 at the probe sites:")
for s in SITES:
    say("  %-26s %s" % ("site %d" % s, " ".join("%+10.3f" % V[L, "pol%d" % s] for L in LEV)))
say("  levels: " + "; ".join("%s %s" % (L, NAME[L]) for L in LEV))
ctrl = {L: all(V[L, x] < 0 for x in K) for L in LEV}
say("\ncontrols (Arg90 alone and a(37) both negative): " + ", ".join("%s %s" % (L, "yes" if ctrl[L] else "NO") for L in LEV))
if not ctrl["L4"]:
    say("the reference level (L4, MP2) fails the controls -> not assessable"); finish()
LOWGAP = {r[0] for r in RUNS if P[r[0]]["gap"] is not None and P[r[0]]["gap"] < GAP_FLAG}
def runs_of(L, x):
    envs = [x] if x.startswith("ALONE_") else ["P%sp" % x[1:], "P%sm" % x[1:]]
    return {"%s_%s_%s" % (L, e, g) for e in envs + ["bare"] for g in ("R", "TS")}
UNRES = {(L, x) for L in LEV for x in ALL11 if runs_of(L, x) & LOWGAP}
def signs_agree(L):
    return [x for x in C + K if abs(V["L4", x]) >= 0.5 and (L, x) not in UNRES and ("L4", x) not in UNRES
            and (V[L, x] > 0) != (V["L4", x] > 0)]
def fit(L):
    x = [V[L, v] for v in ALL11]; y = [V["L4", v] for v in ALL11]
    k = sum(a * b for a, b in zip(x, y)) / sum(a * a for a in x)
    return k, (sum((b - k * a) ** 2 for a, b in zip(x, y)) / len(x)) ** 0.5
say("\nagainst the reference L4 (signs compared where |L4| >= 0.5; fit L4 = k x level over the eleven values):")
for L in ("L1", "L2", "L3", "L5"):
    k, res = fit(L); bad = signs_agree(L)
    say("  %s %-20s k %.3f  residual %.3f  C and K signs %s" % (L, NAME[L], k, res, "all as L4" if not bad else
        "DIFFER for " + ", ".join(label(b) for b in bad)))
bad1 = signs_agree("L1")
if not bad1:
    verdict = "SIGN ROBUST"
    cons = ("the sign beside the carboxylates does not follow the method at this geometry; Stage 5 builds the A matrix at"
            " L1, and the disagreement with Szefczyk et al. is recorded as structural, not resolved")
else:
    verdict = "METHOD-DEPENDENT"
    elig = [L for L in ("L2", "L5") if ctrl[L] and not signs_agree(L)]
    if elig:
        best = min(elig, key=lambda L: fit(L)[1])
        cons = ("L1 is not used to build the A matrix; Stage 5 builds it at %s (%s), the eligible DFT level with the smallest"
                " residual against L4; every Stage 4c and 4d value beside a carboxylate is marked as level-dependent" % (best, NAME[best]))
    else:
        cons = ("L1 is not used to build the A matrix; neither L2 nor L5 keeps L4's signs, so Stage 5 builds it from MP2"
                " densities; every Stage 4c and 4d value beside a carboxylate is marked as level-dependent")
say("\ninformative: HF/6-31G(d) here against Szefczyk et al.'s HF-level values (Table 1, dSCF) for the residues alone:")
for r in RESIDUES:
    v = V["L3", "ALONE_" + r]
    say("  %-7s here %+8.3f   Szefczyk dSCF %+7.2f   %s" % (RNAME[r], v, SZ_SCF[r], "same sign" if (v > 0) == (SZ_SCF[r] > 0) else "opposite sign"))
say("informative: L1's first-order parts against A_v2's a for frame 41786 (aligned frame, in vacuo densities):")
for s in SITES:
    say("  site %3d  L1 a %+7.3f   A_v2 %+7.3f" % (s, V["L1", "a%d" % s], PROBE[s][1]))
flag = sorted(LOWGAP)
say("electronic screen: smallest HOMO-LUMO gap per level (eV): " + ", ".join(
    "%s %.2f" % (L, min(P[r[0]]["gap"] for r in RUNS if r[1] == L)) for L in LEV))
say("  runs with a gap below %.1f eV (values unresolved): %s" % (GAP_FLAG, " ".join(flag) if flag else "none"))
say("  C and K values left out of the sign comparison as unresolved: %s" % (", ".join(
    "%s %s" % (L, label(x)) for (L, x) in sorted(UNRES) if x in C + K and L in ("L1", "L2", "L4", "L5")) or "none"))
say("\nS4e verdict (L1 against L4): %s" % verdict)
say("consequence (STAGE4E_CRITERIA.txt): " + cons)
finish()
