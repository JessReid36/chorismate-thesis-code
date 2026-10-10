#!/usr/bin/env python3
"""
cs4d_analyse.py - Stage 4d analysis (positive and negative controls; STAGE4D_CRITERIA.txt). Reads runs.tsv, lj.tsv,
residues.tsv and each run's job.out (or job.out.xz once committed). Python 3.6, standard library.
Writes STAGE4D_REPORT.txt. Every threshold and reference value below is the one fixed in STAGE4D_CRITERIA.txt.
USAGE  python3 cs4d_analyse.py [cs_stage4d folder]
"""
import lzma, os, re, sys
D = sys.argv[1] if len(sys.argv) > 1 else "/home/18660916/system_development/phase2.2/cs_stage4d"
H = 627.5094740631
E_VAC = {"R": -836.142634941140, "TS": -836.116074935829}    # 05_qmmm/20_invacuo/sp_41786_{R,TS}.out
E_QM_R_P1 = -836.914536233376                                 # 05_qmmm/19_ensemble/frame_41786/reactant_opt.out, QM part
STAB_41786 = -0.806655                                        # ensemble_barriers.tsv, frame 41786
BAND = (-13.3, 3.3)                                           # stab_TS mean +- 2 sd, 43 frames (results 5fa30a2)
GAP_FLAG = 1.0                                                # eV
# Szefczyk et al. 2004, Table 1 (kcal/mol): dV, dV-HF, dEL.MTP(1), dEL(1), d(1), dSCF, dMP2 (Arg7/Arg90: first contact)
SZ = {"R90": (-7.55, -3.89, -5.23, -11.00, -7.18, -11.67, -9.06), "R7": (-8.18, -4.11, -6.75, -8.35, -2.81, -6.66, -5.90),
      "E78": (-3.15, -2.80, -4.67, -1.33, -11.17, -6.12, -3.57), "R116": (-7.32, -0.48, -2.58, -2.47, -2.76, -3.70, -2.45),
      "R63p": (0.76, 0.32, -0.09, -0.37, -1.32, -0.79, -1.40), "K60p": (1.08, 0.38, 1.30, 1.30, 1.30, 1.32, 1.39)}
JUDGED = {"R90": -1, "R7": -1, "R116": -1, "K60p": 1}
NAME = {"R90": "Arg90", "R7": "Arg7", "E78": "Glu78", "R116": "Arg116", "R63p": "Arg63'", "K60p": "Lys60'"}
def text(run):
    p = os.path.join(D, run, "job.out")
    if os.path.exists(p): return open(p, errors="replace").read()
    if os.path.exists(p + ".xz"): return lzma.open(p + ".xz", "rt", errors="replace").read()
    return ""
def parse(run):
    t = text(run)
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
    i = t.rfind("HIRSHFELD ANALYSIS"); hq = []
    for l in t[i:i + 4000].splitlines():
        m = re.match(r"\s*(\d+)\s+[A-Za-z]+\s+(-?\d+\.\d+)\s+-?\d+\.\d+\s*$", l)
        if m: hq.append(float(m.group(2)))
    return {"E": float(e[-1]), "nel": int(nel[-1]) if nel else None, "nbf": int(nbf[0]) if nbf else None,
            "npc": int(npc[-1]) if npc else 0, "gap": gap, "hq": hq[:24]}
RUNS = [l.rstrip("\n").split("\t") for l in open(os.path.join(D, "runs.tsv")) if not l.startswith("#") and l.strip()]
LJ = {}
for l in open(os.path.join(D, "lj.tsv")):
    if not l.startswith("#") and l.strip():
        f = l.split("\t"); LJ[f[0]] = float(f[4])
ENVOF = {}; LJOF = {}
for r in RUNS:
    env = r[0].rsplit("_", 1)[0]; ENVOF[env] = r[1]; LJOF[env] = r[4]
rep = []; say = rep.append
def finish():
    open(os.path.join(D, "STAGE4D_REPORT.txt"), "w").write("\n".join(rep) + "\n"); print("\n".join(rep)); sys.exit(0)
say("Stage 4d report (cs4d_analyse.py) - criteria from STAGE4D_CRITERIA.txt")
P = {r[0]: parse(r[0]) for r in RUNS}
miss = [k for k, v in P.items() if v is None]
if miss:
    say("%d run(s) missing or not terminated normally: %s -> not assessable" % (len(miss), " ".join(miss[:12]))); finish()
ok0 = all(P[r[0]]["nel"] == 118 and P[r[0]]["nbf"] == 264 and P[r[0]]["npc"] == int(r[3]) for r in RUNS)
say("check 0 (%d runs terminated; 118 electrons; 264 basis functions; point charges read = runs.tsv): %s"
    % (len(RUNS), "PASS" if ok0 else "FAIL"))
for r in RUNS:
    p = P[r[0]]
    if not (p["nel"] == 118 and p["nbf"] == 264 and p["npc"] == int(r[3])):
        say("  %s: NEL %s, basis functions %s, point charges %s (expected %s)" % (r[0], p["nel"], p["nbf"], p["npc"], r[3]))
d0 = {g: P["bare_" + g]["E"] - E_VAC[g] for g in ("R", "TS")}
okA = all(abs(v) <= 1e-6 for v in d0.values())
say("R0a bare R, TS against Phase 1's in vacuo single points: %+.2e, %+.2e Eh (tolerance 1e-6) -> %s"
    % (d0["R"], d0["TS"], "PASS" if okA else "FAIL"))
dB = P["ENZ_R"]["E"] - E_QM_R_P1
okB = abs(dB) <= 5e-4
say("R0b ENZ R against Phase 1's QM part of the QM/MM reactant: %+.2e Eh (%+.3f kcal/mol; tolerance 5e-4 Eh) -> %s"
    % (dB, dB * H, "PASS" if okB else "FAIL"))
dN = {g: P["NC0_" + g]["E"] - P["bare_" + g]["E"] for g in ("R", "TS")}
okN = all(abs(v) <= 1e-6 for v in dN.values())
say("NC0 zero charges against bare: %+.2e, %+.2e Eh (tolerance 1e-6) -> %s" % (dN["R"], dN["TS"], "PASS" if okN else "FAIL"))
if not (ok0 and okA and okB and okN):
    say("\nmachinery not verified -> stage not assessable; the cause is found before anything else uses frame 41786's yardstick")
    finish()
def qm(env):
    return ((P[env + "_TS"]["E"] - P["bare_TS"]["E"]) - (P[env + "_R"]["E"] - P["bare_R"]["E"])) * H
def lj(env): return LJ[LJOF[env]]
def tot(env): return qm(env) + lj(env)
envs = [r[0].rsplit("_", 1)[0] for r in RUNS if r[0].endswith("_R") and r[0] != "bare_R"]
say("\nddE by environment (kcal/mol; negative lowers the barrier; total = QM + LJ):")
say("  %-11s %10s %10s %10s" % ("", "QM", "LJ", "total"))
for e in envs:
    say("  %-11s %+10.3f %+10.3f %+10.3f" % (e, qm(e), lj(e), tot(e)))
own_qm = ((P["ENZTS_TS"]["E"] - P["bare_TS"]["E"]) - (P["ENZ_R"]["E"] - P["bare_R"]["E"])) * H
# own environment: R in env_R, TS in env_TS; LJ(TS in env_TS) - LJ(R in env_R) from lj.tsv's absolute values
LJABS = {}
for l in open(os.path.join(D, "lj.tsv")):
    if not l.startswith("#") and l.strip():
        f = l.split("\t"); LJABS[f[0]] = (float(f[2]), float(f[3]))
own_lj = LJABS["envTS_all"][1] - LJABS["envR_all"][0]
say("\nPC1 the whole enzyme, measured with the design yardstick:")
v1 = tot("ENZ"); pc1 = BAND[0] <= v1 <= BAND[1]
say("  ENZ (frozen at the reactant's environment) total %+.3f kcal/mol; band [%.1f, %.1f] -> %s" % (v1, BAND[0], BAND[1], "PASS" if pc1 else "FAIL"))
say("  informative: ENZTS (frozen at the TS's environment) %+.3f; own environment (R in env_R, TS in env_TS) QM %+.3f + LJ %+.3f = %+.3f"
    % (tot("ENZTS"), own_qm, own_lj, own_qm + own_lj))
say("  informative: frame 41786's stab_TS (Phase 1, includes the protein's own MM energy change) %+.3f; the remainder (mainly that MM change) %+.3f;"
    " Claeyssens et al. 2005 average 4.2 kcal/mol of TS stabilisation" % (STAB_41786, STAB_41786 - (own_qm + own_lj)))
odd = (qm("ENZ") - qm("REV")) / 2; even = (qm("ENZ") + qm("REV")) / 2
say("NC1 (informative) reversed charges: QM part odd (first order) %+.3f, even (polarisation) %+.3f kcal/mol" % (odd, even))
say("\nPC2 single residues against Szefczyk et al. 2004 (MP2 and the range over their seven levels, kcal/mol):")
labs = [l.split("\t")[0] for l in open(os.path.join(D, "residues.tsv")) if not l.startswith("#") and l.strip()]
pc2 = True
for lab in labs:
    a = tot("ALONE_" + lab); insitu = qm("ENZ") - qm("KO_" + lab)
    sz = SZ[lab]
    if lab in JUDGED:
        m = (a * JUDGED[lab] > 0) and abs(a) >= 0.05; pc2 = pc2 and m; tag = "sign %s" % ("matches" if m else "DIFFERS")
    else:
        tag = "not judged"
    say("  %-7s alone %+8.3f (QM %+8.3f)  in situ (ENZ - KO) %+8.3f   Szefczyk MP2 %+6.2f, range %+6.2f to %+6.2f   %s"
        % (NAME[lab], a, qm("ALONE_" + lab), insitu, sz[6], min(sz), max(sz), tag))
say("  sum of the six alone %+.3f against ENZ %+.3f (total)" % (sum(tot("ALONE_" + l) for l in labs), tot("ENZ")))
say("  PC2 -> %s" % ("PASS" if pc2 else "FAIL"))
say("\nPC3 a sphere standing in for a residue (QM part; totals in the table above):")
s90, a90 = qm("SPH_R90CZ"), qm("ALONE_R90"); s7, a7 = qm("SPH_R7CZ"), qm("ALONE_R7")
pc3 = (s90 * a90 > 0) and (s7 * a7 > 0)
say("  +1 at Arg90 CZ %+.3f against Arg90 alone %+.3f: r90 = %.3f" % (s90, a90, s90 / a90 if a90 else float("nan")))
say("  +1 at Arg7 CZ  %+.3f against Arg7 alone  %+.3f: r7  = %.3f" % (s7, a7, s7 / a7 if a7 else float("nan")))
say("  +1 at grid site 37 (contact distance) %+.3f against +1 at Arg90 CZ: g = %.3f" % (qm("SPH_G37"), qm("SPH_G37") / s90 if s90 else float("nan")))
say("  Arg90 replaced by a sphere inside the enzyme: SUB_R90 - ENZ = %+.3f" % (qm("SUB_R90") - qm("ENZ")))
say("  PC3 (signs kept for Arg90 and Arg7) -> %s" % ("PASS" if pc3 else "FAIL"))
say("\nelectronic screen, HOMO-LUMO gap (eV) and largest shift of a substrate Hirshfeld charge against bare (e):")
flag = []
for r in RUNS:
    p = P[r[0]]; g = r[0].rsplit("_", 1)[1]; b = P["bare_" + g]["hq"]
    sh = max(abs(x - y) for x, y in zip(p["hq"], b)) if len(p["hq"]) == 24 and len(b) == 24 else float("nan")
    if p["gap"] is not None and p["gap"] < GAP_FLAG: flag.append(r[0])
    say("  %-14s gap %5.2f  max dq %.3f" % (r[0], p["gap"] if p["gap"] is not None else float("nan"), sh))
say("  runs with a gap below %.1f eV (values unresolved): %s" % (GAP_FLAG, " ".join(flag) if flag else "none"))
say("\nverdicts: PC1 %s, PC2 %s, PC3 %s" % ("PASS" if pc1 else "FAIL", "PASS" if pc2 else "FAIL", "PASS" if pc3 else "FAIL"))
say("consequences (STAGE4D_CRITERIA.txt): " + "; ".join([
    "PC1 PASS: the enzyme's frozen value (ENZ total) is the screening reference for designs; the catalytic verdict stays with the relaxed NEB-CI barrier against D5's reference 2"
    if pc1 else "PC1 FAIL: frozen values rank designs only and are never compared with an enzyme number",
    "PC2 PASS" if pc2 else "PC2 FAIL: the design objective is re-examined before the optimiser stage",
    "PC3 PASS: r90, r7, g and SUB_R90 - ENZ are quoted beside every design magnitude" if pc3
    else "PC3 FAIL: the sphere representation is re-examined before Stage 5"]))
finish()
