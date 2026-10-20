"""The active site's charged residues at the reactant and the transition state, extracted exactly as the map's substrate.

The map draws frame 20000's substrate: the reactant is NEB image 0 and the transition state the climbing image
(19_ensemble_barriers/frame_N/reactant_qm.xyz, transition_state_qm.xyz), mapped into the common frame with the frame's
committed transform (phase2.2/aligned/frame_N/transform.txt). The residues come from the same two structures: the NEB's
optimised region at image 0 (neb_MEP.activeRegion_trj.xyz.xz) and at the climbing image
(neb_NEB-CI_converged.activeRegion.xyz.xz), brought down by x89 and checked there (05_qmmm/19_ensemble/active_regions_x89/),
mapped with the same transform. Check: the substrate inside those structures, mapped, must equal the aligned substrate the
map draws (phase2.2/aligned/frame_N/*.xyz) in every one of the 30 band frames.

Residues: the six Stage 4d analyses (cs4d_prepare.py; Szefczyk et al. 2004's charged residues at this site): Arg90, Arg7,
Glu78, Arg116 (chain B 217, 134, 205, 243), Arg63', Lys60' (chain A 63, 60). Drawn: the side chain from CA, every atom,
hydrogens included. Charge centres as Stage 4d: Arg CZ, Lys NZ, Glu CD.
Usage: python3 residues.py <results_repo> <work_dir>   (writes <work_dir>/residues.json and RESIDUES_LOG.txt)
"""
import json, lzma, os, re, sys
import numpy as np

ROOT, HERE = sys.argv[1], sys.argv[2]
P22, E = ROOT + "/phase2.2", ROOT + "/05_qmmm/19_ensemble"
FRAMES = open(P22 + "/cs_stage5/exploratory/dev_A/A_dev.tsv").read().split("# frames (30): ")[1].split("\n")[0].split()
SHOW = "20000"                                   # the frame whose substrate the map draws (the alignment reference)
QM0, NQM = 6207, 24
RES = [("Arg90", 217, "ARG", +1, "CZ", "B"), ("Arg7", 134, "ARG", +1, "CZ", "B"), ("Glu78", 205, "GLU", -1, "CD", "B"),
       ("Arg116", 243, "ARG", +1, "CZ", "B"), ("Arg63′", 63, "ARG", +1, "CZ", "A"), ("Lys60′", 60, "LYS", +1, "NZ", "A")]
BACKBONE = {"N", "H", "C", "O", "HA"}
COV = {"C": 0.76, "H": 0.31, "N": 0.71, "O": 0.66, "S": 1.05}
NAMES = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
out = []
def say(s=""): print(s); out.append(s)
def stop(s): print("STOP: " + s); sys.exit(1)

T = [l.rstrip() for l in open(ROOT + "/03_amber/tleap_build/complex_solvated.prmtop")]
def flag(name):
    i = T.index("%FLAG " + name); fmt = T[i + 1]; j = i + 2; v = []
    w = int(fmt.split("a")[1].split(")")[0]) if "a" in fmt else int(fmt.split("I")[1].split(")")[0]) if "I" in fmt else int(fmt.split("E")[1].split(".")[0])
    while not T[j].startswith("%FLAG"): v += [T[j][k:k + w].strip() for k in range(0, len(T[j]), w)]; j += 1
    return [x for x in v if x]
an, rl = flag("ATOM_NAME"), flag("RESIDUE_LABEL")
rp = [int(x) - 1 for x in flag("RESIDUE_POINTER")] + [len(an)]
qa = [float(x) / 18.2223 for x in flag("CHARGE")]
EL = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S"}
el = [EL.get(int(x), "?") for x in flag("ATOMIC_NUMBER")]
for nm, seq, rn, z, cc, ch in RES:
    idx = range(rp[seq - 1], rp[seq])
    if rl[seq - 1] != rn or abs(sum(qa[i] for i in idx) - z) > 2e-3: stop("%s: residue %d is not %s with charge %+d" % (nm, seq, rn, z))
say("residues 217 134 205 243 (chain B), 63 60 (chain A): ARG ARG GLU ARG ARG LYS, charges +1 +1 -1 +1 +1 +1 (prmtop)")

def structures(path):
    L = lzma.open(path, "rt").read().splitlines(); k = 0; F = []
    while k < len(L) and L[k].strip():
        n = int(L[k].split()[0]); F.append(np.array([[float(x) for x in l.split()[1:4]] for l in L[k + 2:k + 2 + n]])); k += 2 + n
    return F
def xyz(p):
    L = open(p).read().splitlines(); n = int(L[0].split()[0])
    return np.array([[float(x) for x in l.split()[1:4]] for l in L[2:2 + n]])

per = {}        # frame -> state -> aligned coordinates of every system atom in the active region (dict index -> xyz)
worst = 0.0
for f in FRAMES:
    act = [int(x) for x in re.search(r"ActiveAtoms \{([^}]*)\}", open(E + "/frame_%s/neb.inp" % f).read()).group(1).split()]
    M = np.loadtxt(P22 + "/aligned/frame_%s/transform.txt" % f, comments="#"); Rm, t = M[:3], M[3]
    S = {"R": structures(E + "/frame_%s/neb_MEP.activeRegion_trj.xyz.xz" % f)[0],
         "TS": structures(E + "/frame_%s/neb_NEB-CI_converged.activeRegion.xyz.xz" % f)[0]}
    q = act.index(QM0); per[f] = {}
    for st, X in S.items():
        A = X @ Rm + t
        ref = xyz(P22 + "/aligned/frame_%s/%s" % (f, "reactant_qm.xyz" if st == "R" else "transition_state_qm.xyz"))
        d = float(np.abs(A[q:q + NQM] - ref).max()); worst = max(worst, d)
        if d > 2e-5: stop("frame %s %s: the substrate in the NEB structure, mapped, differs from the aligned substrate by %.1e A" % (f, st, d))
        per[f][st] = {i: A[k] for k, i in enumerate(act)}
say("all 30 frames, R and TS: the substrate inside the NEB structures, mapped with the frame's transform, equals the aligned")
say("  substrate the map draws (phase2.2/aligned/frame_N/reactant_qm.xyz, transition_state_qm.xyz) within %.1e A" % worst)

def bonds(names, els, X):
    b = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            if np.linalg.norm(X[i] - X[j]) < COV[els[i]] + COV[els[j]] + 0.40: b.append([i, j])
    return b

rows = []
say("\nresidue  charge centre: 30-frame spread (rms about the mean, A) R / TS;  R->TS movement at frame %s: centre / largest atom;"
    " median over 30 frames" % SHOW)
for nm, seq, rn, z, cc, ch in RES:
    idx = [i for i in range(rp[seq - 1], rp[seq]) if an[i] not in BACKBONE]
    ci = [i for i in idx if an[i] == cc][0]
    if any(i not in per[SHOW]["R"] for i in idx): stop("%s: not every side-chain atom is in the active region" % nm)
    draw = {st: np.array([per[SHOW][st][i] for i in idx]) for st in ("R", "TS")}
    spread = {}
    for st in ("R", "TS"):
        P = np.array([per[f][st][ci] for f in FRAMES]); spread[st] = float(np.sqrt((np.linalg.norm(P - P.mean(0), axis=1) ** 2).mean()))
    mv = [float(np.linalg.norm(per[f]["TS"][ci] - per[f]["R"][ci])) for f in FRAMES]
    big = float(np.linalg.norm(draw["TS"] - draw["R"], axis=1).max())
    sub = {st: np.array([per[SHOW][st][i] for i in range(QM0, QM0 + NQM)]) for st in ("R", "TS")}
    near = {st: (lambda d: (NAMES[int(d.argmin())], float(d.min())))(np.linalg.norm(sub[st] - per[SHOW][st][ci], axis=1)) for st in ("R", "TS")}
    say("%-8s %-2s  spread %.2f / %.2f   frame %s: %.2f / %.2f   median %.2f   nearest substrate atom R %s %.2f, TS %s %.2f"
        % (nm, cc, spread["R"], spread["TS"], SHOW, mv[FRAMES.index(SHOW)], big, float(np.median(mv)), near["R"][0], near["R"][1], near["TS"][0], near["TS"][1]))
    names = [an[i] for i in idx]; els = [el[i] for i in idx]
    rows.append({"lab": nm, "res": rn, "seq": seq, "chain": ch, "q": z, "centre": cc, "names": names, "el": els,
                 "bonds": bonds(names, els, draw["R"]),
                 "R": np.round(draw["R"], 3).tolist(), "TS": np.round(draw["TS"], 3).tolist(),
                 "spreadR": round(spread["R"], 2), "spreadTS": round(spread["TS"], 2), "moveCentre": round(mv[FRAMES.index(SHOW)], 2),
                 "moveMax": round(big, 2), "moveMedian": round(float(np.median(mv)), 2),
                 "nearR": [near["R"][0], round(near["R"][1], 2)], "nearTS": [near["TS"][0], round(near["TS"][1], 2)]})
    if bonds(names, els, draw["R"]) != bonds(names, els, draw["TS"]): stop("%s: bonds differ between R and TS" % nm)
json.dump({"frame": SHOW, "residues": rows, "frames": FRAMES}, open(os.path.join(HERE, "residues.json"), "w"), indent=1)
open(os.path.join(HERE, "RESIDUES_LOG.txt"), "w").write("\n".join(out) + "\n")
