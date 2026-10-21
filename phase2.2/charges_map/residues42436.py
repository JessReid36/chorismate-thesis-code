"""residues42436.py - the active site's six charged residues at frame 42436's reactant and transition state.

Usage: python3 residues42436.py <results_repo> <work_dir>   (writes <work_dir>/residues.json and RESIDUES_LOG.txt)

Why this is not site_effect_map/residues.py. That script takes the residues from the NEB active-region structures of
the 30 BAND frames (05_qmmm/19_ensemble/frame_N/), and frame 42436 is not one of them - it is an SP-package frame.
Here the residues come from the SAME committed archives cs6_prepare.py took the substrate from
(05_qmmm/24_sp_barriers/f42436_r8_a1/hpc_text_outputs_2_{ts,r_irc}.tar.xz), whose final.xyz holds all 55680 system
atoms. Residues and substrate therefore come out of one structure, and no transform is involved because frame
42436's is the identity (STAGE6_CRITERIA.txt, cs6_prepare.py).

Residues, charges and charge centres exactly as cs4d_prepare.py and site_effect_map/residues.py define them, so the
two pages name the same six things the same way: Arg90, Arg7, Glu78, Arg116 (chain B 217, 134, 205, 243), Arg63',
Lys60' (chain A 63, 60); side chain from CA, every atom, hydrogens included; charge centres Arg CZ, Lys NZ, Glu CD.

Every check is a stop:
  (a) each archive's sha256 is the value STAGE6_CRITERIA.txt records (the same constants cs6_prepare.py pins);
  (b) each residue is the expected label in the prmtop and its atomic charges sum to its formal charge;
  (c) the 24 substrate atoms inside each final.xyz are byte-equal in value to the committed frame_42436/*.xyz that
      the page draws - this is what ties the residues and the substrate to one structure;
  (d) every side-chain atom is present, and the side chain's bonds are the same at R and at TS.
"""
import hashlib, json, lzma, os, sys, tarfile
import numpy as np

ROOT, WORK = sys.argv[1], sys.argv[2]
RUN = "05_qmmm/24_sp_barriers/f42436_r8_a1"
ARCHIVES = {"ts": "6f01f6a3efe5c0d0b94c0de9de805a23cc126e7fde7e9f50de25bdbf465c9022",
            "r_irc": "3f9ae488565c11d6b9aa14f98886a5e37e6dfe0b4e314038778a35d359bf354b"}
STATE_OF = {"r_irc": "R", "ts": "TS"}
QM0, NQM = 6207, 24
NAMES = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
RES = [("Arg90", 217, "ARG", +1, "CZ", "B"), ("Arg7", 134, "ARG", +1, "CZ", "B"), ("Glu78", 205, "GLU", -1, "CD", "B"),
       ("Arg116", 243, "ARG", +1, "CZ", "B"), ("Arg63′", 63, "ARG", +1, "CZ", "A"), ("Lys60′", 60, "LYS", +1, "NZ", "A")]
BACKBONE = {"N", "H", "C", "O", "HA"}
COV = {"C": 0.76, "H": 0.31, "N": 0.71, "O": 0.66, "S": 1.05}
# Literature catalytic importance, as site_effect_map/prep.py records it from the papers read in full (kcal/mol):
#   Szefczyk et al. 2004, J. Am. Chem. Soc. 126, 16148, Table 1, MP2 differential TS stabilisation (negative lowers
#   the barrier). Claeyssens et al. 2011, Org. Biomol. Chem. 9, 1578, Table 3, electrostatic TS stabilisation, mean
#   (sd) over 16 B3LYP/CHARMM27 pathways (positive stabilises the TS; Lys60' not reported).
LIT = {"Arg90": (-9.06, 6.0, 0.8), "Arg7": (-5.90, 3.1, 0.5), "Glu78": (-3.57, 2.1, 0.3), "Arg116": (-2.45, 1.0, 0.9),
       "Arg63′": (-1.40, -0.6, 0.5), "Lys60′": (1.39, None, None)}

log = []


def say(s=""):
    print(s)
    log.append(s)


def stop(s):
    sys.exit("STOP: " + s)


def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()


# ---- (a) the archives, and the full-system geometry out of each
X = {}
for key, want in ARCHIVES.items():
    p = os.path.join(ROOT, RUN, "hpc_text_outputs_2_%s.tar.xz" % key)
    got = sha_bytes(open(p, "rb").read())
    if got != want:
        stop("%s sha256 %s, STAGE6_CRITERIA.txt records %s" % (os.path.basename(p), got, want))
    with tarfile.open(p) as t:
        m = [n for n in t.getnames() if n.endswith("/final.xyz")]
        if len(m) != 1:
            stop("%s holds %d final.xyz" % (p, len(m)))
        L = t.extractfile(m[0]).read().decode().splitlines()
    n = int(L[0].split()[0])
    X[STATE_OF[key]] = np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])
    say("%-6s %s  %d atoms  archive sha256 ok" % (STATE_OF[key], m[0], n))
if X["R"].shape != X["TS"].shape:
    stop("R and TS have different atom counts")

# ---- (c) the substrate inside these structures is the one the page draws
for st, fn in (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz")):
    ref = os.path.join(ROOT, "phase2.2", "cs_stage6", "frame_42436", fn)
    Lr = open(ref).read().splitlines()
    A = np.array([[float(v) for v in l.split()[1:4]] for l in Lr[2:2 + NQM]])
    d = float(np.abs(X[st][QM0:QM0 + NQM] - A).max())
    if d > 1e-9:
        stop("%s: the substrate inside final.xyz differs from %s by %.1e A" % (st, fn, d))
    say("%-6s the 24 substrate atoms inside final.xyz equal %s exactly (max |d| %.1e A)" % (st, fn, d))

# ---- (b) the prmtop, and each residue's identity and charge
T = [l.rstrip() for l in open(os.path.join(ROOT, "03_amber/tleap_build/complex_solvated.prmtop"))]


def flag(name):
    i = T.index("%FLAG " + name)
    fmt = T[i + 1]
    j = i + 2
    v = []
    w = (int(fmt.split("a")[1].split(")")[0]) if "a" in fmt else
         int(fmt.split("I")[1].split(")")[0]) if "I" in fmt else
         int(fmt.split("E")[1].split(".")[0]))
    while not T[j].startswith("%FLAG"):
        v += [T[j][k:k + w].strip() for k in range(0, len(T[j]), w)]
        j += 1
    return [x for x in v if x]


an, rl = flag("ATOM_NAME"), flag("RESIDUE_LABEL")
rp = [int(x) - 1 for x in flag("RESIDUE_POINTER")] + [len(an)]
qa = [float(x) / 18.2223 for x in flag("CHARGE")]
EL = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S"}
el = [EL.get(int(x), "?") for x in flag("ATOMIC_NUMBER")]
if len(an) != X["R"].shape[0]:
    stop("the prmtop has %d atoms, final.xyz has %d" % (len(an), X["R"].shape[0]))
for nm, seq, rn, z, cc, ch in RES:
    idx = range(rp[seq - 1], rp[seq])
    if rl[seq - 1] != rn or abs(sum(qa[i] for i in idx) - z) > 2e-3:
        stop("%s: residue %d is not %s with charge %+d" % (nm, seq, rn, z))
say("\nresidues 217 134 205 243 (chain B), 63 60 (chain A): ARG ARG GLU ARG ARG LYS, charges +1 +1 -1 +1 +1 +1 (prmtop)")


def bonds(els, P):
    return [[i, j] for i in range(len(els)) for j in range(i + 1, len(els))
            if np.linalg.norm(P[i] - P[j]) < COV[els[i]] + COV[els[j]] + 0.40]


rows = []
say("\nresidue  charge centre, frame 42436: nearest substrate atom at R and at TS; R->TS movement of the centre and")
say("  of the largest-moving side-chain atom")
for nm, seq, rn, z, cc, ch in RES:
    idx = [i for i in range(rp[seq - 1], rp[seq]) if an[i] not in BACKBONE]
    ci = [i for i in idx if an[i] == cc]
    if not ci:
        stop("%s has no %s atom" % (nm, cc))
    ci = ci[0]
    names, els = [an[i] for i in idx], [el[i] for i in idx]
    draw = {st: X[st][idx] for st in ("R", "TS")}
    if bonds(els, draw["R"]) != bonds(els, draw["TS"]):
        stop("%s: the side chain's bonds differ between R and TS" % nm)
    near = {}
    for st in ("R", "TS"):
        d = np.linalg.norm(X[st][QM0:QM0 + NQM] - X[st][ci], axis=1)
        near[st] = [NAMES[int(d.argmin())], round(float(d.min()), 2)]
    mvc = float(np.linalg.norm(X["TS"][ci] - X["R"][ci]))
    mvx = float(np.linalg.norm(draw["TS"] - draw["R"], axis=1).max())
    s, c, csd = LIT[nm]
    say("%-8s %-2s  R %s %.2f A, TS %s %.2f A   centre moved %.2f A, largest atom %.2f A"
        % (nm, cc, near["R"][0], near["R"][1], near["TS"][0], near["TS"][1], mvc, mvx))
    rows.append({"lab": nm, "res": rn, "seq": seq, "chain": ch, "q": z, "centre": cc,
                 "names": names, "el": els, "bonds": bonds(els, draw["R"]),
                 "ci": idx.index(ci),
                 "R": np.round(draw["R"], 3).tolist(), "TS": np.round(draw["TS"], 3).tolist(),
                 "nearR": near["R"], "nearTS": near["TS"],
                 "moveCentre": round(mvc, 2), "moveMax": round(mvx, 2),
                 "szef": s, "clae": None if c is None else [c, csd]})

order = sorted(LIT, key=lambda k: LIT[k][0])
for r in rows:
    r["rank"] = order.index(r["lab"]) + 1
say("\nliterature rank by Szefczyk 2004 Table 1 (most barrier-lowering first): " + ", ".join(order))
json.dump({"frame": "42436", "residues": rows}, open(os.path.join(WORK, "residues.json"), "w"), indent=1)
open(os.path.join(WORK, "RESIDUES_LOG.txt"), "w").write("\n".join(log) + "\n")
print("\nwrote residues.json (%d residues) and RESIDUES_LOG.txt" % len(rows))
