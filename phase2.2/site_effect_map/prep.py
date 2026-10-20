"""prep.py - the site-effect map's data: sites, first-order effects, frame 20000's substrate (R, TS) and the residues.
Usage: python3 prep.py <results_repo> <work_dir>   (reads <work_dir>/residues.json from residues.py; writes <work_dir>/data.json)
Sites and admitted charges: phase2.2/cs_stage5/exploratory/dev_sites_a2.tsv (Stage 5 Amendment 2); effects a = V_TS - V_R
per unit charge: exploratory/dev_A/A_dev.tsv; substrate: phase2.2/aligned/frame_20000 (the alignment reference)."""
import json, os, sys, numpy as np
ROOT, WORK = sys.argv[1], sys.argv[2]
R = ROOT + "/phase2.2"
E = R + "/cs_stage5/exploratory"
NAMES = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
def xyz(p):
    L = open(p).read().splitlines(); n = int(L[0])
    return [l.split()[0] for l in L[2:2+n]], np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2+n]])
el, XR = xyz(R + "/aligned/frame_20000/reactant_qm.xyz")
_, XT = xyz(R + "/aligned/frame_20000/transition_state_qm.xyz")
_, XP = xyz(R + "/aligned/frame_20000/product_qm.xyz")
assert [n.rstrip("0123456789") for n in NAMES] == el, (el)
cov = {"C": 0.76, "H": 0.31, "O": 0.66}
def bonds(X):
    out = set()
    for i in range(24):
        for j in range(i+1, 24):
            if np.linalg.norm(X[i]-X[j]) < cov[el[i]] + cov[el[j]] + 0.40: out.add((i, j))
    return out
bR, bT, bP = bonds(XR), bonds(XT), bonds(XP)
breaking = sorted(bR - bP); forming = sorted(bP - bR)
print("breaking", [(NAMES[i], NAMES[j], round(float(np.linalg.norm(XT[i]-XT[j])),2)) for i,j in breaking])
print("forming", [(NAMES[i], NAMES[j], round(float(np.linalg.norm(XT[i]-XT[j])),2), round(float(np.linalg.norm(XR[i]-XR[j])),2)) for i,j in forming])
# admitted charges from the Stage 5 Amendment 2 development grid (same sites, positions and order as dev_sites.tsv)
A2 = E + "/dev_sites_a2.tsv"  # committed: results 269f54b
import hashlib; assert hashlib.sha256(open(A2, "rb").read()).hexdigest() == "43be59bc86273c2a7f063965a429dc7dd8aa39957080f7b0de57548561b645e8"
rows = [l.rstrip("\n").split("\t") for l in open(A2) if not l.startswith("#")]
h = rows[0]; rows = rows[1:]
A = np.loadtxt(E + "/dev_A/A_dev.tsv", comments="#")
assert A.shape == (len(rows), 30)
dn = ["d_" + n for n in NAMES]
D = np.array([[float(r[h.index(k)]) for k in dn] for r in rows])
adm = np.array([r[h.index("q_min_pos")] != "none" for r in rows])
m, sd, lo, hi = A.mean(1), A.std(1, ddof=1), A.min(1), A.max(1)
near = D.argmin(1)
sites = []
for i, r in enumerate(rows):
    sites.append([i, 1 if r[1] == "cation" else 0, round(float(r[2]),3), round(float(r[3]),3), round(float(r[4]),3),
                  round(float(m[i]),4), round(float(sd[i]),3), round(float(lo[i]),3), round(float(hi[i]),3), int(near[i]), round(float(D[i].min()),2),
                  r[h.index("q_min_pos")], 1 if m[i] > 0 else -1, 1 if r[h.index("q_min_neg")] == "-1.00" else 0])
idx = np.where(adm)[0]
from collections import Counter
low = idx[np.argsort(m[idx])]; high = idx[np.argsort(-m[idx])]
nlow = int((m[idx] < 0).sum()); nhigh = int((m[idx] > 0).sum())
print("admit +1:", len(idx), "lower (mean a<0):", nlow, "raise:", nhigh)
print("top50 lower nearest:", Counter(NAMES[near[i]] for i in low[:50]).most_common())
print("top50 raise nearest:", Counter(NAMES[near[i]] for i in high[:50]).most_common())
print("most neg %.2f most pos %.2f; median sd (admitting +1) %.2f, all sites %.2f" % (m[idx].min(), m[idx].max(), np.median(sd[idx]), np.median(sd)))
print("frame range of top lowering site", lo[low[0]], hi[low[0]])
RJ = json.load(open(os.path.join(WORK, "residues.json")))
# Literature catalytic importance of the six residues (read in full from the papers; kcal/mol):
#   Szefczyk et al. 2004, J. Am. Chem. Soc. 126, 16148, Table 1, MP2 differential TS stabilisation (negative lowers
#   the barrier; the paper orders residues by it). Claeyssens et al. 2011, Org. Biomol. Chem. 9, 1578, Table 3,
#   electrostatic TS stabilisation, mean (sd) over 16 B3LYP/CHARMM27 pathways (positive stabilises the TS; Lys60 not reported).
LIT = {"Arg90": (-9.06, 6.0, 0.8), "Arg7": (-5.90, 3.1, 0.5), "Glu78": (-3.57, 2.1, 0.3), "Arg116": (-2.45, 1.0, 0.9),
       "Arg63′": (-1.40, -0.6, 0.5), "Lys60′": (1.39, None, None)}
order = sorted(LIT, key=lambda k: LIT[k][0])
cl = [k for k in order if LIT[k][1] is not None]
assert cl == sorted(cl, key=lambda k: -LIT[k][1]), "Claeyssens 2011 orders the residues differently"
for r in RJ["residues"]:
    s, c, csd = LIT[r["lab"]]; r["rank"] = order.index(r["lab"]) + 1; r["szef"] = s; r["clae"] = None if c is None else [c, csd]
print("literature rank:", order)
assert RJ["frame"] == "20000" and RJ["frames"] == open(E + "/dev_A/A_dev.tsv").read().split("# frames (30): ")[1].split("\n")[0].split()
data = {"names": NAMES, "el": el, "R": np.round(XR,3).tolist(), "TS": np.round(XT,3).tolist(),
        "bR": sorted(bR), "bT": sorted(bT & bR & bP), "breaking": breaking, "forming": forming,
        "sites": sites, "nframes": 30,
        "res": RJ["residues"], "res_frame": RJ["frame"],   # from residues.py: side chains at this frame's R and TS, as the substrate
        "frames": open(E + "/dev_A/A_dev.tsv").read().split("# frames (30): ")[1].split("\n")[0].split()}
json.dump(data, open(os.path.join(WORK, "data.json"), "w"), separators=(",", ":"))
print("json bytes", os.path.getsize(os.path.join(WORK, "data.json")))
