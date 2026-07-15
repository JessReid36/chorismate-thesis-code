#!/usr/bin/env python3
# Phase 2 - extract charged active-site residue positions in the TS frame, using attempt_3's own
# active-site mapping (step79a) applied to the step85 NEB-TS full-system PDB. Charged-group center
# found geometrically. (Fixed: exclude self when counting bonded neighbours.)
import sys, numpy as np

PDB = sys.argv[1]
SELECTED = [
    ("Arg7",  +1, (4231,4248)),
    ("Lys60", +1, (3019,3034)),
    ("Arg63", +1, (3067,3084)),
    ("Glu78", -1, (5371,5379)),
    ("Arg90", +1, (5552,5569)),
    ("Arg116",+1, (5989,6006)),
]
def element_of(ln):
    el = ln[76:78].strip()
    if not el:
        nm = ln[12:16].strip(); el = "".join(c for c in nm if c.isalpha())[:1]
    return el.upper()
atoms=[]
for ln in open(PDB):
    if not ln.startswith(("ATOM","HETATM")): continue
    try:
        serial=int(ln[6:11]); x,y,z=float(ln[30:38]),float(ln[38:46]),float(ln[46:54])
    except ValueError:
        continue
    atoms.append((serial-1, element_of(ln), np.array([x,y,z])))
by_idx={a[0]:(a[1],a[2]) for a in atoms}
def charged_center(name, lo, hi):
    idxs=[i for i in range(lo,hi+1) if i in by_idx]
    els={i:by_idx[i][0] for i in idxs}; pos={i:by_idx[i][1] for i in idxs}
    def bonded(i,j): return i!=j and np.linalg.norm(pos[i]-pos[j])<1.8   # exclude self
    if name.startswith("Arg"):
        for i in idxs:
            if els[i]=="C" and sum(1 for j in idxs if els[j]=="N" and bonded(i,j))>=2:
                grp=[pos[i]]+[pos[j] for j in idxs if els[j]=="N" and bonded(i,j)]
                return np.mean(grp,axis=0),"guanidinium"
    elif name.startswith("Lys"):
        for i in idxs:
            if els[i]=="N" and sum(1 for j in idxs if els[j] in("C","N","O") and bonded(i,j))<=1:
                return pos[i],"NZ"
    elif name.startswith(("Glu","Asp")):
        for i in idxs:
            if els[i]=="C" and sum(1 for j in idxs if els[j]=="O" and bonded(i,j))==2:
                grp=[pos[i]]+[pos[j] for j in idxs if els[j]=="O" and bonded(i,j)]
                return np.mean(grp,axis=0),"carboxylate"
    return np.mean([pos[i] for i in idxs],axis=0),"CENTROID-FALLBACK"
rows=[]; fallbacks=0
print(f"{'residue':>8} {'sign':>4} {'center-type':>20}  charge-center xyz")
for nm,s,(lo,hi) in SELECTED:
    cc,ctype=charged_center(nm,lo,hi)
    if "FALLBACK" in ctype: fallbacks+=1
    rows.append((nm,s,cc))
    print(f"{nm:>8} {s:+4d} {ctype:>20}  ({cc[0]:.3f}, {cc[1]:.3f}, {cc[2]:.3f})")
if fallbacks:
    print(f"\n*** WARNING: {fallbacks} residue(s) used centroid fallback ***")
else:
    print(f"\nAll {len(rows)} charged groups identified geometrically (no fallbacks).")
np.savez("enzyme_residues.npz",
         names=np.array([r[0] for r in rows]),
         signs=np.array([r[1] for r in rows]),
         centers=np.array([r[2] for r in rows]),
         dmin=np.array([0.0]*len(rows)))
with open("residue_meta.tsv","w") as fh:
    fh.write("idx\tname\tsign\tx\ty\tz\n")
    for i,(nm,s,cc) in enumerate(rows):
        fh.write(f"{i}\t{nm}\t{int(s)}\t{cc[0]:.4f}\t{cc[1]:.4f}\t{cc[2]:.4f}\n")
print(f"saved enzyme_residues.npz + residue_meta.tsv ({len(rows)} residues, TS frame, canonical IDs)")
