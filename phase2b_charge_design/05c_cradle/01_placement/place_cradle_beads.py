#!/usr/bin/env python3
# 05c cradle Step 1: adaptive placement of NEUTRAL LJ steric beads (q=0) forming a confining pocket
# around the reacting C1-C6 region, WITHOUT clashing the substrate or the fixed charge surrogates.
#
# Why adaptive (not a fixed cage): the substrate ring sits right at the C1-C6 midpoint, so beads dropped
# at fixed offsets clash (verified: symmetric 6-cage had 2 hard clashes). Instead each candidate bead is
# pushed OUTWARD along its direction until it clears all atoms by >= margin, then kept only if it lands in
# a genuinely CONFINING position (outboard of the carbons / off the reacting-plane faces).
#
# Target 6 beads covering 3 opening DOFs:
#   DOF-1 axial separation : 2 beads outboard beyond C1 and beyond C6 (along C1-C6 axis)
#   DOF-2 out-of-plane swing: 2 beads above/below the reacting plane at the midpoint
#   DOF-3 in-plane splay    : 2 beads lateral (in plane, perpendicular to axis)
# Beads that cannot be placed within max_push without clashing are dropped (reported).
#
# Beads are NEUTRAL (q=0) carbon-like LJ (sigma 3.40, eps 0.086) -> a soft hydrophobic wall, like the
# Val/Ile steric pocket in the enzyme. FROZEN. They exert ONLY LJ on the substrate nuclei.

import sys, math
import numpy as np

def load_xyz(path):
    L=open(path).read().splitlines(); n=int(L[0].split()[0])
    els=[L[2+i].split()[0] for i in range(n)]
    xyz=np.array([[float(v) for v in L[2+i].split()[1:4]] for i in range(n)])
    return els,xyz

def main(reactant_xyz, charge_coords_xyz, out_prefix, margin=2.3, max_push=6.0, step=0.1):
    _,R = load_xyz(reactant_xyz)          # 24-atom substrate
    C1,C6,C4 = R[0],R[12],R[8]
    mid=(C1+C6)/2
    axis=(C6-C1)/np.linalg.norm(C6-C1)
    normal=np.cross(C6-C1, C4-C1); normal/=np.linalg.norm(normal)
    inplane=np.cross(normal,axis); inplane/=np.linalg.norm(inplane)

    # obstacles = substrate atoms + charge-surrogate atoms (if a charge design is supplied)
    obstacles=[R[i] for i in range(24)]
    if charge_coords_xyz:
        _,CC=load_xyz(charge_coords_xyz)
        for i in range(24,len(CC)): obstacles.append(CC[i])   # the MM surrogate atoms
    obstacles=np.array(obstacles)

    def clears(p):
        return min(np.linalg.norm(p-o) for o in obstacles) >= margin
    def push_out(origin, direction):
        # start at origin, step outward along direction until clear (or give up at max_push)
        for k in range(1,int(max_push/step)+1):
            p=origin+direction*(k*step)
            if clears(p): return p, k*step
        return None,None

    # candidate beads: (name, origin, outward direction, DOF)
    half=np.linalg.norm(C6-C1)/2
    cands=[
      ("axial_C1", C1, -axis,    "DOF1-axial"),
      ("axial_C6", C6, +axis,    "DOF1-axial"),
      ("oop_up",   mid, +normal, "DOF2-outofplane"),
      ("oop_down", mid, -normal, "DOF2-outofplane"),
      ("lat_A",    mid, +inplane,"DOF3-inplane"),
      ("lat_B",    mid, -inplane,"DOF3-inplane"),
    ]
    placed=[]; dropped=[]
    for name,origin,direction,dof in cands:
        p,dist = push_out(origin, direction)
        if p is None:
            dropped.append((name,dof)); continue
        # confinement strength = distance from the bead to the nearer reacting carbon it guards
        dC1=np.linalg.norm(p-C1); dC6=np.linalg.norm(p-C6)
        guard=min(dC1,dC6)
        # per-bead sigma so the wall REACHES the guarded carbon (guard/sigma ~1.1), capped physical max 5.0
        sigma=min(max(guard/1.1, 3.40), 5.0)   # never smaller than carbon 3.40, never bigger than 5.0
        placed.append((name,p,dof,guard,dist,sigma))

    # write bead xyz (element X = dummy/LJ bead marker; we will map to a real LJ in the runner)
    with open(out_prefix+"_beads.xyz","w") as fh:
        fh.write("%d\ncradle LJ beads (neutral) name x y z sigma eps=0.086\n"%len(placed))
        for name,p,dof,guard,push,sigma in placed:
            fh.write("X  %12.6f %12.6f %12.6f   %.3f\n"%(p[0],p[1],p[2],sigma))
    print("PLACED %d/6 beads (margin %.1f A):"%(len(placed),margin))
    for name,p,dof,guard,push,sigma in placed:
        ratio=guard/sigma
        print("  %-9s %-16s guard %.2f A, sigma %.2f (guard/sigma %.2f) %s"%(
            name,dof,guard,sigma,ratio,"STRONG" if ratio<1.15 else "moderate"))
    if dropped:
        print("DROPPED (could not clear within %.1f A):"%max_push)
        for name,dof in dropped: print("  %-9s %s"%(name,dof))
    # DOF coverage summary
    covered=set(d.split("-")[0] for _,_,d,_,_,_ in placed)
    print("DOF coverage: %s  (%d/3 opening modes blocked)"%(sorted(covered),len(covered)))
    return placed

if __name__=="__main__":
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument("reactant_xyz")
    ap.add_argument("--charges", default=None, help="charge-surrogate coords.xyz (to avoid clashing them)")
    ap.add_argument("--out", default="cradle")
    ap.add_argument("--margin", type=float, default=2.3)
    a=ap.parse_args()
    main(a.reactant_xyz, a.charges, a.out, a.margin)
