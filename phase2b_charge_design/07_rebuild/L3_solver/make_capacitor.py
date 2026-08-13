#!/usr/bin/env python3
"""Parallel-plate CAPACITOR design: tens of charges making a UNIFORM -z field over the whole substrate."""
import os, sys
import numpy as np

HERE  = os.path.dirname(os.path.abspath(__file__))
INP   = os.path.join(HERE, "..", "inputs")
REACT = "reactant.xyz"
FTARGET = float(sys.argv[1]) if len(sys.argv) > 1 else 0.02
NPLATE  = int(sys.argv[2]) if len(sys.argv) > 2 else 5
PLATE_SPAN = float(sys.argv[3]) if len(sys.argv) > 3 else 12.0
BOHR = 0.529177
REACTING = [7, 8, 0, 12]
SIG = {"C": 3.40, "N": 3.25, "O": 2.96, "H": 1.06}
QSIG = 3.40
TWO16 = 2.0 ** (1.0/6.0)
STEP = 0.1


def load_sub():
    L = open(os.path.join(INP, REACT)).read().splitlines()[2:26]
    els = [l.split()[0] for l in L]
    R = np.array([[float(v) for v in l.split()[1:4]] for l in L])
    return els, R


def make_plate(center, span, n, zdir, R, rmin_atom):
    xs = np.linspace(-span/2, span/2, n)
    ys = np.linspace(-span/2, span/2, n)
    d = 0.0
    while True:
        c = center + d*zdir
        pts = np.array([[c[0]+x, c[1]+y, c[2]] for x in xs for y in ys])
        ok = True
        for p in pts:
            if (np.linalg.norm(R - p, axis=1) - rmin_atom < 0).any():
                ok = False; break
        if ok:
            return pts
        d += STEP


def field_z_at(centroid, plus_pts, minus_pts, q):
    F = np.zeros(3)
    for p in plus_pts:
        v=(centroid-p)/BOHR; r=np.linalg.norm(v); F += (+q)*v/r**3
    for p in minus_pts:
        v=(centroid-p)/BOHR; r=np.linalg.norm(v); F += (-q)*v/r**3
    return F


def main():
    els, R = load_sub()
    rmin_atom = np.array([TWO16*0.5*(QSIG+SIG[e]) for e in els])
    centroid = R[REACTING].mean(axis=0)
    zhat = np.array([0,0,1.0])
    plus_pts  = make_plate(centroid, PLATE_SPAN, NPLATE, -zhat, R, rmin_atom)
    minus_pts = make_plate(centroid, PLATE_SPAN, NPLATE, +zhat, R, rmin_atom)
    F1 = field_z_at(centroid, plus_pts, minus_pts, 1.0)
    q = -FTARGET / F1[2]
    Ffinal = field_z_at(centroid, plus_pts, minus_pts, q)
    print("plates: %d + %d = %d charges, span %.1f A, per-charge q=%.4f" %
          (len(plus_pts), len(minus_pts), len(plus_pts)+len(minus_pts), PLATE_SPAN, q))
    print("field at centroid: Fz=%.4f a.u. (target %.3f)" % (Ffinal[2], -FTARGET))
    print("uniformity (-z field at reacting atoms):")
    for a in REACTING:
        Fa = field_z_at(R[a], plus_pts, minus_pts, q)
        print("   atom %2d: Fz=%.4f" % (a, Fa[2]))
    minplate = min(np.linalg.norm(R - p, axis=1).min() for p in np.vstack([plus_pts, minus_pts]))
    print("min plate-substrate distance: %.2f A" % minplate)
    cradle = None
    cp = os.path.join(HERE, "design_cradle2_coords.xyz")
    if os.path.exists(cp):
        L=open(cp).read().splitlines(); n=int(L[0].split()[0])
        cradle=np.array([[float(v) for v in l.split()[1:4]] for l in L[2+24:2+n]])
    sub=open(os.path.join(INP,REACT)).read().splitlines()[2:26]
    tag="cap_f%s_n%d" % (("%.3f"%FTARGET).replace(".","p"), NPLATE)
    ncr=len(cradle) if cradle is not None else 0
    allpts=np.vstack([plus_pts, minus_pts])
    with open(os.path.join(HERE,"design_%s_coords.xyz"%tag),"w") as f:
        f.write("%d\ncapacitor uniform -z field (%dx%d plates) + cradle(%d)\n"%(24+ncr+len(allpts),NPLATE,NPLATE,ncr))
        for l in sub: f.write(l.rstrip()+"\n")
        if cradle is not None:
            for p in cradle: f.write("C   %14.8f %14.8f %14.8f\n"%(p[0],p[1],p[2]))
        for p in allpts: f.write("C   %14.8f %14.8f %14.8f\n"%(p[0],p[1],p[2]))
    with open(os.path.join(HERE,"design_%s_charges.txt"%tag),"w") as f:
        vals=["0.0000"]*24 + (["1.0000"]*ncr) + ["%.4f"%q]*len(plus_pts) + ["%.4f"%(-q)]*len(minus_pts)
        f.write(",".join(vals))
    print("wrote design_%s_* (%d plate charges + %d cradle)"%(tag,len(allpts),ncr))


if __name__=="__main__":
    main()
