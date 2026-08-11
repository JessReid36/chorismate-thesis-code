#!/usr/bin/env python3
"""Hand-built ORIENTED charge pair to make a -z field at the reacting region (OEEF catalysis test)."""
import os, sys
import numpy as np

HERE  = os.path.dirname(os.path.abspath(__file__))
INP   = os.path.join(HERE, "..", "inputs")
REACT = "reactant.xyz"
Q     = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
QSIG  = 3.40
SIG   = {"C": 3.40, "N": 3.25, "O": 2.96, "H": 1.06}
TWO16 = 2.0 ** (1.0 / 6.0)
STEP  = 0.05
REACTING = [7, 8, 0, 12]
CATALYTIC_AXIS = np.array([0.0, 0.0, -1.0])


def load_sub():
    L = open(os.path.join(INP, REACT)).read().splitlines()[2:26]
    els = [l.split()[0] for l in L]
    R = np.array([[float(v) for v in l.split()[1:4]] for l in L])
    return els, R


def march_out(start, direction, R, rmin_atom):
    d = 0.0
    while True:
        p = start + d * direction
        if (np.linalg.norm(R - p, axis=1) - rmin_atom >= 0).all():
            return p, d
        d += STEP


def main():
    els, R = load_sub()
    rmin_atom = np.array([TWO16 * 0.5 * (QSIG + SIG[e]) for e in els])
    centroid = R[REACTING].mean(axis=0)
    zhat = np.array([0.0, 0.0, 1.0])
    p_plus,  d1 = march_out(centroid, -zhat, R, rmin_atom)
    p_minus, d2 = march_out(centroid, +zhat, R, rmin_atom)
    print("reacting centroid: %s" % np.round(centroid, 2))
    print("+%.2f placed %.2f A along -z, min dist to substrate %.2f A" %
          (Q, d1, np.linalg.norm(R - p_plus, axis=1).min()))
    print("-%.2f placed %.2f A along +z, min dist to substrate %.2f A" %
          (Q, d2, np.linalg.norm(R - p_minus, axis=1).min()))
    BOHR = 0.529177
    def fld(qc, pos):
        v = (centroid - pos) / BOHR; r = np.linalg.norm(v); return qc * v / r**3
    F = fld(+Q, p_plus) + fld(-Q, p_minus)
    print("field at centroid: Fz = %.4f a.u. (target ~ -0.02 for ~-10.6 kcal/mol)" % F[2])
    cradle = None
    cp = os.path.join(HERE, "design_cradle2_coords.xyz")
    if os.path.exists(cp):
        L = open(cp).read().splitlines(); n = int(L[0].split()[0])
        cradle = np.array([[float(v) for v in l.split()[1:4]] for l in L[2 + 24:2 + n]])
    sub_lines = open(os.path.join(INP, REACT)).read().splitlines()[2:26]
    tag = "oez%s" % ("%.2f" % Q).replace(".", "p")
    ncr = len(cradle) if cradle is not None else 0
    with open(os.path.join(HERE, "design_%s_coords.xyz" % tag), "w") as f:
        f.write("%d\noriented -z pair (+-%.2f) + cradle(%d) for OEEF catalysis test\n" % (24 + ncr + 2, Q, ncr))
        for l in sub_lines:
            f.write(l.rstrip() + "\n")
        if cradle is not None:
            for p in cradle:
                f.write("C   %14.8f %14.8f %14.8f\n" % (p[0], p[1], p[2]))
        f.write("C   %14.8f %14.8f %14.8f\n" % (p_plus[0], p_plus[1], p_plus[2]))
        f.write("C   %14.8f %14.8f %14.8f\n" % (p_minus[0], p_minus[1], p_minus[2]))
    with open(os.path.join(HERE, "design_%s_charges.txt" % tag), "w") as f:
        vals = ["0.0000"] * 24
        if cradle is not None:
            vals += ["1.0000"] * ncr
        vals += ["%.4f" % Q, "%.4f" % (-Q)]
        f.write(",".join(vals))
    print("wrote design_%s_* (24 substrate + %d cradle + oriented +-%.2f pair)" % (tag, ncr, Q))


if __name__ == "__main__":
    main()
