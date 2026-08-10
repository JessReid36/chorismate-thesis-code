#!/usr/bin/env python3
"""Clean-pipeline neutraliser test: two +1 counter-charges on the carboxylates, LJ-standoff-respecting."""
import os
import numpy as np

HERE   = os.path.dirname(os.path.abspath(__file__))
INP    = os.path.join(HERE, "..", "inputs")
REACT  = "reactant.xyz"
CARBOXYLATES = [(4, 5, 6), (21, 22, 23)]
QSIG   = 3.40
SIG    = {"C": 3.40, "N": 3.25, "O": 2.96, "H": 1.06}
TWO16  = 2.0 ** (1.0 / 6.0)
STEP   = 0.05
START  = 2.0


def load_sub():
    L = open(os.path.join(INP, REACT)).read().splitlines()[2:26]
    els = [l.split()[0] for l in L]
    R = np.array([[float(v) for v in l.split()[1:4]] for l in L])
    return els, R


def main():
    els, R = load_sub()
    rmin_atom = np.array([TWO16 * 0.5 * (QSIG + SIG[e]) for e in els])
    sites = []
    for (cC, oA, oB) in CARBOXYLATES:
        b = (R[oA] - R[cC]) + (R[oB] - R[cC]); b = b / np.linalg.norm(b)
        mid = 0.5 * (R[oA] + R[oB])
        d = START
        while True:
            p = mid + d * b
            if (np.linalg.norm(R - p, axis=1) - rmin_atom >= 0).all():
                break
            d += STEP
        sites.append(p)
        dmin = np.linalg.norm(R - p, axis=1).min()
        print("carboxylate C%d: +1 placed %.2f A along bisector | min dist to substrate %.2f A"
              % (cC, d, dmin))
    sites = np.array(sites)
    sub_lines = open(os.path.join(INP, REACT)).read().splitlines()[2:26]
    with open(os.path.join(HERE, "design_cradle2_coords.xyz"), "w") as f:
        f.write("%d\ncradle: 24 substrate + 2 counter-charges (+1 each, LJ-standoff), 0 catalytic\n" % (24 + 2))
        for l in sub_lines:
            f.write(l.rstrip() + "\n")
        for p in sites:
            f.write("C   %14.8f %14.8f %14.8f\n" % (p[0], p[1], p[2]))
    with open(os.path.join(HERE, "design_cradle2_charges.txt"), "w") as f:
        f.write(",".join(["0.0000"] * 24 + ["1.0000", "1.0000"]))
    print("system net charge: %.1f ; QM charge stays -2" % (-2 + 2))
    print("wrote design_cradle2_{coords.xyz,charges.txt}")


if __name__ == "__main__":
    main()
