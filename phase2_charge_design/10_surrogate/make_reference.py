#!/usr/bin/env python3
# Build the field-off reference for the surrogate calculation.
#
# The catalytic effect is (TS - R) with the design present, minus (TS - R) without it. Using the
# bare 24-atom substrate as the second term works -- the frozen surrogates cancel within the
# design's own TS-minus-R subtraction -- but it leaves basis-set superposition error uncancelled,
# because the two calculations have different basis sizes.
#
# This reference keeps the identical 44 atoms, identical electron count and identical basis, and
# removes only the field, by moving each surrogate far from the substrate. The two are placed
# equidistant from the substrate centroid on opposite sides, so their monopole interactions with
# the -2 substrate cancel by symmetry instead of merely decaying. What survives is a small
# quadrupolar residual, and only the part of it that differs between the reactant and TS
# geometries enters the barrier -- second order on a substrate whose shape barely changes.
#
# Surrogate internal geometry is translated rigidly, not rotated, so their mutual interaction is
# identical in the reactant and TS references and cancels exactly.
#
# Build this from the INPUT geometry, not from a converged output, so the reference and the
# field-on run start from the same substrate and differ only in where the surrogates sit.

import math
import sys

N_SUB = 24
SEP = 15.0                       # A from the substrate centroid, each side
SRC = sys.argv[1] if len(sys.argv) > 1 else "surr_ropt.inp"
OUT = sys.argv[2] if len(sys.argv) > 2 else "surr_ropt_ref.inp"

N_MGDN, N_ACET = 13, 7


def read_inp(path):
    head, atoms, tail, where = [], [], [], "head"
    for line in open(path):
        s = line.rstrip("\n")
        if where == "head":
            head.append(s)
            if s.strip().startswith("*") and "xyz" in s:
                where = "atoms"
            continue
        if where == "atoms":
            if s.strip().startswith("*"):
                tail.append(s); where = "tail"; continue
            p = s.split()
            if len(p) == 4:
                atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
            continue
        tail.append(s)
    return head, atoms, tail


def centroid(atoms):
    n = len(atoms)
    return tuple(sum(a[i+1] for a in atoms)/n for i in range(3))


def translate(group, shift):
    return [(a[0], a[1]+shift[0], a[2]+shift[1], a[3]+shift[2]) for a in group]


def main():
    head, atoms, tail = read_inp(SRC)
    if len(atoms) != N_SUB + N_MGDN + N_ACET:
        print("expected %d atoms, found %d" % (N_SUB+N_MGDN+N_ACET, len(atoms))); sys.exit(1)

    sub = atoms[:N_SUB]
    mgdn = atoms[N_SUB:N_SUB+N_MGDN]
    acet = atoms[N_SUB+N_MGDN:]
    cen = centroid(sub)

    spans = [max(a[i+1] for a in sub) - min(a[i+1] for a in sub) for i in range(3)]
    axis = spans.index(min(spans))
    d = [0.0, 0.0, 0.0]; d[axis] = 1.0
    print("substrate spans %.2f %.2f %.2f -> displacing along axis %d" % (spans[0], spans[1], spans[2], axis))

    mg_c = (mgdn[0][1], mgdn[0][2], mgdn[0][3])
    ac_c = (acet[0][1], acet[0][2], acet[0][3])
    tgt_mg = tuple(cen[i] + SEP*d[i] for i in range(3))
    tgt_ac = tuple(cen[i] - SEP*d[i] for i in range(3))

    mgdn_new = translate(mgdn, tuple(tgt_mg[i]-mg_c[i] for i in range(3)))
    acet_new = translate(acet, tuple(tgt_ac[i]-ac_c[i] for i in range(3)))

    def mindist(g):
        return min(math.sqrt(sum((g[k][i+1]-s[i+1])**2 for i in range(3)))
                   for k in range(len(g)) for s in sub)

    print("  methylguanidinium charge centre moved to %.2f %.2f %.2f, closest substrate atom %.2f A"
          % (tgt_mg[0], tgt_mg[1], tgt_mg[2], mindist(mgdn_new)))
    print("  acetate           charge centre moved to %.2f %.2f %.2f, closest substrate atom %.2f A"
          % (tgt_ac[0], tgt_ac[1], tgt_ac[2], mindist(acet_new)))

    r_mg = math.sqrt(sum((tgt_mg[i]-cen[i])**2 for i in range(3)))
    r_ac = math.sqrt(sum((tgt_ac[i]-cen[i])**2 for i in range(3)))
    resid = 332.06*(-2)*(1.0)/r_mg + 332.06*(-2)*(-1.0)/r_ac
    print("  residual monopole term at the centroid: %+.4f kcal/mol (cancels by symmetry)" % resid)

    fh = open(OUT, "w")
    for s in head:
        fh.write(s + "\n")
    for a in sub + mgdn_new + acet_new:
        fh.write("%-2s %16.8f %16.8f %16.8f\n" % (a[0], a[1], a[2], a[3]))
    for s in tail:
        fh.write(s + "\n")
    fh.close()
    print("\nwrote %s: same %d atoms, same charge, field removed" % (OUT, len(atoms)))


if __name__ == "__main__":
    main()
