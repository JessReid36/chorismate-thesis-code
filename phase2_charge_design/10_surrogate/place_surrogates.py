#!/usr/bin/env python3
# Replace each designed unit charge with the real ionised group it represents: methylguanidinium
# for +1 (an arginine side chain), acetate for -1 (a glutamate side chain).
#
# Bare point charges carry no Pauli repulsion, so a substrate carboxylate oxygen falls onto a
# designed +1 without limit (0.879 A in production, 0.64 A inside its own van der Waals radius).
# Real groups carry core electrons, so the approach stops at a salt-bridge distance instead.
#
# Each surrogate is placed with its formal charge centre on the designed grid point and its
# charged end facing the substrate centroid, which is how a side chain approaches from the
# protein. All surrogate atoms are frozen; only the substrate relaxes.
#
# Pure Python: numpy on the login node needs the thread environment set and this is only
# rotations and distances.

import math
import os
import sys

N_SUB = 24
SUB_INP = "../08_validate_relaxed/ropt_charged.inp"
PC_FILE = "../08_validate_relaxed/design_charges.pc"
OUT = "surr_ropt.inp"

BOND = {"CN": 1.33, "NC": 1.46, "NH": 1.01, "CH": 1.09, "CO": 1.25, "CC": 1.52}
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "H": 1.20}


def unit(v):
    n = math.sqrt(sum(c*c for c in v))
    return tuple(c/n for c in v)


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def rotation_to(a_local, b_target):
    a, b = unit(a_local), unit(b_target)
    v = cross(a, b)
    c = sum(x*y for x, y in zip(a, b))
    s = math.sqrt(sum(x*x for x in v))
    if s < 1e-9:
        return [[1,0,0],[0,1,0],[0,0,1]] if c > 0 else [[-1,0,0],[0,-1,0],[0,0,1]]
    vx = [[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]]
    R = [[0.0]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            vv = sum(vx[i][k]*vx[k][j] for k in range(3))
            R[i][j] = (1.0 if i == j else 0.0) + vx[i][j] + vv*((1-c)/(s*s))
    return R


def apply(R, p):
    return tuple(sum(R[i][j]*p[j] for j in range(3)) for i in range(3))


def methylguanidinium():
    d = BOND["CN"]
    atoms = [("C", (0.0, 0.0, 0.0))]
    n1 = (0.0, 0.0, -d)
    n2 = (d*math.sin(math.radians(60)), 0.0, d*math.cos(math.radians(60)))
    n3 = (-n2[0], 0.0, n2[2])
    atoms += [("N", n1), ("N", n2), ("N", n3)]
    m_dir = (math.sin(math.radians(120)), 0.0, math.cos(math.radians(120)))
    cm = tuple(n1[i] + BOND["NC"]*m_dir[i] for i in range(3))
    atoms.append(("C", cm))
    h_dir = (-m_dir[0], 0.0, m_dir[2])
    atoms.append(("H", tuple(n1[i] + BOND["NH"]*h_dir[i] for i in range(3))))
    for n in (n2, n3):
        dnc = unit((-n[0], 0.0, -n[2]))
        perp = unit((-dnc[2], 0.0, dnc[0]))
        for sgn in (1, -1):
            hd = tuple(-0.5*dnc[i] + sgn*math.sin(math.radians(60))*perp[i] for i in range(3))
            atoms.append(("H", tuple(n[i] + BOND["NH"]*hd[i] for i in range(3))))
    ax = unit(tuple(cm[i]-n1[i] for i in range(3)))
    e1 = unit(cross(ax, (0.0, 1.0, 0.0))); e2 = cross(ax, e1)
    for k in range(3):
        th = 2*math.pi*k/3
        hd = tuple(math.cos(math.radians(70.5))*ax[i]
                   + math.sin(math.radians(70.5))*(math.cos(th)*e1[i]+math.sin(th)*e2[i])
                   for i in range(3))
        atoms.append(("H", tuple(cm[i] + BOND["CH"]*hd[i] for i in range(3))))
    return atoms


def acetate():
    d = BOND["CO"]
    half = math.radians(63.0)
    atoms = [("C", (0.0, 0.0, 0.0))]
    o1 = (d*math.sin(half), 0.0, d*math.cos(half))
    atoms += [("O", o1), ("O", (-o1[0], 0.0, o1[2]))]
    cm = (0.0, 0.0, -BOND["CC"])
    atoms.append(("C", cm))
    ax = (0.0, 0.0, -1.0)
    for k in range(3):
        th = 2*math.pi*k/3
        hd = tuple(math.cos(math.radians(70.5))*ax[i]
                   + math.sin(math.radians(70.5))*(math.cos(th)*(1,0,0)[i]+math.sin(th)*(0,1,0)[i])
                   for i in range(3))
        atoms.append(("H", tuple(cm[i] + BOND["CH"]*hd[i] for i in range(3))))
    return atoms


def read_substrate(path):
    out, started = [], False
    for line in open(path):
        s = line.strip()
        if s.startswith("*") and "xyz" in s:
            started = True
            continue
        if started:
            if s.startswith("*"):
                break
            p = s.split()
            if len(p) == 4:
                out.append((p[0], (float(p[1]), float(p[2]), float(p[3]))))
        if len(out) == N_SUB:
            break
    return out


def read_pc(path):
    out = []
    for line in open(path).read().splitlines()[1:]:
        p = line.split()
        if len(p) == 4:
            out.append((float(p[0]), (float(p[1]), float(p[2]), float(p[3]))))
    return out


def main():
    sub = read_substrate(SUB_INP)
    if len(sub) != N_SUB:
        print("read %d substrate atoms, expected %d" % (len(sub), N_SUB)); sys.exit(1)
    charges = read_pc(PC_FILE)
    cen = tuple(sum(a[1][i] for a in sub)/len(sub) for i in range(3))
    print("substrate %d atoms, centroid %.3f %.3f %.3f" % (len(sub), cen[0], cen[1], cen[2]))

    placed = []
    for q, pos in charges:
        local = methylguanidinium() if q > 0 else acetate()
        label = "methylguanidinium(+1)" if q > 0 else "acetate(-1)"
        R = rotation_to((0.0, 0.0, 1.0), tuple(cen[i]-pos[i] for i in range(3)))
        grp = [(el, tuple(pos[i] + apply(R, p)[i] for i in range(3))) for el, p in local]
        # closest approach to the substrate
        worst = None
        for el, p in grp:
            for sel, sp in sub:
                d = math.sqrt(sum((a-b)**2 for a, b in zip(p, sp)))
                clear = d - VDW.get(el, 1.7) - VDW.get(sel, 1.7)
                if worst is None or clear < worst[0]:
                    worst = (clear, d, el, sel)
        clear, d, el, sel = worst
        flag = "  CLASH" if clear < -0.4 else ("  close" if clear < 0.0 else "")
        print("  %-22s at %.3f %.3f %.3f   %2d atoms   closest %s...%s %.3f A (%+.3f vs vdW)%s"
              % (label, pos[0], pos[1], pos[2], len(grp), el, sel, d, clear, flag))
        placed.append(grp)

    n_surr = sum(len(g) for g in placed)
    first = N_SUB
    ranges = "{C %d:%d C}" % (first, first + n_surr - 1)

    fh = open(OUT, "w")
    fh.write("! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM Opt\n")
    fh.write("%maxcore 3000\n%pal nprocs 8 end\n")
    fh.write("%cpcm epsilon 4.0 end\n")
    fh.write("%scf MaxIter 300 end\n")
    fh.write("%geom MaxIter 200\n  Constraints\n    " + ranges + "\n  end\nend\n")
    fh.write("* xyz -2 1\n")
    for el, p in sub:
        fh.write("%-2s %16.8f %16.8f %16.8f\n" % (el, p[0], p[1], p[2]))
    for g in placed:
        for el, p in g:
            fh.write("%-2s %16.8f %16.8f %16.8f\n" % (el, p[0], p[1], p[2]))
    fh.write("*\n")
    fh.close()
    print("\nwrote %s: %d atoms total, surrogate atoms %d-%d frozen"
          % (OUT, N_SUB + n_surr, first, first + n_surr - 1))
    print("total charge -2 = substrate -2, plus +1 and -1 from the surrogates")


if __name__ == "__main__":
    main()
