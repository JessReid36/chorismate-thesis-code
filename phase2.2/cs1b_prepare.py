#!/usr/bin/env python3
"""
cs1b_prepare.py - Stage 1 extension S1.5: do ORCA's QM/MM electrostatics between a sphere and the substrate equal the
A matrix's first-order terms (orca_vpot on the stored in vacuo density), with no cutoff?
Builds six more QM/MM single points on frame 41786 (same method, substrate geometry and force-field rows as Stage 1)
and one orca_vpot point file. Stage 1 showed the substrate geometry equals 05_qmmm/20_invacuo/41786_R.xyz exactly and
the q = 0 QM energy equals that in vacuo single point to 2e-10 Eh, so sp_41786_R's density is the q -> 0 density.
Sites (Amber N sphere, as Stage 1):
  C (contact) = Stage 1 sphere 1 (40.832, 18.750, 42.551), nearest atom 3.40 A: q = -0.5, -1 (q = 0, +0.5, +1 exist)
  M (mid)     = on the line centroid -> C, nearest atom 9.000 A (7 of 24 atoms beyond 12 A): q = +1, -1
  F (far)     = same line, nearest atom 13.000 A (all 24 atoms beyond 12 A): q = +1, -1
USAGE  python3 cs1b_prepare.py <results repo> <output folder>
"""
import os, sys
import numpy as np
REPO, OUT = sys.argv[1], sys.argv[2]
ANG2BOHR = 1.8897261339          # ORCA 6.0.1's own factor (s15's note: 1/0.5291772083)
S1 = os.path.join(REPO, "phase2.2/cs_stage1/s1a_q1_N")
xyz = open(os.path.join(S1, "combined.xyz")).read().splitlines()
SUB = [l.split() for l in xyz[2:26]]
X = np.array([[float(v) for v in a[1:4]] for a in SUB])
vac = np.array([[float(v) for v in l.split()[1:4]] for l in open(os.path.join(REPO, "05_qmmm/20_invacuo/41786_R.xyz")).read().splitlines()[2:26]])
assert np.abs(X - vac).max() == 0.0, "Stage 1 substrate is not the in vacuo reactant geometry"
P = open(os.path.join(S1, "site.ORCAFF.prms")).read().splitlines()
i = P.index("$atoms")
ROWS = P[i + 2:i + 2 + 24]
C = np.array([40.832, 18.750, 42.551])
u = (C - X.mean(0)) / np.linalg.norm(C - X.mean(0))
def along(target):
    lo, hi = 0.0, 40.0
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if np.linalg.norm(X - (C + m * u), axis=1).min() < target: lo = m
        else: hi = m
    return C + hi * u
M, F = along(9.0), along(13.0)
SITES = {"C": C, "M": M, "F": F}
INP = open(os.path.join(S1, "job.inp")).read()
def write(name, site, q):
    d = os.path.join(OUT, name); os.makedirs(d)
    p = SITES[site]
    with open(os.path.join(d, "combined.xyz"), "w") as f:
        f.write("25\nframe 41786 R + sphere q=%+.1f at site %s\n" % (q, site))
        for a in SUB: f.write("%-2s %15.8f%15.8f%15.8f\n" % (a[0], float(a[1]), float(a[2]), float(a[3])))
        f.write("%-2s %15.8f%15.8f%15.8f\n" % ("N", p[0], p[1], p[2]))
    with open(os.path.join(d, "site.ORCAFF.prms"), "w") as f:
        f.write("$fftype\nAMBER\n$atoms\n25 1 4\n")
        for l in ROWS: f.write(l + "\n")
        f.write("%6d   %-2s%13.6f%13.6f%13.6f%13.6f%13.6f\n" % (25, "N", q, -0.17, 3.648, -0.085, 3.648))
        f.write("$bonds\n0 2 2\n$angles\n0 3 2\n$dihedrals\n0 4 3\n")
    open(os.path.join(d, "job.inp"), "w").write(INP)
write("s1h_C_qm05", "C", -0.5); write("s1i_C_qm1", "C", -1.0)
write("s1j_M_qp1", "M", 1.0);   write("s1k_M_qm1", "M", -1.0)
write("s1l_F_qp1", "F", 1.0);   write("s1m_F_qm1", "F", -1.0)
with open(os.path.join(OUT, "vpot_points_bohr.xyz"), "w") as f:
    f.write("3\n")
    for k in ("C", "M", "F"):
        p = SITES[k] * ANG2BOHR
        f.write("%.10f %.10f %.10f\n" % tuple(p))
with open(os.path.join(OUT, "sites.tsv"), "w") as f:
    f.write("site\tx_A\ty_A\tz_A\tnearest_atom_A\tatoms_beyond_12A\n")
    for k in ("C", "M", "F"):
        d = np.linalg.norm(X - SITES[k], axis=1)
        f.write("%s\t%.6f\t%.6f\t%.6f\t%.4f\t%d\n" % (k, *SITES[k], d.min(), int((d > 12).sum())))
print("6 runs, vpot_points_bohr.xyz and sites.tsv written to", OUT)
