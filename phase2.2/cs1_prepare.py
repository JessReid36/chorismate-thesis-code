#!/usr/bin/env python3
"""
cs1_prepare.py - Stage 1 of the charged-sphere development (implementation checks in ORCA).
Builds seven QM/MM single points (energy + gradient) on frame 41786 from committed inputs only:
  substrate reactant geometry and force-field rows: phase2.2/wp4_s13_rerun/armb_anch/{combined.xyz,site.ORCAFF.prms}
  substrate TS geometry: 05_qmmm/19_ensemble/frame_41786/neb_NEB-CI_converged.QMRegion.xyz (same frame and atom order)
Spheres: Amber N Lennard-Jones (R_min 3.648 A full, eps 0.170 kcal/mol), as WP4 arm b.
  sphere 1 = WP4 arm b start (40.832, 18.750, 42.551): 3.40 A from the nearest substrate atom at R, 3.18 at TS
  sphere 2 = 3.5 A outward from O1 along centroid->O1 at R (44.525, 25.580, 44.863): 3.5 A at R, 3.6 at TS;
             8.10 A from sphere 1, inside both cutoffs (LJ switch 10-12 A, Coulomb 12 A), so the
             sphere-sphere terms are non-zero
Runs: s1a_q1_N, s1b_q1_C, s1c_q1_Ar (element label only differs); s1d_q0_N, s1e_q05_N (charge only differs);
      s1f_two_R, s1g_two_TS (two spheres, +0.5 and -0.5, frozen; substrate at R vs TS).
USAGE  python3 cs1_prepare.py <results repo> <output folder>
"""
import os, subprocess, sys
REPO, OUT = sys.argv[1], sys.argv[2]
W = os.path.join(REPO, "phase2.2/wp4_s13_rerun/armb_anch")
xyzR = open(os.path.join(W, "combined.xyz")).read().splitlines()
R = [l.split() for l in xyzR[2:26]]
TS = [l.split() for l in subprocess.run(["git", "-C", REPO, "show", "HEAD:05_qmmm/19_ensemble/frame_41786/neb_NEB-CI_converged.QMRegion.xyz"],
                                        capture_output=True, text=True, check=True).stdout.splitlines()[2:26]]
assert len(R) == 24 and len(TS) == 24 and [a[0] for a in R] == [a[0] for a in TS]
P = open(os.path.join(W, "site.ORCAFF.prms")).read().splitlines()
i = P.index("$atoms")
SUB = P[i + 2:i + 2 + 24]
assert len(SUB) == 24 and SUB[0].split()[1] == "C" and SUB[23].split()[1] == "O"
S1 = (40.832, 18.750, 42.551)
S2 = (44.525, 25.580, 44.863)
INP = """! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF QMMM EnGrad
%maxcore 3000
%pal nprocs 4 end
%scf MaxIter 250 end
%qmmm
  QMAtoms {0:23} end
  ORCAFFFilename "site.ORCAFF.prms"
end
* xyzfile -2 1 combined.xyz
"""
def row(k, el, q):
    return "%6d   %-2s%13.6f%13.6f%13.6f%13.6f%13.6f" % (k, el, q, -0.17, 3.648, -0.085, 3.648)
def write(name, sub_xyz, spheres, title):
    d = os.path.join(OUT, name); os.makedirs(d)
    n = 24 + len(spheres)
    with open(os.path.join(d, "combined.xyz"), "w") as f:
        f.write("%d\n%s\n" % (n, title))
        for a in sub_xyz: f.write("%-2s %15.8f%15.8f%15.8f\n" % (a[0], float(a[1]), float(a[2]), float(a[3])))
        for el, q, p in spheres: f.write("%-2s %15.8f%15.8f%15.8f\n" % (el, p[0], p[1], p[2]))
    with open(os.path.join(d, "site.ORCAFF.prms"), "w") as f:
        f.write("$fftype\nAMBER\n$atoms\n%d 1 4\n" % n)
        for l in SUB: f.write(l + "\n")
        for k, (el, q, p) in enumerate(spheres, 25): f.write(row(k, el, q) + "\n")
        f.write("$bonds\n0 2 2\n$angles\n0 3 2\n$dihedrals\n0 4 3\n")
    open(os.path.join(d, "job.inp"), "w").write(INP)
write("s1a_q1_N",  R, [("N", 1.0, S1)],  "frame 41786 R + sphere q=+1 label N")
write("s1b_q1_C",  R, [("C", 1.0, S1)],  "frame 41786 R + sphere q=+1 label C")
write("s1c_q1_Ar", R, [("Ar", 1.0, S1)], "frame 41786 R + sphere q=+1 label Ar")
write("s1d_q0_N",  R, [("N", 0.0, S1)],  "frame 41786 R + sphere q=0 label N")
write("s1e_q05_N", R, [("N", 0.5, S1)],  "frame 41786 R + sphere q=+0.5 label N")
write("s1f_two_R", R,  [("N", 0.5, S1), ("N", -0.5, S2)], "frame 41786 R + spheres +0.5 and -0.5")
write("s1g_two_TS", TS, [("N", 0.5, S1), ("N", -0.5, S2)], "frame 41786 TS + spheres +0.5 and -0.5")
print("7 runs written to", OUT)
