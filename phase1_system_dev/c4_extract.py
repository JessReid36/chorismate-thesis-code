#!/usr/bin/env python3
"""
c4_extract.py RUNDIR - C4 in-job extraction (C4_CRITERIA.txt). Called at the end of each run's reopt.pbs on hpc1, after
ORCA. Reads RUNDIR/reopt.xyz (ORCA's last geometry, full system), RUNDIR/start.xyz and RUNDIR/../active_atoms.txt;
checks that every atom outside the active region is where start.xyz put it (the active-region displacement A3 is taken
without superposition, so the fixed atoms must not have moved); writes RUNDIR/final_active.xyz: the active atoms in
active_atoms.txt's order, coordinates copied verbatim, 0-based atom index in a 5th column. The full-system file stays
on hpc1 (the results commit records its sha256); the analysis reads final_active.xyz.
Python 3.6, standard library; deterministic. Exit 1 if reopt.xyz is missing, not the full system, or a fixed atom moved.
"""
import glob, os, sys
D = sys.argv[1]
NATOMS = 55680
def rows(p):
    L = open(p).read().splitlines()
    n = int(L[0].split()[0]); R = [l.split() for l in L[2:2 + n]]
    if n != NATOMS or len(R) != NATOMS or any(len(r) < 4 for r in R):
        sys.exit("c4_extract: %s holds %d atoms, expected %d" % (p, n, NATOMS))
    return R
p = os.path.join(D, "reopt.xyz")
if not os.path.isfile(p):
    sys.exit("c4_extract: no reopt.xyz in %s (xyz files present: %s)" % (D, " ".join(sorted(os.path.basename(x) for x in glob.glob(os.path.join(D, "*.xyz"))))))
R = rows(p); S = rows(os.path.join(D, "start.xyz"))
ACT = [int(x) for x in open(os.path.join(D, "..", "active_atoms.txt")).read().split()]
AS = set(ACT)
moved = max(abs(float(R[i][k]) - float(S[i][k])) for i in range(NATOMS) if i not in AS for k in (1, 2, 3))
if moved > 2e-6: sys.exit("c4_extract: an atom outside the active region moved by %.1e A" % moved)
open(os.path.join(D, "final_active.xyz"), "w").write(
    "%d\n%s: last geometry (reopt.xyz), active atoms; 5th column: 0-based atom index\n" % (len(ACT), os.path.basename(os.path.abspath(D)))
    + "".join("%s %s %s %s %d\n" % (R[i][0], R[i][1], R[i][2], R[i][3], i) for i in ACT))
