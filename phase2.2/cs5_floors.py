#!/usr/bin/env python3
"""
cs5_floors.py - Stage 5 of the charged-sphere development: how close a sphere of charge q may sit to each substrate
atom (its FLOOR), by charge level. The grid (cs5_grid.py) places sites on these floors and gives every site the charges
that may sit there, so that no sphere is placed nearer an atom than any resting sphere of its charge was seen to come
(rule B, Jessica, 19 October 2026). STAGE2_PROTOCOL.txt: "where a sign of charge has no contact minimum, the floor for
it comes from the Lennard-Jones pair minimum and Stage 3".

Floor of atom J for a sphere of charge q (Amber N sphere, R* 1.824 A, eps 0.170 kcal/mol):

  q >= +0.25   per tabulated level q_k = +0.25, +0.50, +0.75, +1.00:
                 floor_J(q_k) = min( LJ_J , closest approach to J of a Stage 2 minimum at q_k )
               over the minima of that level (four oxygen starts x 60 geometries): with minima=all (the default),
               every minimum, own basin and migrated, since each is a place a sphere of that charge comes to rest
               (gradient below 1e-5 kcal/mol/A); with minima=own, own-basin minima only. A floor is one atom's closest
               approach, applied around that atom in every direction. A charge between two levels takes the floors of
               the lower level (farther out). The levels are made monotone (a larger charge never has a larger floor):
               floor_J(q_k) := max over levels >= q_k, so a site cleared for one level is cleared for every higher one.
               H8 (GAFF ho) has no LJ parameters, so for positive spheres its floor is the closest approach of a sphere
               of that level; these spheres rest at oxygens (O4 and the carboxylates), so the H8 floors (3.68-4.02 A)
               are distances observed from there, not a limit set by H8.
  0 < q < 0.25 no Stage 2 data; the LJ pair minima, with H8 at the q = +0.25 level (farther out).
  q <= 0       negative spheres have no contact minimum at any q (Stage 2: all 1920 starts escaped), so the floor is
               the LJ pair minimum, where a sphere held near the substrate meets its repulsive wall; for H8 (no LJ)
               it is Stage 3's smallest H8-sphere distance, 2.331 A (n1_41786_qm1: a -1 sphere in QM/MM, no artefact).
  LJ_J         R_min,ij = (R_min,sphere + R_min,J)/2 on the full R_min, ORCA's rule (Stage 1 S1.4).

Stage 2's minima are recomputed here with Stage 2's own code (cs2_contact_scan.py's starts and options, WP5's energy
model and minimiser, imported unchanged), because cs2_contacts.tsv kept distances, not positions. Every recomputed
minimum must reproduce its cs2_contacts.tsv row: the same outcome (own basin or migrated), contact atom and nearest
atom, and the contact, nearest and H8 distances within 1e-3 A. That is Stage 2's reproducibility criterion
(STAGE2_REPRODUCIBILITY.txt, x46) extended from own-basin to migrated minima, because migrated minima set floors here.
Otherwise the script stops. The header records the largest differences and the versions.

USAGE  python3 cs5_floors.py <results repo> <out floors.tsv> [minima=all|own]   (needs numpy and scipy, as Stage 2)
"""
import hashlib, importlib.util, sys
import scipy
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, HERE / fn)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
wp5 = load("wp5", "wp5_rigid_lj_scan.py")
cs2 = load("cs2", "cs2_contact_scan.py")
NAMES = wp5.NAMES
H8 = NAMES.index("H8")
STAGE3_H8 = 2.331          # cs_stage3/STAGE3_REPORT.txt, n1_41786_qm1: H8..sphere smallest along the path 2.331 A
LEVELS = [0.25, 0.5, 0.75, 1.0]
COLS = ["neg", "pos0", "p025", "p050", "p075", "p100"]
MINIMA = "all"

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    global MINIMA
    if len(sys.argv) not in (3, 4) or (len(sys.argv) == 4 and sys.argv[3] not in ("minima=all", "minima=own")): sys.exit(__doc__)
    MINIMA = sys.argv[3].split("=")[1] if len(sys.argv) == 4 else "all"
    R = Path(sys.argv[1]); out = Path(sys.argv[2])
    if out.exists(): sys.exit("STOP: %s exists" % out)
    el, q, eps, rc = wp5.load_substrate(R)
    ref = {}
    tsv = R / "phase2.2/cs_stage2/cs2_contacts.tsv"
    for l in tsv.read_text().splitlines()[1:]:
        t = l.split("\t")
        if t[2] == "N" and t[3] in ("+0.25", "+0.50", "+0.75", "+1.00"): ref[(t[0], t[1], float(t[3]), t[4])] = t
    frames = sorted(set(k[0] for k in ref))
    if len(frames) != 30 or len(ref) != 960: sys.exit("STOP: cs2_contacts.tsv does not hold the 960 N-sphere positive rows")
    rs, ep = 1.824, 0.170
    closest = {k: np.full(24, np.inf) for k in LEVELS}; where = {k: [""] * 24 for k in LEVELS}; n = 0
    dmax = {"own": 0.0, "migrated": 0.0}
    for fr in frames:
        for st in ("R", "TS"):
            _, X = wp5.read_xyz(R / "05_qmmm/20_invacuo/{}_{}.xyz".format(fr, st))
            for qk in LEVELS:
                fE = wp5.make_energy(X, q, eps, rc, 2 * rs, ep, q_s=qk)
                for gname, members, anchor in wp5.GROUPS:
                    r = minimize(lambda p: fE(p), cs2.o_start(X, members, anchor), jac=True, method="BFGS", options=cs2.OPT)
                    _, g = fE(r.x)
                    dd = np.linalg.norm(X - r.x, axis=1)
                    jn = min(wp5.O_IDX, key=lambda i: dd[i]); ja = int(dd.argmin())
                    t = ref[(fr, st, qk, gname)]
                    own = "own" if NAMES[jn] in members else "migrated"
                    diff = max(abs(dd[jn] - float(t[7])), abs(dd[ja] - float(t[9])), abs(dd[H8] - float(t[10])))
                    if (np.linalg.norm(g) > 1e-5 or own != t[5] or NAMES[jn] != t[6] or NAMES[ja] != t[8] or diff > 1e-3):
                        sys.exit("STOP: {} {} q={} {}: the recomputed minimum does not reproduce cs2_contacts.tsv".format(fr, st, qk, gname))
                    dmax[own] = max(dmax[own], diff)
                    n += 1
                    if MINIMA == "own" and own != "own": continue
                    for j in range(24):
                        if dd[j] < closest[qk][j]:
                            closest[qk][j] = dd[j]; where[qk][j] = "{} {} {} start ({})".format(fr, st, gname, own)
    lj = 0.5 * (2 * rs + rc)
    F = {}
    for qk in LEVELS:
        f = np.minimum(lj, closest[qk]); f[H8] = closest[qk][H8]
        F[qk] = f
    raw = {qk: F[qk].copy() for qk in LEVELS}
    for i in range(len(LEVELS) - 2, -1, -1):              # monotone: a smaller charge never sits closer than a larger one
        F[LEVELS[i]] = np.maximum(F[LEVELS[i]], F[LEVELS[i + 1]])
    neg = lj.copy(); neg[H8] = STAGE3_H8
    pos0 = lj.copy(); pos0[H8] = F[0.25][H8]
    cols = {"neg": neg, "pos0": pos0, "p025": F[0.25], "p050": F[0.5], "p075": F[0.75], "p100": F[1.0]}
    mono = sum(int(np.sum(F[qk] > raw[qk] + 1e-12)) for qk in LEVELS)
    hdr = ["# grid floors by charge level, cs5_floors.py (Stage 5; sha256 %s): A, Amber N sphere (R* %.3f A, eps %.3f kcal/mol)" % (sha(__file__), rs, ep),
           "#   neg   q <= 0: LJ pair minimum; H8 Stage 3 (%.3f A)" % STAGE3_H8,
           "#   pos0  0 < q < 0.25: LJ pair minimum; H8 at the +0.25 level",
           "#   pNNN  q >= level: min(LJ pair minimum, closest approach of a Stage 2 minimum at that level); H8 its closest approach;",
           "#         minima=%s (%s)" % (MINIMA, "own basin and migrated" if MINIMA == "all" else "own basin only"),
           "#         monotone in q (%d entries raised to the next level's floor)" % mono,
           "# Stage 2 minima recomputed and checked against cs2_contacts.tsv (sha256 %s): %d of 960 reproduced (same outcome" % (sha(tsv), n),
           "#   and atoms; distances within 1e-3 A); largest distance difference: own basin %.1e A, migrated %.1e A" % (
               dmax["own"], dmax["migrated"]),
           "# python %s, numpy %s, scipy %s" % (sys.version.split()[0], np.__version__, scipy.__version__),
           "# force field 05_qmmm/13_bridge/complex_solvated.ORCAFF.prms sha256 %s" % sha(R / "05_qmmm/13_bridge/complex_solvated.ORCAFF.prms"),
           "atom\t" + "\t".join(COLS) + "\tLJ_pair_min\tclosest_p100_at"]
    rows = ["%s\t%s\t%.3f\t%s" % (NAMES[j], "\t".join("%.3f" % cols[c][j] for c in COLS), lj[j], where[1.0][j]) for j in range(24)]
    out.write_text("\n".join(hdr + rows) + "\n")
    print("\n".join(hdr + rows))

if __name__ == "__main__":
    main()
