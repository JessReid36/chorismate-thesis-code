#!/usr/bin/env python3
"""
s1_signcheck.py - sign-convention check for the difference potential Dv = V_TS - V_R.

Two modes, called in sequence by s1_signcheck.pbs:

  points <geom_dir> <provenance.tsv> <out_bohr.xyz> <out_labels.tsv>
      Validate the reactant and TS geometries against the committed
      provenance table, confirm they share a coordinate frame, and write the
      two probe points in Bohr for orca_vpot.

  assemble <labels.tsv> <vR.out> <vTS.out> <out_table.tsv>
      Read the two potentials, form Dv = V_TS - V_R, and decide PASS / FAIL.

Probe points
------------
  A  2.0 A beyond the ether oxygen O3, along C4 -> O3.
     Arg90 occupies this region in BsCM. Negative charge develops on O3 as
     the C4-O3 bond breaks, so a +1 cation here must be stabilising: Dv < 0.

  B  2.0 A beyond C6, along C1 -> C6.
     Charge is depleted in the forming-bond region, so the sign must invert:
     Dv > 0. This is the negative control. Without it, a Dv that was
     uniformly negative everywhere - a global offset carrying no
     reaction-coordinate information - would still satisfy point A.

Pass criterion: Dv(A) < 0 AND Dv(B) > 0. Both are required.

Scope
-----
Both probes sit inside the candidate-charge envelope (roughly 0.1-0.3 A from
the substrate vdW surface, against a grid minimum standoff near 2.0 A). They
are placed close in deliberately, to maximise signal for a yes/no convention
test. The MAGNITUDES here are therefore diagnostic only and must not be
reported as stabilisation energies; the physically meaningful quantities are
the Dv values on the candidate grid, where every site is a position a charge
could actually occupy.

Nothing in this script is hardcoded to a particular run. Expected bond lengths
are read from the committed provenance table, so the guards track the
geometry of record rather than a value typed in here.
"""

import math
import sys
from pathlib import Path

ANG2BOHR = 1.8897259886
HARTREE2KCAL = 627.5094740631

# 0-based indices, per 00_admin/reacting_atoms.tsv
IDX = {"C1": 0, "O3": 7, "C4": 8, "C6": 12}

# Guard tolerances
BOND_TOL = 0.02        # A, on each reacting bond vs the provenance table
CENTROID_TOL = 0.50    # A, max centroid separation for a shared frame
ROTATION_TOL = 0.10    # A, max RMSD improvement from rigid-body alignment


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    return [tuple(float(x) for x in line.split()[1:4]) for line in lines[2:2 + n]]


def read_provenance(path):
    """Expected reacting-bond lengths, from the committed provenance table."""
    rows = Path(path).read_text().splitlines()
    header = rows[0].split("\t")
    try:
        li, bi, fi = (header.index("label"),
                      header.index("break_expect"),
                      header.index("form_expect"))
    except ValueError:
        sys.exit(f"FAIL: {path} lacks label/break_expect/form_expect columns")
    out = {}
    for row in rows[1:]:
        if not row.strip():
            continue
        p = row.split("\t")
        out[p[li]] = (float(p[bi]), float(p[fi]))
    return out


def sub(a, b):
    return tuple(a[i] - b[i] for i in range(3))


def add(a, b):
    return tuple(a[i] + b[i] for i in range(3))


def scale(a, s):
    return tuple(x * s for x in a)


def norm(a):
    return math.sqrt(sum(x * x for x in a))


def centroid(pts):
    n = len(pts)
    return tuple(sum(p[i] for p in pts) / n for i in range(3))


def bonds(xyz):
    return (norm(sub(xyz[IDX["O3"]], xyz[IDX["C4"]])),
            norm(sub(xyz[IDX["C1"]], xyz[IDX["C6"]])))


def rmsd_as_stored(P, Q):
    return math.sqrt(sum(sum((P[i][k] - Q[i][k]) ** 2 for k in range(3))
                         for i in range(len(P))) / len(P))


def rmsd_aligned(P, Q):
    """
    Kabsch RMSD after optimal rigid-body superposition, via SVD.

    Returns None if numpy is unavailable, in which case the caller falls back
    to the centroid test alone. (The PBS wrapper sets OPENBLAS_NUM_THREADS=1
    and OMP_NUM_THREADS=1 before invoking this, which is what numpy needs to
    be well behaved here.)
    """
    try:
        import numpy as np
    except ImportError:
        return None

    A = np.asarray(P, float)
    B = np.asarray(Q, float)
    A = A - A.mean(0)
    B = B - B.mean(0)

    V, _, Wt = np.linalg.svd(A.T @ B)
    d = np.sign(np.linalg.det(V @ Wt))          # guard against a reflection
    R = V @ np.diag([1.0, 1.0, d]) @ Wt

    return float(np.sqrt(((A @ R - B) ** 2).sum(1).mean()))


def check_geometry(label, xyz, expect):
    """Confirm a geometry is the one of record."""
    if len(xyz) != 24:
        sys.exit(f"FAIL: expected 24 atoms in the {label}, found {len(xyz)}")
    d_break, d_form = bonds(xyz)
    e_break, e_form = expect
    print(f"  {label:9s} break O3-C4 = {d_break:.3f} A (expect {e_break:.3f})   "
          f"form C1-C6 = {d_form:.3f} A (expect {e_form:.3f})")
    if abs(d_break - e_break) > BOND_TOL or abs(d_form - e_form) > BOND_TOL:
        sys.exit(f"FAIL: {label} does not match the provenance table "
                 f"(tolerance {BOND_TOL} A) - wrong file?")


def check_common_frame(R, T):
    """
    Dv is evaluated at points fixed in the laboratory frame, so V_R and V_TS
    are only comparable if the two geometries are expressed in the same frame.
    A relative translation or rotation would make the difference meaningless
    while still producing plausible-looking numbers.
    """
    sep = norm(sub(centroid(R), centroid(T)))
    stored = rmsd_as_stored(R, T)
    aligned = rmsd_aligned(R, T)

    print(f"  centroid separation R vs TS      = {sep:.3f} A  (tolerance {CENTROID_TOL})")
    if sep > CENTROID_TOL:
        sys.exit("FAIL: reactant and TS centroids are too far apart to share a frame.")

    if aligned is None:
        print("  WARNING: numpy unavailable, rotation check SKIPPED (centroid test only)")
    else:
        gain = stored - aligned
        print(f"  RMSD as stored / rigidly aligned = {stored:.3f} / {aligned:.3f} A")
        print(f"  improvement from alignment       = {gain:.3f} A  (tolerance {ROTATION_TOL})")
        if gain > ROTATION_TOL:
            sys.exit("FAIL: rigid-body alignment substantially reduces the RMSD, so the two\n"
                     "      geometries are rotated relative to one another. Re-express them\n"
                     "      in a common frame before evaluating Dv at fixed points.")
    print("  frame check OK - the two densities are directly comparable at fixed points")


def build_points(geom_dir, provenance, out_bohr, out_labels):
    geom_dir = Path(geom_dir)
    expect = read_provenance(provenance)
    R = read_xyz(geom_dir / "reactant.xyz")
    T = read_xyz(geom_dir / "ts.xyz")

    print("=== geometry guards (expectations read from the provenance table) ===")
    check_geometry("reactant", R, expect["reactant"])
    check_geometry("TS", T, expect["ts"])
    print()
    print("=== common-frame check ===")
    check_common_frame(R, T)
    print()

    O3, C4 = R[IDX["O3"]], R[IDX["C4"]]
    C1, C6 = R[IDX["C1"]], R[IDX["C6"]]

    u = scale(sub(O3, C4), 1.0 / norm(sub(O3, C4)))
    pA = add(O3, scale(u, 2.0))

    u2 = scale(sub(C6, C1), 1.0 / norm(sub(C6, C1)))
    pB = add(C6, scale(u2, 2.0))

    pts = [("A_etherO_breaking", pA), ("B_formingC6", pB)]

    with open(out_bohr, "w") as fh:
        fh.write(f"{len(pts)}\n")
        for _, p in pts:
            fh.write("  ".join(f"{c * ANG2BOHR:.10f}" for c in p) + "\n")

    with open(out_labels, "w") as fh:
        fh.write("label\tx_ang\ty_ang\tz_ang\n")
        for name, p in pts:
            fh.write(f"{name}\t{p[0]:.4f}\t{p[1]:.4f}\t{p[2]:.4f}\n")

    print(f"wrote {len(pts)} probe points to {out_bohr} (Bohr) and {out_labels} (Angstrom)")


def read_vpot(path):
    vals = []
    for line in Path(path).read_text().splitlines():
        p = line.split()
        if len(p) >= 4:
            try:
                vals.append(float(p[-1]))
            except ValueError:
                pass
    return vals


def assemble(labels_path, vR_path, vTS_path, out_table):
    labels = [l.split("\t")[0]
              for l in Path(labels_path).read_text().splitlines()[1:]]
    vR, vTS = read_vpot(vR_path), read_vpot(vTS_path)

    if not (len(vR) == len(vTS) == len(labels)):
        sys.exit(f"FAIL: length mismatch - vR={len(vR)} vTS={len(vTS)} labels={len(labels)}")

    dv = {}
    print(f"{'probe':22s} {'V_R (Eh)':>12s} {'V_TS (Eh)':>12s} "
          f"{'Dv (Eh)':>12s} {'Dv (kcal/mol per +1e)':>24s}")
    print("-" * 88)
    rows = []
    for lab, a, b in zip(labels, vR, vTS):
        d = b - a
        dv[lab[0]] = d
        print(f"{lab:22s} {a:12.6f} {b:12.6f} {d:+12.6f} "
              f"{d * HARTREE2KCAL:+18.3f}   "
              f"{'STABILISING' if d < 0 else 'destabilising'}")
        rows.append((lab, a, b, d))

    with open(out_table, "w") as fh:
        fh.write("label\tV_R_Eh\tV_TS_Eh\tDv_Eh\tDv_kcal_per_e\n")
        for lab, a, b, d in rows:
            fh.write(f"{lab}\t{a:.6f}\t{b:.6f}\t{d:.6f}\t{d * HARTREE2KCAL:.4f}\n")

    for key in ("A", "B"):
        if key not in dv:
            sys.exit(f"FAIL: probe {key} missing from the results")

    ok_A, ok_B = dv["A"] < 0, dv["B"] > 0
    print()
    print(f"  A (ether O, breaking bond) Dv < 0 : {'PASS' if ok_A else 'FAIL'}")
    print(f"  B (forming bond)           Dv > 0 : {'PASS' if ok_B else 'FAIL'}")
    print()

    if ok_A and ok_B:
        print("SIGN CHECK PASS.")
        print("A +1 cation beyond the ether oxygen lowers the barrier, and the sign")
        print("inverts in the forming-bond region. Dv therefore carries genuine")
        print("reaction-coordinate information rather than a uniform offset, and the")
        print("convention used by the charge-selection objective is correct.")
        print()
        print("Magnitudes above are diagnostic only: both probes sit inside the")
        print("candidate-charge envelope and must not be quoted as stabilisation")
        print("energies. Report signs here; report magnitudes from the grid Dv map.")
        return 0

    print("SIGN CHECK FAIL - do not proceed to charge selection.")
    print("Investigate the subtraction order and the units before anything downstream.")
    return 1


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    mode = sys.argv[1]
    if mode == "points":
        build_points(*sys.argv[2:6])
    elif mode == "assemble":
        sys.exit(assemble(*sys.argv[2:6]))
    else:
        sys.exit(f"unknown mode: {mode}")


if __name__ == "__main__":
    main()
