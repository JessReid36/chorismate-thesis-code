#!/usr/bin/env python3
"""
s16_align_frames.py - superimpose the ensemble onto a common reference frame.

WHY THIS IS NEEDED, AND WHAT IT FIXES
The grid and the difference-potential matrix are built in a FIXED Cartesian frame. The
ensemble geometries come straight out of molecular dynamics, so each sits wherever the
substrate happened to be in the simulation box. Measured on the 30 post-cut frames:

    centroid spread          x 10.50   y 13.59   z 25.06  Angstrom
    centroid-to-centroid     min 1.16, mean 10.79, max 25.25 Angstrom

The substrate is about 6 Angstrom across, so the frames are separated by more than the
molecule's own size. A fixed grid point is therefore in a COMPLETELY DIFFERENT PLACE
relative to each frame's substrate, and an a-matrix built that way is not thirty
measurements of one site but thirty measurements of thirty unrelated positions.

The symptom, from the first unaligned run: only 4 of 580 sites had all frames agreeing on
the SIGN of the difference potential (0.7 per cent), and the per-frame stabilising
fraction ranged from 19 to 88 per cent. Those look like findings about conformational
diversity. They are not. They are a coordinate-frame error.

[DTHESIS] Ch. 7.3 states the requirement plainly: the surface is built over ALIGNED path
frames.

WHAT IS ALIGNED ON, AND WHY
Measured per-atom deviation after a whole-molecule fit across the post-cut frames:

    ring carbons and the ether oxygen      0.05 - 0.11 A    rigid
    carboxylate oxygens (5, 6, 22, 23)     0.29 - 0.47 A    freely rotating
    hydroxyl and methylene hydrogens       0.19 - 0.27 A

Fitting on all 24 atoms lets the rotating carboxylates drag the superposition, which
moves the reacting centre for no good reason. The default alignment set is therefore the
RIGID CORE: the ring carbons, the ether oxygen and the atoms directly bonded into the
ring. Pass align_set=all to fit on everything instead and compare.

THE TRANSFORM IS DERIVED ONCE PER FRAME AND APPLIED TO EVERY STATE
The rotation and translation are computed from the REACTANT geometry and then applied
unchanged to that frame's transition state and product. Aligning each state separately
would superimpose the transition states on each other and destroy the reactant-to-TS
motion, which is the entire quantity being designed against.

USAGE
    python3 s16_align_frames.py <ensemble_dir> <out_dir> [key=value ...]
    python3 s16_align_frames.py 05_qmmm/19_ensemble_barriers phase2.2/aligned \\
        frames=20000,21634,... align_set=core
"""
import sys
from pathlib import Path

import numpy as np

CONFIG = {
    "align_set": "core",     # core | all | heavy
    "reference": "",         # frame to align onto; empty = the first requested
    "frames": "",            # comma-separated; empty = every frame found
    "states": "R,TS,P",
}

# 0-based indices into the 24-atom QM region, in the committed atom order:
#   0:C 1:H 2:H 3:C 4:C 5:O 6:O 7:O 8:C 9:H 10:C 11:H 12:C 13:C 14:H 15:C 16:H
#   17:C 18:H 19:O 20:H 21:C 22:O 23:O
# Index 7 is the ether oxygen O3, the atom the design targets.
CORE = [3, 7, 8, 10, 12, 13, 15, 17]        # ring carbons plus the ether oxygen
CARBOXYLATE_O = [5, 6, 22, 23]              # free to rotate; excluded from the fit

SRC_STATE = {"R": "reactant_qm.xyz",
             "TS": "transition_state_qm.xyz",
             "P": "product_qm.xyz"}
ALT_STATE = {"R": "_R.xyz", "TS": "_TS.xyz", "P": "_P.xyz"}


def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    els = [l.split()[0] for l in L[2:2 + n]]
    xyz = np.array([[float(x) for x in l.split()[1:4]] for l in L[2:2 + n]])
    return els, xyz


def write_xyz(p, els, xyz, comment):
    with open(p, "w") as fh:
        fh.write(f"{len(els)}\n{comment}\n")
        for e, (x, y, z) in zip(els, xyz):
            fh.write(f"{e:<3}{x:15.8f}{y:15.8f}{z:15.8f}\n")


def kabsch(P, Q, idx):
    """Rotation and translation taking P onto Q, fitted on the atoms in idx.

    Returns (R, t) such that P @ R + t is the aligned structure. The reflection guard
    on det(V @ W) is essential: without it the fit can return an improper rotation,
    which mirrors the molecule and silently inverts its chirality.
    """
    Pc, Qc = P[idx].mean(0), Q[idx].mean(0)
    A = (P[idx] - Pc).T @ (Q[idx] - Qc)
    V, S, W = np.linalg.svd(A)
    d = np.sign(np.linalg.det(V @ W))
    R = V @ np.diag([1.0, 1.0, d]) @ W
    return R, Qc - Pc @ R


def find_geom(ens, alt, fr, state):
    f = ens / f"frame_{fr}" / SRC_STATE[state]
    if f.exists():
        return f
    f = alt / f"{fr}{ALT_STATE[state]}"
    return f if f.exists() else None


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    ens, out = Path(sys.argv[1]), Path(sys.argv[2])
    cfg = dict(CONFIG)
    for kv in sys.argv[3:]:
        k, v = kv.split("=", 1)
        if k not in cfg:
            sys.exit(f"unknown parameter {k!r}. Known: {', '.join(sorted(cfg))}")
        cfg[k] = v
    alt = ens.parent / "20_invacuo"
    states = [s for s in cfg["states"].split(",") if s]
    frames = [f for f in cfg["frames"].split(",") if f]
    if not frames:
        frames = sorted(d.name.replace("frame_", "") for d in ens.glob("frame_*"))

    # the reference frame must have a reactant geometry
    ref_name = cfg["reference"] or frames[0]
    rf = find_geom(ens, alt, ref_name, "R")
    if rf is None:
        sys.exit(f"reference frame {ref_name} has no reactant geometry")
    ref_els, ref_xyz = read_xyz(rf)
    n = len(ref_els)

    if cfg["align_set"] == "core":
        idx = CORE
        what = "rigid core: ring carbons and the ether oxygen"
    elif cfg["align_set"] == "heavy":
        idx = [i for i, e in enumerate(ref_els) if e.upper() != "H"]
        what = "all heavy atoms"
    elif cfg["align_set"] == "all":
        idx = list(range(n))
        what = "all atoms"
    else:
        sys.exit("align_set must be core, heavy or all")

    out.mkdir(parents=True, exist_ok=True)
    print(f"reference: frame {ref_name}, reactant")
    print(f"aligning on {len(idx)} atoms ({what})")
    print(f"excluded from the fit: carboxylate oxygens {CARBOXYLATE_O} rotate freely "
          f"and would drag the superposition" if cfg["align_set"] == "core" else "")

    rows, done, skipped = [], 0, []
    for fr in frames:
        f_r = find_geom(ens, alt, fr, "R")
        if f_r is None:
            skipped.append(fr)
            continue
        els, xyz = read_xyz(f_r)
        if len(els) != n or els != ref_els:
            skipped.append(f"{fr}(atom order differs)")
            continue
        # ONE transform per frame, from the reactant, applied to every state
        R, t = kabsch(xyz, ref_xyz, idx)
        d = out / f"frame_{fr}"
        d.mkdir(exist_ok=True)
        # The transform is SAVED because the wavefunctions were computed on the
        # ORIGINAL geometries and still live in those coordinates. To evaluate a
        # frame's potential at a point p of the common aligned grid, that point must be
        # mapped BACK into the frame's own frame of reference:
        #     p_orig = (p - t) @ R^T
        # s15_dv_matrix.py reads this file and does exactly that, which is why no SCF
        # has to be repeated.
        np.savetxt(d / "transform.txt", np.vstack([R, t.reshape(1, 3)]),
                   fmt="%.12f",
                   header="rows 0-2: rotation R; row 3: translation t\n"
                          "aligned = original @ R + t ; original = (aligned - t) @ R.T")
        pre = np.linalg.norm(xyz[idx] - ref_xyz[idx], axis=1)
        post_xyz = xyz @ R + t
        post = np.linalg.norm(post_xyz[idx] - ref_xyz[idx], axis=1)
        rows.append((fr, float(np.sqrt((pre ** 2).mean())),
                     float(np.sqrt((post ** 2).mean()))))
        for s in states:
            fs = find_geom(ens, alt, fr, s)
            if fs is None:
                continue
            e2, x2 = read_xyz(fs)
            write_xyz(d / SRC_STATE[s], e2, x2 @ R + t,
                      f"frame {fr} {s}, aligned onto {ref_name} on {len(idx)} atoms")
        done += 1

    print(f"\naligned {done} frames, {len(states)} states each")
    if skipped:
        print(f"  SKIPPED {len(skipped)}: {' '.join(skipped)}")
    print(f"\n  {'frame':>8}{'RMSD before':>14}{'RMSD after':>13}")
    for fr, a, b in rows:
        print(f"  {fr:>8}{a:>13.2f} A{b:>12.3f} A")
    after = [b for _, _, b in rows]
    print(f"\n  RMSD on the fitted atoms after alignment: "
          f"min {min(after):.3f}, mean {sum(after)/len(after):.3f}, "
          f"max {max(after):.3f} A")

    cent = []
    for fr, _, _ in rows:
        _, x = read_xyz(out / f"frame_{fr}" / SRC_STATE["R"])
        cent.append(x.mean(0))
    C = np.array(cent)
    print(f"  centroid spread now: x {np.ptp(C[:,0]):.2f}  y {np.ptp(C[:,1]):.2f}  "
          f"z {np.ptp(C[:,2]):.2f} A   (was 10.50 / 13.59 / 25.06)")

    (out / "ALIGNMENT.txt").write_text(
        f"reference: frame {ref_name}, reactant geometry\n"
        f"alignment set: {cfg['align_set']} ({len(idx)} atoms) - {what}\n"
        f"atom indices: {idx}\n"
        f"frames aligned: {done}\n"
        f"states: {','.join(states)}\n"
        "\nThe transform was derived ONCE per frame from the reactant and applied\n"
        "unchanged to that frame's other states, so the reactant-to-transition-state\n"
        "motion is preserved. Aligning each state separately would superimpose the\n"
        "transition states on one another and destroy the quantity being designed for.\n"
        f"\nRMSD on the fitted atoms after alignment: min {min(after):.3f}, "
        f"mean {sum(after)/len(after):.3f}, max {max(after):.3f} A\n")
    print(f"\nwrote {out}/ and {out}/ALIGNMENT.txt")
    print("\nnow rebuild the grid and the a-matrix on THIS directory, not the original.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
