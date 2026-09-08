#!/usr/bin/env python3
"""
step_2_5_grid_v2.py - candidate-position grid, second generation.

What changed from the v1 pipeline (step_2_1 .. step_2_4) and why
----------------------------------------------------------------
v1 used a single parameter, r_min = 1.5 A, to do two unrelated jobs: it set how
finely a charge could be positioned, AND it set how close two placed charges
could be. Those are different requirements and are now separated.

  RESOLUTION      is a property of the grid. Finer is better, limited only by
                  problem size. Set here by --r-min (default 1.0 A, ~1100 pts).

  SEPARATION      is a property of the physical species being placed - about
                  3.5 A centre-to-centre for formate, 4.5 A for guanidinium -
                  and is now enforced as a constraint in the optimiser, not
                  baked into the grid. This grid does not decide it.

Measured cost of the v1 conflation: an arbitrary position on the shell sat on
average 0.909 A from the nearest available site (max 1.800 A). Against the
measured Dv gradient (mean 0.279, p95 0.888 kcal/mol/A per +1e) that is a
typical placement penalty of 0.25 kcal/mol and up to ~1.6 kcal/mol in the steep
regions near the reacting bonds - against a best single site of -5.27. At
r_min = 1.0 the mean error falls to 0.527 A.

Shells
------
Default 3.0 / 4.0 / 5.0 A, replacing v1's 2.0 / 3.0 / 4.0. The offsets are
measured from the substrate vdW surface to a POINT. A molecular surrogate
carries its own radius, so a guanidinium carbon at 2 A standoff puts its
hydrogens through the substrate surface. The 2 A shell is retained as an option
(--shells) for point-charge work but is not the default.

Per-site metadata
-----------------
Each site carries:
  shell     nominal shell it was sampled from
  standoff  measured distance to the substrate vdW surface (union of R/TS/P)
  normal    outward unit normal of the shell at that point

The normal is what makes oriented surrogates tractable. A charged group in an
active site points its charged face at the substrate, so aligning the surrogate
axis to the local normal removes two of the three rotational degrees of
freedom; the residual spin about the normal is weak for near-symmetric groups
such as guanidinium. Orientation is therefore deferred to validation, and the
screen treats each site as a monopole - an approximation worth roughly
1 kcal/mol, which is stated rather than hidden.

Usage
-----
  python3 step_2_5_grid_v2.py [--r-min 1.0] [--shells 3.0 4.0 5.0]
                              [--density 6.0] [--seed 0] [--out grid_v2]

Requires sdf_grid.npz and shell_*.npz from steps 2.1 and 2.2. Shells not yet
extracted are built here from the SDF directly.
"""

import argparse
import itertools
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from skimage import measure


def load_sdf(path="sdf_grid.npz"):
    z = np.load(path, allow_pickle=True)
    return (z["sdf"], z["origin"].astype(float), float(z["voxel"]),
            z["xyz"].astype(float), z["radii"].astype(float))


def extract_shell(sdf, origin, voxel, level):
    verts, faces, normals, _ = measure.marching_cubes(sdf, level=level)
    return origin + verts * voxel, faces, normals


def sample_shell(verts, faces, normals, density, seed):
    mesh = trimesh.Trimesh(vertices=verts, faces=faces,
                           vertex_normals=normals, process=False)
    n = max(4, int(round(mesh.area * density)))
    pts, face_idx = trimesh.sample.sample_surface(mesh, n, seed=seed)
    # face normal at each sample; marching_cubes normals point along -grad SDF,
    # so flip to get the outward (away from substrate) direction
    fn = mesh.face_normals[face_idx]
    return np.asarray(pts), np.asarray(fn), mesh.area


def thin_poisson(pts, r_min, seed=0):
    """Global blue-noise thinning with a guaranteed floor. Grid-bucketed dart
    elimination; deterministic given the seed."""
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(pts))
    origin = pts.min(0)
    buckets, accepted, r2 = {}, [], r_min * r_min
    for idx in order:
        p = pts[idx]
        k = tuple(((p - origin) // r_min).astype(int))
        ok = True
        for d in itertools.product((-1, 0, 1), repeat=3):
            for q in buckets.get((k[0] + d[0], k[1] + d[1], k[2] + d[2]), ()):
                if ((pts[q] - p) ** 2).sum() < r2:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            accepted.append(idx)
            buckets.setdefault(k, []).append(idx)
    return np.array(accepted)


def standoff(pts, atoms, radii):
    return (np.sqrt(((pts[:, None, :] - atoms[None, :, :]) ** 2).sum(-1))
            - radii[None, :]).min(1)


def orient_outward(pts, normals, atoms):
    """Ensure normals point away from the substrate centroid."""
    c = atoms.mean(0)
    flip = ((pts - c) * normals).sum(1) < 0
    normals[flip] *= -1.0
    return normals


def capacity(pts, sep):
    """Greedy count of mutually compatible sites at a given separation -
    how many surrogates of that footprint could coexist."""
    tree = cKDTree(pts)
    used = np.zeros(len(pts), bool)
    n = 0
    for i in range(len(pts)):
        if used[i]:
            continue
        n += 1
        used[tree.query_ball_point(pts[i], sep)] = True
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r-min", type=float, default=1.0,
                    help="grid RESOLUTION, not placement separation")
    ap.add_argument("--shells", nargs="+", type=float, default=[3.0, 4.0, 5.0])
    ap.add_argument("--density", type=float, default=6.0,
                    help="dense pre-thinning samples per A^2")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="grid_v2")
    args = ap.parse_args()

    sdf, origin, voxel, atoms, radii = load_sdf()
    print(f"SDF {sdf.shape}, voxel {voxel} A, range [{sdf.min():.3f}, {sdf.max():.3f}] A")
    print(f"resolution r_min = {args.r_min} A   shells = {args.shells} A\n")

    pooled, tags, norms = [], [], []
    for d in args.shells:
        if not (sdf.min() < d < sdf.max()):
            print(f"shell {d}: SKIP (outside SDF range - increase margin in step 2.1)")
            continue
        verts, faces, vn = extract_shell(sdf, origin, voxel, d)
        pts, fn, area = sample_shell(verts, faces, vn, args.density, args.seed)
        V, F = len(verts), len(faces)
        chi = V - (F * 3 // 2) + F
        print(f"shell {d:.1f} A: area {area:7.1f} A^2, chi={chi}, "
              f"{len(pts):6d} dense samples")
        pooled.append(pts)
        norms.append(fn)
        tags += [d] * len(pts)

    pooled = np.vstack(pooled)
    norms = np.vstack(norms)
    tags = np.array(tags, float)

    keep = thin_poisson(pooled, args.r_min, seed=args.seed)
    pts, nrm, tag = pooled[keep], norms[keep], tags[keep]
    nrm = orient_outward(pts, nrm, atoms)
    so = standoff(pts, atoms, radii)

    d, _ = cKDTree(pts).query(pts, k=2)
    nn = d[:, 1]

    print(f"\nGRID v2: {len(pts)} candidate positions")
    for s in args.shells:
        m = tag == s
        if m.sum():
            print(f"  shell {s:.1f} A: {m.sum():5d} sites   "
                  f"standoff {so[m].min():.3f} - {so[m].max():.3f} A")
    print(f"  NN spacing  min {nn.min():.4f}  mean {nn.mean():.4f} "
          f"+/- {nn.std():.4f}  (CV {nn.std()/nn.mean():.4f}) A")
    print(f"  standoff overall {so.min():.3f} - {so.max():.3f} A")

    print("\nplacement capacity (independent sites at a given species separation):")
    for sep, who in [(3.0, ""), (3.5, "formate"), (4.5, "guanidinium"), (6.0, "")]:
        print(f"  sep {sep:.1f} A {who:14s}: {capacity(pts, sep):4d} sites")

    np.savez_compressed(f"{args.out}.npz", pts=pts, shell=tag, standoff=so,
                        normal=nrm, r_min=args.r_min, shells=np.array(args.shells),
                        density=args.density, seed=args.seed)
    with open(f"{args.out}.xyz", "w") as fh:
        fh.write(f"{len(pts)}\ncandidate grid v2 r_min={args.r_min} "
                 f"shells={args.shells} (resolution only; separation is an "
                 f"optimiser constraint)\n")
        for (x, y, z), s, o in zip(pts, tag, so):
            fh.write(f"X {x:.4f} {y:.4f} {z:.4f}   shell={s:.1f} standoff={o:.3f}\n")
    with open(f"{args.out}_sites.tsv", "w") as fh:
        fh.write("idx\tx\ty\tz\tshell\tstandoff\tnx\tny\tnz\n")
        for i, ((x, y, z), s, o, (a, b, c)) in enumerate(zip(pts, tag, so, nrm)):
            fh.write(f"{i}\t{x:.4f}\t{y:.4f}\t{z:.4f}\t{s:.1f}\t{o:.3f}"
                     f"\t{a:.4f}\t{b:.4f}\t{c:.4f}\n")
    print(f"\nsaved {args.out}.npz, {args.out}.xyz, {args.out}_sites.tsv")


if __name__ == "__main__":
    main()
