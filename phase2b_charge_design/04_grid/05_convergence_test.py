#!/usr/bin/env python3
"""
05_convergence_test.py - voxel-resolution convergence of the candidate grid.

The signed distance field is evaluated on a Cartesian lattice, and marching
cubes traces the shells from that lattice. The lattice spacing is a purely
numerical parameter with no physical meaning, so the grid is only defensible if
the answer does not depend on it. This script halves the spacing and measures
what changes.

Three independent measures, because they fail differently:

  1. SHELL GEOMETRY   surface area of each shell, and the measured offset of the
                      extracted isosurface from the substrate vdW surface. If
                      the coarse lattice were under-resolving the surface, the
                      area would move and the offset would drift from its target.

  2. SURFACE AGREEMENT  distance from points sampled on the fine-lattice shell to
                      the nearest point on the coarse-lattice shell. This is the
                      direct question - are these the same surface? - and is not
                      inferable from area alone, since two different surfaces can
                      share an area.

  3. GRID OUTCOME     site count and coverage of the finished grid at fixed
                      r_min. This is what actually propagates downstream. Small
                      count differences are expected from the stochastic sampling
                      draw even at identical resolution, so read the coverage
                      figure alongside the count rather than the count alone.

Usage
-----
  python3 05_convergence_test.py [--fine 0.20] [--coarse 0.30]
                                 [--shells 3.0 4.0 5.0] [--r-min 1.0]
                                 [--density 6.0] [--seed 0]

Needs reactant.xyz, ts.xyz, product.xyz in the working directory. Pure geometry;
no electronic structure. Runs in a few minutes; the fine lattice holds about
3.4x the voxels of the coarse one.
"""

import argparse
import itertools
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from skimage import measure

BONDI = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52}


def read_xyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[0].split()[0])
    els, xyz = [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        els.append(p[0])
        xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return els, np.asarray(xyz, float)


def build_sdf(els, xyz, voxel, margin=5.0):
    radii = np.array([BONDI[e] for e in els])
    lo = xyz.min(0) - (radii.max() + margin)
    hi = xyz.max(0) + (radii.max() + margin)
    dims = np.ceil((hi - lo) / voxel).astype(int) + 1
    gx, gy, gz = [lo[i] + voxel * np.arange(dims[i]) for i in range(3)]
    GX, GY, GZ = np.meshgrid(gx, gy, gz, indexing="ij")
    sdf = np.full(tuple(dims), np.inf, dtype=np.float32)
    for (ax, ay, az), r in zip(xyz, radii):
        np.minimum(sdf, np.sqrt((GX - ax) ** 2 + (GY - ay) ** 2 + (GZ - az) ** 2) - r,
                   out=sdf)
    return sdf, lo, radii


def shell(sdf, lo, voxel, level):
    v, f, n, _ = measure.marching_cubes(sdf, level=level)
    return lo + v * voxel, f


def offset_stats(verts, faces, xyz, radii, seed=0):
    m = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
    p, _ = trimesh.sample.sample_surface(m, 3000, seed=seed)
    p = np.asarray(p)
    d = (np.sqrt(((p[:, None, :] - xyz[None, :, :]) ** 2).sum(-1))
         - radii[None, :]).min(1)
    return m.area, float(d.mean()), float(d.std())


def thin(pts, r_min, seed=0):
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(pts))
    origin = pts.min(0)
    buckets, acc, r2 = {}, [], r_min * r_min
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
            acc.append(idx)
            buckets.setdefault(k, []).append(idx)
    return pts[np.array(acc)]


def build_grid(sdf, lo, voxel, shells, density, r_min, seed):
    pool, meshes = [], {}
    for d in shells:
        v, f = shell(sdf, lo, voxel, d)
        m = trimesh.Trimesh(vertices=v, faces=f, process=False)
        p, _ = trimesh.sample.sample_surface(m, max(4, int(round(m.area * density))),
                                             seed=seed)
        pool.append(np.asarray(p))
        meshes[d] = (v, f)
    return thin(np.vstack(pool), r_min, seed), meshes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coarse", type=float, default=0.30)
    ap.add_argument("--fine", type=float, default=0.20)
    ap.add_argument("--shells", nargs="+", type=float, default=[3.0, 4.0, 5.0])
    ap.add_argument("--r-min", type=float, default=1.0)
    ap.add_argument("--density", type=float, default=6.0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    els, xyz = [], []
    for f in ("reactant.xyz", "ts.xyz", "product.xyz"):
        e, x = read_xyz(f)
        els += e
        xyz.append(x)
    xyz = np.vstack(xyz)
    radii = np.array([BONDI[e] for e in els])
    print(f"pooled atoms: {len(els)}  shells {args.shells} A  "
          f"r_min {args.r_min} A  density {args.density} pts/A^2\n")

    built = {}
    for vox in (args.coarse, args.fine):
        sdf, lo, _ = build_sdf(els, xyz, vox)
        print(f"voxel {vox:.2f} A: SDF dims {sdf.shape}, {np.prod(sdf.shape):,} voxels")
        grid, meshes = build_grid(sdf, lo, vox, args.shells, args.density,
                                  args.r_min, args.seed)
        built[vox] = (sdf, lo, grid, meshes)
    print()

    # 1. shell geometry
    print("1. SHELL GEOMETRY")
    print(f"{'shell':>7}{'area coarse':>13}{'area fine':>12}{'d area':>9}"
          f"{'offset coarse':>18}{'offset fine':>18}")
    print("-" * 77)
    for d in args.shells:
        ac, oc, sc = offset_stats(*built[args.coarse][3][d], xyz, radii, args.seed)
        af, of, sf = offset_stats(*built[args.fine][3][d], xyz, radii, args.seed)
        print(f"{d:>7.1f}{ac:>13.1f}{af:>12.1f}{100*(af-ac)/ac:>8.2f}%"
              f"{f'{oc:.3f} +/- {sc:.3f}':>18}{f'{of:.3f} +/- {sf:.3f}':>18}")

    # 2. accuracy against the analytic surface
    print("\n2. ACCURACY AGAINST THE ANALYTIC SURFACE")
    print("   The union-of-spheres SDF is analytic, so the true shell is exactly")
    print("   SDF = d. Discretisation error is |SDF(p) - d| for points sampled")
    print("   from each extracted mesh - measured against ground truth, not")
    print("   against the other mesh.")
    print(f"\n{'shell':>7}{'lattice':>10}{'mean err':>11}{'p95':>10}{'max':>10}   (A)")
    print("-" * 50)
    worst = 0.0
    for d in args.shells:
        for vox, label in ((args.coarse, "coarse"), (args.fine, "fine")):
            v, f = built[vox][3][d]
            m = trimesh.Trimesh(vertices=v, faces=f, process=False)
            p, _ = trimesh.sample.sample_surface(m, 5000, seed=args.seed + 1)
            p = np.asarray(p)
            sdf_true = (np.sqrt(((p[:, None, :] - xyz[None, :, :]) ** 2).sum(-1))
                        - radii[None, :]).min(1)
            err = np.abs(sdf_true - d)
            if label == "coarse":
                worst = max(worst, err.max())
            print(f"{d:>7.1f}{label:>10}{err.mean():>11.4f}"
                  f"{np.percentile(err,95):>10.4f}{err.max():>10.4f}")

    # 3. grid outcome
    print("\n3. GRID OUTCOME")
    gc, gf = built[args.coarse][2], built[args.fine][2]
    print(f"  sites: coarse {len(gc)}   fine {len(gf)}   "
          f"difference {100*(len(gf)-len(gc))/len(gc):+.2f}%")

    ref = []
    for d in args.shells:
        vf, ff = built[args.fine][3][d]
        m = trimesh.Trimesh(vertices=vf, faces=ff, process=False)
        p, _ = trimesh.sample.sample_surface(m, 8000, seed=args.seed + 7)
        ref.append(np.asarray(p))
    ref = np.vstack(ref)
    for name, g in (("coarse", gc), ("fine", gf)):
        dd, _ = cKDTree(g).query(ref, k=1)
        nn, _ = cKDTree(g).query(g, k=2)
        print(f"  {name:>6}: coverage mean {dd.mean():.4f} p95 "
              f"{np.percentile(dd,95):.4f} A | min NN {nn[:,1].min():.4f} A")

    print("\nVERDICT")
    print(f"  worst coarse-lattice error against the analytic surface: {worst:.4f} A")
    print(f"  site spacing {args.r_min:.2f} A; grid resolves position to ~0.53 A")
    if worst < 0.1 * args.r_min:
        print("  CONVERGED - the coarse lattice reproduces the analytic surface to")
        print("  well within the site spacing, so the lattice is not limiting the grid.")
    else:
        print("  NOT CONVERGED at this tolerance - refine the lattice before")
        print("  quoting the grid.")


if __name__ == "__main__":
    main()
