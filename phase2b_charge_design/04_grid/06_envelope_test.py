#!/usr/bin/env python3
"""
06_envelope_test.py - justification for building the exclusion envelope from the
union of the reactant, transition-state and product geometries.

The designed charges are static while the substrate reacts, so every candidate
position must clear the substrate at every point along the reaction path. An
envelope built from one geometry alone will contain positions the molecule
subsequently moves into.

The script builds a grid from each single geometry and from the union, all at
identical settings, then asks of each single-geometry grid: how many of its
sites fall below the intended standoff once the other two geometries are taken
into account?

It also reports which geometry actually sets the envelope at each site of the
union grid. If one geometry dominated, the union would be unnecessary.

Usage
-----
  python3 06_envelope_test.py [--shells 3.0 4.0 5.0] [--r-min 1.0]
                              [--density 6.0] [--voxel 0.30] [--seed 0]

Defaults are the production settings. Needs reactant.xyz, ts.xyz, product.xyz.
Pure geometry; no electronic structure.
"""

import argparse
import itertools
import numpy as np
import trimesh
from skimage import measure

BONDI = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52}
NAMES = ["reactant", "ts", "product"]


def read_xyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[0].split()[0])
    els, xyz = [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        els.append(p[0])
        xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return els, np.asarray(xyz, float)


def sdf_field(els, xyz, voxel, margin=5.0):
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
    return sdf, lo


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


def build(els, xyz, shells, voxel, density, r_min, seed):
    sdf, lo = sdf_field(els, xyz, voxel)
    pool = []
    for d in shells:
        if not (sdf.min() < d < sdf.max()):
            raise SystemExit(f"shell {d} outside SDF range")
        v, f, _, _ = measure.marching_cubes(sdf, level=d)
        m = trimesh.Trimesh(vertices=lo + v * voxel, faces=f, process=False)
        p, _ = trimesh.sample.sample_surface(m, max(4, int(round(m.area * density))),
                                             seed=seed)
        pool.append(np.asarray(p))
    return thin(np.vstack(pool), r_min, seed)


def standoff(pts, xyz, radii):
    return (np.sqrt(((pts[:, None, :] - xyz[None, :, :]) ** 2).sum(-1))
            - radii[None, :]).min(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shells", nargs="+", type=float, default=[3.0, 4.0, 5.0])
    ap.add_argument("--r-min", type=float, default=1.0)
    ap.add_argument("--density", type=float, default=6.0)
    ap.add_argument("--voxel", type=float, default=0.30)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    floor = min(args.shells)
    G = {}
    for name in NAMES:
        e, x = read_xyz(f"{name}.xyz")
        G[name] = (np.array([BONDI[q] for q in e]), x, e)

    print(f"shells {args.shells} A   r_min {args.r_min} A   "
          f"density {args.density} pts/A^2   voxel {args.voxel} A")
    print(f"intended minimum standoff: {floor:.1f} A\n")

    grids = {}
    for keys, label in [(["reactant"], "reactant only"),
                        (["ts"], "TS only"),
                        (["product"], "product only"),
                        (NAMES, "union R+TS+P")]:
        els = sum([G[k][2] for k in keys], [])
        xyz = np.vstack([G[k][1] for k in keys])
        grids[label] = (build(els, xyz, args.shells, args.voxel,
                              args.density, args.r_min, args.seed), keys)

    print("1. PATH CLEARANCE")
    print("   For each envelope, the standoff of its own sites measured against")
    print("   the geometries it did not include.\n")
    print(f"{'envelope':<18}{'sites':>7}{'min standoff':>15}{'below floor':>16}"
          f"{'worst':>9}")
    print("-" * 66)
    for label, (P, keys) in grids.items():
        others = [k for k in NAMES if k not in keys]
        if not others:
            print(f"{label:<18}{len(P):>7}{'-':>15}"
                  f"{'0 (by construction)':>16}{'-':>9}")
            continue
        mins = np.min([standoff(P, G[k][1], G[k][0]) for k in others], axis=0)
        bad = int((mins < floor - 1e-6).sum())
        print(f"{label:<18}{len(P):>7}{mins.min():>15.3f}"
              f"{f'{bad} ({100*bad/len(P):.1f}%)':>16}{mins.min():>9.3f}")

    print("\n2. WHICH GEOMETRY SETS THE ENVELOPE")
    print("   At each site of the union grid, which geometry is closest. If one")
    print("   dominated, the union would be unnecessary.\n")
    U = grids["union R+TS+P"][0]
    S = np.array([standoff(U, G[k][1], G[k][0]) for k in NAMES])
    who = np.argmin(S, axis=0)
    for i, n in enumerate(NAMES):
        c = int((who == i).sum())
        print(f"   {n:<10} binds at {c:>5} / {len(U)} sites  ({100*c/len(U):.1f}%)")

    print("\n3. SUBSTRATE MOTION ALONG THE PATH")
    R, T, P = G["reactant"][1], G["ts"][1], G["product"][1]
    for a, b, lab in ((R, T, "R -> TS"), (R, P, "R -> P")):
        d = np.linalg.norm(a - b, axis=1)
        print(f"   {lab:<9} per-atom displacement: mean {d.mean():.3f}  "
              f"max {d.max():.3f} A")

    print("\n4. COST OF THE UNION")
    nu = len(grids["union R+TS+P"][0])
    nr = len(grids["reactant only"][0])
    print(f"   union {nu} sites vs reactant-only {nr}: "
          f"{100*(nu-nr)/nr:+.1f}% in site count")
    Pr, _ = grids["reactant only"]
    mins = np.min([standoff(Pr, G[k][1], G[k][0]) for k in ("ts", "product")], axis=0)
    nbreach = int((mins < floor - 1e-6).sum())
    print(f"   reactant-only admits {nbreach} sites ({100*nbreach/nr:.1f}%) that "
          f"breach the {floor:.1f} A floor later on the path")

    # Sites present on the reactant-only grid but absent from the union grid.
    # This is a different and larger quantity than the count above: a site is
    # lost either because the union field excludes it, or because the union
    # grid's stochastic thinning did not select that position. Both are
    # reported so the write-up can quote the measure it means.
    from scipy.spatial import cKDTree
    Pu = grids["union R+TS+P"][0]
    tol = 1e-3
    dmatch, _ = cKDTree(Pu).query(Pr, k=1)
    nlost = int((dmatch > tol).sum())
    print(f"   reactant-only sites with no counterpart on the union grid "
          f"(within {tol} A): {nlost} ({100*nlost/nr:.1f}%)")
    if nlost:
        so_r = standoff(Pr, G["reactant"][1], G["reactant"][0])
        lost = so_r[dmatch > tol]
        print(f"   those sites stand {lost.min():.3f} to {lost.max():.3f} A from "
              f"the reactant surface")
    print("\n   The union is the cheaper of the two errors: a small loss of")
    print("   candidate positions against a guarantee of path-wide clearance.")


if __name__ == "__main__":
    main()
