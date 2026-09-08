#!/usr/bin/env python3
"""
03_grid_cvt.py - candidate-position grid by restricted centroidal Voronoi
tessellation (CVT) with Lloyd relaxation.

This is the alternative to 02_grid_poisson.py, provided so the two placement
strategies can be compared on identical inputs. Both scripts read the same
sdf_grid.npz, extract the same shells, and sample the same dense cloud with the
same seed; only the point-selection step differs.

Method
------
For each shell, N seeds are drawn from that shell's dense cloud by
farthest-point sampling (deterministic given the seed), then relaxed:

  1. assign every dense-cloud point to its nearest seed  (Voronoi partition,
     restricted to the surface because the cloud lies on it)
  2. move each seed to the centroid of the points assigned to it
  3. project the centroid back onto the shell mesh, since the centroid of a
     patch of curved surface lies slightly inside it
  4. repeat until the mean seed displacement falls below --tol

That is Lloyd's algorithm restricted to a surface. It converges to a
centroidal tessellation in which each seed sits at the centre of mass of its
own cell and the cells have near-equal area, giving a regular, near-hexagonal
arrangement.

Two modes
---------
  per_shell  each shell relaxed independently (the textbook restricted CVT)
  global     seeds partition the POOLED cloud from all shells, so a seed
             competes with seeds on neighbouring shells as well as its own.
             Each seed is still projected back onto its own shell, preserving
             the per-shell counts.

The global mode exists because the shells sit ~1 A apart radially while
in-shell spacing is larger, so for many sites the nearest neighbour is on an
adjacent shell - something per-shell relaxation cannot see.

What CVT does and does not give you
-----------------------------------
Lloyd minimises the VARIANCE of the spacing. It places no lower bound on any
individual pair distance, because it moves seeds to centroids and a centroid
has no notion of a minimum distance. Poisson-disk elimination does the
opposite: a hard floor by construction, with irregular spacing above it.

Usage
-----
  python3 03_grid_cvt.py --counts 365 479 604 --shells 3.0 4.0 5.0 --mode global

Pass --counts to match a Poisson grid exactly (read them from its output), or
--total to have them apportioned by shell area.
"""

import argparse
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from skimage import measure


def load_sdf(path="sdf_grid.npz"):
    z = np.load(path, allow_pickle=True)
    return (z["sdf"], z["origin"].astype(float), float(z["voxel"]),
            z["xyz"].astype(float), z["radii"].astype(float))


def shell_mesh(sdf, origin, voxel, level):
    verts, faces, normals, _ = measure.marching_cubes(sdf, level=level)
    return origin + verts * voxel, faces, normals


def dense_cloud(verts, faces, normals, density, seed):
    """Identical sampling to 02_grid_poisson.py, so both methods see the same
    candidate population."""
    mesh = trimesh.Trimesh(vertices=verts, faces=faces,
                           vertex_normals=normals, process=False)
    n = max(4, int(round(mesh.area * density)))
    pts, face_idx = trimesh.sample.sample_surface(mesh, n, seed=seed)
    return np.asarray(pts), np.asarray(mesh.face_normals[face_idx]), mesh.area


def farthest_point_seeds(pts, n, seed=0):
    """Deterministic, well-spread starting seeds. A better Lloyd initialisation
    than a random draw, and reproducible."""
    rng = np.random.default_rng(seed)
    idx = [int(rng.integers(len(pts)))]
    d2 = ((pts - pts[idx[0]]) ** 2).sum(1)
    for _ in range(n - 1):
        nxt = int(np.argmax(d2))
        idx.append(nxt)
        d2 = np.minimum(d2, ((pts - pts[nxt]) ** 2).sum(1))
    return pts[np.array(idx)].copy()


def project(pts, verts, vtree):
    """Snap points back onto the shell. Vertex-nearest projection suffices:
    mesh vertex spacing is ~0.25 A, well below the site spacing."""
    _, j = vtree.query(pts, k=1)
    return verts[j]


def lloyd_per_shell(cloud, verts, n, iters, tol, seed):
    vtree = cKDTree(verts)
    seeds = project(farthest_point_seeds(cloud, n, seed), verts, vtree)
    for it in range(iters):
        _, owner = cKDTree(seeds).query(cloud, k=1)
        new = seeds.copy()
        for k in range(n):
            m = cloud[owner == k]
            if len(m):
                new[k] = m.mean(0)
        new = project(new, verts, vtree)
        shift = float(np.linalg.norm(new - seeds, axis=1).mean())
        seeds = new
        if shift < tol:
            break
    return seeds, it + 1


def lloyd_global(clouds, shells, verts_by_shell, counts, iters, tol, seed):
    pooled = np.vstack(clouds)
    cloud_shell = np.concatenate([[s] * len(c) for s, c in zip(shells, clouds)])
    vtrees = {s: cKDTree(v) for s, v in verts_by_shell.items()}

    seeds, seed_shell = [], []
    for s, c, n in zip(shells, clouds, counts):
        seeds.append(project(farthest_point_seeds(c, n, seed),
                             verts_by_shell[s], vtrees[s]))
        seed_shell += [s] * n
    seeds = np.vstack(seeds)
    seed_shell = np.array(seed_shell)

    for it in range(iters):
        _, owner = cKDTree(seeds).query(pooled, k=1)
        new = seeds.copy()
        for k in range(len(seeds)):
            sel = owner == k
            if not sel.any():
                continue
            m = pooled[sel]
            # keep the centroid on the seed's own shell
            same = m[cloud_shell[sel] == seed_shell[k]]
            use = same if len(same) else m
            new[k] = use.mean(0)
        for s in shells:
            msk = seed_shell == s
            new[msk] = project(new[msk], verts_by_shell[s], vtrees[s])
        shift = float(np.linalg.norm(new - seeds, axis=1).mean())
        seeds = new
        if shift < tol:
            break
    return seeds, seed_shell, it + 1


def standoff(pts, atoms, radii):
    return (np.sqrt(((pts[:, None, :] - atoms[None, :, :]) ** 2).sum(-1))
            - radii[None, :]).min(1)


def nearest_normal(pts, cloud, cloud_normals):
    _, j = cKDTree(cloud).query(pts, k=1)
    return cloud_normals[j]


def orient_outward(pts, normals, atoms):
    c = atoms.mean(0)
    flip = ((pts - c) * normals).sum(1) < 0
    normals[flip] *= -1.0
    return normals


def capacity(pts, sep):
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
    ap.add_argument("--shells", nargs="+", type=float, default=[3.0, 4.0, 5.0])
    ap.add_argument("--counts", nargs="+", type=int, default=None,
                    help="points per shell; match the Poisson grid for a fair test")
    ap.add_argument("--total", type=int, default=None,
                    help="total points, apportioned by shell area (if --counts absent)")
    ap.add_argument("--density", type=float, default=6.0)
    ap.add_argument("--mode", choices=["per_shell", "global"], default="global")
    ap.add_argument("--iters", type=int, default=300)
    ap.add_argument("--tol", type=float, default=1e-4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    sdf, origin, voxel, atoms, radii = load_sdf()
    base = args.out or f"grid_cvt_{args.mode}"
    print(f"restricted CVT (Lloyd), mode={args.mode}, seed={args.seed}, "
          f"tol={args.tol}, max {args.iters} iters")
    print(f"shells = {args.shells} A, cloud density = {args.density} pts/A^2\n")

    clouds, cnormals, verts_by_shell, areas = [], [], {}, []
    for d in args.shells:
        if not (sdf.min() < d < sdf.max()):
            raise SystemExit(f"shell {d} A lies outside the SDF range "
                             f"[{sdf.min():.2f}, {sdf.max():.2f}] - "
                             f"increase the margin in 01_build_sdf.py")
        v, f, vn = shell_mesh(sdf, origin, voxel, d)
        c, cn, area = dense_cloud(v, f, vn, args.density, args.seed)
        V, F = len(v), len(f)
        print(f"shell {d:.1f} A: area {area:7.1f} A^2, chi={V - (F*3//2) + F}, "
              f"{len(c):6d} dense samples")
        clouds.append(c)
        cnormals.append(cn)
        verts_by_shell[d] = v
        areas.append(area)

    if args.counts:
        counts = args.counts
        if len(counts) != len(args.shells):
            raise SystemExit("--counts needs one entry per shell")
    elif args.total:
        w = np.array(areas) / sum(areas)
        counts = np.maximum(4, np.round(w * args.total).astype(int)).tolist()
    else:
        raise SystemExit("give --counts (preferred, to match the Poisson grid) "
                         "or --total")
    print(f"\ntarget counts per shell: {counts}  (total {sum(counts)})\n")

    if args.mode == "per_shell":
        pts_list, tags = [], []
        for d, n, cloud in zip(args.shells, counts, clouds):
            seeds, used = lloyd_per_shell(cloud, verts_by_shell[d], n,
                                          args.iters, args.tol, args.seed)
            dd, _ = cKDTree(seeds).query(seeds, k=2)
            print(f"shell {d:.1f} A: {n:5d} points, {used:3d} iters | "
                  f"within-shell NN min {dd[:,1].min():.3f} mean {dd[:,1].mean():.3f}")
            pts_list.append(seeds)
            tags += [d] * n
        pts, tag = np.vstack(pts_list), np.array(tags, float)
    else:
        pts, tag, used = lloyd_global(clouds, args.shells, verts_by_shell,
                                      counts, args.iters, args.tol, args.seed)
        print(f"pooled relaxation: {len(pts)} points, {used} iters")
        for d in args.shells:
            m = tag == d
            dd, _ = cKDTree(pts[m]).query(pts[m], k=2)
            print(f"  shell {d:.1f} A: {m.sum():5d} points | "
                  f"within-shell NN min {dd[:,1].min():.3f} "
                  f"mean {dd[:,1].mean():.3f}")

    allcloud = np.vstack(clouds)
    allnorm = np.vstack(cnormals)
    nrm = orient_outward(pts, nearest_normal(pts, allcloud, allnorm), atoms)
    so = standoff(pts, atoms, radii)
    d, _ = cKDTree(pts).query(pts, k=2)
    nn = d[:, 1]

    print(f"\nCVT GRID ({args.mode}): {len(pts)} candidate positions")
    print(f"  NN spacing  min {nn.min():.4f}  mean {nn.mean():.4f} "
          f"+/- {nn.std():.4f}  (CV {nn.std()/nn.mean():.4f}) A")
    print(f"  standoff {so.min():.3f} - {so.max():.3f} A")
    print("\nplacement capacity (independent sites at a given species separation):")
    for sep, who in [(3.0, ""), (3.5, "formate"), (4.5, "guanidinium"), (6.0, "")]:
        print(f"  sep {sep:.1f} A {who:14s}: {capacity(pts, sep):4d} sites")

    np.savez_compressed(f"{base}.npz", pts=pts, shell=tag, standoff=so,
                        normal=nrm, method=f"restricted_CVT_Lloyd_{args.mode}",
                        shells=np.array(args.shells), counts=np.array(counts),
                        density=args.density, seed=args.seed)
    with open(f"{base}.xyz", "w") as fh:
        fh.write(f"{len(pts)}\ncandidate grid (restricted CVT, {args.mode}) "
                 f"shells={args.shells} counts={counts}\n")
        for (x, y, z), s, o in zip(pts, tag, so):
            fh.write(f"X {x:.4f} {y:.4f} {z:.4f}   shell={s:.1f} standoff={o:.3f}\n")
    with open(f"{base}_sites.tsv", "w") as fh:
        fh.write("idx\tx\ty\tz\tshell\tstandoff\tnx\tny\tnz\n")
        for i, ((x, y, z), s, o, (a, b, c)) in enumerate(zip(pts, tag, so, nrm)):
            fh.write(f"{i}\t{x:.4f}\t{y:.4f}\t{z:.4f}\t{s:.1f}\t{o:.3f}"
                     f"\t{a:.4f}\t{b:.4f}\t{c:.4f}\n")
    print(f"\nsaved {base}.npz, {base}.xyz, {base}_sites.tsv")


if __name__ == "__main__":
    main()
