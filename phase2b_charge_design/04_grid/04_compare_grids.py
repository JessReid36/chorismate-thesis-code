#!/usr/bin/env python3
"""
04_compare_grids.py - head-to-head comparison of candidate grids.

Reports, for each grid passed on the command line:

  N                 number of candidate positions
  min NN            smallest nearest-neighbour distance. For Poisson this is a
                    guaranteed floor; for CVT there is no floor at all.
  mean NN, CV       spacing and its coefficient of variation. CV is the metric
                    CVT exists to minimise, so it is the fair test of that
                    method's own claim.
  Coulomb cond      condition number of the inter-site 1/r matrix. This is the
                    quadratic term in the charge-selection problem; a large
                    value means near-degenerate sites the optimiser cannot tell
                    apart.
  max 1/r           worst-case inter-site Coulomb term.
  resolution        mean / p95 / max distance from an arbitrary position on the
                    shells to the nearest available site. This is the quantity
                    that limits how well a charge can be positioned, and it is
                    measured against a dense independent reference set.
  Dv penalty        resolution error multiplied by the measured Dv gradient,
                    i.e. the placement error expressed in kcal/mol.
  capacity          how many mutually compatible sites exist at a given species
                    separation - the number of surrogates that could coexist.

Usage
-----
  python3 04_compare_grids.py grid_poisson.npz grid_cvt_global.npz [...]

Optionally pass --dv <dv_grid.tsv> to calibrate the Dv gradient from a real
difference-potential map rather than the built-in default.
"""

import argparse
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from skimage import measure

HARTREE2KCAL = 627.5094740631
# Measured on the production Dv map (1448 sites, shells 3/4/5 A).
# Overridden by --dv, which is preferable: pass the map you are actually using.
DEFAULT_GRAD_MEAN, DEFAULT_GRAD_P95 = 0.185, 0.586


def reference_set(shells, n_per_shell=8000, seed=7):
    """Dense, independent sample of the shells: the 'true' positions a design
    might want. Deliberately seeded differently from the grids."""
    z = np.load("sdf_grid.npz", allow_pickle=True)
    sdf, origin, voxel = z["sdf"], z["origin"].astype(float), float(z["voxel"])
    out = []
    for d in shells:
        v, f, _, _ = measure.marching_cubes(sdf, level=d)
        m = trimesh.Trimesh(vertices=origin + v * voxel, faces=f, process=False)
        p, _ = trimesh.sample.sample_surface(m, n_per_shell, seed=seed)
        out.append(np.asarray(p))
    return np.vstack(out)


def dv_gradient(path):
    rows = [l.split("\t") for l in open(path).read().splitlines()[1:]]
    P = np.array([[float(r[1]), float(r[2]), float(r[3])] for r in rows])
    dv = np.array([float(r[7]) for r in rows]) * HARTREE2KCAL
    d, j = cKDTree(P).query(P, k=2)
    g = np.abs(dv - dv[j[:, 1]]) / d[:, 1]
    return float(g.mean()), float(np.percentile(g, 95))


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
    ap.add_argument("grids", nargs="+")
    ap.add_argument("--dv", default=None, help="dv_grid.tsv to calibrate the gradient")
    args = ap.parse_args()

    gmean, gp95 = (dv_gradient(args.dv) if args.dv
                   else (DEFAULT_GRAD_MEAN, DEFAULT_GRAD_P95))
    src = args.dv if args.dv else "built-in default (production map)"
    print(f"Dv gradient: mean {gmean:.3f}, p95 {gp95:.3f} kcal/mol/A per +1e "
          f"[{src}]\n")

    loaded = []
    for g in args.grids:
        z = np.load(g, allow_pickle=True)
        loaded.append((g.replace(".npz", ""), z["pts"], z["shell"]))

    shells = sorted(set(np.concatenate([s for _, _, s in loaded]).tolist()))
    print(f"reference set: dense independent sampling of shells {shells}\n")
    ref = reference_set(shells)

    hdr = (f"{'grid':<24}{'N':>6}{'minNN':>8}{'meanNN':>8}{'CV':>7}"
           f"{'cond':>10}{'max1/r':>8}{'res.mean':>10}{'res.p95':>9}{'res.max':>9}")
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for name, P, S in loaded:
        d, _ = cKDTree(P).query(P, k=2)
        nn = d[:, 1]
        D = np.sqrt(((P[:, None, :] - P[None, :, :]) ** 2).sum(-1))
        np.fill_diagonal(D, np.inf)
        C = 1.0 / D
        ev = np.abs(np.linalg.eigvalsh(C))
        ev = ev[ev > 1e-9]
        cond = ev.max() / ev.min()
        r, _ = cKDTree(P).query(ref, k=1)
        print(f"{name:<24}{len(P):>6}{nn.min():>8.3f}{nn.mean():>8.3f}"
              f"{nn.std()/nn.mean():>7.3f}{cond:>10.0f}{C.max():>8.3f}"
              f"{r.mean():>10.3f}{np.percentile(r,95):>9.3f}{r.max():>9.3f}")
        rows.append((name, P, S, nn, r))

    print(f"\n{'grid':<24}{'Dv penalty typ.':>18}{'steep-region':>15}"
          f"{'pairs<0.8A':>12}{'pairs<1.0A':>12}{'min NN':>9}{'p1 NN':>9}")
    print("-" * 99)
    for name, P, S, nn, r in rows:
        print(f"{name:<24}{gmean*r.mean():>18.3f}"
              f"{gp95*np.percentile(r,95):>15.3f}"
              f"{len(cKDTree(P).query_pairs(0.8)):>12,}"
              f"{len(cKDTree(P).query_pairs(1.0)):>12,}"
              f"{nn.min():>9.3f}{np.percentile(nn,1):>9.3f}")

    print(f"\n{'grid':<24}" + "".join(f"{f'sep {s}A':>12}" for s in [3.0, 3.5, 4.5, 6.0]))
    print("-" * 72)
    for name, P, S, nn, r in rows:
        print(f"{name:<24}" + "".join(f"{capacity(P, s):>12}"
                                      for s in [3.0, 3.5, 4.5, 6.0]))

    print(f"\n{'grid':<24}{'per-shell counts':>40}")
    print("-" * 64)
    for name, P, S, nn, r in rows:
        counts = " / ".join(f"{int((S==s).sum())}" for s in shells)
        print(f"{name:<24}{counts:>40}")

    print("\nNotes")
    print("  min NN is a guaranteed floor for Poisson-disk; CVT imposes none.")
    print("  CV is the metric CVT minimises, so compare it directly.")
    print("  resolution is what limits how well a charge can be positioned.")
    print("  cond is the quadratic term the optimiser sees; lower is better.")
    print("  cond depends on the smallest eigenvalue and so is sensitive to a")
    print("  few near-degenerate pairs. The pairs<0.8A / pairs<1.0A counts are")
    print("  the robust companions: quote those alongside it, not cond alone.")


if __name__ == "__main__":
    main()
