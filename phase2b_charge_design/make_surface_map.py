#!/usr/bin/env python3
"""
make_surface_map.py - the difference potential as an interpolated surface over
the unrolled shell, in filled-contour and three-dimensional forms.

RELATION TO make_unrolled_map.py
    That script marks each site individually, which is the honest representation
    of where values were computed. This one interpolates between sites to produce
    a continuous field, which reads more naturally as a landscape but shows
    values at positions where nothing was calculated.

    Interpolation is defensible here because the sampling is dense relative to
    the structure being shown: 365 sites on the innermost shell at a guaranteed
    minimum separation of 1.0 A, against a potential that varies smoothly over
    several angstroms. The site positions are overlaid on the contour figure so a
    reader can see the sampling for themselves, and any caption should say the
    surface is interpolated.

THE SEAM AT THE WRAP
    The azimuth is periodic at plus and minus 180 degrees. An interpolator given
    the raw angles produces a discontinuity there, because it does not know the
    left and right edges are adjacent. The data are therefore tiled once to each
    side before gridding and the result cropped back, so the field is continuous
    across the wrap.

    The polar angle is not periodic and needs no such treatment, but the sampling
    thins towards the poles in this projection, so the interpolated surface is
    least reliable at the top and bottom of the figure.

WHAT IS PLOTTED
    One figure per shell. The filled contours show the ensemble mean across
    frames, with the site positions marked and the ten most stabilising sites of
    the whole map circled. A second figure shows the same field as a
    three-dimensional surface, with the potential as height, so that stabilising
    regions appear as basins and destabilising ones as peaks.

Usage:
    python3 make_surface_map.py --sites <tsv> --vpot-dir <dir> --out <prefix>
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import griddata

H = 627.5094740631


def read_sites(path):
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    hdr = [h.strip() for h in lines[0].lstrip("#").split("\t")]
    col = {h: i for i, h in enumerate(hdr)}
    for n in ("idx", "x", "y", "z", "shell"):
        if n not in col:
            raise SystemExit(f"{path}: no column '{n}'; found {hdr}")
    out = []
    for l in lines[1:]:
        if l.startswith("#"):
            continue
        p = l.split("\t")
        out.append((int(p[col["idx"]]), float(p[col["x"]]), float(p[col["y"]]),
                    float(p[col["z"]]), float(p[col["shell"]])))
    return out


def read_vpot(path):
    out = []
    for l in open(path):
        s = l.split()
        if len(s) == 1:
            continue
        try:
            out.append(float(s[-1]))
        except ValueError:
            continue
    return out


def grid_periodic(phi, theta, val, nphi=360, ntheta=180):
    """Interpolate onto a regular grid, tiling in azimuth so the field is
    continuous across the wrap at plus and minus 180 degrees."""
    p = np.concatenate([phi - 360.0, phi, phi + 360.0])
    t = np.concatenate([theta, theta, theta])
    v = np.concatenate([val, val, val])
    gp = np.linspace(-180, 180, nphi)
    gt = np.linspace(0, 180, ntheta)
    P, T = np.meshgrid(gp, gt)
    Z = griddata((p, t), v, (P, T), method="cubic")
    # cubic leaves gaps at the convex hull edges; fill them from a linear pass
    if np.isnan(Z).any():
        Zl = griddata((p, t), v, (P, T), method="linear")
        Z = np.where(np.isnan(Z), Zl, Z)
    if np.isnan(Z).any():
        Zn = griddata((p, t), v, (P, T), method="nearest")
        Z = np.where(np.isnan(Z), Zn, Z)
    return P, T, Z


def main():
    home = os.path.expanduser("~")
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", default="grid_poisson_sites.tsv")
    ap.add_argument("--vpot-dir", default=".")
    ap.add_argument("--out", default="fig_surface")
    ap.add_argument("--frames", default="02450 04085 09025 09900 10775 "
                                        "11630 12485 14155 17505 19185")
    ap.add_argument("--elev", type=float, default=38.0,
                    help="viewing elevation for the three-dimensional figure")
    ap.add_argument("--azim", type=float, default=-125.0)
    args = ap.parse_args()

    sites = read_sites(args.sites)
    n = len(sites)
    print(f"{n} sites from {os.path.basename(args.sites)}")

    maps = {}
    for f in args.frames.split():
        r = f"{args.vpot_dir}/vpot_{f}_R.out"
        t = f"{args.vpot_dir}/vpot_{f}_TS.out"
        if not (os.path.exists(r) and os.path.exists(t)):
            print(f"  {f}: missing, skipped"); continue
        vr, vt = read_vpot(r), read_vpot(t)
        if len(vr) != n or len(vt) != n:
            print(f"  {f}: {len(vr)} and {len(vt)} against {n} sites, skipped")
            continue
        maps[f] = np.array([(vt[i] - vr[i]) * H for i in range(n)])
    if len(maps) < 2:
        raise SystemExit("need at least two frames")
    print(f"{len(maps)} frames: {' '.join(sorted(maps))}")

    M = np.vstack([maps[f] for f in sorted(maps)])
    mean = M.mean(0)

    xyz = np.array([[s[1], s[2], s[3]] for s in sites])
    shell = np.array([s[4] for s in sites])
    c = xyz[shell == shell.min()].mean(0)
    v = xyz - c
    rr = np.linalg.norm(v, axis=1)
    theta = np.degrees(np.arccos(np.clip(v[:, 2] / rr, -1, 1)))
    phi = np.degrees(np.arctan2(v[:, 1], v[:, 0]))

    order = np.argsort(mean)
    top10 = set(order[:10].tolist())

    for s in sorted(set(shell)):
        m = shell == s
        P, T, Z = grid_periodic(phi[m], theta[m], mean[m])
        lim = np.nanmax(np.abs(Z))

        # ---------------------------------------------- filled contours
        fig, ax = plt.subplots(figsize=(11, 5.6))
        lv = np.linspace(-lim, lim, 25)
        cf = ax.contourf(P, T, Z, levels=lv, cmap="RdBu_r", extend="both")
        ax.contour(P, T, Z, levels=lv[::4], colors="k", linewidths=0.35, alpha=0.5)
        ax.scatter(phi[m], theta[m], s=3, c="k", alpha=0.35, linewidths=0)
        tp = [i for i in top10 if m[i]]
        if tp:
            ax.scatter(phi[tp], theta[tp], s=140, facecolors="none",
                       edgecolors="k", linewidths=1.4)
        ax.set_xlim(-180, 180); ax.set_ylim(180, 0)
        ax.set_xticks([-180, -90, 0, 90, 180]); ax.set_yticks([0, 45, 90, 135, 180])
        ax.set_xlabel(r"azimuth $\phi$ / degrees")
        ax.set_ylabel(r"polar angle $\theta$ / degrees")
        ax.set_title(f"difference potential, {s:.0f} Å shell, "
                     f"mean of {len(maps)} frames")
        cb = fig.colorbar(cf, ax=ax)
        cb.set_label(r"$\Delta v$ / kcal mol$^{-1}$ per $e$")
        fig.tight_layout()
        o = f"{args.out}_contour_shell{s:.0f}.png"
        fig.savefig(o, dpi=200); plt.close(fig)
        print(f"  wrote {o}")

        # ------------------------------------------ three-dimensional surface
        fig = plt.figure(figsize=(11, 7))
        ax = fig.add_subplot(111, projection="3d")
        surf = ax.plot_surface(P, T, Z, cmap="RdBu_r", vmin=-lim, vmax=lim,
                               rstride=2, cstride=2, linewidth=0,
                               antialiased=True, alpha=0.97)
        ax.contour(P, T, Z, levels=12, zdir="z",
                   offset=float(np.nanmin(Z)) - 0.4, cmap="RdBu_r",
                   linewidths=0.6)
        ax.set_xlabel(r"azimuth $\phi$ / degrees")
        ax.set_ylabel(r"polar angle $\theta$ / degrees")
        ax.set_zlabel(r"$\Delta v$ / kcal mol$^{-1}$ per $e$")
        ax.set_title(f"difference potential as a surface, {s:.0f} Å shell, "
                     f"mean of {len(maps)} frames")
        ax.view_init(elev=args.elev, azim=args.azim)
        ax.set_zlim(float(np.nanmin(Z)) - 0.4, float(np.nanmax(Z)) + 0.2)
        fig.colorbar(surf, ax=ax, shrink=0.6, pad=0.08,
                     label=r"$\Delta v$ / kcal mol$^{-1}$ per $e$")
        fig.tight_layout()
        o = f"{args.out}_3d_shell{s:.0f}.png"
        fig.savefig(o, dpi=200); plt.close(fig)
        print(f"  wrote {o}   range {np.nanmin(Z):.2f} to {np.nanmax(Z):.2f}")

    print()
    print("The surfaces are interpolated between sites; the contour figures show")
    print("the site positions so the sampling is visible. Basins are positions")
    print("where a positive charge lowers the barrier, peaks where it raises it.")
    print("The azimuth is periodic and has been tiled before gridding, so there")
    print("is no seam at the edges; the polar sampling thins towards the poles,")
    print("so the top and bottom of each figure are the least reliable.")


if __name__ == "__main__":
    main()
