#!/usr/bin/env python3
"""
make_consensus_map.py - one panel carrying both halves of the frame-dependence
argument.

WHAT THE PER-FRAME MAPS SHOWED
    Across ten frames at the 3 A shell, every deepest basin falls in the same
    region: the polar angle spans 23.7 to 39.1 degrees and the azimuth about 80
    degrees of one quadrant, adjacent to the ether oxygen O3 at -168.6, 25.3.
    That is the atom on which negative charge accumulates as the C4-O3 bond
    breaks, so the landscape agrees with the chemistry established independently
    by the sign-validation probes.

    What the frames do not agree on is the magnitude. The basin depth ranges from
    -3.42 to -5.56 kcal/mol, a variation of 62 per cent, and the site-wise
    standard deviation at the ten most stabilising sites is 0.66 to 0.72
    kcal/mol, which is above the resolution limit of the design objective.

    The case for an ensemble-averaged map is therefore not that single frames
    point in the wrong direction. It is that a single frame gives a misleading
    magnitude, and that within the broad favourable region the fine ranking
    cannot be resolved from one structure.

WHAT THIS FIGURE SHOWS
    Filled contours are the mean across frames. Overlaid line contours are the
    standard deviation across frames at the same positions: where those lines
    crowd over a basin, the frames disagree about a place a design would choose.
    Each frame's own deepest basin is marked, so the clustering is visible
    directly. The reacting atoms are labelled for orientation.

Usage:
    python3 make_consensus_map.py --shell 3.0 --geom <reactant.xyz>
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

H = 627.5094740631


def read_sites(path):
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    hdr = [h.strip() for h in lines[0].lstrip("#").split("\t")]
    col = {h: i for i, h in enumerate(hdr)}
    for n in ("idx", "x", "y", "z", "shell"):
        if n not in col:
            raise SystemExit(f"{path}: no column '{n}'")
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


def grid_periodic(phi, theta, val, nphi=320, ntheta=160):
    p = np.concatenate([phi - 360.0, phi, phi + 360.0])
    t = np.concatenate([theta, theta, theta])
    v = np.concatenate([val, val, val])
    P, T = np.meshgrid(np.linspace(-180, 180, nphi),
                       np.linspace(0, 180, ntheta))
    Z = griddata((p, t), v, (P, T), method="cubic")
    for meth in ("linear", "nearest"):
        if not np.isnan(Z).any():
            break
        Zf = griddata((p, t), v, (P, T), method=meth)
        Z = np.where(np.isnan(Z), Zf, Z)
    return P, T, Z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", default="grid_poisson_sites.tsv")
    ap.add_argument("--vpot-dir", default=".")
    ap.add_argument("--out", default="fig_consensus")
    ap.add_argument("--shell", type=float, default=3.0)
    ap.add_argument("--geom", default="")
    ap.add_argument("--frames", default="02450 04085 09025 09900 10775 "
                                        "11630 12485 14155 17505 19185")
    args = ap.parse_args()

    sites = read_sites(args.sites)
    n = len(sites)
    maps = {}
    for f in args.frames.split():
        r = f"{args.vpot_dir}/vpot_{f}_R.out"
        t = f"{args.vpot_dir}/vpot_{f}_TS.out"
        if not (os.path.exists(r) and os.path.exists(t)):
            continue
        vr, vt = read_vpot(r), read_vpot(t)
        if len(vr) != n or len(vt) != n:
            continue
        maps[f] = np.array([(vt[i] - vr[i]) * H for i in range(n)])
    names = sorted(maps)
    print(f"{len(names)} frames at the {args.shell:.0f} A shell")

    shell = np.array([s[4] for s in sites])
    m = shell == args.shell
    xyz = np.array([[s[1], s[2], s[3]] for s in sites])
    c = xyz[shell == shell.min()].mean(0)
    v = xyz - c
    rr = np.linalg.norm(v, axis=1)
    theta = np.degrees(np.arccos(np.clip(v[:, 2] / rr, -1, 1)))[m]
    phi = np.degrees(np.arctan2(v[:, 1], v[:, 0]))[m]

    M = np.vstack([maps[f][m] for f in names])
    mean = M.mean(0)
    sd = M.std(0, ddof=1)

    P, T, Zm = grid_periodic(phi, theta, mean)
    _, _, Zs = grid_periodic(phi, theta, sd)

    # each frame's own deepest basin, on the interpolated field
    basins = []
    for f in names:
        _, _, Z = grid_periodic(phi, theta, maps[f][m])
        b = np.unravel_index(np.nanargmin(Z), Z.shape)
        basins.append((f, float(P[b]), float(T[b]), float(np.nanmin(Z))))

    marks = []
    if args.geom and os.path.exists(args.geom):
        L = open(args.geom).read().splitlines()
        nat = int(L[0].split()[0])
        g = np.array([[float(x) for x in l.split()[1:4]] for l in L[2:2+nat]])
        for i, lbl in ((7, "O3"), (8, "C4"), (0, "C1"), (12, "C6")):
            w = g[i] - c
            rw = np.linalg.norm(w)
            marks.append((lbl,
                          float(np.degrees(np.arctan2(w[1], w[0]))),
                          float(np.degrees(np.arccos(np.clip(w[2]/rw, -1, 1))))))

    lim = float(np.nanmax(np.abs(Zm)))
    fig, ax = plt.subplots(figsize=(12, 6.0))
    lv = np.linspace(-lim, lim, 25)
    cf = ax.contourf(P, T, Zm, levels=lv, cmap="RdBu_r", extend="both")

    # the spread, as line contours over the mean
    sl = np.linspace(np.nanmin(Zs), np.nanmax(Zs), 7)[1:]
    cs = ax.contour(P, T, Zs, levels=sl, colors="k", linewidths=0.9,
                    linestyles="--", alpha=0.75)
    ax.clabel(cs, fmt="%.2f", fontsize=7, inline=True)

    # every frame's basin, then the mean's
    for f, ph, th, val in basins:
        ax.plot(ph, th, "o", ms=7, mfc="none", mec="k", mew=1.3, zorder=5)
    bm = np.unravel_index(np.nanargmin(Zm), Zm.shape)
    ax.plot(P[bm], T[bm], "*", ms=20, mfc="yellow", mec="k", mew=1.4, zorder=6)

    for lbl, ph, th in marks:
        ax.plot(ph, th, "s", ms=7, mfc="w", mec="k", mew=1.3, zorder=7)
        ax.annotate(lbl, (ph, th), textcoords="offset points", xytext=(7, 5),
                    fontsize=10, fontweight="bold", zorder=8,
                    bbox=dict(boxstyle="round,pad=0.2", fc="w", ec="none",
                              alpha=0.85))

    ax.set_xlim(-180, 180); ax.set_ylim(180, 0)
    ax.set_xticks([-180, -90, 0, 90, 180]); ax.set_yticks([0, 45, 90, 135, 180])
    ax.set_xlabel(r"azimuth $\phi$ / degrees")
    ax.set_ylabel(r"polar angle $\theta$ / degrees")
    depths = [b[3] for b in basins]
    ax.set_title(
        f"Difference potential on the {args.shell:.0f} Å shell, mean of "
        f"{len(names)} frames.\n"
        f"Dashed contours are the standard deviation across frames. "
        f"Circles mark each frame's deepest basin, the star the mean's. "
        f"Depths span {min(depths):.2f} to {max(depths):.2f} kcal mol$^{{-1}}$.",
        fontsize=11)
    cb = fig.colorbar(cf, ax=ax)
    cb.set_label(r"mean $\Delta v$ / kcal mol$^{-1}$ per $e$")
    fig.tight_layout()
    o = f"{args.out}_shell{args.shell:.0f}.png"
    fig.savefig(o, dpi=200); plt.close(fig)
    print(f"  wrote {o}")

    # ---------------------------------------------------- the numbers to quote
    ph = np.array([b[1] for b in basins])
    th = np.array([b[2] for b in basins])
    # azimuth is periodic; take the spread on the circle
    a = np.radians(ph)
    R = np.hypot(np.cos(a).mean(), np.sin(a).mean())
    circ = np.degrees(np.sqrt(-2 * np.log(R))) if R > 0 else float("nan")
    print()
    print(f"  basin positions: theta {th.min():.1f} to {th.max():.1f} degrees, "
          f"circular s.d. in azimuth {circ:.1f} degrees")
    print(f"  basin depths:    {min(depths):.2f} to {max(depths):.2f} kcal/mol, "
          f"a spread of {100*(max(depths)-min(depths))/abs(np.mean(depths)):.0f}%")
    if marks:
        o3 = [mk for mk in marks if mk[0] == "O3"]
        if o3:
            _, p3, t3 = o3[0]
            d = np.degrees(np.arccos(np.clip(
                np.sin(np.radians(th))*np.sin(np.radians(t3)) *
                np.cos(np.radians(ph - p3)) +
                np.cos(np.radians(th))*np.cos(np.radians(t3)), -1, 1)))
            print(f"  angular distance of each basin from O3: "
                  f"{d.min():.0f} to {d.max():.0f} degrees")
    print()
    print("  The frames agree on where a charge helps and disagree on how much.")
    print("  That is the case for designing against the mean with its spread")
    print("  attached, rather than against any one structure.")


if __name__ == "__main__":
    main()
