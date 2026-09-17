#!/usr/bin/env python3
"""
make_perframe_maps.py - the difference potential at one shell, drawn separately
for each frame, to show how much the landscape depends on the conformation.

WHY
    A single frame's map is what a design would be built on if one structure were
    chosen. Test B6 found that the site ranking depends on which frame is used:
    Spearman correlations of 0.87 to 0.98 against the reference, with only three
    to ten of the ten most stabilising sites common between frames, and per-site
    differences up to 3.2 kcal/mol on a map whose whole range is 8.25.

    Those are numbers. This shows the same thing as a picture: the same shell,
    the same colour scale, one panel per frame, so the differences can be seen
    rather than read off a table. The final panel is the mean, which is what an
    ensemble-averaged design would select against.

WHAT TO LOOK FOR
    Whether the deep basin sits in the same place in every panel. Where it moves,
    a design built on one frame would place its charge somewhere another frame
    says is worse. Where every panel agrees, the position is robust to
    conformation and is a safer target.

CONSTRUCTION
    The shell is unrolled onto spherical angles about the substrate and the
    values interpolated onto a regular grid, with the azimuth tiled before
    gridding so the field is continuous across the periodic boundary. All panels
    share one colour scale, set by the largest absolute value across every frame,
    so the panels are directly comparable; a per-panel scale would hide exactly
    the differences this figure exists to show.

    The surface is interpolated between sites. Site positions are drawn on the
    mean panel so the sampling is visible.

Usage:
    python3 make_perframe_maps.py [--shell 3.0] [--out fig_perframe]
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


def grid_periodic(phi, theta, val, nphi=260, ntheta=130):
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
    ap.add_argument("--out", default="fig_perframe")
    ap.add_argument("--shell", type=float, default=3.0)
    ap.add_argument("--frames", default="02450 04085 09025 09900 10775 "
                                        "11630 12485 14155 17505 19185")
    ap.add_argument("--geom", default="",
                    help="substrate geometry used only to mark the "
                         "reacting atoms; any of the three works, since "
                         "they share a coordinate frame")
    ap.add_argument("--elev", type=float, default=42.0)
    ap.add_argument("--azim", type=float, default=-125.0)
    args = ap.parse_args()

    sites = read_sites(args.sites)
    n = len(sites)
    maps = {}
    for f in args.frames.split():
        r = f"{args.vpot_dir}/vpot_{f}_R.out"
        t = f"{args.vpot_dir}/vpot_{f}_TS.out"
        if not (os.path.exists(r) and os.path.exists(t)):
            print(f"  {f}: missing, skipped"); continue
        vr, vt = read_vpot(r), read_vpot(t)
        if len(vr) != n or len(vt) != n:
            print(f"  {f}: length mismatch, skipped"); continue
        maps[f] = np.array([(vt[i] - vr[i]) * H for i in range(n)])
    names = sorted(maps)
    if len(names) < 2:
        raise SystemExit("need at least two frames")
    print(f"{len(names)} frames at the {args.shell:.0f} A shell")

    shell = np.array([s[4] for s in sites])
    m = shell == args.shell
    if not m.any():
        raise SystemExit(f"no sites on the {args.shell} A shell; "
                         f"available: {sorted(set(shell))}")
    xyz = np.array([[s[1], s[2], s[3]] for s in sites])
    c = xyz[shell == shell.min()].mean(0)
    v = xyz - c
    rr = np.linalg.norm(v, axis=1)
    theta = np.degrees(np.arccos(np.clip(v[:, 2] / rr, -1, 1)))[m]
    phi = np.degrees(np.arctan2(v[:, 1], v[:, 0]))[m]

    # ------------------------------------- the four reacting atoms, projected
    # Delta v is a difference between two states, so it belongs to neither
    # geometry. The reacting atoms are marked instead, in the same angular frame
    # as the sites, because the landscape tracks where density moves between the
    # states: a basin sits where the transition state holds more negative charge
    # than the reactant. The atoms lie inside the shell, so each is placed at the
    # angle it would occupy if pushed radially outward.
    marks = []
    if args.geom and os.path.exists(args.geom):
        L = open(args.geom).read().splitlines()
        nat = int(L[0].split()[0])
        g = np.array([[float(x) for x in l.split()[1:4]] for l in L[2:2 + nat]])
        for i, lbl in ((7, "O3"), (8, "C4"), (0, "C1"), (12, "C6")):
            w = g[i] - c
            rw = np.linalg.norm(w)
            th = np.degrees(np.arccos(np.clip(w[2] / rw, -1, 1)))
            ph = np.degrees(np.arctan2(w[1], w[0]))
            marks.append((lbl, ph, th))
        print("  reacting atoms projected from %s:" % os.path.basename(args.geom))
        for lbl, ph, th in marks:
            print("    %-3s phi %+7.1f  theta %6.1f" % (lbl, ph, th))
    elif args.geom:
        print("  NOTE geometry not found at %s; atoms not marked" % args.geom)

    def annotate(ax):
        for lbl, ph, th in marks:
            ax.plot(ph, th, "o", ms=5, mfc="w", mec="k", mew=1.2, zorder=5)
            ax.annotate(lbl, (ph, th), textcoords="offset points",
                        xytext=(5, 4), fontsize=7, fontweight="bold", zorder=6,
                        bbox=dict(boxstyle="round,pad=0.15", fc="w", ec="none",
                                  alpha=0.75))

    M = np.vstack([maps[f][m] for f in names])
    mean = M.mean(0)
    lim = np.abs(M).max()
    print(f"  {m.sum()} sites, shared colour scale +/- {lim:.2f} kcal/mol")

    fields = []
    for f in names:
        P, T, Z = grid_periodic(phi, theta, maps[f][m])
        fields.append((f, Z))
    P, T, Zm = grid_periodic(phi, theta, mean)

    ncol = 4
    nrow = int(np.ceil((len(names) + 1) / ncol))

    # ------------------------------------------------ contour panel grid
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0*ncol, 2.6*nrow),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    lv = np.linspace(-lim, lim, 21)
    for ax, (f, Z) in zip(axes, fields):
        ax.contourf(P, T, Z, levels=lv, cmap="RdBu_r", extend="both")
        b = np.unravel_index(np.nanargmin(Z), Z.shape)
        ax.plot(P[b], T[b], "k*", ms=11, mfc="none", mew=1.4)
        annotate(ax)
        ax.set_title(f"frame {f}", fontsize=10)
        ax.set_xlim(-180, 180); ax.set_ylim(180, 0)
    axm = axes[len(fields)]
    cf = axm.contourf(P, T, Zm, levels=lv, cmap="RdBu_r", extend="both")
    axm.scatter(phi, theta, s=2, c="k", alpha=0.3, linewidths=0)
    b = np.unravel_index(np.nanargmin(Zm), Zm.shape)
    axm.plot(P[b], T[b], "k*", ms=11, mfc="none", mew=1.4)
    annotate(axm)
    axm.set_title(f"mean of {len(names)}", fontsize=10, fontweight="bold")
    axm.set_xlim(-180, 180); axm.set_ylim(180, 0)
    for ax in axes[len(fields)+1:]:
        ax.axis("off")
    for ax in axes[:len(fields)+1]:
        ax.set_xticks([-180, 0, 180]); ax.set_yticks([0, 90, 180])
    fig.supxlabel(r"azimuth $\phi$ / degrees")
    fig.supylabel(r"polar angle $\theta$ / degrees")
    sub = ("Difference potential on the %.0f Å shell, one panel per frame. "
           "The star marks the deepest basin." % args.shell)
    if marks:
        sub += (" Labelled points are the reacting atoms, projected radially "
                "onto the shell from inside it.")
    fig.suptitle(sub, fontsize=11)
    cb = fig.colorbar(cf, ax=axes.tolist(), shrink=0.8, pad=0.02)
    cb.set_label(r"$\Delta v$ / kcal mol$^{-1}$ per $e$")
    o = f"{args.out}_contour_shell{args.shell:.0f}.png"
    fig.savefig(o, dpi=190, bbox_inches="tight"); plt.close(fig)
    print(f"  wrote {o}")

    def mark3d(ax, Z, floor):
        """Mark the deepest basin and the ether oxygen on a surface panel.

        A point drawn at the surface height disappears whenever something in
        front of it is higher, so the basin is shown as a stem from the floor up
        to the surface, plus a marker on the floor itself. The floor marker is
        never occluded; the stem reads from most angles."""
        b = np.unravel_index(np.nanargmin(Z), Z.shape)
        bp, bt, bz = float(P[b]), float(T[b]), float(Z[b])
        ax.plot([bp, bp], [bt, bt], [floor, bz], color="k", lw=1.4, zorder=10)
        ax.scatter([bp], [bt], [bz], c="yellow", edgecolors="k",
                   s=55, linewidths=1.2, depthshade=False, zorder=11)
        ax.scatter([bp], [bt], [floor], c="k", marker="x", s=32,
                   depthshade=False, zorder=9)
        for lbl, ph, th in marks:
            if lbl != "O3":
                continue
            ax.plot([ph, ph], [th, th], [floor, 0.0], color="0.25", lw=1.0,
                    ls=":", zorder=9)
            ax.scatter([ph], [th], [floor], c="w", edgecolors="k", marker="s",
                       s=40, linewidths=1.1, depthshade=False, zorder=10)
            ax.text(ph, th, floor, "  O3", fontsize=8, fontweight="bold",
                    zorder=12)
        return bz

    # ------------------------------------------------ surface panel grid
    fig = plt.figure(figsize=(4.2*ncol, 3.1*nrow))
    for k, (f, Z) in enumerate(fields):
        ax = fig.add_subplot(nrow, ncol, k+1, projection="3d")
        ax.plot_surface(P, T, Z, cmap="RdBu_r", vmin=-lim, vmax=lim,
                        rstride=3, cstride=3, linewidth=0, antialiased=True)
        mark3d(ax, Z, -lim)
        ax.set_title(f"frame {f}", fontsize=10)
        ax.set_zlim(-lim, lim)
        ax.view_init(elev=args.elev, azim=args.azim)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_zlabel("", fontsize=7); ax.tick_params(labelsize=7)
    ax = fig.add_subplot(nrow, ncol, len(fields)+1, projection="3d")
    s = ax.plot_surface(P, T, Zm, cmap="RdBu_r", vmin=-lim, vmax=lim,
                        rstride=3, cstride=3, linewidth=0, antialiased=True)
    mark3d(ax, Zm, -lim)
    ax.set_title(f"mean of {len(names)}", fontsize=10, fontweight="bold")
    ax.set_zlim(-lim, lim)
    ax.view_init(elev=args.elev, azim=args.azim)
    ax.set_xticks([]); ax.set_yticks([]); ax.tick_params(labelsize=7)
    st = (f"Difference potential as a surface, {args.shell:.0f} Å shell. "
          f"Basins stabilise the transition state; peaks destabilise it. "
          f"All panels share one scale.")
    if marks:
        st += ("\nThe stem marks each panel's deepest basin; the square on the "
               "floor is the ether oxygen, projected radially onto the shell.")
    fig.suptitle(st, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    o = f"{args.out}_3d_shell{args.shell:.0f}.png"
    fig.savefig(o, dpi=170); plt.close(fig)
    print(f"  wrote {o}")

    # ------------------------------------------------ where is the best site?
    print()
    print("  deepest basin in each frame, as an angle and a value:")
    for f, Z in fields:
        b = np.unravel_index(np.nanargmin(Z), Z.shape)
        print(f"    frame {f}   phi {P[b]:+7.1f}  theta {T[b]:6.1f}   "
              f"{np.nanmin(Z):+.2f} kcal/mol")
    b = np.unravel_index(np.nanargmin(Zm), Zm.shape)
    print(f"    mean       phi {P[b]:+7.1f}  theta {T[b]:6.1f}   "
          f"{np.nanmin(Zm):+.2f} kcal/mol")
    print()
    print("  Where those angles differ between frames, a design built on one")
    print("  structure would place its charge where another structure says the")
    print("  effect is smaller. That is the case for an ensemble-averaged map.")


if __name__ == "__main__":
    main()
