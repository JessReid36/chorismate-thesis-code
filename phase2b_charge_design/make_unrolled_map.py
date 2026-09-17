#!/usr/bin/env python3
"""
make_unrolled_map.py - unrolled difference-potential maps, per shell, with the
across-frame spread.

WHAT IS PLOTTED
    The candidate grid is a set of closed shells around the substrate, so it has
    no natural two-dimensional layout. Each shell is therefore unrolled onto
    spherical angles about the substrate centroid: the azimuth phi runs from -180
    to 180 degrees across the figure and the polar angle theta from 0 to 180 down
    it, with the difference potential as colour.

    Two panels per shell. The left shows the ensemble mean over the frames, which
    is the quantity a design would select against. The right shows the standard
    deviation across frames at the same sites, which is how much the value at a
    given position depends on which conformation the enzyme is in.

WHY UNROLL RATHER THAN PROJECT
    An equirectangular unrolling distorts area towards the poles, as any flat map
    of a sphere must. It is used here because it preserves the ordering of sites
    in both angles and is simple to read against a coordinate frame, not because
    it is area-true. Site positions are marked individually rather than
    interpolated into a continuous field, so no site appears where none was
    computed.

WHAT THE SPREAD MEANS
    The standard deviation is taken across frames at each site independently. It
    describes how much that site's value varies between conformations. It is not
    the standard error on the mean, which is smaller by the square root of the
    number of frames, and it says nothing about whether sites vary together:
    neighbouring sites almost certainly do, and that covariance is not shown
    here.

Usage:
    python3 make_unrolled_map.py --dv <mean map> --sites <sites tsv> \
        --vpot-dir <directory of vpot_*.out> --out <prefix>
"""
import argparse
import math
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = 627.5094740631


def read_sites(path):
    """Site index, coordinates, shell and standoff, by header name."""
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    hdr = [h.strip() for h in lines[0].lstrip("#").split("\t")]
    col = {h: i for i, h in enumerate(hdr)}
    need = ("idx", "x", "y", "z", "shell")
    for n in need:
        if n not in col:
            raise SystemExit(f"{path}: no column '{n}'; found {hdr}")
    out = []
    for l in lines[1:]:
        if l.startswith("#"):
            continue
        p = l.split("\t")
        out.append((int(p[col["idx"]]),
                    float(p[col["x"]]), float(p[col["y"]]), float(p[col["z"]]),
                    float(p[col["shell"]])))
    return out


def read_vpot(path):
    """One value per point, in the order the point list was submitted."""
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


def main():
    home = os.path.expanduser("~")
    root = f"{home}/system_development/phase2b_charge_design"
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", default=f"{root}/04_grid/grid_poisson_sites.tsv")
    ap.add_argument("--vpot-dir", default=f"{root}/08_frame_sensitivity")
    ap.add_argument("--out", default=f"{root}/08_frame_sensitivity/fig_unrolled")
    ap.add_argument("--frames", default="02450 04085 09025 09900 10775 "
                                        "11630 12485 14155 17505 19185")
    args = ap.parse_args()

    sites = read_sites(args.sites)
    n = len(sites)
    print(f"{n} sites from {os.path.basename(args.sites)}")

    frames = args.frames.split()
    maps = {}
    for f in frames:
        r = f"{args.vpot_dir}/vpot_{f}_R.out"
        t = f"{args.vpot_dir}/vpot_{f}_TS.out"
        if not (os.path.exists(r) and os.path.exists(t)):
            print(f"  {f}: missing, skipped"); continue
        vr, vt = read_vpot(r), read_vpot(t)
        if len(vr) != n or len(vt) != n:
            print(f"  {f}: {len(vr)} and {len(vt)} values against {n} sites, skipped")
            continue
        maps[f] = np.array([(vt[i] - vr[i]) * H for i in range(n)])
    if len(maps) < 2:
        raise SystemExit("need at least two frames")
    print(f"{len(maps)} frames: {' '.join(sorted(maps))}")

    M = np.vstack([maps[f] for f in sorted(maps)])
    mean = M.mean(0)
    sd = M.std(0, ddof=1)

    # angles about the centroid of the substrate, taken as the centroid of the
    # innermost shell's sites, which surrounds it
    xyz = np.array([[s[1], s[2], s[3]] for s in sites])
    shell = np.array([s[4] for s in sites])
    inner = xyz[shell == shell.min()]
    c = inner.mean(0)
    v = xyz - c
    r = np.linalg.norm(v, axis=1)
    theta = np.degrees(np.arccos(np.clip(v[:, 2] / r, -1, 1)))
    phi = np.degrees(np.arctan2(v[:, 1], v[:, 0]))

    shells = sorted(set(shell))
    print(f"shells: {shells}")

    for s in shells:
        m = shell == s
        fig, ax = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
        lim = max(abs(mean[m].min()), abs(mean[m].max()))
        a = ax[0].scatter(phi[m], theta[m], c=mean[m], s=42,
                          cmap="RdBu_r", vmin=-lim, vmax=lim,
                          edgecolors="none")
        ax[0].set_title(f"ensemble mean, {s:.0f} Å shell, {len(maps)} frames")
        cb = fig.colorbar(a, ax=ax[0]); cb.set_label(r"$\Delta v$ / kcal mol$^{-1}$ per $e$")

        b = ax[1].scatter(phi[m], theta[m], c=sd[m], s=42,
                          cmap="viridis", vmin=0, vmax=sd.max(),
                          edgecolors="none")
        ax[1].set_title(f"standard deviation across frames, {s:.0f} Å shell")
        cb = fig.colorbar(b, ax=ax[1]); cb.set_label(r"s.d. / kcal mol$^{-1}$")

        for a_ in ax:
            a_.set_xlabel(r"azimuth $\phi$ / degrees")
            a_.set_xlim(-180, 180); a_.set_ylim(180, 0)
            a_.set_xticks([-180, -90, 0, 90, 180])
            a_.set_yticks([0, 45, 90, 135, 180])
        ax[0].set_ylabel(r"polar angle $\theta$ / degrees")

        # mark the ten most stabilising sites of the mean map that lie on this
        # shell, since those are what a design would select from
        order = np.argsort(mean)
        top = [i for i in order[:10] if m[i]]
        if top:
            ax[0].scatter(phi[top], theta[top], s=150, facecolors="none",
                          edgecolors="k", linewidths=1.2)

        fig.tight_layout()
        out = f"{args.out}_shell{s:.0f}.png"
        fig.savefig(out, dpi=200)
        plt.close(fig)
        print(f"  wrote {out}   {m.sum()} sites, "
              f"mean {mean[m].min():.2f} to {mean[m].max():.2f}, "
              f"s.d. up to {sd[m].max():.2f}")

    with open(f"{args.out}_data.tsv", "w") as fh:
        fh.write(f"# unrolled difference potential, {len(maps)} frames\n")
        fh.write(f"# frames: {' '.join(sorted(maps))}\n")
        fh.write("# angles about the centroid of the innermost shell\n")
        fh.write("idx\tshell\ttheta_deg\tphi_deg\tdv_mean\tdv_sd\n")
        for i, s in enumerate(sites):
            fh.write(f"{s[0]}\t{s[4]:.1f}\t{theta[i]:.2f}\t{phi[i]:.2f}"
                     f"\t{mean[i]:.4f}\t{sd[i]:.4f}\n")
    print(f"  wrote {args.out}_data.tsv")
    print()
    print("The circled sites on the mean panel are the ten most stabilising of")
    print("the whole map that fall on that shell. Compare their positions with")
    print("the right panel: where the spread there is comparable to the")
    print("separation between candidates, the frames do not agree on the order.")


if __name__ == "__main__":
    main()
