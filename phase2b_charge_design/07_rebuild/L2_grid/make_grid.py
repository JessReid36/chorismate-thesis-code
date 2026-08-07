#!/usr/bin/env python3
"""Layer 2 -- candidate GRID with a principled LJ-equilibrium standoff."""
import os
import numpy as np

HERE   = os.path.dirname(os.path.abspath(__file__))
INP    = os.path.join(HERE, "..", "inputs")
REACT  = "reactant.xyz"
SPACING = 0.7
OUTER   = 9.0
QSIG    = 3.40
SIG = {"C": 3.40, "N": 3.25, "O": 2.96, "H": 1.06}
TWO16 = 2.0 ** (1.0 / 6.0)


def load_substrate():
    L = open(os.path.join(INP, REACT)).read().splitlines()[2:26]
    els = [l.split()[0] for l in L]
    R = np.array([[float(v) for v in l.split()[1:4]] for l in L])
    return els, R


def main():
    els, R = load_substrate()
    rmin_atom = np.array([TWO16 * 0.5 * (QSIG + SIG[e]) for e in els])

    lo = R.min(0) - (OUTER + 2.0); hi = R.max(0) + (OUTER + 2.0)
    xs = np.arange(lo[0], hi[0], SPACING); ys = np.arange(lo[1], hi[1], SPACING); zs = np.arange(lo[2], hi[2], SPACING)
    grid = np.array(np.meshgrid(xs, ys, zs)).reshape(3, -1).T

    keep_xyz, keep_md, keep_el = [], [], []
    B = 20000
    for s in range(0, len(grid), B):
        chunk = grid[s:s + B]
        d = np.linalg.norm(chunk[:, None, :] - R[None, :, :], axis=2)
        outside = (d - rmin_atom[None, :]).min(1) >= 0.0
        md = d.min(1); near = md <= OUTER
        nearest = d.argmin(1)
        m = outside & near
        for p, mdv, ne in zip(chunk[m], md[m], nearest[m]):
            keep_xyz.append(p); keep_md.append(mdv); keep_el.append(els[ne])
    keep_xyz = np.array(keep_xyz)

    with open(os.path.join(HERE, "grid_final.tsv"), "w") as f:
        f.write("idx\tx\ty\tz\tmin_dist\tnearest_el\n")
        for i, (p, md, e) in enumerate(zip(keep_xyz, keep_md, keep_el)):
            f.write("%d\t%.4f\t%.4f\t%.4f\t%.3f\t%s\n" % (i, p[0], p[1], p[2], md, e))
    with open(os.path.join(HERE, "grid_final.xyz"), "w") as f:
        f.write("%d\ncandidate grid (LJ-min standoff to OUTER=%.1f A)\n" % (len(keep_xyz), OUTER))
        for p in keep_xyz:
            f.write("X %.4f %.4f %.4f\n" % (p[0], p[1], p[2]))

    md = np.array(keep_md)
    print("substrate: 24 atoms | charge sigma %.2f A | outer bound %.1f A" % (QSIG, OUTER))
    print("accepted sites: %d" % len(keep_xyz))
    print("min-dist-to-substrate: min %.2f  median %.2f  max %.2f A" % (md.min(), np.median(md), md.max()))
    edges = np.arange(np.floor(md.min()), np.ceil(md.max()) + 1, 1.0)
    h, _ = np.histogram(md, bins=edges)
    print("min-dist histogram:", " ".join("%.0f-%.0fA:%d" % (edges[i], edges[i+1], h[i]) for i in range(len(h))))
    print("wrote grid_final.tsv + grid_final.xyz")


if __name__ == "__main__":
    main()
