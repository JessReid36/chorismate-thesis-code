#!/usr/bin/env python3
"""
check_top_site_degeneracy.py - reconcile two statements that appear to conflict.

THE APPARENT CONFLICT
    Test B6 reported that between some pairs of frames only three of the ten most
    stabilising sites are common. Taken alone that suggests the frames disagree
    profoundly about where to place a charge.

    The per-frame maps show the opposite: every frame's deepest basin lies within
    9 to 30 degrees of the ether oxygen, with a circular standard deviation in
    azimuth of 32 degrees, and the depth varying by 49 per cent.

    Both cannot be the whole story. The likely reconciliation is that the sites
    within the favourable region are nearly degenerate, so their identities
    reshuffle between frames while the region itself does not move. That is the
    same phenomenon the single-charge tests measured as a resolution limit: below
    a certain separation in Delta v, the objective cannot order candidates.

WHAT THIS CHECKS
    1. How far apart the ten most stabilising sites are, in Delta v, within each
       frame. If they differ by less than the across-frame spread, their ordering
       is not determined and reshuffling is expected.

    2. How far apart those sites are in space. If they cluster, the reshuffling
       moves the choice by a small distance; if they scatter, it does not.

    3. Whether the union of each frame's top ten is much larger than ten. A small
       union means the frames draw from the same pool in a different order; a
       large one means they genuinely look elsewhere.

    4. The same test at a looser cut, the top twenty and top fifty, to see
       whether agreement improves once near-degenerate sites are grouped rather
       than ranked.

Usage:
    python3 check_top_site_degeneracy.py [--shell 3.0]
"""
import argparse
import itertools
import math
import os

import numpy as np

H = 627.5094740631


def read_sites(path):
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    hdr = [h.strip() for h in lines[0].lstrip("#").split("\t")]
    col = {h: i for i, h in enumerate(hdr)}
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", default="grid_poisson_sites.tsv")
    ap.add_argument("--vpot-dir", default=".")
    ap.add_argument("--shell", type=float, default=0.0,
                    help="restrict to one shell; 0 uses the whole grid, which is "
                         "what test B6 used")
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
    if len(names) < 2:
        raise SystemExit("need at least two frames")

    shell = np.array([s[4] for s in sites])
    sel = np.ones(n, bool) if args.shell == 0 else (shell == args.shell)
    xyz = np.array([[s[1], s[2], s[3]] for s in sites])
    where = "the whole grid" if args.shell == 0 else f"the {args.shell:.0f} A shell"
    print(f"{len(names)} frames, {sel.sum()} sites on {where}\n")

    M = np.vstack([maps[f][sel] for f in names])
    pos = xyz[sel]
    mean = M.mean(0)
    sd = M.std(0, ddof=1)

    # ---------------------------------------- 1. separation within each frame
    print("=" * 74)
    print("How well separated are the top ten sites, within a single frame?")
    print("=" * 74)
    print(f"{'frame':<9}{'best':>9}{'10th':>9}{'gap':>9}"
          f"{'median step':>14}")
    print("-" * 74)
    gaps = []
    for k, f in enumerate(names):
        v = M[k]
        o = np.argsort(v)[:10]
        vals = v[o]
        steps = np.diff(vals)
        gaps.append(vals[-1] - vals[0])
        print(f"{f:<9}{vals[0]:>9.3f}{vals[-1]:>9.3f}{vals[-1]-vals[0]:>9.3f}"
              f"{np.median(steps):>14.3f}")
    sd_top = sd[np.argsort(mean)[:10]]
    print()
    print(f"  mean spread across the top ten within a frame: {np.mean(gaps):.3f} kcal/mol")
    print(f"  across-frame standard deviation at those sites: "
          f"{sd_top.min():.3f} to {sd_top.max():.3f}")
    if np.mean(gaps) < sd_top.mean() * 2:
        print()
        print("  The ten best sites within a frame span less than about twice the")
        print("  amount by which any one of them varies between frames. Their")
        print("  ordering is therefore not determined by the data, and reshuffling")
        print("  between frames is expected rather than evidence of disagreement.")
    else:
        print()
        print("  The ten best sites are separated by more than the across-frame")
        print("  variation, so their ordering is meaningful and reshuffling would")
        print("  indicate genuine disagreement.")

    # ---------------------------------------- 2. how far apart are they
    print()
    print("=" * 74)
    print("How far apart are those sites in space?")
    print("=" * 74)
    for k, f in enumerate(names[:3]):
        o = np.argsort(M[k])[:10]
        d = [np.linalg.norm(pos[a] - pos[b]) for a, b in itertools.combinations(o, 2)]
        print(f"  frame {f}: top ten span {min(d):.2f} to {max(d):.2f} A, "
              f"median {np.median(d):.2f}")
    om = np.argsort(mean)[:10]
    d = [np.linalg.norm(pos[a] - pos[b]) for a, b in itertools.combinations(om, 2)]
    print(f"  mean map: top ten span {min(d):.2f} to {max(d):.2f} A, "
          f"median {np.median(d):.2f}")

    # ---------------------------------------- 3. the union
    print()
    print("=" * 74)
    print("Do the frames draw from the same pool?")
    print("=" * 74)
    for k in (10, 20, 50):
        tops = [set(np.argsort(M[i])[:k].tolist()) for i in range(len(names))]
        union = set().union(*tops)
        inter = set.intersection(*tops)
        pair = [len(a & b) for a, b in itertools.combinations(tops, 2)]
        print(f"  top {k:<3}  union {len(union):<4} "
              f"({len(union)/k:.1f} times {k})   "
              f"common to all {len(inter):<3}   "
              f"pairwise overlap {min(pair)}/{k} to {max(pair)}/{k}, "
              f"median {int(np.median(pair))}/{k}")
    print()
    print("  A union near k means the frames rank the same sites differently.")
    print("  A union far above k means they look in different places.")

    # ---------------------------------------- 4. grouping instead of ranking
    print()
    print("=" * 74)
    print("Does agreement improve if near-degenerate sites are grouped?")
    print("=" * 74)
    thr = float(sd_top.mean())
    print(f"  grouping threshold: {thr:.3f} kcal/mol, the mean across-frame")
    print("  standard deviation at the ten most stabilising sites")
    print()
    agree = []
    for a, b in itertools.combinations(range(len(names)), 2):
        va, vb = M[a], M[b]
        ta = set(np.argsort(va)[:10].tolist())
        # sites in frame b within the threshold of b's tenth best
        cut = np.sort(vb)[9] + thr
        tb_loose = set(np.where(vb <= cut)[0].tolist())
        agree.append(len(ta & tb_loose))
    print(f"  frame A's top ten, found within {thr:.2f} kcal/mol of frame B's")
    print(f"  tenth best: {min(agree)}/10 to {max(agree)}/10, "
          f"median {int(np.median(agree))}/10")
    print()
    print("  If this is much higher than the strict top-ten overlap, the")
    print("  disagreement is about ordering within a degenerate set and not")
    print("  about which region of the surface is favourable.")


if __name__ == "__main__":
    main()
