#!/usr/bin/env python3
# Add a measured standoff column to the Dv grid.
#
# The sweep so far restricted charge placement by the grid's "shell" label, which is a
# construction artefact rather than a measured quantity. The physical constraint is distance from
# the substrate, so this computes, for every grid point, the minimum distance to any substrate
# atom across the pooled reactant + TS + product geometries -- the same union the grid's original
# 1.7 A floor was applied against, so the two are on one footing.
#
# It also re-measures the active-site charged residues on that same basis. Those coordinates are
# TS-frame, so they are additionally reported against the TS geometry alone, which is their native
# frame; the earlier reactant-frame figures were approximate.

import os
import numpy as np

GRID = "dv_grid.tsv"
OUT = "dv_grid_dist.tsv"
COL_X, COL_SHELL, COL_DV = 1, 4, 7
N_SUBSTRATE = 24                  # CHA

# charge-centre coordinates, TS frame, as used in the field visualisations
RESIDUES = {
    "Arg90":  (22.022, 35.773, 55.567),
    "Arg7":   (22.749, 39.245, 51.757),
    "Glu78":  (25.518, 34.447, 56.683),
    "Lys60'": (18.868, 42.675, 61.459),
    "Arg63'": (18.797, 45.286, 58.422),
}

SEARCH_DIRS = [
    ".", "..", "../01_substrate", "../..",
    os.path.expanduser("~/Desktop/phase2_local"),
    os.path.expanduser("~/Desktop/phase2_local/01_substrate"),
    os.path.expanduser("~/Desktop/hpc1_home_18660916/attempts/attempt_3"),
    os.path.expanduser("~/system_dev_offline"),
]
WANTED = ("reactant.xyz", "ts.xyz", "product.xyz")


def find_geometries():
    found = {}
    for d in SEARCH_DIRS:
        if not os.path.isdir(d):
            continue
        for root, _dirs, files in os.walk(d):
            for w in WANTED:
                if w in files and w not in found:
                    found[w] = os.path.join(root, w)
            if len(found) == len(WANTED):
                return found
    return found


def read_xyz(path, nmax):
    lines = open(path).read().splitlines()
    coords = []
    for ln in lines[2:]:
        parts = ln.split()
        if len(parts) < 4:
            continue
        try:
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        except ValueError:
            continue
        if len(coords) == nmax:
            break
    return np.array(coords)


def min_dist(points, atoms):
    out = np.empty(len(points))
    for i, p in enumerate(points):
        out[i] = np.sqrt(np.sum((atoms - p) ** 2, axis=1)).min()
    return out


def main():
    geo = find_geometries()
    for w in WANTED:
        print("%-14s %s" % (w, geo.get(w, "NOT FOUND")))
    if len(geo) < len(WANTED):
        print("\nmissing geometries; scp them from HPC 01_substrate before continuing")
        return

    frames = {}
    for w in WANTED:
        frames[w] = read_xyz(geo[w], N_SUBSTRATE)
        print("  %-14s %d atoms" % (w, len(frames[w])))
    pooled = np.vstack([frames[w] for w in WANTED])

    xyz, shell, dv = [], [], []
    for line in open(GRID):
        parts = line.split()
        if len(parts) <= COL_DV:
            continue
        try:
            xyz.append([float(parts[COL_X]), float(parts[COL_X + 1]), float(parts[COL_X + 2])])
            shell.append(int(float(parts[COL_SHELL])))
            dv.append(float(parts[COL_DV]))
        except ValueError:
            continue
    xyz = np.array(xyz); shell = np.array(shell); dv = np.array(dv)

    dist = min_dist(xyz, pooled)
    print("\ngrid standoff (min distance to pooled R+TS+P substrate):")
    print("  range %.2f - %.2f A, mean %.2f" % (dist.min(), dist.max(), dist.mean()))
    print("  grid floor was 1.7 A; points below that: %d" % int((dist < 1.7).sum()))

    print("\nshell label vs measured standoff:")
    for s in sorted(set(shell)):
        sel = shell == s
        print("  shell %d  %3d points   %.2f - %.2f A   mean %.2f   max|Dv| %.6f"
              % (s, sel.sum(), dist[sel].min(), dist[sel].max(), dist[sel].mean(),
                 np.abs(dv[sel]).max()))

    print("\nresidue standoff, measured two ways:")
    print("  residue   vs TS frame   vs pooled R+TS+P")
    ts_only = frames["ts.xyz"]
    band = []
    for n, p in RESIDUES.items():
        p = np.array(p)
        d_ts = np.sqrt(np.sum((ts_only - p) ** 2, axis=1)).min()
        d_pool = np.sqrt(np.sum((pooled - p) ** 2, axis=1)).min()
        band.append(d_pool)
        print("  %-7s   %5.2f         %5.2f" % (n, d_ts, d_pool))
    band = np.array(band)
    print("\n  enzyme band (pooled): %.2f - %.2f A" % (band.min(), band.max()))
    inband = (dist >= band.min()) & (dist <= band.max())
    print("  grid points inside that band: %d of %d" % (int(inband.sum()), len(dist)))
    if inband.sum():
        print("  their Dv range: %+.6f to %+.6f, max|Dv| %.6f"
              % (dv[inband].min(), dv[inband].max(), np.abs(dv[inband]).max()))

    fh = open(OUT, "w")
    fh.write("idx\tx_ang\ty_ang\tz_ang\tshell\tDv\tstandoff\n")
    for i in range(len(dv)):
        fh.write("%d\t%.4f\t%.4f\t%.4f\t%d\t%.8f\t%.4f\n"
                 % (i, xyz[i][0], xyz[i][1], xyz[i][2], shell[i], dv[i], dist[i]))
    fh.close()
    print("\nwrote %s" % OUT)


if __name__ == "__main__":
    main()
