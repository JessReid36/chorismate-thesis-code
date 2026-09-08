#!/usr/bin/env python3
"""
s3_dv_on_grid_v2.py - evaluate the difference potential Dv = V_TS - V_R at
every site of the candidate grid.

Two modes, called in sequence by s3_dv_on_grid_v2.pbs:

  prepare <sites.tsv> <out_bohr.xyz> <out_order.tsv>
      Read the grid, convert Angstrom -> Bohr, and write the points file for
      orca_vpot. Row order is preserved and recorded separately, because
      orca_vpot returns values positionally with no site identifiers - if the
      order were lost the entire map would be silently scrambled.

  assemble <order.tsv> <vR.out> <vTS.out> <out_dv.tsv> [--compare old_dv.tsv]
      Difference the two potentials site by site, write the objective vector,
      and report the statistics needed for the methods write-up.

Sign convention (validated separately by the sign check): Dv < 0 means a
positive charge at that site stabilises the transition state relative to the
reactant, and so lowers the barrier.

No new SCF is performed anywhere. The reactant and TS densities are the
committed single points, and orca_vpot only evaluates their potentials at new
positions, so the resulting map is directly comparable to the barrier those
same densities define.
"""

import argparse
import sys
from pathlib import Path

ANG2BOHR = 1.8897259886
HARTREE2KCAL = 627.5094740631


def read_sites(path):
    """Accept either the *_sites.tsv written by the grid builder or a plain
    grid .xyz. Returns (idx, x, y, z, shell, standoff)."""
    lines = Path(path).read_text().splitlines()
    rows = []

    if lines[0].startswith("idx\t") or lines[0].split("\t")[0] == "idx":
        header = lines[0].split("\t")
        ci = {name: header.index(name) for name in ("x", "y", "z", "shell")}
        so_i = header.index("standoff") if "standoff" in header else None
        for n, line in enumerate(lines[1:]):
            if not line.strip():
                continue
            p = line.split("\t")
            rows.append((n,
                         float(p[ci["x"]]), float(p[ci["y"]]), float(p[ci["z"]]),
                         float(p[ci["shell"]]),
                         float(p[so_i]) if so_i is not None else float("nan")))
    else:
        n_expected = int(lines[0].split()[0])
        for n, line in enumerate(lines[2:2 + n_expected]):
            p = line.split()
            shell = float(line.split("shell=")[1].split()[0]) if "shell=" in line else float("nan")
            so = float(line.split("standoff=")[1].split()[0]) if "standoff=" in line else float("nan")
            rows.append((n, float(p[1]), float(p[2]), float(p[3]), shell, so))
        if len(rows) != n_expected:
            sys.exit(f"FAIL: header says {n_expected} sites, read {len(rows)}")

    if not rows:
        sys.exit(f"FAIL: no sites read from {path}")
    return rows


def prepare(sites_path, out_bohr, out_order):
    rows = read_sites(sites_path)

    from collections import Counter
    per_shell = Counter(r[4] for r in rows)
    print(f"read {len(rows)} sites from {sites_path}")
    for s in sorted(per_shell):
        print(f"  shell {s:.1f} A: {per_shell[s]} sites")

    so = [r[5] for r in rows if r[5] == r[5]]          # drop NaN
    if so:
        print(f"  standoff {min(so):.3f} - {max(so):.3f} A")
        if min(so) < 2.5:
            print(f"  NOTE: closest site is {min(so):.3f} A off the vdW surface. "
                  f"Molecular surrogates need roughly 3 A; this grid may be "
                  f"point-charge only.")

    with open(out_bohr, "w") as fh:
        fh.write(f"{len(rows)}\n")
        for _, x, y, z, _, _ in rows:
            fh.write(f"{x*ANG2BOHR:.10f}  {y*ANG2BOHR:.10f}  {z*ANG2BOHR:.10f}\n")

    with open(out_order, "w") as fh:
        fh.write("row\tidx\tx\ty\tz\tshell\tstandoff\n")
        for n, (i, x, y, z, s, o) in enumerate(rows):
            fh.write(f"{n}\t{i}\t{x:.4f}\t{y:.4f}\t{z:.4f}\t{s:.1f}\t{o:.3f}\n")

    print(f"\nwrote {len(rows)} points to {out_bohr} (Bohr) "
          f"and the row order to {out_order}")


def read_vpot(path):
    vals = []
    for line in Path(path).read_text().splitlines():
        p = line.split()
        if len(p) >= 4:
            try:
                vals.append(float(p[-1]))
            except ValueError:
                pass
    return vals


def assemble(order_path, vR_path, vTS_path, out_path, compare=None):
    order = [l.split("\t") for l in Path(order_path).read_text().splitlines()[1:] if l.strip()]
    vR, vTS = read_vpot(vR_path), read_vpot(vTS_path)

    if not (len(vR) == len(vTS) == len(order)):
        sys.exit(f"FAIL: length mismatch - order={len(order)} "
                 f"vR={len(vR)} vTS={len(vTS)}. The row order and the potentials "
                 f"must correspond one-to-one; do not proceed.")

    recs = []
    for row, a, b in zip(order, vR, vTS):
        idx, x, y, z, shell, so = row[1], row[2], row[3], row[4], float(row[5]), row[6]
        recs.append((idx, x, y, z, shell, so, a, b, b - a))

    with open(out_path, "w") as fh:
        fh.write("idx\tx\ty\tz\tshell\tstandoff\tV_R_Eh\tV_TS_Eh\tdv_Eh\tdv_kcal_per_e\n")
        for idx, x, y, z, shell, so, a, b, d in recs:
            fh.write(f"{idx}\t{x}\t{y}\t{z}\t{shell:.1f}\t{so}\t"
                     f"{a:.9f}\t{b:.9f}\t{d:.9f}\t{d*HARTREE2KCAL:.4f}\n")

    dv = [r[8] for r in recs]
    shells = sorted({r[4] for r in recs})
    n = len(dv)
    mean = sum(dv) / n
    var = sum((d - mean) ** 2 for d in dv) / n
    stab = sum(1 for d in dv if d < 0)

    print(f"\nDv MAP: {n} sites")
    print(f"  range   {min(dv):+.6f} to {max(dv):+.6f} Eh")
    print(f"          {min(dv)*HARTREE2KCAL:+.3f} to {max(dv)*HARTREE2KCAL:+.3f} kcal/mol per +1e")
    print(f"  mean    {mean:+.6f} Eh   sd {var**0.5:.6f} Eh")
    print(f"  stabilising (Dv < 0): {stab}/{n} = {100*stab/n:.1f}%")

    print(f"\n{'shell':>8}{'sites':>8}{'min Dv':>12}{'max Dv':>12}"
          f"{'max |Dv| kcal':>16}{'stabilising':>14}")
    print("-" * 70)
    for s in shells:
        sub = [r[8] for r in recs if r[4] == s]
        ns = sum(1 for d in sub if d < 0)
        print(f"{s:>8.1f}{len(sub):>8}{min(sub):>12.6f}{max(sub):>12.6f}"
              f"{max(abs(min(sub)), abs(max(sub)))*HARTREE2KCAL:>16.3f}"
              f"{f'{ns} ({100*ns/len(sub):.1f}%)':>14}")

    best = min(recs, key=lambda r: r[8])
    worst = max(recs, key=lambda r: r[8])
    print(f"\n  most stabilising: idx {best[0]}, shell {best[4]:.1f} A, "
          f"({best[1]}, {best[2]}, {best[3]}), "
          f"{best[8]:+.6f} Eh = {best[8]*HARTREE2KCAL:+.3f} kcal/mol per +1e")
    print(f"  most destabilising: idx {worst[0]}, shell {worst[4]:.1f} A, "
          f"{worst[8]*HARTREE2KCAL:+.3f} kcal/mol per +1e")
    print(f"\n  best achievable single charge: "
          f"{min(min(dv), -max(dv))*HARTREE2KCAL:+.3f} kcal/mol")

    # local gradient - sets the placement-error penalty quoted in the methods
    try:
        import numpy as np
        from scipy.spatial import cKDTree
        P = np.array([[float(r[1]), float(r[2]), float(r[3])] for r in recs])
        D = np.array(dv) * HARTREE2KCAL
        d, j = cKDTree(P).query(P, k=2)
        g = np.abs(D - D[j[:, 1]]) / d[:, 1]
        print(f"\n  local Dv gradient: mean {g.mean():.3f}, p95 "
              f"{np.percentile(g,95):.3f} kcal/mol/A per +1e")
        print(f"  (multiply by the grid placement error to get the positioning penalty)")
    except ImportError:
        pass

    if compare:
        try:
            old = [l.split("\t") for l in Path(compare).read_text().splitlines()[1:] if l.strip()]
            odv = [float(r[7]) for r in old]
            ob = min(min(odv), -max(odv)) * HARTREE2KCAL
            nb = min(min(dv), -max(dv)) * HARTREE2KCAL
            print(f"\n  previous grid ({len(odv)} sites): best single charge {ob:+.3f} kcal/mol")
            print(f"  this grid     ({n} sites): best single charge {nb:+.3f} kcal/mol")
            print(f"  difference {nb-ob:+.3f} kcal/mol "
                  f"({100*(abs(nb)-abs(ob))/abs(ob):+.0f}% of the achievable single-charge effect)")
            print(f"  This is the cost of the shell change, and should be reported as such.")
        except Exception as e:
            print(f"\n  (comparison against {compare} failed: {e})")

    print(f"\nwrote {out_path}")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    mode = sys.argv[1]
    if mode == "prepare":
        prepare(*sys.argv[2:5])
    elif mode == "assemble":
        ap = argparse.ArgumentParser()
        ap.add_argument("order")
        ap.add_argument("vR")
        ap.add_argument("vTS")
        ap.add_argument("out")
        ap.add_argument("--compare", default=None)
        a = ap.parse_args(sys.argv[2:])
        assemble(a.order, a.vR, a.vTS, a.out, a.compare)
    else:
        sys.exit(f"unknown mode: {mode}")


if __name__ == "__main__":
    main()
