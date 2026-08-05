#!/usr/bin/env python3
# Sanity gate: build a 6-site FRACTIONAL, net-neutral charge design at a SAFE (intact) distance
# and feed it to the PINNED gocat_screen.py reactant relaxation. Tests ONLY whether fractional
# point charges behind the LJ wall run cleanly under CPCM (SCF converges, ether intact, finite RMSG)
# -- the go/no-go for the abstract-fractional-charge representation.
# Pure Python (no numpy) so it runs on the login node.
import os

HERE = os.path.dirname(os.path.abspath(__file__))
INP  = os.path.join(HERE, "..", "inputs")
GRID = os.path.join(INP, "dv_grid_ext.tsv")   # shells 2-8 A (full inner+mid grid with V_R/V_TS/Dv)
NSEL = 6
QMAG = 0.5                                     # fractional magnitude; net-neutral (3x -0.5, 3x +0.5)
SHELL_LO, SHELL_HI = 7.0, 8.0                  # safe intact band (7-15 A stayed intact in the ladder)

def dist(a, b):
    return sum((a[k]-b[k])**2 for k in range(3))**0.5

# --- substrate (24 atoms) verbatim from reactant.xyz ---
rlines = open(os.path.join(INP, "reactant.xyz")).read().splitlines()
sub = rlines[2:26]
assert len(sub) == 24, "reactant.xyz must have 24 substrate atoms"
sub_xyz = [tuple(float(v) for v in l.split()[1:4]) for l in sub]

# --- grid sites (idx x y z shell V_R V_TS Dv) ---
sites = []
for ln in open(GRID).read().splitlines()[1:]:
    p = ln.split("\t")
    if len(p) < 8:
        continue
    sites.append({"idx": int(p[0]), "x": float(p[1]), "y": float(p[2]), "z": float(p[3]),
                  "shell": float(p[4]), "dv": float(p[7])})

# --- select the safe band, then the NSEL sites of largest |Dv| ---
band = [s for s in sites if SHELL_LO <= s["shell"] <= SHELL_HI]
if len(band) < NSEL:
    band = sorted(sites, key=lambda s: abs(s["shell"] - 7.5))[:40]
band.sort(key=lambda s: -abs(s["dv"]))
sel = band[:NSEL]

# --- net-neutral fractional charges: -QMAG to the highest-Dv half, +QMAG to the lowest ---
sel.sort(key=lambda s: -s["dv"])
for i, s in enumerate(sel):
    s["q"] = -QMAG if i < NSEL // 2 else +QMAG

# --- write coords (24 substrate placeholder + 6 C quasi-atoms) and charges ---
n = 24 + NSEL
with open(os.path.join(HERE, "sanity6_coords.xyz"), "w") as f:
    f.write("%d\nsanity gate: 24 substrate (from reactant.xyz) + %d fractional C quasi-atoms\n" % (n, NSEL))
    for l in sub:
        f.write(l.rstrip() + "\n")
    for s in sel:
        f.write("C   %14.8f %14.8f %14.8f\n" % (s["x"], s["y"], s["z"]))
with open(os.path.join(HERE, "sanity6_charges.txt"), "w") as f:
    f.write(",".join("%.4f" % s["q"] for s in sel))

# --- report ---
shells = [s["shell"] for s in sel]
print("selected %d sites | shells %.1f-%.1f A | net charge %+.3f" %
      (NSEL, min(shells), max(shells), sum(s["q"] for s in sel)))
print("idx     shell    Dv          q      min_dist_to_substrate(A)")
for s in sel:
    md = min(dist((s["x"], s["y"], s["z"]), a) for a in sub_xyz)
    print("%-6d  %5.1f  %+10.6f  %+.2f   %6.2f" % (s["idx"], s["shell"], s["dv"], s["q"], md))
print("\nwrote sanity6_coords.xyz (%d atoms) and sanity6_charges.txt" % n)
