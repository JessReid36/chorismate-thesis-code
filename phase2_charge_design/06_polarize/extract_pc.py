#!/usr/bin/env python3
# Extract a named discrete design from sol_discrete.npz and write an ORCA .pc point-charge file.
import sys, numpy as np
npz, design, out = sys.argv[1], sys.argv[2], sys.argv[3]
z = np.load(npz, allow_pickle=True)
q = z[f"{design}_q"]; xyz = z["xyz"]; idx = z["idx"]
sel = np.where(q != 0)[0]
with open(out, "w") as fh:
    fh.write(f"{len(sel)}\n")
    for s in sel:
        x, y, z_ = xyz[s]
        fh.write(f"{int(q[s]):+d}   {x:.6f}   {y:.6f}   {z_:.6f}\n")
print(f"design {design}: {len(sel)} charges -> {out}")
for s in sel:
    print(f"  grid idx {int(idx[s])}: q={int(q[s]):+d}  xyz=({xyz[s][0]:.3f},{xyz[s][1]:.3f},{xyz[s][2]:.3f})")
print(f"net charge = {int(q[sel].sum()):+d}")
