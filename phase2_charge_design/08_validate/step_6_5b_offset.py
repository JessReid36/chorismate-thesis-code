#!/usr/bin/env python3
# Phase 2 step 6.5b - residue electrostatic-optimality OFFSET map. For each active-site residue, find
# the nearest grid point whose Dv is of the sign that HELPS at that residue (min Dv where the field
# wants a cation; max Dv where it wants an anion), within a local search radius. Report offset distance
# + direction: how far, and which way, toward the local electrostatic optimum on the sampled envelope.
# COARSE: the grid samples only the outer envelope, so offsets are RELATIVE/DIRECTIONAL, not precise.
import sys, numpy as np

def load_grid_dv(tsv):
    idx,xyz,dv=[],[],[]
    for ln in open(tsv).read().splitlines()[1:]:
        p=ln.split("\t")
        idx.append(int(p[0])); xyz.append([float(p[1]),float(p[2]),float(p[3])]); dv.append(float(p[7]))
    return np.array(idx),np.array(xyz),np.array(dv)

def residue_offset(res_pos, wants_sign, gxyz, gdv, search_radius=6.0):
    d = np.linalg.norm(gxyz - res_pos, axis=1)
    local = d <= search_radius
    if local.sum()==0: return None
    lx, ld = gxyz[local], gdv[local]
    best = np.argmin(ld) if wants_sign>0 else np.argmax(ld)
    opt = lx[best]
    vec = opt - res_pos
    return dict(opt=opt, vec=vec, dist=np.linalg.norm(vec), dv_at_opt=ld[best],
                nearest_grid_d=d[local].min())

if __name__ == "__main__":
    grid_tsv = sys.argv[1] if len(sys.argv)>1 else "dv_grid.tsv"
    res_npz  = sys.argv[2] if len(sys.argv)>2 else "residue_dv.npz"
    sr = float(sys.argv[3]) if len(sys.argv)>3 else 6.0
    gidx,gxyz,gdv = load_grid_dv(grid_tsv)
    z=np.load(res_npz, allow_pickle=True)
    names=z["names"]; signs=z["signs"]; centers=z["centers"]; dv=z["dv"]
    # Szefczyk DTSS for context
    szef={"Arg90":-9.06,"Arg7":-5.90,"Glu78":-3.57,"Arg116":-2.45,"Arg63":-1.40,"Lys60":+1.39}
    print(f"grid Dv points: {len(gdv)} (Dv range {gdv.min():+.5f} to {gdv.max():+.5f}) | "
          f"residues: {len(names)} | search radius {sr} A")
    print("COARSE offset map: grid samples the outer envelope, so magnitudes are RELATIVE/DIRECTIONAL.\n")
    print(f"{'residue':>8} {'Dv@res':>10} {'wants':>7} {'offset':>7} {'Dv@opt':>10} {'near-grid':>9}  {'direction (dx,dy,dz)':>22}")
    for nm,s,c,d in zip(names,signs,centers,dv):
        key=str(nm)
        wants = +1 if d<0 else -1   # what the field wants THERE
        r = residue_offset(np.array(c), wants, gxyz, gdv, sr)
        wsym = "cation" if wants>0 else "anion"
        if r is None:
            print(f"{key:>8} {d:+10.6f} {wsym:>7}   (no grid points within {sr} A)")
            continue
        v=r["vec"]
        print(f"{key:>8} {d:+10.6f} {wsym:>7} {r['dist']:7.2f} {r['dv_at_opt']:+10.6f} {r['nearest_grid_d']:9.2f}  "
              f"({v[0]:+.2f},{v[1]:+.2f},{v[2]:+.2f})")
    print("\nInterpretation (relative): smaller offset + smaller near-grid distance = residue sits closer to")
    print("the favorable-field region (well-placed). Larger = its local optimum is displaced (investigate).")
    print("Note: 'near-grid' = distance to closest grid point at all; large values mean the residue sits")
    print("far from the sampled envelope, so its offset is less reliable (grid doesn't sample near it).")
