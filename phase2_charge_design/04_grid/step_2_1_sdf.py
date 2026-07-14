#!/usr/bin/env python3
# Phase 2 step 2.1 - analytic union-of-spheres signed distance field (SDF) around the substrate.
# SDF(r) = min_atoms(|r - r_atom| - R_vdw(atom))  (Bondi radii): 0 on the vdW surface, <0 inside,
# >0 outside. Evaluated on a voxel grid enclosing the molecule with margin for the outer shells.
import sys
import numpy as np

BONDI = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52}

def read_xyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[0].split()[0])
    els, xyz = [], []
    for ln in lines[2:2+n]:
        p = ln.split()
        els.append(p[0]); xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return els, np.asarray(xyz, float)

def build_sdf(els, xyz, voxel=0.30, margin=5.0):
    radii = np.array([BONDI[e] for e in els])
    rmax = radii.max()
    lo = xyz.min(0) - (rmax + margin)
    hi = xyz.max(0) + (rmax + margin)
    nx, ny, nz = np.ceil((hi - lo) / voxel).astype(int) + 1
    gx = lo[0] + voxel*np.arange(nx)
    gy = lo[1] + voxel*np.arange(ny)
    gz = lo[2] + voxel*np.arange(nz)
    sdf = np.full((nx, ny, nz), np.inf, dtype=np.float32)
    GX, GY, GZ = np.meshgrid(gx, gy, gz, indexing="ij")
    for (ax, ay, az), r in zip(xyz, radii):
        d = np.sqrt((GX-ax)**2 + (GY-ay)**2 + (GZ-az)**2) - r
        np.minimum(sdf, d, out=sdf)
    return sdf, (lo, voxel, (nx, ny, nz)), radii

def verify(sdf, grid, els, xyz, radii):
    lo, voxel, dims = grid
    checks = {}
    centre_vals = []
    for (ax,ay,az) in xyz:
        i = int(round((ax-lo[0])/voxel)); j = int(round((ay-lo[1])/voxel)); k = int(round((az-lo[2])/voxel))
        centre_vals.append(float(sdf[i,j,k]))
    checks["max_sdf_at_atom_centres"] = max(centre_vals)
    checks["sdf_far_corner"] = float(sdf[0,0,0])
    checks["sdf_min"] = float(sdf.min()); checks["sdf_max"] = float(sdf.max())
    checks["has_zero_crossing"] = bool(sdf.min() < 0 < sdf.max())
    ix = int(np.argmax(xyz[:,0])); r = radii[ix]
    probe = xyz[ix] + np.array([r, 0, 0])
    i = int(round((probe[0]-lo[0])/voxel)); j = int(round((probe[1]-lo[1])/voxel)); k = int(round((probe[2]-lo[2])/voxel))
    checks["surface_probe_value_near_zero"] = float(sdf[i,j,k])
    return checks

if __name__ == "__main__":
    xyzpath = sys.argv[1] if len(sys.argv) > 1 else "reactant.xyz"
    voxel = float(sys.argv[2]) if len(sys.argv) > 2 else 0.30
    els, xyz = read_xyz(xyzpath)
    print(f"atoms: {len(els)}  formula elements: {sorted(set(els))}")
    sdf, grid, radii = build_sdf(els, xyz, voxel=voxel)
    lo, vox, dims = grid
    print(f"voxel: {vox} A   grid dims: {dims}   total voxels: {np.prod(dims):,}")
    print(f"grid origin (A): {lo.round(3)}")
    c = verify(sdf, grid, els, xyz, radii)
    print("\nverification:")
    for k,v in c.items():
        print(f"  {k}: {v}")
    np.savez_compressed("sdf_grid.npz", sdf=sdf, origin=lo, voxel=vox, dims=np.array(dims),
                        xyz=xyz, radii=radii, els=np.array(els))
    print("\nsaved sdf_grid.npz")
    ok = (c["max_sdf_at_atom_centres"] < 0 and c["sdf_far_corner"] > 0 and c["has_zero_crossing"]
          and abs(c["surface_probe_value_near_zero"]) < 2*vox)
    print("SDF VALID" if ok else "SDF CHECK FAILED")
