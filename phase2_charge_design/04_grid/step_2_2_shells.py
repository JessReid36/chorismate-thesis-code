#!/usr/bin/env python3
# Phase 2 step 2.2 - extract constant-distance shells as triangular meshes via marching cubes.
# Each shell is the isosurface SDF = d (additive, absolute offset in Angstrom).
import sys, numpy as np
from skimage import measure

def load_sdf(npz="sdf_grid.npz"):
    z = np.load(npz, allow_pickle=True)
    return z["sdf"], z["origin"].astype(float), float(z["voxel"]), z["xyz"].astype(float), z["radii"].astype(float)

def extract_shell(sdf, origin, voxel, level):
    verts, faces, normals, _ = measure.marching_cubes(sdf, level=level)
    verts_ang = origin + verts * voxel
    return verts_ang, faces, normals

def mesh_area(verts, faces):
    v = verts[faces]
    tri = np.cross(v[:,1]-v[:,0], v[:,2]-v[:,0])
    return 0.5*np.linalg.norm(tri, axis=1).sum()

def min_surface_distance(pts, xyz, radii):
    d = np.sqrt(((pts[:,None,:]-xyz[None,:,:])**2).sum(-1)) - radii[None,:]
    return d.min(1)

if __name__ == "__main__":
    shells = [float(x) for x in (sys.argv[1:] or ["2.0","3.0","4.0"])]
    sdf, origin, voxel, xyz, radii = load_sdf()
    print(f"loaded SDF grid {sdf.shape}, voxel {voxel} A, {len(xyz)} atoms")
    print(f"SDF range [{sdf.min():.3f}, {sdf.max():.3f}] A\n")
    for d in shells:
        if not (sdf.min() < d < sdf.max()):
            print(f"shell {d} A: SKIP (outside SDF range - increase margin in 2.1)"); continue
        verts, faces, normals = extract_shell(sdf, origin, voxel, d)
        area = mesh_area(verts, faces)
        samp = verts[np.random.default_rng(0).choice(len(verts), size=min(500,len(verts)), replace=False)]
        off = min_surface_distance(samp, xyz, radii)
        print(f"shell d={d:.1f} A: {len(verts):6d} verts, {len(faces):6d} faces, area {area:8.1f} A^2 | "
              f"measured offset mean {off.mean():.3f} +/- {off.std():.3f} A (target {d:.1f})")
        np.savez_compressed(f"shell_{d:.1f}.npz", verts=verts, faces=faces, normals=normals, level=d)
    print("\nsaved shell_*.npz meshes")
    print("\nEuler characteristic check (closed surface => 2 per component):")
    for d in shells:
        try:
            z = np.load(f"shell_{d:.1f}.npz")
            V, F = len(z["verts"]), len(z["faces"])
            E = F*3//2
            print(f"  shell {d:.1f}: V={V} E={E} F={F}  chi = {V-E+F}")
        except FileNotFoundError:
            pass
