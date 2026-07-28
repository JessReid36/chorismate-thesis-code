#!/usr/bin/env python3
# Rebuild the union-of-spheres SDF with a LARGER margin so outer shells (5-8 A) fit.
# Same construction as 04_grid/step_2_1_sdf.py (union envelope over R, TS, P) but margin 9.0.
# Needs R/TS/P xyz (from 01_geometry). Outputs sdf_grid_ext.npz.
import sys, numpy as np
VDW={"H":1.10,"C":1.70,"N":1.55,"O":1.52}  # standard vdW radii (A) - match original
def load(p):
    L=open(p).read().splitlines(); n=int(L[0].split()[0])
    els=[L[2+i].split()[0] for i in range(n)]
    xyz=np.array([[float(v) for v in L[2+i].split()[1:4]] for i in range(n)])
    return els,xyz
def build(els_list, xyz_list, voxel=0.30, margin=9.0):
    allx=np.vstack(xyz_list)
    radii_list=[np.array([VDW.get(e,1.7) for e in els]) for els in els_list]
    rmax=max(r.max() for r in radii_list)
    lo=allx.min(0)-(rmax+margin); hi=allx.max(0)+(rmax+margin)
    nx,ny,nz=[int(np.ceil((hi[i]-lo[i])/voxel))+1 for i in range(3)]
    gx=lo[0]+voxel*np.arange(nx); gy=lo[1]+voxel*np.arange(ny); gz=lo[2]+voxel*np.arange(nz)
    X,Y,Z=np.meshgrid(gx,gy,gz,indexing='ij')
    P=np.stack([X,Y,Z],-1).reshape(-1,3)
    sdf=np.full(len(P), 1e9)
    for xyz,radii in zip(xyz_list,radii_list):
        for a in range(len(xyz)):
            d=np.linalg.norm(P-xyz[a],axis=1)-radii[a]
            sdf=np.minimum(sdf,d)   # union of spheres = min over atoms
    sdf=sdf.reshape(nx,ny,nz)
    np.savez_compressed("sdf_grid_ext.npz", sdf=sdf, origin=lo, voxel=voxel,
                        xyz=allx, radii=np.concatenate(radii_list))
    print("SDF ext: shape %s, range [%.2f, %.2f] A, margin %.1f"%(sdf.shape,sdf.min(),sdf.max(),margin))
if __name__=="__main__":
    # args: reactant.xyz ts.xyz product.xyz
    els=[]; xyz=[]
    for p in sys.argv[1:4]:
        e,x=load(p); els.append(e); xyz.append(x)
    build(els,xyz)
