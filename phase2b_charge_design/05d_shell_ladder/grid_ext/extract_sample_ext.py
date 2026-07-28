#!/usr/bin/env python3
# Extract outer shells from the extended SDF + sample+thin to candidate points. Mirrors 04_grid steps 2.2-2.4.
import sys, numpy as np
from skimage import measure
def load_sdf(npz="sdf_grid_ext.npz"):
    z=np.load(npz,allow_pickle=True)
    return z["sdf"], z["origin"].astype(float), float(z["voxel"]), z["xyz"].astype(float), z["radii"].astype(float)
def min_surf_dist(pts,xyz,radii):
    d=np.sqrt(((pts[:,None,:]-xyz[None,:,:])**2).sum(-1))-radii[None,:]
    return d.min(1)
def thin(pts, r=1.5):
    kept=[]
    for p in pts:
        if all(np.linalg.norm(p-q)>=r for q in kept): kept.append(p)
    return np.array(kept)
if __name__=="__main__":
    shells=[float(x) for x in (sys.argv[1:] or ["5.0","6.0","7.0","8.0"])]
    sdf,origin,voxel,xyz,radii=load_sdf()
    print("SDF range [%.2f,%.2f]"%(sdf.min(),sdf.max()))
    allpts=[]; allsh=[]
    rng=np.random.default_rng(0)
    for d in shells:
        if not (sdf.min()<d<sdf.max()):
            print("shell %.1f SKIP (outside SDF)"%d); continue
        verts,faces,normals,_=measure.marching_cubes(sdf,level=d)
        va=origin+verts*voxel
        # sample down then thin to ~1.5A spacing (match original grid density)
        idx=rng.choice(len(va),size=min(3000,len(va)),replace=False)
        kept=thin(va[idx],1.5)
        allpts.append(kept); allsh+=[d]*len(kept)
        off=min_surf_dist(kept,xyz,radii)
        print("shell %.1f A: %d candidate pts, offset %.2f+/-%.2f"%(d,len(kept),off.mean(),off.std()))
    P=np.vstack(allpts)
    with open("outer_shell_points.tsv","w") as fh:
        fh.write("idx\tx_ang\ty_ang\tz_ang\tshell\n")
        for i,(p,s) in enumerate(zip(P,allsh)):
            fh.write("%d\t%.4f\t%.4f\t%.4f\t%.1f\n"%(1000+i,p[0],p[1],p[2],s))  # idx 1000+ to not collide with existing
    print("wrote outer_shell_points.tsv (%d points, idx 1000+)"%len(P))
