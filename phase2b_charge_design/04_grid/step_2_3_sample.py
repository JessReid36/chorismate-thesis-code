#!/usr/bin/env python3
# Phase 2 step 2.3 - area-proportional sampling of each shell mesh into a dense point cloud.
import sys, numpy as np, trimesh
from scipy.spatial import cKDTree

def load_shell(d):
    z = np.load(f"shell_{d:.1f}.npz")
    return z["verts"], z["faces"], float(z["level"])

def sample_shell(verts, faces, density):
    mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
    area = mesh.area
    n = max(4, int(round(area * density)))
    pts, _ = trimesh.sample.sample_surface(mesh, n, seed=0)
    return np.asarray(pts), area, n

if __name__ == "__main__":
    shells = [float(x) for x in (sys.argv[1:] or ["2.0","3.0","4.0"])]
    density = 2.0  # dense pre-thinning samples per A^2
    allpts = {}
    for d in shells:
        try:
            verts, faces, lvl = load_shell(d)
        except FileNotFoundError:
            print(f"shell {d}: missing shell_{d:.1f}.npz (run 2.2)"); continue
        pts, area, n = sample_shell(verts, faces, density)
        dists, _ = cKDTree(pts).query(pts, k=2); nn = dists[:,1]
        print(f"shell {d:.1f} A: area {area:7.1f} A^2 -> {n:5d} samples "
              f"(density {density}/A^2) | NN spacing mean {nn.mean():.3f} min {nn.min():.3f} A")
        np.savez_compressed(f"cloud_{d:.1f}.npz", pts=pts, level=d, area=area)
        allpts[d] = pts
    print(f"\ntotal dense samples: {sum(len(p) for p in allpts.values())}")
    print("saved cloud_*.npz")
