#!/usr/bin/env python3
# Phase 2 step 2.4 - GLOBAL Poisson-disk / blue-noise thinning to a guaranteed min spacing.
# All shells' dense clouds are POOLED and thinned once in 3D, so the spacing guarantee holds
# across shells as well as within them. Grid-bucketed dart-elimination, pure numpy/scipy,
# deterministic given the seed. Per-point shell tags retained for later per-shell analysis.
import sys, itertools, numpy as np
from scipy.spatial import cKDTree

def thin_poisson_fast(pts, r_min, seed=0):
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(pts))
    origin = pts.min(0)
    def key(p): return tuple(((p - origin)//r_min).astype(int))
    buckets = {}; accepted = []; r2 = r_min*r_min
    for idx in order:
        p = pts[idx]; kx,ky,kz = key(p); ok = True
        for dx,dy,dz in itertools.product((-1,0,1), repeat=3):
            for q in buckets.get((kx+dx,ky+dy,kz+dz), ()):
                if ((pts[q]-p)**2).sum() < r2:
                    ok = False; break
            if not ok: break
        if ok:
            accepted.append(idx); buckets.setdefault((kx,ky,kz), []).append(idx)
    return np.array(accepted)

def coherence_report(pts):
    tree = cKDTree(pts); d, _ = tree.query(pts, k=2)
    min_sp = float(d[:,1].min())
    D = np.sqrt(((pts[:,None,:]-pts[None,:,:])**2).sum(-1))
    np.fill_diagonal(D, np.inf); C = 1.0/D
    ev = np.abs(np.linalg.eigvalsh(C)); ev = ev[ev > 1e-9]
    cond = float(ev.max()/ev.min()) if len(ev) else float('nan')
    return min_sp, cond, len(pts)

if __name__ == "__main__":
    r_min  = float(sys.argv[1]) if len(sys.argv) > 1 else 1.5
    shells = [float(x) for x in (sys.argv[2:] or ["2.0","3.0","4.0"])]
    pooled, tags = [], []
    for d in shells:
        try:
            pts = np.load(f"cloud_{d:.1f}.npz")["pts"]
        except FileNotFoundError:
            print(f"shell {d}: missing cloud_{d:.1f}.npz (run 2.3)"); continue
        pooled.append(pts); tags += [d]*len(pts)
    pooled = np.vstack(pooled); tags = np.array(tags)
    print(f"pooled dense cloud: {len(pooled)} points across {len(shells)} shells")
    keep = thin_poisson_fast(pooled, r_min, seed=0)
    grid = pooled[keep]; grid_tags = tags[keep]
    gmin, gcond, gn = coherence_report(grid)
    print("retained per shell (global thinning):")
    for d in shells:
        print(f"  shell {d:.1f} A: {(grid_tags==d).sum():4d} points")
    print(f"\nFINAL GRID: {gn} candidate points")
    print(f"  global min spacing {gmin:.3f} A (target >= {r_min})")
    print(f"  global Coulomb-matrix condition number {gcond:.1f}")
    np.savez_compressed("grid_final.npz", pts=grid, shell=grid_tags, r_min=r_min)
    with open("grid_final.xyz","w") as fh:
        fh.write(f"{len(grid)}\ncandidate charge grid r_min={r_min} shells={shells}\n")
        for (x,y,z),s in zip(grid, grid_tags):
            fh.write(f"X {x:.4f} {y:.4f} {z:.4f}   shell={s:.1f}\n")
    print("saved grid_final.npz + grid_final.xyz")
    print("GRID OK" if gmin >= r_min - 1e-6 else "GRID FAILED min-spacing")
