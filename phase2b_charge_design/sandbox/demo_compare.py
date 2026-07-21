#!/usr/bin/env python3
# Reproduces the sandbox findings on synthetic fields (magnitudes illustrative; the STRUCTURAL
# results are properties of the math). Swap in real loaders when the consistent-region grid exists.
from synth_field import single_image_grid, path_field
from formulations import certified_milp, distributed_minmax, wholepath_minmax, true_path_barrier_from_singlepoint
BARE=18.549
print("=== certified max-lowering MILP (single-image) — K>=3 drives barrier negative (the pathology) ===")
grid=single_image_grid()
print("%3s %14s %14s %8s"%("K","ddE(kcal/mol)","barrier","gap"))
for K in (2,4,6,8):
    ddE,g=certified_milp(grid,K); print("%3d %14.3f %14.3f %8.4f"%(K,ddE,BARE+ddE,g))
print("\n=== distributed min-max (total lowering held in enzyme-target band [-13.5,-9.1]) — max|contrib| falls as charges share the load ===")
print("%3s %16s %8s"%("K","max|contrib|(kcal)","gap"))
for K in (2,4,6,8):
    mx,g=distributed_minmax(grid,K,band=(-(BARE-5.0), -9.1)); print("%3d %16.3f %8.4f"%(K,mx,g))
print("\n=== single-point (deceptive) vs whole-path (honest, bounded, certified) ===")
pts,b,dv=path_field()
print("bare path barrier: +%.1f kcal/mol"%max(b))
print("%3s | %12s %14s | %12s %8s"%("K","TS-img val","TRUE barrier","whole-path","gap"))
for K in (2,4,6,8):
    tsv,tb=true_path_barrier_from_singlepoint(pts,b,dv,K); wp,g=wholepath_minmax(pts,b,dv,K)
    print("%3d | %12.2f %14.2f | %12.2f %8.4f"%(K,tsv,tb,wp,g))
