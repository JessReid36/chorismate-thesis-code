#!/usr/bin/env bash
# Build both candidate grids and compare them.
# Run from a directory containing reactant.xyz, ts.xyz, product.xyz.
set -euo pipefail

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1   # numpy needs this on hpc1

# Record the environment. scikit-image selects between two marching-cubes
# implementations by version, so the citation in the methods depends on this.
echo "############ 0. environment"
python3 - <<'VER'
import numpy, scipy, skimage, trimesh, sys
print("python    ", sys.version.split()[0])
print("numpy     ", numpy.__version__)
print("scipy     ", scipy.__version__)
print("scikit-image", skimage.__version__)
print("trimesh   ", trimesh.__version__)
VER
echo

SHELLS="3.0 4.0 5.0"
RMIN=1.0
DENSITY=6.0
SEED=0

echo "############ 1. signed distance field (union of R + TS + P)"
python3 01_build_sdf.py reactant.xyz ts.xyz product.xyz 0.30

echo
echo "############ 2. Poisson-disk grid (ours)"
python3 02_grid_poisson.py --r-min $RMIN --shells $SHELLS \
        --density $DENSITY --seed $SEED --out grid_poisson

# Match the CVT point count to the Poisson grid, shell by shell, so the
# comparison isolates the placement algorithm.
COUNTS=$(python3 -c "
import numpy as np
z=np.load('grid_poisson.npz'); s=z['shell']
print(' '.join(str(int((s==v).sum())) for v in [3.0,4.0,5.0]))")
echo
echo "matched per-shell counts for CVT: $COUNTS"

echo
echo "############ 3. CVT grid (colleague's method), global mode"
python3 03_grid_cvt.py --shells $SHELLS --counts $COUNTS \
        --density $DENSITY --seed $SEED --mode global --out grid_cvt_global

echo
echo "############ 4. CVT grid, per-shell mode (textbook restricted CVT)"
python3 03_grid_cvt.py --shells $SHELLS --counts $COUNTS \
        --density $DENSITY --seed $SEED --mode per_shell --out grid_cvt_per_shell

echo
echo "############ 5. comparison"
DV=""
[ -f dv_grid_v2.tsv ] && DV="--dv dv_grid_v2.tsv"
python3 04_compare_grids.py grid_poisson.npz grid_cvt_global.npz \
        grid_cvt_per_shell.npz $DV | tee grid_comparison.txt

echo
echo "done. outputs:"
ls -1 grid_*.xyz grid_*_sites.tsv grid_comparison.txt

echo
echo "############ 6. voxel convergence"
python3 05_convergence_test.py 2>&1 | tee convergence.txt

echo
echo "############ 7. envelope justification"
python3 06_envelope_test.py 2>&1 | tee envelope.txt
