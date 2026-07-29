#!/bin/bash
# Build + submit TS-under-field barrier jobs for ALL shells (coarse 2-15A + fine 7.0-9.0 0.1A).
# Run from 05d_shell_ladder/relax. Needs each shell's <tag>_coords.xyz + <tag>_charges.txt + ts.xyz here,
# and the completed reactant ladder job in run_<tag>/ for E_react.
set -e
cp ../../01_geometry/ts.xyz . 2>/dev/null || true

# coarse whole-number shells 2-15
COARSE="2p0 3p0 4p0 5p0 6p0 7p0 8p0 9p0 10p0 11p0 12p0 13p0 14p0 15p0"
# fine sweep 7.0-9.0 in 0.1 (skip 7p0,8p0,9p0 already in COARSE)
FINE="7p1 7p2 7p3 7p4 7p5 7p6 7p7 7p8 7p9 8p1 8p2 8p3 8p4 8p5 8p6 8p7 8p8 8p9"

for tag in $COARSE $FINE; do
  dtag="shell$tag"
  if [ ! -f "${dtag}_coords.xyz" ]; then echo "SKIP $dtag (no design coords here)"; continue; fi
  # build TS-under-field design
  python3 make_ts_any.py $dtag >/dev/null
  # PBS
  cat > run_barr_${tag}.pbs <<PBS
#!/bin/bash
#PBS -N barr_${tag}
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -l walltime=48:00:00
#PBS -M 18660916@sun.ac.za -m ae
#PBS -o barr_${tag}.stdout
#PBS -e barr_${tag}.stderr
cd \$PBS_O_WORKDIR
source /apps/mambaforge/etc/profile.d/conda.sh
conda activate \$HOME/envs/ash
export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PATH=/apps/openmpi/4.1.1/bin:\$PATH
export LD_LIBRARY_PATH=/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:\$LD_LIBRARY_PATH
python3 relax_barrier.py ${dtag}ts run > barr_${tag}_relax.log 2>&1
PBS
  # own dir, submit
  mkdir -p run_barr_$tag
  cp relax_barrier.py ${dtag}ts_coords.xyz ${dtag}ts_charges.txt run_barr_${tag}.pbs run_barr_$tag/
  ( cd run_barr_$tag && qsub run_barr_${tag}.pbs )
  echo "submitted barr_$tag"
done
echo "ALL barrier jobs submitted."
