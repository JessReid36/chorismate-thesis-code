#!/bin/bash
# GOCAT Stage A: screen endpoints (R and P) for all intact designs. Relax under field + get RMSG/energy.
# Run from 05d_shell_ladder/relax. Produces shell<tag>_{reactant,product}_screen.txt per design.
set -e
INTACT="2p0 7p0 7p1 7p2 7p3 7p4 7p5 7p6 7p7 7p8 7p9 8p0 8p1 8p2 8p3 8p5 8p6 8p7 8p8 8p9 9p0 10p0 11p0 12p0 13p0 15p0"
cp ../../01_geometry/reactant.xyz ../../01_geometry/product.xyz . 2>/dev/null || true
for tag in $INTACT; do
  dtag="shell$tag"
  if [ ! -f "${dtag}_coords.xyz" ]; then
    [ -f "run_${dtag}/${dtag}_coords.xyz" ] && cp run_${dtag}/${dtag}_coords.xyz run_${dtag}/${dtag}_charges.txt . 2>/dev/null || { echo "SKIP $dtag"; continue; }
  fi
  for ep in reactant product; do
    mkdir -p run_screen_${tag}_${ep}
    cp gocat_screen.py ${dtag}_coords.xyz ${dtag}_charges.txt reactant.xyz product.xyz run_screen_${tag}_${ep}/
    { cat <<PBS
#!/bin/bash
#PBS -N scr_${tag}_${ep}
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -l walltime=48:00:00
#PBS -M 18660916@sun.ac.za -m ae
#PBS -o scr_${tag}_${ep}.stdout
#PBS -e scr_${tag}_${ep}.stderr
cd \$PBS_O_WORKDIR
source /apps/mambaforge/etc/profile.d/conda.sh
conda activate \$HOME/envs/ash
export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PATH=/apps/openmpi/4.1.1/bin:\$PATH
export LD_LIBRARY_PATH=/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:\$LD_LIBRARY_PATH
python3 gocat_screen.py $dtag $ep run > scr_${tag}_${ep}.log 2>&1
PBS
    } > run_screen_${tag}_${ep}/run.pbs
    ( cd run_screen_${tag}_${ep} && qsub run.pbs )
  done
  echo "submitted screen $tag (R+P)"
done
echo "Stage A screening submitted (2 jobs per design)."
