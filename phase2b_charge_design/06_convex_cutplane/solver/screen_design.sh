#!/bin/bash
# Screen one design through the pinned oracle: relax reactant AND product under its field.
# Usage (on the HPC):  bash screen_design.sh <name> <explore_tag>
#   e.g.  bash screen_design.sh s3    s3        (control: close, expect fragment)
#         bash screen_design.sh m4    m4_ge7    (intact-band test)
# Expects explore_<tag>_coords.xyz / _charges.txt in this dir (scp'd from the PC).
set -e
NAME="$1"; TAG="$2"
[ -n "$NAME" ] && [ -n "$TAG" ] || { echo "usage: bash screen_design.sh <name> <explore_tag>"; exit 1; }
[ -f "explore_${TAG}_coords.xyz" ] || { echo "missing explore_${TAG}_coords.xyz here"; exit 1; }

for EP in reactant product; do
  RUN="screen_${NAME}_${EP}"
  mkdir -p "$RUN"
  cp ../oracle/gocat_screen.py "$RUN"/
  cp ../inputs/reactant.xyz ../inputs/ts.xyz ../inputs/product.xyz "$RUN"/
  cp "explore_${TAG}_coords.xyz"  "$RUN/${NAME}_coords.xyz"
  cp "explore_${TAG}_charges.txt" "$RUN/${NAME}_charges.txt"
  {
    printf '%s\n' '#!/bin/bash'
    printf '%s\n' "#PBS -N scr_${NAME}_${EP}"
    printf '%s\n' '#PBS -l select=1:ncpus=8:mem=32gb'
    printf '%s\n' '#PBS -l walltime=24:00:00'
    printf '%s\n' '#PBS -M 18660916@sun.ac.za -m ae'
    printf '%s\n' "#PBS -o scr_${NAME}_${EP}.stdout"
    printf '%s\n' "#PBS -e scr_${NAME}_${EP}.stderr"
    printf '%s\n' 'cd $PBS_O_WORKDIR'
    printf '%s\n' 'source /apps/mambaforge/etc/profile.d/conda.sh'
    printf '%s\n' 'conda activate $HOME/envs/ash'
    printf '%s\n' 'export PYTHONNOUSERSITE=1'
    printf '%s\n' 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1'
    printf '%s\n' 'export PATH=/apps/openmpi/4.1.1/bin:$PATH'
    printf '%s\n' 'export LD_LIBRARY_PATH=/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:$LD_LIBRARY_PATH'
    printf '%s\n' "python3 gocat_screen.py ${NAME} ${EP} run > scr_${NAME}_${EP}.log 2>&1"
  } > "$RUN"/run.pbs
  ( cd "$RUN" && qsub run.pbs )
  echo "submitted ${NAME} ${EP}  (dir $RUN)"
done
