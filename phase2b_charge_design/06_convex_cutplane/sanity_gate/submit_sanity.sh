#!/bin/bash
# Sanity gate: relax the reactant under the 6-site FRACTIONAL net-neutral field via the PINNED oracle.
# Run from 06_convex_cutplane/sanity_gate/ AFTER build_sanity.py. Sets up run dir, writes PBS, qsubs.
set -e
[ -f sanity6_coords.xyz ] && [ -f sanity6_charges.txt ] || { echo "run build_sanity.py first"; exit 1; }
RUN=run_sanity
mkdir -p "$RUN"
cp ../oracle/gocat_screen.py "$RUN"/
cp ../inputs/reactant.xyz ../inputs/product.xyz ../inputs/ts.xyz "$RUN"/
cp sanity6_coords.xyz sanity6_charges.txt "$RUN"/

# PBS built with printf (no heredoc), env + resources mirrored from the committed submit_screen_batch.sh
{
printf '%s\n' '#!/bin/bash'
printf '%s\n' '#PBS -N sanity6_frac'
printf '%s\n' '#PBS -l select=1:ncpus=8:mem=32gb'
printf '%s\n' '#PBS -l walltime=24:00:00'
printf '%s\n' '#PBS -M 18660916@sun.ac.za -m ae'
printf '%s\n' '#PBS -o sanity6_frac.stdout'
printf '%s\n' '#PBS -e sanity6_frac.stderr'
printf '%s\n' 'cd $PBS_O_WORKDIR'
printf '%s\n' 'source /apps/mambaforge/etc/profile.d/conda.sh'
printf '%s\n' 'conda activate $HOME/envs/ash'
printf '%s\n' 'export PYTHONNOUSERSITE=1'
printf '%s\n' 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1'
printf '%s\n' 'export PATH=/apps/openmpi/4.1.1/bin:$PATH'
printf '%s\n' 'export LD_LIBRARY_PATH=/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:$LD_LIBRARY_PATH'
printf '%s\n' 'python3 gocat_screen.py sanity6 reactant run > sanity6_frac.log 2>&1'
} > "$RUN"/run.pbs

( cd "$RUN" && qsub run.pbs )
echo "submitted sanity gate: reactant relaxation under 6-site fractional net-neutral field"
