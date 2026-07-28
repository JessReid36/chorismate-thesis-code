# 05d_shell_ladder — WORKFLOW

## Goal
Find the CLOSEST vdW shell at which K=10 discrete +-1 charges catalyse WITHOUT rupturing the substrate.
Read PLAN.md for the interpretive caveat (equidistant +-1 = uniform strength, unlike GOCAT's graded [-1,1];
a null result motivates but does not preempt the fractional optimiser).

## Two parts
### PART 1 - existing shells (2/3/4 A): RUN NOW (grid already exists)
Designs already optimised (K=10, certified): shell2p0 (proxy -11.22), shell3p0 (-4.74), shell4p0 (-2.70).
coords+charges are in relax/. Just relax them:
  cd relax
  for s in 2p0 3p0 4p0; do mkdir -p run_shell$s; cp relax_shell.py shell${s}_coords.xyz shell${s}_charges.txt run_shell${s}.pbs run_shell$s/; done
  # dry-run on a compute node, then:
  cd run_shell2p0 && qsub run_shell2p0.pbs && cd ..   (etc for 3p0,4p0)
Analyse: python3 <path>/analyse_relaxed.py run_shell$s/shell${s}_relaxed.xyz

### PART 2 - outer shells (5/6/7/8 A): need grid extension + QM vpot first
1. grid_ext:  python3 rebuild_sdf_ext.py <R.xyz> <TS.xyz> <P.xyz>   (margin 9)
              python3 extract_sample_ext.py 5.0 6.0 7.0 8.0          -> outer_shell_points.tsv
2. dvpot_ext: edit run_vpot_ext.sh paths to your 02_singlepoints .gbw/.densities, then run on HPC
              -> dv_grid_ext.tsv ; then merge: cat existing dv_grid.tsv + dv_grid_ext.tsv -> dv_grid_full.tsv
   (REQUIRES the committed R and TS .gbw/.densities from 02_singlepoints - if discarded, re-run those SPs)
3. optimise:  for s in 5.0 6.0 7.0 8.0; do python3 optimise_shell.py dv_grid_full.tsv $s <R.xyz> --K 10; done
4. relax:     same as Part 1, for shell5p0..8p0 (PBS wrappers already built)

## Run ALL in parallel
Each shell = its own PBS job (distinct name ladder_<shell>, own run dir). Submit all 7 (2-8A) once their
coords exist. K=10 per shell, own scratch dir (collision-safe).

## Decision
Closest INTACT shell = usable discrete sweet spot. None intact even at 8A -> uniform +-1 cannot
catalyse-without-rupture at any distance -> fractional/graded field (certified QP) required.
