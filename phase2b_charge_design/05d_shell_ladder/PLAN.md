# 05d_shell_ladder — discrete +-1 charges, incremental vdW shells: at what distance do they catalyse without rupturing?

## Question (user)
Our optimiser clusters +-1 charges on the CLOSEST vdW shell (strongest, harshest) -> fragments. Move
OUTWARD incrementally, re-optimising K=10 discrete +-1 per shell, on the vdW surface (our approach, NOT
GOCAT's enclosing sphere). Find the CLOSEST shell at which discrete charges catalyse WITHOUT rupturing the
substrate = the usable sweet spot for the discrete model.

## CRITICAL INTERPRETIVE CAVEAT (user-flagged)
On an EQUIDISTANT shell with DISCRETE +-1, every charge has the SAME magnitude (1.0) AND the same distance
-> the field is UNIFORM-STRENGTH, tunable only by sign + angular position. GOCAT's advantage was CONTINUOUS
q in [-1,+1]: it could make some charges strong, some weak -> a GRADED field. So this ladder tests a MORE
CONSTRAINED design than GOCAT. If it FAILS at every shell, that does NOT rule out that FRACTIONAL charges at
the same distances would work - the failure could be "uniform +-1 is too blunt", not "distance is wrong".
=> the shell ladder is a clean test of DISTANCE FOR THE DISCRETE MODEL; a null result motivates (does not
preempt) the fractional-charge optimiser. Read the result with this confound in mind.

## Proxy already known (existing 2/3/4 A shells, K=6): catalysis FADES outward
  2A -6.19 | 3A -0.04 | 4A +3.47 (bare +17.47). Catalysis is STRONGEST close; what changes outward is
  RUPTURE, not catalysis. So the target is the CLOSEST shell that stays INTACT (max usable catalysis).

## Design: K=10 discrete +-1 (net-free), re-optimised per shell, shells 2,3,4,5,6,7,8 A (vdW-surface offset)
Existing grid stops at 4A -> BUILD outer shells (5-8A): rebuild SDF (margin 5->9), extract shells, sample,
then orca_vpot on the committed R and TS wavefunctions (02_singlepoints .gbw/.densities) to get V_R/V_TS/Dv
at the new points (E3 convention in TIER2_KNOWN_ISSUES). Then optimise K=10 per shell + relax each.

## Steps
  grid_ext/  : rebuild_sdf_ext.py, extract_shells_ext.py, sample_ext.py -> outer-shell candidate points
  dvpot_ext/ : run_vpot_ext.sh (orca_vpot on R,TS wavefns at new points) -> dv_grid_ext.tsv (all shells)
  optimise/  : optimise_shell.py (K=10 +-1 net-free per shell, HiGHS certified)
  relax/     : relax_shell.py + per-shell coords/charges + PBS (own dir per job, distinct names)

## Decision
Closest INTACT shell = usable discrete sweet spot. If NONE intact even at 8A -> discrete uniform +-1 cannot
catalyse-without-rupture at any distance -> the fractional/graded field (certified QP) is REQUIRED (and the
caveat above means this is the honest conclusion, not "distance failed").

## Dependency
Requires committed R and TS .gbw/.densities from 02_singlepoints on the HPC. If discarded, re-run those two
single points first (geometries known; cheap).
