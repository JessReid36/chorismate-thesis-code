# Step 3 higher-K: does more charge (GOCAT regime) recover near-attack geometry?

Tests K=6 and K=10 x {a_pot, a_force}, mu=100, SURROGATE-SIZED placement (min_dist 4.7, clash-free).
Same protocol as K2. GOCAT context: working designs used 81 charges (GA benchmarks 5/10/20); K6/K10
bracket their GA range. Question: does the RELAXED C1-C6 enter the near-attack window [2.6-3.5] with
more charges, where K2 could not (a_pot 5.28, a_force 4.06)?

## The four designs (all certified gap 0.000, steric-verified as molecular surrogates)
  K6_apot    : 72 atoms, MM net +2, surr-surr 3.50 surr-sub 2.49
  K6_aforce  : 66 atoms, MM net  0 (3+/3-, most enzyme-like balanced cavity), surr-surr 3.50 surr-sub 2.58
  K10_apot   :100 atoms, MM net +2, surr-surr 3.04 surr-sub 2.44
  K10_aforce : 88 atoms, MM net -2, surr-surr 2.97 surr-sub 2.53
Net charges modest (no dangerous same-sign pile-up); over-polarisation risk manageable. QM stays -2.

## Certified proxy-barrier trend (a_force mu=100): K2 +6.5 -> K6 -5.5 -> K10 -12.0. Proxy != geometry.

## Run
Dry-run each: python3 relax_higherK.py <name>       (expect DRY RUN OK ... LJ pairs, frozen)
Submit:       qsub run_<name>.pbs                    (one line each, four independent jobs)
Analyse:      python3 ../analyse_relaxed.py <name>_relaxed.xyz

## Decision
- Any design brings C1-C6 into [2.6-3.5] -> MORE CHARGES recover geometry; identify the charge count
  needed (between 4 and this K). Compare a_pot vs a_force (expect a_force >= a_pot, per K2).
- Even K10 stays open (~4-5) -> charges alone insufficient at realizable K on this grid -> pivot to the
  neutral cradle (secondary layer) with a fully-earned mandate.
Watch: connectivity (single molecule) + net-charge over-polarisation on the higher-net designs.

## Files
relax_higherK.py (one parameterized runner), <name>_coords.xyz + <name>_charges.txt (x4), run_<name>.pbs (x4)
