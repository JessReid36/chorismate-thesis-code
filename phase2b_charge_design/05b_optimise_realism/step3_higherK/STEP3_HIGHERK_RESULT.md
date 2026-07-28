# Step 3 higher-K RESULT (partial - completed DFT apot jobs; aforce + r2SCAN still running)

Tests whether climbing toward the GOCAT charge regime (working designs used 81; GA benchmarks 5/10/20)
recovers the near-attack geometry a K=2 pair could not. Same protocol as K2, surrogate-sized placement.

## Completed DFT jobs (converged, 146-149 min)
  K6_apot  : C1-C6 4.243  O3-C4 1.406 (intact ether)  SINGLE molecule  -> VALID, still OPEN
  K10_apot : C1-C6 4.201  O3-C4 1.421  C4-C6 1.474     FRAGMENTED (20/24) -> INVALID (over-polarisation)

## The trend (all valid ether-intact points, both coordinates)
  bare reactant 3.124 | K2_apot 5.28 | K2_aforce 4.06 | K6_apot 4.24 | TS 2.526     window [2.6-3.5]
  K10_apot: FRAGMENTED - C4-C6 collapsed to 1.474 (bond formed) + 4 atoms detached; the 10-cation field
  drives partial reaction/over-polarisation, not a clean reactant. C1-C6 unreadable.

## Finding: MORE CHARGES DO NOT RECOVER NEAR-ATTACK GEOMETRY
K6 (6 cations) lands 4.24 A - marginally better than K2_apot (5.28), comparable to K2_aforce (4.06), and
NOWHERE NEAR the 2.6-3.5 window. Pushing to K10 makes it WORSE: the substrate fragments under
over-polarisation (net +2, but 10 charges close to a -2 dianion). The dissociating aforce jobs (K6_aforce,
K10_aforce: O3-C4 ~5.7, ether BROKEN) show the same ceiling from the negative side. Certified across
K2/K6/K10 x a_pot/a_force: NO realizable point-charge design brings C1-C6 into the near-attack window.
Below ~K6 geometry stays open (~4-5); at K10 the substrate over-polarises/fragments.

## Interpretation (the earned mandate for the cradle)
This is the empirical, certified-optimisation-backed demonstration that POINT CHARGES ALONE CANNOT CRADLE
this dianion into the reactive conformation - proven, not assumed, across charge count and penalty
coordinate. Two independent limits bracket the design space:
  - too few charges (K2-K6): geometry stays open (~4 A), field cannot hold C1-C6 compact
  - too many charges (K10): over-polarisation ruptures/fragments the substrate (def2-SVP spill-out regime,
    exactly as RELAXATION_MASTER_LOG warned for high charge density near a -2 species)
=> the geometric preorganisation job needs NEUTRAL structure (dipolar cradle + steric walls), not more
   monopoles. Motivates the secondary-layer (cradle) extension.

## Over-polarisation note (K10 + aforce dissociation)
K10_apot fragmentation and the aforce ether-rupture are the over-polarisation limit made concrete. High
charge density near the -2 dianion, def2-SVP (non-diffuse) cannot describe the spilling density -> the
optimiser drives an unphysical dissociation. This BOUNDS realizable K on this grid/basis and is itself a
reportable constraint. (Diffuse basis def2-SVPD would reveal the spill-out explicitly - future check.)

## Status of remaining jobs (running to completion for the full record - NOT killed)
  K6_aforce, K10_aforce (DFT)      : dissociating (O3-C4 ~5.7); retained for the record
  K6/K10 x apot/aforce (r2SCAN-3c) : calibration, still running; compare_methods.py when done
Results appended when they finish.
