# Placement-strategy gap: point-spacing vs molecular-surrogate excluded volume (fix)

## The bug (why K3/K4 designs clashed)
Candidate sites (04_grid) were Poisson-thinned to r_min = 1.5 A between POINTS. The Tier-1 MILP
(tier1_sweep.py / _preorg / _preorg_v2) then forbade two SELECTED charges closer than min_dist = 2.5 A -
again a POINT criterion. But when a design is realised as MOLECULAR SURROGATES, each charge becomes a body:
guanidinium radius ~1.33 A, formate ~1.25 A (centre atom to farthest atom). Two surrogate CENTRES therefore
need to be >= ~3.7 A apart (radii sum + ~1.0 A vdW gap) or their atoms overlap. At 2.5-2.74 A centres the
bodies clashed (observed atom-atom down to 0.33 A on K3_aforce). So the point-charge placement is correct
FOR POINTS but has no notion of surrogate excluded volume - and this only gets worse at higher K.

## Required clearances (from templates, +1.0 A gap)
guan-guan 3.66 | guan-form 3.58 | form-form 3.50  -> conservative single threshold 3.7 A covers all.

## The fix (new script, old ones untouched)
tier1_sweep_preorg_surrogate.py : identical certified sweep (net-free, bounded +-1, CPCM-screened Dv,
linear preorg penalty mu*|D|, a_pot/a_force) but the excluded-volume constraint uses --surr_mindist
(default 3.7 A) between surrogate CENTRES instead of 2.5 A between points. Encoded as the SAME linear
pairwise no-overlap rows -> HiGHS gap = 0.000 preserved (certificate intact; ~3010 rows added, still solves).
This is exactly the report's Step-1 excluded-volume realism ingredient, now made surrogate-correct.

## Verified
K4 a_force mu=100 (new): sites 150/188/209/210, min centre distance 4.11 A -> CLASH-FREE. Certified.
Scales cleanly to higher K (all gap 0.000): certified proxy barrier at mu=100, near-zero distorting drive:
  K3 +2.68 | K4 -0.85 | K6 -8.00 | K8 -13.14   (proxy only - relaxed geometry still to be tested)
NOTE: new designs are DIFFERENT sites from the old (clashing) ones - they are the correctly-constrained
certified optimum, which is the right thing to relax.

## GOCAT context (why higher K is appropriate)
GOCAT (Dittner 2019) GA-benchmarked at N_Ch = 5,10,20 and its WORKING Diels-Alder catalyst designs used
N_Ch = 81 point charges on a sphere. Behrens 2024 is explicitly about finding the MINIMUM required number.
Our K=1-4 is ~20x sparser than any working GOCAT design; K2 provably cannot preorganise this dianion
(a_pot 5.28, a_force 4.06, D->0 achievable yet still 4.06). Moving to higher K with SURROGATE-CORRECT
spacing is the appropriate next test of whether more (realizable) charges recover near-attack geometry.

## Use for new K
Generate designs with tier1_sweep_preorg_surrogate.py (a_force, the resolved better coordinate) at the
higher K values, steric-verify (should now pass by construction), relax, measure C1-C6 vs the K2 4.06 anchor.
