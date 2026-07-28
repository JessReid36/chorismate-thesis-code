# Charge-placement excluded-volume: point-spacing bug + surrogate-sized fix + higher-K (GOCAT regime)

## The bug (flagged by K3/K4 steric checks)
The MILP's no-overlap constraint used min_dist=2.5 A - the center-center spacing appropriate for POINT
charges (the original Tier-1 abstraction, 05_optimise/tier1_sweep.py). But Phase-2b places MOLECULAR
SURROGATES (guanidinium/formate). Surrogate body radius (center->farthest atom) = 1.33 A (guan), 1.25 A
(formate). Two surrogate CENTERS at 2.5 A => their ~1.3 A bodies interpenetrate. Empirically, several
K3/K4 mu=100 designs built as molecular surrogates CLASHED:
  K3_aforce (centers 2.74 A) -> surr-surr 0.33 A  (severe overlap, not buildable)
  K3_apot, K4_apot           -> surr-substrate 1.98 A (surrogate atom inside substrate vdW)
Only K4_aforce happened to be realizable. The grid itself has global min spacing 1.5 A (04_grid), so
min_dist is the ONLY guard - and 2.5 was sized for points, not bodies.

## The fix (new script, old one retained)
tier1_sweep_preorg_surrsized.py (NEW; tier1_sweep_preorg.py/_v2.py unchanged). Sets the no-overlap
distance to the SURROGATE BODIES, orientation-independent:
  safe center-center = r_body_i + r_body_j + clearance(2.0 A)
  worst case guan-guan = 1.33+1.33+2.0 = 4.66 -> DEFAULT min_dist = 4.7 A
This GUARANTEES >=2.0 A atom-atom for any pair in any orientation. Verified: new K4/K6 a_force mu=100
designs give surr-surr 4.78/3.50 A (clash-free). --min_dist is a CLI flag (4.7 surrogate / 2.5 point).
Certificate preserved: it is still a linear no-overlap constraint, gap 0.000 at every K.
[A per-pair-type min_dist (guan-guan 4.66, guan-form 4.58, form-form 4.50) would allow slightly tighter
formate packing; deferred - the single conservative value is simpler and the grid stays usable.]

## Higher K - the GOCAT regime (why K=2-4 was always going to be too sparse)
GOCAT primary sources (Dittner 2019, Behrens 2024): GA benchmarks at N_Ch = 5, 10, 20; WORKING
Diels-Alder catalyst designs used **81 point charges on a sphere** (Figs 7.11-7.13, incl. the
concerted->zwitterionic mechanism switch). Behrens' whole follow-on is about finding the *minimum*
required charge count - i.e. "how many charges" is treated as open and non-trivial, and their empirical
working answer was ~81. Our K=1-4 is ~20x sparser than any working GOCAT design and below even their
smallest GA benchmark (5). This CONTEXTUALISES the K2 result: a 2-charge pair failing to preorganise the
dianion (a_pot 5.28, a_force 4.06, both certified, both outside the 2.6-3.5 window) is EXPECTED, not a
failure of the method - it is the sparse-charge limit. The new script accepts arbitrary --K so we can
climb toward the GOCAT regime and test whether more charges recover near-attack geometry.

## Certified higher-K trend (proxy barrier; a_force mu=100, min_dist 4.7)
  K2 +6.52 | K4 +0.76 | K6 -5.45 | K8 -8.92 | K10 -11.97 | K12 -14.37   (all gap 0.000)
Proxy barrier drops steadily with K - MORE charges lower it further while keeping distortion drive ~0.
BUT proxy != geometry (the whole Step-3 lesson). Whether the RELAXED C1-C6 enters the window with higher K
is the open test -> relax a higher-K design and measure (next).

## Feasibility
At 4.7 A spacing the 331-site grid supports at least K=15 (certified). Max-K-that-fits is bounded by grid
extent; if we need GOCAT-scale (tens), the grid/shell density (04_grid) may need revisiting - flagged.

## Recommended new K values to test (given GOCAT)
Skip incremental K3/K5. Test K=6 and K=10 (a_force, surrogate-sized) as the "meaningfully more charges"
points bracketing the GOCAT GA-benchmark range (5/10/20), relax them, and see if C1-C6 reaches 2.6-3.5.
If even K=10 stays open -> pivot to the neutral cradle (charges alone insufficient at any realizable K).
