# The optimiser stage's code (development set)

Two scripts, written against the decisions of 20 October 2026 (`CHARGE_MODEL_DECISIONS.md`, "Update: the optimiser
stage decided", O1–O10; the reasoning is in `OPTIMISER_DECISIONS_20261020.md`):

| Script | What it does |
|---|---|
| `cs_opt_milp.py` | the MILP: the score, the rules, the pool of distinct designs, and the timing table |
| `cs_opt_null.py` | the null distribution: random designs under the same rules (O10) |

**Nothing either script produces on the development set is a result.** The sites, L and A are the development ones
(`cs_stage5/exploratory/dev_sites_a2.tsv`, `dev_L.tsv`, `dev_A/A_dev.tsv`), not Stage 5 production, and every run
says so in its own output.

## Running them

```
python3 cs_opt_milp.py selftest
python3 cs_opt_milp.py solve <sites.tsv> <L.tsv> <A.tsv> [key=value ...]
python3 cs_opt_milp.py time  <sites.tsv> <L.tsv> <A.tsv> [key=value ...]
python3 cs_opt_null.py selftest
python3 cs_opt_null.py run   <sites.tsv> <L.tsv> <A.tsv> [key=value ...]
```

Settings are `key=value` on the command line: `n`, `lam`, `budget`, `r_min`, `tangents`, `pool`, `distinct_r`,
`mip_gap`, `time_limit`, `seed`, and for the null script `draws` and `kind`. Every run prints the sha256 of each
input file and the settings it used.

## How each decision appears in the code

| Decision | Where |
|---|---|
| **O1** polarisation out of the score | the score is linear in q: `d_f = sum_i a_if q_i + sum_i L_if z_i`. Nothing quadratic enters the objective. |
| **O3** budget `sum q^2 <= B` | `s_i` per site with 21 tangent rows; `sum_i s_i <= B`. The true `sum q^2` is recomputed and printed for every design. |
| **O5** mean + λ × semi-deviation | `dbar` and `u_f >= d_f - dbar`, `u_f >= 0`; objective `dbar + (lam/F) sum_f u_f`. |
| **O6** n and the stopping rule | `n` is a setting; `time` mode solves n = 3 upward and prints the time, the solver's gap and the certified gap at each n. The stopping rule itself is applied to results, not inside the solver. |
| **O7** pool and distinct designs | `Problem.pool`: after each design, a cut forbidding designs whose spheres all lie within `distinct_r` of that design's. |
| **O8** "optimal" | every solve reports its gap, and the word "certified" appears nowhere. |
| **O10** null distribution | `cs_opt_null.py`, two kinds, with the sampling rule stated in its docstring. |
| Stage 5's two layers | `read_sites` reads `q_min_neg` and `q_min_pos` into one or two admitted ranges per site; at most one range per site may be chosen. |
| Stage 5's 3.25 Å spacing | close pairs grouped into cliques, one row each. |

**O2** (the reactant-gradient term) and the re-scoring step are deliberately **not** here: they act after the pool is
chosen, on full single points with gradients, and their threshold and weight come from Stage 6.

## Reaching n beyond about nine: the LP-guided search

Measured on the development grid at B = 2.0, the full search over all 2769 sites solves to a zero gap up to n = 9
(31 s at n = 3, 110 s at 8, 197 s at 9) and then returns **nothing at all** inside 300 s at n = 10. The cause is not
a weak bound and not infeasibility:

- designs at n = 10, 14 and 20 were constructed directly, respecting the spacing, so all of them are feasible;
- the relaxation is nearly exact — at n = 6 its value equals the integer optimum to four decimals — and it solves in
  4 to 7 s at every n.

The difficulty is finding the matching design among thousands of near-equivalent sites, which is the degeneracy the
records already describe (`HANDOFF.md` section 5: the frames agree on a favourable region of about 20 sites but
cannot order the sites within it).

`strategy=lp` therefore: relaxes the whole grid for a bound, takes the sites that relaxation uses (plus anything
within `grow`, 1.0 Å by default), and searches only those. What comes back is **a design of the full problem** — the
restricted problem has the full problem's rules on fewer sites, and the induced cliques still cover every close pair
inside the subset — and the bound is the **full** problem's, so the difference between them is a true optimality
gap: no design anywhere on the grid beats the one returned by more than that. Measured at B = 2.0:

| n | bound (whole grid) | design | certified gap | candidates | time |
|---|---|---|---|---|---|
| 10 | −16.1736 | −16.0455 | 0.128 (0.80 %) | 106 | 11 s |
| 14 | −16.8576 | −16.3454 | 0.512 (3.13 %) | 166 | 14 s |
| 20 | −17.5942 | −16.9124 | 0.682 (4.03 %) | 212 | 34 s |

**How far the candidate set should reach** (`grow`), measured at B = 2.0 on the development grid:

| `grow` | n = 10 | n = 14 | n = 20 |
|---|---|---|---|
| 0.0 Å (only the sites the relaxation used) | 0.80 %, 11 s | 3.13 %, 14 s | 4.03 %, 35 s |
| **1.0 Å (the default)** | 0.64 %, 47 s | 1.69 %, 111 s | 2.89 %, 421 s |
| 2.0 Å | 0.64 %, 79 s | 1.69 %, 418 s | no design in 420 s |

1.0 Å is the default because it tightens the gap where 0.0 Å is loosest and 2.0 Å buys nothing beyond it: at n = 10
and 14 the 2.0 Å search returns the same design and the same gap for two to four times the time, and at n = 20 it
returns none. Set `grow=0.0` when speed matters more than the last percent of the gap. Enlarging the candidate set
further - to a few hundred sites chosen by their mean effect - reinstates the original difficulty, which is why the
candidates come from the relaxation rather than from a score.

`strategy=auto`, the default, gives the full search `probe` seconds (120 by default). If it finishes, its design is
returned. If not, the LP-guided search runs too, and the better of the two designs is kept, with the tighter of the two
bounds (both are valid for the whole grid: the relaxation's, and whatever the full search had proved when it stopped).
Self-test 11 forces each branch. The switch is decided by what the solver does
on the machine running it, not by a chosen n, because how far the full search reaches depends on the HiGHS build: the
same n = 7, 8 and 9 solved exactly in 39, 102 and 194 s under scipy 1.18.1 and returned nothing in 600 s each under
scipy 1.11.4. On the older build the fallback at n = 7, 8 and 9 returned the same objectives as the exact search on
the newer one (−15.2136, −15.5435, −15.8251), with certified gaps of 0.012, 0.042 and 0.086 kcal/mol: the LP-guided
designs were in fact the best ones there, but the search could not prove it to zero, so the row says so.
`strategy=full` forces the full search at every n, and reports honestly when it finds nothing; `strategy=lp` forces
the LP-guided search. Decision O8 is unaffected: every solve still reports a gap as a property of this score and these rules, and
nothing is called a certified optimal catalyst.

## Two things worth knowing about the formulation

**The tangent constants are tied to the indicator.** Each tangent is written `s_i >= 2c q_i - c^2 z_i` rather than
`s_i >= 2c q_i - c^2`. For `z_i` in {0, 1} the two are the same rule, so the designs allowed are exactly those the
decisions describe and the 0.0025 e² per sphere bound is unchanged. The difference is in the solver's relaxations,
where `z_i` may be fractional: without the indicator, a hundredth of a sphere buys a hundredth of a charge for a
ten-thousandth of the budget, so the relaxed bound collapses. Measured on the development grid at n = 4: the plain
form relaxes to −19.32 kcal/mol and HiGHS returns no design at all inside 120 s, while the tied form relaxes to
−13.54 and the solve finishes with a zero gap. Self-test 9 checks that both forms give the same optimum.

**Tangents that cannot bind are not written.** For a site whose admitted charges lie in [lo, hi], the tangent at c
can only bind where `hi > c/2` (c > 0) or `lo < c/2` (c < 0); at c = 0 it never binds, because `s_i >= 0` already.
Omitting the rest leaves the feasible set untouched and removes about a third of the rows on this grid (40 930
instead of roughly 55 000).

## What the self-tests check

`cs_opt_milp.py selftest` — ten checks, each against something known independently of the solver:

1. **one sphere, closed form.** With n = 1, λ = 0 and a budget that does not bind, the best design is computable by
   hand; the MILP must match it, site and charge.
2. **every rule honoured** on a crowded instance using both layers: cardinality, charge ranges, spacing, and the
   arithmetic of `d_f`, the mean, the semi-deviation and the objective, all recomputed from the returned numbers.
3. **the q² stand-in stays inside its bound**: the true `sum q^2` never exceeds B by more than 0.0025 e² per sphere.
4. **charges proportional to −a**: with λ = 0 and no range bound active, the closed form is
   `q = -sqrt(B) abar / ||abar||`; checked against it.
5. **the optimality conditions with a bound active**: the budget's multiplier µ is recovered from an interior site
   and every other site is checked against `q = -a/(2µ)`, with bounded sites taking the sign that lowers the score.
6. **cliques equal the pairwise rules**: every close pair is covered; every pair inside a clique really is closer
   than r_min, so no valid design is excluded; and the two forms give the same optimum on a small instance.
7. **distinct designs**: supports differ, each design places at least one sphere beyond `distinct_r` of every sphere
   of each earlier design, and the pool is in order.
8. **L on its own**: (a) with A = 0 the design takes the n most negative mean-L sites; (b) with B = 0 the charges sit
   at the stand-in's own bound of 0.05 e rather than at zero, which is the slack of the tangent step expressed as a
   charge, and the design matches the closed form `mean L - 0.05 |mean a|` per site.
9. **the two tangent forms agree**: plain and indicator-tied give the same optimum, and the tied one relaxes tighter.
10. **the LP-guided search is sound**: its design obeys every rule of the full problem and its arithmetic re-checks;
    it is never better than the full optimum; its bound really does bound the full optimum; and the gap it reports
    covers the real shortfall against the full optimum.

`cs_opt_null.py selftest` — five checks: the charge draw is uniform inside the budget ball (checked against the
median radius a uniform ball must have), random sites always respect the spacing, out-of-range charges are moved to
the nearest admitted value, the two kinds agree with the MILP's own arithmetic on forced sites and kind (ii) beats
kind (i) on average, and percentiles read the right way round.

## The budget and the cation layer

Not a fault in the code, but a property of the grid and the budget that shapes both the designs and the null
distribution, and that the calibration of B (decision O3) should have in view.

**1228 of the 2769 development sites cannot carry zero charge.** They are the whole cation layer, which admits
positive charges only, from each site's own floor (Stage 5): 1064 sites admit only |q| = 1, 72 at least 0.75, 82 at
least 0.5 and 10 at least 0.25. A sphere on one of them costs the budget at least its floor squared, whatever charge
the design would prefer.

Two consequences:

1. **The budget caps the number of cation-layer spheres.** At most about B of them can sit on full-charge sites: two
   at B = 2.0, the value used here only as a code test. Every further sphere must sit on the anion/neutral layer,
   the only layer where small and zero charges are admitted. The n = 12 development design shows it: one sphere at
   +1.00, the other eleven between 0.15 and 0.45 in size. Stage 2 put the resting cations, and the strongest
   barrier-lowering sites of the map, on the cation layer, so the budget bites hardest exactly there. The
   stress-test designs that calibrate B ("the whole budget as close to the substrate as the grid allows") will meet
   this directly.
2. **Most random site sets cannot be charged inside a tight budget at all.** On their floors alone, before any
   charge is chosen:

   | | B = 2.0 | B = 5.0 | B = 10.0 |
   |---|---|---|---|
   | n = 4 | 22.6 % break it | 0 % | 0 % |
   | n = 10 | 88.8 % | 22.5 % | 0 % |
   | n = 20 | 99.9 % | 92.3 % | 17.6 % |

   So the null distribution applies decision O10's rule as written - "sites whose minimum charge alone exceeds the
   remaining budget are skipped" - while the sites are drawn, not afterwards. The first version of cs_opt_null.py
   drew the sites first and rejected afterwards; at n = 20 it would have rejected almost every draw, and whatever it
   did return would have been conditioned on lucky draws without saying so. Found in the second rehearsal, fixed
   before anything was committed, and covered by a self-test in each direction.

### What the cation layer does to the cost of the null distribution

Skipping budget-breaking sites (above) makes every drawn SET OF SITES chargeable, but decision O10's charge rule then
draws the charges uniformly inside the whole budget and moves each to the nearest value its site admits. On a
cation-layer site that move is up to the site's floor, so the floors and the drawn charges are both paid from the
same B, and the redraw rate grows with n. Measured on the development grid (kind (i), 100 designs each):

| | B = 2.0 | B = 5.0 |
|---|---|---|
| n = 4 | acceptance 0.36 | 1.00 |
| n = 10 | 0.020 | 0.28 |
| n = 20 | 0.0004 (490 s per 100 designs) | not measured |

So at a small B and large n, O10's 10 000 designs per n would take about 14 hours per n for kind (i) alone. Every
design it does draw obeys the rules, and the acceptance rate is reported, so nothing is wrong with the numbers; the
issue is cost. Decision O10 already says the exact sampling rule "is a pre-registered choice. It shapes the
distribution, so it is stated", so this belongs in the optimiser-stage pre-registration, which is written before B
is calibrated and before n is chosen. Two candidate rules for it:

- **(a) the rule as written, kept** (what the code does now). Affordable wherever the calibrated B is a few e² or
  more, or the stopping rule ends in the single digits.
- **(b) floors first.** Each sphere's floor is paid first; the remaining budget, B minus the floors, is what the
  charges above the floors are drawn uniformly inside. Every draw is then accepted by construction. It is a
  different distribution from (a), so it would be stated as the rule, not as a speed-up.

Claude's assessment of which is more likely to succeed: (a) is the decided rule, and whether it is affordable
depends on two things not yet known (the calibrated B and where the stopping rule ends). Fixing (a) with a stated
switch to (b) if (a)'s acceptance at the calibrated B and a carried-forward n falls below a pre-registered rate is
the choice that does not depend on any design's score, because acceptance depends only on the grid and the budget.
(b) is not built yet.

## Versions, and what reproduces

Both scripts need only numpy and scipy (HiGHS through `scipy.optimize.milp`). All 78 self-test checks pass on
python 3.12.3 + numpy 1.26.4 + scipy 1.11.4 (Jessica's PC) and on python 3.13.16 + numpy 2.5.3 + scipy 1.18.1.

**The designs reproduce on the numbers, on both searches.** The same development-set instance solved under both
versions returned the same answer in every digit printed:

| | n = 5, full search | n = 12, LP-guided |
|---|---|---|
| objective | −14.2489 both | −16.3531 both |
| bound and gap | exact, 0 | −16.5430, gap 0.1899 both |
| sites and charges | identical | identical (12 sites, 12 charges) |
| wall time | 16.6 s / 15.0 s | 78.4 s / 377.0 s |

So a rebuild of these outputs can be checked against the committed numbers directly, which is the check this project
uses everywhere else. The one thing that does **not** carry across is the time: the same LP-guided solve took 78 s
under scipy 1.18.1 and 377 s under 1.11.4, nearly five times longer, so every timing in this project's records is
labelled with the version that produced it and timings must not be compared across versions.

Two cautions, because neither was proved and both would be easy to assume:

- this was checked on two instances, one per search, not on every n. Where many designs tie exactly, nothing in
  either HiGHS build guarantees the same one is returned, so a future difference in the *design* at an equal
  objective would not be a fault.
- the LP-guided search depends on which sites the relaxation uses. That came out the same here (159 sites both
  times), but a degenerate relaxation could in principle be solved at a different vertex by a different build, which
  would change the candidate set. The bound and the design it certifies would still be valid; only the particular
  design could differ.
