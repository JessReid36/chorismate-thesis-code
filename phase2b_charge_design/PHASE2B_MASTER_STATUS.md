# Phase 2b — Master Status & Open-Threads List
(chorismate dianion external-charge catalyst design; de novo reverse design, enzyme = yardstick)

## A. SETTLED RESULTS (committed or ready to commit)
1. **Charges alone (discrete ±1) CANNOT recover near-attack geometry** — certified across K2/K6/K10 ×
   {a_pot, a_force}. K2af 4.06, K6ap 4.24 (valid, OPEN); K10 fragments; a_force ruptures at K6.
   Two limits bracket the space: too few → open ~4Å; too many → over-polarise/rupture.
2. **Neutral steric cradle alone FAILS** — heuristic 6-bead cage: K2af 5.31 (opened MORE than bare),
   K6ap 5.93 FRAGMENTED. Frozen LJ walls leak; substrate opens along routes the cage doesn't intercept.
3. **Certified set-cover cage machinery built** (Option A) — union grid (charge lattice + inner shell,
   435 pts). All 3 K get full 6-target certified cage (gap 0). Earlier "K10 uncoverable" RETRACTED
   (grid artifact). NOT relaxed (heuristic stronger-wall cage already failed).
4. **r2SCAN-3c FAILS calibration** — K6_apot: DFT valid (4.24, intact) vs r2SCAN FRAGMENTED (5.86).
   SIE-driven; no cheap fix (dispersion already in 3c; only exact-exchange hybrids fix it, at hybrid cost).
   → DFT-only. [Full compare_methods table pending remaining r2SCAN jobs.]

## B. KEY MECHANISTIC INSIGHTS (the physics we've established)
- **Linear objective → charges RAIL to ±1** (LP vertex solutions) → harsh field → fragmentation. The
  certificate (linearity) CAUSES the physical failure. Fundamental tension.
- **Continuous q + linear objective still rails** (LP theory: box + linear = vertex). Fractional needs
  a NONLINEAR (e.g. quadratic) objective.
- **Quadratic penalty (convex QP via HiGHS passHessian) → gentle fractional charges** (max|q| 0.06–0.29,
  distributed, certified convex-global). Certificate AND gentleness — but plain L2 is modest novelty.
- **Distance tradeoff (inner vs outer shell):** inner 2Å shell = strongest catalysis (|Dv| 0.0092) but
  harshest; outer 4Å = gentle but weak (|Dv| 0.0043). Forcing K6 to 4Å: barrier +3.47 (vs −6.93 all-shell)
  = catalysis gutted. Distance reduces harshness AND catalysis (coupled). Can't fix by distance alone.

## C. OPEN THREADS / THINGS TO TRY (the actual to-do list)
### C1. Cheap physical controls (do first — cheap, rule out simple explanations)
- [ ] **Outer-shell ±1 relaxation test** (user's "too close?" hypothesis): force K6/K10 to 3–4Å shells,
      relax, see if it fragments LESS than inner-shell. Confirms/refutes distance as a rupture contributor.
      Cheap (1–2 jobs). Bridges to the proximity-weighted penalty.

### C2. Novel certified optimiser (the main thrust — from the research report, ranked)
- [ ] **Stage 1 (low risk): convex QP with PHYSICAL regulariser** — proximity-weighted Σ w_i q_i²
      (w heavier near breaking C–O bond) OR RESP hyperbolic restraint. Certified convex-global, gentle,
      fractional. Solve HiGHS/OSQP/Clarabel. → relax → does C1–C6 hit window without rupture?
- [ ] **Stage 2 (novel core): MISOCP best-subset** — binary site-select + continuous fractional q +
      SOC field-energy budget ||Wq||₂≤τ. Certified gap. Sweep k,τ → Pareto front of TS-stabilisation vs
      over-polarisation (operationalises our empirical bracket). Gurobi(academic)/SCIP/Pajarito.
- [ ] **Stage 3 (aspirational): Lasserre moment/SOS or spatial B&B** on a fitted non-convex barrier
      surrogate over a REDUCED site set (10–30 from Stage 2) — certified-global of the true-ish objective.
      High risk (SDP scaling); proof-of-concept figure only.
- [ ] **Outer QM-validation loop** for all stages: solve certified proxy → validate few in ORCA →
      reweight penalty → repeat. Certificate on PROXY, physics validated externally (state explicitly).

### C3. Fractional-charge relaxation modelling (needed once QP gives fractional designs)
- [ ] Decide relaxation charge-model for fractional q: (a) bare SCALED point charges (simplest; gentle
      charges less collapse-prone), or (b) scaled partial-charge surrogates. Lean (a) first.

### C4. If external charges fundamentally can't (fallback)
- [ ] DIPOLAR/active cradle (oriented net-neutral dipoles, not passive LJ walls) — option a from report.
- [ ] Or accept + publish the CERTIFIED NEGATIVE result: external charge/steric design hits a wall for
      this dianion (comprehensive, certified — itself a thesis contribution).

### C5. Protocol / completeness (carry-forward)
- [ ] For EVERY K: test BOTH a_pot AND a_force (winner is K-dependent: K2→a_force, K6→a_pot).
- [ ] Add K20 to the charge-count ladder (toward GOCAT's 81).
- [ ] Full compare_methods.py table once all r2SCAN jobs done → commit calibration verdict.

## D. FIXED CONTEXT (don't re-derive)
- Substrate: chorismate dianion, 24 atoms, −2, singlet. C1=0,C6=12,O3=7,C4=8. Reactant C1-C6 3.124,
  O3-C4 1.472. TS C1-C6 2.526. Window [2.6,3.5]. Bare barrier +17.47.
- Level: B3LYP-D3BJ/def2-SVP/def2-J/RIJCOSX/CPCM(ε=4), ORCA 6.0.1. r2SCAN-3c FAILED (SIE).
- Grid: dv_grid 331 sites, shells 2/3/4Å (actual 3.2–6.1Å to substrate atoms), bounded ±1, NET-FREE, our grid.
- Solver: HiGHS (LP/MILP/QP). GOCAT = GA + continuous q∈[-1,1] + net-neutral + 81 charges on sphere (the foil).
- min_dist 4.7Å (surrogate-body-sized). MM surrogates: guanidinium+1/formate−1 (unit charges only).

## E. THESIS FRAMING (the through-line)
Novelty = DETERMINISTIC, PRINCIPLED optimiser with an OPTIMALITY GUARANTEE for external-charge catalyst
design, vs GOCAT's stochastic GA (no certificate). "Certified MILP" was ONE instance; the deeper claim is
certified/deterministic design. The TS-stabilisation-vs-over-polarisation tradeoff (our empirical finding)
becomes the OBJECTIVE (regularised convex/conic program) — turning the finding into the method.
Literature check (research report): NO prior certified math-programming optimiser for external-charge/field
catalyst design exists → the whole direction is novel.
