# Decision record: after the QMAX=0.3 barrier -- optimizer vs proxy-ceiling

## Context
Certified fields fragment the bare dianion at every K (carboxylate 12,21). The cradle (2x+1
counter-charges) holds the substrate (both endpoints VALID). Catalytic +-1 fragments the ring;
magnitude sweep found QMAX<=0.3 holds both endpoints. Barrier on cradleK4q0p30:
  - relaxed scan gave +29.43 kcal/mol but DERAILED past the peak (67 kcal/mol cliff-drop at
    O3-C4=2.40 -- lost the reaction path), so NOT a clean number.
  - NEB-TS running as the trustworthy cross-check (the decision hinges on it).

## Two SEPARATE problems (do not conflate)
1. Bound-slamming: all catalytic charges land at +-QMAX. This is a THEOREM (linear objective optimizes
   at a polytope vertex), not a bug. Fix = convex quadratic objective -> varied interior charges,
   certificate preserved.
2. Proxy ceiling: Sum q.dV on FROZEN geometries is the first-order "extreme approximation" limit of
   DTSS (Beker & Sokalski JCTC 2020). It locates WHERE charges go (a field map), not the relaxed-path
   barrier. GOCAT's own fixed-path version "allowed for only small catalytic effects"; they needed
   on-the-fly path optimization for real catalysis. The +29-vs-17.5 result is the textbook symptom.
VARIED MAGNITUDES FIX #1, NOT #2. They are orthogonal. Solve both; don't expect the first to solve
the second.

## Staged plan (in order; each gates the next)
1. READ THE NEB. If barrier < 17.47 -> the scan derailed and we HAVE catalysis; premise flips, stop
   here and celebrate/verify TS. If >= 17.47 -> uncatalytic, proceed.
2. BARRIER-CORRELATION TEST (the decision-maker, BEFORE any new optimizer): take designs we ALREADY
   have (K-ladder K6/10/15/20, cradle, cradleK4, sweep q0p2/q0p3) whose Sum q.dV objectives are known,
   compute their TRUE relaxed barriers (scan/NEB), plot true-barrier vs objective. Poor correlation
   (likely, given +29-vs-17.5) => the proxy has a ceiling for this substrate; reweighting it (varied
   magnitudes included) will NOT catalyse.
3. IF proxy correlates (surprising): build Tier-1 convex QP -- minimize dV^T q + lambda q^T W q
   (+ mu||q||_1 for sparsity), W = reacting-bond field Gram matrix (penalize the field that tears
   bonds, not raw magnitude), Sum q = 0, |q|<=QMAX. Solve on CLARABEL (real primal-dual certificate;
   HiGHS-QP as cross-check). Gives varied, buildable, integrity-aware charges with certificate intact.
   For exactly-K: Tier-2 MIQP on SCIP (open) / Gurobi (free academic, auditable gap); HiGHS CANNOT do
   MIQP.
4. IF proxy does NOT correlate (likely): STOP reweighting. Fork:
   (a) Response-based objective: minimize -(dMu)^T C q s.t. field trust-region ||Bq||<=Fmax (convex
       SOCP/QP, certifiable on Clarabel; Shaik OEEF logic; barrier-correlated via the real TS-reactant
       dipole difference). Requires computing dMu, dAlpha for chorismate. Convexity caveat: full
       quadratic -1/2 q^T(C^T dAlpha C)q is only convex if -dAlpha is PSD (TSs often MORE polarizable
       -> indefinite -> lose certificate; then restrict to the linear-dipole + trust-region form).
   (b) Neutral STERIC cradle: geometric preorganization (near-attack conformation), a different
       (non-electrostatic) catalytic channel that sidesteps the tearing-vs-catalysis vise. Burschowsky
       PNAS 2014 caveat: Arg90->citrulline "poor catalyst even though it preorganizes" -- so geometry
       ALONE may also be limited; the enzyme needs BOTH a placed charge and geometry.

## Novelty reframing (the north star, regardless of outcome)
Strongest form of "beats GOCAT" is NOT "certified optimum of Sum q.dV" (a proxy that may not track the
barrier) but "certified optimum of a physically-grounded CONVEX SURROGATE of the relaxed-path barrier"
-- same physical target as GOCAT, provable convex inner problem instead of a genetic algorithm. The
METHOD is the contribution; it survives even if this substrate proves uncatalytic by point charges.

## Solver facts (for the write-up)
- HiGHS: LP + convex QP + integer-LINEAR MIP only. CANNOT do MIQP. Use native highspy (scipy-HiGHS
  segfault is scipy issue #15456, not a HiGHS defect).
- Clarabel: interior-point, genuine primal-dual + infeasibility certificates, native QP. Certificate-
  of-record for continuous QP.
- SCIP: open MIQP/MINLP, certified gap, inspectable. Gurobi: free academic, fastest MIQP, auditable
  gap (closed core -- cross-check against SCIP).
- OSQP/SCS: first-order ADMM, tolerance-based residuals -- pre-screen ONLY, not certificate-of-record.

## Key citations to verify in primary PDFs before formal use
- Beker & Sokalski, JCTC 2020, 16, 3420 (DOI 10.1021/acs.jctc.0c00139) -- point-charge Sum q.dV as the
  extreme-approximation DTSS limit.
- Dittner & Hartke, JCP 2020, 152, 114106 (10.1063/5.0023434) -- GOCAT fixed-path "small effects only".
- Besalu-Sala et al., ACS Catal. 2021, 11, 14467 (10.1021/acscatal.1c04247) -- FDBbeta response model.
- Bofill et al., JCP 2023, 159, 114112 (10.1063/5.0167749) -- PMED, sign-equivalent.
- Burschowsky et al., PNAS 2014, 111, 17516 (10.1073/pnas.1408512111) -- single placed charge (Arg90)
  dominates; preorganization "contributes relatively little".
