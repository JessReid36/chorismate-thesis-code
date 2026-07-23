# Tier-2 Endpoint-Relaxation — Master Attempt Log (what we tried, what we learned, what NOT to redo)

Goal: re-optimise the bare 24-atom chorismate substrate under a fixed external ±1 point-charge
design (the certified Tier-1 field), as the first step of Tier-2 relaxed-path validation. This
log records every approach tried so we never repeat a dead end. Evidence lives in the results
repo (07_stage0/ for native; relaxation_attempts/ash_guarded/ for ASH).

Reacting-atom indices (0-based): C1=0, ether O3=7, C4=8, C6=12. Level (fixed throughout):
B3LYP-D3BJ/def2-SVP/def2-J/RIJCOSX/CPCM(eps=4), charge -2, mult 1. Designs: max-lowering K1-K4;
the +1 sits 3.84 A from O3 (the Arg90 valley) in K2/K3/K4 and is the implosion driver. K1 = lone
-1 at 8.40 A, no +1 -> the only benign design.

## Attempt 1 — Native ORCA bare %pointcharges + %geom Opt (coordsys cartesian).  CLOSED.
8 jobs (K1-4 x reactant/product). Result: 0/8 converged (72-cycle cap); 6/8 imploded (K2/K3/K4
both endpoints) — a nucleus pulled to 0.01-1.06 A of a charge from a >=3.5 A start; O3 torn from
C4 at products (O3-C4 5.3-5.9 A). K4_product diverged (MAX grad ~4767). K1 intact but did NOT
converge (MAX grad plateau ~0.0022 vs 3e-4 tol). SCF converged throughout (failure is geometric,
not electronic). Evidence: results 07_stage0/runs/*, STAGE0_results_table.md, STAGE0_FINDINGS.md.
DO NOT RETRY: bare point charges provide no wall against nuclear approach; the +1 pulls electron-
rich atoms in directly.

## Attempt 2 — Native bare + Calc_Hess true + Recalc_Hess 5 (break the plateau).  CLOSED.
Same 8. Result: 0/8 converged at ~4.5 h/job (Recalc_Hess ~10x slower for zero gain). A full
analytic Hessian recomputed every 5 cycles did NOT fix the K1 plateau: MAX grad oscillates
6-30x over tol with no downward trend. Evidence: results 07_stage0/stage1/. 
KEY INSIGHT: the non-convergence is a FRUSTRATED, shallow, field-distorted surface, not a
curvature deficit -> more Hessian is useless. Loose reaction-relevant tolerance also ruled out
(K1 bounces above even 1e-3). DO NOT RETRY Calc_Hess/looser-tol on native bare charges.
COORDINATE-SYSTEM NOTE: in a FIXED external field the net rigid-body force/torque is physical;
redundant/TRIC internals are translation/rotation-invariant and cannot relax it -> contraindicated.
Use Cartesian or HDLC.

=> DECISION after Attempts 1-2: native %geom Opt + bare point charges is exhausted for these
designs. Both failure modes (implosion; frustrated-surface non-convergence) are outside what
%geom can fix. Shift to ASH (QM/MM: ORCA QM + OpenMM MM charge sites + geomeTRIC opt) so we get
(a) a distance guard and (b) a restraint-capable optimiser. See DESIGN_DECISIONS.md.

## ASH environment (works; reproducible).  ash_guarded/ash_environment.yml
Env at $HOME/envs/ash: python 3.11, OpenMM 8.5.2, geomeTRIC 1.1.1, ASH 0.95 (from
git+https://github.com/RagnarB83/ash.git — NOT a PyPI package named ash-multiscale).
BUILD GOTCHAS (don't relearn): system python is 3.6 (too old); /apps/mambaforge/envs is
read-only (use --prefix $HOME/envs/ash); the conda solve OOMs on the login node -> build in a
PBS job. RUNTIME (every ASH job): PYTHONNOUSERSITE=1 (a stray ~/.local numpy shadows the env) +
OPENBLAS/OMP/MKL/NUMEXPR_NUM_THREADS=1 (login-node OpenBLAS OOM) + OpenMPI on PATH
(/apps/openmpi/4.1.1/bin) or ORCATheory init aborts. ASH prints a banner to STDOUT -> never parse
`python -c "...print(path)"` for a path; use a shell glob.

## ASH architecture (verified from source + docs).
- ORCATheory has NO point-charge argument -> external charges MUST be an OpenMM MM region fed via
  QMMMTheory(qm_theory, mm_theory, qmatoms, charges=[...], embedding="elstat"). QM charges auto-
  zeroed; MM +/-1 live in the field; no covalent boundary/link atoms (charge sites are isolated).
- Optimizer(...) constraints are FIXED-VALUE only: constraints={'bond':[[i,j,value]]},
  constrainvalue=True. No one-sided/flat-bottom option in the geom module.
- One-sided wall => inject an OpenMM CustomBondForce into mm.system BEFORE the QMMMTheory wrap.
- geomeTRIC: use HDLC (coordsystem='hdlc') for these systems; TRIC warned against.
Charge sites are added to the Fragment as He placeholders (element irrelevant for MM points),
indices 24 (-1) and 25 (+1).

## Attempt 3 — ASH QM/MM, fixed pin O3(7)->+1(25) at 3.84 A (job 395168).  INFORMATIVE FAIL.
Result: the pin holds O3 cleanly (3.82 A; O3-C4 1.50; C1-C6 3.14) BUT the implosion RELOCATES —
carboxylate oxygen O5 collapses onto the +1 at 0.955 A. Evidence: results
relaxation_attempts/ash_guarded/pin_run/.
KEY INSIGHT: the +1 is a NON-SPECIFIC electron-density attractor; pinning one pair just frees the
next electron-rich atom. A per-pair pin is the wrong shape of guard. DO NOT RETRY single-pair pins.
=> the guard must be GLOBAL (every substrate atom vs every charge).

## Attempt 4 — ASH QM/MM, GLOBAL one-sided floor, k=100 kcal/mol/A^2, floor 2.2 A (wall run).
CustomBondForce("step(r0-r)*0.5*k*(r-r0)^2") over all 24x2 atom-charge pairs. Result: 81 steps /
73 min, NO convergence at maxiter=80. The wall ENGAGED (stopped the 0.02 A total collapse) but was
PENETRATED: an atom (C3) settles at ~1.26 A from the +1 in a STABLE force balance (min atom-charge
1.16-1.27 A across steps 30-80), O3 pinned at the floor (2.14 A). Evidence: results
relaxation_attempts/ash_guarded/wall_run_k100/.
KEY INSIGHT: k=100/floor=2.2 is too SOFT and too LATE — at 2.2 A a -2-density atom is already deep
in the +1 attractive well, and a soft quadratic can't hold it. The wall works in principle; the
parameters don't.

## Standing conclusions / DO-NOT-REDO
1. Native bare %pointcharges + %geom Opt: CLOSED (implosion + frustrated-surface non-convergence).
2. Calc_Hess / Recalc_Hess / looser tolerance on native bare: useless, expensive.
3. Redundant/TRIC internals in a fixed field: contraindicated (can't relax rigid-body field force).
4. Single fixed-pair pin guard: relocates the implosion. Guard must be global.
5. Soft global floor k=100/2.2 A: penetrated at ~1.26 A. Too soft/late.
The ASH QM/MM pipeline ITSELF is validated (constructs, embeds, runs, wall injects) — the open
问题 is guard parameters and whether a strong field is stabilizable at all.

## Next (planned, NOT yet done)
- Wall iteration 2: stiffen + widen -> k~1500 kcal/mol/A^2, floor ~2.7 A (intercept before the
  deep well), on K2. This is the "can ANY reasonable guard hold K2?" test.
- Run K1 (benign, never implodes) through the SAME ASH guarded pipeline to confirm it converges on
  an easy design before declaring strong designs intractable (isolates guard-tuning from
  field-too-strong).
- Emerging hypothesis to test, not assume: K2's certified-optimal-on-proxy field may be physically
  too aggressive to admit a relaxed substrate — which would itself be a Tier-2 result
  (proxy-optimal =/= physically realisable), consistent with the frozen-field proxy's known limits.
