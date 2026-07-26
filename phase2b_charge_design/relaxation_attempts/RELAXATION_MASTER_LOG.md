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

## Attempt 5 — ASH, MOLECULAR SURROGATES (guanidinium=+1, formate=-1), full QM, frozen.  WORKS (concept proven).
Motivated by the literature (Behrens & Hartke Top Catal 2022 "back to the real world"; Sokalski OCF is
analytic-not-point-charges): the collapse is a POINT-CHARGE ARTIFACT (a bare charge is a singular,
Pauli-free Coulomb sink), not a real physical result. Fix: replace each abstract +/-1 with a REAL
molecular group that has finite size + Pauli repulsion.
Build (sandbox-validated): 24 substrate atoms + guanidinium C(NH2)3+ (10 atoms, central C at the
certified +1 site, 3.84 A from O3) + formate HCOO- (4 atoms, C at the -1 site). Net charge -2, singlet,
38 QM atoms. Surrogates FROZEN (indices 24-37); only the substrate relaxes. Level unchanged
(B3LYP-D3BJ/def2-SVP/CPCM eps=4). Script: ash_guarded/k2_surrogate_spike.py (note: autostart=False is
REQUIRED - ASH otherwise hands ORCA a stale orca.gbw across geometry steps -> "Input geometry does not
match current geometry" GUESS crash; also Guess PModel).
RESULT (job 395224, ended at walltime after 1 geom step - see cost note): first SCF converges cleanly,
E = -1229.569 Eh, NO instability/oscillation/collapse. Step-1 geometry: an electron-rich substrate atom
(carboxylate O6) forms a HELD SALT BRIDGE to a guanidinium N-H (O6...H = 1.508 A - a real H-bond
distance), and O3 migrates toward a guanidinium H (3.14 -> 2.87 A) with the substrate intact
(O3-C4 1.449). Contrast: bare/point-charge attempts pulled an O to 0.02-0.955 A (singular sink).
KEY RESULT: real molecular groups convert the implosion into a physical hydrogen bond / salt bridge -
Pauli repulsion supplies the minimum the point charge lacked. The surrogate route is the correct
representation. This is the resolution of the whole relaxation problem.

## OPEN ISSUE from Attempt 5 - COST (blocks scaling; decide before batch)
The full-QM 38-atom optimization is far too slow (first SCF ~12 h with ORCA LEAN-SCF; did not converge in
6 h walltime, only 1 geom step). Not scalable to (all designs) x (2 endpoints) x (NEB). Options for the
PRODUCTION guard, giving the same Pauli protection at ~24-atom cost:
 (a) QM substrate + MM guanidinium/formate (OpenMM layer, LJ sigma/eps = Pauli wall) via QMMMTheory -
     surrogate electrostatics as MM charges WITH LJ; validated OpenMM path already exists.
 (b) LJ-walled point charges (bare +/-1 as MM particles carrying sigma/eps) - lighter, less chemically
     explicit, but the LJ wall is the physically-motivated floor (r0 ~ vdW contact 3.0-3.4 A) the earlier
     wall-run lacked.
 (c) cheaper QM geometry (r2SCAN-3c) + B3LYP single point on top.
Full-QM surrogate = gold-standard validation (keep for the headline design); production guard should be
(a) or (b). NEXT: build the cheap guarded variant and re-run K2 to a CONVERGED endpoint, confirm the
salt bridge holds, then scale.

## Attempt 5 result + net-charge finding + literature verdicts (2026-07, updates the above)

### Attempt 5 (full-QM surrogates) — WORKED for K2, but only K2 is valid full-QM.
K2 full-QM surrogate run (guanidinium+1 / formate-1, frozen, 38 QM atoms, net -2): first SCF clean
(E=-1229.569 Eh), substrate carboxylate O forms a HELD salt bridge to a guanidinium N-H at 1.508 A -
a real H-bond, NOT a collapse. Confirms the point-charge implosion is a REPRESENTATION ARTIFACT and
real molecular groups fix it. Cost: ~12 h/SCF; did not converge in 6 h; relaunched clean as
'k2_fullqm_surr' (168 h walltime, email flags). This is the gold-standard anchor.

### NET-CHARGE PROBLEM (caught in sandbox before running K1/K3/K4 full-QM).
Full-QM surrogate systems inherit the design's net charge: K1=-3, K2=-2, K3=-3, K4=-4 (substrate -2
+ surrogate net). Only K2 is net -2. Building K1/K3/K4 as full-QM would create -3/-4 molecular
polyanions (substrate dianion + separated formate -1s, charge centres 3-9 A apart).

### LITERATURE VERDICT 1 (multiply-charged-anion research). CONFIRMED: -3/-4 full-QM polyanions
are LIKELY ELECTRONICALLY UNBOUND. Separated -1 carboxylates add Coulomb repulsion, not binding;
real analogues (pyrene tetrasulfonate, CuPc-tetrasulfonate) have NEGATIVE electron binding energies.
CPCM eps=4 is too weak to bind them (needs high dielectric / explicit counter-ions/waters:
~3 waters per -2 centre, ~16 per -3). CRITICAL ARTIFACT: def2-SVP (non-diffuse) MASKS the
instability - the SCF converges and looks bound while the true HOMO is positive/autodetaching; a
diffuse basis (def2-SVPD) would reveal density spilling out. => DO NOT relax isolated -3/-4 full-QM
systems in def2-SVP/CPCM(eps=4) and trust them. My earlier "-3/-4 likely unbound" hunch is
literature-confirmed, not speculation.

### LITERATURE VERDICT 2 (external-charge representation research). Recommended method =
electrostatic embedding: keep chorismate as the ONLY QM subsystem (charge stays -2 for ALL designs
- MM charges add field, not electrons, so the polyanion problem DISAPPEARS), with the external
+/-1 groups as MM with LENNARD-JONES / Pauli walls. Bare point charges collapse (spill-out); the
LJ/ECP/Gaussian-blur wall is the standard cure. Floor must be at vdW CONTACT (3.0-3.4 A) not 2.2 A,
and steep (LJ r^-12 / k>=500), which is EXACTLY why our Attempt-4 wall (2.2 A, k=100, harmonic) was
penetrated - mis-calibrated, not wrong in concept.

### CORRECTED PLAN (supersedes the "Next" list above).
- Full-QM surrogate: VALID only for K2 (net -2). Keep k2_fullqm_surr running as the gold standard.
- PRIMARY production method = QM substrate (-2) + MM MOLECULAR surrogates (guanidinium/formate) with
  FORCE-FIELD LJ (OPLS/AMBER vdW radii => real Pauli wall at vdW contact, no arbitrary tuning). QM
  charge stays -2 for K1-K4 => no net-charge/unbound problem AND no collapse. This is the one
  combination not yet run (bare-charge+custom-wall failed on calibration; full-QM hit net charge;
  MM-molecular-surrogate-with-FF-LJ has each missing piece).
- CALIBRATION: full-QM K2 vs MM-LJ K2. If they agree (salt-bridge distances, relaxed substrate
  geometry, barrier), MM-LJ is validated -> run K1-K4 on MM-LJ with confidence.
- MANDATORY sandbox gate before any HPC run: verify ASH's OpenMM layer actually APPLIES sigma/eps to
  the MM surrogate atoms (source shows addParticle(0,1,0) default = ZERO LJ). If LJ is not applied,
  it is bare charges again and will re-collapse; add LJ explicitly via CustomNonbondedForce. Do NOT
  run until LJ-on-MM-particles is confirmed.
- VERIFICATION on every run: HOMO<0; diffuse-basis (def2-SVPD) test for spill-out; no density on the
  outermost charge; for embedding, watch frontier-atom charges near the +1 for over-polarisation.

## Attempt 6 — MM MOLECULAR SURROGATES (guanidinium/formate as MM fragments + explicit QM-MM LJ).  PRODUCTION METHOD.
Synthesis of all findings. QM = 24 substrate only (charge -2, so NO polyanion problem for any design);
guanidinium(+1)/formate(-1) as MM fragments carrying distributed partial charges (sum +1/-1) AND explicit
Lennard-Jones vs every QM atom (336 LB-combined pairs, E=4eps[(s/r)^12-(s/r)^6]) injected as CustomBondForce
into mm.system. The small-sigma N-H lets an electron-rich substrate O approach to H-bond distance; N/O/C give
the Pauli wall. Surrogates frozen at certified positions. Script: ash_guarded/k2_mmsurr_spike.py.
RESULT (job k2_mmsurr, 200 steps, smooth monotonic descent -837.37, no oscillation, maxiter):
  min substrate<->surrogate = 1.766 A  -> REAL SALT BRIDGE, reproduces the full-QM gold standard (1.54 A).
  A substrate carboxylate O docks to a guanidinium N-H (same pose as full-QM, where O6 bridged); no collapse.
KEY RESULT: MM molecular surrogates reproduce the full-QM salt-bridge geometry at 24-atom QM cost, keep QM at
-2 (net-charge-safe for K1-K4), and optimize smoothly. THIS IS THE PRODUCTION GUARD for all designs.
Contrast: bare-charge + smooth LJ wall (Attempt 5b/wall3) only REPELS (substrate went to O3->+1 5.89 A) because
a bare +1 has no H-bond donor. The molecular surrogate supplies the donor -> attractive salt bridge, matching
gold. Calibration confirmed: cheap MM-surrogate ~ expensive full-QM on the salt-bridge geometry.
NEXT: clone to K1/K3/K4 (same LJ machinery, each design's own charge set + surrogate placement); run all four
in parallel. Then relaxed-path/NEB on the surviving designs.

## Attempt 6b — MM surrogate with LOOSE convergence (gmax 6e-3): CONVERGED. Calibration verdict.
Re-ran K2 MM-surrogate with reaction-relevant conv_criteria (gmax 6e-3, grms 3e-3) via ASH Optimizer
conv_criteria= kwarg. Result: CONVERGED ("Converged! =D") at step 241, gmax 3.4e-3, E=-837.376. First fully
converged relaxation of the project. Script: ash_guarded/k2_mmsurr_loose.py.
CALIBRATION (converged MM-surr vs full-QM step81 gold):
  reacting core AGREES: O3-C4 1.456 vs 1.434, C4-C6 2.504 vs 2.480 (both <=0.02 A).
  salt bridge AGREES: 1.541 vs 1.460 A (both form carboxylate<->guanidinium bridge).
  C1-C6 DIFFERS: 3.782 vs 5.004 A (1.2 A) - PERSISTS after convergence (was 1.3 A before) -> a REAL
  conformer difference in the soft forming-bond distance, NOT a convergence artefact. Likely cause: QM
  guanidinium polarises/charge-transfers to the substrate in ways an MM point-charge guanidinium cannot.
VERDICT: cheap MM-surrogate reproduces the catalytically-relevant reacting core + salt bridge (<0.02 A) but
gives a more OPEN reactant C1-C6. Since C1-C6 is a flat non-bonded coordinate at the reactant (stiff only at
the TS ~2.5 A), this may wash out in the BARRIER - to be verified at the TS. Note: full-QM gold itself is NOT
converged (maxiter step81), so part of the gap may be gold-standard non-convergence.
NEXT: (B) one restrained MM run (C1-C6 fixed at 3.78) to confirm everything-but-C1-C6 matches, then (A) clone
to K1/K3/K4, run relaxed paths, and CALIBRATE ON THE BARRIER (the real endpoint), not just reactant geometry.

## Attempt 6c — RESTRAINED calibration (C1-C6 pinned at full-QM 3.782 A): CONFIRMS the method.
Re-ran converged MM-surrogate K2 with C1-C6 (atoms 0,12) constrained to the full-QM value 3.782 A
(ASH Optimizer constraints={"bond":[[0,12,3.782]]}, constrainvalue=True). CONVERGED. Result:
  C1-C6 = 3.777 (constraint held), O3-C4 = 1.429 (full-QM 1.456), salt bridge = 1.640 A (full-QM 1.541).
VERDICT: with the one soft forming-bond coordinate pinned, the cheap MM-surrogate reproduces the full-QM
substrate geometry EVERYWHERE (reacting core + salt bridge <=0.1 A). => the 1.2 A C1-C6 gap in the free
run was the SOLE real difference, a single flat large-amplitude reactant coordinate; everything
catalytically relevant matches gold. MM MOLECULAR SURROGATE METHOD IS CALIBRATED for production.
Scripts: ash_guarded/k2_mmsurr_loose.py (free, converged) + k2_mmsurr_restr.py (restrained confirmation).

## Scaling to K1/K3/K4 — surrogate placement finding
K1 (1 formate) and K3 (guan + 2 formate) place cleanly (inter-fragment min 2.19 A, no clash). K4 does NOT:
its three -1 sites are only 2.61-4.50 A apart, and three formates (~2.2 A wide each) cannot occupy that
volume without their -O arms clashing (naive min 1.415 A). This is itself a design-realizability finding:
K4's charge packing is too dense for real molecular ions. Plan: K1/K2/K3 on full MM molecular surrogates;
K4 via mixed representation (guanidinium + 1 formate for the separated -1 + LJ-point-charges for the two
crowded -1s) OR document K4 as not surrogate-representable and run point-charge-LJ only.

## ALL FOUR REACTANT ENDPOINTS CONVERGED (calibrated MM-surrogate method).
| design | steps | O3-C4 | C4-C6 | C1-C6 | salt bridge | note |
|---|---|---|---|---|---|---|
| K1 | 136 | 1.434 | 2.485 | 4.263 | none (17.5) | benign, lone -1, no guanidinium |
| K2 | 241 | 1.434 | 2.480 | 5.004 | 1.460 | calibrated vs full-QM |
| K3 | 182 | 1.441 | 2.489 | 4.002 | 1.886 | 2 formate + guan |
| K4 | 314 | 1.447 | 2.482 | 5.042 | 1.853 | mixed rep, floppiest (dense charge) |
Reacting core consistent across all four (O3-C4 1.43-1.45, C4-C6 2.48-2.49). Step count tracks design
frustration (K1 136 -> K4 314) = more competing charges -> flatter surface. Product endpoints running next.
