# Relaxation Attempts — Job & File Manifest (phase2b endpoint relaxation)

Single source of truth mapping every HPC job -> script -> working dir -> outcome -> keep/archive.
Working tree root (HPC): ~/system_development/phase2b_charge_design/07_stage0/stage1/
Companion: RELAXATION_MASTER_LOG.md (insight/do-not-redo narrative). This file is the index.

## Job ledger
| Job ID | Name | Wall | Exit | Method / script | Dir | Outcome | Keep? |
|---|---|---|---|---|---|---|---|
| array | s1_conv | ~4.5h ea | maxiter | native bare + Calc_Hess (Stage 1) | stage1/runs/ | 0/8 conv; implosion K2-K4 | committed (results 07_stage0/stage1) |
| 395168 | k2_guarded | - | - | ASH fixed-pin O3->+1 | stage1/ | pin holds O3, O5 collapses 0.955A | KEEP pin_run |
| wall-k100 | k2_wall | - | maxiter | ASH harmonic wall 2.2A/k100 | stage1/ | penetrated, C3 @1.26A | KEEP wall_run_k100 |
| 395254 | k2_fullqm_surr | 4:49 | maxiter | full-QM surrogates guan+formate | surrogate/ | 81 steps, salt bridge O6..H 1.54A HELD; net -2 (K2 only) | KEEP GOLD STANDARD |
| 395255 | k2w2_dry | 8s | 0 | step-wall dry-run | wall2/ | ok(dry) | archive |
| 395256 | k2_mmlj_wall2 | 1m | fail | step-wall (step()) | wall2/ | early fail dup | archive |
| 395257 | k2_mmlj_wall2 | 2m | killed | step-wall dup | wall2/ | qdel | archive |
| 395258 | k2_mmlj_wall2 | 37s | killed | step-wall dup | wall2/ | qdel | archive |
| 395259 | k2_mmlj_wall2 | 2:16 | maxiter | step-wall floor3.2 CLIFF | wall2/ | OSCILLATED across step() discontinuity; held 3.199A but un-optimizable | KEEP instructive |
| 395364 | k2w3_dry | 5s | 0 | smooth-wall dry-run | wall3/ | ok(dry) | archive |
| 395365 | k2_mmlj_wall3 | 3:19 | maxiter | SMOOTH wall eps*(sigma/r)^12 sigma3.8 | wall3/ | 200 steps NO oscillation, held 3.41A; substrate relaxed AWAY O3->+1 5.89A | KEEP |

## Scripts (HPC stage1/, committed copies in relaxation_attempts/ash_guarded/)
- surrogate/k2_surrogate_spike.py -- full-QM surrogates (frozen, net -2). GOLD STANDARD.
- wall2/k2_wall2_spike.py         -- step-wall. SUPERSEDED (discontinuity -> oscillation).
- wall3/k2_wall3_spike.py         -- smooth wall eps*(sigma/r)^12 sigma=3.8. Collapse-safe; bare +1 = no H-bond donor.

## KEY CALIBRATION FINDING (2026-07-24)
Full-QM surrogate (gold) vs cheap smooth-wall bare-point-charge, K2 reactant:
- Full-QM: substrate moves TOWARD the +1, forms SALT BRIDGE (carboxylate O..guanidinium N-H 1.54 A).
- Cheap wall (bare +1 + LJ): substrate moves AWAY (O3->+1 5.89 A); a bare +1 offers no H-bond donor.
=> They DISAGREE qualitatively. Bare-charge+wall is collapse-safe and net-charge-safe but models the WRONG
   interaction (repulsive standoff, not attractive salt bridge). Too crude to reproduce the gold standard.
=> NEXT: MM MOLECULAR surrogates (guanidinium/formate as MM fragments WITH LJ) - keeps QM at -2 (cheap,
   net-charge-safe) AND supplies the N-H donor geometry -> should reproduce the full-QM salt bridge.

## STANDING TECHNICAL NOTE
All three guard methods HOLD geometry but none trip tight default geomeTRIC convergence in the step budget
(exit 1 = GeomOptNotConverged each) -> field-distorted surface is hard to converge tightly. Lever = convergence
criterion (raise maxiter / loosen to reaction-relevant), NOT more wall tuning. For calibration compare
best-reached relaxed geometry (as done for full-QM step-81).
