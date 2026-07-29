#!/bin/bash
# RECONCILE_AND_COMMIT.sh — pull everything from HPC into the repo and stage it.
# Run on the Ubuntu PC from the repo root (~/Desktop/chorismate_thesis_code), branch tier1-realism.
# It rsyncs the HPC phase2b tree (scripts, logs, relaxed geoms, outputs) into the local repo, EXCLUDING
# only ORCA binary scratch (.gbw/.tmp/.densities) which are huge and regenerable. Keeps everything else
# incl. logs and _relaxed.xyz so results can be regenerated + rationale reconstructed.
set -e
HPC=18660916@hpc1.sun.ac.za
HROOT=/home/18660916/system_development/phase2b_charge_design
LROOT=phase2b_charge_design   # local repo path

echo "=== 1. pull the NEW step directories from HPC (scripts+logs+xyz+outputs, skip binary scratch) ==="
for d in 05b_optimise_realism 05c_cradle 05d_shell_ladder; do
  echo "--- $d ---"
  mkdir -p $LROOT/$d
  rsync -av \
    --exclude='*.gbw' --exclude='*.tmp' --exclude='*.densities' --exclude='*.densitiesinfo' \
    --exclude='*.cpcm' --exclude='*.cpcm_corr' --exclude='__pycache__' \
    --exclude='*.bibtex' --exclude='*.property.txt' \
    $HPC:$HROOT/$d/ $LROOT/$d/
done

echo "=== 2. pull updated 04_grid (extended-grid scripts + dv_grid_ext/fine + points) ==="
rsync -av \
  --include='*_ext.py' --include='*.tsv' --include='grid_final.xyz' --include='run_vpot*.sh' \
  --include='*.md' --exclude='*.npz' --exclude='*' \
  $HPC:$HROOT/04_grid/ $LROOT/04_grid/

echo "=== 3. pull any updated TIER2_KNOWN_ISSUES / DESIGN_DECISIONS from HPC (if edited there) ==="
for f in TIER2_KNOWN_ISSUES.md DESIGN_DECISIONS.md; do
  rsync -av $HPC:$HROOT/$f $LROOT/$f 2>/dev/null || true
done

echo "=== 4. drop the session notes + reconciliation doc into place ==="
cp PHASE2B_SESSION_RECONCILIATION.md $LROOT/
# (notes_this_session/*.md are already the per-step READMEs; place any that aren't in their step dirs)

echo "=== 5. stage everything ==="
git add $LROOT/
echo ""
echo "=== staged. review with: git status | head -50 ==="
echo "then: git commit -m '...' (see suggested message in PHASE2B_SESSION_RECONCILIATION.md) and git push"
