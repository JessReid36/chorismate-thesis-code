#!/usr/bin/env bash
# Phase 2 step 0.1 - scaffold the charge-design tree and define the single substrate
# swap point. Idempotent: only regenerates its own bookkeeping, never deletes stage work.
set -euo pipefail

phase2_root="$HOME/system_development/phase2_charge_design"
admindir="$phase2_root/00_admin"
substrate_dir="$phase2_root/01_substrate"

parent="$HOME/system_development"
[[ -d "$parent" ]] || { echo "FAIL: expected Phase 1 parent tree at $parent (edit phase2_root if your layout differs)"; exit 1; }

subdirs=(00_admin 01_substrate 02_baseline 03_dvpot 04_grid 05_select 06_polarize 07_polish 08_validate 09_benchmark local_workstation)
for d in "${subdirs[@]}"; do
  mkdir -p "$phase2_root/$d"
done

paths_file="$admindir/phase2_paths.sh"
metadata_template="$substrate_dir/substrate_metadata.tsv"
layout_doc="$admindir/phase2_layout.md"
audit="$admindir/step00_1_scaffold_audit.txt"
checksum="$admindir/sha256_step00_1.txt"

rm -f "$paths_file" "$layout_doc" "$audit" "$checksum"

cat > "$paths_file" <<PATHS
# Phase 2 canonical paths - source at the top of every Phase 2 script.
# Changing the reaction (placeholder -> final Phase 1 NEB triple) = overwrite the
# three files in \$SUBSTRATE_DIR and re-run. No downstream script hard-codes a geometry path.
export PHASE2_ROOT="$phase2_root"
export SUBSTRATE_DIR="$phase2_root/01_substrate"
export SUBSTRATE_REACTANT="reactant.xyz"
export SUBSTRATE_TS="ts.xyz"
export SUBSTRATE_PRODUCT="product.xyz"
export SUBSTRATE_CHARGE="-2"
export SUBSTRATE_MULT="1"
PATHS

# printf (not a heredoc) so literal tabs are generated locally and survive any paste.
if [[ ! -s "$metadata_template" ]]; then
  {
    printf '# Substrate geometry provenance. One row per endpoint; fill in at step 0.2.\n'
    printf '# r = antisymmetric reaction coordinate d(C4-O3 breaking) - d(C6-C1 forming), Angstrom.\n'
    printf 'endpoint\tfilename\tsource\tcharge\tmult\tr_angstrom\tnotes\n'
    printf 'reactant\treactant.xyz\tPLACEHOLDER_cp2k\t-2\t1\tNA\tplaceholder-geometry\n'
    printf 'ts\tts.xyz\tPLACEHOLDER_cp2k\t-2\t1\tNA\tplaceholder-geometry\n'
    printf 'product\tproduct.xyz\tPLACEHOLDER_cp2k\t-2\t1\tNA\tplaceholder-geometry\n'
  } > "$metadata_template"
fi

cat > "$layout_doc" <<DOC
# Phase 2 layout - charge-design pipeline
Root: \`$phase2_root\` (sibling of phase1_system_dev)

## Swap point
\`01_substrate/\` holds reactant/TS/product geometries - the only reaction-specific input.
Paths/filenames defined once in \`00_admin/phase2_paths.sh\`, sourced by every downstream
script. Finalising Phase 1 = overwrite the three .xyz files and re-run; no code edits.

## Stage directories
00_admin, 01_substrate (SWAP POINT), 02_baseline, 03_dvpot, 04_grid, 05_select,
06_polarize, 07_polish, 08_validate, 09_benchmark, local_workstation.
Placeholder-derived results carry [placeholder-geometry] until the final Phase 1 triple lands.
DOC

{
  echo "phase2 root: $phase2_root"
  echo "subdirs verified: ${#subdirs[@]}"
  echo "swap point: $substrate_dir"
  echo "generated (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$audit"

sha256sum "$paths_file" "$metadata_template" "$layout_doc" > "$checksum"
echo "step 0.1 complete"
