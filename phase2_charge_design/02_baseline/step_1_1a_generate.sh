#!/usr/bin/env bash
# Phase 2 step 1.1a - generate reactant single-point inputs: vacuum diagnostic + CPCM(eps=4)
# production. Vacuum run is expected to expose the metastable-dianion pathology (positive/
# near-zero HOMO); CPCM eps=4 (Phase 1 QM-cluster convention) binds the density. KeepDens
# retains .scfp for orca_vpot. Identical CPCM at TS => continuum cancels in Dv = V_TS - V_R.
set -euo pipefail
source "$HOME/system_development/phase2_charge_design/00_admin/phase2_paths.sh"

stage="$PHASE2_ROOT/02_baseline"
mkdir -p "$stage"
react="$SUBSTRATE_DIR/$SUBSTRATE_REACTANT"
[[ -s "$react" ]] || { echo "FAIL: reactant geometry missing: $react"; exit 1; }

coords="$(tail -n +3 "$react")"
[[ -n "$coords" ]] || { echo "FAIL: no coordinates parsed from $react"; exit 1; }

cat > "$stage/react_vac.inp" <<INP
! B3LYP D3BJ def2-SVP def2/J RIJCOSX KeepDens
%maxcore 3000
%pal nprocs 8 end
%scf MaxIter 300 end
* xyz $SUBSTRATE_CHARGE $SUBSTRATE_MULT
$coords
*
INP

cat > "$stage/react_cpcm.inp" <<INP
! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM KeepDens
%maxcore 3000
%pal nprocs 8 end
%cpcm epsilon 4.0 end
%scf MaxIter 300 end
* xyz $SUBSTRATE_CHARGE $SUBSTRATE_MULT
$coords
*
INP

echo "wrote:"
echo "  $stage/react_vac.inp"
echo "  $stage/react_cpcm.inp"
