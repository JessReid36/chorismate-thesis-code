#!/usr/bin/env bash
# Phase 2 step 1.1b - production single-points for the Dv engine: reactant AND TS at
# CPCM(eps=4), identical settings, KeepDens for orca_vpot. Reference-state decision (CPCM
# eps=4 binds the dianion; vacuum leaves it metastable) was established and committed in 1.1a.
# HOMO checked on both; the TS is more charge-separated, so we verify it is also bound.
set -euo pipefail
source "$HOME/system_development/phase2_charge_design/00_admin/phase2_paths.sh"

stage="$PHASE2_ROOT/03_dvpot"
mkdir -p "$stage"

emit () {  # $1 = geometry file, $2 = output basename
  local geo="$1" base="$2"
  [[ -s "$geo" ]] || { echo "FAIL: geometry missing: $geo"; exit 1; }
  local coords; coords="$(tail -n +3 "$geo")"
  [[ -n "$coords" ]] || { echo "FAIL: no coordinates parsed from $geo"; exit 1; }
  cat > "$stage/${base}.inp" <<INP
! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM KeepDens
%maxcore 3000
%pal nprocs 8 end
%cpcm epsilon 4.0 end
%scf MaxIter 300 end
* xyz $SUBSTRATE_CHARGE $SUBSTRATE_MULT
$coords
*
INP
  echo "  wrote $stage/${base}.inp  (from $(basename "$geo"))"
}

echo "generating production single-point inputs:"
emit "$SUBSTRATE_DIR/$SUBSTRATE_REACTANT" "sp_reactant"
emit "$SUBSTRATE_DIR/$SUBSTRATE_TS"       "sp_ts"
