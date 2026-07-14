#!/usr/bin/env bash
# Phase 2 step 4.1 - polarization test for the K=2 net-neutral design. Four single-points:
# {reactant, TS} x {bare, with-charges}, B3LYP-D3BJ/def2-SVP, CPCM eps=4, KeepDens. Designed charges
# enter via %pointcharges (.pc, Angstrom), inside the eps=4 cavity with the substrate. Polarized
# barrier change = difference-of-differences (charge-charge self-energy cancels, positions fixed
# at R and TS); compared to the frozen-density q.Dv prediction.
set -euo pipefail
source "$HOME/system_development/phase2_charge_design/00_admin/phase2_paths.sh"
stage="$PHASE2_ROOT/06_polarize"
mkdir -p "$stage"; cd "$stage"
pc="design_K2_neutral.pc"
[[ -s "$pc" ]] || { echo "FAIL: $pc not found"; exit 1; }

emit () {
  local geo="$1" base="$2" withpc="$3"
  local coords; coords="$(tail -n +3 "$geo")"
  {
    echo "! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM KeepDens"
    echo "%maxcore 3000"
    echo "%pal nprocs 8 end"
    echo "%cpcm epsilon 4.0 end"
    echo "%scf MaxIter 300 end"
    [[ "$withpc" == "yes" ]] && echo "%pointcharges \"$pc\""
    echo "* xyz $SUBSTRATE_CHARGE $SUBSTRATE_MULT"
    echo "$coords"
    echo "*"
  } > "$stage/${base}.inp"
  echo "  wrote ${base}.inp (with_pc=$withpc)"
}
echo "generating four single-point inputs:"
emit "$SUBSTRATE_DIR/$SUBSTRATE_REACTANT" pol_R_bare no
emit "$SUBSTRATE_DIR/$SUBSTRATE_TS"       pol_TS_bare no
emit "$SUBSTRATE_DIR/$SUBSTRATE_REACTANT" pol_R_pc  yes
emit "$SUBSTRATE_DIR/$SUBSTRATE_TS"       pol_TS_pc yes
