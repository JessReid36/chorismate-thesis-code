#!/usr/bin/env bash
# Phase 2 step 4.1 (net-free) - polarization test for the K=2 net-FREE design (two -1 anions, net -2).
# Same four-single-point scheme as net-neutral. This design placed anions on the extreme-positive-Dv
# innermost-shell points; the hypothesis is that its frozen-density edge leans on the model's least
# reliable regime and will deflate more than the net-neutral design under self-consistent relaxation.
set -euo pipefail
source "$HOME/system_development/phase2_charge_design/00_admin/phase2_paths.sh"
stage="$PHASE2_ROOT/06_polarize"
mkdir -p "$stage"; cd "$stage"
pc="design_K2_free.pc"
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
echo "generating four single-point inputs (net-free):"
emit "$SUBSTRATE_DIR/$SUBSTRATE_REACTANT" polf_R_bare no
emit "$SUBSTRATE_DIR/$SUBSTRATE_TS"       polf_TS_bare no
emit "$SUBSTRATE_DIR/$SUBSTRATE_REACTANT" polf_R_pc  yes
emit "$SUBSTRATE_DIR/$SUBSTRATE_TS"       polf_TS_pc yes
