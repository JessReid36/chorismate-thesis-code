#!/bin/bash
# get the reactant-under-8A-field final energy from the completed ladder job log
LOG=~/system_development/phase2b_charge_design/05d_shell_ladder/relax/run_shell8p0/ladder_8p0_relax.log
echo "=== final energy lines in shell8p0 (reactant under 8A field) ==="
grep -iE "FINAL SINGLE POINT ENERGY|Current energy|final energy|Energy:" "$LOG" | tail -5
echo "(use the LAST converged 'FINAL SINGLE POINT ENERGY' as E_reactant)"
