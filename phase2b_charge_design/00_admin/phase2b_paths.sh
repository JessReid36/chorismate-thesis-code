#!/usr/bin/env bash
# Central path config for phase2b. Source before any stage script.
# HPC working tree; geometry sourced from the validated Phase 1 reduced-region path.
export PHASE2B_ROOT="$HOME/system_development/phase2b_charge_design"
export QMMM_ROOT="$HOME/system_development/05_qmmm"

# Validated consistent-region geometry sources (last frame of each QM-region trajectory):
export TS_SOURCE="$QMMM_ROOT/18c_reduced_region/ts_reduced.QMRegion_trj.xyz"
export REACTANT_SOURCE="$QMMM_ROOT/18e_reactant_reduced/reactant_reduced.QMRegion_trj.xyz"
export PRODUCT_SOURCE="$QMMM_ROOT/18f_product_reduced/product_reduced.QMRegion_trj.xyz"

# Extracted substrate geometries (written by 01_geometry):
export GEOM_DIR="$PHASE2B_ROOT/01_geometry"
export SUBSTRATE_REACTANT="reactant.xyz"
export SUBSTRATE_TS="ts.xyz"
export SUBSTRATE_PRODUCT="product.xyz"

# Reacting-atom map (0-based): C1=0 O3=7 C4=8 C6=12
export REACTING_MAP="$PHASE2B_ROOT/00_admin/reacting_atoms.tsv"
