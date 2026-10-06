#!/usr/bin/env bash
# x25_c5_prepare.sh - RUN ON hpc1 (login node: file copies and text edits only, no ORCA).
#
# PHASE1_AUDIT_CHECKLIST.md item C5. The eight post-cut frames that ran NEB-TS took their barrier
# from a band halted at the NEB-TS thresholds (climbing image 2.0e-3 Eh/bohr). This prepares one
# NEB-CI restart per frame, from that frame's final band (neb_MEP.allxyz), converged to the NEB-CI
# defaults used by the other 22 post-cut frames (climbing image 5.0e-4 Eh/bohr).
#
# Everything is written to a NEW stage folder, 05_qmmm/22_c5_nebci/frame_<N>/. The original
# 19_ensemble folders are read, never written. Nothing is submitted: the script prints the qsub
# commands to run after you have looked at one input.
#
# Each input is the frame's own neb.inp (same active region, endpoints and settings) with exactly
# two changes: NEB-TS -> NEB-CI, and a Restart_ALLXYZFile line. The script then checks that the
# result, with the restart line removed and the active list masked, is byte-identical to the
# NEB-CI template frame's neb.inp - i.e. that these restarts run under exactly the settings of the
# 22 NEB-CI frames.
#
# No heredocs: they fail on hpc1. Files are built with sed and printf.
set -euo pipefail
SD=/home/18660916/system_development/05_qmmm
ENS=$SD/19_ensemble
OUT=$SD/22_c5_nebci
TEMPLATE=36665                       # an NEB-CI frame whose neb.inp defines the settings to match
FRAMES="20000 21634 23268 24883 26495 33320 34991 43738"
PRMS_SHA=3181d624af818b668375a8dc4dd137708eeb8f80cc2ae616767334ae1ace4980   # committed 13_bridge copy
NATOM=55680
stop() { echo "STOP: $*"; exit 1; }

[ -e "$OUT" ] && stop "$OUT exists; remove it first"
mask() { sed -E 's/ActiveAtoms \{[^}]*\}/ActiveAtoms {..}/' "$1"; }
grep -q ' NEB-CI ' "$ENS/frame_$TEMPLATE/neb.inp" || stop "template $TEMPLATE is not an NEB-CI input"

# ---- pre-flight: every frame must have what the restart needs, before anything is written
for f in $FRAMES; do
  F=$ENS/frame_$f
  grep -q ' NEB-TS ' "$F/neb.inp" || stop "frame $f: neb.inp is not NEB-TS"
  for x in neb.inp neb.pbs reactant.pdb product.pdb tsguess.pdb complex_solvated.ORCAFF.prms neb_MEP.allxyz; do
    [ -s "$F/$x" ] || stop "frame $f: missing $x"
  done
  [ "$(sha256sum "$F/complex_solvated.ORCAFF.prms" | cut -c1-64)" = "$PRMS_SHA" ] || stop "frame $f: force field differs from the committed one"
  n1=$(head -1 "$F/neb_MEP.allxyz" | tr -d ' ')
  [ "$n1" = "$NATOM" ] || stop "frame $f: neb_MEP.allxyz first block has $n1 atoms, expected $NATOM"
  nb=$(( $(grep -c '^>' "$F/neb_MEP.allxyz") + 1 ))
  [ "$nb" -eq 10 ] || stop "frame $f: neb_MEP.allxyz holds $nb structures, expected 10 (8 images + 2 end points)"
done
echo "pre-flight passed for all eight frames"

mkdir -p "$OUT"
for f in $FRAMES; do
  F=$ENS/frame_$f
  D=$OUT/frame_$f
  mkdir "$D"
  cp -p "$F/reactant.pdb" "$F/product.pdb" "$F/tsguess.pdb" "$F/complex_solvated.ORCAFF.prms" "$D/"
  cp -p "$F/neb_MEP.allxyz" "$D/restart.allxyz"
  # input: the frame's own neb.inp, NEB-TS -> NEB-CI, plus the restart line inside %neb
  sed -e '1s/ NEB-TS / NEB-CI /' -e '/^%neb$/a\  Restart_ALLXYZFile "restart.allxyz"' "$F/neb.inp" > "$D/neb.inp"
  [ "$(grep -c 'Restart_ALLXYZFile' "$D/neb.inp")" -eq 1 ] || stop "frame $f: restart line not inserted exactly once"
  grep -q ' NEB-TS ' "$D/neb.inp" && stop "frame $f: NEB-TS still present"
  # same settings as the NEB-CI frames: identical once the restart line is removed and the list masked
  cmp -s <(grep -v 'Restart_ALLXYZFile' "$D/neb.inp" | sed -E 's/ActiveAtoms \{[^}]*\}/ActiveAtoms {..}/') <(mask "$ENS/frame_$TEMPLATE/neb.inp") \
    || { diff <(grep -v 'Restart_ALLXYZFile' "$D/neb.inp" | sed -E 's/ActiveAtoms \{[^}]*\}/ActiveAtoms {..}/') <(mask "$ENS/frame_$TEMPLATE/neb.inp"); stop "frame $f: settings differ from NEB-CI frame $TEMPLATE"; }
  # and the same active region as this frame's original run
  cmp -s <(grep -o 'ActiveAtoms {[^}]*}' "$D/neb.inp") <(grep -o 'ActiveAtoms {[^}]*}' "$F/neb.inp") || stop "frame $f: active region changed"
  # job script: the frame's own neb.pbs with a new name, folder, walltime and restart-file check
  sed -e "s|^#PBS -N cm19_n$f\$|#PBS -N cm22_c5_$f|" \
      -e 's|^#PBS -l walltime=.*$|#PBS -l walltime=168:00:00|' \
      -e "s|/19_ensemble/frame_$f|/22_c5_nebci/frame_$f|g" \
      -e 's|for f in reactant.pdb product.pdb; do|for f in reactant.pdb product.pdb restart.allxyz; do|' \
      "$F/neb.pbs" > "$D/neb.pbs"
  grep -q "^#PBS -N cm22_c5_$f\$" "$D/neb.pbs" || stop "frame $f: job name not set"
  grep -q "19_ensemble" "$D/neb.pbs" && stop "frame $f: job script still points into 19_ensemble"
  [ "$(grep -c "/22_c5_nebci/frame_$f" "$D/neb.pbs")" -eq 2 ] || stop "frame $f: job script should name its new folder twice (-o and cd)"
  grep -q 'restart.allxyz' "$D/neb.pbs" || stop "frame $f: restart-file check not added"
done

# ---- the pass criteria, written before anything runs (PHASE1_AUDIT_CHECKLIST.md C5)
{
  printf '%s\n' "C5 pass criteria, fixed $(date -Iseconds) before submission"
  printf '%s\n' "1. Final CI-NEB convergence table meets the NEB-CI defaults: MAX(|FCI|) <= 5.0e-4, RMS(FCI) <= 2.5e-4,"
  printf '%s\n' "   MAX(|Fp|) <= 5.0e-3, RMS(Fp) <= 2.5e-3 Eh/bohr, and THE NEB OPTIMIZATION HAS CONVERGED is printed."
  printf '%s\n' "2. The barrier varies by at most 0.1 kcal/mol over the last 5 iterations of neb.NEB.log."
  printf '%s\n' "3. Expected: the barrier falls by 0.2-0.8 kcal/mol from its current value. A fall above 2 kcal/mol is"
  printf '%s\n' "   inspected for environment rearrangement (C4) before it is accepted; a rise is inspected too."
  printf '%s\n' "4. The new climbing image lies within about 0.05 A RMSD (QM region) of the current one."
  printf '%s\n' "5. Image 0 and the last image keep their energies (fixed end points): within 1e-5 Eh."
} > "$OUT/C5_CRITERIA.txt"

P=$SD/PROVENANCE.log
[ -f "$P" ] || stop "no $P"
echo "$(date -Iseconds)  22_c5_nebci  prepared NEB-CI restarts (C5) for $FRAMES from each frame's neb_MEP.allxyz; settings identical to NEB-CI frame $TEMPLATE apart from the restart line (x25)" >> "$P"

echo "== written:"
ls "$OUT" "$OUT/frame_20000"
echo "== frame 20000 input, active list masked:"
mask "$OUT/frame_20000/neb.inp"
echo "== frame 20000 job script:"
cat "$OUT/frame_20000/neb.pbs"
echo "== restart band comment lines, frame 20000 (energies of the 10 structures, if ORCA wrote them):"
awk 'p==1{print; p=0} /^ *55680 *$/{p=1}' "$OUT/frame_20000/restart.allxyz" | cut -c1-100
echo
echo "Look these over. To submit all eight:"
echo "  for f in $FRAMES; do qsub $OUT/frame_\$f/neb.pbs; done"
