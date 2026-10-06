#!/usr/bin/env bash
# wp4d_cartesian.sh - RUN ON hpc1 (login node: copies and edits only). WP4 addendum.
#
# WHY: in the free variants ORCA optimises in redundant internal coordinates, which carry no overall
# translation or rotation. Whether the substrate can still move towards a fixed site then depends on
# whether ORCA links that site into its coordinate set: it did in s13 (site 3.2 A from O3: "Will
# constrain atom 24", B(Na 24, O 7)), it did not in WP4's armb_free (site 3.55 A away), where the
# substrate's centroid stayed put to 1e-6 A. So "free" was decided by a connectivity heuristic.
# The ORCA 6 manual: "coordsys ... can be set to cartesian or redundant"; "ProjectTR false" is the
# default and "MUST be false for redundant internals". In Cartesian coordinates with ProjectTR false,
# translation and rotation are ordinary degrees of freedom, for every arm alike.
#
# WHAT: three new folders, armX_freecart for X = a, b, c, next to the existing ones (nothing existing is
# touched). Each copies its armX_free combined.xyz and site.ORCAFF.prms byte for byte, and its job.inp
# with one change - the %geom block gains "coordsys cartesian" and "ProjectTR false", and MaxIter rises
# to 1000 because Cartesian steps can need more cycles (it does not change a converged result).
# Job names cm24_wp4_armX_freecart. Nothing is submitted.
set -euo pipefail
W=${W:-/home/18660916/system_development/phase2.2/wp4_s13_rerun}
stop() { echo "STOP: $*"; exit 1; }
[ -d "$W" ] || stop "$W not found"
for X in a b c; do
  S=$W/arm${X}_free; D=$W/arm${X}_freecart
  [ -d "$S" ] || stop "$S not found"
  [ -e "$D" ] && stop "$D exists; remove it first"
done
for X in a b c; do
  S=$W/arm${X}_free; D=$W/arm${X}_freecart
  mkdir "$D"
  cp -p "$S/combined.xyz" "$S/site.ORCAFF.prms" "$D/"
  cmp -s "$S/combined.xyz" "$D/combined.xyz" && cmp -s "$S/site.ORCAFF.prms" "$D/site.ORCAFF.prms" || stop "copy of arm $X inputs differs"
  grep -q '^%geom MaxIter 500 end$' "$S/job.inp" || stop "arm $X job.inp has no single-line %geom block to edit"
  sed 's/^%geom MaxIter 500 end$/%geom\n  MaxIter 1000\n  coordsys cartesian\n  ProjectTR false\nend/' "$S/job.inp" > "$D/job.inp"
  # the only difference must be that block
  n=$(diff "$S/job.inp" "$D/job.inp" | grep -c '^[<>]' || true)
  [ "$n" -eq 6 ] || { diff "$S/job.inp" "$D/job.inp"; stop "arm $X: unexpected input differences"; }
  sed -e "s|^#PBS -N cm24_wp4_arm${X}_free\$|#PBS -N cm24_wp4_arm${X}_freecart|" \
      -e "s|/arm${X}_free/|/arm${X}_freecart/|g" -e "s|/arm${X}_free\$|/arm${X}_freecart|g" "$S/job.pbs" > "$D/job.pbs"
  grep -q "^#PBS -N cm24_wp4_arm${X}_freecart\$" "$D/job.pbs" || stop "arm $X: job name not set"
  [ "$(grep -c "/arm${X}_freecart" "$D/job.pbs")" -eq 2 ] || stop "arm $X: job script should name its new folder twice (-o and cd)"
  grep -q -E "/arm${X}_free(/|\$)" "$D/job.pbs" && stop "arm $X: job script still points at the internal-coordinate run"
done

cat > "$W/WP4_CRITERIA_CARTESIAN.txt" <<'EOF'
WP4 addendum, fixed before submission: the free variants rerun in Cartesian coordinates.

The armX_free runs optimised in redundant internal coordinates, which contain no overall translation or
rotation; in armb_free the point site was not linked into the coordinate set and the substrate's
centroid did not move (1e-6 A). Those runs are kept and reported as "internal motion only". The
armX_freecart runs start from byte-identical files and differ only in the optimisation coordinates
(Cartesian, ProjectTR false, MaxIter 1000), so the substrate may translate and rotate in every arm.
WP4's conclusions on the representation rest on the anchored variants; these runs show what a lone
site or group does to a substrate that is free to move as a whole.

Expectations (same rigid model as WP4_CRITERIA.txt, which allows exactly this whole-body motion):
  arma_freecart  bare +1: collapse onto a substrate oxygen (closest O below 1.2 A).
  armb_freecart  N-sized LJ site: 2.584 A to O2 with O3 at 2.78 A; expected 2.28-2.88 A to O2 or O3.
                 Compared with armb_free (2.563 A to O2, O3 at 3.46 A, centroid fixed).
  armc_freecart  methylguanidinium: informative. With the group fully frozen the rigid model turns the
                 substrate by about 230 degrees and 5.5 A to a double-carboxylate bridge; here the group
                 can turn about its attachment. Reported: binding mode, contacts against the D-WP4
                 windows, centroid displacement and rotation.
Every run: converged, or the contact distance moves by at most 0.05 A over the last third of the cycles.
EOF
echo "== written:"; ls -d "$W"/arm*_freecart
echo "== arm b input, original against Cartesian:"; diff "$W/armb_free/job.inp" "$W/armb_freecart/job.inp" || true
echo "== arm b job script header:"; sed -n '2p;8,9p' "$W/armb_freecart/job.pbs"
echo "Submit with qsub after checking the queue (see the chat for the chained commands)."
