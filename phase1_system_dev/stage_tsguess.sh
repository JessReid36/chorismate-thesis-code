#!/usr/bin/env bash
# stage_tsguess.sh - copy the highest-energy window of a completed scan to
# tsguess.pdb, for use as the transition-state guess of the band calculation.
#
# WHY THIS IS A SEPARATE SCRIPT
#   The band needs three structures: the two optimised endpoints and a guess at
#   the saddle. The scan produces all three, and its maximum is where the
#   reference calculation's guess came from: 18_nebts/tsguess.pdb is
#   byte-identical to 16_scan/win_11.pdb, and window 11 is that scan's energy
#   maximum.
#
#   Keeping this in its own file, copied into each frame directory the way
#   run_scan.sh already is, means there is one definition of how the guess is
#   chosen rather than one per generated job. The scan stage did not drift
#   between the reference and the ensemble because it copies its driver; the
#   band stage did drift because its input was written from scratch.
#
# USAGE
#   Run from a frame's scan directory, after the scan has completed:
#       bash stage_tsguess.sh [output_path]
#   The default output is ../tsguess.pdb, which is where the band input expects
#   it.
#
# EXIT STATUS
#   0 on success, 1 if no maximum can be identified or the window file is
#   missing. The job that calls this should stop on a non-zero status rather
#   than leave the band to start from nothing.

set -uo pipefail
OUT="${1:-../tsguess.pdb}"

if [[ ! -s scan_energies.tsv ]]; then
  echo "TSGUESS_MISSING: no scan_energies.tsv"
  exit 1
fi

# scan_energies.tsv columns: window, achieved reaction coordinate, energy in Eh.
# The maximum energy is the least negative, so a plain numeric sort on column
# three and taking the last row gives the window at the top of the profile.
read -r W E < <(sort -k3 -g scan_energies.tsv | tail -1 | awk '{print $1, $3}')

if [[ -z "${W:-}" ]]; then
  echo "TSGUESS_MISSING: could not identify a maximum in scan_energies.tsv"
  exit 1
fi

if [[ ! -s "win_${W}.pdb" ]]; then
  echo "TSGUESS_MISSING: window $W has no win_${W}.pdb"
  exit 1
fi

cp "win_${W}.pdb" "$OUT" || { echo "TSGUESS_MISSING: copy failed"; exit 1; }

# Report the guess against the first window, which is the reactant end of the
# scan. A guess that sits only a few kcal/mol above the reactant is a poor
# starting point for the band, and one frame of the first ensemble whose guess
# was 7.18 kcal/mol above its reactant produced a band that relaxed to a
# monotonically decreasing profile and found no saddle.
E1=$(sort -k1 -n scan_energies.tsv | head -1 | awk '{print $3}')
if [[ -n "${E1:-}" ]]; then
  REL=$(awk -v a="$E" -v b="$E1" 'BEGIN{printf "%.2f", (a-b)*627.5094740631}')
  echo "TSGUESS window $W, $REL kcal/mol above the first scan window"
  LOW=$(awk -v r="$REL" 'BEGIN{print (r<8.0) ? 1 : 0}')
  if [[ "$LOW" == "1" ]]; then
    echo "TSGUESS_LOW: the guess is less than 8 kcal/mol above the scan start."
    echo "  The band may relax back into the reactant basin rather than climbing."
    echo "  Check the converged path for a strictly decreasing profile."
  fi
else
  echo "TSGUESS window $W"
fi
exit 0
