#!/usr/bin/env bash
# run_ensemble_batched.sh - run the ensemble in stages, keeping disk use bounded.
#
# NOTHING ABOUT THE METHOD CHANGES. The same generators, the same inputs, the
# same convergence criteria. Only the order and concurrency of submission, and a
# harvest-then-clean step between batches of the expensive stage.
#
# WHERE THE SPACE GOES
#   A finished frame holds about 5.7 GB, of which 4.78 is neb.appr.hess and 0.65
#   is neb_MEP_ALL_trj.xyz. Everything before the band stage is small: the
#   reactant optimisation, the scan and the product optimisation together come to
#   roughly 400 MB.
#
#   The Hessian is built during the transition-state optimisation that NEB-TS
#   performs after the band has converged. In frame 02450 the band converged at
#   line 7902 of the output and the transition-state optimisation began at 18630,
#   so more than half the run happens after the barrier is already determined.
#   The barrier is read from the path summary written at convergence.
#
#   Peak disk is therefore set by how many band calculations run at once, not by
#   how many frames there are in total.
#
# THE STAGES
#   1  reactant optimisation   all frames at once, about 12 GB for thirty
#   2  scan                    all frames at once, needs stage 1 per frame
#   3  product optimisation    all frames at once, needs stage 2 per frame
#   4  band                    in batches, about 5.7 GB per concurrent frame
#
#   After each batch of stage 4: harvest, then remove the heavy files from the
#   frames just harvested. Their path summary, barrier and converged geometries
#   are already copied into 19_ensemble_barriers by the harvest, at a few MB per
#   frame.
#
# WHAT IS REMOVED AFTER HARVEST
#   neb.appr.hess, neb.gu.tmp, and the full-system trajectories. The
#   quantum-region trajectory of the same path is kept, as are every converged
#   geometry, the path summary, and every input and job script.
#
# Usage:
#   bash run_ensemble_batched.sh status
#   bash run_ensemble_batched.sh stage1
#   bash run_ensemble_batched.sh stage2
#   bash run_ensemble_batched.sh stage3
#   bash run_ensemble_batched.sh stage4 [batch size, default 4]
#   bash run_ensemble_batched.sh harvest     # harvest and clean what has finished

set -uo pipefail
ROOT="$HOME/system_development/05_qmmm"
ENS="$ROOT/19_ensemble"
MAN="$ROOT/12_frame_selection/selection_manifest.tsv"
cd "$ROOT" || exit 1

frames_new() { awk -F'\t' 'NR>1 && $2>=20000 {printf "%05d\n",$2}' "$MAN"; }
running()    { qstat -u "$USER" 2>/dev/null | grep -c "cm19_" ; }

# Is a job for this frame and stage already in the queue? The queue is the only
# reliable test: output files appear when a stage finishes, so a running job and
# an unstarted one look identical on disk. Job names are cm19_<letter><frame>
# and qstat truncates them, so match on the leading characters.
in_queue() {  # in_queue <frame> <letter>
  local want="cm19_$2${1#0}"
  # qstat truncates the name to ten characters including a trailing asterisk,
  # so the field it prints is a PREFIX of the real name, not the whole of it.
  # Test that the expected name starts with the printed field, not the reverse.
  local j
  for j in $(qstat -u "$USER" 2>/dev/null | awk 'NR>5{print $4}' | sed 's/\*$//'); do
    [[ -n "$j" && "$want" == "$j"* ]] && return 0
  done
  return 1
}

have() {  # have <frame> <stage>
  local d="$ENS/frame_$1"
  case "$2" in
    reactant) grep -q "The minimization has converged" "$d/reactant_opt.out" 2>/dev/null ;;
    scan)     [[ -s "$d/scan/win_20.pdb" ]] ;;
    product)  grep -q "The minimization has converged" "$d/product_opt.out" 2>/dev/null ;;
    band)     grep -q "THE NEB OPTIMIZATION HAS CONVERGED" "$d/neb.out" 2>/dev/null ;;
    harvested) [[ -s "$ROOT/19_ensemble_barriers/frame_$1/path_summary.txt" ]] ;;
  esac
}

case "${1:-status}" in

status)
  printf "%-8s %-9s %-6s %-8s %-6s %-10s %s\n" frame reactant scan product band harvested size
  printf -- "------------------------------------------------------------------\n"
  for f in $(frames_new); do
    d="$ENS/frame_$f"
    [[ -d "$d" ]] || { printf "%-8s not prepared\n" "$f"; continue; }
    r=$(have "$f" reactant  && echo yes || echo "-")
    s=$(have "$f" scan      && echo yes || echo "-")
    p=$(have "$f" product   && echo yes || echo "-")
    b=$(have "$f" band      && echo yes || echo "-")
    h=$(have "$f" harvested && echo yes || echo "-")
    sz=$(du -sh "$d" 2>/dev/null | cut -f1)
    printf "%-8s %-9s %-6s %-8s %-6s %-10s %s\n" "$f" "$r" "$s" "$p" "$b" "$h" "$sz"
  done
  echo
  echo "jobs in the queue: $(running)"
  echo "05_qmmm total:     $(du -sh "$ROOT" 2>/dev/null | cut -f1)"
  echo "filesystem free:   $(df -h "$HOME" | awk 'END{print $4}')"
  ;;

stage1)
  n=0
  for f in $(frames_new); do
    have "$f" reactant && continue
    in_queue "$f" r && { echo "  $f already queued"; continue; }
    [[ -s "$ENS/frame_$f/reactant_opt.pbs" ]] || { echo "  $f: no job script"; continue; }
    (cd "$ENS/frame_$f" && qsub reactant_opt.pbs >/dev/null) && { echo "  submitted $f"; n=$((n+1)); }
  done
  echo "stage 1: $n reactant optimisations submitted"
  echo "these are small on disk; the next stage needs each frame's own result"
  ;;

stage2)
  n=0; w=0
  for f in $(frames_new); do
    have "$f" scan && continue
    if ! have "$f" reactant; then w=$((w+1)); continue; fi
    in_queue "$f" s && { echo "  $f already running a scan"; continue; }
    [[ -s "$ENS/frame_$f/scan.pbs" ]] || { echo "  $f: run step19c first"; continue; }
    (cd "$ENS/frame_$f" && qsub scan.pbs >/dev/null) && { echo "  submitted $f"; n=$((n+1)); }
  done
  echo "stage 2: $n scans submitted, $w still waiting on their reactant"
  ;;

stage3)
  n=0; w=0
  for f in $(frames_new); do
    have "$f" product && continue
    if ! have "$f" scan; then w=$((w+1)); continue; fi
    in_queue "$f" p && { echo "  $f already running a product optimisation"; continue; }
    [[ -s "$ENS/frame_$f/product_opt.pbs" ]] || { echo "  $f: no job script"; continue; }
    (cd "$ENS/frame_$f" && qsub product_opt.pbs >/dev/null) && { echo "  submitted $f"; n=$((n+1)); }
  done
  echo "stage 3: $n product optimisations submitted, $w waiting on their scan"
  ;;

stage4)
  size="${2:-4}"
  live=$(running)
  room=$(( size - live ))
  if [[ $room -le 0 ]]; then
    echo "$live already running, batch size $size; nothing submitted"
    echo "wait for these to finish, then: bash $0 harvest"
    exit 0
  fi
  n=0
  for f in $(frames_new); do
    [[ $n -ge $room ]] && break
    have "$f" band && continue
    have "$f" harvested && continue
    have "$f" product || continue
    in_queue "$f" n && { echo "  $f already running a band"; continue; }
    [[ -s "$ENS/frame_$f/tsguess.pdb" ]] || { echo "  $f: no transition-state guess"; continue; }
    (cd "$ENS/frame_$f" && qsub neb.pbs >/dev/null) && { echo "  submitted $f"; n=$((n+1)); }
  done
  echo "stage 4: $n bands submitted, $((live+n)) now running"
  echo "each holds about 5.7 GB while it runs, mostly the approximate Hessian"
  echo "when they finish: bash $0 harvest"
  ;;

harvest)
  echo "harvesting everything that has converged ..."
  OPENBLAS_NUM_THREADS=1 python3 step19e_harvest_barriers.py 2>&1 | tail -6
  echo
  echo "removing the heavy files from harvested frames whose job has ended:"
  freed=0
  for f in $(frames_new); do
    have "$f" harvested || continue
    qstat -u "$USER" 2>/dev/null | grep -q "cm19_n${f:1}" && { echo "  $f still running, left alone"; continue; }
    d="$ENS/frame_$f"
    before=$(du -sm "$d" 2>/dev/null | cut -f1)
    rm -f "$d"/neb.appr.hess "$d"/neb.gu.tmp "$d"/neb_MEP_ALL_trj.xyz \
          "$d"/neb_MEP_trj.xyz "$d"/neb_initial_path_trj.xyz \
          "$d"/neb_MEP_ALL.activeRegion_trj.xyz "$d"/neb_im*.gbw \
          "$d"/neb*.tmp 2>/dev/null
    after=$(du -sm "$d" 2>/dev/null | cut -f1)
    freed=$(( freed + before - after ))
    echo "  $f: ${before} MB -> ${after} MB"
  done
  echo
  echo "freed ${freed} MB"
  echo "05_qmmm total: $(du -sh "$ROOT" 2>/dev/null | cut -f1)"
  echo
  echo "The path summary, barrier and converged geometries are in"
  echo "19_ensemble_barriers and are what the analysis reads."
  ;;

*)
  echo "usage: $0 {status|stage1|stage2|stage3|stage4 [n]|harvest}"
  exit 1
  ;;
esac
