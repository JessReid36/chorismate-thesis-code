#!/usr/bin/env bash
# State of every new-ensemble frame: band, barrier, in vacuo, size, exit.
# Reads only. Falls back to neb.out when PBS did not stage neb.pbs.out back
# (seen on frame 41786).
ENS="$HOME/system_development/05_qmmm/19_ensemble"
BAR="$HOME/system_development/05_qmmm/19_ensemble_barriers"
Q=$(qstat -u "$USER" 2>/dev/null)
printf '%-8s %-9s %-9s %-9s %-8s %-7s %s\n' frame state barrier vac size runtime exit
printf -- '---------------------------------------------------------------------------------\n'
run=0; done_=0; todo=0
for d in "$ENS"/frame_*; do
  f=$(basename "$d"); f=${f#frame_}
  [ "$f" -ge 20000 ] 2>/dev/null || continue
  jid=$(printf '%s\n' "$Q" | awk -v p="cm19_n${f}" '
    $1 ~ /^[0-9]+\./ { n=$4; sub(/\*$/,"",n); if (n ~ /^cm19_n/ && index(p,n)==1) { print $1; exit } }')
  bar="-"; vac="-"; rt="-"; ex="-"; st="not started"
  if [ -f "$d/neb.out" ]; then
    if grep -q 'THE NEB OPTIMIZATION HAS CONVERGED' "$d/neb.out"; then
      bar=$(awk '/THE NEB OPTIMIZATION HAS CONVERGED/{c=1} c&&/<= CI/{print $4; exit}' "$d/neb.out")
      st="converged"
    else
      st="band running"
    fi
    rt=$(grep -m1 'TOTAL RUN TIME' "$d/neb.out" | sed 's/.*TIME: *//; s/ days /d/; s/ hours /h/; s/ minutes.*//')
    if grep -q 'ORCA TERMINATED NORMALLY' "$d/neb.out"; then ex="normal"; fi
    if [ -f "$d/neb.pbs.out" ] && grep -q 'walltime' "$d/neb.pbs.out"; then ex="WALLTIME"; fi
  fi
  [ -n "$jid" ] && st="RUNNING ${jid%%.*}"
  v=$(awk -F'\t' -v fr="$f" '!/^#/ && $1==fr {print $7}' "$BAR/ensemble_barriers.tsv" 2>/dev/null)
  [ -n "$v" ] && vac=$(printf '%.2f' "$v")
  sz=$(du -sh "$d" 2>/dev/null | cut -f1)
  case "$st" in RUNNING*) run=$((run+1));; converged) done_=$((done_+1));; *) todo=$((todo+1));; esac
  printf '%-8s %-9s %-9s %-9s %-8s %-7s %s\n' "$f" "${st:0:9}" "$bar" "$vac" "$sz" "$rt" "$ex"
done
echo
echo "running $run   converged $done_   not started $todo"
echo
echo "converged but NO in vacuo yet (next s8_invacuo_new.pbs batch):"
awk -F'\t' '!/^#/ && $1+0>=20000 && $6!="" && $6+0!=0 && $7=="" {printf " %s", $1}' "$BAR/ensemble_barriers.tsv" 2>/dev/null
echo
