#!/bin/bash
# Collect ALL barriers once jobs finish. barrier = (E_TS - E_react)*627.509. Flags fragmented reactants.
echo "shell | E_react(Eh)   E_TS(Eh)     | barrier(kcal/mol) | vs bare +17.47 | note"
for tag in 2p0 3p0 4p0 5p0 6p0 7p0 7p1 7p2 7p3 7p4 7p5 7p6 7p7 7p8 7p9 8p0 8p1 8p2 8p3 8p4 8p5 8p6 8p7 8p8 8p9 9p0 10p0 11p0 12p0 13p0 14p0 15p0; do
  rlog=run_shell$tag/ladder_${tag}_relax.log
  tsE=run_barr_$tag/shell${tag}ts_energy.txt
  [ -f "$rlog" ] || continue
  er=$(grep "Final optimized energy" "$rlog" | tail -1 | grep -oE '\-?[0-9]+\.[0-9]+')
  # fragmented reactant? check O3-C4 of the relaxed reactant
  rx=run_shell$tag/shell${tag}_relaxed.xyz
  frag=""
  if [ -f "$rx" ]; then
    frag=$(python3 -c "
import math
L=open('$rx').read().splitlines();a=[[float(v) for v in L[2+i].split()[1:4]] for i in range(24)]
d=math.sqrt(sum((a[7][k]-a[8][k])**2 for k in range(3)))
print('FRAG' if d>2.2 else '')
")
  fi
  if [ -f "$tsE" ] && [ -n "$er" ]; then
    ets=$(cat "$tsE")
    python3 -c "
er=float('$er'); ets=float('$ets'); b=(ets-er)*627.509
note='$frag' or ('catalytic' if b<17.47 else 'anti/non-cat')
print('%-5s | %12.6f %12.6f | %+8.2f        | %+6.2f        | %s'%('$tag',er,ets,b,b-17.47,note))
"
  else
    echo "$(printf '%-5s' $tag) | (waiting: E_react=$([ -n "$er" ]&&echo ok||echo -) E_TS=$([ -f "$tsE" ]&&echo ok||echo -)) $frag"
  fi
done
