#!/bin/bash
# Read NEB barriers for PASS designs from ASH_NEB.result / the final Image energy table.
echo "design | status | barrier dE(kcal/mol) | vs bare +17.47"
for tag in 2p0 7p1 7p3 7p4 7p6 7p8 8p1 8p2 8p5 8p6 8p9 12p0 13p0 15p0; do
  d=run_nebP_$tag; log=$d/nebP_$tag.log; res=$d/ASH_NEB.result
  [ -f "$log" ] || { printf '%-6s| no log\n' "$tag"; continue; }
  if [ -f "$res" ] && grep -q '(fail)' "$res"; then
    printf '%-6s| NEB failed to converge\n' "$tag"
  elif grep -qE "^>> NEB DONE" "$log"; then
    # pull the HEI dE (kcal/mol) from the final Image energy table
    b=$(grep -A20 "Image  Energy(Eh)" "$log" | grep -E "^ *[0-9]+ " | awk '{print $3}' | sort -g | tail -1)
    [ -n "$b" ] && python3 -c "b=float('$b'); print('%-6s| DONE   | %+8.2f            | %+.2f'%('$tag',b,b-17.47))" || printf '%-6s| DONE but no barrier parsed\n' "$tag"
  else
    n=$(grep -c 'ORCA energy:' "$log")
    printf '%-6s| running (%d img E)\n' "$tag" "$n"
  fi
done
