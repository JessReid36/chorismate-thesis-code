#!/bin/bash
cd ~/system_development/05_qmmm
todo=$(for f in $(awk -F'\t' 'NR>1 && $2>=20000 {printf "%05d\n",$2}' 12_frame_selection/selection_manifest.tsv); do
  grep -q "The minimization has converged" 19_ensemble/frame_$f/reactant_opt.out 2>/dev/null || continue
  [[ -s 19_ensemble/frame_$f/scan/targets.txt ]] && continue
  echo -n "$f "
done)
[ -n "$todo" ] && { echo "preparing: $todo"; bash step19c_scan_product.sh $todo >/dev/null 2>&1; }
bash run_ensemble_batched.sh stage2 | tail -1
bash run_ensemble_batched.sh stage3 | tail -1
qstat -u 18660916 | awk 'NR>5{print substr($4,1,7), $10}' | sort | uniq -c
