#!/usr/bin/env bash
# Assessment priorities 2-5 (PROTOCOL 7): iterates listed in results/assess_list.txt as "run iterate reason".
ROOT=C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad/study3
export OMP_NUM_THREADS=${1:-12}
A=$ROOT/assess; S=$A/snap; mkdir -p $S
while read -r r k why; do
  [ -z "$r" ] && continue
  tag=${r}@$k
  grep -q "^$tag," $A/iters.csv 2>/dev/null && continue
  sed -n "$((k+1))p" $ROOT/runs/$r/ys.txt | tr ' ' '\n' | sed '/^$/d' > $S/y_${r}_it$k.txt
  $ROOT/driver/s3.exe assess --y $S/y_${r}_it$k.txt --seed 1 --tag $tag --out $A/iters.csv
done < $ROOT/results/assess_list.txt 2>&1 | tee -a $ROOT/logs/assess_iters.log
