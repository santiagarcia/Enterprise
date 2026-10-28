#!/usr/bin/env bash
# Study-4 assessment (PROTOCOL-v3 section 7). Code: s4 assess = s3 assess verbatim.
#   assess.sh finals RUN...      : run finals, --full 1, seeds 1 and 2 -> assess/assess.csv
#   assess.sh starts             : start designs, seeds 1 and 2         -> assess/starts.csv
#   assess.sh iters LISTFILE     : "RUN IT" lines, seed 1               -> assess/iters.csv
ROOT=C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad
S3=$ROOT/study3/runs; S4=$ROOT/study4; X=$S4/driver/s4.exe; A=$S4/assess; mkdir -p $A/snap
export OMP_NUM_THREADS=${THREADS:-12}
mode=$1; shift
case $mode in
finals)
  for r in "$@"; do
    grep -q "^$r," $A/assess.csv 2>/dev/null && continue
    for s in 1 2; do $X assess --y $S4/runs/$r/y_final.txt --seed $s --tag $r --full 1 --out $A/assess.csv; done
  done ;;
starts)
  for t in A_4_L2_it84:$S4/starts/A_4_L2_it84.txt A_6_L2_it90:$S4/starts/A_6_L2_it90.txt \
           A_4_L2:$S3/A_4_L2/y_final.txt A_6_L2:$S3/A_6_L2/y_final.txt B_4_L2:$S3/B_4_L2/y_final.txt \
           B_6_L2:$S3/B_6_L2/y_final.txt A_4_L1:$S3/A_4_L1/y_final.txt A_6_L1:$S3/A_6_L1/y_final.txt \
           B_4_L1:$S3/B_4_L1/y_final.txt B_6_L1:$S3/B_6_L1/y_final.txt; do
    tag=${t%%:*}; f=${t#*:}
    grep -q "^start_$tag," $A/starts.csv 2>/dev/null && continue
    for s in 1 2; do $X assess --y $f --seed $s --tag start_$tag --full 1 --out $A/starts.csv; done
  done ;;
iters)
  while read -r r k; do
    [ -z "$r" ] && continue
    tag=${r}@$k
    grep -q "^$tag," $A/iters.csv 2>/dev/null && continue
    sed -n "$((k+1))p" $S4/runs/$r/ys.txt | tr ' ' '\n' | sed '/^$/d' > $A/snap/y_${r}_it$k.txt
    $X assess --y $A/snap/y_${r}_it$k.txt --seed 1 --tag $tag --out $A/iters.csv
  done < "$1" ;;
esac 2>&1 | tee -a $S4/logs/assess.log
