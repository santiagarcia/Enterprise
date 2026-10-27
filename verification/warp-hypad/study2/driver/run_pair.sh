#!/usr/bin/env bash
# Run the matched A/B pair for one (density, start). args: density start_label start_name threads
ROOT=C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad/study2
for F in A B; do
  d=$ROOT/runs/${F}_$1_$2; mkdir -p $d
  ( s=$(date +%s); OMP_NUM_THREADS=$4 $ROOT/driver/s2.exe run --form $F --ws $1:1 --start $3 --kmax 150 --maxfev 1500 --out $d > $d/log.txt 2>&1
    echo "${F}_$1_$2 rc=$? wall=$(( $(date +%s) - s ))s start=$(date -d @$s +%T)" >> $ROOT/runs/runs.txt ) &
done
wait
