#!/usr/bin/env bash
# One lane (PROTOCOL 5): form F, density D, threads T; runs lambda1 then lambda2 from S1.
# args: F D T GLOBAL_DEADLINE_EPOCH L1 L2
ROOT=C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad/study3
F=$1; D=$2; T=$3; G=$4; L1=$5; L2=$6
k=1
for L in $L1 $L2; do
  now=$(date +%s)
  if (( now >= G )); then echo "${F}_${D}_L$k NOT_STARTED (global deadline)" >> $ROOT/runs/runs.txt; k=$((k+1)); continue; fi
  if (( k == 1 )); then dl=$(( now + 8100 < G ? now + 8100 : G )); else dl=$G; fi
  d=$ROOT/runs/${F}_${D}_L$k; mkdir -p $d
  s=$(date +%s)
  OMP_NUM_THREADS=$T $ROOT/driver/s3.exe run --form $F --ws $D:1 --start baseline_zero --lambda $L --kmax 600 \
      --maxfev 3000 --deadline $dl --diag-every 1 --out $d > $d/log.txt 2>&1
  echo "${F}_${D}_L$k rc=$? wall=$(( $(date +%s) - s ))s start=$(date -d @$s +%T) deadline=$(date -d @$dl +%T) threads=$T lambda=$L $(cat $d/termination.txt)" >> $ROOT/runs/runs.txt
  k=$((k+1))
done
