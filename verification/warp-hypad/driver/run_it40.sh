#!/usr/bin/env bash
# Zero-start runs with per-iterate validation tracking (wh_probe4). args: name ws maxacc threads
ROOT=C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad
d=$ROOT/results/refine/$1; mkdir -p $d; cd $d; s=$(date +%s)
OMP_NUM_THREADS=$4 $ROOT/driver/wh_probe4.exe o2 --ws $2 --vs 8:1:3.5:12.25 --maxacc $3 --only baseline_zero --out . > o2.log 2>&1
echo "$1 ws=$2 maxacc=$3 rc=$? wall=$(( $(date +%s) - s ))s" | tee -a $ROOT/results/refine/runs.txt
