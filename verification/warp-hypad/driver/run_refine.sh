#!/usr/bin/env bash
# Re-run the O2 zero start on coarser/finer WORKING (optimizer) sets, and a
# longer run at the original set. Validation always at the original 8:1 set;
# finer validation is done afterwards by `wh_probe eval`.
set -u
ROOT=C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad
EXE=$ROOT/driver/wh_probe2.exe
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-12}
run() { # name ws maxacc
  d=$ROOT/results/refine/$1; mkdir -p $d; cd $d
  s=$(date +%s)
  $EXE o2 --ws $2 --vs 8:1:3.5:12.25 --maxacc $3 --only baseline_zero --out . > o2.log 2>&1
  echo "$1 ws=$2 maxacc=$3 rc=$? wall=$(( $(date +%s) - s ))s" | tee -a $ROOT/results/refine/runs.txt
}
run w3s1_it12 3:1:1.7:2.6 12
run w4s2_it12 4:2:1.7:2.6 12
run w8s1_it12 8:1:1.7:2.6 12
run w4s1_it60 4:1:1.7:2.6 60
run w2s1_it12 2:1:1.7:2.6 12
run w8s2_it12 8:2:1.7:2.6 12
