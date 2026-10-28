#!/usr/bin/env bash
# Assessment priority 1 (PROTOCOL 7): baseline (2 seeds, full), floor, every final (2 seeds, full), axi, fields.
ROOT=C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad/study3
export OMP_NUM_THREADS=${1:-12}
A=$ROOT/assess; mkdir -p $A $ROOT/fields
X=$ROOT/driver/s3.exe
{
$X floor
for s in 1 2; do $X assess --y BASELINE --seed $s --tag baseline --full 1 --out $A/assess.csv; done
$X axi --y BASELINE
$X field --y BASELINE --out $ROOT/fields/baseline.txt
for d in $ROOT/runs/*_L*/; do
  r=$(basename $d)
  [ -f $d/y_final.txt ] || { echo "$r: no y_final"; continue; }
  for s in 1 2; do $X assess --y $d/y_final.txt --seed $s --tag $r --full 1 --out $A/assess.csv; done
  $X axi --y $d/y_final.txt
  $X field --y $d/y_final.txt --out $ROOT/fields/$r.txt
done
} 2>&1 | tee -a $ROOT/logs/assess_finals.log
