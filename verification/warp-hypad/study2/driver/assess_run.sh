#!/usr/bin/env bash
# Independent assessment of one finished run: final design (seeds 1 and 2), snapshots every 25
# accepted iterates (seed 1), meridian field of the final design, axisymmetry check. args: run threads
ROOT=C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad/study2
r=$1; T=${2:-24}; d=$ROOT/runs/$r; S=$ROOT/assess/snap/$r; mkdir -p $S $ROOT/fields
export OMP_NUM_THREADS=$T
$ROOT/driver/s2.exe assess --y $d/y_final.txt --seed 1 --tag $r --out $ROOT/assess/assess.csv
$ROOT/driver/s2.exe assess --y $d/y_final.txt --seed 2 --tag $r --out $ROOT/assess/assess.csv
$ROOT/driver/s2.exe axi --y $d/y_final.txt
$ROOT/driver/s2.exe field --y $d/y_final.txt --out $ROOT/fields/$r.txt
n=$(wc -l < $d/ys.txt)
for ((k=25; k<n; k+=25)); do
  sed -n "$((k+1))p" $d/ys.txt | tr ' ' '\n' | sed '/^$/d' > $S/y_it$k.txt
  $ROOT/driver/s2.exe assess --y $S/y_it$k.txt --seed 1 --tag ${r}@$k --out $ROOT/assess/snapshots.csv
done
