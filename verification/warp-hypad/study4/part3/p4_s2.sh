#!/usr/bin/env bash
# Part 3: S2 (v4 s7, descriptive) at the strict P4 finals; same steps as Part 2 (B on its final set).
ROOT=C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad; S4=$ROOT/study4; X=$S4/driver/s4.exe; P=$S4/driver/hessian_passes.txt
PY="C:/Users/vfn333/AppData/Local/anaconda3/python.exe -I"; export OMP_NUM_THREADS=12; L1=3.75e-8
{
for spec in A4L1f:4 B4L1f:4 A6L1f:6 B6L1f:6; do
  r=${spec%%:*}; d=${spec#*:}; extra=""; case $r in B*) extra="--set $S4/runs/$r/set_final.txt";; esac
  $PY $S4/driver/s2_newton.py $S4/hessians/final_$r $S4/hessians/final_$r/y_newton.txt
  [ -f $S4/hessians/s2_$r/meta.txt ] || $X hessian --y $S4/hessians/final_$r/y_newton.txt --ws $d:1 $extra --lambda $L1 --passes $P --out $S4/hessians/s2_$r --tag s2_$r
  $PY $S4/driver/hess_analyze.py $S4 design $S4/hessians/s2_$r > $S4/hessians/s2_$r/analysis.json
  echo "$(date -Iseconds) s2_$r done"
done
echo "$(date -Iseconds) P4 S2 DONE"
} 2>&1 | tee -a $S4/logs/post_p4.log
