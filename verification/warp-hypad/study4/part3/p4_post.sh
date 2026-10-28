#!/usr/bin/env bash
# Part 3 (P4) post-run work, writing only to new paths so Part 2 artifacts stay byte-identical.
# Same code and arguments as driver/assess.sh and the Part 2 post-run steps; outputs under part3/ and new hessians/ dirs.
ROOT=C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad; S4=$ROOT/study4; X=$S4/driver/s4.exe; P=$S4/driver/hessian_passes.txt
M=$S4/part3; PY="C:/Users/vfn333/AppData/Local/anaconda3/python.exe -I"
L1=3.75e-8; P4="A4L1f B4L1f A6L1f B6L1f"
export OMP_NUM_THREADS=12
dens() { case $1 in A4*|B4*) echo 4;; *) echo 6;; esac; }
{
echo "$(date -Iseconds) part3 mirror build"
mkdir -p $M/runs $M/assess/snap $M/results $M/hessians
for d in $S4/runs/*/; do r=$(basename $d); mkdir -p $M/runs/$r
  for f in trace.csv termination.txt final_diag.csv ys.txt y_final.txt outer.csv; do cp $d/$f $M/runs/$r/; done; done
for f in assess.csv starts.csv iters.csv; do cp $S4/assess/$f $M/assess/$f; done
for d in $S4/hessians/final_*/; do mkdir -p $M/hessians/$(basename $d); cp $d/summary.csv $M/hessians/$(basename $d)/; done
echo "$(date -Iseconds) P4 converged-final Hessians"
for r in $P4; do
  st=$(cut -d' ' -f1 $S4/runs/$r/termination.txt)
  case $st in CONVERGED|CONVERGED_LS_LIMITED) ;; *) continue;; esac
  extra=""; case $r in B*) extra="--set $S4/runs/$r/set_final.txt";; esac
  [ -f $S4/hessians/final_$r/meta.txt ] || $X hessian --y $S4/runs/$r/y_final.txt --ws $(dens $r):1 $extra --lambda $L1 --passes $P --out $S4/hessians/final_$r --tag final_$r
  $PY $S4/driver/hess_analyze.py $S4 design $S4/hessians/final_$r > $S4/hessians/final_$r/analysis.json
  mkdir -p $M/hessians/final_$r; cp $S4/hessians/final_$r/summary.csv $M/hessians/final_$r/
  echo "$(date -Iseconds) final_$r done"
done
echo "$(date -Iseconds) P4 finals assessment (seeds 1, 2, --full 1)"
for r in $P4; do for s in 1 2; do $X assess --y $S4/runs/$r/y_final.txt --seed $s --tag $r --full 1 --out $M/assess/assess.csv; done; done
echo "$(date -Iseconds) P4 trajectory assessment"
$PY $S4/driver/analyze_runs.py $M select
tr -d '\r' < $M/results/assess_list.txt | grep -E "^(A4L1f|B4L1f|A6L1f|B6L1f) " > $M/results/assess_list_p4.txt
while read -r r k; do
  [ -z "$r" ] && continue
  sed -n "$((k+1))p" $S4/runs/$r/ys.txt | tr ' ' '\n' | sed '/^$/d' > $M/assess/snap/y_${r}_it$k.txt
  $X assess --y $M/assess/snap/y_${r}_it$k.txt --seed 1 --tag ${r}@$k --out $M/assess/iters.csv
done < $M/results/assess_list_p4.txt
echo "$(date -Iseconds) X4L1fromL2 iterate 1 (Part 2 deviation 4; Vera recommendation)"
sed -n "2p" $S4/runs/X4L1fromL2/ys.txt | tr ' ' '\n' | sed '/^$/d' > $M/assess/snap/y_X4L1fromL2_it1.txt
$X assess --y $M/assess/snap/y_X4L1fromL2_it1.txt --seed 1 --tag X4L1fromL2@1 --out $M/assess/iters.csv
echo "$(date -Iseconds) PART3 POST DONE"
} 2>&1 | tee -a $S4/logs/post_p4.log
