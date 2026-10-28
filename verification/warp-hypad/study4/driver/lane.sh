#!/usr/bin/env bash
# One study-4 run (PROTOCOL-v3 sections 5, 6). args: ID FORM DENS THREADS LAMBDA K Y0 HDIR WALLCAP_S [RECON_DIR]
ROOT=C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad; S4=$ROOT/study4
ID=$1; F=$2; D=$3; T=$4; L=$5; K=$6; Y0=$7; H=$8; CAP=$9; RC=${10}
d=$S4/runs/$ID; mkdir -p $d; s=$(date +%s); dl=$((s + CAP))
extra=""; [ -n "$RC" ] && extra="--recon $RC"
OMP_NUM_THREADS=$T $S4/driver/s4.exe run --form $F --ws $D:1 --y $Y0 --lambda $L --kmax $K --hessian $H --passes $S4/driver/hessian_passes.txt --deadline $dl --diag-every 1 --out $d $extra > $d/log.txt 2>&1
echo "$ID rc=$? wall=$(( $(date +%s) - s ))s start=$(date -d @$s +%T) cap=${CAP}s threads=$T lambda=$L K=$K $(cat $d/termination.txt)" >> $S4/runs/runs.txt
