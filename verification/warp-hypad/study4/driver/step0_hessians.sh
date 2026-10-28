#!/usr/bin/env bash
# Step 0 (PROTOCOL-v3 2.3): OTI Hessians at the 10 designs, 12 threads, one design at a time.
ROOT=C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad
S3=$ROOT/study3/runs; S4=$ROOT/study4; X=$S4/driver/s4.exe; P=$S4/driver/hessian_passes.txt
export OMP_NUM_THREADS=12
L1=3.75e-8; L2=3.75e-9
h() { local tag=$1; shift; [ -f $S4/hessians/$tag/meta.txt ] && { echo "$tag exists"; return; }
      $X hessian --passes $P --out $S4/hessians/$tag --tag $tag "$@"; }
{
h A_4_L2_it84 --y $S4/starts/A_4_L2_it84.txt --ws 4:1 --lambda $L2
h B_4_L2 --y $S3/B_4_L2/y_final.txt --ws 4:1 --recon $S3/B_4_L2 --lambda $L2
h A_6_L2_it90 --y $S4/starts/A_6_L2_it90.txt --ws 6:1 --lambda $L2
h B_6_L2 --y $S3/B_6_L2/y_final.txt --ws 6:1 --recon $S3/B_6_L2 --lambda $L2
h A_4_L2 --y $S3/A_4_L2/y_final.txt --ws 4:1 --lambda $L2
h A_6_L2 --y $S3/A_6_L2/y_final.txt --ws 6:1 --lambda $L2
h A_4_L1 --y $S3/A_4_L1/y_final.txt --ws 4:1 --lambda $L1
h A_6_L1 --y $S3/A_6_L1/y_final.txt --ws 6:1 --lambda $L1
h B_4_L1 --y $S3/B_4_L1/y_final.txt --ws 4:1 --recon $S3/B_4_L1 --lambda $L1
h B_6_L1 --y $S3/B_6_L1/y_final.txt --ws 6:1 --recon $S3/B_6_L1 --lambda $L1
echo STEP0_DONE $(date)
} 2>&1 | tee -a $S4/logs/step0_hessians.log
