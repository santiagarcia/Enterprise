#!/usr/bin/env bash
ROOT=C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad; S3=$ROOT/study3/runs; S4=$ROOT/study4; F=$S4/hessians/fd
export OMP_NUM_THREADS=12
{
$S4/driver/s4.exe fd --y $S3/A_4_L2/y_final.txt --ws 4:1 --dirs $F/A_4_L2.dirs --out $F/A_4_L2.out
$S4/driver/s4.exe fd --y $S3/A_4_L1/y_final.txt --ws 4:1 --dirs $F/A_4_L1.dirs --out $F/A_4_L1.out
$S4/driver/s4.exe fd --y $S3/B_4_L2/y_final.txt --ws 4:1 --recon $S3/B_4_L2 --dirs $F/B_4_L2.dirs --out $F/B_4_L2.out
$S4/driver/s4.exe fd --y $S3/A_6_L2/y_final.txt --ws 6:1 --dirs $F/A_6_L2.dirs --out $F/A_6_L2.out
$S4/driver/s4.exe fd --y $S3/B_6_L2/y_final.txt --ws 6:1 --recon $S3/B_6_L2 --dirs $F/B_6_L2.dirs --out $F/B_6_L2.out
echo FD_DONE $(date)
} 2>&1 | tee -a $S4/logs/fd.log
