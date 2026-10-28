#!/usr/bin/env bash
# Build s4.exe (and the V2 s3 harness) with the study-3 flags (study3/driver/build.sh), unmodified libraries.
set -e
ROOT=C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad
E=$ROOT/Enterprise/warp_hypad; B=$ROOT/build-oti; D=$ROOT/study4/driver
FL="-Wno-error=literal-suffix -Wno-error=unused-but-set-variable -Wno-error=unused-variable -Wno-error=unused-function -Wno-error=misleading-indentation -O2 -g -DNDEBUG -std=c++17 -Wall -Wextra -Wpedantic -Werror -Wshadow -Wconversion -ffp-contract=off -fopenmp -DWARP_HYPAD_HAVE_OTILIB=1"
INC="-I$E/tests -I$E/include -I$E/external/otilib/include -I$E/external/otilib/src/c"
C:/Users/vfn333/AppData/Local/anaconda3/python.exe -I $D/assemble_s4.py $ROOT
g++ $FL $INC -c $D/s4.cpp -o $D/s4.o
g++ $FL $INC -c $E/tests/o0d_convergence_checks.cpp -o $D/o0d.o
g++ -fopenmp -O2 $D/s4.o $D/o0d.o $B/libwarp_hypad_core.a $B/libwarp_hypad_otilib.a $B/libwarp_hypad_core.a -o $D/s4.exe
g++ $FL $INC -c $D/v2_s3harness.cpp -o $D/v2h.o
g++ -fopenmp -O2 $D/v2h.o $D/o0d.o $B/libwarp_hypad_core.a $B/libwarp_hypad_otilib.a $B/libwarp_hypad_core.a -o $D/v2_s3harness.exe
g++ --version | head -1
echo built
