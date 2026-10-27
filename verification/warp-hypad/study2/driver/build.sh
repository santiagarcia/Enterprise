#!/usr/bin/env bash
# Build s2.exe against the unmodified repo libraries (study-1 flags + -fopenmp).
set -e
ROOT=C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad
E=$ROOT/Enterprise/warp_hypad; B=$ROOT/build-oti
FL="-Wno-error=literal-suffix -Wno-error=unused-but-set-variable -Wno-error=unused-variable -Wno-error=unused-function -Wno-error=misleading-indentation -O2 -g -DNDEBUG -std=c++17 -Wall -Wextra -Wpedantic -Werror -Wshadow -Wconversion -ffp-contract=off -fopenmp -DWARP_HYPAD_HAVE_OTILIB=1"
INC="-I$E/tests -I$E/include -I$E/external/otilib/include -I$E/external/otilib/src/c"
g++ $FL $INC -c $ROOT/study2/driver/s2.cpp -o $ROOT/study2/driver/s2.o
g++ $FL $INC -c $E/tests/o0d_convergence_checks.cpp -o $ROOT/study2/driver/o0d.o
g++ -fopenmp -O2 $ROOT/study2/driver/s2.o $ROOT/study2/driver/o0d.o $B/libwarp_hypad_core.a $B/libwarp_hypad_otilib.a $B/libwarp_hypad_core.a -o $ROOT/study2/driver/s2.exe
echo built
