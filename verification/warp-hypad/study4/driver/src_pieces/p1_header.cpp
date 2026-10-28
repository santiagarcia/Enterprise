// s4: study-4 driver (study4/PROTOCOL-v3.md). Builds against the UNMODIFIED Enterprise/warp_hypad headers +
// libraries at 477141d with the flags of study3/driver/build.sh.
//
// Assembled from study3/driver/s3.cpp: every block not listed below is copied verbatim (node sets, rho_pt,
// node_rho, residual, compass polish, in-loop search, 8:1 diagnostics, traj, Nelder-Mead, box_search, assess,
// floor, axi, field, timing, check). Changed/new blocks (PROTOCOL-v3 sections 1-4, 9):
//   * gradient(): forward sum unchanged (bit-identical, gate V2); per-node contributions are also summed in
//     reverse order, with Neumaier compensation and in long double (H1b/H1c, v3 N3).
//   * refine_B(): split decisions verbatim; templated on the rho evaluator; returns an outcome record; no
//     persistent `capped` flag (v2 R1).
//   * run: dense BFGS initialised from the measured Hessian, refresh rule, B nested loop (sections 3, 4).
//   * hessian / fd / recon / v2dump / selftest modes (sections 2, 9).
//
// Per-node values are computed in parallel and summed serially in node order (deterministic).

#include <omp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <ctime>
#include <deque>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

