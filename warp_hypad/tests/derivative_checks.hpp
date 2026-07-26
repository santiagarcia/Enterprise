#pragma once

#include <string>

#include "warp_hypad/core/verification.hpp"

namespace warp_hypad {
namespace checks {

// Milestone A2. Writes convergence.csv, profile_x.csv and error_map_xy.csv
// into out_dir.
void run_derivative_checks(VerificationReport& report, const std::string& out_dir);

}  // namespace checks
}  // namespace warp_hypad
