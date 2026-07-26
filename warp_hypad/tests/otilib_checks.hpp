#pragma once

#include <string>

#include "warp_hypad/core/verification.hpp"

namespace warp_hypad {
namespace checks {

// Analytic and identity-based verification of the OTIlib adapter.
// Writes analytic_derivatives.csv into out_dir.
void run_otilib_checks(VerificationReport& report, const std::string& out_dir);

// OTI against centered finite differences across a range of step sizes.
// Writes finite_difference_comparison.csv into out_dir.
void run_otilib_fd_checks(VerificationReport& report, const std::string& out_dir);

}  // namespace checks
}  // namespace warp_hypad
