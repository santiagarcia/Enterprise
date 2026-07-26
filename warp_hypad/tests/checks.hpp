#pragma once

#include <string>

#include "warp_hypad/core/verification.hpp"

// Milestone 1 verification groups. Each function appends its cases to the
// shared report; none of them terminate the process on failure, because the
// gate is meant to show the full picture of what is broken, not just the first
// thing that broke.

namespace warp_hypad {
namespace checks {

void run_grid_checks(VerificationReport& report);
void run_layout_checks(VerificationReport& report);

// Also writes inverse_residuals.csv into out_dir for the Python plots.
void run_algebra_checks(VerificationReport& report, const std::string& out_dir);

void run_contraction_checks(VerificationReport& report);
void run_genericity_checks(VerificationReport& report);
void run_config_checks(VerificationReport& report);

}  // namespace checks
}  // namespace warp_hypad
