#pragma once

#include <string>

#include "warp_hypad/core/verification.hpp"

namespace program_b {
namespace checks {

// Program B, Milestone B1: metric <-> constitutive-tensor mapping.
// Writes warp_medium_profile.csv and realizability_map.csv into out_dir.
void run_b1_checks(warp_hypad::VerificationReport& report,
                   const std::string& out_dir);

// HyPAD sensitivities through the mapping. Writes medium_sensitivity.csv.
// Records a single skipped case when the OTIlib backend is absent.
void run_b1_sensitivity_checks(warp_hypad::VerificationReport& report,
                               const std::string& out_dir);

}  // namespace checks
}  // namespace program_b
