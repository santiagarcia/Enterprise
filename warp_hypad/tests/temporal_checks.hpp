#pragma once

#include <string>

#include "warp_hypad/core/verification.hpp"

namespace warp_hypad {
namespace checks {

// Milestone A2T: temporal derivative operators, and the rigid-translation
// identity the time-dependence decision rests on.
// Writes temporal_convergence.csv and identity_profile.csv into out_dir.
void run_temporal_checks(VerificationReport& report, const std::string& out_dir);

}  // namespace checks
}  // namespace warp_hypad
