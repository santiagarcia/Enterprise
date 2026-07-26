#pragma once

#include <string>

#include "warp_hypad/core/verification.hpp"

namespace warp_hypad {
namespace checks {

// A3 structural checks. Writes shape_profiles.csv, metric_slice_xy.csv and
// centerline_x.csv into out_dir.
void run_metric_checks(VerificationReport& report, const std::string& out_dir);

// A3 HyPAD parameter sensitivities. Writes hypad_vs_analytic.csv,
// hypad_vs_fd.csv and sensitivity_slice_xy.csv. Records a single skipped case
// when the OTIlib backend is not compiled in.
void run_metric_sensitivity_checks(VerificationReport& report,
                                   const std::string& out_dir);

}  // namespace checks
}  // namespace warp_hypad
