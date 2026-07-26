#pragma once

#include <string>

#include "warp_hypad/core/verification.hpp"

namespace program_b {
namespace checks {

// Quantities the B3 gate is required to report.
struct B3Metrics {
    std::string fourier_convention;
    std::string model_family;
    std::string pole_locations;

    double min_loss_eigenvalue = 0.0;
    double max_passivity_violation = 0.0;
    double worst_passivity_frequency = 0.0;

    double centre_mismatch = 0.0;
    double mean_mismatch = 0.0;
    double max_mismatch = 0.0;
    double worst_mismatch_frequency = 0.0;

    double hard_bound_violation = 0.0;
    double soft_bound_violation = 0.0;

    double worst_gradient_analytic = 0.0;
    double worst_gradient_fd = 0.0;
    double worst_gradient_complex_step = 0.0;
    double worst_gradient_oti = 0.0;
    int oti_batches = 0;

    std::string termination_reason;
    double worst_parameter_recovery_error = 0.0;
    double one_pole_recovery_error = 0.0;
    double two_pole_recovery_error = 0.0;
    double restricted_residual = 0.0;
    double worst_coordinate_change = 0.0;
    double real_limit_residual = 0.0;
    bool oti_available = false;
};

void run_b3_checks(warp_hypad::VerificationReport& report,
                   const std::string& out_dir, B3Metrics& metrics);

}  // namespace checks
}  // namespace program_b
