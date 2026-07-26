#pragma once

#include <string>

#include "warp_hypad/core/verification.hpp"

namespace program_b {
namespace checks {

// Quantities the B2 gate is required to report, filled in as the cases run.
struct B2Metrics {
    double worst_parameter_roundtrip = 0.0;
    double worst_symmetry_residual = 0.0;
    double worst_spatiality_residual = 0.0;
    double worst_objective_at_target = 0.0;
    double worst_gradient_analytic = 0.0;
    double worst_gradient_fd = 0.0;
    double worst_gradient_oti = 0.0;
    double worst_batch_consistency = 0.0;
    int oti_batches = 0;
    int oti_evaluations = 0;

    double isotropic_residual = 0.0;
    double anisotropic_residual = 0.0;
    double magnetoelectric_residual = 0.0;
    double metric_induced_residual = 0.0;
    double restricted_family_residual = 0.0;
    double restricted_vs_least_squares = 0.0;
    double bounded_residual = 0.0;
    double worst_coordinate_change = 0.0;
    bool oti_available = false;
};

void run_b2_checks(warp_hypad::VerificationReport& report,
                   const std::string& out_dir, B2Metrics& metrics);

}  // namespace checks
}  // namespace program_b
