#pragma once

#include <string>

#include "warp_hypad/core/verification.hpp"

namespace program_b {
namespace checks {

struct B4CMetrics {
    std::string old_metricity_claim;
    std::string corrected_metricity_claim;

    double fresnel_equality_under_axion = 0.0;
    double constitutive_difference_under_axion = 0.0;
    double zeta_block_difference_under_axion = 0.0;
    double axion_coefficient_recovered = 0.0;
    double conformal_fresnel_zero_set_shift = 0.0;

    double worst_repeated_root_residual = 0.0;
    double metric_double_root_separation = 0.0;
    double anisotropic_split = 0.0;
    double coefficient_reconstruction_error = 0.0;

    double worst_stack_gradient_analytic = 0.0;
    double worst_stack_gradient_fd = 0.0;
    double magnetoelectric_sensitivity = 0.0;
    double magnetoelectric_sensitivity_restricted = 0.0;

    double energy_balance_regression = 0.0;
    std::string exclusions = "none";
    bool oti_available = false;
};

void run_b4c_checks(warp_hypad::VerificationReport& report,
                    const std::string& out_dir, B4CMetrics& metrics);

}  // namespace checks
}  // namespace program_b
