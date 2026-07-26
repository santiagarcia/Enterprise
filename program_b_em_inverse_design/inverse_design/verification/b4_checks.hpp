#pragma once

#include <string>

#include "warp_hypad/core/verification.hpp"

namespace program_b {
namespace checks {

struct B4Metrics {
    std::string fourier_convention;
    std::string maxwell_curl_signs;

    double worst_plane_wave_residual = 0.0;
    double worst_transversality_residual = 0.0;
    double fresnel_factorization_error = 0.0;
    double generic_birefringence_split = 0.0;
    double group_velocity_error = 0.0;

    double coarse_min_passivity = 0.0;
    double adaptive_min_passivity = 0.0;
    double adaptive_improvement = 0.0;
    int adaptive_samples = 0;

    double worst_energy_balance = 0.0;
    double reflection = 0.0;
    double transmission = 0.0;
    double absorption = 0.0;

    double phase_mismatch = 0.0;
    double group_delay_mismatch = 0.0;
    double exact_vs_restricted_error = 0.0;

    double hypad_gradient_error = 0.0;
    std::string branch_warnings = "none";
    bool oti_available = false;
};

void run_b4_checks(warp_hypad::VerificationReport& report,
                   const std::string& out_dir, B4Metrics& metrics);

}  // namespace checks
}  // namespace program_b
