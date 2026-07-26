#pragma once
#include <string>
#include "warp_hypad/core/verification.hpp"

namespace program_b {
namespace checks {

// B1C: correction gate. Separates invariant metric realizability from
// observer-dependent 3+1 decomposition. Writes ergosurface_conditioning.csv
// and corrected_status_map.csv into out_dir.
// B1C step 1: four-dimensional Hodge duals, verified in isolation. Every
// later B1C case (Bel extraction A7, magnetoelectric B17, Z reconstruction
// A11) depends on these, so they run first.
void run_dual_checks(warp_hypad::VerificationReport& report);

// B1C steps 2 and 3: Bel extraction from (A7a)-(A7d), and the spatial
// projector with block spatiality. Gate B1C-S23.
void run_s23_checks(warp_hypad::VerificationReport& report);

// B1C step 5: pin the magnetoelectric normalisation via (B17).
void run_s5_checks(warp_hypad::VerificationReport& report);

// B1C steps 4, 5-remainder and 6. Gate B1C-S456.
//
// The quantities the gate is required to report, filled in as the checks run so
// the driver can print them without recomputing anything. Every one of them is
// also the recorded value of at least one case, so the summary cannot drift
// from the verdicts.
struct S456Metrics {
    double worst_direct_vs_block_H = 0.0;
    double worst_F_reconstruction = 0.0;
    double b17_zeta_error = 0.0;
    double b17_zeta_adjoint_error = 0.0;
    std::string handedness_result;
    double worst_A11_roundtrip = 0.0;
    double worst_observer_change = 0.0;
    double worst_coordinate_change = 0.0;
};

void run_s456_checks(warp_hypad::VerificationReport& report, S456Metrics& metrics);

// The rebuilt velocity-sensitivity cases, in their own translation unit because
// they carry an independent reference implementation of the metric that must not
// be visible to anything else.
struct VelocityMetrics {
    double worst_analytic = 0.0;      // OTI versus the analytic closed form
    double worst_fd = 0.0;            // analytic versus centred finite differences
    double worst_chain_rule = 0.0;    // the trajectory-term identity
};

void run_velocity_checks(warp_hypad::VerificationReport& report,
                         VelocityMetrics& metrics);

void run_b1c_checks(warp_hypad::VerificationReport& report,
                    const std::string& out_dir, VelocityMetrics& velocity);

}  // namespace checks
}  // namespace program_b
