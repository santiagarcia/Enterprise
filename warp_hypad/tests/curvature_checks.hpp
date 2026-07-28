#pragma once

#include <string>

#include "warp_hypad/core/verification.hpp"

namespace warp_hypad {
namespace checks {

// Milestone A4, group "connection": Christoffel assembly, metric compatibility,
// and the index symmetries of the Riemann tensor. Writes
// flat_spherical_christoffel_error.csv into out_dir.
void run_connection_checks(VerificationReport& report, const std::string& out_dir);

// Group "analytic": the closed-form spacetimes. Writes
// schwarzschild_kretschmann_profile.csv and desitter_invariants.csv.
void run_analytic_curvature_checks(VerificationReport& report,
                                   const std::string& out_dir);

// Group "temporal": what the curvature pipeline refuses, and the exact
// rigid-translation identity relating the resampled temporal derivative to the
// spatial one.
void run_curvature_temporal_checks(VerificationReport& report,
                                   const std::string& out_dir);

// Group "alcubierre": the time-dependent case, against the independently
// generated symbolic reference. Writes alcubierre_curvature_slice_xy.csv and
// alcubierre_centerline_x.csv.
void run_alcubierre_curvature_checks(VerificationReport& report,
                                     const std::string& out_dir);

// Group "convergence": grid-refinement study, interior and boundary separately.
// Writes curvature_convergence.csv.
void run_curvature_convergence_checks(VerificationReport& report,
                                      const std::string& out_dir);

}  // namespace checks
}  // namespace warp_hypad
