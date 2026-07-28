#pragma once

#include <string>

#include "warp_hypad/core/verification.hpp"

namespace warp_hypad {
namespace checks {

// Group "frames": ADM recovery, tetrad orthonormality, handedness, continuity.
// Writes tetrad_verification.csv.
void run_frame_checks(VerificationReport& report, const std::string& out_dir);

// Group "observers": the four observer families and their admissibility rules,
// including the coordinate-stationary observer's failure at and beyond
// g_00 = 0. Writes observer_comparison.csv.
void run_observer_checks(VerificationReport& report, const std::string& out_dir);

// Group "fixtures": direct orthonormal-frame stress tensors with known
// energy-condition status, independent of the metric pipeline. Writes
// observer_reference_tensors.csv.
void run_fixture_checks(VerificationReport& report, const std::string& out_dir);

// Group "nec": the exact global minimiser on the sphere and its three
// independent cross-checks. Writes null_direction_samples.csv and
// worst_null_directions.csv.
void run_nec_checks(VerificationReport& report, const std::string& out_dir);

// Group "alcubierre": the Eulerian energy density against the closed form, and
// the energy-condition fields. Writes alcubierre_energy_slice_xy.csv,
// alcubierre_energy_centerline.csv, alcubierre_nec_slice_xy.csv,
// principal_stresses.csv and energy_condition_status.csv.
void run_alcubierre_energy_checks(VerificationReport& report,
                                  const std::string& out_dir);

// Group "integrals": E+, E-, E_abs, E_net and J_NEC, exact and smooth, with
// resolution and domain convergence. Writes integrated_energy_vs_velocity.csv,
// energy_condition_convergence.csv and domain_convergence.csv.
void run_integral_checks(VerificationReport& report, const std::string& out_dir);

}  // namespace checks
}  // namespace warp_hypad
