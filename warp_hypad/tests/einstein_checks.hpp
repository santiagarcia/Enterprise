#pragma once

#include <string>

#include "warp_hypad/core/verification.hpp"

namespace warp_hypad {
namespace checks {

// Group "reference": the Einstein tensor of geometries that must have none.
// Writes einstein_minkowski.csv, einstein_flat_spherical.csv and
// einstein_schwarzschild.csv.
void run_einstein_reference_checks(VerificationReport& report,
                                   const std::string& out_dir);

// Group "desitter": the one case where G != 0 while the spacetime is still a
// vacuum, so it is the case that separates the geometric Einstein tensor from
// matter stress-energy. Writes einstein_desitter.csv.
void run_einstein_desitter_checks(VerificationReport& report,
                                  const std::string& out_dir);

// Group "alcubierre": selected components against the independent symbolic
// reference. Writes alcubierre_einstein_slice_xy.csv and
// alcubierre_einstein_centerline.csv.
void run_einstein_alcubierre_checks(VerificationReport& report,
                                    const std::string& out_dir);

// Group "identities": raw versus symmetric projection, the trace identity
// G^mu_mu = -R, and that neither tensor is silently replaced by the other.
void run_einstein_identity_checks(VerificationReport& report,
                                  const std::string& out_dir);

// Group "bianchi": the contracted Bianchi residual, absolute and relative.
void run_bianchi_checks(VerificationReport& report, const std::string& out_dir);

// Group "convergence": antisymmetry and Bianchi under refinement. Writes
// einstein_convergence.csv and bianchi_convergence.csv.
void run_einstein_convergence_checks(VerificationReport& report,
                                     const std::string& out_dir);

// Group "units": geometrized and SI conversion, and that the two never mix.
// Writes units_and_scaling.json.
void run_units_checks(VerificationReport& report, const std::string& out_dir);

// Group "policy": what A5 refuses, and that its outputs are deterministic.
void run_einstein_policy_checks(VerificationReport& report,
                                const std::string& out_dir);

}  // namespace checks
}  // namespace warp_hypad
