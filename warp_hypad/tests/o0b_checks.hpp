#pragma once

#include <string>
#include <vector>

#include "warp_hypad/core/verification.hpp"
#include "warp_hypad/design/coefficient_spaces.hpp"
#include "warp_hypad/design/evaluation_routes.hpp"
#include "warp_hypad/design/shift_family.hpp"

// O0B verification groups. Each appends its cases to the shared report and
// writes its records into out_dir; none terminates the process.

namespace warp_hypad {
namespace checks {

// The canonical O0B configuration, defined once so every group, record and
// figure reports the same dimensions. The basis is O0A's, unchanged.
inline design::ShiftFamilyConfig o0b_config() {
    design::ShiftFamilyConfig c;
    c.velocity = 0.6;
    c.bubble_radius = 2.0;
    c.wall_sigma = 1.0;
    c.scale_R0 = 2.0;
    c.bubble_center_x0 = 0.0;
    return c;
}

inline design::DesignRegion o0b_region() { return design::DesignRegion(); }

inline design::AxisymmetricBasis o0b_basis() {
    return design::make_f1_basis(o0b_region(), o0b_config().scale_R0,
                                 design::FieldFamily::Scalar);
}

// The same canonical constraint set O0A used: moment, cabin mean, and one
// planted dependency. Rank 2 of 3, so 34 free coordinates and the established
// twelve-batch OTI schedule.
design::ConstraintSystem o0b_constraints(const design::AxisymmetricBasis& b,
                                         const design::DesignQuadrature& q);

inline std::vector<double> o0b_coefficients(int n, double amplitude) {
    std::vector<double> c(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        c[static_cast<std::size_t>(i)] =
            amplitude * std::sin(2.1 * static_cast<double>(i) + 0.4);
    return c;
}

// Cases 1-26 and 33-36: coefficient spaces, ADM assembly, references, jets.
void run_o0b_metric_checks(VerificationReport& report,
                           const std::string& out_dir);

// Cases 27-32 and 37-44: OTI coefficient gradients, rejections, routes,
// readiness and determinism.
void run_o0b_route_checks(VerificationReport& report,
                          const std::string& out_dir);

void run_o0bc_checks(VerificationReport& report, const std::string& out_dir);

void write_o0b_frame_records(const std::string& out_dir);
void write_o0b_readiness(const VerificationReport& report,
                         const std::string& out_dir);
void write_o0b_manifest(const VerificationReport& report,
                        const std::string& out_dir);
void write_o0a_consistency_record(const std::string& out_dir);

}  // namespace checks
}  // namespace warp_hypad
