#pragma once

#include <string>
#include <vector>

#include "warp_hypad/core/verification.hpp"
#include "warp_hypad/design/constraint_projection.hpp"

// O0A verification groups. Each appends its cases to the shared report and
// writes its diagnostic records into out_dir; none terminates the process, so a
// failing gate shows everything that is broken rather than the first thing.

namespace warp_hypad {
namespace checks {

// --- The canonical O0A configuration ----------------------------------------
//
// Defined ONCE and shared by every group, so the raw, active, fixed and free
// coefficient counts that appear in the cases, the records, the readiness report
// and the figures are the same numbers. A group that built its own basis would
// report a self-consistent set of dimensions that quietly disagreed with the
// others.

inline constexpr double kO0aR0 = 2.0;

inline design::DesignRegion canonical_region() { return design::DesignRegion(); }

inline design::AxisymmetricBasis canonical_basis(
    design::FieldFamily family = design::FieldFamily::Scalar,
    bool protect_cabin = false) {
    return design::make_f1_basis(canonical_region(), kO0aR0, family, 10, 8,
                                 protect_cabin);
}

// A weighted integral of the field over a sub-region, assembled as constraint
// weights from the quadrature. One helper covers the first moment, the cabin
// mean and any global normalisation.
inline std::vector<double> functional_weights(
    const design::AxisymmetricBasis& b, const design::DesignQuadrature& q,
    double s_upper, bool weight_by_xi) {
    std::vector<double> w(static_cast<std::size_t>(b.coefficient_count()), 0.0);
    for (std::size_t p = 0; p < q.point.size(); ++p) {
        if (q.point[p].radial_squared > s_upper) continue;
        const design::TensorBasisEvaluation e = b.evaluate(q.point[p]);
        const double factor = weight_by_xi ? q.point[p].xi : 1.0;
        for (std::size_t k = 0; k < e.size(); ++k)
            w[static_cast<std::size_t>(e.flat_index[k])] +=
                q.weight[p] * factor * e.value[k];
    }
    return w;
}

// The canonical constraint set: two INDEPENDENT constraints plus one
// deliberately dependent row.
//
//   zero first xi-moment      removes a correction that is a pure translation
//   zero protected-cabin mean the cabin's mean value must not move
//   -3 x the moment row       a planted dependency, so rank detection has
//                             something real to find
//
// Rank 2 of 3 rows, so the free space is coefficient_count() - 2. With the
// mandated 6 x 6 = 36 active tensor products that is 34 independent optimization
// coordinates, which is what fixes the OTI gradient batching schedule.
inline design::ConstraintSystem canonical_constraints(
    const design::AxisymmetricBasis& basis,
    const design::DesignQuadrature& quad) {
    const design::DesignRegion region = canonical_region();
    design::ConstraintSystem sys(basis.coefficient_count());
    const std::vector<double> moment =
        functional_weights(basis, quad, region.s_max, /*weight_by_xi=*/true);
    sys.add_functional(moment, 0.0, "moment_constraint",
                       "zero first xi-moment over the region: removes a "
                       "correction that merely translates the bubble");
    sys.add_functional(
        functional_weights(basis, quad, region.protected_cabin_radius, false),
        0.0, "fixed_cabin_value", "zero change in the protected-cabin mean");
    std::vector<double> dependent = moment;
    for (double& v : dependent) v *= -3.0;
    sys.add_functional(dependent, 0.0, "moment_constraint",
                       "-3 times the first moment row: planted dependency, so "
                       "rank detection is exercised rather than assumed");
    return sys;
}

// Deterministic coefficients. No randomness anywhere in O0A: two runs must agree
// bit for bit, and a seeded generator would make that depend on the library.
inline std::vector<double> deterministic_coefficients(int n, double phase) {
    std::vector<double> c(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        c[static_cast<std::size_t>(i)] =
            0.37 * std::sin(2.11 * static_cast<double>(i) + phase) +
            0.11 * std::cos(0.83 * static_cast<double>(i));
    return c;
}

// --- Groups ------------------------------------------------------------------

// Cases 1-10: the univariate quintic kernel.
void run_o0a_kernel_checks(VerificationReport& report, const std::string& out_dir);

// Cases 11-20: the tensor product, the regular field families, the collars.
void run_o0a_tensor_checks(VerificationReport& report, const std::string& out_dir);

// Cases 21-28: constraints, projection, Gram conditioning, whitening.
void run_o0a_constraint_checks(VerificationReport& report,
                               const std::string& out_dir);

// Cases 29-35: sparse stencils, the F0 baseline projection, resolution.
void run_o0a_projection_checks(VerificationReport& report,
                               const std::string& out_dir);

// Cases 36-40: coefficient derivatives (OTI and finite difference), determinism,
// and the output checksums.
void run_o0a_derivative_checks(VerificationReport& report,
                               const std::string& out_dir);

// Written by the driver once every group has run, because readiness depends on
// results from several of them.
void write_o0a_readiness(const VerificationReport& report,
                         const std::string& out_dir);

// Animation frame records: production computes every frame so the renderer never
// re-implements the projection or the basis evaluation.
void write_o0a_frame_records(const std::string& out_dir);

// Written last: the checksum manifest, which must exclude itself from its own
// input set.
void write_o0a_manifest(const VerificationReport& report,
                        const std::string& out_dir);

}  // namespace checks
}  // namespace warp_hypad
