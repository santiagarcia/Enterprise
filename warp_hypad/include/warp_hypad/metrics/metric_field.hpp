#pragma once

#include <cmath>
#include <cstdint>
#include <vector>

#include "warp_hypad/core/field.hpp"
#include "warp_hypad/core/grid.hpp"
#include "warp_hypad/geometry/signature.hpp"
#include "warp_hypad/metrics/metric.hpp"

// ---------------------------------------------------------------------------
// Building a metric field over the grid, and checking that what was built is
// actually a metric.
//
// The validity mask is a first-class output rather than an assertion. A
// parameter sweep will inevitably wander into combinations that produce a
// degenerate or wrong-signature metric, and the useful response is a map
// showing exactly where -- not a crash, and certainly not silently continuing
// into a curvature calculation that means nothing.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace metrics {

// Bit flags, so one integer per point records every way a point can be bad.
enum ValidityFlag : std::uint32_t {
    kValid = 0u,
    kNonFinite = 1u << 0,
    kAsymmetric = 1u << 1,
    kSingular = 1u << 2,
    kIllConditioned = 1u << 3,
    kBadSignature = 1u << 4,
    kInverseResidual = 1u << 5,
    kEigenSolverFailed = 1u << 6,
};

struct ValidityReport {
    std::vector<std::uint32_t> flags;  // one per grid point

    std::size_t num_invalid = 0;
    std::size_t num_nonfinite = 0;
    std::size_t num_asymmetric = 0;
    std::size_t num_singular = 0;
    std::size_t num_bad_signature = 0;
    std::size_t num_ill_conditioned = 0;

    double worst_symmetry_residual = 0.0;
    double worst_inverse_residual = 0.0;
    double worst_det_error = 0.0;      // vs the expected determinant
    double min_abs_determinant = 1e300;

    bool all_valid() const { return num_invalid == 0; }
};

// Fill a rank-2 field from any metric exposing evaluate(t,x,y,z,params).
template <typename Scalar, typename MetricT>
Rank2Field<Scalar> build_metric_field(
    const Grid& grid, const MetricT& metric,
    const typename MetricT::Parameters& params, double t = 0.0) {
    Rank2Field<Scalar> field(grid);
    for (int i = 0; i < grid.nx(); ++i)
        for (int j = 0; j < grid.ny(); ++j)
            for (int k = 0; k < grid.nz(); ++k) {
                const std::size_t p = grid.flat(i, j, k);
                field.set_matrix(
                    p, metric.evaluate(t, grid.x(i), grid.y(j), grid.z(k), params));
            }
    return field;
}

// Also fill the inverse metric, numerically, via the shared Mat4 inverse.
template <typename Scalar>
Rank2Field<Scalar> build_inverse_field(const Grid& grid,
                                       const Rank2Field<Scalar>& g,
                                       ValidityReport* report = nullptr) {
    Rank2Field<Scalar> out(grid);
    const std::size_t n = grid.num_points();
    if (report) report->flags.assign(n, kValid);

    for (std::size_t p = 0; p < n; ++p) {
        const Mat4<Scalar> m = g.matrix_at(p);
        const InverseResult<Scalar> inv = inverse(m);
        if (inv.ok) out.set_matrix(p, inv.inverse);

        if (report) {
            std::uint32_t f = kValid;
            if (!inv.ok) f |= kSingular;
            if (inv.ill_conditioned) f |= kIllConditioned;
            report->flags[p] |= f;
        }
    }
    return out;
}

// Full validity sweep. `expected_determinant` is checked when finite; for the
// Alcubierre metric it is exactly -1 for every parameter combination, which
// makes it a strong pointwise test.
template <typename Scalar>
ValidityReport check_metric_validity(const Grid& grid,
                                     const Rank2Field<Scalar>& g,
                                     double expected_determinant,
                                     double inverse_tolerance = 1e-10) {
    ValidityReport rep;
    const std::size_t n = grid.num_points();
    rep.flags.assign(n, kValid);

    for (std::size_t p = 0; p < n; ++p) {
        const Mat4<Scalar> m = g.matrix_at(p);
        std::uint32_t flag = kValid;

        Mat4<double> mr;
        bool finite = true;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) {
                const double v = real_part(m(a, b));
                mr(a, b) = v;
                if (!std::isfinite(v)) finite = false;
            }
        if (!finite) flag |= kNonFinite;

        const double sym = symmetry_residual(m);
        rep.worst_symmetry_residual = std::fmax(rep.worst_symmetry_residual, sym);
        if (sym > 1e-14) flag |= kAsymmetric;

        const InverseResult<Scalar> inv = inverse(m);
        if (!inv.ok) {
            flag |= kSingular;
        } else {
            const double res = identity_residual(m, inv.inverse);
            rep.worst_inverse_residual = std::fmax(rep.worst_inverse_residual, res);
            if (res > inverse_tolerance) flag |= kInverseResidual;
            if (inv.ill_conditioned) flag |= kIllConditioned;
        }

        const double det_abs = std::fabs(real_part(det(m)));
        rep.min_abs_determinant = std::fmin(rep.min_abs_determinant, det_abs);
        if (std::isfinite(expected_determinant)) {
            rep.worst_det_error = std::fmax(
                rep.worst_det_error,
                std::fabs(real_part(det(m)) - expected_determinant));
        }

        if (finite) {
            const geometry::SignatureResult sig = geometry::analyse_signature(mr);
            if (!sig.converged) flag |= kEigenSolverFailed;
            if (!sig.lorentzian) flag |= kBadSignature;
        }

        rep.flags[p] = flag;
        if (flag != kValid) {
            ++rep.num_invalid;
            if (flag & kNonFinite) ++rep.num_nonfinite;
            if (flag & kAsymmetric) ++rep.num_asymmetric;
            if (flag & kSingular) ++rep.num_singular;
            if (flag & kBadSignature) ++rep.num_bad_signature;
            if (flag & kIllConditioned) ++rep.num_ill_conditioned;
        }
    }
    return rep;
}

}  // namespace metrics
}  // namespace warp_hypad
