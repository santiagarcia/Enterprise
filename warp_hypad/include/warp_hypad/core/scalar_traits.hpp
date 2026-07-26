#pragma once

#include <cmath>
#include <limits>

// ---------------------------------------------------------------------------
// The scalar seam.
//
// Every numerical kernel in this project is templated on Scalar. Kernels are
// permitted to use only:
//
//   * the arithmetic operators + - * / and unary -,
//   * ADL-visible sqrt/fabs/exp/log/tanh/... where the physics needs them,
//   * default construction and copy,
//   * ScalarTraits<Scalar>, declared below.
//
// Kernels are NOT permitted to:
//
//   * name `double` for any quantity that participates in the physics,
//   * compare two Scalars with < > <= >=,
//   * branch on a Scalar value.
//
// The last two matter for hypercomplex automatic differentiation: a comparison
// on a hypercomplex number is either ill-defined or silently discards the
// derivative channels, and a branch introduces a derivative discontinuity that
// the AD type cannot see. Wherever an ordering decision is genuinely required
// (singularity detection, tolerance checks, reporting) the code goes through
// real_part() and makes the decision on the real component only, explicitly.
//
// Milestone 10 specialises this trait for HyPADScalar. Nothing else changes.
// ---------------------------------------------------------------------------

namespace warp_hypad {

template <typename Scalar>
struct ScalarTraits {
    // The real (derivative-free) component. Used only for comparisons,
    // tolerances, diagnostics and I/O -- never inside a physics expression.
    static double real(const Scalar& x) { return static_cast<double>(x); }

    // Lift a compile-time or configuration constant into the scalar type.
    // Constants carry no parameter sensitivity, so this seeds zero derivative.
    static Scalar from_double(double x) { return static_cast<Scalar>(x); }

    // True when the type carries derivative channels alongside its real part.
    static constexpr bool is_differential = false;

    // Number of independent derivative channels (0 for plain arithmetic).
    static constexpr int derivative_channels = 0;
};

template <typename Scalar>
inline double real_part(const Scalar& x) {
    return ScalarTraits<Scalar>::real(x);
}

template <typename Scalar>
inline Scalar scalar_cast(double x) {
    return ScalarTraits<Scalar>::from_double(x);
}

// Magnitude of the real part. The natural quantity for residual reporting:
// a residual is a diagnostic about the primal solution, not a differentiable
// output, so taking |real| here is correct rather than lossy.
template <typename Scalar>
inline double real_abs(const Scalar& x) {
    const double r = real_part(x);
    return r < 0.0 ? -r : r;
}

}  // namespace warp_hypad
