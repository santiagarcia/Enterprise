#pragma once

#include <cassert>
#include <cmath>
#include <string>

#include "warp_hypad/core/mat4.hpp"
#include "warp_hypad/core/scalar_traits.hpp"

// ---------------------------------------------------------------------------
// Metric evaluation.
//
// Every metric exposes the same call:
//
//     Mat4<Scalar> evaluate(double t, double x, double y, double z,
//                           const Parameters& p) const;
//
// The time argument is present even where it is unused. That is deliberate: it
// is the seam through which a moving bubble is added later without changing a
// single call site, and it forces the frozen-time approximation to be visible
// in the signature rather than hidden in a comment.
//
// COORDINATES ARE double, PARAMETERS ARE Scalar. Grid coordinates carry no
// sensitivity -- differentiating with respect to x is a job for the
// finite-difference operators, not the AD scalar. The metric parameters (v, R,
// sigma) are the things HyPAD seeds. Mixing the two would make every grid
// coordinate a derivative channel and blow up the cost for nothing.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace metrics {

// How the metric treats time. Recorded in run metadata so no results file can
// be misread as a full spacetime calculation.
enum class TimeTreatment {
    // DIAGNOSTIC ONLY. The bubble centre is held fixed and d_t g_{mu nu} is
    // DISCARDED rather than computed.
    //
    // This is a spatial slice of a moving spacetime, and it is not a solution
    // of anything: imposing d_t g = 0 on a metric that genuinely depends on t
    // does not yield the metric of some other, static spacetime -- it yields an
    // inconsistent input to the Christoffel formula.
    //
    // Legitimate: metric components, determinant, signature, validity masks,
    // parameter sensitivities, regression-testing the spatial operators.
    // ILLEGITIMATE: any Christoffel symbol, curvature, stress-energy or energy
    // condition. The curvature code must REFUSE this mode, not warn about it.
    //
    // See docs/time_dependence_decision.md section 3.
    FrozenTimeDiagnostic,

    // RECOMMENDED. Inertial Cartesian coordinates with x_s(t) = v_s t.
    //
    // Coordinate observers are asymptotically the inertial rest observers, the
    // metric approaches eta at infinity, and d_t g_{mu nu} is genuinely
    // nonzero. For constant v_s it obeys the exact rigid-translation identity
    // d_t g = -v_s d_x g, so it costs no more than a spatial derivative.
    // See docs/time_dependence_decision.md sections 0 and 1.
    InertialTimeDependent,

    // Declared, not implemented. Co-moving coordinates x' = x - v_s t, in which
    // the metric is manifestly stationary but does NOT approach eta at infinity
    // and coordinate observers move with the bubble. Retained as the planned
    // Milestone 5 cross-check: curvature INVARIANTS must agree with
    // InertialTimeDependent while components differ.
    CoMoving
};

const char* to_string(TimeTreatment t);

// WHICH derivative d/dv_s means. The two differ whenever t != 0 and the
// difference is not small, so the choice is named rather than implied.
//
// See docs/time_dependence_decision.md and Program B's claims register.
enum class MetricSensitivityMode {
    // d/dv_s holding the bubble centre x_s fixed: the SHAPE sensitivity.
    // The trajectory is treated as a coordinate and carries no channel.
    // Appropriate when optimizing the bubble configuration itself.
    FixedBubbleCenter,

    // d/dv_s at fixed laboratory coordinates (t, x): the TOTAL derivative,
    //
    //   dg/dv|_{t,x} = dg/dv|_{r_s} + (dg/dr_s)(dr_s/dx_s)(dx_s/dv_s)
    //
    // with dx_s/dv_s = t. Appropriate for anything anchored to the laboratory:
    // fixed sensors, fixed hardware, location-dependent constraints.
    // This is the default for Program B measurement objectives.
    FixedLaboratoryCoordinates
};

const char* to_string(MetricSensitivityMode m);

// ---------------------------------------------------------------------------
// Minkowski
// ---------------------------------------------------------------------------

template <typename Scalar>
class MinkowskiMetric {
public:
    struct Parameters {};  // none

    static constexpr TimeTreatment time_treatment = TimeTreatment::FrozenTimeDiagnostic;
    static constexpr bool exact_time_independence = true;  // genuinely static

    Mat4<Scalar> evaluate(double t, double x, double y, double z,
                          const Parameters& p) const {
        (void)t; (void)x; (void)y; (void)z; (void)p;
        return Mat4<Scalar>::minkowski();
    }
};

// ---------------------------------------------------------------------------
// Alcubierre-type warp metric
// ---------------------------------------------------------------------------

template <typename Scalar>
struct AlcubierreParameters {
    Scalar velocity;  // v, bubble velocity (units of c)
    Scalar radius;    // R, bubble radius
    Scalar sigma;     // wall steepness

    AlcubierreParameters() = default;
    AlcubierreParameters(const Scalar& v, const Scalar& r, const Scalar& s)
        : velocity(v), radius(r), sigma(s) {}
};

// Alcubierre's original shape function:
//
//            tanh(sigma (r_s + R)) - tanh(sigma (r_s - R))
//   f(r_s) = ---------------------------------------------
//                        2 tanh(sigma R)
//
// f(0) = 1 exactly, and f -> 0 as r_s -> infinity, so the interior is a flat
// boosted region and the exterior is asymptotically flat. Smooth everywhere in
// the parameters, which is what makes it safe to differentiate with OTI.
//
// r_s is a double: it is built from grid coordinates and the (fixed) bubble
// centre, neither of which carries a derivative channel.
template <typename Scalar>
Scalar shape_function(double r_s, const AlcubierreParameters<Scalar>& p) {
    using std::tanh;  // ADL picks the OTIScalar/TestScalar overloads

    const Scalar r = scalar_cast<Scalar>(r_s);
    const Scalar two = scalar_cast<Scalar>(2.0);

    const Scalar a = tanh(p.sigma * (r + p.radius));
    const Scalar b = tanh(p.sigma * (r - p.radius));
    const Scalar d = two * tanh(p.sigma * p.radius);

    return (a - b) / d;
}

// Scalar-radius overload. Needed by MetricSensitivityMode::FixedLaboratoryCoordinates,
// where x_s(t) = v_s t carries a derivative channel and therefore so does r_s.
//
// AT r_s = 0 the radius is not differentiable (d(sqrt)/du diverges), but f is
// EVEN in r and f(0) = 1 identically for every R and sigma, so every derivative
// channel of f vanishes there. The guard returns the constant branch rather
// than propagating 0 * inf = NaN. This is exact, not a regularisation.
template <typename Scalar>
Scalar shape_function_scalar_radius(const Scalar& r_squared,
                                    const AlcubierreParameters<Scalar>& p) {
    using std::sqrt;
    using std::tanh;

    if (real_part(r_squared) <= 1e-24) {
        // f(0) = 1 exactly, with zero derivative in every channel.
        return scalar_cast<Scalar>(1.0);
    }

    const Scalar r = sqrt(r_squared);
    const Scalar two = scalar_cast<Scalar>(2.0);
    const Scalar a = tanh(p.sigma * (r + p.radius));
    const Scalar b = tanh(p.sigma * (r - p.radius));
    const Scalar d = two * tanh(p.sigma * p.radius);
    return (a - b) / d;
}

// df/dr_s, analytic. With sech^2(u) = 1 - tanh^2(u):
//
//     df/dr = sigma [ sech^2(sigma(r+R)) - sech^2(sigma(r-R)) ] / (2 tanh(sigma R))
//
// Note f is even in r and df/dr -> 0 as r -> 0, which is what makes the
// temporal derivative well defined at the bubble centre where dr/dt is
// undefined.
template <typename Scalar>
Scalar shape_function_radial_derivative(double r_s,
                                        const AlcubierreParameters<Scalar>& p) {
    using std::tanh;

    const Scalar r = scalar_cast<Scalar>(r_s);
    const Scalar one = scalar_cast<Scalar>(1.0);
    const Scalar two = scalar_cast<Scalar>(2.0);

    const Scalar ta = tanh(p.sigma * (r + p.radius));
    const Scalar tb = tanh(p.sigma * (r - p.radius));
    const Scalar d = two * tanh(p.sigma * p.radius);

    return p.sigma * ((one - ta * ta) - (one - tb * tb)) / d;
}

class AlcubierreMetricBase {
public:
    // Bubble centre at t = 0. A fixed coordinate, NOT a differentiated
    // parameter: seeding it would need a fourth OTI direction and is not part
    // of the A3 parameter set.
    double bubble_center_x = 0.0;

    TimeTreatment time_treatment = TimeTreatment::FrozenTimeDiagnostic;

    // Which d/dv_s the OTI channels represent. Defaults to the shape
    // sensitivity, which is what A3 verified; laboratory-anchored callers must
    // ask for the other one explicitly.
    MetricSensitivityMode sensitivity_mode =
        MetricSensitivityMode::FixedBubbleCenter;

    // True when this configuration may be fed to curvature code. Curvature
    // built from a discarded temporal derivative is meaningless, so the
    // Christoffel entry point will assert on this rather than warn.
    bool valid_for_curvature() const {
        return time_treatment == TimeTreatment::InertialTimeDependent;
    }
};

template <typename Scalar>
class AlcubierreMetric : public AlcubierreMetricBase {
public:
    using Parameters = AlcubierreParameters<Scalar>;

    AlcubierreMetric() = default;

    explicit AlcubierreMetric(double center_x,
                              TimeTreatment tt = TimeTreatment::FrozenTimeDiagnostic) {
        bubble_center_x = center_x;
        time_treatment = tt;
    }

    // Bubble centre at time t. This one function is the entire difference
    // between the time treatments, which is why it is factored out: a reader
    // can see the whole of the frozen-time approximation in three lines.
    double center_at(double t, const Parameters& p) const {
        switch (time_treatment) {
            case TimeTreatment::FrozenTimeDiagnostic:
                return bubble_center_x;  // t discarded -- diagnostic only
            case TimeTreatment::InertialTimeDependent:
                // x_s(t) = x_s(0) + v_s t, constant velocity.
                // NOTE this returns a double and therefore drops the channel.
                // FixedLaboratoryCoordinates does NOT route through here --
                // see center_at_scalar() and evaluate().
                //
                // real_part, not the Scalar velocity: the bubble centre is a
                // coordinate, and coordinates carry no derivative channels.
                // The v_s dependence of the METRIC still flows through
                // p.velocity in the components below; what is dropped here is
                // only d(x_s)/dv_s = t, which would be a derivative of the
                // coordinate system rather than of the field. Milestone 5 must
                // revisit this if sensitivities at fixed comoving position are
                // ever wanted.
                return bubble_center_x + real_part(p.velocity) * t;
            case TimeTreatment::CoMoving:
                return bubble_center_x;  // handled by CoMovingAlcubierreMetric
        }
        return bubble_center_x;
    }

    // ds^2 = -dt^2 + (dx - v f dt)^2 + dy^2 + dz^2
    //
    // expanded in the (-,+,+,+) signature:
    //
    //   g_00 = -1 + v^2 f^2
    //   g_01 = g_10 = -v f
    //   g_11 = g_22 = g_33 = 1
    //
    // det g = -1 identically for ANY v and f: the 2x2 (t,x) block has
    // determinant (-1 + v^2 f^2) - (v f)^2 = -1. That makes the determinant a
    // parameter-independent exact test rather than a tolerance one.
    // Bubble centre carrying the velocity derivative channel, for the
    // FixedLaboratoryCoordinates mode. Here dx_s/dv_s = t is propagated rather
    // than discarded.
    Scalar center_at_scalar(double t, const Parameters& p) const {
        if (time_treatment == TimeTreatment::InertialTimeDependent)
            return scalar_cast<Scalar>(bubble_center_x) +
                   p.velocity * scalar_cast<Scalar>(t);
        return scalar_cast<Scalar>(bubble_center_x);
    }

    Mat4<Scalar> evaluate(double t, double x, double y, double z,
                          const Parameters& p) const {
        assert(time_treatment != TimeTreatment::CoMoving &&
               "CoMoving is declared but not implemented; see "
               "docs/time_dependence_decision.md section 2");

        Scalar f;
        if (sensitivity_mode == MetricSensitivityMode::FixedLaboratoryCoordinates) {
            // r_s carries the channel, so the chain rule supplies the
            // trajectory term automatically -- no separate correction is added
            // by hand, which is the point of doing it this way.
            const Scalar dxs = scalar_cast<Scalar>(x) - center_at_scalar(t, p);
            const Scalar r2 = dxs * dxs + scalar_cast<Scalar>(y * y + z * z);
            f = shape_function_scalar_radius(r2, p);
        } else {
            const double dx = x - center_at(t, p);
            const double r_s = std::sqrt(dx * dx + y * y + z * z);
            f = shape_function(r_s, p);
        }

        const Scalar vf = p.velocity * f;

        Mat4<Scalar> g;
        g(0, 0) = scalar_cast<Scalar>(-1.0) + vf * vf;
        g(0, 1) = -vf;
        g(1, 0) = -vf;
        g(1, 1) = scalar_cast<Scalar>(1.0);
        g(2, 2) = scalar_cast<Scalar>(1.0);
        g(3, 3) = scalar_cast<Scalar>(1.0);
        return g;
    }

    // Analytic inverse, available in closed form because det g = -1:
    //
    //   g^00 = -1,  g^01 = -v f,  g^11 = 1 - v^2 f^2,  g^22 = g^33 = 1
    //
    // Provided as an INDEPENDENT reference for the numerical inverse, not as a
    // shortcut used by the pipeline. The pipeline inverts numerically so that
    // the same code path serves every metric.
    Mat4<Scalar> analytic_inverse(double t, double x, double y, double z,
                                  const Parameters& p) const {
        const double dx = x - center_at(t, p);
        const double r_s = std::sqrt(dx * dx + y * y + z * z);
        const Scalar f = shape_function(r_s, p);
        const Scalar vf = p.velocity * f;

        Mat4<Scalar> gi;
        gi(0, 0) = scalar_cast<Scalar>(-1.0);
        gi(0, 1) = -vf;
        gi(1, 0) = -vf;
        gi(1, 1) = scalar_cast<Scalar>(1.0) - vf * vf;
        gi(2, 2) = scalar_cast<Scalar>(1.0);
        gi(3, 3) = scalar_cast<Scalar>(1.0);
        return gi;
    }

    double radial_distance(double t, double x, double y, double z,
                           const Parameters& p) const {
        const double dx = x - center_at(t, p);
        return std::sqrt(dx * dx + y * y + z * z);
    }

    // Analytic d_t g_{mu nu}, exact.
    //
    // For a rigidly translating bubble every component depends on t and x only
    // through u = x - v_s t, so by the chain rule
    //
    //     d_t g = -v_s d_x g          (EXACT, constant v_s)
    //
    // This is why InertialTimeDependent costs nothing beyond a spatial
    // derivative, and it is the identity verified at Level A in the A2T suite
    // against an independently computed numerical d_x. See
    // docs/time_dependence_decision.md section 0.
    //
    // Implemented here by analytic differentiation of the closed form rather
    // than by invoking that identity, so the identity check compares two
    // genuinely independent routes.
    Mat4<Scalar> analytic_time_derivative(double t, double x, double y, double z,
                                          const Parameters& p) const {
        Mat4<Scalar> dg;  // zero
        if (time_treatment == TimeTreatment::FrozenTimeDiagnostic) return dg;

        using std::tanh;
        const double vs = real_part(p.velocity);
        const double dx = x - center_at(t, p);
        const double r_s = std::sqrt(dx * dx + y * y + z * z);

        // dr_s/dt = -v_s * dx / r_s. At r_s = 0 the radius is not
        // differentiable, but f'(0) = 0 there for this shape function (it is
        // even in r), so the product is zero and the limit is well defined.
        const Scalar f = shape_function(r_s, p);
        const Scalar dfdr = shape_function_radial_derivative(r_s, p);
        const double drdt = (r_s > 0.0) ? (-vs * dx / r_s) : 0.0;
        const Scalar dfdt = dfdr * scalar_cast<Scalar>(drdt);

        // g_00 = -1 + v^2 f^2  ->  d_t = 2 v^2 f df/dt
        dg(0, 0) = scalar_cast<Scalar>(2.0) * p.velocity * p.velocity * f * dfdt;
        // g_0x = -v f          ->  d_t = -v df/dt
        dg(0, 1) = -p.velocity * dfdt;
        dg(1, 0) = dg(0, 1);
        // spatial block is constant
        return dg;
    }
};

// --- Analytic references, double only --------------------------------------
//
// Closed-form parameter derivatives of the shape function, for verifying the
// HyPAD gradient. Written out independently from the expressions above rather
// than derived from them by any shared helper, so the comparison is meaningful.

struct ShapeDerivatives {
    double f;
    double df_dR;
    double df_dsigma;
};

ShapeDerivatives shape_function_derivatives(double r_s, double R, double sigma);

}  // namespace metrics
}  // namespace warp_hypad
