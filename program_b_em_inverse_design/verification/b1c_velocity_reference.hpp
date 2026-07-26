#pragma once

#include <cmath>

#include "warp_hypad/core/mat4.hpp"

// ---------------------------------------------------------------------------
// Independent reference implementations for the velocity-sensitivity tests.
//
// NOTHING HERE MAY CALL PRODUCTION CODE. No AlcubierreMetric, no
// shape_function, no shape_function_scalar_radius, no center_at, no
// center_at_scalar, no shape_function_radial_derivative, no OTI. The whole
// point is that these are written from the equations, so that agreement with
// the production evaluator is evidence rather than tautology. The only shared
// symbol is Mat4, which is a container.
//
// Metric, written out (Alcubierre, (-,+,+,+), c = 1):
//
//     g_00 = -1 + v^2 f^2       g_01 = g_10 = -v f
//     g_11 = g_22 = g_33 = 1    everything else 0
//
//     f(r_s) = [ tanh(sigma (r_s + R)) - tanh(sigma (r_s - R)) ]
//              / [ 2 tanh(sigma R) ]
//
//     r_s = sqrt( (x - x_s)^2 + y^2 + z^2 )
//
// R and sigma are fixed constants here, not arguments: the tests differentiate
// with respect to v only, and a signature carrying parameters that never vary
// invites them to be varied by accident.
//
// TWO DISTINCT DERIVATIVES, which is the whole subject:
//
//   FIXED CENTRE      x_s is held at a given value while v varies. Note that
//                     evaluate_fixed_center_reference() takes t and never uses
//                     it -- that absence IS the fixed-centre condition, and it
//                     is why a finite-difference sweep in v cannot move the
//                     bubble.
//
//   FIXED LABORATORY  x_s(v) = x_s0 + v t is recomputed inside every call, so a
//                     finite-difference sweep in v moves the bubble exactly as
//                     the physical trajectory does.
// ---------------------------------------------------------------------------

namespace program_b {
namespace checks {

inline constexpr double kRefR = 3.0;
inline constexpr double kRefSigma = 2.0;

// f(r)
inline double reference_shape(double r) {
    const double a = std::tanh(kRefSigma * (r + kRefR));
    const double b = std::tanh(kRefSigma * (r - kRefR));
    return (a - b) / (2.0 * std::tanh(kRefSigma * kRefR));
}

// df/dr, using sech^2(u) = 1 - tanh^2(u). Even f, so this vanishes at r = 0
// exactly -- not approximately, which is what makes the bubble centre a legal
// test point.
inline double reference_shape_derivative(double r) {
    const double a = std::tanh(kRefSigma * (r + kRefR));
    const double b = std::tanh(kRefSigma * (r - kRefR));
    return kRefSigma * ((1.0 - a * a) - (1.0 - b * b)) /
           (2.0 * std::tanh(kRefSigma * kRefR));
}

inline warp_hypad::Mat4<double> reference_metric_from_shape(double v, double f) {
    warp_hypad::Mat4<double> g;
    g(0, 0) = -1.0 + v * v * f * f;
    g(0, 1) = -v * f;
    g(1, 0) = -v * f;
    g(1, 1) = 1.0;
    g(2, 2) = 1.0;
    g(3, 3) = 1.0;
    return g;
}

// The bubble centre is GIVEN and does not move when v is perturbed. `t` is
// accepted so the two reference signatures are parallel, and is deliberately
// unused.
inline warp_hypad::Mat4<double> evaluate_fixed_center_reference(
    double v, double fixed_xs, double t, double x, double y, double z) {
    (void)t;
    const double dx = x - fixed_xs;
    const double r = std::sqrt(dx * dx + y * y + z * z);
    return reference_metric_from_shape(v, reference_shape(r));
}

// The bubble centre is RECOMPUTED from v on every call.
inline warp_hypad::Mat4<double> evaluate_fixed_lab_reference(double v, double xs0,
                                                            double t, double x,
                                                            double y, double z) {
    const double xs = xs0 + v * t;
    const double dx = x - xs;
    const double r = std::sqrt(dx * dx + y * y + z * z);
    return reference_metric_from_shape(v, reference_shape(r));
}

// --- Analytic derivatives ----------------------------------------------------
//
//   dg_00/df = 2 v^2 f        dg_01/df = dg_10/df = -v
//   dg_00/dv |_f = 2 v f^2    dg_01/dv |_f = dg_10/dv |_f = -f
//
// every other component is a constant and has zero derivative in v, which the
// gate asserts rather than assumes.

inline warp_hypad::Mat4<double> analytic_dgdv_fixed_center(double v, double fixed_xs,
                                                          double x, double y,
                                                          double z) {
    const double dx = x - fixed_xs;
    const double r = std::sqrt(dx * dx + y * y + z * z);
    const double f = reference_shape(r);
    warp_hypad::Mat4<double> d;
    d(0, 0) = 2.0 * v * f * f;
    d(0, 1) = -f;
    d(1, 0) = -f;
    return d;
}

// dr_s/dv = (dr_s/dx_s)(dx_s/dv) = ( -(x - x_s)/r_s )( t ).
//
// At r_s = 0 the radius is not differentiable, but the TRAJECTORY TERM is still
// exactly zero there: it always carries a factor df/dr, which vanishes at
// r = 0 because f is even, while |(x - x_s)/r_s| <= 1 stays bounded. Returning
// zero is therefore the value of the limit, not a regularisation.
inline double analytic_drs_dv(double xs, double t, double x, double y, double z) {
    const double dx = x - xs;
    const double r = std::sqrt(dx * dx + y * y + z * z);
    if (r < 1e-14) return 0.0;
    return -(dx / r) * t;
}

// The trajectory contribution, assembled component by component as
// (dg/df)(df/dr_s)(dr_s/dv).
inline warp_hypad::Mat4<double> analytic_trajectory_term(double v, double xs,
                                                         double t, double x,
                                                         double y, double z) {
    const double dx = x - xs;
    const double r = std::sqrt(dx * dx + y * y + z * z);
    warp_hypad::Mat4<double> d;
    if (r < 1e-14) return d;   // df/dr = 0 exactly at the centre

    const double f = reference_shape(r);
    const double dfdr = reference_shape_derivative(r);
    const double drdv = -(dx / r) * t;
    const double dfdv = dfdr * drdv;

    d(0, 0) = 2.0 * v * v * f * dfdv;   // (dg_00/df) dfdv
    d(0, 1) = -v * dfdv;                // (dg_01/df) dfdv
    d(1, 0) = -v * dfdv;
    return d;
}

inline warp_hypad::Mat4<double> analytic_dgdv_fixed_lab(double v, double xs0,
                                                       double t, double x,
                                                       double y, double z) {
    const double xs = xs0 + v * t;
    const warp_hypad::Mat4<double> p =
        analytic_dgdv_fixed_center(v, xs, x, y, z);
    const warp_hypad::Mat4<double> q = analytic_trajectory_term(v, xs, t, x, y, z);
    warp_hypad::Mat4<double> d;
    for (int i = 0; i < 16; ++i) d.m[i] = p.m[i] + q.m[i];
    return d;
}

}  // namespace checks
}  // namespace program_b
