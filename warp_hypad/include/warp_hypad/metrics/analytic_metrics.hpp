#pragma once

#include <cmath>

#include "warp_hypad/core/mat4.hpp"
#include "warp_hypad/core/scalar_traits.hpp"
#include "warp_hypad/core/tensor.hpp"
#include "warp_hypad/geometry/chart.hpp"

// ---------------------------------------------------------------------------
// Closed-form reference spacetimes for Milestone A4.
//
// These exist to verify the curvature pipeline, not to be optimized. Each one
// pins down a different failure mode, and the set is chosen so that no single
// bug can pass all four:
//
//   Minkowski, Cartesian      Gamma = 0 and Riemann = 0. Catches almost nothing
//                             on its own -- a pipeline that returns zero always
//                             passes it -- which is exactly why the next case
//                             exists.
//
//   flat space, spherical     Gamma != 0 and Riemann = 0. Separates "the
//                             connection vanishes" from "the curvature
//                             vanishes". This is the case that catches a sign
//                             error in the Gamma-Gamma terms of the Riemann
//                             assembly, because those terms are individually
//                             large and must cancel against the dGamma terms.
//
//   Schwarzschild             R_{mu nu} = 0 and R = 0 while the Kretschmann
//                             invariant is 48 M^2 / r^6. Separates "Ricci flat"
//                             from "flat", and supplies a nonzero closed form
//                             with a steep radial profile for the convergence
//                             study.
//
//   static de Sitter          R_{mu nu} = (3/L^2) g_{mu nu}, R = 12/L^2,
//                             K = 24/L^4. The only case here with nonzero Ricci,
//                             so it fixes the sign of the Ricci contraction and
//                             therefore the sign of the energy density that
//                             Milestone A6 will report. Under the MTW (+,+,+)
//                             conventions of docs/conventions.md section 6, R is
//                             POSITIVE for de Sitter; a pipeline with the
//                             opposite Riemann sign returns -12/L^2 and would
//                             otherwise sail through every identity check.
//
// Schwarzschild appears twice, in Schwarzschild and in isotropic coordinates, so
// that a scalar invariant can be compared between two charts describing the same
// curved geometry. Chart independence of K is the strongest statement available
// about the pipeline short of cross-code validation.
//
// COORDINATE SINGULARITIES. The spherical and isotropic charts are singular at
// r = 0 and on the polar axis sin(theta) = 0, and Schwarzschild additionally at
// the horizon r = 2M. Those are properties of the chart, not of the spacetime,
// and a grid that straddles one produces large finite numbers that look like
// physics. Every domain used with these metrics must exclude them; the A4 suite
// states its domains explicitly and the chart records that it has an axis.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace metrics {

// --- Flat space in spherical coordinates -------------------------------------
//
//   ds^2 = -dt^2 + dr^2 + r^2 dtheta^2 + r^2 sin^2(theta) dphi^2
//
// Grid axes are (r, theta, phi).

template <typename Scalar>
class FlatSphericalMetric {
public:
    struct Parameters {};

    static constexpr geometry::ChartKind chart = geometry::ChartKind::SphericalTRThetaPhi;

    Mat4<Scalar> evaluate(double t, double r, double theta, double phi,
                          const Parameters& p) const {
        (void)t; (void)phi; (void)p;
        const double s = std::sin(theta);
        Mat4<Scalar> g;
        g(0, 0) = scalar_cast<Scalar>(-1.0);
        g(1, 1) = scalar_cast<Scalar>(1.0);
        g(2, 2) = scalar_cast<Scalar>(r * r);
        g(3, 3) = scalar_cast<Scalar>(r * r * s * s);
        return g;
    }
};

// Exact Christoffel symbols of the flat spherical chart, written out from the
// standard result rather than derived from the code under test. Index map
// 0 = t, 1 = r, 2 = theta, 3 = phi.
//
//   Gamma^r_{theta theta}  = -r
//   Gamma^r_{phi phi}      = -r sin^2(theta)
//   Gamma^theta_{r theta}  = 1/r
//   Gamma^theta_{phi phi}  = -sin(theta) cos(theta)
//   Gamma^phi_{r phi}      = 1/r
//   Gamma^phi_{theta phi}  = cot(theta)
//
// every other component zero.
inline Rank3<double> flat_spherical_christoffel(double r, double theta) {
    const double s = std::sin(theta);
    const double c = std::cos(theta);
    Rank3<double> G;
    G(1, 2, 2) = -r;
    G(1, 3, 3) = -r * s * s;
    G(2, 1, 2) = 1.0 / r;
    G(2, 2, 1) = 1.0 / r;
    G(2, 3, 3) = -s * c;
    G(3, 1, 3) = 1.0 / r;
    G(3, 3, 1) = 1.0 / r;
    G(3, 2, 3) = c / s;
    G(3, 3, 2) = c / s;
    return G;
}

// --- Schwarzschild, Schwarzschild coordinates --------------------------------
//
//   ds^2 = -(1 - 2M/r) dt^2 + (1 - 2M/r)^-1 dr^2
//          + r^2 dtheta^2 + r^2 sin^2(theta) dphi^2
//
// Grid axes are (r, theta, phi). Valid for r > 2M only.

template <typename Scalar>
struct SchwarzschildParameters {
    Scalar mass;
    SchwarzschildParameters() = default;
    explicit SchwarzschildParameters(const Scalar& m) : mass(m) {}
};

template <typename Scalar>
class SchwarzschildMetric {
public:
    using Parameters = SchwarzschildParameters<Scalar>;

    static constexpr geometry::ChartKind chart = geometry::ChartKind::SphericalTRThetaPhi;

    Mat4<Scalar> evaluate(double t, double r, double theta, double phi,
                          const Parameters& p) const {
        (void)t; (void)phi;
        const double s = std::sin(theta);
        const Scalar one = scalar_cast<Scalar>(1.0);
        const Scalar two_over_r = scalar_cast<Scalar>(2.0 / r);
        const Scalar lapse2 = one - two_over_r * p.mass;  // 1 - 2M/r

        Mat4<Scalar> g;
        g(0, 0) = -lapse2;
        g(1, 1) = one / lapse2;
        g(2, 2) = scalar_cast<Scalar>(r * r);
        g(3, 3) = scalar_cast<Scalar>(r * r * s * s);
        return g;
    }

    // Kretschmann invariant, exact: K = 48 M^2 / r^6.
    static double kretschmann(double r, double mass) {
        const double r3 = r * r * r;
        return 48.0 * mass * mass / (r3 * r3);
    }
};

// --- Schwarzschild, isotropic coordinates ------------------------------------
//
//   ds^2 = -((1 - M/2rho)/(1 + M/2rho))^2 dt^2
//          + (1 + M/2rho)^4 ( drho^2 + rho^2 dtheta^2 + rho^2 sin^2(theta) dphi^2 )
//
// The areal radius is r = rho (1 + M/2rho)^2, so the SAME geometry is described
// with different components. Valid for rho > M/2, which is the horizon.

template <typename Scalar>
class SchwarzschildIsotropicMetric {
public:
    using Parameters = SchwarzschildParameters<Scalar>;

    static constexpr geometry::ChartKind chart = geometry::ChartKind::IsotropicTRhoThetaPhi;

    Mat4<Scalar> evaluate(double t, double rho, double theta, double phi,
                          const Parameters& p) const {
        (void)t; (void)phi;
        const double s = std::sin(theta);
        const Scalar one = scalar_cast<Scalar>(1.0);
        const Scalar half_over_rho = scalar_cast<Scalar>(0.5 / rho);
        const Scalar u = half_over_rho * p.mass;  // M / (2 rho)

        const Scalar minus = one - u;
        const Scalar plus = one + u;
        const Scalar a = minus / plus;      // lapse
        const Scalar b2 = plus * plus;      // (1 + M/2rho)^2
        const Scalar b4 = b2 * b2;          // conformal factor

        Mat4<Scalar> g;
        g(0, 0) = -(a * a);
        g(1, 1) = b4;
        g(2, 2) = b4 * scalar_cast<Scalar>(rho * rho);
        g(3, 3) = b4 * scalar_cast<Scalar>(rho * rho * s * s);
        return g;
    }

    // Areal radius corresponding to isotropic rho.
    static double areal_radius(double rho, double mass) {
        const double f = 1.0 + mass / (2.0 * rho);
        return rho * f * f;
    }

    // K expressed in the isotropic coordinate, via the areal radius.
    static double kretschmann(double rho, double mass) {
        const double r = areal_radius(rho, mass);
        const double r3 = r * r * r;
        return 48.0 * mass * mass / (r3 * r3);
    }
};

// --- Static de Sitter --------------------------------------------------------
//
//   ds^2 = -(1 - r^2/L^2) dt^2 + (1 - r^2/L^2)^-1 dr^2
//          + r^2 dtheta^2 + r^2 sin^2(theta) dphi^2
//
// Grid axes are (r, theta, phi). Valid for 0 < r < L; r = L is the cosmological
// horizon, a coordinate singularity of this chart.
//
// A maximally symmetric spacetime, so every curvature scalar is a constant. That
// makes it the one reference here whose expected values do not vary across the
// grid, which in turn makes a radial trend in the measured error unambiguously a
// property of the discretisation rather than of the solution.

template <typename Scalar>
struct DeSitterParameters {
    Scalar hubble_length;  // L, with Lambda = 3/L^2
    DeSitterParameters() = default;
    explicit DeSitterParameters(const Scalar& l) : hubble_length(l) {}
};

template <typename Scalar>
class DeSitterStaticMetric {
public:
    using Parameters = DeSitterParameters<Scalar>;

    static constexpr geometry::ChartKind chart = geometry::ChartKind::SphericalTRThetaPhi;

    Mat4<Scalar> evaluate(double t, double r, double theta, double phi,
                          const Parameters& p) const {
        (void)t; (void)phi;
        const double s = std::sin(theta);
        const Scalar one = scalar_cast<Scalar>(1.0);
        const Scalar r2 = scalar_cast<Scalar>(r * r);
        const Scalar l2 = p.hubble_length * p.hubble_length;
        const Scalar lapse2 = one - r2 / l2;

        Mat4<Scalar> g;
        g(0, 0) = -lapse2;
        g(1, 1) = one / lapse2;
        g(2, 2) = r2;
        g(3, 3) = scalar_cast<Scalar>(r * r * s * s);
        return g;
    }

    // Exact invariants, MTW (+,+,+) conventions with signature (-,+,+,+):
    //
    //   R_{mu nu} = (3/L^2) g_{mu nu},  R = 12/L^2,
    //   R_{mu nu} R^{mu nu} = 36/L^4,   K = 24/L^4.
    static double ricci_scalar(double l) { return 12.0 / (l * l); }
    static double ricci_square(double l) {
        const double l2 = l * l;
        return 36.0 / (l2 * l2);
    }
    static double kretschmann(double l) {
        const double l2 = l * l;
        return 24.0 / (l2 * l2);
    }
    static double lambda(double l) { return 3.0 / (l * l); }
};

}  // namespace metrics
}  // namespace warp_hypad
