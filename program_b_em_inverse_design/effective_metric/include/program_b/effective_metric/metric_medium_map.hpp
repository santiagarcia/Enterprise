#pragma once

#include "program_b/constitutive/constitutive.hpp"

// ---------------------------------------------------------------------------
// The metric <-> medium correspondence.
//
// FORWARD (Plebanski-Tamm). For a spacetime metric g_{mu nu} with g_00 < 0,
// Maxwell's equations in that spacetime take the flat-space form with
//
//     eps^{ij} = mu^{ij} = - sqrt(-g) g^{ij} / g_00
//     w_i      = - g_{0i} / g_00
//     D = eps E + w x H
//     B = mu  H - w x E
//
// so that xi = cross_matrix(w) and zeta = -cross_matrix(w).
//
// Reference: J. Plebanski, "Electromagnetic waves in gravitational fields",
// Phys. Rev. 118, 1396 (1960), building on earlier work by Tamm; the pairing is
// commonly called the Plebanski-Tamm constitutive relations. The form is also
// reproduced in U. Leonhardt and T. G. Philbin, "General relativity in
// electrical engineering", New J. Phys. 8, 247 (2006).
//
// ON THE SIGN CONVENTION. Sources differ on the sign of w and on whether the
// magnetoelectric term is written with a cross product or an explicit
// Levi-Civita contraction, and the bibliographic details above are recorded
// from memory and should be checked against the originals before being cited in
// a paper. Rather than inherit a sign this project cannot verify, the
// convention above is FIXED HERE and pinned by three checks that do not depend
// on any external convention:
//
//   1. Minkowski must map to exact vacuum (eps = mu = I, w = 0). This fixes the
//      overall normalisation and the sign of eps.
//   2. The map must be invariant under g -> Omega^2 g. This is a property of
//      the formula itself and holds for either sign of w, so it does not fix
//      the sign but does confirm the sqrt(-g) power.
//   3. Round-tripping a metric with g_{0i} = 0 must return the conformal class
//      of the original.
//
// What the sign of w DOES affect is the direction of the magnetoelectric drag.
// Until a propagation test is implemented (Milestone B2), this module makes no
// claim about the physical handedness of that term; see
// docs/effective_metric_theory.md section 7 and the limitations list.
//
// CONFORMAL AMBIGUITY. The forward map is invariant under g -> Omega^2 g:
//
//     sqrt(-g) -> Omega^4 sqrt(-g),  g^{ij} -> Omega^-2 g^{ij},
//     g_00 -> Omega^2 g_00
//
// so eps^{ij} -> Omega^{4-2-2} eps^{ij} = eps^{ij}, and w_i is a ratio of two
// factors of Omega^2 and is likewise unchanged. A medium therefore determines
// its metric only up to a conformal factor. The inverse map below returns an
// explicit representative together with the gauge used to pick it, rather than
// silently normalising.
// ---------------------------------------------------------------------------

namespace program_b {

// Which representative of the conformal class the inverse map returns.
enum class ConformalGauge {
    // g_00 = -1 at the evaluation point. Always available when the medium came
    // from a metric with a timelike coordinate direction.
    UnitLapse,
};

const char* to_string(ConformalGauge g);

// --- Forward -----------------------------------------------------------------

template <typename Scalar>
struct ForwardResult {
    ConstitutiveResponse<Scalar> response;
    Scalar w[3];

    // False when g_00 >= 0, i.e. the coordinate time direction is not timelike
    // and no medium in this frame reproduces the metric. For the Alcubierre
    // metric this happens exactly where v_s f > 1 -- the same condition that
    // costs Program A its coordinate observers.
    bool ok = false;
    double g00_real = 0.0;
    std::string obstruction;
};

template <typename Scalar>
ForwardResult<Scalar> metric_to_medium(const Mat4<Scalar>& g) {
    ForwardResult<Scalar> out;
    out.g00_real = real_part(g(0, 0));

    if (!(out.g00_real < 0.0)) {
        out.obstruction =
            "g_00 >= 0: the coordinate time direction is not timelike here, so "
            "there is no equivalent medium in this frame. For a warp metric "
            "this is the region where v_s f > 1.";
        return out;
    }

    const warp_hypad::InverseResult<Scalar> gi = warp_hypad::inverse(g);
    if (!gi.ok) {
        out.obstruction = "the metric is singular or too ill-conditioned to invert here.";
        return out;
    }

    const Scalar detg = gi.determinant;
    if (!(real_part(detg) < 0.0)) {
        out.obstruction =
            "det g >= 0: not a Lorentzian metric in the (-,+,+,+) signature.";
        return out;
    }

    using std::sqrt;
    const Scalar root = sqrt(-detg);            // sqrt(-g)
    const Scalar factor = -root / g(0, 0);      // -sqrt(-g)/g_00, positive

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            // g^{ij} is the SPATIAL block of the inverse metric, indices 1..3.
            const Scalar v = factor * gi.inverse(i + 1, j + 1);
            out.response.epsilon(i, j) = v;
            out.response.mu(i, j) = v;
        }

    for (int i = 0; i < 3; ++i) out.w[i] = -g(0, i + 1) / g(0, 0);

    out.response.xi = cross_matrix<Scalar>(out.w);
    out.response.zeta = mat3_scale(out.response.xi, scalar_cast<Scalar>(-1.0));

    out.ok = true;
    return out;
}

// --- Inverse -----------------------------------------------------------------

enum class InverseStatus {
    Ok,
    NotMetricRepresentable,   // classify() rejected it
    MagnetoelectricUnsupported,  // w != 0: valid but not implemented in B1
    Degenerate,               // eps not invertible
};

const char* to_string(InverseStatus s);

template <typename Scalar>
struct InverseResultMetric {
    Mat4<Scalar> metric;
    InverseStatus status = InverseStatus::NotMetricRepresentable;
    ConformalGauge gauge = ConformalGauge::UnitLapse;
    ResponseClassification classification;
    std::string message;

    bool ok() const { return status == InverseStatus::Ok; }
};

// Reconstruct a representative metric from a medium.
//
// DERIVATION for the implemented case (w = 0, gauge g_00 = -1). The forward map
// with g_00 = -1 and g_{0i} = 0 gives, writing h for the spatial block,
//
//     sqrt(-g) = sqrt(det h),   g^{ij} = (h^{-1})^{ij}
//     eps = sqrt(det h) h^{-1}
//
// Taking determinants of both sides:
//
//     det eps = (det h)^{3/2} / det h = (det h)^{1/2}
//
// so sqrt(det h) = det eps, and substituting back
//
//     h = det(eps) * eps^{-1}
//
// which is exact and closed form. Checks: eps = I gives h = I, recovering
// Minkowski; eps = n^2 I gives h = n^4 I, whose null cones travel at
// 1/n^2 = 1/sqrt(eps mu) as they must for an impedance-matched medium with
// eps = mu = n^2.
//
// The returned metric is ONE REPRESENTATIVE of the conformal class. Any
// Omega^2 g maps to the same medium, so componentwise comparison against an
// original metric is the wrong test; compare conformal invariants instead.
template <typename Scalar>
InverseResultMetric<Scalar> medium_to_metric(
    const ConstitutiveResponse<Scalar>& r, double tolerance = 1e-10) {
    InverseResultMetric<Scalar> out;
    out.classification = classify(r, tolerance);

    if (!out.classification.metric_representable) {
        out.status = InverseStatus::NotMetricRepresentable;
        out.message = out.classification.obstruction;
        return out;
    }

    if (mat3_max_abs(r.xi) > tolerance) {
        out.status = InverseStatus::MagnetoelectricUnsupported;
        out.message =
            "the medium has a nonzero magnetoelectric coupling (w != 0). This "
            "IS metric-representable -- it corresponds to g_{0i} != 0 -- but "
            "the closed-form reconstruction implemented in B1 covers only "
            "w = 0. Reporting rather than approximating.";
        return out;
    }

    const Mat3Inverse<Scalar> ei = mat3_inverse(r.epsilon);
    if (!ei.ok) {
        out.status = InverseStatus::Degenerate;
        out.message = "eps is singular; no metric representative exists.";
        return out;
    }

    if (!(ei.det_real > 0.0)) {
        out.status = InverseStatus::Degenerate;
        out.message =
            "det eps <= 0: the spatial block would not be positive definite, so "
            "the result would not have signature (-,+,+,+).";
        return out;
    }

    const Mat3<Scalar> h = mat3_scale(ei.inverse, ei.determinant);

    out.metric = Mat4<Scalar>();
    out.metric(0, 0) = scalar_cast<Scalar>(-1.0);  // the UnitLapse gauge
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) out.metric(i + 1, j + 1) = h(i, j);

    out.status = InverseStatus::Ok;
    out.gauge = ConformalGauge::UnitLapse;
    out.message = "representative in the g_00 = -1 gauge; the medium determines "
                  "the metric only up to a conformal factor.";
    return out;
}

// --- Conformal comparison ----------------------------------------------------
//
// Two metrics are conformally equivalent iff g2 = Omega^2 g1 for some positive
// Omega. Componentwise equality is the WRONG test after a round trip, because
// the forward map discards exactly that factor.
//
// The estimator below recovers Omega^2 from the ratio of the largest-magnitude
// component pair and then reports the worst residual of g2 - Omega^2 g1. It
// returns the factor as well, so a caller can see what was removed rather than
// having it silently divided out.

template <typename Scalar>
struct ConformalComparison {
    double omega_squared = 0.0;
    double residual = 1e300;
    bool conformally_equal = false;
};

template <typename Scalar>
ConformalComparison<Scalar> compare_up_to_conformal(const Mat4<Scalar>& a,
                                                    const Mat4<Scalar>& b,
                                                    double tolerance = 1e-10) {
    ConformalComparison<Scalar> out;

    int best = -1;
    double best_mag = 0.0;
    for (int i = 0; i < 16; ++i) {
        const double m = std::fabs(real_part(a.m[i]));
        if (m > best_mag) { best_mag = m; best = i; }
    }
    if (best < 0 || best_mag <= 0.0) return out;

    out.omega_squared = real_part(b.m[best]) / real_part(a.m[best]);
    if (!(out.omega_squared > 0.0)) return out;  // orientation-reversing: not conformal

    double worst = 0.0;
    double scale = 0.0;
    for (int i = 0; i < 16; ++i) {
        const double predicted = out.omega_squared * real_part(a.m[i]);
        worst = std::fmax(worst, std::fabs(real_part(b.m[i]) - predicted));
        scale = std::fmax(scale, std::fabs(real_part(b.m[i])));
    }
    out.residual = worst / std::fmax(scale, 1e-300);
    out.conformally_equal = out.residual <= tolerance;
    return out;
}

// --- Optical (Gordon) metric --------------------------------------------------
//
// For an isotropic medium of refractive index n at rest in flat spacetime, ray
// propagation follows the null geodesics of the Gordon metric
//
//     gbar_{mu nu} = g_{mu nu} + (1 - 1/n^2) u_mu u_nu
//
// with u the medium's 4-velocity. Reference: W. Gordon, "Zur
// Lichtfortpflanzung nach der Relativitatstheorie", Ann. Phys. 377, 421 (1923).
// (Bibliographic details recorded from memory; verify before citing.)
//
// For a static medium in Minkowski, u_mu = (-1,0,0,0) and
//
//     gbar = diag(-1/n^2, 1, 1, 1)
//
// whose null condition gives coordinate speed 1/n.
//
// IMPORTANT DISTINCTION. The Gordon metric describes RAYS in a medium that need
// not be impedance matched. It is NOT the same statement as the Plebanski map,
// which requires eps = mu and reproduces the FULL Maxwell system. A plain
// dielectric (eps = n^2, mu = 1) has a Gordon metric but is not
// metric-representable in the Plebanski sense. Conflating the two is the most
// common error in this area; see docs/effective_metric_theory.md section 4.
template <typename Scalar>
Mat4<Scalar> gordon_metric_static_isotropic(const Scalar& n) {
    Mat4<Scalar> g;
    const Scalar one = scalar_cast<Scalar>(1.0);
    g(0, 0) = -one / (n * n);
    g(1, 1) = one;
    g(2, 2) = one;
    g(3, 3) = one;
    return g;
}

}  // namespace program_b
