#pragma once

#include <cmath>
#include <string>

#include "warp_hypad/core/mat4.hpp"
#include "warp_hypad/core/scalar_traits.hpp"

// ---------------------------------------------------------------------------
// Program B, Milestone B1: linear electromagnetic constitutive response.
//
// The bianisotropic constitutive relations, in the 3+1 form this module uses:
//
//     D^i = eps^{ij} E_j + xi^{ij}  H_j
//     B^i = zeta^{ij} E_j + mu^{ij} H_j
//
// UNITS: Heaviside-Lorentz with c = 1, so vacuum is eps = mu = identity and
// xi = zeta = 0. This matches the geometrized units used by Program A and is
// the convention in which the Plebanski map takes its simplest form. SI
// conversion is a single overall factor per block and is NOT applied anywhere
// in this module.
//
// SIGN CONVENTION for the magnetoelectric blocks: see metric_to_medium.hpp.
// Sources differ, so this project fixes its own and verifies it rather than
// inheriting one.
//
// Everything is templated on Scalar so that OTI parameter sensitivities flow
// through the constitutive mapping exactly as they do through the metric.
// ---------------------------------------------------------------------------

namespace program_b {

using warp_hypad::Mat4;
using warp_hypad::real_part;
using warp_hypad::scalar_cast;

// Fixed-size 3x3, row-major: m(row, col) at offset 3*row + col.
template <typename Scalar>
struct Mat3 {
    Scalar m[9];

    Mat3() {
        const Scalar z = scalar_cast<Scalar>(0.0);
        for (int i = 0; i < 9; ++i) m[i] = z;
    }

    Scalar& operator()(int r, int c) { return m[3 * r + c]; }
    const Scalar& operator()(int r, int c) const { return m[3 * r + c]; }

    static Mat3 zero() { return Mat3(); }

    static Mat3 identity() {
        Mat3 r;
        const Scalar one = scalar_cast<Scalar>(1.0);
        for (int i = 0; i < 3; ++i) r(i, i) = one;
        return r;
    }

    static Mat3 diagonal(const Scalar& a, const Scalar& b, const Scalar& c) {
        Mat3 r;
        r(0, 0) = a; r(1, 1) = b; r(2, 2) = c;
        return r;
    }
};

template <typename Scalar>
Mat3<Scalar> mat3_mul(const Mat3<Scalar>& a, const Mat3<Scalar>& b) {
    Mat3<Scalar> c;
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k) {
            Scalar s = scalar_cast<Scalar>(0.0);
            for (int j = 0; j < 3; ++j) s = s + a(i, j) * b(j, k);
            c(i, k) = s;
        }
    return c;
}

template <typename Scalar>
Mat3<Scalar> mat3_scale(const Mat3<Scalar>& a, const Scalar& s) {
    Mat3<Scalar> r;
    for (int i = 0; i < 9; ++i) r.m[i] = a.m[i] * s;
    return r;
}

template <typename Scalar>
Mat3<Scalar> mat3_sub(const Mat3<Scalar>& a, const Mat3<Scalar>& b) {
    Mat3<Scalar> r;
    for (int i = 0; i < 9; ++i) r.m[i] = a.m[i] - b.m[i];
    return r;
}

template <typename Scalar>
Scalar mat3_det(const Mat3<Scalar>& a) {
    return a(0, 0) * (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)) -
           a(0, 1) * (a(1, 0) * a(2, 2) - a(1, 2) * a(2, 0)) +
           a(0, 2) * (a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0));
}

// Adjugate-based inverse, branch-free for the same reason as Mat4: pivoting
// would branch on a Scalar and break hypercomplex AD.
template <typename Scalar>
struct Mat3Inverse {
    Mat3<Scalar> inverse;
    Scalar determinant;
    double det_real = 0.0;
    bool ok = false;
};

template <typename Scalar>
Mat3Inverse<Scalar> mat3_inverse(const Mat3<Scalar>& a) {
    Mat3Inverse<Scalar> out;
    out.determinant = mat3_det(a);
    out.det_real = real_part(out.determinant);

    double scale = 0.0;
    for (int i = 0; i < 9; ++i) scale = std::fmax(scale, std::fabs(real_part(a.m[i])));
    const double cube = scale * scale * scale;
    if (!std::isfinite(out.det_real) ||
        std::fabs(out.det_real) <= 1e-13 * std::fmax(cube, 1e-300)) {
        return out;  // ok stays false
    }

    const Scalar inv = scalar_cast<Scalar>(1.0) / out.determinant;
    Mat3<Scalar>& b = out.inverse;
    b(0, 0) = (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)) * inv;
    b(0, 1) = (a(0, 2) * a(2, 1) - a(0, 1) * a(2, 2)) * inv;
    b(0, 2) = (a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1)) * inv;
    b(1, 0) = (a(1, 2) * a(2, 0) - a(1, 0) * a(2, 2)) * inv;
    b(1, 1) = (a(0, 0) * a(2, 2) - a(0, 2) * a(2, 0)) * inv;
    b(1, 2) = (a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2)) * inv;
    b(2, 0) = (a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0)) * inv;
    b(2, 1) = (a(0, 1) * a(2, 0) - a(0, 0) * a(2, 1)) * inv;
    b(2, 2) = (a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0)) * inv;
    out.ok = true;
    return out;
}

// Skew matrix with cross_matrix(w) * H == w x H, i.e. M_{ik} = eps_{ijk} w_j.
template <typename Scalar>
Mat3<Scalar> cross_matrix(const Scalar w[3]) {
    Mat3<Scalar> m;
    const Scalar z = scalar_cast<Scalar>(0.0);
    m(0, 0) = z;      m(0, 1) = -w[2];  m(0, 2) = w[1];
    m(1, 0) = w[2];   m(1, 1) = z;      m(1, 2) = -w[0];
    m(2, 0) = -w[1];  m(2, 1) = w[0];   m(2, 2) = z;
    return m;
}

template <typename Scalar>
double mat3_max_abs_diff(const Mat3<Scalar>& a, const Mat3<Scalar>& b) {
    double worst = 0.0;
    for (int i = 0; i < 9; ++i)
        worst = std::fmax(worst, std::fabs(real_part(a.m[i]) - real_part(b.m[i])));
    return worst;
}

template <typename Scalar>
double mat3_max_abs(const Mat3<Scalar>& a) {
    double worst = 0.0;
    for (int i = 0; i < 9; ++i)
        worst = std::fmax(worst, std::fabs(real_part(a.m[i])));
    return worst;
}

template <typename Scalar>
double mat3_symmetry_residual(const Mat3<Scalar>& a) {
    double worst = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j)
            worst = std::fmax(worst,
                              std::fabs(real_part(a(i, j)) - real_part(a(j, i))));
    return worst;
}

// ---------------------------------------------------------------------------
// The full linear response.
// ---------------------------------------------------------------------------

template <typename Scalar>
struct ConstitutiveResponse {
    Mat3<Scalar> epsilon;  // D = eps E + xi  H
    Mat3<Scalar> mu;       // B = zeta E + mu H
    Mat3<Scalar> xi;
    Mat3<Scalar> zeta;

    static ConstitutiveResponse vacuum() {
        ConstitutiveResponse r;
        r.epsilon = Mat3<Scalar>::identity();
        r.mu = Mat3<Scalar>::identity();
        return r;  // xi, zeta zero
    }
};

// --- Structural classification ---------------------------------------------
//
// Which of these hold decides whether a single pseudo-Riemannian metric can
// reproduce the response at all. See docs/effective_metric_theory.md.

struct ResponseClassification {
    // eps == mu (in these units). Plebanski media are impedance matched; a
    // medium that is not cannot be a metric in the full-Maxwell sense.
    bool impedance_matched = false;
    double impedance_mismatch = 0.0;

    // xi == -zeta, both antisymmetric. A metric produces exactly this form.
    // Departures correspond to the skewon and axion parts of the general
    // constitutive tensor, neither of which a metric can generate.
    bool magnetoelectric_is_metric_form = false;
    double magnetoelectric_residual = 0.0;

    bool epsilon_symmetric = false;
    bool mu_symmetric = false;
    double symmetry_residual = 0.0;

    // Necessary conditions above, all satisfied.
    bool metric_representable = false;

    // Human-readable reason when metric_representable is false.
    std::string obstruction;
};

template <typename Scalar>
ResponseClassification classify(const ConstitutiveResponse<Scalar>& r,
                                double tolerance = 1e-10) {
    ResponseClassification c;

    c.impedance_mismatch = mat3_max_abs_diff(r.epsilon, r.mu);
    c.impedance_matched = c.impedance_mismatch <= tolerance;

    // xi + zeta must vanish, and xi must be antisymmetric.
    Mat3<Scalar> sum;
    for (int i = 0; i < 9; ++i) sum.m[i] = r.xi.m[i] + r.zeta.m[i];
    double anti = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            anti = std::fmax(anti, std::fabs(real_part(r.xi(i, j)) +
                                             real_part(r.xi(j, i))));
    c.magnetoelectric_residual = std::fmax(mat3_max_abs(sum), anti);
    c.magnetoelectric_is_metric_form =
        c.magnetoelectric_residual <= tolerance;

    const double es = mat3_symmetry_residual(r.epsilon);
    const double ms = mat3_symmetry_residual(r.mu);
    c.symmetry_residual = std::fmax(es, ms);
    c.epsilon_symmetric = es <= tolerance;
    c.mu_symmetric = ms <= tolerance;

    c.metric_representable = c.impedance_matched &&
                             c.magnetoelectric_is_metric_form &&
                             c.epsilon_symmetric && c.mu_symmetric;

    if (!c.metric_representable) {
        if (!c.impedance_matched) {
            c.obstruction =
                "eps != mu: the medium is not impedance matched, so it is "
                "birefringent in the general case and no single "
                "pseudo-Riemannian metric reproduces full Maxwell propagation "
                "in it. An optical (Gordon) metric may still describe rays.";
        } else if (!c.magnetoelectric_is_metric_form) {
            c.obstruction =
                "the magnetoelectric blocks are not of the form xi = -zeta "
                "antisymmetric: the response carries a skewon and/or axion "
                "part, neither of which a metric can generate.";
        } else {
            c.obstruction =
                "eps or mu is not symmetric, which a metric-induced response "
                "always is.";
        }
    }
    return c;
}

}  // namespace program_b
