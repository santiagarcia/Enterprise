#pragma once

#include <array>
#include <cmath>
#include <string>

#include "warp_hypad/core/mat4.hpp"
#include "warp_hypad/core/scalar_traits.hpp"

// ---------------------------------------------------------------------------
// O0B: typed ADM design fields, structurally valid metric assembly, and the
// analytic inverse.
//
// WHY THE PUBLIC INTERFACE IS NOT A Mat4.
//
// A bare 4x4 of numbers can hold anything: a metric with a positive g_00, a
// spatial block that is not positive definite, a signature that is not
// Lorentzian. If the design layer accepted one, then every downstream
// guarantee -- that the lapse is positive, that gamma is a metric, that the
// signature is (-,+,+,+) -- would depend on whoever filled the array in. So the
// only public route to a metric is through ADMDesignFields, which cannot
// represent an invalid one:
//
//     alpha = exp(a)          positive for every finite a, by construction
//     gamma = L L^T           positive definite for every finite L with
//                             positive diagonal, by construction
//
// and the assembled g is then Lorentzian for free. These are STRUCTURAL
// guarantees, not penalties an optimizer is asked to respect: there is no
// coefficient vector, valid or not, that produces a metric with the wrong
// signature, because the parameterisation has no such point.
//
// THE SIGN CONVENTION IS READ FROM THE REPOSITORY, NOT FROM MEMORY.
//
// metrics::AlcubierreMetric::evaluate writes, in the (-,+,+,+) signature,
//
//     g_00 = -1 + (v f)^2       g_01 = g_10 = -v f       g_ij = delta_ij
//
// The ADM form of the same line element is
//
//     g_00 = -alpha^2 + gamma_ij beta^i beta^j    g_0i = gamma_ij beta^j
//
// so with alpha = 1 and gamma = I the two agree if and only if
//
//     beta^x = -v f
//
// NEGATIVE. That sign is not a convention this file is free to choose; it is
// fixed by the existing production metric, and O0B's central claim -- that a
// zero correction reproduces the Alcubierre metric exactly -- is what tests it.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace design {

// Spacetime index order, recorded in every output that carries a derivative.
// 0 = t, 1 = x, 2 = y, 3 = z. Stated as data rather than as a comment, because
// a jet whose index order is only documented is a jet whose index order can be
// wrong in one file and right in another.
inline constexpr int kTimeIndex = 0;
inline constexpr int kIndexX = 1;
inline constexpr int kIndexY = 2;
inline constexpr int kIndexZ = 3;
inline const char* index_ordering() { return "0=t, 1=x, 2=y, 3=z"; }

// --- Symmetric 3x3 -----------------------------------------------------------

template <typename Scalar>
struct SymmetricMatrix3 {
    // Upper-triangular storage: (00, 01, 02, 11, 12, 22). Symmetry is a
    // property of the STORAGE, so a caller cannot make it asymmetric.
    std::array<Scalar, 6> m{};

    static int slot(int i, int j) {
        if (i > j) { const int t = i; i = j; j = t; }
        static const int table[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
        return table[i][j];
    }
    Scalar& at(int i, int j) { return m[static_cast<std::size_t>(slot(i, j))]; }
    const Scalar& at(int i, int j) const {
        return m[static_cast<std::size_t>(slot(i, j))];
    }

    static SymmetricMatrix3 identity() {
        SymmetricMatrix3 s;
        for (int i = 0; i < 6; ++i) s.m[static_cast<std::size_t>(i)] =
            scalar_cast<Scalar>(0.0);
        s.at(0, 0) = scalar_cast<Scalar>(1.0);
        s.at(1, 1) = scalar_cast<Scalar>(1.0);
        s.at(2, 2) = scalar_cast<Scalar>(1.0);
        return s;
    }
    static SymmetricMatrix3 zero() {
        SymmetricMatrix3 s;
        for (int i = 0; i < 6; ++i) s.m[static_cast<std::size_t>(i)] =
            scalar_cast<Scalar>(0.0);
        return s;
    }
};

// det of a symmetric 3x3.
template <typename Scalar>
Scalar determinant3(const SymmetricMatrix3<Scalar>& g) {
    const Scalar a = g.at(0, 0), b = g.at(0, 1), c = g.at(0, 2);
    const Scalar d = g.at(1, 1), e = g.at(1, 2), f = g.at(2, 2);
    return a * (d * f - e * e) - b * (b * f - e * c) + c * (b * e - d * c);
}

// The inverse of a symmetric positive-definite 3x3, by cofactors.
template <typename Scalar>
SymmetricMatrix3<Scalar> inverse3(const SymmetricMatrix3<Scalar>& g,
                                  bool& ok) {
    const Scalar det = determinant3(g);
    SymmetricMatrix3<Scalar> out = SymmetricMatrix3<Scalar>::zero();
    if (real_abs(det) < 1e-300) { ok = false; return out; }
    ok = true;
    const Scalar a = g.at(0, 0), b = g.at(0, 1), c = g.at(0, 2);
    const Scalar d = g.at(1, 1), e = g.at(1, 2), f = g.at(2, 2);
    const Scalar inv = scalar_cast<Scalar>(1.0) / det;
    out.at(0, 0) = (d * f - e * e) * inv;
    out.at(0, 1) = (c * e - b * f) * inv;
    out.at(0, 2) = (b * e - c * d) * inv;
    out.at(1, 1) = (a * f - c * c) * inv;
    out.at(1, 2) = (b * c - a * e) * inv;
    out.at(2, 2) = (a * d - b * b) * inv;
    return out;
}

// --- The design fields -------------------------------------------------------

template <typename Scalar>
struct ADMDesignFields {
    Scalar lapse = scalar_cast<Scalar>(1.0);           // alpha > 0
    std::array<Scalar, 3> shift_contravariant{};       // beta^i, INDEX UP
    SymmetricMatrix3<Scalar> spatial_metric =
        SymmetricMatrix3<Scalar>::identity();          // gamma_ij, INDICES DOWN

    static ADMDesignFields minkowski() {
        ADMDesignFields a;
        a.lapse = scalar_cast<Scalar>(1.0);
        for (int i = 0; i < 3; ++i)
            a.shift_contravariant[static_cast<std::size_t>(i)] =
                scalar_cast<Scalar>(0.0);
        a.spatial_metric = SymmetricMatrix3<Scalar>::identity();
        return a;
    }
};

// The two metric roles, distinguished by TYPE and not by a naming convention.
// Passing a contravariant metric where a covariant one is wanted is the kind of
// error that produces a plausible, wrong answer; here it does not compile.
template <typename Scalar>
struct CovariantMetric {
    Mat4<Scalar> m;
};
template <typename Scalar>
struct ContravariantMetric {
    Mat4<Scalar> m;
};

// --- Structural parameterisations --------------------------------------------

// alpha = exp(a). Positive for every finite a; the caller supplies a, never
// alpha, so a non-positive lapse has no representation.
template <typename Scalar>
Scalar lapse_from_log(const Scalar& a) {
    using std::exp;
    return exp(a);
}

// The largest |a| whose exponential is still comfortably finite. Rejected
// BEFORE a metric is produced, so an overflowing lapse never reaches a
// downstream consumer as an infinity.
inline constexpr double kMaxLogLapse = 300.0;

// gamma = L L^T with L lower triangular, exp on the diagonal. Positive definite
// for every finite (l1, l2, l3, l21, l31, l32), by construction: det gamma =
// exp(2(l1+l2+l3)) > 0 and every leading minor is likewise a product of
// exponentials.
template <typename Scalar>
struct CholeskyFactorL {
    Scalar l1 = scalar_cast<Scalar>(0.0);   // log of L(0,0)
    Scalar l2 = scalar_cast<Scalar>(0.0);   // log of L(1,1)
    Scalar l3 = scalar_cast<Scalar>(0.0);   // log of L(2,2)
    Scalar l21 = scalar_cast<Scalar>(0.0);
    Scalar l31 = scalar_cast<Scalar>(0.0);
    Scalar l32 = scalar_cast<Scalar>(0.0);

    static CholeskyFactorL identity() { return CholeskyFactorL(); }
};

template <typename Scalar>
SymmetricMatrix3<Scalar> spatial_metric_from_cholesky(
    const CholeskyFactorL<Scalar>& L) {
    using std::exp;
    const Scalar e1 = exp(L.l1), e2 = exp(L.l2), e3 = exp(L.l3);
    // L =  [ e1   0   0 ]
    //      [ l21  e2  0 ]
    //      [ l31 l32  e3 ]
    SymmetricMatrix3<Scalar> g;
    g.at(0, 0) = e1 * e1;
    g.at(0, 1) = L.l21 * e1;
    g.at(0, 2) = L.l31 * e1;
    g.at(1, 1) = L.l21 * L.l21 + e2 * e2;
    g.at(1, 2) = L.l31 * L.l21 + L.l32 * e2;
    g.at(2, 2) = L.l31 * L.l31 + L.l32 * L.l32 + e3 * e3;
    return g;
}

// d(gamma)/dp for one parameter direction, given dL/dp. The product rule on
// gamma = L L^T, written out so the F2/F3 spatial-metric families inherit exact
// derivatives rather than acquiring finite differences later. Exercised with
// nonzero synthetic L even though O0B-F1 uses L = I, so a dormant path is not
// silently wrong when it is switched on.
template <typename Scalar>
SymmetricMatrix3<Scalar> spatial_metric_derivative(
    const CholeskyFactorL<Scalar>& L, const CholeskyFactorL<Scalar>& dL) {
    using std::exp;
    const Scalar e1 = exp(L.l1), e2 = exp(L.l2), e3 = exp(L.l3);
    const Scalar de1 = e1 * dL.l1, de2 = e2 * dL.l2, de3 = e3 * dL.l3;
    SymmetricMatrix3<Scalar> g;
    g.at(0, 0) = scalar_cast<Scalar>(2.0) * e1 * de1;
    g.at(0, 1) = dL.l21 * e1 + L.l21 * de1;
    g.at(0, 2) = dL.l31 * e1 + L.l31 * de1;
    g.at(1, 1) = scalar_cast<Scalar>(2.0) * (L.l21 * dL.l21 + e2 * de2);
    g.at(1, 2) = dL.l31 * L.l21 + L.l31 * dL.l21 + dL.l32 * e2 + L.l32 * de2;
    g.at(2, 2) = scalar_cast<Scalar>(2.0) *
                 (L.l31 * dL.l31 + L.l32 * dL.l32 + e3 * de3);
    return g;
}

// --- Assembly ----------------------------------------------------------------

// g_00 = -alpha^2 + gamma_ij beta^i beta^j
// g_0i = gamma_ij beta^j
// g_ij = gamma_ij
template <typename Scalar>
CovariantMetric<Scalar> assemble_covariant(const ADMDesignFields<Scalar>& a) {
    CovariantMetric<Scalar> out;
    Mat4<Scalar>& g = out.m;

    // beta_i = gamma_ij beta^j, the shift with its index DOWN. Distinguished
    // from beta^i throughout: they differ whenever gamma is not the identity,
    // and conflating them is the classic ADM sign error.
    Scalar beta_down[3];
    for (int i = 0; i < 3; ++i) {
        Scalar acc = scalar_cast<Scalar>(0.0);
        for (int j = 0; j < 3; ++j)
            acc = acc + a.spatial_metric.at(i, j) *
                            a.shift_contravariant[static_cast<std::size_t>(j)];
        beta_down[i] = acc;
    }

    Scalar beta_sq = scalar_cast<Scalar>(0.0);
    for (int i = 0; i < 3; ++i)
        beta_sq = beta_sq +
                  beta_down[i] * a.shift_contravariant[static_cast<std::size_t>(i)];

    g(0, 0) = -(a.lapse * a.lapse) + beta_sq;
    for (int i = 0; i < 3; ++i) {
        g(0, i + 1) = beta_down[i];
        g(i + 1, 0) = beta_down[i];
    }
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) g(i + 1, j + 1) = a.spatial_metric.at(i, j);
    return out;
}

// g^00 = -1/alpha^2
// g^0i = beta^i/alpha^2
// g^ij = gamma^ij - beta^i beta^j / alpha^2
//
// Assembled ANALYTICALLY rather than by inverting the 4x4. The inverse of the
// assembled matrix is then an independent check with a known answer, which a
// numerical inversion of the same object could not provide.
template <typename Scalar>
ContravariantMetric<Scalar> assemble_contravariant(
    const ADMDesignFields<Scalar>& a, bool& ok) {
    ContravariantMetric<Scalar> out;
    Mat4<Scalar>& g = out.m;
    const Scalar a2 = a.lapse * a.lapse;
    if (real_abs(a2) < 1e-300) { ok = false; return out; }
    const SymmetricMatrix3<Scalar> gamma_inv = inverse3(a.spatial_metric, ok);
    if (!ok) return out;
    const Scalar inv_a2 = scalar_cast<Scalar>(1.0) / a2;

    g(0, 0) = -inv_a2;
    for (int i = 0; i < 3; ++i) {
        const Scalar v =
            a.shift_contravariant[static_cast<std::size_t>(i)] * inv_a2;
        g(0, i + 1) = v;
        g(i + 1, 0) = v;
    }
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            g(i + 1, j + 1) =
                gamma_inv.at(i, j) -
                a.shift_contravariant[static_cast<std::size_t>(i)] *
                    a.shift_contravariant[static_cast<std::size_t>(j)] * inv_a2;
    ok = true;
    return out;
}

// det g = -alpha^2 det gamma. An identity of the ADM decomposition, so it is an
// exact test rather than a tolerance one: for F1 (alpha = 1, gamma = I) it must
// give exactly -1 for EVERY finite shift, including where g_00 = 0.
template <typename Scalar>
Scalar determinant_from_adm(const ADMDesignFields<Scalar>& a) {
    return -(a.lapse * a.lapse) * determinant3(a.spatial_metric);
}

// --- Validity ----------------------------------------------------------------

enum class MetricRejection {
    None,
    NonFiniteCoefficient,
    LapseOverflow,
    NonPositiveLapse,
    SpatialMetricNotPositiveDefinite,
    SingularSpatialMetric,
    DimensionMismatch,
};

inline const char* to_string(MetricRejection r) {
    switch (r) {
        case MetricRejection::None: return "none";
        case MetricRejection::NonFiniteCoefficient:
            return "nonfinite_coefficient";
        case MetricRejection::LapseOverflow: return "lapse_overflow";
        case MetricRejection::NonPositiveLapse: return "non_positive_lapse";
        case MetricRejection::SpatialMetricNotPositiveDefinite:
            return "spatial_metric_not_positive_definite";
        case MetricRejection::SingularSpatialMetric:
            return "singular_spatial_metric";
        case MetricRejection::DimensionMismatch: return "dimension_mismatch";
    }
    return "unknown";
}

// Sylvester's criterion on the leading minors: positive definite iff all three
// are positive. Checked on the ASSEMBLED gamma rather than trusted from the
// parameterisation, so a future family that builds gamma some other way cannot
// bypass it.
template <typename Scalar>
bool spatial_metric_positive_definite(const SymmetricMatrix3<Scalar>& g) {
    const double m1 = real_part(g.at(0, 0));
    const double m2 = real_part(g.at(0, 0) * g.at(1, 1) - g.at(0, 1) * g.at(0, 1));
    const double m3 = real_part(determinant3(g));
    return m1 > 0.0 && m2 > 0.0 && m3 > 0.0;
}

// The signature of the assembled 4x4, by the sign of the eigenvalues of the
// symmetric matrix. Reported as counts so "(-,+,+,+)" is a measurement.
struct SignatureCount {
    int negative = 0;
    int positive = 0;
    int zero = 0;
    bool lorentzian() const {
        return negative == 1 && positive == 3 && zero == 0;
    }
};

// --- The metric jet ----------------------------------------------------------

template <typename Scalar>
struct MetricJet {
    CovariantMetric<Scalar> g_cov;
    ContravariantMetric<Scalar> g_contra;

    // Index order is kTimeIndex..kIndexZ throughout, recorded by
    // index_ordering() in every emitted record.
    std::array<CovariantMetric<Scalar>, 4> dg{};
    std::array<std::array<CovariantMetric<Scalar>, 4>, 4> d2g{};
    std::array<std::array<std::array<CovariantMetric<Scalar>, 4>, 4>, 4> d3g{};

    bool valid = false;
    MetricRejection rejection = MetricRejection::None;
    std::string index_order = index_ordering();
};

// max |d_mu d_nu g - d_nu d_mu g| over the jet: mixed partials of a smooth
// field commute, so this is zero by construction if the chain rule was applied
// correctly and nonzero if a channel was assembled in the wrong order.
template <typename Scalar>
double second_derivative_symmetry_residual(const MetricJet<Scalar>& j) {
    double worst = 0.0;
    for (int m = 0; m < 4; ++m)
        for (int n = 0; n < 4; ++n)
            for (int A = 0; A < 4; ++A)
                for (int B = 0; B < 4; ++B)
                    worst = std::fmax(
                        worst, real_abs(j.d2g[static_cast<std::size_t>(m)]
                                            [static_cast<std::size_t>(n)].m(A, B) -
                                        j.d2g[static_cast<std::size_t>(n)]
                                            [static_cast<std::size_t>(m)].m(A, B)));
    return worst;
}

// The same for every permutation of the three third-derivative indices.
template <typename Scalar>
double third_derivative_symmetry_residual(const MetricJet<Scalar>& j) {
    const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                            {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    double worst = 0.0;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c) {
                const int idx[3] = {a, b, c};
                for (int A = 0; A < 4; ++A)
                    for (int B = 0; B < 4; ++B) {
                        const double ref =
                            real_part(j.d3g[static_cast<std::size_t>(a)]
                                          [static_cast<std::size_t>(b)]
                                          [static_cast<std::size_t>(c)].m(A, B));
                        for (int p = 0; p < 6; ++p) {
                            const double v = real_part(
                                j.d3g[static_cast<std::size_t>(idx[perm[p][0]])]
                                     [static_cast<std::size_t>(idx[perm[p][1]])]
                                     [static_cast<std::size_t>(idx[perm[p][2]])]
                                         .m(A, B));
                            worst = std::fmax(worst, std::fabs(v - ref));
                        }
                    }
            }
    return worst;
}

// max |g_{mu lambda} g^{lambda nu} - delta| for the assembled pair.
template <typename Scalar>
double inverse_residual(const CovariantMetric<Scalar>& cov,
                        const ContravariantMetric<Scalar>& con) {
    double worst = 0.0;
    for (int m = 0; m < 4; ++m)
        for (int n = 0; n < 4; ++n) {
            Scalar acc = scalar_cast<Scalar>(0.0);
            for (int l = 0; l < 4; ++l) acc = acc + cov.m(m, l) * con.m(l, n);
            worst = std::fmax(worst,
                              std::fabs(real_part(acc) - (m == n ? 1.0 : 0.0)));
        }
    return worst;
}

}  // namespace design
}  // namespace warp_hypad
