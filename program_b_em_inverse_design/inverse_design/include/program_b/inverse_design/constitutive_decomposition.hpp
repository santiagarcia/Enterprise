#pragma once

#include "program_b/inverse_design/maxwell.hpp"
#include "program_b/inverse_design/polynomial_roots.hpp"

// ---------------------------------------------------------------------------
// The premetric decomposition of a constitutive tensor, reported piece by piece.
//
// WHAT THIS REPLACES. B4 carried a single boolean called `metric_induced`,
// inferred from whether the Fresnel quartic factorized into a doubled quadric.
// That is not a metricity test. It conflates several independent questions, and
// the direction it appears to establish is the one it cannot:
//
//   * a metric-induced Z DOES give a repeated quadratic cone -- provable;
//   * a Z with a repeated quadratic cone need NOT be metric-induced -- because
//     the axion piece is invisible to the Fresnel polynomial entirely, and
//     because the cone only ever determines a conformal class.
//
// So light-cone factorization is a NECESSARY condition read off a cubic
// invariant, not a sufficient one, and the flags below are named for what they
// actually measure.
//
// THE PIECES. Following the premetric decomposition (Hehl & Obukhov), a tensor
// with the antisymmetries of Z splits into a principal part, a skewon part and
// an axion part. For the reciprocal tensors this project works with, pair
// symmetry Z^{abcd} = Z^{cdab} kills the skewon piece, leaving
//
//     21 = 20 (principal) + 1 (axion)
//
// and the axion part is
//
//     Z_axion^{abcd} = alpha eps^{abcd}
//
// with eps^{abcd} the Levi-Civita TENSOR of g_lab, in this project's own
// convention eps^{0123} = -1/sqrt(-det g_lab). It is totally antisymmetric,
// hence antisymmetric in each index pair and symmetric under pair exchange:
// it is a legitimate reciprocal constitutive tensor.
//
// WHY THE AXION IS INVISIBLE TO PROPAGATION. The Fresnel construction goes
// through
//
//     W^{bd} = Z^{abcd} k_a k_c
//
// and for the axion piece W_axion^{bd} = alpha eps^{abcd} k_a k_c = 0
// IDENTICALLY, because a totally antisymmetric symbol contracted with the
// symmetric k_a k_c vanishes. Not approximately, not up to a conformal factor:
// the axion contributes NOTHING to W, so nothing to its adjugate, so nothing to
// the light cone. That is the cleanest possible demonstration that
//
//     same Fresnel cone  =/=>  identical constitutive tensor.
//
// The axion does change Z, and it does change the observer blocks -- it appears
// as the TRACE part of a symmetric zeta, since
// zeta^{bd} = 2 V_a V_c (Z*)^{abcd} = 2 alpha h^{bd} for the axion piece alone.
// ---------------------------------------------------------------------------

namespace program_b {
namespace inverse_design {

// The axion tensor alpha eps^{abcd}, in this project's Levi-Civita convention.
template <typename Scalar>
Rank4<Scalar> axion_tensor(double alpha, const LaboratoryMetric<Scalar>& g_lab) {
    const Rank4<Scalar> e = levi_civita_upper(g_lab.value);
    Rank4<Scalar> Z;
    const Scalar a = scalar_cast<Scalar>(alpha);
    for (int i = 0; i < 256; ++i) Z.t[i] = a * e.t[i];
    return Z;
}

// The axion coefficient of a given Z: alpha = Z^{abcd} eps_{abcd} / (eps.eps).
// A projection, so it is exact for a pure axion tensor and zero for a pure
// principal one.
template <typename Scalar>
double axion_coefficient(const Rank4<Scalar>& Z,
                         const LaboratoryMetric<Scalar>& g_lab) {
    const Rank4<Scalar> eu = levi_civita_upper(g_lab.value);
    const Rank4<Scalar> el = levi_civita_lower(g_lab.value);
    double num = 0.0, den = 0.0;
    for (int i = 0; i < 256; ++i) {
        num += real_part(Z.t[i]) * real_part(el.t[i]);
        den += real_part(eu.t[i]) * real_part(el.t[i]);
    }
    return (std::fabs(den) > 0.0) ? num / den : 0.0;
}

// The skewon part is the piece that breaks pair symmetry:
//     Z_skewon^{abcd} = (Z^{abcd} - Z^{cdab}) / 2
// Its magnitude is exactly the reciprocity residual up to a factor of two, and
// it is reported as its own number rather than folded into a verdict.
template <typename Scalar>
double skewon_magnitude(const Rank4<Scalar>& Z) {
    double worst = 0.0;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d)
                    worst = std::fmax(worst,
                                      0.5 * std::fabs(real_part(Z(a, b, c, d)) -
                                                      real_part(Z(c, d, a, b))));
    return worst;
}

// --- The report ---------------------------------------------------------------
//
// Deliberately NOT a single `metric_induced` boolean. Each name says what it
// measures and nothing more.

struct ConstitutiveClassification {
    // Structural symmetries.
    double pair_antisymmetry_residual = 0.0;   // Z^{(ab)cd}, Z^{ab(cd)}
    double reciprocity_residual = 0.0;         // Z^{abcd} - Z^{cdab}
    bool pair_antisymmetric = false;
    bool reciprocal = false;

    // Premetric pieces.
    double axion_coefficient = 0.0;
    double skewon_magnitude = 0.0;
    double principal_magnitude = 0.0;
    bool contains_axion = false;
    bool contains_skewon = false;

    // What the Fresnel polynomial actually says.
    bool fresnel_single_cone = false;          // quartic is a doubled quadric
    double fresnel_double_root_residual = 0.0;
    double fresnel_branch_split = 0.0;

    // The most that can be concluded from the cone alone.
    bool principal_part_metric_compatible = false;

    std::string interpretation;
};

// Classify. The wording of `interpretation` is the corrected claim.
template <typename Scalar>
ConstitutiveClassification classify_constitutive(
    const Rank4<Scalar>& Z, const LaboratoryMetric<Scalar>& g_lab,
    double omega, const double direction[3]) {
    ConstitutiveClassification c;

    // --- symmetries ---------------------------------------------------------
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int d = 0; d < 4; ++d)
                for (int e = 0; e < 4; ++e) {
                    c.pair_antisymmetry_residual = std::fmax(
                        c.pair_antisymmetry_residual,
                        std::fabs(real_part(Z(a, b, d, e)) + real_part(Z(b, a, d, e))));
                    c.pair_antisymmetry_residual = std::fmax(
                        c.pair_antisymmetry_residual,
                        std::fabs(real_part(Z(a, b, d, e)) + real_part(Z(a, b, e, d))));
                }
    c.reciprocity_residual = pair_symmetry_residual(Z);

    double scale = 0.0;
    for (int i = 0; i < 256; ++i)
        scale = std::fmax(scale, std::fabs(real_part(Z.t[i])));
    const double rel = std::fmax(scale, 1e-300);
    c.pair_antisymmetric = c.pair_antisymmetry_residual <= 1e-12 * rel;
    c.reciprocal = c.reciprocity_residual <= 1e-12 * rel;

    // --- premetric pieces ---------------------------------------------------
    c.axion_coefficient = axion_coefficient(Z, g_lab);
    c.skewon_magnitude = skewon_magnitude(Z);

    const Rank4<Scalar> ax = axion_tensor<Scalar>(c.axion_coefficient, g_lab);
    for (int i = 0; i < 256; ++i)
        c.principal_magnitude = std::fmax(
            c.principal_magnitude,
            std::fabs(real_part(Z.t[i]) - real_part(ax.t[i])));

    c.contains_axion = std::fabs(c.axion_coefficient) > 1e-12 * rel;
    c.contains_skewon = c.skewon_magnitude > 1e-12 * rel;

    // --- what the Fresnel polynomial says -----------------------------------
    double coeff[5];
    fresnel_quartic_coefficients(Z, omega, direction, coeff);
    const RootSet rs = quartic_roots(coeff);
    c.fresnel_double_root_residual = rs.worst_residual;

    // A doubled quadric means every root has multiplicity two, i.e. exactly two
    // distinct roots for a quartic.
    int max_mult = 1;
    for (const PolynomialRoot& r : rs.roots)
        max_mult = std::max(max_mult, r.multiplicity);
    c.fresnel_single_cone = (rs.distinct_count == 2) && (max_mult == 2);

    double lo = 1e300, hi = -1e300;
    for (const PolynomialRoot& r : rs.real_roots)
        if (r.real > 0.0) { lo = std::fmin(lo, r.real); hi = std::fmax(hi, r.real); }
    c.fresnel_branch_split =
        (hi > 0.0 && lo < 1e300) ? (hi - lo) / std::fmax(hi, 1e-300) : 0.0;

    // The most the cone alone licenses.
    c.principal_part_metric_compatible = c.fresnel_single_cone;

    std::string m;
    m += c.reciprocal ? "reciprocal; " : "NOT reciprocal; ";
    m += c.contains_skewon ? "contains skewon; " : "no skewon; ";
    m += c.contains_axion ? "contains axion; " : "no axion; ";
    m += c.fresnel_single_cone ? "single Fresnel cone (principal part "
                                "metric-COMPATIBLE, which does NOT establish "
                                "metricity: the axion piece is invisible to the "
                                "cone and the cone fixes only a conformal class)"
                              : "birefringent: two distinct Fresnel cones, so the "
                                "response is NOT metric-induced";
    c.interpretation = m;
    return c;
}

}  // namespace inverse_design
}  // namespace program_b
