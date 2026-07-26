#pragma once

#include <string>

#include "program_b/inverse_design/design_space.hpp"

// ---------------------------------------------------------------------------
// B2 targets, and their classification.
//
// Two entry points: an effective metric, mapped through the verified covariant
// correspondence, or a constitutive tensor supplied directly. The second is what
// makes manufactured tests and non-metric targets possible.
//
// Every target carries a classification rather than an assumption. A target that
// is not reciprocal is REPORTED, not silently fitted: the B2 design space is
// reciprocal by construction, so the achievable residual against a
// skewon-carrying target is bounded below by the target's own skewon part, and a
// gate that did not say so would present that floor as an optimizer failure.
// ---------------------------------------------------------------------------

namespace program_b {
namespace inverse_design {

template <typename Scalar>
struct ConstitutiveTarget {
    Rank4<Scalar> Z;
    LaboratoryMetric<Scalar> g_lab;
    Scalar V[4];

    // Classification. Each is measured, not assumed.
    bool metric_induced = false;
    bool reciprocal = false;              // Z^{abcd} = Z^{cdab}
    bool pair_antisymmetric = false;      // Z^{(ab)cd} = Z^{ab(cd)} = 0
    bool finite = false;
    bool observer_decomposable = false;
    bool inside_real_lossless_model = false;

    double reciprocity_residual = 0.0;
    double antisymmetry_residual = 0.0;
    double max_component = 0.0;

    // The target's observer blocks, computed ONCE. They depend only on Z, g_lab
    // and V, all of which are fixed for the life of the target, so recomputing
    // them inside the optimizer's inner loop would be pure waste -- and
    // decompose_A7 is by far the most expensive thing in an evaluation.
    ObserverDecomposition4<Scalar> blocks;

    std::string origin;
    std::string message;

    bool usable() const {
        return finite && pair_antisymmetric && observer_decomposable;
    }
};

// max |Z^{(ab)cd}| and |Z^{ab(cd)}|: the symmetric part of either index pair,
// which must vanish by Schuster & Visser (12).
template <typename Scalar>
double pair_antisymmetry_residual(const Rank4<Scalar>& Z) {
    double worst = 0.0;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d) {
                    worst = std::fmax(worst, std::fabs(real_part(Z(a, b, c, d)) +
                                                       real_part(Z(b, a, c, d))));
                    worst = std::fmax(worst, std::fabs(real_part(Z(a, b, c, d)) +
                                                       real_part(Z(a, b, d, c))));
                }
    return worst;
}

// Shared classification, applied whichever entry point was used.
template <typename Scalar>
void classify_target(ConstitutiveTarget<Scalar>& t,
                     const ObserverFrame<Scalar>& frame) {
    t.reciprocity_residual = pair_symmetry_residual(t.Z);
    t.antisymmetry_residual = pair_antisymmetry_residual(t.Z);

    t.finite = true;
    t.max_component = 0.0;
    for (int i = 0; i < 256; ++i) {
        const double z = real_part(t.Z.t[i]);
        if (!std::isfinite(z)) t.finite = false;
        t.max_component = std::fmax(t.max_component, std::fabs(z));
    }

    const double scale = std::fmax(t.max_component, 1e-300);
    t.reciprocal = t.reciprocity_residual <= 1e-12 * scale;
    t.pair_antisymmetric = t.antisymmetry_residual <= 1e-12 * scale;
    t.observer_decomposable = frame.ok && t.finite;
    if (t.observer_decomposable)
        t.blocks = decompose_A7(t.Z, t.g_lab.value, frame.V);

    // The supported model: real, lossless, non-dispersive, local, reciprocal.
    // Nothing here can represent a complex or frequency-dependent response, so
    // this flag records the scope rather than testing for it.
    t.inside_real_lossless_model = t.finite && t.pair_antisymmetric;

    std::string m;
    m += t.reciprocal ? "reciprocal; " : "NOT reciprocal (skewon part present); ";
    m += t.pair_antisymmetric ? "pair-antisymmetric; "
                              : "NOT pair-antisymmetric; ";
    m += t.finite ? "finite; " : "NOT finite; ";
    m += t.observer_decomposable ? "observer-decomposable; "
                                 : "NOT observer-decomposable; ";
    m += t.metric_induced ? "metric-induced." : "not metric-induced.";
    t.message = m;
}

// Entry point 1: an effective metric, through Schuster & Visser eq (6).
template <typename Scalar>
ConstitutiveTarget<Scalar> target_from_effective_metric(
    const EffectiveMetric<Scalar>& g_target, const LaboratoryMetric<Scalar>& g_lab,
    const ObserverFrame<Scalar>& frame, const std::string& origin) {
    ConstitutiveTarget<Scalar> t;
    t.g_lab = g_lab;
    for (int a = 0; a < 4; ++a) t.V[a] = frame.V[a];
    t.origin = origin;

    const CovariantConstitutive<Scalar> cc =
        covariant_constitutive(g_target, g_lab);
    if (!cc.ok) {
        t.message = "no covariant constitutive tensor: " + cc.message;
        return t;
    }
    t.Z = cc.Z;
    t.metric_induced = true;
    classify_target(t, frame);
    return t;
}

// Entry point 2: a constitutive tensor supplied directly.
template <typename Scalar>
ConstitutiveTarget<Scalar> target_from_tensor(const Rank4<Scalar>& Z,
                                             const LaboratoryMetric<Scalar>& g_lab,
                                             const ObserverFrame<Scalar>& frame,
                                             const std::string& origin) {
    ConstitutiveTarget<Scalar> t;
    t.Z = Z;
    t.g_lab = g_lab;
    for (int a = 0; a < 4; ++a) t.V[a] = frame.V[a];
    t.origin = origin;
    t.metric_induced = false;
    classify_target(t, frame);
    return t;
}

// A target built from a known design, which is what every exact-recovery problem
// needs: the answer is known by construction and the optimizer has to find it.
template <typename Scalar>
ConstitutiveTarget<Scalar> target_from_design(
    const ConstitutiveDesignVariables<Scalar>& v,
    const LaboratoryMetric<Scalar>& g_lab, const ObserverFrame<Scalar>& frame,
    const std::string& origin) {
    return target_from_tensor(reconstruct_design_tensor(v, g_lab, frame), g_lab,
                              frame, origin);
}

}  // namespace inverse_design
}  // namespace program_b
