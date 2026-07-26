#pragma once

#include <limits>

#include "program_b/effective_metric/bel_extraction.hpp"
#include "program_b/effective_metric/covariant_constitutive.hpp"

// ---------------------------------------------------------------------------
// The five-way status report, on the PRODUCTION (A7) decomposition.
//
// This replaces legacy_unverified::analyse(), which routed through the guessed
// tetrad decomposition. The classification logic is unchanged -- it is the
// correction B1C exists to make, and it was never the part that was wrong --
// but the constitutive blocks it inspects now come from Schuster & Visser
// (A7a)-(A7d) via decompose_A7().
//
// TWO OBSERVERS, deliberately separate arguments:
//
//   lab_observer         a physical observer of g_lab. The constitutive blocks
//                        belong to this one. Whether they are finite and usable
//                        is a statement about the MEDIUM.
//
//   effective_observer   a diagnostic observer of g_eff. Whether it is timelike
//                        is a statement about the ANALOGUE GEOMETRY -- it is
//                        what fails on an ergosurface, and it says nothing
//                        about the medium or about other observers.
//
// Collapsing those two questions into one verdict is exactly the error B1 made.
// The types make the two arguments impossible to swap.
// ---------------------------------------------------------------------------

namespace program_b {

struct MetricMediumReport {
    bool metric_lorentzian_nondegenerate = false;
    bool covariant_tensor_exists = false;
    bool effective_observer_admissible = false;
    bool constitutive_finite = false;

    MetricMediumStatus status =
        MetricMediumStatus::MetricDegenerateOrNonLorentzian;

    // Magnitude diagnostics on the laboratory observer's blocks. Reported
    // rather than reduced to a single conditioning number: a block can scale
    // uniformly towards infinity with its condition number pinned at 1, which
    // is precisely why the old cond(eps) divergence test could not fire.
    double epsilon_max_component = 0.0;
    double mu_inverse_max_component = 0.0;
    double zeta_max_component = 0.0;

    std::string message;
};

template <typename Scalar>
MetricMediumReport analyse_metric_medium(
    const EffectiveMetric<Scalar>& g_eff, const LaboratoryMetric<Scalar>& g_lab,
    const Observer<Scalar>& lab_observer,
    const EffectiveGeometryObserver<Scalar>& effective_observer,
    double magnitude_limit = 1e12) {
    MetricMediumReport r;

    const CovariantConstitutive<Scalar> cc =
        covariant_constitutive(g_eff, g_lab);
    r.metric_lorentzian_nondegenerate = cc.ok;
    r.covariant_tensor_exists = cc.ok;
    if (!cc.ok) {
        r.status = MetricMediumStatus::MetricDegenerateOrNonLorentzian;
        r.message = cc.message;
        return r;
    }

    r.effective_observer_admissible = effective_observer.admissible;
    if (!effective_observer.admissible) {
        r.status = MetricMediumStatus::SelectedObserverDecompositionInvalid;
        r.message = effective_observer.reason;
        return r;
    }

    if (!lab_observer.admissible) {
        r.status = MetricMediumStatus::SelectedObserverDecompositionInvalid;
        r.message = "the LABORATORY observer is not timelike in g_lab: " +
                    lab_observer.reason;
        return r;
    }

    const ObserverDecomposition4<Scalar> d =
        decompose_A7(cc.Z, g_lab.value, lab_observer.V);

    bool finite = true;
    for (int i = 0; i < 16; ++i) {
        const double e = std::fabs(real_part(d.epsilon.m[i]));
        const double m = std::fabs(real_part(d.mu_inverse.m[i]));
        const double z = std::fabs(real_part(d.zeta.m[i]));
        const double za = std::fabs(real_part(d.zeta_adjoint.m[i]));
        if (!std::isfinite(e) || !std::isfinite(m) || !std::isfinite(z) ||
            !std::isfinite(za))
            finite = false;
        r.epsilon_max_component = std::fmax(r.epsilon_max_component, e);
        r.mu_inverse_max_component = std::fmax(r.mu_inverse_max_component, m);
        r.zeta_max_component = std::fmax(r.zeta_max_component, z);
    }

    const double biggest =
        std::fmax(r.epsilon_max_component,
                  std::fmax(r.mu_inverse_max_component, r.zeta_max_component));

    if (!finite || biggest > magnitude_limit) {
        r.status = MetricMediumStatus::ConstitutiveParametersSingular;
        r.message =
            "the covariant tensor is finite and both observers are admissible, "
            "but the laboratory observer's constitutive blocks have run away in "
            "magnitude. A different observer may still work.";
        return r;
    }

    r.constitutive_finite = true;
    r.status = MetricMediumStatus::CovariantMetricMediumValid;
    r.message =
        "covariant tensor valid, both observers admissible, laboratory blocks "
        "finite.";
    return r;
}

}  // namespace program_b
