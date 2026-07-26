#pragma once

// ---------------------------------------------------------------------------
// The HyPAD-seeded gradient. Separate header because it is the ONLY part of B2
// that names the sensitivity backend, and it exists only when the backend does.
//
// The objective itself is NOT reimplemented here. evaluate_design() is
// instantiated at Scalar = OTIScalar and called unchanged; this file only seeds
// the active parameters and reads the channels back.
// ---------------------------------------------------------------------------

#ifdef WARP_HYPAD_HAVE_OTILIB

#include "program_b/inverse_design/design_gradient.hpp"
#include "warp_hypad/sensitivity/otilib_traits.hpp"

namespace program_b {
namespace inverse_design {

namespace sens = warp_hypad::sensitivity;
using sens::OTIScalar;

// How many parameters one evaluation can differentiate. Three, for onumm3n2.
inline constexpr int kGradientBatchSize = sens::kOtiDirections;

// One batch: seed the listed parameters into distinct OTI directions, evaluate,
// read the channels back.
//
// `active` must hold at most kGradientBatchSize indices. Passing more is a
// programming error and is reported rather than silently truncated.
inline GradientBatchResult evaluate_gradient_batch(
    const ConstitutiveDesignVariables<double>& variables,
    const ConstitutiveTarget<double>& target,
    const DesignConstraints<double>& constraints,
    const ObjectiveWeights<double>& weights, const int* active, int count) {
    GradientBatchResult r;
    if (count < 0 || count > kGradientBatchSize) {
        r.message = "batch larger than the backend's derivative directions (" +
                    std::to_string(kGradientBatchSize) + ")";
        return r;
    }

    // Promote the design, target and constraints to the AD scalar. The frame is
    // rebuilt at the AD scalar too, but it depends only on g_lab and the
    // observer, so it carries no channel -- which is the point.
    ConstitutiveDesignVariables<OTIScalar> v;
    v.family = variables.family;
    v.positive_definite_floor = variables.positive_definite_floor;
    for (int i = 0; i < kMaxDesignParameters; ++i)
        v.q[i] = OTIScalar(variables.q[i]);
    for (int b = 0; b < count; ++b) {
        r.active[b] = active[b];
        v.q[active[b]] = sens::seed_parameter(variables.q[active[b]], b);
    }
    r.count = count;

    Mat4<OTIScalar> lab_o;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            lab_o(a, b) = OTIScalar(target.g_lab.value(a, b));
    const LaboratoryMetric<OTIScalar> lab(lab_o);

    OTIScalar V_o[4];
    for (int a = 0; a < 4; ++a) V_o[a] = OTIScalar(target.V[a]);
    const Observer<OTIScalar> obs =
        make_laboratory_observer(lab, V_o, "laboratory");
    const ObserverFrame<OTIScalar> frame = build_observer_frame(lab, obs);
    if (!frame.ok) {
        r.message = "could not build the observer frame at the AD scalar: " +
                    frame.message;
        return r;
    }

    ConstitutiveTarget<OTIScalar> t;
    for (int i = 0; i < 256; ++i) t.Z.t[i] = OTIScalar(target.Z.t[i]);
    t.g_lab = lab;
    for (int a = 0; a < 4; ++a) t.V[a] = V_o[a];
    t.origin = target.origin;

    DesignConstraints<OTIScalar> c;
    c.active = constraints.active;
    c.max_epsilon_entry = constraints.max_epsilon_entry;
    c.max_mu_inverse_entry = constraints.max_mu_inverse_entry;
    c.max_zeta_entry = constraints.max_zeta_entry;
    c.max_epsilon_frobenius = constraints.max_epsilon_frobenius;
    c.max_mu_inverse_frobenius = constraints.max_mu_inverse_frobenius;
    c.max_zeta_frobenius = constraints.max_zeta_frobenius;
    c.enforce_determinant_floor = constraints.enforce_determinant_floor;
    c.min_epsilon_determinant = constraints.min_epsilon_determinant;
    c.min_mu_inverse_determinant = constraints.min_mu_inverse_determinant;
    c.penalty_weight = constraints.penalty_weight;
    c.hinge_smoothing = constraints.hinge_smoothing;

    ObjectiveWeights<OTIScalar> w;
    w.primary = weights.primary;
    w.norm = weights.norm;
    w.w_epsilon = weights.w_epsilon;
    w.w_mu_inverse = weights.w_mu_inverse;
    w.w_zeta = weights.w_zeta;
    w.delta_Z = weights.delta_Z;
    w.delta_block = weights.delta_block;

    // The SAME production objective, at a different scalar.
    const DesignEvaluation<OTIScalar> e =
        evaluate_design(v, t, c, frame, w);

    r.objective = real_part(e.J_total);
    for (int b = 0; b < count; ++b)
        r.partial[b] = sens::extract_gradient(e.J_total, b);
    r.ok = std::isfinite(r.objective);
    r.message = r.ok ? "seeded " + std::to_string(count) + " of " +
                           std::to_string(kGradientBatchSize) + " directions"
                     : "non-finite objective";
    return r;
}

// The full gradient, assembled across batches. Deterministic partition: indices
// in order, kGradientBatchSize at a time.
inline GradientResult oti_gradient(
    const ConstitutiveDesignVariables<double>& variables,
    const ConstitutiveTarget<double>& target,
    const DesignConstraints<double>& constraints,
    const ObjectiveWeights<double>& weights) {
    GradientResult r;
    r.size = variables.size();
    r.ok = true;

    for (int start = 0; start < r.size; start += kGradientBatchSize) {
        int active[kGradientBatchSize];
        int count = 0;
        for (int b = 0; b < kGradientBatchSize && start + b < r.size; ++b)
            active[count++] = start + b;

        const GradientBatchResult batch = evaluate_gradient_batch(
            variables, target, constraints, weights, active, count);
        r.batches.push_back(batch);
        ++r.evaluations;
        if (!batch.ok) {
            r.ok = false;
            r.message = batch.message;
            return r;
        }
        r.objective = batch.objective;
        for (int b = 0; b < count; ++b) r.g[active[b]] = batch.partial[b];
    }
    r.message = "assembled from " + std::to_string(r.evaluations) +
                " batches of at most " + std::to_string(kGradientBatchSize);
    return r;
}

// The same partial derivative, obtained with the parameter seeded ALONE. If a
// batch of three and a batch of one disagree, the directions are leaking into
// each other.
inline double oti_partial_alone(
    const ConstitutiveDesignVariables<double>& variables,
    const ConstitutiveTarget<double>& target,
    const DesignConstraints<double>& constraints,
    const ObjectiveWeights<double>& weights, int index) {
    const int active[1] = {index};
    const GradientBatchResult b = evaluate_gradient_batch(
        variables, target, constraints, weights, active, 1);
    return b.ok ? b.partial[0] : std::numeric_limits<double>::quiet_NaN();
}

}  // namespace inverse_design
}  // namespace program_b

#endif  // WARP_HYPAD_HAVE_OTILIB
