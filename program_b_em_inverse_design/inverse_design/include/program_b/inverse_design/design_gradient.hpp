#pragma once

#include <vector>

#include "program_b/inverse_design/design_objective.hpp"

// ---------------------------------------------------------------------------
// Gradients of the B2 objective with respect to the design vector.
//
// BATCHING IS NOT OPTIONAL. The available HyPAD scalar is onumm3n2: three
// derivative directions, order two. It cannot seed 21 variables at once, and an
// interface that pretended otherwise would be a lie about the backend. So the
// unit of work is a BATCH of at most kOtiDirections active parameters, and the
// full gradient is assembled across ceil(n/3) = 7 batches for the 21-parameter
// family.
//
// There is no separate OTI implementation of the objective. evaluate_design() is
// scalar-generic and the same function serves double, the test scalar and the
// HyPAD scalar; only the seeding differs.
//
// ONE EXACTNESS PROPERTY IS EXPLOITED. reconstruct_Z_A11 is homogeneous linear in
// the blocks and the blocks are linear in q for every family except
// PositiveDefinite, so J_Z is an exact QUADRATIC in q. Centred differences of a
// quadratic carry zero truncation error, so the finite-difference gradient is
// exact to roundoff and the comparison can be tight rather than step-limited.
// The gate asserts the quadraticity rather than assuming it.
// ---------------------------------------------------------------------------

namespace program_b {
namespace inverse_design {

struct GradientBatchResult {
    // Which parameters this batch differentiated, and their partials.
    int active[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
    double partial[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    int count = 0;
    double objective = 0.0;
    bool ok = false;
    std::string message;
};

struct GradientResult {
    double g[kMaxDesignParameters] = {};
    int size = 0;
    double objective = 0.0;
    int evaluations = 0;
    std::vector<GradientBatchResult> batches;
    bool ok = false;
    std::string message;

    double norm() const {
        double s = 0.0;
        for (int i = 0; i < size; ++i) s += g[i] * g[i];
        return std::sqrt(s);
    }
};

// --- Finite-difference gradient ----------------------------------------------
//
// Available in every build, with or without the HyPAD backend. This is what
// keeps the optimizer and every inverse-design demonstration runnable when
// OTIlib is disabled.

template <typename Scalar>
GradientResult finite_difference_gradient(
    const ConstitutiveDesignVariables<Scalar>& variables,
    const ConstitutiveTarget<Scalar>& target,
    const DesignConstraints<Scalar>& constraints,
    const ObserverFrame<Scalar>& frame, const ObjectiveWeights<Scalar>& weights,
    double step) {
    GradientResult r;
    r.size = variables.size();
    const DesignEvaluation<Scalar> base =
        evaluate_design(variables, target, constraints, frame, weights);
    r.objective = real_part(base.J_total);
    r.evaluations = 1;

    for (int k = 0; k < r.size; ++k) {
        ConstitutiveDesignVariables<Scalar> hi = variables, lo = variables;
        hi.q[k] = hi.q[k] + scalar_cast<Scalar>(step);
        lo.q[k] = lo.q[k] - scalar_cast<Scalar>(step);
        const double Jhi = real_part(
            evaluate_design(hi, target, constraints, frame, weights).J_total);
        const double Jlo = real_part(
            evaluate_design(lo, target, constraints, frame, weights).J_total);
        r.g[k] = (Jhi - Jlo) / (2.0 * step);
        r.evaluations += 2;
    }
    r.ok = std::isfinite(r.objective);
    r.message = r.ok ? "centred finite differences" : "non-finite objective";
    return r;
}

// --- Manufactured quadratic, with a hand-derived gradient --------------------
//
// J(q) = sum_k c_k (q_k - a_k)^2,   dJ/dq_k = 2 c_k (q_k - a_k)
//
// Deliberately independent of the constitutive objective: it tests the seeding
// and batch-assembly machinery in isolation, so a failure names the AD plumbing
// rather than the physics.

struct ManufacturedQuadratic {
    double c[kMaxDesignParameters];
    double a[kMaxDesignParameters];
    int size = 0;

    ManufacturedQuadratic() {
        for (int i = 0; i < kMaxDesignParameters; ++i) {
            c[i] = 0.0;
            a[i] = 0.0;
        }
    }
};

// Deterministic, reproducible coefficients: no random engine, no seed to lose.
inline ManufacturedQuadratic make_manufactured_quadratic(int n) {
    ManufacturedQuadratic m;
    m.size = n;
    for (int i = 0; i < n; ++i) {
        m.c[i] = 0.5 + 0.25 * static_cast<double>(i % 5);
        m.a[i] = -0.7 + 0.13 * static_cast<double>(i);
    }
    return m;
}

template <typename Scalar>
Scalar manufactured_objective(const ManufacturedQuadratic& m, const Scalar q[]) {
    Scalar acc = scalar_cast<Scalar>(0.0);
    for (int i = 0; i < m.size; ++i) {
        const Scalar d = q[i] - scalar_cast<Scalar>(m.a[i]);
        acc = acc + scalar_cast<Scalar>(m.c[i]) * d * d;
    }
    return acc;
}

inline void manufactured_gradient(const ManufacturedQuadratic& m, const double q[],
                                  double g[]) {
    for (int i = 0; i < m.size; ++i) g[i] = 2.0 * m.c[i] * (q[i] - m.a[i]);
}

}  // namespace inverse_design
}  // namespace program_b
