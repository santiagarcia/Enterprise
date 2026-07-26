#pragma once

#include "program_b/inverse_design/design_target.hpp"

// ---------------------------------------------------------------------------
// B2 objectives, norms and constraints.
//
// TWO NORMS, AND WHICH IS WHICH.
//
//   COORDINATE-BASIS norm, labelled as such because it is NOT invariant:
//
//       ||A||^2_coord = sum_{abcd} (A^{abcd})^2
//
//   OBSERVER-REFERRED norm, which is a genuine scalar:
//
//       ||A||^2_V = A^{abcd} A^{efgh} G_ae G_bf G_cg G_dh
//       G_ab = g_lab_{ab} + 2 V_a V_b
//
// G is used rather than g_lab itself for a mathematical reason, not a stylistic
// one: contracting with a Lorentzian metric is INDEFINITE, so ||A||^2_g can
// vanish or go negative for A != 0 and is unusable as an objective. G = h + V(x)V
// is g_lab with the observer's timelike direction flipped, hence positive
// definite, so ||A||^2_V is a positive definite quadratic form and a scalar under
// coordinate change, provided V transforms with the coordinates.
//
// The two coincide EXACTLY when g_lab = eta and V = (1,0,0,0), because there
// G = diag(1,1,1,1). That is the analogue configuration, which is why the
// coordinate-basis norm is usable there -- and why calling it invariant would
// still be wrong. Both facts are asserted in the gate.
//
// CONSTRAINTS ARE SMOOTH. No branch is taken on a differentiable scalar inside
// the objective. One-sided penalties use
//
//     hinge(x) = (x + sqrt(x^2 + tau^2))/2 - tau/2
//
// which is C-infinity, needs only sqrt, and requires no comparison on a Scalar.
// Eigenvalues and condition numbers need an iterative algorithm with
// comparisons, so they are REPORTED as real-part diagnostics and never used as
// penalties; where a bound on them is wanted the AD-safe surrogate is used and
// labelled as a bound.
// ---------------------------------------------------------------------------

namespace program_b {
namespace inverse_design {

// --- Norms -------------------------------------------------------------------

// G_ab = g_lab_{ab} + 2 V_a V_b. Positive definite whenever V is unit timelike.
template <typename Scalar>
Mat4<Scalar> observer_gram(const LaboratoryMetric<Scalar>& g_lab,
                           const Scalar V_up[4]) {
    Scalar V_dn[4];
    lower_vector(g_lab.value, V_up, V_dn);
    Mat4<Scalar> G;
    const Scalar two = scalar_cast<Scalar>(2.0);
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            G(a, b) = g_lab.value(a, b) + two * V_dn[a] * V_dn[b];
    return G;
}

enum class ObjectiveNorm {
    CoordinateBasis,     // sum of squared components; NOT invariant
    ObserverReferred,    // contracted with G = g_lab + 2 V(x)V; a scalar
};

const char* to_string(ObjectiveNorm n);
const char* describe(ObjectiveNorm n);

// ||A||^2 for a rank-4, in whichever norm.
template <typename Scalar>
Scalar rank4_norm_squared(const Rank4<Scalar>& A, const Mat4<Scalar>& G,
                          ObjectiveNorm norm) {
    Scalar acc = scalar_cast<Scalar>(0.0);
    if (norm == ObjectiveNorm::CoordinateBasis) {
        for (int i = 0; i < 256; ++i) acc = acc + A.t[i] * A.t[i];
        return acc;
    }
    // A^{abcd} A^{efgh} G_ae G_bf G_cg G_dh, contracted one index pair at a
    // time so the cost is 4 x 4^5 rather than 4^8.
    Rank4<Scalar> t1, t2, t3, t4;
    for (int e = 0; e < 4; ++e)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d) {
                    Scalar s = scalar_cast<Scalar>(0.0);
                    for (int a = 0; a < 4; ++a) s = s + A(a, b, c, d) * G(a, e);
                    t1(e, b, c, d) = s;
                }
    for (int e = 0; e < 4; ++e)
        for (int f = 0; f < 4; ++f)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d) {
                    Scalar s = scalar_cast<Scalar>(0.0);
                    for (int b = 0; b < 4; ++b) s = s + t1(e, b, c, d) * G(b, f);
                    t2(e, f, c, d) = s;
                }
    for (int e = 0; e < 4; ++e)
        for (int f = 0; f < 4; ++f)
            for (int g = 0; g < 4; ++g)
                for (int d = 0; d < 4; ++d) {
                    Scalar s = scalar_cast<Scalar>(0.0);
                    for (int c = 0; c < 4; ++c) s = s + t2(e, f, c, d) * G(c, g);
                    t3(e, f, g, d) = s;
                }
    for (int e = 0; e < 4; ++e)
        for (int f = 0; f < 4; ++f)
            for (int g = 0; g < 4; ++g)
                for (int h = 0; h < 4; ++h) {
                    Scalar s = scalar_cast<Scalar>(0.0);
                    for (int d = 0; d < 4; ++d) s = s + t3(e, f, g, d) * G(d, h);
                    t4(e, f, g, h) = s;
                }
    for (int i = 0; i < 256; ++i) acc = acc + A.t[i] * t4.t[i];
    return acc;
}

template <typename Scalar>
Scalar mat4_norm_squared(const Mat4<Scalar>& A, const Mat4<Scalar>& G,
                         ObjectiveNorm norm) {
    Scalar acc = scalar_cast<Scalar>(0.0);
    if (norm == ObjectiveNorm::CoordinateBasis) {
        for (int i = 0; i < 16; ++i) acc = acc + A.m[i] * A.m[i];
        return acc;
    }
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int e = 0; e < 4; ++e)
                for (int f = 0; f < 4; ++f)
                    acc = acc + A(a, b) * A(e, f) * G(a, e) * G(b, f);
    return acc;
}

// --- Smooth one-sided penalty ------------------------------------------------

// hinge(x) ~ max(x, 0), smooth everywhere, no comparison on a Scalar.
// hinge(0) = 0 exactly, and hinge(x) -> x for x >> tau.
template <typename Scalar>
Scalar smooth_hinge(const Scalar& x, double tau) {
    using std::sqrt;
    const Scalar t2 = scalar_cast<Scalar>(tau * tau);
    const Scalar half = scalar_cast<Scalar>(0.5);
    return half * (x + sqrt(x * x + t2)) - scalar_cast<Scalar>(0.5 * tau);
}

// --- Constraints -------------------------------------------------------------

template <typename Scalar>
struct DesignConstraints {
    bool active = false;

    // Per-entry magnitude bounds on the 3x3 coefficient matrices.
    double max_epsilon_entry = 1e6;
    double max_mu_inverse_entry = 1e6;
    double max_zeta_entry = 1e6;

    // Frobenius-norm bounds. An upper bound on the spectral radius, AD-safe
    // where an eigenvalue bound is not.
    double max_epsilon_frobenius = 1e6;
    double max_mu_inverse_frobenius = 1e6;
    double max_zeta_frobenius = 1e6;

    // Determinant floors, from the closed-form 3x3 determinant: a lower bound on
    // the product of eigenvalues, again AD-safe.
    bool enforce_determinant_floor = false;
    double min_epsilon_determinant = 1e-6;
    double min_mu_inverse_determinant = 1e-6;

    double penalty_weight = 1.0;
    double hinge_smoothing = 1e-6;
};

template <typename Scalar>
Scalar determinant3(const Scalar m[3][3]) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
           m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

template <typename Scalar>
Scalar frobenius_squared3(const Scalar m[3][3]) {
    Scalar acc = scalar_cast<Scalar>(0.0);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) acc = acc + m[i][j] * m[i][j];
    return acc;
}

// --- Diagnostics -------------------------------------------------------------
//
// Real parts only, and deliberately so: these are eigenvalue and conditioning
// comparisons, which is exactly the family of quantity that must not be inside a
// differentiated objective.

struct BlockDiagnostics {
    double eig_min_abs = 0.0;
    double eig_max_abs = 0.0;
    double determinant = 0.0;
    double max_component = 0.0;
    double frobenius = 0.0;
    // ||A||_F ||A^-1||_F, an UPPER BOUND on the condition number, not cond.
    double condition_bound = 0.0;
};

template <typename Scalar>
BlockDiagnostics diagnose_block(const Scalar m[3][3]) {
    BlockDiagnostics d;
    double a[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            a[i][j] = real_part(m[i][j]);
            d.max_component = std::fmax(d.max_component, std::fabs(a[i][j]));
            d.frobenius += a[i][j] * a[i][j];
        }
    d.frobenius = std::sqrt(d.frobenius);

    double eig[3];
    symmetric_eigenvalues_3x3(a, eig);
    d.eig_min_abs = std::fmin(std::fabs(eig[0]),
                              std::fmin(std::fabs(eig[1]), std::fabs(eig[2])));
    d.eig_max_abs = std::fmax(std::fabs(eig[0]),
                              std::fmax(std::fabs(eig[1]), std::fabs(eig[2])));
    d.determinant = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
                    a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
                    a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);

    if (std::fabs(d.determinant) > 0.0) {
        // ||adj A||_F / |det A| = ||A^-1||_F
        double adj = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                const int i1 = (i + 1) % 3, i2 = (i + 2) % 3;
                const int j1 = (j + 1) % 3, j2 = (j + 2) % 3;
                const double c = a[i1][j1] * a[i2][j2] - a[i1][j2] * a[i2][j1];
                adj += c * c;
            }
        d.condition_bound = d.frobenius * std::sqrt(adj) / std::fabs(d.determinant);
    } else {
        d.condition_bound = std::numeric_limits<double>::infinity();
    }
    return d;
}

// --- The evaluation ----------------------------------------------------------

enum class PrimaryObjective {
    CovariantTensorMismatch,   // J_Z
    BlockMismatch,             // J_blocks
};

const char* to_string(PrimaryObjective p);

template <typename Scalar>
struct ObjectiveWeights {
    PrimaryObjective primary = PrimaryObjective::CovariantTensorMismatch;
    ObjectiveNorm norm = ObjectiveNorm::CoordinateBasis;
    double w_epsilon = 1.0;
    double w_mu_inverse = 1.0;
    double w_zeta = 1.0;
    double delta_Z = 1e-12;      // floor on the target scale in J_Z
    double delta_block = 1e-12;
};

template <typename Scalar>
struct DesignEvaluation {
    // Reported separately. Never collapsed into one number by the evaluator.
    Scalar J_total;
    Scalar J_Z;
    Scalar J_blocks;
    Scalar J_epsilon;
    Scalar J_mu_inverse;
    Scalar J_zeta;
    Scalar J_reciprocity;   // diagnostic; not in J_total
    Scalar J_penalty;

    Rank4<Scalar> Z_design;
    ObserverDecomposition4<Scalar> blocks;
    BlockCoefficients<Scalar> coefficients;

    BlockDiagnostics epsilon_diagnostics;
    BlockDiagnostics mu_inverse_diagnostics;
    BlockDiagnostics zeta_diagnostics;

    double spatiality_residual = 0.0;
    double pair_symmetry_residual_design = 0.0;
    bool finite = false;

    // Which bounds are pressing, for the constraint-activity record.
    bool epsilon_entry_active = false;
    bool mu_inverse_entry_active = false;
    bool zeta_entry_active = false;
    bool determinant_floor_active = false;
};

template <typename Scalar>
DesignEvaluation<Scalar> evaluate_design(
    const ConstitutiveDesignVariables<Scalar>& variables,
    const ConstitutiveTarget<Scalar>& target,
    const DesignConstraints<Scalar>& constraints,
    const ObserverFrame<Scalar>& frame,
    const ObjectiveWeights<Scalar>& weights) {
    DesignEvaluation<Scalar> out;
    const Scalar zero = scalar_cast<Scalar>(0.0);
    out.J_total = zero;
    out.J_Z = zero;
    out.J_blocks = zero;
    out.J_epsilon = zero;
    out.J_mu_inverse = zero;
    out.J_zeta = zero;
    out.J_reciprocity = zero;
    out.J_penalty = zero;

    out.coefficients = coefficients_from_parameters(variables);
    out.blocks = blocks_from_parameters(variables, frame);
    out.Z_design = reconstruct_Z_A11(out.blocks, target.g_lab, frame.V);

    const Mat4<Scalar> G = observer_gram(target.g_lab, frame.V);

    // --- J_Z ----------------------------------------------------------------
    {
        Rank4<Scalar> diff;
        for (int i = 0; i < 256; ++i)
            diff.t[i] = out.Z_design.t[i] - target.Z.t[i];
        const Scalar num = rank4_norm_squared(diff, G, weights.norm);
        const Scalar den0 = rank4_norm_squared(target.Z, G, weights.norm);
        const Scalar den = den0 + scalar_cast<Scalar>(weights.delta_Z);
        out.J_Z = num / den;
    }

    // --- J_blocks -----------------------------------------------------------
    {
        const ObserverDecomposition4<Scalar>& tb = target.blocks;
        auto term = [&](const Mat4<Scalar>& d, const Mat4<Scalar>& t) {
            Mat4<Scalar> diff;
            for (int i = 0; i < 16; ++i) diff.m[i] = d.m[i] - t.m[i];
            const Scalar num = mat4_norm_squared(diff, G, weights.norm);
            const Scalar den = mat4_norm_squared(t, G, weights.norm) +
                               scalar_cast<Scalar>(weights.delta_block);
            return num / den;
        };
        out.J_epsilon = term(out.blocks.epsilon, tb.epsilon);
        out.J_mu_inverse = term(out.blocks.mu_inverse, tb.mu_inverse);
        out.J_zeta = term(out.blocks.zeta, tb.zeta);
        out.J_blocks = scalar_cast<Scalar>(weights.w_epsilon) * out.J_epsilon +
                       scalar_cast<Scalar>(weights.w_mu_inverse) * out.J_mu_inverse +
                       scalar_cast<Scalar>(weights.w_zeta) * out.J_zeta;
    }

    // --- J_reciprocity, a diagnostic ---------------------------------------
    {
        Rank4<Scalar> diff;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int c = 0; c < 4; ++c)
                    for (int d = 0; d < 4; ++d)
                        diff(a, b, c, d) =
                            out.Z_design(a, b, c, d) - out.Z_design(c, d, a, b);
        out.J_reciprocity = rank4_norm_squared(diff, G, weights.norm);
    }

    // --- Penalties ----------------------------------------------------------
    if (constraints.active) {
        const double tau = constraints.hinge_smoothing;
        Scalar pen = zero;

        auto entry_penalty = [&](const Scalar m[3][3], double bound) {
            Scalar acc = zero;
            const Scalar b = scalar_cast<Scalar>(bound);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    // Two-sided: |m| <= bound, written without fabs on a Scalar.
                    const Scalar hi = smooth_hinge(m[i][j] - b, tau);
                    const Scalar lo = smooth_hinge(-m[i][j] - b, tau);
                    acc = acc + hi * hi + lo * lo;
                }
            return acc;
        };
        auto frobenius_penalty = [&](const Scalar m[3][3], double bound) {
            const Scalar over = frobenius_squared3(m) -
                               scalar_cast<Scalar>(bound * bound);
            const Scalar h = smooth_hinge(over, tau);
            return h * h;
        };

        pen = pen + entry_penalty(out.coefficients.epsilon,
                                  constraints.max_epsilon_entry);
        pen = pen + entry_penalty(out.coefficients.mu_inverse,
                                  constraints.max_mu_inverse_entry);
        pen = pen + entry_penalty(out.coefficients.zeta,
                                  constraints.max_zeta_entry);
        pen = pen + frobenius_penalty(out.coefficients.epsilon,
                                      constraints.max_epsilon_frobenius);
        pen = pen + frobenius_penalty(out.coefficients.mu_inverse,
                                      constraints.max_mu_inverse_frobenius);
        pen = pen + frobenius_penalty(out.coefficients.zeta,
                                      constraints.max_zeta_frobenius);

        if (constraints.enforce_determinant_floor) {
            const Scalar de = scalar_cast<Scalar>(constraints.min_epsilon_determinant) -
                              determinant3(out.coefficients.epsilon);
            const Scalar dm = scalar_cast<Scalar>(constraints.min_mu_inverse_determinant) -
                              determinant3(out.coefficients.mu_inverse);
            const Scalar he = smooth_hinge(de, tau);
            const Scalar hm = smooth_hinge(dm, tau);
            pen = pen + he * he + hm * hm;
        }

        out.J_penalty = scalar_cast<Scalar>(constraints.penalty_weight) * pen;
    }

    out.J_total = (weights.primary == PrimaryObjective::CovariantTensorMismatch)
                      ? out.J_Z + out.J_penalty
                      : out.J_blocks + out.J_penalty;

    // --- Diagnostics --------------------------------------------------------
    out.epsilon_diagnostics = diagnose_block(out.coefficients.epsilon);
    out.mu_inverse_diagnostics = diagnose_block(out.coefficients.mu_inverse);
    out.zeta_diagnostics = diagnose_block(out.coefficients.zeta);

    if (constraints.active) {
        out.epsilon_entry_active = out.epsilon_diagnostics.max_component >=
                                   0.999 * constraints.max_epsilon_entry;
        out.mu_inverse_entry_active = out.mu_inverse_diagnostics.max_component >=
                                      0.999 * constraints.max_mu_inverse_entry;
        out.zeta_entry_active = out.zeta_diagnostics.max_component >=
                                0.999 * constraints.max_zeta_entry;
        out.determinant_floor_active =
            constraints.enforce_determinant_floor &&
            (out.epsilon_diagnostics.determinant <=
                 1.001 * constraints.min_epsilon_determinant ||
             out.mu_inverse_diagnostics.determinant <=
                 1.001 * constraints.min_mu_inverse_determinant);
    }

    {
        double s1 = 0.0, s2 = 0.0, worst = 0.0;
        const Mat4<Scalar>* blocks[4] = {&out.blocks.epsilon, &out.blocks.mu_inverse,
                                         &out.blocks.zeta, &out.blocks.zeta_adjoint};
        for (int k = 0; k < 4; ++k) {
            spatiality_residuals(*blocks[k], target.g_lab.value, frame.V, s1, s2);
            worst = std::fmax(worst, std::fmax(s1, s2));
        }
        out.spatiality_residual = worst;
    }
    out.pair_symmetry_residual_design = pair_symmetry_residual(out.Z_design);

    out.finite = std::isfinite(real_part(out.J_total)) &&
                 std::isfinite(real_part(out.J_Z)) &&
                 std::isfinite(real_part(out.J_blocks));
    return out;
}

}  // namespace inverse_design
}  // namespace program_b
