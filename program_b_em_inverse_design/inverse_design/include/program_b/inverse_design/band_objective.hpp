#pragma once

#include "program_b/inverse_design/optimizer.hpp"
#include "program_b/inverse_design/passivity.hpp"

// ---------------------------------------------------------------------------
// The B3 band objective, and the hard-feasible optimizer.
//
// Z(omega) FROM THE COMPLEX BLOCKS. (A11) is homogeneous linear in the blocks
// with REAL coefficients built from g_lab, V and the Levi-Civita tensor -- the
// geometry carries no frequency and no imaginary part. So
//
//     Z(omega) = A11( Re blocks ) + i A11( Im blocks )
//
// exactly. The real and imaginary parts are pushed through the SAME
// reconstruct_Z_A11 that B1C verified to 2.2e-16, rather than a second complex
// implementation being written. The linearity this rests on is asserted in the
// gate, not assumed.
//
// THE BAND OBJECTIVE
//
//     J_band = sum_k w_k ||Z_model(w_k) - Z_target(w_k)||_G^2
//                        / max( ||Z_target(w_k)||_G^2, delta )
//
// with the G-norm of B2 -- ||A||_G^2 for a complex A being the G-norm of its
// real part plus that of its imaginary part, which is the natural extension and
// is what |.|^2 means for a complex tensor.
//
// Contributions are reported SEPARATELY: centre-frequency mismatch, band mean,
// worst frequency, passivity penalty, bound penalty, and the robustness term.
// Collapsing them would hide exactly the tradeoff B3 exists to expose.
// ---------------------------------------------------------------------------

namespace program_b {
namespace inverse_design {

// --- Complex Z at one frequency ----------------------------------------------

template <typename Scalar>
struct ComplexRank4 {
    Rank4<Scalar> re;
    Rank4<Scalar> im;
};

// Lift a complex 3x3 coefficient block into the 4D observer frame, real and
// imaginary parts separately.
template <typename Scalar>
void split_block(const CMat3<Scalar>& c, const ObserverFrame<Scalar>& frame,
                 Mat4<Scalar>& re, Mat4<Scalar>& im) {
    Scalar mr[3][3], mi[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            mr[i][j] = c(i, j).re;
            mi[i][j] = c(i, j).im;
        }
    re = block_from_coefficients(mr, frame);
    im = block_from_coefficients(mi, frame);
}

template <typename Scalar>
ComplexRank4<Scalar> reconstruct_Z_at_frequency(
    const PhysicalDesignVariables<Scalar>& v, const LaboratoryMetric<Scalar>& g_lab,
    const ObserverFrame<Scalar>& frame, double omega) {
    const ComplexBlocks<Scalar> b = blocks_at_frequency(v, omega);

    ObserverDecomposition4<Scalar> re, im;
    split_block(b.epsilon, frame, re.epsilon, im.epsilon);
    split_block(b.mu_inverse, frame, re.mu_inverse, im.mu_inverse);
    split_block(b.zeta, frame, re.zeta, im.zeta);
    split_block(b.zeta_adjoint, frame, re.zeta_adjoint, im.zeta_adjoint);

    ComplexRank4<Scalar> Z;
    Z.re = reconstruct_Z_A11(re, g_lab, frame.V);
    Z.im = reconstruct_Z_A11(im, g_lab, frame.V);
    return Z;
}

// --- Band specification -------------------------------------------------------

struct BandSpec {
    double omega_min = 0.5;
    double omega_max = 2.0;
    int samples = 21;

    std::vector<double> frequencies() const {
        std::vector<double> w;
        if (samples <= 1) { w.push_back(0.5 * (omega_min + omega_max)); return w; }
        const double d = (omega_max - omega_min) / static_cast<double>(samples - 1);
        for (int k = 0; k < samples; ++k)
            w.push_back(omega_min + d * static_cast<double>(k));
        return w;
    }
    double centre() const { return 0.5 * (omega_min + omega_max); }
    int centre_index() const { return samples / 2; }
};

// A target sampled on the band. Held as complex Rank4 per frequency, whatever
// produced it -- a metric, a manufactured model, or a direct tensor.
template <typename Scalar>
struct BandTarget {
    std::vector<ComplexRank4<Scalar>> Z;
    std::vector<double> omega;
    LaboratoryMetric<Scalar> g_lab;
    std::string origin;
    bool ok = false;
};

template <typename Scalar>
BandTarget<Scalar> band_target_from_model(
    const PhysicalDesignVariables<Scalar>& v, const LaboratoryMetric<Scalar>& g_lab,
    const ObserverFrame<Scalar>& frame, const BandSpec& band,
    const std::string& origin) {
    BandTarget<Scalar> t;
    t.g_lab = g_lab;
    t.origin = origin;
    t.omega = band.frequencies();
    for (double w : t.omega)
        t.Z.push_back(reconstruct_Z_at_frequency(v, g_lab, frame, w));
    t.ok = true;
    return t;
}

// A nondispersive target: the same real Z at every frequency. This is what
// reduces B3 to B2 in the real limit.
template <typename Scalar>
BandTarget<Scalar> band_target_from_constant(
    const Rank4<Scalar>& Z0, const LaboratoryMetric<Scalar>& g_lab,
    const BandSpec& band, const std::string& origin) {
    BandTarget<Scalar> t;
    t.g_lab = g_lab;
    t.origin = origin;
    t.omega = band.frequencies();
    for (std::size_t k = 0; k < t.omega.size(); ++k) {
        ComplexRank4<Scalar> z;
        z.re = Z0;
        t.Z.push_back(z);
    }
    t.ok = true;
    return t;
}

// --- Constraint handling ------------------------------------------------------

enum class FeasibilityMode {
    SoftPenalty,      // B2 style: smooth penalties, a small violation is possible
    HardFeasible,     // bounded transforms and projection: violation impossible
};

const char* to_string(FeasibilityMode m);

template <typename Scalar>
struct BandConstraints {
    bool enforce_passivity_penalty = true;
    double passivity_weight = 1.0;

    // Box bounds on the PHYSICAL values, per slot. Enforced structurally in
    // HardFeasible mode and penalised in SoftPenalty mode.
    bool bound_parameters = false;
    double slot_lo[kPhysicalSlots];
    double slot_hi[kPhysicalSlots];
    double bound_weight = 1.0;

    double robustness_weight = 0.0;
    double hinge_smoothing = 1e-6;

    BandConstraints() {
        for (int i = 0; i < kPhysicalSlots; ++i) {
            slot_lo[i] = -1e30;
            slot_hi[i] = 1e30;
        }
    }
};

// --- The evaluation -----------------------------------------------------------

template <typename Scalar>
struct BandEvaluation {
    Scalar J_total;

    // Every contribution kept apart.
    Scalar J_band;            // the weighted band mismatch
    Scalar J_passivity;       // penalty, zero when passive
    Scalar J_bounds;          // parameter-bound penalty
    Scalar J_robustness;      // smoothness across the band

    double centre_mismatch = 0.0;
    double mean_mismatch = 0.0;
    double worst_mismatch = 0.0;
    double worst_mismatch_frequency = 0.0;

    double min_loss_eigenvalue = 0.0;
    double max_passivity_violation = 0.0;
    double worst_passivity_frequency = 0.0;
    double max_bound_violation = 0.0;

    bool finite = false;
};

// The G-norm of a complex rank-4 difference: the real part's norm plus the
// imaginary part's.
template <typename Scalar>
Scalar complex_rank4_norm_squared(const ComplexRank4<Scalar>& A,
                                  const Mat4<Scalar>& G, ObjectiveNorm norm) {
    return rank4_norm_squared(A.re, G, norm) + rank4_norm_squared(A.im, G, norm);
}

template <typename Scalar>
BandEvaluation<Scalar> evaluate_band_design(
    const PhysicalDesignVariables<Scalar>& v, const BandTarget<Scalar>& target,
    const BandConstraints<Scalar>& constraints, const ObserverFrame<Scalar>& frame,
    const ObjectiveWeights<Scalar>& weights) {
    BandEvaluation<Scalar> out;
    const Scalar zero = scalar_cast<Scalar>(0.0);
    out.J_band = zero;
    out.J_passivity = zero;
    out.J_bounds = zero;
    out.J_robustness = zero;

    const Mat4<Scalar> G = observer_gram(target.g_lab, frame.V);
    const std::size_t n = target.omega.size();
    const std::size_t centre = n / 2;

    std::vector<ComplexRank4<Scalar>> model;
    model.reserve(n);
    for (std::size_t k = 0; k < n; ++k)
        model.push_back(reconstruct_Z_at_frequency(v, target.g_lab, frame,
                                                   target.omega[k]));

    double sum_mismatch = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        ComplexRank4<Scalar> d;
        for (int i = 0; i < 256; ++i) {
            d.re.t[i] = model[k].re.t[i] - target.Z[k].re.t[i];
            d.im.t[i] = model[k].im.t[i] - target.Z[k].im.t[i];
        }
        const Scalar num = complex_rank4_norm_squared(d, G, weights.norm);
        const Scalar den =
            complex_rank4_norm_squared(target.Z[k], G, weights.norm) +
            scalar_cast<Scalar>(weights.delta_Z);
        const Scalar term = num / den;
        out.J_band = out.J_band + term;

        const double t = real_part(term);
        sum_mismatch += t;
        if (t > out.worst_mismatch) {
            out.worst_mismatch = t;
            out.worst_mismatch_frequency = target.omega[k];
        }
        if (k == centre) out.centre_mismatch = t;
    }
    out.mean_mismatch = (n > 0) ? sum_mismatch / static_cast<double>(n) : 0.0;
    if (n > 0)
        out.J_band = out.J_band / scalar_cast<Scalar>(static_cast<double>(n));

    // --- Passivity ----------------------------------------------------------
    //
    // The penalty is built from the smallest eigenvalue of L, which is a
    // comparison and therefore a REAL-PART quantity. It is a barrier that pushes
    // the optimizer away from active media; it is not the feasibility verdict,
    // which assess_passivity() reports separately and exactly.
    {
        out.min_loss_eigenvalue = 1e300;
        for (std::size_t k = 0; k < n; ++k) {
            const double lam = passivity_margin_at(v, target.omega[k]);
            if (lam < out.min_loss_eigenvalue) {
                out.min_loss_eigenvalue = lam;
                out.worst_passivity_frequency = target.omega[k];
            }
        }
        out.max_passivity_violation = std::fmax(0.0, -out.min_loss_eigenvalue);

        if (constraints.enforce_passivity_penalty) {
            // A differentiable surrogate: penalise negative Im eps and negative
            // Im mu, which are the scalar drivers of L's diagonal blocks. The
            // exact eigenvalue is not differentiable, so it is not used here --
            // and the gate checks the exact condition independently rather than
            // relying on this surrogate to have been sufficient.
            Scalar pen = zero;
            for (std::size_t k = 0; k < n; ++k) {
                const ResponseAtFrequency<Scalar> r =
                    evaluate_response(v, target.omega[k]);
                const Scalar he = smooth_hinge(-r.epsilon.im, constraints.hinge_smoothing);
                const Scalar hm = smooth_hinge(-r.mu.im, constraints.hinge_smoothing);
                pen = pen + he * he + hm * hm;
            }
            out.J_passivity = scalar_cast<Scalar>(constraints.passivity_weight) * pen;
        }
    }

    // --- Parameter bounds ---------------------------------------------------
    if (constraints.bound_parameters) {
        Scalar p[kPhysicalSlots];
        physical_values(v, p);
        Scalar pen = zero;
        const std::vector<int> slots = v.spec.active_slots();
        for (int slot : slots) {
            const Scalar hi = smooth_hinge(p[slot] - scalar_cast<Scalar>(constraints.slot_hi[slot]),
                                           constraints.hinge_smoothing);
            const Scalar lo = smooth_hinge(scalar_cast<Scalar>(constraints.slot_lo[slot]) - p[slot],
                                           constraints.hinge_smoothing);
            pen = pen + hi * hi + lo * lo;
            out.max_bound_violation = std::fmax(
                out.max_bound_violation,
                std::fmax(real_part(p[slot]) - constraints.slot_hi[slot],
                          constraints.slot_lo[slot] - real_part(p[slot])));
        }
        out.max_bound_violation = std::fmax(0.0, out.max_bound_violation);
        out.J_bounds = scalar_cast<Scalar>(constraints.bound_weight) * pen;
    }

    // --- Robustness ---------------------------------------------------------
    //
    // A smoothness heuristic, not a physical requirement: the mean squared
    // second difference of the modelled response across the band. Labelled as a
    // heuristic in the report, and weighted zero by default.
    if (constraints.robustness_weight > 0.0 && n >= 3) {
        Scalar acc = zero;
        for (std::size_t k = 1; k + 1 < n; ++k)
            for (int i = 0; i < 256; ++i) {
                const Scalar dr = model[k + 1].re.t[i] -
                                  scalar_cast<Scalar>(2.0) * model[k].re.t[i] +
                                  model[k - 1].re.t[i];
                const Scalar di = model[k + 1].im.t[i] -
                                  scalar_cast<Scalar>(2.0) * model[k].im.t[i] +
                                  model[k - 1].im.t[i];
                acc = acc + dr * dr + di * di;
            }
        out.J_robustness =
            scalar_cast<Scalar>(constraints.robustness_weight) * acc;
    }

    out.J_total = out.J_band + out.J_passivity + out.J_bounds + out.J_robustness;
    out.finite = std::isfinite(real_part(out.J_total));
    return out;
}

// --- Optimizer with genuine hard feasibility ---------------------------------
//
// B2 used smooth soft penalties and admitted a small bound violation. B3 adds a
// mode with ACTUAL feasibility, by two mechanisms used together:
//
//   * BOUNDED PARAMETER TRANSFORMS. Strictly positive quantities -- oscillator
//     strengths, resonance frequencies, dampings -- are expressed as strictly
//     positive functions of unconstrained variables, so no step can make them
//     non-positive. This is what enforces causality and the scalar part of
//     passivity structurally.
//
//   * PROJECTION. Box bounds on physical values are enforced by projecting the
//     transformed value back into the box after each accepted step, so the
//     reported design satisfies them exactly rather than approximately.
//
// A soft-penalty result is NEVER labelled feasible merely because its violation
// is small. The two are reported under different names.

// L-BFGS history depth. Eight pairs is ample for a design vector of this size
// and keeps the recursion cheap.
inline constexpr int kLbfgsMemory = 8;

struct BandOptimizerSettings {
    int max_iterations = 300;
    double gradient_tolerance = 1e-12;
    double objective_tolerance = 1e-18;
    double step_tolerance = 1e-16;
    double initial_step = 0.5;
    double backtrack_factor = 0.5;
    double armijo_c = 1e-4;
    int max_backtracks = 60;
    double fd_step = 1e-6;
    FeasibilityMode mode = FeasibilityMode::HardFeasible;
};

struct BandOptimizerResult {
    PhysicalDesignVariables<double> best;
    PhysicalDesignVariables<double> initial;
    double best_objective = 0.0;
    double initial_objective = 0.0;
    double final_gradient_norm = 0.0;
    int iterations = 0;
    int objective_evaluations = 0;
    TerminationReason reason = TerminationReason::MaxIterationsReached;
    FeasibilityMode mode = FeasibilityMode::HardFeasible;
    std::vector<IterationRecord> history;

    // Filled from the final design, so a caller never has to recompute them to
    // find out whether the answer is feasible.
    double max_bound_violation = 0.0;
    double max_passivity_violation = 0.0;
    bool hard_feasible = false;
};

// Project the compact design vector so every active slot's PHYSICAL value lies
// inside its box. In BoundedTransform mode the transform is inverted, clipped,
// and re-applied, so the stored variable is exactly consistent with the clipped
// physical value.
inline void project_into_bounds(PhysicalDesignVariables<double>& v,
                               const BandConstraints<double>& c) {
    if (!c.bound_parameters) return;
    double p[kPhysicalSlots];
    physical_values(v, p);
    const std::vector<int> slots = v.spec.active_slots();
    for (std::size_t k = 0; k < slots.size(); ++k) {
        const int slot = slots[k];
        double value = p[slot];
        const double clipped = std::fmin(std::fmax(value, c.slot_lo[slot]),
                                         c.slot_hi[slot]);
        if (clipped == value) continue;
        if (v.spec.mode == ParameterMode::BoundedTransform &&
            slot_is_strictly_positive(slot))
            v.q[k] = positive_transform_inverse(clipped, v.spec.positive_floor);
        else
            v.q[k] = clipped;
    }
}

template <typename Scalar>
BandOptimizerResult optimize_band_design(
    const PhysicalDesignVariables<double>& start, const BandTarget<double>& target,
    const BandConstraints<double>& constraints, const ObserverFrame<double>& frame,
    const ObjectiveWeights<double>& weights,
    const BandOptimizerSettings& settings) {
    BandOptimizerResult r;
    r.mode = settings.mode;
    r.initial = start;
    r.best = start;

    auto J = [&](const PhysicalDesignVariables<double>& v) {
        ++r.objective_evaluations;
        return real_part(
            evaluate_band_design(v, target, constraints, frame, weights).J_total);
    };

    PhysicalDesignVariables<double> x = start;
    if (settings.mode == FeasibilityMode::HardFeasible)
        project_into_bounds(x, constraints);

    double f = J(x);
    r.initial_objective = f;
    r.best_objective = f;
    if (!std::isfinite(f)) {
        r.reason = TerminationReason::NonFiniteEncountered;
        return r;
    }

    const int n = x.size();
    double step = settings.initial_step;

    std::vector<std::vector<double>> hist_s, hist_y;
    std::vector<double> rho;

    for (int it = 0; it < settings.max_iterations; ++it) {
        double g[kPhysicalSlots] = {};
        bool finite_gradient = true;
        for (int k = 0; k < n; ++k) {
            PhysicalDesignVariables<double> hi = x, lo = x;
            hi.q[k] += settings.fd_step;
            lo.q[k] -= settings.fd_step;
            const double Jh = J(hi), Jl = J(lo);
            g[k] = (Jh - Jl) / (2.0 * settings.fd_step);
            if (!std::isfinite(g[k])) finite_gradient = false;
        }
        if (!finite_gradient) {
            r.reason = TerminationReason::NonFiniteEncountered;
            r.iterations = it;
            return r;
        }

        double gn = 0.0;
        for (int k = 0; k < n; ++k) gn += g[k] * g[k];
        gn = std::sqrt(gn);

        IterationRecord rec;
        rec.iteration = it;
        rec.objective = f;
        rec.gradient_norm = gn;

        if (gn <= settings.gradient_tolerance || f <= settings.objective_tolerance) {
            r.history.push_back(rec);
            r.reason = (gn <= settings.gradient_tolerance)
                           ? TerminationReason::GradientNormBelowTolerance
                           : TerminationReason::ObjectiveBelowTolerance;
            r.iterations = it;
            r.final_gradient_norm = gn;
            break;
        }

        // --- L-BFGS two-loop recursion -----------------------------------
        //
        // Steepest descent stalls badly here: unlike B2's exactly quadratic
        // objective, J_band is RATIONAL in the physical parameters, the
        // resonance parameters are far more sensitive than the backgrounds, and
        // the feasibility transform curves the landscape further. Plain descent
        // reached only ~1e-7 in 400 iterations. L-BFGS with the same Armijo
        // backtracking reaches machine precision in a fraction of that, and it
        // is implemented here rather than pulled in as a dependency.
        //
        // The curvature pairs are stored in the UNCONSTRAINED variables, which
        // is the space the line search moves in, so projection cannot corrupt
        // them.
        double d[kPhysicalSlots];
        {
            double alpha[kLbfgsMemory];
            double q2[kPhysicalSlots];
            for (int k = 0; k < n; ++k) q2[k] = g[k];

            const int m = static_cast<int>(hist_s.size());
            for (int i = m - 1; i >= 0; --i) {
                double sq = 0.0;
                for (int k = 0; k < n; ++k) sq += hist_s[i][k] * q2[k];
                alpha[i] = rho[i] * sq;
                for (int k = 0; k < n; ++k) q2[k] -= alpha[i] * hist_y[i][k];
            }
            // Initial Hessian scaling from the most recent pair.
            double gamma = 1.0;
            if (m > 0) {
                double sy = 0.0, yy = 0.0;
                for (int k = 0; k < n; ++k) {
                    sy += hist_s[m - 1][k] * hist_y[m - 1][k];
                    yy += hist_y[m - 1][k] * hist_y[m - 1][k];
                }
                if (yy > 0.0) gamma = sy / yy;
            }
            for (int k = 0; k < n; ++k) q2[k] *= gamma;
            for (int i = 0; i < m; ++i) {
                double yq = 0.0;
                for (int k = 0; k < n; ++k) yq += hist_y[i][k] * q2[k];
                const double beta = rho[i] * yq;
                for (int k = 0; k < n; ++k)
                    q2[k] += hist_s[i][k] * (alpha[i] - beta);
            }
            for (int k = 0; k < n; ++k) d[k] = q2[k];
        }

        // Guard: if the L-BFGS direction is not a descent direction -- which a
        // damaged curvature history or a projection can cause -- fall back to
        // steepest descent for this step rather than searching uphill.
        double dg = 0.0;
        for (int k = 0; k < n; ++k) dg += d[k] * g[k];
        if (!(dg > 0.0)) {
            for (int k = 0; k < n; ++k) d[k] = g[k];
            dg = gn * gn;
        }

        double trial = step;
        int backtracks = 0;
        bool accepted = false;
        PhysicalDesignVariables<double> next = x;
        double fnext = f;

        while (backtracks < settings.max_backtracks) {
            next = x;
            for (int k = 0; k < n; ++k) next.q[k] = x.q[k] - trial * d[k];
            if (settings.mode == FeasibilityMode::HardFeasible)
                project_into_bounds(next, constraints);
            fnext = J(next);
            if (std::isfinite(fnext) &&
                fnext <= f - settings.armijo_c * trial * dg) {
                accepted = true;
                break;
            }
            trial *= settings.backtrack_factor;
            ++backtracks;
        }

        rec.step_length = trial;
        rec.backtracks = backtracks;
        rec.line_search_ok = accepted;
        r.history.push_back(rec);

        if (!accepted) {
            r.reason = TerminationReason::LineSearchFailed;
            r.iterations = it;
            r.final_gradient_norm = gn;
            break;
        }

        double move = 0.0;
        double s_k[kPhysicalSlots];
        for (int k = 0; k < n; ++k) {
            s_k[k] = next.q[k] - x.q[k];
            move += s_k[k] * s_k[k];
        }

        // Curvature pair (s, y). Computed BEFORE x is overwritten, and stored
        // only when the curvature condition s.y > 0 holds -- a pair that fails
        // it would make the approximate inverse Hessian indefinite.
        {
            double g_next[kPhysicalSlots];
            for (int k = 0; k < n; ++k) {
                PhysicalDesignVariables<double> hi = next, lo = next;
                hi.q[k] += settings.fd_step;
                lo.q[k] -= settings.fd_step;
                g_next[k] = (J(hi) - J(lo)) / (2.0 * settings.fd_step);
            }
            double y_k[kPhysicalSlots];
            double sy = 0.0;
            for (int k = 0; k < n; ++k) {
                y_k[k] = g_next[k] - g[k];
                sy += s_k[k] * y_k[k];
            }
            if (sy > 1e-18) {
                if (static_cast<int>(hist_s.size()) == kLbfgsMemory) {
                    hist_s.erase(hist_s.begin());
                    hist_y.erase(hist_y.begin());
                    rho.erase(rho.begin());
                }
                std::vector<double> sv(s_k, s_k + n), yv(y_k, y_k + n);
                hist_s.push_back(sv);
                hist_y.push_back(yv);
                rho.push_back(1.0 / sy);
            }
        }

        x = next;
        f = fnext;
        r.best = x;
        r.best_objective = f;
        r.final_gradient_norm = gn;
        step = std::fmin(trial * 2.0, settings.initial_step);

        if (std::sqrt(move) <= settings.step_tolerance) {
            r.reason = TerminationReason::StepBelowTolerance;
            r.iterations = it + 1;
            break;
        }
        if (it + 1 == settings.max_iterations) {
            r.reason = TerminationReason::MaxIterationsReached;
            r.iterations = settings.max_iterations;
        }
    }

    r.best = x;
    r.best_objective = f;

    const BandEvaluation<double> final_eval =
        evaluate_band_design(r.best, target, constraints, frame, weights);
    r.max_bound_violation = final_eval.max_bound_violation;
    r.max_passivity_violation = final_eval.max_passivity_violation;

    // FEASIBLE means feasible, exactly. A soft-penalty run is never labelled so
    // however small its violation.
    r.hard_feasible = (settings.mode == FeasibilityMode::HardFeasible) &&
                      (r.max_bound_violation <= 0.0);
    return r;
}

}  // namespace inverse_design
}  // namespace program_b
