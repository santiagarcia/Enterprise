#pragma once

#include <vector>

#include "program_b/inverse_design/design_gradient.hpp"

// ---------------------------------------------------------------------------
// A deliberately small deterministic optimizer: gradient descent with a
// backtracking Armijo line search. No external dependency, nothing stochastic,
// nothing adaptive that could not be reproduced from the recorded inputs.
//
// It runs on the FINITE-DIFFERENCE gradient by default, so every inverse-design
// demonstration is runnable with the HyPAD backend disabled. That is not a
// weakness of the gradient work: the HyPAD gradient is verified against this one
// and against the analytic case, and either can drive the optimizer.
//
// Everything a run needs to be reproducible or diagnosed after the fact is
// recorded: the initialization, the objective and gradient-norm histories, the
// per-iteration line-search outcome, and a typed termination reason.
// ---------------------------------------------------------------------------

namespace program_b {
namespace inverse_design {

enum class TerminationReason {
    GradientNormBelowTolerance,
    ObjectiveBelowTolerance,
    StepBelowTolerance,
    LineSearchFailed,
    MaxIterationsReached,
    NonFiniteEncountered,
};

const char* to_string(TerminationReason r);

struct OptimizerSettings {
    int max_iterations = 400;
    double gradient_tolerance = 1e-10;
    double objective_tolerance = 1e-16;
    double step_tolerance = 1e-16;

    double initial_step = 1.0;
    double backtrack_factor = 0.5;
    double armijo_c = 1e-4;
    int max_backtracks = 60;

    // Step used by the finite-difference gradient. Chosen for a quadratic
    // objective, where centred differences carry no truncation error.
    double fd_step = 1e-5;
};

struct IterationRecord {
    int iteration = 0;
    double objective = 0.0;
    double gradient_norm = 0.0;
    double step_length = 0.0;
    int backtracks = 0;
    bool line_search_ok = false;
};

struct OptimizerResult {
    ConstitutiveDesignVariables<double> best;
    double best_objective = 0.0;
    double final_gradient_norm = 0.0;
    int iterations = 0;
    int objective_evaluations = 0;
    TerminationReason reason = TerminationReason::MaxIterationsReached;
    std::vector<IterationRecord> history;

    // The initialization, recorded so a run can be replayed exactly.
    ConstitutiveDesignVariables<double> initial;
    double initial_objective = 0.0;
};

// Deterministic perturbation of a known design, used to start the recovery
// problems away from the answer. A fixed integer recurrence rather than a random
// engine: no seed to lose, and identical on every platform.
inline ConstitutiveDesignVariables<double> perturb_design(
    const ConstitutiveDesignVariables<double>& v, double amplitude, int seed) {
    ConstitutiveDesignVariables<double> out = v;
    unsigned int s = static_cast<unsigned int>(seed) * 2654435761u + 1u;
    for (int i = 0; i < v.size(); ++i) {
        s = s * 1664525u + 1013904223u;
        const double u = static_cast<double>((s >> 8) & 0xFFFFu) / 65535.0;
        out.q[i] = v.q[i] + amplitude * (2.0 * u - 1.0);
    }
    return out;
}

inline OptimizerResult optimize_design(
    const ConstitutiveDesignVariables<double>& start,
    const ConstitutiveTarget<double>& target,
    const DesignConstraints<double>& constraints,
    const ObserverFrame<double>& frame, const ObjectiveWeights<double>& weights,
    const OptimizerSettings& settings) {
    OptimizerResult r;
    r.initial = start;
    r.best = start;

    auto J = [&](const ConstitutiveDesignVariables<double>& v) {
        ++r.objective_evaluations;
        return real_part(
            evaluate_design(v, target, constraints, frame, weights).J_total);
    };

    ConstitutiveDesignVariables<double> x = start;
    double f = J(x);
    r.initial_objective = f;
    r.best_objective = f;

    if (!std::isfinite(f)) {
        r.reason = TerminationReason::NonFiniteEncountered;
        r.final_gradient_norm = 0.0;
        return r;
    }

    const int n = x.size();
    double step = settings.initial_step;

    for (int it = 0; it < settings.max_iterations; ++it) {
        const GradientResult g = finite_difference_gradient(
            x, target, constraints, frame, weights, settings.fd_step);
        r.objective_evaluations += g.evaluations;

        bool finite_gradient = g.ok;
        for (int k = 0; k < n && finite_gradient; ++k)
            if (!std::isfinite(g.g[k])) finite_gradient = false;
        if (!finite_gradient) {
            r.reason = TerminationReason::NonFiniteEncountered;
            r.iterations = it;
            return r;
        }

        const double gn = g.norm();
        IterationRecord rec;
        rec.iteration = it;
        rec.objective = f;
        rec.gradient_norm = gn;

        if (gn <= settings.gradient_tolerance) {
            r.history.push_back(rec);
            r.reason = TerminationReason::GradientNormBelowTolerance;
            r.iterations = it;
            r.final_gradient_norm = gn;
            r.best = x;
            r.best_objective = f;
            return r;
        }
        if (f <= settings.objective_tolerance) {
            r.history.push_back(rec);
            r.reason = TerminationReason::ObjectiveBelowTolerance;
            r.iterations = it;
            r.final_gradient_norm = gn;
            r.best = x;
            r.best_objective = f;
            return r;
        }

        // Backtracking Armijo search along -g.
        double trial = step;
        int backtracks = 0;
        bool accepted = false;
        ConstitutiveDesignVariables<double> next = x;
        double fnext = f;

        while (backtracks < settings.max_backtracks) {
            next = x;
            for (int k = 0; k < n; ++k) next.q[k] = x.q[k] - trial * g.g[k];
            fnext = J(next);
            if (std::isfinite(fnext) &&
                fnext <= f - settings.armijo_c * trial * gn * gn) {
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
            r.best = x;
            r.best_objective = f;
            return r;
        }

        double move = 0.0;
        for (int k = 0; k < n; ++k) {
            const double d = next.q[k] - x.q[k];
            move += d * d;
        }
        move = std::sqrt(move);

        x = next;
        f = fnext;
        r.best = x;
        r.best_objective = f;
        r.final_gradient_norm = gn;

        // Grow the trial step slightly on success, so a run that starts far away
        // does not stay stuck at a tiny step. Bounded, and deterministic.
        step = std::fmin(trial * 2.0, settings.initial_step);

        if (move <= settings.step_tolerance) {
            r.reason = TerminationReason::StepBelowTolerance;
            r.iterations = it + 1;
            return r;
        }
    }

    r.reason = TerminationReason::MaxIterationsReached;
    r.iterations = settings.max_iterations;
    return r;
}

}  // namespace inverse_design
}  // namespace program_b
