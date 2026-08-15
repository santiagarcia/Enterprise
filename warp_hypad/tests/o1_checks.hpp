#pragma once

#include <string>
#include <vector>

#include "o0d_checks.hpp"

namespace warp_hypad {
namespace checks {

// ---------------------------------------------------------------------------
// O1: safeguarded optimization of the inferred SOURCE REQUIREMENT for a fixed
// warp-bubble kinematic family.
//
// THE THREE NAMES FOR THE SURROGATE, KEPT DISTINCT EVERYWHERE.
//
//     raw_E_neg_tau(y; D)  = integral softplus_tau(-rho) dV
//     renormalized         = raw - tau log2 V(D)
//     delta                = raw(y) - raw(0)      on the same D
//
// O0D showed softplus_tau(0) = tau log 2 > 0, so the raw integral carries a
// y-independent vacuum floor proportional to the domain volume. The floor
// cancels from renormalized and from delta identically, and all three have the
// SAME gradient in y -- which is why the verified OTI gradient is used
// unaltered. The optimizer minimises the RENORMALIZED surrogate; the raw value
// is recorded but never called an energy.
//
// The exact reporting objective is E_neg = integral max(-rho, 0) dV, and it is
// the only quantity a reduction claim may refer to.
// ---------------------------------------------------------------------------

struct OptimizerConfig {
    int lbfgs_memory = 8;
    double c1_armijo = 1.0e-4;
    double trust_radius_initial = 0.25;
    double trust_radius_min = 1.0e-6;
    double trust_radius_max = 4.0;
    double trust_expand = 2.0;
    double trust_shrink = 0.5;
    double backtrack = 0.5;
    int max_line_search = 12;
    int max_accepted = 40;
    int max_objective_evals = 250;
    double gradient_norm_tol = 1.0e-8;
    double step_norm_tol = 1.0e-8;
    int stall_window = 5;
};

// One accepted design, with every diagnostic Phase 4 requires.
struct Iterate {
    int iteration = 0;
    int objective_evals = 0;
    int gradient_evals = 0;

    double raw_working = 0.0, renorm_working = 0.0, delta_working = 0.0;
    double raw_validation = 0.0, renorm_validation = 0.0,
           delta_validation = 0.0;

    double E_neg = 0.0, E_plus = 0.0, E_net = 0.0, E_abs = 0.0;
    double J_NEC = 0.0, N_min = 0.0, K_max = 0.0;   // on the WORKING set

    double grad_norm = 0.0, y_norm = 0.0, step_norm = 0.0;
    double trust_radius = 0.0, alpha = 0.0;
    int line_search_evals = 0;

    double max_correction = 0.0;
    double max_d1 = 0.0, max_d2 = 0.0, max_d3 = 0.0;
    std::size_t invalid_points = 0, observer_failures = 0;

    std::vector<double> y;
    std::vector<double> gradient;
};

struct RejectedTrial {
    int iteration = 0;
    int attempt = 0;
    double alpha = 0.0;
    double trust_radius = 0.0;
    double step_norm = 0.0;
    double raw_working = 0.0, renorm_working = 0.0;
    double raw_validation = 0.0, renorm_validation = 0.0;
    double E_neg = 0.0;
    double drift = 0.0;
    double interior = 0.0, collar = 0.0;
    std::string classification;   // comparable | unresolved | behavior_changing | invalid
    std::string reason;           // machine-readable
};

struct LineSearchEntry {
    int iteration = 0;
    int attempt = 0;
    double alpha = 0.0;
    double f_trial = 0.0;
    double f_current = 0.0;
    double armijo_target = 0.0;
    bool armijo_satisfied = false;
    double directional_derivative = 0.0;
};

struct OptimizerResult {
    bool ran = false;
    std::string skip_reason;

    std::vector<Iterate> accepted;
    std::vector<RejectedTrial> rejected;
    std::vector<LineSearchEntry> line_search;

    int objective_evals = 0;
    int gradient_evals = 0;
    int lbfgs_fallbacks = 0;
    int curvature_pairs_accepted = 0;
    int curvature_pairs_rejected = 0;
    double trust_radius_min_seen = 0.0;
    double trust_radius_max_seen = 0.0;

    std::string termination_reason;
    bool terminated_by_limit = false;
    bool converged = false;

    // Baseline, on both evaluators.
    double baseline_raw_working = 0.0, baseline_renorm_working = 0.0;
    double baseline_raw_validation = 0.0, baseline_renorm_validation = 0.0;
    double baseline_E_neg = 0.0;
};

// The engine. Deterministic: no randomness, no time-dependence, no adaptive
// tolerance that depends on what it has seen.
OptimizerResult run_o1_optimization(const F1wSystem& sys,
                                    const OptimizerConfig& cfg,
                                    const std::vector<double>& y_start,
                                    int max_accepted_override,
                                    std::string* timing_csv);

// Records and gate.
void write_o1_records(const OptimizerResult& r, const F1wSystem& sys,
                      const std::string& out_dir);
void run_o1_checks(VerificationReport& report, const std::string& out_dir);

inline const char* o1_scope() {
    return "O1 optimizes the inferred source requirement for a FIXED warp-bubble "
           "kinematic family: the F1W shift-only wall correction at fixed "
           "velocity, fixed protected-cabin behavior, the ADM Eulerian observer "
           "and the total-effective source convention with Lambda = 0. It does "
           "not demonstrate faster-than-light travel, a causal shortcut, "
           "reduced travel time, horizon removal, physical realizability, "
           "dynamic stability, satisfaction of quantum energy inequalities, or "
           "a global optimum.";
}

}  // namespace checks
}  // namespace warp_hypad
