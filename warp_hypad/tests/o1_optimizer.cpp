#include "o1_checks.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <deque>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#if defined(WARP_HYPAD_HAVE_OTILIB)
#include "warp_hypad/sensitivity/otilib_extract.hpp"
#include "warp_hypad/sensitivity/otilib_math.hpp"
#include "warp_hypad/sensitivity/otilib_scalar.hpp"
#include "warp_hypad/sensitivity/otilib_seed.hpp"
#endif

// ---------------------------------------------------------------------------
// The O1 engine: limited-memory BFGS in whitened y, safeguarded by an adaptive
// trust radius, an Armijo line search on the WORKING renormalized surrogate,
// and a twelve-condition acceptance test against the VALIDATION evaluator.
//
// WHAT L-BFGS IS AND IS NOT. It builds an inverse-Hessian APPROXIMATION from
// differences of gradients. It never forms, stores or claims the true Hessian,
// and it cannot: O0B-C established that mixed second derivatives between
// coordinates in different OTI batches are not available at all. The two-loop
// recursion below uses only gradient differences, which are available exactly.
//
// DETERMINISM IS A DESIGN CONSTRAINT, NOT AN ASPIRATION. No randomness, no
// wall-clock input, no tolerance that adapts to what the run has already seen.
// Every branch is a function of the iterate history alone, so two runs from the
// same start produce the same sequence bit for bit -- which Phase 8 checks.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {
namespace {

using design::IntegrationSet;
using design::ObjectiveValues;
using design::ShiftFamilyConfig;
using design::SourceSemantics;

inline constexpr int kOtiDirections = 3;
const double kLog2 = 0.69314718055994530942;

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}
double nrm(const std::vector<double>& a) { return std::sqrt(dot(a, a)); }

// The protected-region lock, measured on a fixed deterministic probe set.
struct BehaviorLock {
    double interior = 0.0;
    double collar = 0.0;
    double max_correction = 0.0;
    double d1 = 0.0, d2 = 0.0, d3 = 0.0;
};

BehaviorLock measure_lock(const ShiftFamilyConfig& cfg,
                          const design::AxisymmetricBasis& basis,
                          const std::vector<double>& c) {
    BehaviorLock b;
    for (double q : c) b.max_correction = std::fmax(b.max_correction, std::fabs(q));
    for (int i = 0; i <= 12; ++i)
        for (int j = 0; j <= 12; ++j) {
            const double xi = -1.7 + 3.4 * i / 12.0;
            const double s = 2.6 * j / 12.0;
            const double w = xi * xi + s;
            const bool inside = w <= o0c_wall_region().w_cabin;
            const bool outside = w >= o0c_wall_region().w_collar;
            const double x = xi * cfg.scale_R0;
            const double y = cfg.scale_R0 * std::sqrt(s);
            const design::ShiftJet a =
                design::shift_jet(cfg, basis, c, 0.0, x, y, 0.0, true);
            const design::ShiftJet z =
                design::shift_jet(cfg, basis, c, 0.0, x, y, 0.0, false);
            double d0 = std::fabs(a.value - z.value), e1 = 0.0, e2 = 0.0,
                   e3 = 0.0;
            for (int m = 0; m < 4; ++m) {
                e1 = std::fmax(e1, std::fabs(a.d1[m] - z.d1[m]));
                for (int n = 0; n < 4; ++n) {
                    e2 = std::fmax(e2, std::fabs(a.d2[m][n] - z.d2[m][n]));
                    for (int l = 0; l < 4; ++l)
                        e3 = std::fmax(e3,
                                       std::fabs(a.d3[m][n][l] - z.d3[m][n][l]));
                }
            }
            b.d1 = std::fmax(b.d1, e1);
            b.d2 = std::fmax(b.d2, e2);
            b.d3 = std::fmax(b.d3, e3);
            const double worst =
                std::fmax(d0, std::fmax(e1, std::fmax(e2, e3)));
            if (inside) b.interior = std::fmax(b.interior, worst);
            else if (outside) b.collar = std::fmax(b.collar, worst);
        }
    return b;
}

}  // namespace

OptimizerResult run_o1_optimization(const F1wSystem& sys,
                                    const OptimizerConfig& ocfg,
                                    const std::vector<double>& y_start,
                                    int max_accepted_override,
                                    std::string* timing_csv) {
    OptimizerResult R;
    const CertifiedSettings& cs = o0d_settings();
    const ShiftFamilyConfig cfg = o0c_config();
    const int n = sys.whitened;
    const int max_accepted = max_accepted_override > 0 ? max_accepted_override
                                                       : ocfg.max_accepted;

#if !defined(WARP_HYPAD_HAVE_OTILIB)
    (void)cs; (void)cfg; (void)n; (void)max_accepted; (void)y_start;
    (void)timing_csv;
    R.ran = false;
    R.skip_reason =
        "the OTIlib sensitivity backend is not compiled into this "
        "configuration (WARP_HYPAD_WITH_OTILIB=OFF). O1 requires the VERIFIED "
        "batched OTI gradient; substituting a finite-difference production "
        "gradient is explicitly forbidden, so the optimizer does not run.";
    R.termination_reason = "not_run_no_gradient_backend";
    return R;
#else
    using sensitivity::OTIScalar;

    const IntegrationSet wset = o0d_working_set(sys);
    const IntegrationSet vset = o0d_validation_set(sys);

    // The vacuum floors, one per evaluator domain. Fixed for the whole run.
    double vol_w = 0.0, vol_v = 0.0;
    for (double q : wset.weight) vol_w += q;
    for (double q : vset.weight) vol_v += q;
    const double floor_w = cs.tau * kLog2 * vol_w;
    const double floor_v = cs.tau * kLog2 * vol_v;

    std::ostringstream timing;
    timing.precision(6);
    timing << "# wall clock; NON-REPRODUCIBLE by nature and excluded from every "
              "scientific record\n";
    timing << "stage,iteration,seconds\n";

    auto raw_working = [&](const std::vector<double>& c) {
        ++R.objective_evals;
        return design::smooth_negative_energy<double>(
            cfg, sys.parent, c, wset, cs.tau, 0.0,
            SourceSemantics::TotalEffective);
    };
    auto raw_validation = [&](const std::vector<double>& c) {
        return design::smooth_negative_energy<double>(
            cfg, sys.parent, c, vset, cs.tau, 0.0,
            SourceSemantics::TotalEffective);
    };

    // The complete gradient, through twelve deterministic OTI batches in y.
    auto gradient_at = [&](const std::vector<double>& y) {
        ++R.gradient_evals;
        std::vector<double> g(static_cast<std::size_t>(n), 0.0);
        for (int s0 = 0; s0 < n; s0 += kOtiDirections) {
            const int count = std::min(kOtiDirections, n - s0);
            std::vector<OTIScalar> yy(static_cast<std::size_t>(n));
            for (int j = 0; j < n; ++j)
                yy[static_cast<std::size_t>(j)] = sensitivity::make_constant(
                    y[static_cast<std::size_t>(j)]);
            for (int d = 0; d < count; ++d)
                yy[static_cast<std::size_t>(s0 + d)] =
                    sensitivity::seed_parameter(
                        y[static_cast<std::size_t>(s0 + d)], d);
            std::vector<OTIScalar> reduced(static_cast<std::size_t>(sys.active));
            for (int i = 0; i < sys.active; ++i) {
                OTIScalar a = sensitivity::make_constant(
                    sys.projection.c0[static_cast<std::size_t>(i)]);
                for (int j = 0; j < n; ++j)
                    a = a + yy[static_cast<std::size_t>(j)] *
                                OTIScalar(sys.PW(i, j));
                reduced[static_cast<std::size_t>(i)] = a;
            }
            const std::vector<OTIScalar> c = sys.family.prolong<OTIScalar>(
                reduced, sys.parent.coefficient_count());
            const OTIScalar E = design::smooth_negative_energy<OTIScalar>(
                cfg, sys.parent, c, wset, cs.tau, 0.0,
                SourceSemantics::TotalEffective);
            for (int d = 0; d < count; ++d)
                g[static_cast<std::size_t>(s0 + d)] =
                    sensitivity::extract_gradient(E, d);
        }
        return g;
    };

    // --- baseline ------------------------------------------------------------
    const std::vector<double> y_zero(static_cast<std::size_t>(n), 0.0);
    const std::vector<double> c_zero = f1w_coefficients(sys, y_zero);
    R.baseline_raw_working = raw_working(c_zero);
    R.baseline_renorm_working = R.baseline_raw_working - floor_w;
    R.baseline_raw_validation = raw_validation(c_zero);
    R.baseline_renorm_validation = R.baseline_raw_validation - floor_v;
    {
        const ObjectiveValues v = design::evaluate_objectives(
            cfg, sys.parent, c_zero, vset, 0.0, SourceSemantics::TotalEffective,
            false, false);
        R.baseline_E_neg = v.E_neg;
    }

    // --- the first iterate ---------------------------------------------------
    std::vector<double> y = y_start;
    std::vector<double> c = f1w_coefficients(sys, y);
    double f_w = raw_working(c) - floor_w;
    double f_v = raw_validation(c) - floor_v;
    std::vector<double> g = gradient_at(y);

    auto snapshot = [&](int iter, double alpha, double trust, double step,
                        int ls_evals, const std::vector<double>& yy,
                        const std::vector<double>& cc,
                        const std::vector<double>& gg, double fw, double fv) {
        Iterate it;
        it.iteration = iter;
        it.objective_evals = R.objective_evals;
        it.gradient_evals = R.gradient_evals;
        it.renorm_working = fw;
        it.raw_working = fw + floor_w;
        it.delta_working = it.raw_working - R.baseline_raw_working;
        it.renorm_validation = fv;
        it.raw_validation = fv + floor_v;
        it.delta_validation = it.raw_validation - R.baseline_raw_validation;

        const ObjectiveValues ex = design::evaluate_objectives(
            cfg, sys.parent, cc, vset, 0.0, SourceSemantics::TotalEffective,
            false, true);
        it.E_neg = ex.E_neg;
        it.E_plus = ex.E_plus;
        it.E_net = ex.E_net;
        it.E_abs = ex.E_abs;
        it.K_max = ex.K_max;
        it.invalid_points = ex.invalid_points;

        // J_NEC and N_min are the expensive diagnostics: the global null
        // minimisation runs a trust-region solve per point. They are taken on
        // the WORKING set, and the column names say so, rather than being
        // silently reported at a resolution the run cannot afford.
        const ObjectiveValues nec = design::evaluate_objectives(
            cfg, sys.parent, cc, wset, 0.0, SourceSemantics::TotalEffective,
            true, false);
        it.J_NEC = nec.J_NEC;
        it.N_min = nec.N_min_global;
        it.observer_failures = nec.nec_failures;

        const BehaviorLock b = measure_lock(cfg, sys.parent, cc);
        it.max_correction = b.max_correction;
        it.max_d1 = b.d1; it.max_d2 = b.d2; it.max_d3 = b.d3;

        it.grad_norm = nrm(gg);
        it.y_norm = nrm(yy);
        it.step_norm = step;
        it.trust_radius = trust;
        it.alpha = alpha;
        it.line_search_evals = ls_evals;
        it.y = yy;
        it.gradient = gg;
        return it;
    };

    {
        const std::clock_t t0 = std::clock();
        R.accepted.push_back(snapshot(0, 0.0, ocfg.trust_radius_initial, 0.0, 0,
                                      y, c, g, f_w, f_v));
        timing << "iterate,0,"
               << static_cast<double>(std::clock() - t0) / CLOCKS_PER_SEC
               << "\n";
    }

    double trust = ocfg.trust_radius_initial;
    R.trust_radius_min_seen = trust;
    R.trust_radius_max_seen = trust;
    std::deque<std::vector<double>> S, Y;   // L-BFGS history
    int stall = 0;
    R.ran = true;

    for (int iter = 1; iter <= max_accepted; ++iter) {
        if (R.objective_evals >= ocfg.max_objective_evals) {
            R.termination_reason = "objective_evaluation_limit";
            R.terminated_by_limit = true;
            break;
        }
        const double gnorm = nrm(g);
        if (gnorm < ocfg.gradient_norm_tol) {
            R.termination_reason = "certified_gradient_norm_below_threshold";
            R.converged = true;
            break;
        }

        // --- L-BFGS two-loop recursion --------------------------------------
        std::vector<double> p(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) p[static_cast<std::size_t>(i)] = -g[static_cast<std::size_t>(i)];
        bool fell_back = false;
        if (!S.empty()) {
            std::vector<double> q = g;
            std::vector<double> alphas(S.size(), 0.0);
            for (int k = static_cast<int>(S.size()) - 1; k >= 0; --k) {
                const double rho = 1.0 / dot(Y[static_cast<std::size_t>(k)],
                                             S[static_cast<std::size_t>(k)]);
                const double a = rho * dot(S[static_cast<std::size_t>(k)], q);
                alphas[static_cast<std::size_t>(k)] = a;
                for (int i = 0; i < n; ++i)
                    q[static_cast<std::size_t>(i)] -=
                        a * Y[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)];
            }
            const std::size_t last = S.size() - 1;
            const double gamma = dot(S[last], Y[last]) / dot(Y[last], Y[last]);
            for (int i = 0; i < n; ++i) q[static_cast<std::size_t>(i)] *= gamma;
            for (std::size_t k = 0; k < S.size(); ++k) {
                const double rho = 1.0 / dot(Y[k], S[k]);
                const double b = rho * dot(Y[k], q);
                for (int i = 0; i < n; ++i)
                    q[static_cast<std::size_t>(i)] +=
                        (alphas[k] - b) * S[k][static_cast<std::size_t>(i)];
            }
            for (int i = 0; i < n; ++i)
                p[static_cast<std::size_t>(i)] = -q[static_cast<std::size_t>(i)];
        }
        double dd = dot(p, g);
        if (!(dd < 0.0)) {
            // Not a descent direction: discard the incompatible history and
            // fall back to steepest descent, which always is one.
            S.clear(); Y.clear();
            for (int i = 0; i < n; ++i)
                p[static_cast<std::size_t>(i)] = -g[static_cast<std::size_t>(i)];
            dd = dot(p, g);
            fell_back = true;
            ++R.lbfgs_fallbacks;
        }

        // --- Armijo backtracking inside the trust radius --------------------
        const double pn = nrm(p);
        double alpha = pn > 0.0 ? std::fmin(1.0, trust / pn) : 0.0;
        bool accepted = false;
        int attempt = 0;
        std::vector<double> y_new, c_new, g_new;
        double fw_new = 0.0, fv_new = 0.0, step_norm = 0.0;

        for (; attempt < ocfg.max_line_search; ++attempt) {
            y_new = y;
            for (int i = 0; i < n; ++i)
                y_new[static_cast<std::size_t>(i)] +=
                    alpha * p[static_cast<std::size_t>(i)];
            std::vector<double> d(static_cast<std::size_t>(n));
            for (int i = 0; i < n; ++i)
                d[static_cast<std::size_t>(i)] =
                    y_new[static_cast<std::size_t>(i)] - y[static_cast<std::size_t>(i)];
            step_norm = nrm(d);

            RejectedTrial rt;
            rt.iteration = iter; rt.attempt = attempt; rt.alpha = alpha;
            rt.trust_radius = trust; rt.step_norm = step_norm;

            bool finite_y = true;
            for (double q : y_new) if (!std::isfinite(q)) finite_y = false;
            c_new = f1w_coefficients(sys, y_new);
            bool finite_c = true;
            for (double q : c_new) if (!std::isfinite(q)) finite_c = false;
            if (!finite_y || !finite_c) {
                rt.classification = "invalid";
                rt.reason = finite_y ? "non_finite_coefficient"
                                     : "non_finite_coordinate";
                R.rejected.push_back(rt);
                alpha *= ocfg.backtrack; continue;
            }

            const double raw_w = raw_working(c_new);
            fw_new = raw_w - floor_w;
            rt.raw_working = raw_w; rt.renorm_working = fw_new;

            LineSearchEntry ls;
            ls.iteration = iter; ls.attempt = attempt; ls.alpha = alpha;
            ls.f_trial = fw_new; ls.f_current = f_w;
            ls.directional_derivative = dd;
            ls.armijo_target = f_w + ocfg.c1_armijo * alpha * dd;
            ls.armijo_satisfied = std::isfinite(fw_new) &&
                                  fw_new <= ls.armijo_target;
            R.line_search.push_back(ls);

            if (!ls.armijo_satisfied) {
                rt.classification = "armijo_failed";
                rt.reason = "working_renormalized_surrogate_did_not_satisfy_"
                            "sufficient_decrease";
                R.rejected.push_back(rt);
                alpha *= ocfg.backtrack; continue;
            }

            // --- validation ---------------------------------------------------
            const ObjectiveValues ov = design::evaluate_objectives(
                cfg, sys.parent, c_new, vset, 0.0,
                SourceSemantics::TotalEffective, false, false);
            const double raw_v = raw_validation(c_new);
            fv_new = raw_v - floor_v;
            rt.raw_validation = raw_v; rt.renorm_validation = fv_new;
            rt.E_neg = ov.E_neg;

            const BehaviorLock lock = measure_lock(cfg, sys.parent, c_new);
            rt.interior = lock.interior; rt.collar = lock.collar;

            if (lock.interior != 0.0 || lock.collar != 0.0) {
                rt.classification = "behavior_changing";
                rt.reason = "protected_region_correction_nonzero";
                R.rejected.push_back(rt);
                alpha *= ocfg.backtrack; continue;
            }
            if (ov.invalid_points > 0 || !std::isfinite(ov.E_neg) ||
                !std::isfinite(fv_new) || ov.sqrt_gamma_min <= 0.0) {
                rt.classification = "invalid";
                rt.reason = "invalid_points_or_non_finite_validation_objective";
                R.rejected.push_back(rt);
                alpha *= ocfg.backtrack; continue;
            }
            // The difference of differences: the y-independent domain offset
            // between the evaluators cancels, so this measures resolution.
            const double dw = fw_new - f_w;
            const double dv = fv_new - f_v;
            // THE DRIFT CRITERION IS MIXED, AND THE 2 PER CENT IS MEASURED.
            //
            // A purely relative test on the difference of differences fails for
            // the same reason every other one in this project has: dv shrinks
            // with the step, so the ratio measures the evaluators' disagreement
            // floor rather than the candidate. Worse, the working set is
            // Gauss-Legendre order 4 while the validation set is order 8, so
            // even a perfectly ordinary candidate shows a consistent ~0.3 per
            // cent order difference that has nothing to do with resolution of
            // FEATURES -- measured across 240 trials at step norms spanning
            // seven orders of magnitude, the drift sat at 1.8e-3 to 5.4e-3 and
            // did NOT shrink with the step, which is the signature of an order
            // difference rather than an exploit.
            //
            // So: an absolute floor at the certified eps_abs, and a relative
            // band of 2 per cent, which is six times the measured order-4
            // versus order-8 discrepancy on this integrand. A candidate that
            // exploits quadrature nodes shows drift orders of magnitude larger
            // than that -- O0D's 1e6-magnitude probe drifted by 1.5e-2 with a
            // change of the same size, i.e. order 1 relative.
            rt.drift = std::fabs(dw - dv) /
                       (cs.eps_abs + 0.02 * std::fabs(dv));
            if (rt.drift > 1.0) {
                rt.classification = "unresolved";
                rt.reason = "working_validation_change_disagreement";
                R.rejected.push_back(rt);
                alpha *= ocfg.backtrack; continue;
            }
            if (!(fv_new < f_v)) {
                rt.classification = "validation_not_decreasing";
                rt.reason = "validation_renormalized_surrogate_did_not_decrease";
                R.rejected.push_back(rt);
                alpha *= ocfg.backtrack; continue;
            }
            if (ov.E_neg > R.accepted.back().E_neg +
                               o0d_state().numerical_uncertainty) {
                rt.classification = "exact_energy_increase";
                rt.reason = "exact_E_neg_increased_beyond_numerical_uncertainty";
                R.rejected.push_back(rt);
                alpha *= ocfg.backtrack; continue;
            }
            accepted = true;
            break;
        }

        if (!accepted) {
            trust *= ocfg.trust_shrink;
            R.trust_radius_min_seen = std::fmin(R.trust_radius_min_seen, trust);
            if (trust < ocfg.trust_radius_min) {
                R.termination_reason = "trust_radius_below_minimum";
                break;
            }
            continue;
        }

        g_new = gradient_at(y_new);

        // --- L-BFGS curvature pair, accepted only if positive ---------------
        std::vector<double> s_vec(static_cast<std::size_t>(n)),
            y_vec(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            s_vec[static_cast<std::size_t>(i)] =
                y_new[static_cast<std::size_t>(i)] - y[static_cast<std::size_t>(i)];
            y_vec[static_cast<std::size_t>(i)] =
                g_new[static_cast<std::size_t>(i)] - g[static_cast<std::size_t>(i)];
        }
        const double sy = dot(s_vec, y_vec);
        if (sy > 1e-14 * nrm(s_vec) * nrm(y_vec)) {
            S.push_back(s_vec); Y.push_back(y_vec);
            ++R.curvature_pairs_accepted;
            while (static_cast<int>(S.size()) > ocfg.lbfgs_memory) {
                S.pop_front(); Y.pop_front();
            }
        } else {
            ++R.curvature_pairs_rejected;
        }

        const double prev_E_neg = R.accepted.back().E_neg;
        {
            const std::clock_t t0 = std::clock();
            R.accepted.push_back(snapshot(iter, alpha, trust, step_norm,
                                          attempt + 1, y_new, c_new, g_new,
                                          fw_new, fv_new));
            timing << "iterate," << iter << ","
                   << static_cast<double>(std::clock() - t0) / CLOCKS_PER_SEC
                   << "\n";
        }
        (void)fell_back;

        y = y_new; c = c_new; g = g_new; f_w = fw_new; f_v = fv_new;

        // A full step that succeeded first try earns a larger radius.
        if (attempt == 0) {
            trust = std::fmin(trust * ocfg.trust_expand, ocfg.trust_radius_max);
            R.trust_radius_max_seen =
                std::fmax(R.trust_radius_max_seen, trust);
        }

        if (step_norm < ocfg.step_norm_tol) {
            R.termination_reason = "accepted_step_norm_below_threshold";
            R.converged = true;
            break;
        }
        if (std::fabs(R.accepted.back().E_neg - prev_E_neg) <
            o0d_state().numerical_uncertainty) {
            if (++stall >= ocfg.stall_window) {
                R.termination_reason =
                    "exact_objective_change_below_uncertainty_for_" +
                    std::to_string(ocfg.stall_window) +
                    "_consecutive_accepted_iterations";
                R.converged = true;
                break;
            }
        } else {
            stall = 0;
        }
        if (static_cast<int>(R.accepted.size()) - 1 >= max_accepted) {
            R.termination_reason = "accepted_iteration_limit";
            R.terminated_by_limit = true;
            break;
        }
    }

    if (R.termination_reason.empty()) {
        R.termination_reason = "accepted_iteration_limit";
        R.terminated_by_limit = true;
    }
    if (timing_csv) *timing_csv = timing.str();
    return R;
#endif
}

}  // namespace checks
}  // namespace warp_hypad
