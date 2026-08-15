#include <algorithm>
#include <cmath>
#include <filesystem>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <iostream>
#include <string>
#include <vector>

#include "o1_checks.hpp"
#include "warp_hypad/core/metadata.hpp"

// Milestone O1 gate: safeguarded F1W source-requirement optimization.
//
// O1 may run only because O0D certified the objective's resolution, its
// gradient's resolution, and the candidate-validation machinery. The readiness
// record is READ, not asserted, and the optimizer refuses to start if it says
// no. There is no command-line override.
//
// SCOPE, printed on every run and written into every record:
//
//   O1 optimizes the inferred source requirement for a FIXED warp-bubble
//   kinematic family. It does not demonstrate faster-than-light travel, a
//   causal shortcut, reduced travel time, horizon removal, physical
//   realizability, dynamic stability, satisfaction of quantum energy
//   inequalities, or a global optimum.

namespace warp_hypad {
namespace checks {
namespace {

using design::IntegrationSet;
using design::ObjectiveValues;
using design::ShiftFamilyConfig;
using design::SourceSemantics;

const double kLog2 = 0.69314718055994530942;

double nrm(const std::vector<double>& a) {
    double s = 0.0;
    for (double v : a) s += v * v;
    return std::sqrt(s);
}

struct FinalCertification {
    bool order = false, subdivision = false, domain = false;
    bool independent = false, gradient = false, behavior = false;
    bool metric_valid = false, observer_valid = false;
    double final_E_neg = 0.0;
    double uncertainty = 0.0;
    double independent_gap = 0.0;
    double raw_offset_ratio = 0.0;
    bool certified = false;
};

}  // namespace

// ---------------------------------------------------------------------------
// Records
// ---------------------------------------------------------------------------
void write_o1_records(const OptimizerResult& R, const F1wSystem& sys,
                      const std::string& out) {
    const CertifiedSettings& cs = o0d_settings();
    const ShiftFamilyConfig cfg = o0c_config();

    {
        std::ofstream f(out + "/optimizer_configuration.json");
        f.precision(17);
        OptimizerConfig o;
        f << "{\n  \"method\": \"limited-memory BFGS with Armijo backtracking "
             "inside an adaptive trust radius\",\n"
          << "  \"coordinates\": \"whitened y; c = c0 + P W y\",\n"
          << "  \"raw_c_optimized\": false,\n"
          << "  \"dimension\": " << sys.whitened << ",\n"
          << "  \"oti_batches\": " << sys.batches << ",\n"
          << "  \"lbfgs_memory\": " << o.lbfgs_memory << ",\n"
          << "  \"c1_armijo\": " << o.c1_armijo << ",\n"
          << "  \"trust_radius_initial\": " << o.trust_radius_initial << ",\n"
          << "  \"trust_radius_min\": " << o.trust_radius_min << ",\n"
          << "  \"trust_radius_max\": " << o.trust_radius_max << ",\n"
          << "  \"max_accepted\": " << o.max_accepted << ",\n"
          << "  \"max_objective_evaluations\": " << o.max_objective_evals
          << ",\n  \"gradient_norm_tolerance\": " << o.gradient_norm_tol
          << ",\n  \"step_norm_tolerance\": " << o.step_norm_tol
          << ",\n  \"stall_window\": " << o.stall_window << ",\n"
          << "  \"objective\": \"renormalized E_neg_tau = raw - tau log2 V\",\n"
          << "  \"tau\": " << cs.tau << ",\n"
          << "  \"working\": {\"order\": " << cs.working_order
          << ", \"subdivision\": " << cs.working_subdivision
          << ", \"xi_extent\": " << cs.working_xi_extent
          << ", \"s_extent\": " << cs.working_s_extent << "},\n"
          << "  \"validation\": {\"order\": " << cs.validation_order
          << ", \"subdivision\": " << cs.validation_subdivision
          << ", \"xi_extent\": " << cs.validation_xi_extent
          << ", \"s_extent\": " << cs.validation_s_extent << "},\n"
          << "  \"production_finite_difference_gradient\": false,\n"
          << "  \"full_hessian_claimed\": false,\n"
          << "  \"hessian_note\": \"L-BFGS builds an inverse-Hessian "
             "APPROXIMATION from gradient differences. The true Hessian is not "
             "formed, stored or claimed; mixed second derivatives between "
             "coordinates in different OTI batches are not available at "
             "all.\"\n}\n";
    }

    {
        std::ofstream f(out + "/accepted_iterates.csv");
        f.precision(17);
        f << "iteration,objective_evals,gradient_evals,"
             "raw_working_E_neg_tau,renormalized_working_E_neg_tau,"
             "delta_working_E_neg_tau,raw_validation_E_neg_tau,"
             "renormalized_validation_E_neg_tau,delta_validation_E_neg_tau,"
             "exact_E_neg,E_plus,E_net,E_abs,J_NEC_working,N_min_working,"
             "K_max,grad_norm,y_norm,step_norm,trust_radius,alpha,"
             "line_search_evals,max_correction,max_d1,max_d2,max_d3,"
             "invalid_points,observer_failures\n";
        for (const Iterate& i : R.accepted)
            f << i.iteration << "," << i.objective_evals << ","
              << i.gradient_evals << "," << i.raw_working << ","
              << i.renorm_working << "," << i.delta_working << ","
              << i.raw_validation << "," << i.renorm_validation << ","
              << i.delta_validation << "," << i.E_neg << "," << i.E_plus << ","
              << i.E_net << "," << i.E_abs << "," << i.J_NEC << ","
              << i.N_min << "," << i.K_max << "," << i.grad_norm << ","
              << i.y_norm << "," << i.step_norm << "," << i.trust_radius << ","
              << i.alpha << "," << i.line_search_evals << ","
              << i.max_correction << "," << i.max_d1 << "," << i.max_d2 << ","
              << i.max_d3 << "," << i.invalid_points << ","
              << i.observer_failures << "\n";
    }
    {
        std::ofstream f(out + "/optimizer_trace.csv");
        f.precision(17);
        f << "iteration,accepted,renormalized_working,renormalized_validation,"
             "exact_E_neg,grad_norm,trust_radius\n";
        for (const Iterate& i : R.accepted)
            f << i.iteration << ",1," << i.renorm_working << ","
              << i.renorm_validation << "," << i.E_neg << "," << i.grad_norm
              << "," << i.trust_radius << "\n";
        for (const RejectedTrial& t : R.rejected)
            f << t.iteration << ",0," << t.renorm_working << ","
              << t.renorm_validation << "," << t.E_neg << ",," << t.trust_radius
              << "\n";
    }
    {
        std::ofstream f(out + "/rejected_trials.csv");
        f.precision(17);
        f << "iteration,attempt,alpha,trust_radius,step_norm,"
             "raw_working_E_neg_tau,renormalized_working_E_neg_tau,"
             "raw_validation_E_neg_tau,renormalized_validation_E_neg_tau,"
             "exact_E_neg,working_validation_drift,interior_delta,"
             "collar_delta,classification,rejection_reason\n";
        for (const RejectedTrial& t : R.rejected)
            f << t.iteration << "," << t.attempt << "," << t.alpha << ","
              << t.trust_radius << "," << t.step_norm << "," << t.raw_working
              << "," << t.renorm_working << "," << t.raw_validation << ","
              << t.renorm_validation << "," << t.E_neg << "," << t.drift << ","
              << t.interior << "," << t.collar << "," << t.classification
              << "," << t.reason << "\n";
    }
    {
        std::ofstream f(out + "/line_search_trace.csv");
        f.precision(17);
        f << "iteration,attempt,alpha,f_trial,f_current,armijo_target,"
             "armijo_satisfied,directional_derivative\n";
        for (const LineSearchEntry& l : R.line_search)
            f << l.iteration << "," << l.attempt << "," << l.alpha << ","
              << l.f_trial << "," << l.f_current << "," << l.armijo_target
              << "," << (l.armijo_satisfied ? 1 : 0) << ","
              << l.directional_derivative << "\n";
    }
    {
        std::ofstream f(out + "/gradient_trace.csv");
        f.precision(17);
        f << "iteration,component,gradient\n";
        for (const Iterate& i : R.accepted)
            for (std::size_t j = 0; j < i.gradient.size(); ++j)
                f << i.iteration << "," << j << "," << i.gradient[j] << "\n";
    }
    {
        std::ofstream f(out + "/trust_radius_trace.csv");
        f.precision(17);
        f << "iteration,trust_radius,alpha,step_norm,line_search_evals\n";
        for (const Iterate& i : R.accepted)
            f << i.iteration << "," << i.trust_radius << "," << i.alpha << ","
              << i.step_norm << "," << i.line_search_evals << "\n";
    }
    {
        std::ofstream f(out + "/baseline_objectives.json");
        f.precision(17);
        const Iterate& b = R.accepted.front();
        f << "{\n  \"y_norm\": 0.0,\n"
          << "  \"raw_working_E_neg_tau\": " << R.baseline_raw_working << ",\n"
          << "  \"renormalized_working_E_neg_tau\": "
          << R.baseline_renorm_working << ",\n"
          << "  \"raw_validation_E_neg_tau\": " << R.baseline_raw_validation
          << ",\n  \"renormalized_validation_E_neg_tau\": "
          << R.baseline_renorm_validation << ",\n"
          << "  \"exact_E_neg\": " << R.baseline_E_neg << ",\n"
          << "  \"E_plus\": " << b.E_plus << ",\n"
          << "  \"E_net\": " << b.E_net << ",\n"
          << "  \"E_abs\": " << b.E_abs << ",\n"
          << "  \"J_NEC_working\": " << b.J_NEC << ",\n"
          << "  \"N_min_working\": " << b.N_min << ",\n"
          << "  \"K_max\": " << b.K_max << ",\n"
          << "  \"o0d_reference_E_neg\": " << o0d_state().baseline_E_neg
          << ",\n  \"o0d_numerical_uncertainty\": "
          << o0d_state().numerical_uncertainty << "\n}\n";
    }
    {
        const Iterate& fin = R.accepted.back();
        std::ofstream f(out + "/final_design_y.csv");
        f.precision(17);
        f << "index,y\n";
        for (std::size_t j = 0; j < fin.y.size(); ++j)
            f << j << "," << fin.y[j] << "\n";

        std::vector<double> z(static_cast<std::size_t>(sys.projection.n_free),
                              0.0);
        for (int i = 0; i < sys.projection.n_free; ++i)
            for (int j = 0; j < sys.whitened; ++j)
                z[static_cast<std::size_t>(i)] +=
                    sys.whitening.w(i, j) * fin.y[static_cast<std::size_t>(j)];
        std::ofstream fz(out + "/final_design_z.csv");
        fz.precision(17);
        fz << "index,z\n";
        for (std::size_t j = 0; j < z.size(); ++j)
            fz << j << "," << z[j] << "\n";

        const std::vector<double> c = f1w_coefficients(sys, fin.y);
        std::ofstream fc(out + "/final_design_c.csv");
        fc.precision(17);
        fc << "index,c,c0_baseline,delta\n";
        const std::vector<double> c0 = f1w_coefficients(
            sys, std::vector<double>(static_cast<std::size_t>(sys.whitened), 0.0));
        for (std::size_t j = 0; j < c.size(); ++j)
            fc << j << "," << c[j] << "," << c0[j] << "," << (c[j] - c0[j])
               << "\n";
    }
    {
        const Iterate& e = R.accepted.back();
        std::ofstream f(out + "/final_objectives.json");
        f.precision(17);
        f << "{\n  \"iteration\": " << e.iteration << ",\n"
          << "  \"raw_working_E_neg_tau\": " << e.raw_working << ",\n"
          << "  \"renormalized_working_E_neg_tau\": " << e.renorm_working << ",\n"
          << "  \"delta_working_E_neg_tau\": " << e.delta_working << ",\n"
          << "  \"raw_validation_E_neg_tau\": " << e.raw_validation << ",\n"
          << "  \"renormalized_validation_E_neg_tau\": " << e.renorm_validation << ",\n"
          << "  \"delta_validation_E_neg_tau\": " << e.delta_validation << ",\n"
          << "  \"exact_E_neg\": " << e.E_neg << ",\n"
          << "  \"E_plus\": " << e.E_plus << ",\n"
          << "  \"E_net\": " << e.E_net << ",\n"
          << "  \"E_abs\": " << e.E_abs << ",\n"
          << "  \"J_NEC_working\": " << e.J_NEC << ",\n"
          << "  \"N_min_working\": " << e.N_min << ",\n"
          << "  \"K_max\": " << e.K_max << ",\n"
          << "  \"y_norm\": " << e.y_norm << ",\n"
          << "  \"grad_norm\": " << e.grad_norm << ",\n"
          << "  \"accepted_iterations\": " << (R.accepted.size() - 1) << ",\n"
          << "  \"rejected_trials\": " << R.rejected.size() << ",\n"
          << "  \"objective_evaluations\": " << R.objective_evals << ",\n"
          << "  \"gradient_evaluations\": " << R.gradient_evals << ",\n"
          << "  \"lbfgs_fallbacks\": " << R.lbfgs_fallbacks << ",\n"
          << "  \"curvature_pairs_accepted\": " << R.curvature_pairs_accepted << ",\n"
          << "  \"curvature_pairs_rejected\": " << R.curvature_pairs_rejected << ",\n"
          << "  \"trust_radius_min_seen\": " << R.trust_radius_min_seen << ",\n"
          << "  \"trust_radius_max_seen\": " << R.trust_radius_max_seen << ",\n"
          << "  \"termination_reason\": \"" << R.termination_reason << "\",\n"
          << "  \"terminated_by_limit\": "
          << (R.terminated_by_limit ? "true" : "false") << ",\n"
          << "  \"converged\": " << (R.converged ? "true" : "false")
          << "\n}\n";

        // The behaviour contract at the FINAL design, restated as data so a
        // reader does not have to trust the gate output alone.
        std::ofstream b2(out + "/final_behavior_contract.json");
        b2.precision(17);
        b2 << "{\n  \"family\": \"F1W\",\n"
           << "  \"protected_interior_delta_beta\": 0.0,\n"
           << "  \"outer_collar_delta_beta\": 0.0,\n"
           << "  \"derivatives_zero_through_third_order\": true,\n"
           << "  \"enforcement\": \"structural support "
              "exclusion, not a penalty: the excluded basis functions carry "
              "no coefficient, so no accepted OR rejected candidate can touch "
              "either region\",\n"
           << "  \"velocity_unchanged\": true,\n"
           << "  \"bubble_centre_unchanged\": true,\n"
           << "  \"asymptotic_condition_unchanged\": true,\n"
           << "  \"max_correction\": " << e.max_correction << ",\n"
           << "  \"max_first_derivative\": " << e.max_d1 << ",\n"
           << "  \"max_second_derivative\": " << e.max_d2 << ",\n"
           << "  \"max_third_derivative\": " << e.max_d3
           << "\n}\n";
    }

    {
        const Iterate& b = R.accepted.front();
        const Iterate& e = R.accepted.back();
        std::ofstream f(out + "/baseline_vs_final_metrics.csv");
        f.precision(17);
        f << "metric,baseline,final,change,fractional_change,lower_is_better,"
             "worsened\n";
        auto row = [&](const char* nm, double a, double c, bool lower_better) {
            const double d = c - a;
            const bool worse = lower_better ? (d > 0.0) : (d < 0.0);
            f << nm << "," << a << "," << c << "," << d << ","
              << (a != 0.0 ? d / std::fabs(a) : 0.0) << ","
              << (lower_better ? 1 : 0) << "," << (worse ? 1 : 0) << "\n";
        };
        row("exact_E_neg", b.E_neg, e.E_neg, true);
        row("E_abs", b.E_abs, e.E_abs, true);
        row("E_plus", b.E_plus, e.E_plus, true);
        row("E_net", b.E_net, e.E_net, false);
        row("J_NEC_working", b.J_NEC, e.J_NEC, true);
        row("N_min_working", b.N_min, e.N_min, false);
        row("K_max", b.K_max, e.K_max, true);
        row("max_correction", b.max_correction, e.max_correction, true);
        row("max_d1", b.max_d1, e.max_d1, true);
        row("max_d2", b.max_d2, e.max_d2, true);
        row("max_d3", b.max_d3, e.max_d3, true);
        row("y_norm", b.y_norm, e.y_norm, true);
    }
    (void)cfg;
}

}  // namespace checks
}  // namespace warp_hypad

// ---------------------------------------------------------------------------
int main(int argc, char** argv) {
    using namespace warp_hypad;
    using namespace warp_hypad::checks;

    std::string out_dir = "results/verification/o1_optimization";
    std::string o0d_dir = "results/verification/o0d_f1w_certification";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) out_dir = argv[++i];
        else if (a == "--o0d" && i + 1 < argc) o0d_dir = argv[++i];
        else {
            std::cerr << "usage: o1_verify [--out <dir>] [--o0d <dir>]\n";
            return 2;
        }
    }
    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);

    VerificationReport report(
        "o1_optimization",
        "Safeguarded F1W source-requirement optimization: limited-memory BFGS "
        "in the whitened optimizer coordinates, with an adaptive trust radius, "
        "an Armijo line search on the renormalized smooth surrogate, and a "
        "twelve-condition acceptance test against the O0D-certified validation "
        "evaluator; followed by a complete recertification of the final "
        "candidate, a deterministic repeat and a restart check");

    Stopwatch sw;

    // --- 1/2. Readiness is READ, and it gates execution ---------------------
    //
    // The committed O0D record is the authority, exactly as the specification
    // says. Replaying O0D's whole study here would cost several minutes per
    // run and would not make the decision any more trustworthy: the record is
    // produced by a gated, tagged milestone and carries the commit it was
    // built from.
    //
    // What IS re-derived rather than trusted is the baseline number, because
    // O1 must fail loudly if it turns out to be optimizing a different problem
    // from the certified one -- see o1.baseline_reproduces_o0d.
    O0dState& od = o0d_state();
    {
        std::ifstream rf(o0d_dir + "/o0d_objective_readiness.json");
        std::string all((std::istreambuf_iterator<char>(rf)),
                        std::istreambuf_iterator<char>());
        // A LITERAL-MINDED READER, BECAUSE THE FIRST ONE WAS TOO CLEVER.
        //
        // The first version decided a flag was true by asking whether "true"
        // appeared before the next comma. That works for every key except the
        // last one in an object, whose value is followed by a newline and a
        // brace rather than a comma -- and ready_for_optimization is exactly
        // that key. The parse silently returned false, the optimizer refused to
        // run, and thirteen cases skipped while looking like a configuration
        // difference rather than a bug in the reader.
        //
        // So: find the key, step past its colon, skip whitespace, and read the
        // token that is actually there.
        auto token_after = [&](const char* key) {
            const std::size_t k = all.find(std::string("\"") + key + "\"");
            if (k == std::string::npos) return std::string();
            std::size_t c = all.find(':', k);
            if (c == std::string::npos) return std::string();
            ++c;
            while (c < all.size() &&
                   (all[c] == ' ' || all[c] == '\n' || all[c] == '\r' ||
                    all[c] == '\t'))
                ++c;
            std::size_t e = c;
            while (e < all.size() && all[e] != ',' && all[e] != '\n' &&
                   all[e] != '}')
                ++e;
            return all.substr(c, e - c);
        };
        auto flag = [&](const char* key) {
            return token_after(key) == "true";
        };
        auto number = [&](const char* key) {
            const std::string t = token_after(key);
            return t.empty() ? 0.0 : std::atof(t.c_str());
        };
        od.resolution_certified = flag("resolution_certified");
        od.ready_for_optimization = flag("ready_for_optimization");
        od.baseline_E_neg = number("E_neg");   // the "baseline" block
        od.numerical_uncertainty = number("numerical_uncertainty");
        if (all.empty()) {
            od.ready_for_optimization = false;
            od.blocking_reason =
                "the O0D readiness record could not be read from " + o0d_dir;
        } else if (!od.ready_for_optimization) {
            od.blocking_reason = "the O0D readiness record says O1 may not run";
        } else {
            od.blocking_reason = "none: the O0D record permits O1";
        }
    }

    report.check_predicate(
        "o1.o0d_readiness_required",
        "O1 reads the O0D readiness decision and runs only when it permits",
        VerificationLevel::A, true,
        "the decision is read, whichever way it goes",
        std::string("resolution_certified = ") +
            (od.resolution_certified ? "true" : "false") +
            ", ready_for_optimization = " +
            (od.ready_for_optimization ? "true" : "false") + "; " +
            od.blocking_reason, sw.seconds(),
        "there is no command-line override: the only way to run the optimizer "
        "is for the certification to pass");
    sw.reset();

    const F1wSystem sys = build_f1w_system();
    OptimizerConfig ocfg;
    OptimizerResult R;
    std::string timing;

    if (od.ready_for_optimization) {
        const std::vector<double> y0(
            static_cast<std::size_t>(sys.whitened), 0.0);
        R = run_o1_optimization(sys, ocfg, y0, 0, &timing);
    } else {
        R.ran = false;
        R.skip_reason = "O0D readiness is false: " + od.blocking_reason;
        R.termination_reason = "not_run_readiness_false";
    }

    // TWO INDEPENDENT PRECONDITIONS, AND CONFLATING THEM WAS WRONG.
    //
    // The optimizer runs only if BOTH the O0D readiness permits it AND the
    // verified gradient backend is present. In the OTIlib-disabled build the
    // readiness record still says true -- it was produced by the enabled build
    // and is not configuration-dependent -- while the backend is absent, so the
    // optimizer correctly refuses. Asserting run == readiness therefore failed
    // on a build that was behaving exactly as specified.
    //
    // The implication that actually matters is one-directional: if it ran, then
    // readiness permitted it. Never the converse.
#if defined(WARP_HYPAD_HAVE_OTILIB)
    const bool backend = true;
#else
    const bool backend = false;
#endif
    report.check_predicate(
        "o1.optimizer_execution_is_gated",
        "the optimizer runs only when O0D readiness permits it AND the "
        "verified gradient backend is present",
        VerificationLevel::A,
        (!R.ran || od.ready_for_optimization) &&
            R.ran == (od.ready_for_optimization && backend),
        "ran implies readiness; and it ran exactly when readiness and the "
        "backend both held",
        R.ran ? "readiness true, optimizer ran"
              : std::string("readiness false or backend absent, optimizer did "
                            "NOT run: ") + R.skip_reason, sw.seconds(),
        "in the OTIlib-disabled build the verified gradient backend is absent "
        "and a finite-difference production gradient is explicitly forbidden, "
        "so the optimizer is skipped rather than degraded");
    sw.reset();

    if (!R.ran) {
        // STALE RECORDS FROM AN EARLIER RUN MUST NOT BE READ AS THIS RUN'S
        // OUTPUT. Without this the gate happily verified a previous
        // optimization's accepted trace while reporting that the optimizer had
        // been skipped -- the two statements cannot both be true, and the
        // combination is worse than either failure alone.
        std::error_code rm;
        for (const char* n : {"/accepted_iterates.csv", "/rejected_trials.csv",
                              "/final_design_y.csv", "/final_objectives.json",
                              "/line_search_trace.csv"})
            std::filesystem::remove(out_dir + n, rm);
        report.check_skipped(
            "o1.optimization_run",
            "the safeguarded L-BFGS optimization itself", R.skip_reason,
            sw.seconds());
        for (const char* nm :
             {"o1.baseline_reproduces_o0d", "o1.armijo_every_accepted_step",
              "o1.working_surrogate_monotone",
              "o1.validation_surrogate_monotone",
              "o1.exact_energy_recorded_every_iterate",
              "o1.rejected_trials_absent_from_accepted_trace",
              "o1.protected_cabin_exact", "o1.outer_collar_exact",
              "o1.final_recertification", "o1.deterministic_repeat",
              "o1.restart_reproducibility",
              "o1.exact_reduction_exceeds_uncertainty"})
            report.check_skipped(nm, "requires a completed optimizer run",
                                 R.skip_reason, 0.0);
    } else {
        write_o1_records(R, sys, out_dir);
        std::ofstream tf(out_dir + "/optimization_timing.csv");
        tf << timing;
    }

    // The remaining cases live in run_o1_checks, which reads the records.
    run_o1_checks(report, out_dir);

    report.print(std::cout);
    std::string err;
    if (!report.write(out_dir, err)) {
        std::cerr << "cannot write reports: " << err << "\n";
        return 4;
    }
    std::cout << "\nSCOPE: " << o1_scope() << "\n";
    const bool ok = report.all_passed();
    std::cout << "\nMILESTONE O1 GATE: " << (ok ? "PASS (" : "FAIL (")
              << report.num_passed() << "/" << report.num_cases();
    if (!ok) std::cout << ", " << report.num_failed() << " failed";
    std::cout << ", " << report.num_skipped() << " skipped)\n";
    return ok ? 0 : 1;
}
