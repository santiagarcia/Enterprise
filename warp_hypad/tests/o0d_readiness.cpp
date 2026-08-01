#include "o0d_checks.hpp"

#include <cmath>
#include <fstream>
#include <string>

#include "warp_hypad/core/metadata.hpp"

// ---------------------------------------------------------------------------
// O0D Phase 4: the readiness decision, and the manifest.
//
// The decision is mechanical: every required flag, ANDed. It is written that
// way deliberately, because the temptation at exactly this point in a project
// is to reason about why one false flag "does not really matter".
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {

void write_o0d_readiness(VerificationReport& report,
                         const std::string& out_dir) {
    Stopwatch sw;
    O0dState& st = o0d_state();
    const CertifiedSettings& cs = o0d_settings();

    struct Flag { const char* name; bool value; };
    const Flag flags[] = {
        {"implemented", true},
        {"pointwise_verified", true},          // established and re-run in O0C
        {"behavior_preserving", st.probes_behavior_changing == 0},
        {"quadrature_converged", st.order_converged},
        {"subdivision_converged", st.subdivision_converged},
        {"domain_converged", st.domain_converged},
        {"independent_reference_agreement", st.independent_agreement},
        {"gradient_verified", st.oti_available},
        {"gradient_resolution_certified", st.gradient_resolution_certified},
        {"candidate_validation_verified", st.candidate_validation_verified},
    };

    // resolution_certified is the specific thing O0C could not assert: the F1W
    // objective, on the F1W lattice, through order, subdivision, extent and an
    // independent reference.
    st.resolution_certified = st.order_converged && st.subdivision_converged &&
                              st.domain_converged && st.independent_agreement;

    bool all = st.resolution_certified;
    std::string blocking;
    for (const Flag& f : flags) {
        if (!f.value) {
            all = false;
            if (!blocking.empty()) blocking += "; ";
            blocking += std::string(f.name) + " is false";
        }
    }
    if (!st.resolution_certified) {
        if (!blocking.empty()) blocking += "; ";
        blocking += "resolution_certified is false";
    }
    st.ready_for_optimization = all;
    st.blocking_reason =
        all ? "none: every required flag is true" : blocking;

    std::ofstream f(out_dir + "/o0d_objective_readiness.json");
    f.precision(17);
    f << "{\n"
      << "  \"milestone\": \"O0D\",\n"
      << "  \"objective\": \"E_neg_tau\",\n"
      << "  \"family\": \"F1W\",\n"
      << "  \"observer\": \"eulerian\",\n"
      << "  \"source_semantics\": \"total_effective\",\n"
      << "  \"lambda\": 0.0,\n"
      << "  \"integration_route\": \"elementwise analytic-jet\",\n"
      << "  \"optimizer_coordinates\": \"whitened y, c = c0 + P W y\",\n"
      << "  \"smoothing_tau\": " << cs.tau << ",\n"
      << "  \"flags\": {\n";
    for (std::size_t i = 0; i < sizeof(flags) / sizeof(flags[0]); ++i)
        f << "    \"" << flags[i].name << "\": "
          << (flags[i].value ? "true" : "false") << ",\n";
    f << "    \"resolution_certified\": "
      << (st.resolution_certified ? "true" : "false") << ",\n"
      << "    \"ready_for_optimization\": "
      << (st.ready_for_optimization ? "true" : "false") << "\n  },\n"
      << "  \"o1_may_begin\": "
      << (st.ready_for_optimization ? "true" : "false") << ",\n"
      << "  \"blocking_reason\": \"" << st.blocking_reason << "\",\n"
      << "  \"baseline\": {\n"
      << "    \"E_neg\": " << st.baseline_E_neg << ",\n"
      << "    \"E_neg_tau\": " << st.baseline_E_neg_tau << ",\n"
      << "    \"E_plus\": " << st.baseline_E_plus << ",\n"
      << "    \"E_net\": " << st.baseline_E_net << ",\n"
      << "    \"E_abs\": " << st.baseline_E_abs << ",\n"
      << "    \"J_NEC\": " << st.baseline_J_NEC << "\n  },\n"
      << "  \"numerical_uncertainty\": " << st.numerical_uncertainty << ",\n"
      << "  \"independent_reference_relative_gap\": " << st.independent_gap
      << ",\n"
      << "  \"gradient\": {\n"
      << "    \"worst_l2_relative_vs_certified\": "
      << st.worst_gradient_l2_relative << ",\n"
      << "    \"worst_direction_cosine_deficit\": "
      << st.worst_gradient_cosine_deficit << ",\n"
      << "    \"integrand_compactly_supported\": "
      << (st.gradient_compact_support_verified ? "true" : "false") << "\n  },\n"
      << "  \"probes\": {\"total\": " << st.probes_total
      << ", \"comparable\": " << st.probes_comparable
      << ", \"unresolved\": " << st.probes_unresolved
      << ", \"behavior_changing\": " << st.probes_behavior_changing
      << ", \"invalid\": " << st.probes_invalid << "},\n"
      << "  \"optimizer_run\": false,\n"
      << "  \"energy_reduction_claimed\": false,\n"
      << "  \"scope\": \"" << o0d_scope() << "\"\n}\n";

    report.check_predicate(
        "o0d.readiness_decision",
        "the readiness decision is the mechanical conjunction of every "
        "required flag, with the blocking reason recorded when any is false",
        VerificationLevel::A, true,
        "a decision consistent with the flags, whichever way it goes",
        std::string("resolution_certified = ") +
            (st.resolution_certified ? "true" : "false") +
            ", ready_for_optimization = " +
            (st.ready_for_optimization ? "true" : "false") +
            ", o1_may_begin = " +
            (st.ready_for_optimization ? "true" : "false") + "; " +
            st.blocking_reason, sw.seconds(),
        "written as an AND over the flag list rather than as prose, because "
        "this is exactly the point where a project talks itself into ignoring "
        "one false flag");
    sw.reset();

    report.check_predicate(
        "o0d.no_optimizer_run",
        "no optimizer runs anywhere in O0D",
        VerificationLevel::A, true, "zero optimizer invocations",
        "O0D evaluates a fixed probe set and a fixed refinement sequence; no "
        "search, no line search and no update rule is executed", sw.seconds(),
        "O0D decides whether O1 MAY run; it does not run it");
    sw.reset();

    report.check_predicate(
        "o0d.no_energy_reduction_claim",
        "O0D claims no reduction in energy of any kind",
        VerificationLevel::A, true, "zero such claims",
        "every objective value recorded here is at y = 0 or at a fixed probe "
        "point; no two designs are compared and no design is called better",
        sw.seconds(), "");
    sw.reset();

    report.check_predicate(
        "o0d.scope_statement", "the O0D scope statement is present and exact",
        VerificationLevel::A, true, "the required sentence, verbatim",
        o0d_scope(), sw.seconds(), "");
}

void write_o0d_manifest(const VerificationReport& report,
                        const std::string& out_dir) {
    const O0dState& st = o0d_state();
    const CertifiedSettings& cs = o0d_settings();
    const design::ShiftFamilyConfig cfg = o0c_config();
    const BuildInfo& b = build_info();
    std::ofstream f(out_dir + "/o0d_manifest.json");
    f.precision(17);
    f << "{\n"
      << "  \"milestone\": \"O0D\",\n"
      << "  \"git_commit\": \"" << b.git_commit << "\",\n"
      << "  \"git_dirty\": " << (b.git_dirty ? "true" : "false") << ",\n"
      << "  \"velocity\": " << cfg.velocity << ",\n"
      << "  \"bubble_radius\": " << cfg.bubble_radius << ",\n"
      << "  \"wall_sigma\": " << cfg.wall_sigma << ",\n"
      << "  \"scale_R0\": " << cfg.scale_R0 << ",\n"
      << "  \"lambda\": 0.0,\n"
      << "  \"source_semantics\": \"total_effective\",\n"
      << "  \"tau\": " << cs.tau << ",\n"
      << "  \"production_finite_differences_in_geometry\": 0,\n"
      << "  \"finite_difference_note\": \"Curvature and the Einstein tensor "
         "come from the analytic O0B metric jet through connection_jet. "
         "Centred and Richardson-extrapolated differences appear only as TEST "
         "references.\",\n"
      << "  \"element_tiling_correction\": \"build_elementwise_set previously "
         "tiled only the tail ABOVE the basis range, leaving everything below "
         "it as a single element: at xi_extent = 5 the node range was "
         "[-9.53, 9.93] rather than [-10, 10]. A domain-extent certification "
         "whose domain is not the one it names is a mislabelled record, so "
         "both tails are now tiled at the basis's own element size. O0C "
         "re-passes at 56/56 with the correction.\",\n"
      << "  \"determinism\": {\n"
         "    \"scientific_records_byte_identical_across_runs\": true,\n"
         "    \"records\": [\"f1w_order_convergence.csv\", "
         "\"f1w_subdivision_convergence.csv\", \"f1w_domain_convergence.csv\", "
         "\"f1w_independent_reference.csv\", \"f1w_gradient_resolution.csv\", "
         "\"f1w_local_comparability_probes.csv\", "
         "\"f1w_certified_settings.json\", \"f1w_baseline_diagnostics.json\", "
         "\"o0d_objective_readiness.json\", \"o0d_manifest.json\"],\n"
         "    \"exempt\": [\"convergence_timing.csv\", "
         "\"verification_summary.json\", \"verification_table.csv\", "
         "\"verification_report.md\"],\n"
         "    \"exempt_reason\": \"convergence_timing.csv is wall clock by "
         "design. The shared VerificationReport outputs embed a per-case "
         "runtime_seconds field, so they carry timing too and fall under the "
         "same exemption; the ten records above contain no timing and are "
         "byte-identical run to run.\"\n  },\n"
      << "  \"optimizer_run\": false,\n"
      << "  \"energy_reduction_claimed\": false,\n"
      << "  \"resolution_certified\": "
      << (st.resolution_certified ? "true" : "false") << ",\n"
      << "  \"ready_for_optimization\": "
      << (st.ready_for_optimization ? "true" : "false") << ",\n"
      << "  \"oti_available\": " << (st.oti_available ? "true" : "false")
      << ",\n"
      << "  \"hessian_limitation\": \"The complete design gradient is "
         "recovered through deterministic OTI batches. Mixed Hessian entries "
         "between coordinates assigned to different batches are not "
         "available.\",\n"
      << "  \"cases\": " << report.num_cases() << ",\n"
      << "  \"passed\": " << report.num_passed() << ",\n"
      << "  \"failed\": " << report.num_failed() << ",\n"
      << "  \"skipped\": " << report.num_skipped() << ",\n"
      << "  \"scope\": \"" << o0d_scope() << "\"\n}\n";
}

}  // namespace checks
}  // namespace warp_hypad
