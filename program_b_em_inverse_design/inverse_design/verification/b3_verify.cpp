#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>

#include "b3_checks.hpp"
#include "program_b/inverse_design/band_objective.hpp"
#include "warp_hypad/core/metadata.hpp"

// Gate B3-PHYSICAL-ADMISSIBILITY.
//
// SCOPE. A passing result is a passive, causal parametric constitutive model
// that approximates the target response over a stated frequency band. It is NOT
// a fabricated medium or an experimental device. There is no unit-cell geometry,
// no hardware control model, and no Maxwell propagation anywhere in B3.

int main(int argc, char** argv) {
    using namespace warp_hypad;
    namespace id = program_b::inverse_design;

    std::string out_dir = "results/verification/program_b3";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            std::cerr << "usage: b3_verify [--out <dir>]\n";
            return 2;
        }
    }

    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if (ec) {
        std::cerr << "cannot create " << out_dir << ": " << ec.message() << "\n";
        return 3;
    }

    VerificationReport report(
        "program_b3",
        "B3-PHYSICAL-ADMISSIBILITY: frequency-dependent constitutive models. "
        "q_physical -> eps(w), mu^-1(w), zeta(w) -> Z(w) -> J_band, with "
        "structural causality and a derived passivity condition");

    program_b::checks::B3Metrics m;
    program_b::checks::run_b3_checks(report, out_dir, m);

    report.print(std::cout);

    std::string error;
    if (!report.write(out_dir, error)) {
        std::cerr << "failed to write reports: " << error << "\n";
        return 3;
    }

    {
        const id::FrequencyConventions conv;
        JsonValue meta = JsonValue::object();
        meta.set("milestone", JsonValue("B3"));
        meta.set("scope",
                 JsonValue("A passing result is a passive, causal parametric "
                           "constitutive model that approximates the target "
                           "response over a stated frequency band. It is NOT a "
                           "fabricated medium or an experimental device. No "
                           "unit-cell geometry, no hardware control model, no "
                           "Maxwell propagation."));
        meta.set("fourier_kernel", JsonValue(conv.fourier_kernel));
        meta.set("time_derivative", JsonValue(conv.time_derivative));
        meta.set("constitutive_relation", JsonValue(conv.constitutive_relation));
        meta.set("reciprocity", JsonValue(conv.reciprocity));
        meta.set("dissipated_power", JsonValue(conv.dissipated_power));
        meta.set("poynting", JsonValue(conv.poynting));
        meta.set("passive_sign", JsonValue(conv.passive_sign));
        meta.set("causal_half_plane", JsonValue(conv.causal_half_plane));
        meta.set("reality_condition", JsonValue(conv.reality));
        meta.set("model_family", JsonValue(m.model_family));
        meta.set("pole_locations", JsonValue(m.pole_locations));
        meta.set("causality_enforcement",
                 JsonValue("structural: damping is a strictly positive function "
                           "of an unconstrained variable, so every pole is in the "
                           "lower half omega-plane for every parameter value. "
                           "Kramers-Kronig is a diagnostic only."));
        meta.set("hard_feasibility",
                 JsonValue("bounded parameter transforms plus projection. A "
                           "soft-penalty result is never labelled feasible."));
        meta.set("build", build_info_json());
        meta.set("timestamp_utc", JsonValue(utc_timestamp()));
        if (!json_write_file(out_dir + "/b3_metadata.json", meta, error)) return 3;
    }

    std::cout << "\nB3 REQUIRED QUANTITIES\n";
    std::cout << "  Fourier convention                 " << m.fourier_convention << "\n";
    std::cout << "  model family                       " << m.model_family << "\n";
    std::cout << "  pole locations                     " << m.pole_locations << "\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  total / passed / failed / skipped  " << report.num_cases() << " / "
              << report.num_passed() << " / " << report.num_failed() << " / "
              << report.num_skipped() << "\n";
    std::cout << "  min loss eigenvalue (passive case) " << m.min_loss_eigenvalue << "\n";
    std::cout << "  max passivity violation            " << m.max_passivity_violation << "\n";
    std::cout << "  worst passivity frequency          " << m.worst_passivity_frequency << "\n";
    std::cout << "  centre mismatch                    " << m.centre_mismatch << "\n";
    std::cout << "  mean mismatch                      " << m.mean_mismatch << "\n";
    std::cout << "  max mismatch                       " << m.max_mismatch << "\n";
    std::cout << "  worst mismatch frequency           " << m.worst_mismatch_frequency << "\n";
    std::cout << "  hard-mode bound violation          " << m.hard_bound_violation
              << "   (must be exactly 0)\n";
    std::cout << "  soft-mode bound violation          " << m.soft_bound_violation
              << "   (reported, never called feasible)\n";
    std::cout << "  gradient vs analytic               " << m.worst_gradient_analytic << "\n";
    std::cout << "  gradient vs complex step           " << m.worst_gradient_complex_step << "\n";
    if (m.oti_available) {
        std::cout << "  batched HyPAD gradient             " << m.worst_gradient_oti
                  << "  (" << m.oti_batches << " batches)\n";
    } else {
        std::cout << "  batched HyPAD gradient             skipped (backend disabled)\n";
    }
    std::cout << "  one-pole recovery J_band           " << m.one_pole_recovery_error << "\n";
    std::cout << "  two-pole recovery J_band           " << m.two_pole_recovery_error << "\n";
    std::cout << "  restricted one-pole residual       " << m.restricted_residual
              << "   (must be nonzero)\n";
    std::cout << "  worst parameter recovery error     " << m.worst_parameter_recovery_error << "\n";
    std::cout << "  coordinate-change residual         " << m.worst_coordinate_change << "\n";
    std::cout << "  real-limit residual vs B2          " << m.real_limit_residual << "\n";
    std::cout << "  termination reason                 " << m.termination_reason << "\n";

    std::cout << "\nSCOPE: a passing result is a passive, causal parametric "
                 "constitutive model\n"
                 "  approximating the target over a stated band. It is NOT a "
                 "fabricated medium\n"
                 "  or an experimental device.\n";

    if (!report.all_passed()) {
        std::cout << "\nB3 GATE: FAIL (" << report.num_failed() << " of "
                  << (report.num_cases() - report.num_skipped())
                  << " cases failed)\n";
        return 1;
    }

    std::cout << "\nB3 GATE: PASS (" << report.num_passed() << "/"
              << (report.num_cases() - report.num_skipped()) << ")";
    if (report.num_skipped() > 0)
        std::cout << ", " << report.num_skipped()
                  << " skipped (HyPAD gradient cases only)";
    std::cout << "\n";
    return 0;
}
