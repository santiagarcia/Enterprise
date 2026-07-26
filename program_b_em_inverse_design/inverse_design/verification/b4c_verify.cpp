#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>

#include "b4c_checks.hpp"
#include "warp_hypad/core/metadata.hpp"

// Gate B4C-PROPAGATION-CORRECTION.
//
// Three corrections to B4: an overstated metricity claim, a solver that could
// not be differentiated because it was hard-wired to double, and a root finder
// blind to roots of even multiplicity.

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/program_b4c";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            std::cerr << "usage: b4c_verify [--out <dir>]\n";
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
        "program_b4c",
        "B4C-PROPAGATION-CORRECTION: the metricity claim withdrawn and replaced, "
        "propagation-invisible constitutive pieces demonstrated, the premetric "
        "classification reported piece by piece, the layered solver made "
        "scalar-generic and differentiated through, and a root finder that sees "
        "repeated and complex roots");

    program_b::checks::B4CMetrics m;
    program_b::checks::run_b4c_checks(report, out_dir, m);

    report.print(std::cout);

    std::string error;
    if (!report.write(out_dir, error)) {
        std::cerr << "failed to write reports: " << error << "\n";
        return 3;
    }

    {
        JsonValue meta = JsonValue::object();
        meta.set("milestone", JsonValue("B4C"));
        meta.set("withdrawn_claim", JsonValue(m.old_metricity_claim));
        meta.set("corrected_claim", JsonValue(m.corrected_metricity_claim));
        meta.set("fresnel_does_not_determine_Z",
                 JsonValue("W^{bd} = Z^{abcd}k_a k_c annihilates the axion piece "
                           "identically, and the cone is conformally invariant, so "
                           "the Fresnel polynomial does NOT uniquely reconstruct "
                           "Z^{abcd}."));
        meta.set("classification_flags",
                 JsonValue("fresnel_single_cone, principal_part_metric_compatible, "
                           "contains_axion, contains_skewon, reciprocal, "
                           "pair_antisymmetric -- reported separately, replacing "
                           "B4's single metric_induced boolean"));
        meta.set("root_finder",
                 JsonValue("Durand-Kerner with Newton polish, clustering for "
                           "multiplicity, per-root residual and conditioning "
                           "warning. Replaces a sign-change scan that was blind "
                           "to roots of even multiplicity."));
        meta.set("differentiation_exclusions", JsonValue(m.exclusions));
        meta.set("build", build_info_json());
        meta.set("timestamp_utc", JsonValue(utc_timestamp()));
        if (!json_write_file(out_dir + "/b4c_metadata.json", meta, error)) return 3;
    }

    std::cout << "\nB4C REQUIRED QUANTITIES\n";
    std::cout << "  WITHDRAWN claim                    " << m.old_metricity_claim << "\n";
    std::cout << "  CORRECTED claim                    " << m.corrected_metricity_claim << "\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  total / passed / failed / skipped  " << report.num_cases() << " / "
              << report.num_passed() << " / " << report.num_failed() << " / "
              << report.num_skipped() << "\n";
    std::cout << "  Fresnel equality under axion       " << m.fresnel_equality_under_axion
              << "   (must be ~0)\n";
    std::cout << "  constitutive difference            " << m.constitutive_difference_under_axion
              << "   (must be nonzero)\n";
    std::cout << "  zeta block difference              " << m.zeta_block_difference_under_axion << "\n";
    std::cout << "  axion coefficient recovered        " << m.axion_coefficient_recovered << "\n";
    std::cout << "  conformal Z difference             " << m.conformal_fresnel_zero_set_shift
              << "   (must be ~0)\n";
    std::cout << "  worst repeated-root residual       " << m.worst_repeated_root_residual << "\n";
    std::cout << "  metric double-root separation      " << m.metric_double_root_separation << "\n";
    std::cout << "  anisotropic relative split         " << m.anisotropic_split << "\n";
    std::cout << "  coefficient reconstruction error   " << m.coefficient_reconstruction_error << "\n";
    std::cout << "  worst stack gradient vs analytic   " << m.worst_stack_gradient_analytic << "\n";
    std::cout << "  worst stack gradient FD spread     " << m.worst_stack_gradient_fd << "\n";
    std::cout << "  magnetoelectric sensitivity        " << m.magnetoelectric_sensitivity
              << " (exact), " << m.magnetoelectric_sensitivity_restricted
              << " (restricted)\n";
    std::cout << "  energy-balance regression          " << m.energy_balance_regression << "\n";
    std::cout << "\n  EXCLUSIONS: " << m.exclusions << "\n";

    if (!report.all_passed()) {
        std::cout << "\nB4C GATE: FAIL (" << report.num_failed() << " of "
                  << (report.num_cases() - report.num_skipped())
                  << " cases failed)\n";
        return 1;
    }

    std::cout << "\nB4C GATE: PASS (" << report.num_passed() << "/"
              << (report.num_cases() - report.num_skipped()) << ")";
    if (report.num_skipped() > 0)
        std::cout << ", " << report.num_skipped()
                  << " skipped (HyPAD comparisons only)";
    std::cout << "\n";
    return 0;
}
