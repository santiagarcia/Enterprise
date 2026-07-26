#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>

#include "b2_checks.hpp"
#include "program_b/inverse_design/design_objective.hpp"
#include "program_b/inverse_design/optimizer.hpp"
#include "warp_hypad/core/metadata.hpp"

#ifdef WARP_HYPAD_HAVE_OTILIB
#include "program_b/inverse_design/design_gradient_oti.hpp"
#endif

// Gate B2-ABSTRACT-INVERSE.
//
// The first inverse-design layer, entirely in an abstract constitutive-design
// space. Depends on B1C through CTest: the reconstruction it optimizes over is
// the A11 implementation verified there, unchanged.
//
// ACCEPTANCE: 0 failed, 0 unresolved mathematical skips. With OTIlib disabled
// only the HyPAD gradient cases may skip; every optimization and
// double-precision inverse-design case must still run.

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/program_b2";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            std::cerr << "usage: b2_verify [--out <dir>]\n";
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
        "program_b2",
        "B2-ABSTRACT-INVERSE: abstract constitutive inverse design. "
        "g_eff^target -> Z_target -> observer blocks -> nearest admissible "
        "constitutive design, over the A11 reconstruction verified at B1C");

    program_b::checks::B2Metrics m;
    program_b::checks::run_b2_checks(report, out_dir, m);

    report.print(std::cout);

    std::string error;
    if (!report.write(out_dir, error)) {
        std::cerr << "failed to write reports: " << error << "\n";
        return 3;
    }

    {
        namespace id = program_b::inverse_design;
        JsonValue meta = JsonValue::object();
        meta.set("milestone", JsonValue("B2"));
        meta.set("scope",
                 JsonValue("Abstract constitutive-design space only. No Maxwell "
                           "solver, no coils or antennas, no metamaterial unit "
                           "cells, no dispersion fitting, no experimental control "
                           "model. Nothing here claims physical realizability of "
                           "any computed design, and nothing generates measurable "
                           "gravitational curvature."));
        meta.set("laboratory_metric", JsonValue("eta = diag(-1,1,1,1)"));
        meta.set("observer", JsonValue("laboratory rest observer V = (1,0,0,0), "
                                      "normalised with g_lab"));
        meta.set("design_family_default",
                 JsonValue(id::to_string(id::DesignFamily::UnrestrictedReciprocal)));
        meta.set("design_family_description",
                 JsonValue(id::describe(id::DesignFamily::UnrestrictedReciprocal)));
        meta.set("parameterization",
                 JsonValue("21 parameters in an orthonormal spatial triad adapted "
                           "to V: 6 symmetric eps, 6 symmetric mu^-1, 9 general "
                           "zeta, with zeta_adjoint = zeta^T structural"));
        meta.set("objective_norm",
                 JsonValue(id::to_string(id::ObjectiveNorm::CoordinateBasis)));
        meta.set("objective_norm_description",
                 JsonValue(id::describe(id::ObjectiveNorm::CoordinateBasis)));
        meta.set("reconstruction",
                 JsonValue("Schuster & Visser (A11), the implementation verified at "
                           "B1C-S456 to 2.2e-16 on random reciprocal media"));
#ifdef WARP_HYPAD_HAVE_OTILIB
        meta.set("oti_backend", JsonValue("enabled"));
        meta.set("oti_directions", JsonValue(static_cast<double>(
                                       id::kGradientBatchSize)));
        meta.set("oti_batches", JsonValue(static_cast<double>(m.oti_batches)));
#else
        meta.set("oti_backend", JsonValue("disabled"));
#endif
        meta.set("build", build_info_json());
        meta.set("timestamp_utc", JsonValue(utc_timestamp()));
        if (!json_write_file(out_dir + "/b2_metadata.json", meta, error)) return 3;
    }

    std::cout << "\nB2 REQUIRED QUANTITIES\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  total cases                        " << report.num_cases() << "\n";
    std::cout << "  passed                             " << report.num_passed() << "\n";
    std::cout << "  failed                             " << report.num_failed() << "\n";
    std::cout << "  skipped                            " << report.num_skipped() << "\n";
    std::cout << "  worst parameter round trip         " << m.worst_parameter_roundtrip << "\n";
    std::cout << "  worst symmetry residual            " << m.worst_symmetry_residual << "\n";
    std::cout << "  worst spatiality residual          " << m.worst_spatiality_residual << "\n";
    std::cout << "  worst objective at exact target    " << m.worst_objective_at_target << "\n";
    std::cout << "  worst gradient vs analytic         " << m.worst_gradient_analytic << "\n";
    std::cout << "  worst gradient vs finite diff      " << m.worst_gradient_fd << "\n";
    if (m.oti_available) {
        std::cout << "  worst batched HyPAD gradient       " << m.worst_gradient_oti << "\n";
        std::cout << "  worst batch independence           " << m.worst_batch_consistency << "\n";
        std::cout << "  HyPAD batches / evaluations        " << m.oti_batches << " / "
                  << m.oti_evaluations << "\n";
    } else {
        std::cout << "  batched HyPAD gradient             skipped (backend disabled)\n";
    }
    std::cout << "  isotropic recovery residual        " << m.isotropic_residual << "\n";
    std::cout << "  anisotropic recovery residual      " << m.anisotropic_residual << "\n";
    std::cout << "  magnetoelectric recovery residual  " << m.magnetoelectric_residual << "\n";
    std::cout << "  metric-induced recovery residual   " << m.metric_induced_residual << "\n";
    std::cout << "  restricted-family residual         " << m.restricted_family_residual
              << "  (must be nonzero)\n";
    std::cout << "  restricted vs least squares        " << m.restricted_vs_least_squares << "\n";
    std::cout << "  infeasible bounded residual        " << m.bounded_residual
              << "  (must be nonzero)\n";
    std::cout << "  worst coordinate-change residual   " << m.worst_coordinate_change << "\n";

    if (!report.all_passed()) {
        std::cout << "\nB2 GATE: FAIL (" << report.num_failed() << " of "
                  << (report.num_cases() - report.num_skipped())
                  << " cases failed)\n";
        return 1;
    }

    std::cout << "\nB2 GATE: PASS (" << report.num_passed() << "/"
              << (report.num_cases() - report.num_skipped()) << ")";
    if (report.num_skipped() > 0)
        std::cout << ", " << report.num_skipped()
                  << " skipped (HyPAD gradient cases only; the backend is "
                     "disabled)";
    std::cout << "\n";
    return 0;
}
