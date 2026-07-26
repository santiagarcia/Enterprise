#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>

#include "b4_checks.hpp"
#include "program_b/inverse_design/maxwell.hpp"
#include "warp_hypad/core/metadata.hpp"

// Gate B4-MAXWELL-PROPAGATION.
//
// SCOPE. A passing result demonstrates that the B3 constitutive model reproduces
// selected electromagnetic propagation observables of the target effective
// metric. It still does not establish a fabricated metamaterial, unit cell,
// antenna, coil, or experimental device.

int main(int argc, char** argv) {
    using namespace warp_hypad;
    namespace id = program_b::inverse_design;

    std::string out_dir = "results/verification/program_b4";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            std::cerr << "usage: b4_verify [--out <dir>]\n";
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
        "program_b4",
        "B4-MAXWELL-PROPAGATION: does the B3 constitutive model produce the "
        "electromagnetic propagation the target effective metric predicts? "
        "Plane-wave operator, Fresnel polynomial, layered propagation and energy "
        "balance");

    program_b::checks::B4Metrics m;
    program_b::checks::run_b4_checks(report, out_dir, m);

    report.print(std::cout);

    std::string error;
    if (!report.write(out_dir, error)) {
        std::cerr << "failed to write reports: " << error << "\n";
        return 3;
    }

    {
        const id::MaxwellConventions mx;
        JsonValue meta = JsonValue::object();
        meta.set("milestone", JsonValue("B4"));
        meta.set("scope",
                 JsonValue("A passing result demonstrates that the B3 "
                           "constitutive model reproduces selected "
                           "electromagnetic propagation observables of the target "
                           "effective metric. It still does not establish a "
                           "fabricated metamaterial, unit cell, antenna, coil, or "
                           "experimental device."));
        meta.set("fourier_kernel", JsonValue(mx.fourier_kernel));
        meta.set("faraday", JsonValue(mx.faraday));
        meta.set("ampere", JsonValue(mx.ampere));
        meta.set("gauss_magnetic", JsonValue(mx.gauss_magnetic));
        meta.set("gauss_electric", JsonValue(mx.gauss_electric));
        meta.set("plane_wave", JsonValue(mx.plane_wave));
        meta.set("poynting", JsonValue(mx.poynting));
        meta.set("sign_provenance", JsonValue(mx.note));
        meta.set("solver",
                 JsonValue("1D frequency-domain layered medium, scattering "
                           "matrices with Redheffer star product; only decaying "
                           "exponentials appear, so thick evanescent stacks do "
                           "not overflow"));
        meta.set("fresnel",
                 JsonValue("adj(W) = lambda k (x) k with W^{bd} = Z^{abcd}k_a k_c; "
                           "lambda is the quartic Fresnel scalar, and det W "
                           "vanishes identically"));
        meta.set("passivity_scope",
                 JsonValue("Causality is global and structural. Passivity is "
                           "verified over the tested band only, now narrowed by "
                           "the adaptive scan. Global passivity is NOT claimed."));
        meta.set("build", build_info_json());
        meta.set("timestamp_utc", JsonValue(utc_timestamp()));
        if (!json_write_file(out_dir + "/b4_metadata.json", meta, error)) return 3;
    }

    std::cout << "\nB4 REQUIRED QUANTITIES\n";
    std::cout << "  Fourier convention                 " << m.fourier_convention << "\n";
    std::cout << "  Maxwell curl signs                 " << m.maxwell_curl_signs << "\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  total / passed / failed / skipped  " << report.num_cases() << " / "
              << report.num_passed() << " / " << report.num_failed() << " / "
              << report.num_skipped() << "\n";
    std::cout << "  worst plane-wave residual          " << m.worst_plane_wave_residual << "\n";
    std::cout << "  worst transversality residual      " << m.worst_transversality_residual << "\n";
    std::cout << "  Fresnel factorization error        " << m.fresnel_factorization_error << "\n";
    std::cout << "  generic birefringence split        " << m.generic_birefringence_split
              << "   (must be nonzero)\n";
    std::cout << "  group-velocity error               " << m.group_velocity_error << "\n";
    std::cout << "  coarse min passivity margin        " << m.coarse_min_passivity << "\n";
    std::cout << "  ADAPTIVE min passivity margin      " << m.adaptive_min_passivity
              << "   (" << m.adaptive_samples << " samples)\n";
    std::cout << "  adaptive improvement over coarse   " << m.adaptive_improvement << "\n";
    std::cout << "  worst energy-balance residual      " << m.worst_energy_balance << "\n";
    std::cout << "  R / T / A (lossy slab)             " << m.reflection << " / "
              << m.transmission << " / " << m.absorption << "\n";
    std::cout << "  phase mismatch, exact vs zeta=0    " << m.phase_mismatch << "\n";
    std::cout << "  group-delay mismatch               " << m.group_delay_mismatch << "\n";
    std::cout << "  exact vs restricted |dT|           " << m.exact_vs_restricted_error << "\n";
    if (m.oti_available)
        std::cout << "  HyPAD gradient error               " << m.hypad_gradient_error << "\n";
    else
        std::cout << "  HyPAD gradient                     skipped (backend disabled)\n";
    std::cout << "  branch / degeneracy report         " << m.branch_warnings << "\n";

    std::cout << "\nSCOPE: a passing result shows the B3 constitutive model "
                 "reproduces selected\n"
                 "  propagation observables of the target effective metric. It is "
                 "NOT a fabricated\n"
                 "  metamaterial, unit cell, antenna, coil, or experimental "
                 "device.\n"
                 "  Causality is global; passivity is verified over the tested "
                 "band only.\n";

    if (!report.all_passed()) {
        std::cout << "\nB4 GATE: FAIL (" << report.num_failed() << " of "
                  << (report.num_cases() - report.num_skipped())
                  << " cases failed)\n";
        return 1;
    }

    std::cout << "\nB4 GATE: PASS (" << report.num_passed() << "/"
              << (report.num_cases() - report.num_skipped()) << ")";
    if (report.num_skipped() > 0)
        std::cout << ", " << report.num_skipped()
                  << " skipped (HyPAD sensitivity cases only)";
    std::cout << "\n";
    return 0;
}
