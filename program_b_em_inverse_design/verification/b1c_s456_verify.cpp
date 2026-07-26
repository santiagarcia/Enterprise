#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>

#include "b1c_checks.hpp"

// Gate B1C-S456.
//
//   step 4  H^{ab} = Z^{abcd}F_{cd} by two independent routes
//   step 5  zeta_adjoint, reciprocity, and the handedness benchmark
//   step 6  reconstruction of Z from the observer blocks via (A11), with
//           observer-change and coordinate-change consistency
//
// Depends on B1C-S23, which must already pass 42/42: this gate takes the Hodge
// duals, the (A7) extraction and the pinned zeta normalisation as given.
//
// ACCEPTANCE: zero failures and zero skips. A skip here would mean an
// unresolved normalisation or sign, which is exactly the state this gate exists
// to end.

int main(int argc, char** argv) {
    std::string out_dir = "results/verification/program_b1c_s456";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            std::cerr << "usage: b1c_s456_verify [--out <dir>]\n";
            return 2;
        }
    }

    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if (ec) {
        std::cerr << "cannot create " << out_dir << ": " << ec.message() << "\n";
        return 3;
    }

    warp_hypad::VerificationReport report(
        "program_b1c_s456",
        "B1C-S456: full electromagnetic field reconstruction H^{ab} = "
        "Z^{abcd}F_{cd}, the magnetoelectric remainder (zeta_adjoint, "
        "reciprocity, handedness), and reconstruction of Z from the observer "
        "blocks via Schuster & Visser (A11)");

    program_b::checks::S456Metrics m;
    program_b::checks::run_s456_checks(report, m);

    report.print(std::cout);

    std::string error;
    if (!report.write(out_dir, error)) {
        std::cerr << "failed to write reports: " << error << "\n";
        return 3;
    }

    std::cout << "\nB1C-S456 REQUIRED QUANTITIES\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  total cases                        " << report.num_cases() << "\n";
    std::cout << "  passed                             " << report.num_passed() << "\n";
    std::cout << "  failed                             " << report.num_failed() << "\n";
    std::cout << "  skipped                            " << report.num_skipped() << "\n";
    std::cout << "  worst direct vs block H^{ab}       " << m.worst_direct_vs_block_H << "\n";
    std::cout << "  worst F reconstruction             " << m.worst_F_reconstruction << "\n";
    std::cout << "  B17 zeta error                     " << m.b17_zeta_error << "\n";
    std::cout << "  B17 zeta_adjoint error             " << m.b17_zeta_adjoint_error << "\n";
    std::cout << "  worst A11 Z round trip             " << m.worst_A11_roundtrip << "\n";
    std::cout << "  worst observer-change error        " << m.worst_observer_change << "\n";
    std::cout << "  worst coordinate-change error      " << m.worst_coordinate_change << "\n";
    std::cout << "  handedness                         " << m.handedness_result << "\n";

    if (!report.all_passed()) {
        std::cout << "\nB1C-S456 GATE: FAIL (" << report.num_failed() << " of "
                  << (report.num_cases() - report.num_skipped())
                  << " cases failed)\n";
        return 1;
    }
    if (report.num_skipped() > 0) {
        std::cout << "\nB1C-S456 GATE: FAIL (" << report.num_skipped()
                  << " skipped; this gate does not permit unresolved "
                     "normalisation or sign cases)\n";
        return 1;
    }

    std::cout << "\nB1C-S456 GATE: PASS (" << report.num_passed() << "/"
              << report.num_cases() << "), 0 failures, 0 skips\n";
    return 0;
}
