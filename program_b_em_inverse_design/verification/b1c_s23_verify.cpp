#include <filesystem>
#include <iostream>
#include <string>

#include "b1c_checks.hpp"

// Gate B1C-S23: Hodge duals (step 1) + Bel extraction (step 2) + spatial
// projector and block spatiality (step 3).
//
// An intermediate gate. Full F->H reconstruction, the (B17) magnetoelectric
// benchmark and the (A11) Z reconstruction must not begin until this passes.

int main(int argc, char** argv) {
    std::string out_dir = "results/verification/program_b1c_s23";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            std::cerr << "usage: b1c_s23_verify [--out <dir>]\n";
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
        "program_b1c_s23",
        "B1C-S23: four-dimensional Hodge duals, Bel extraction from Schuster & "
        "Visser (A7a)-(A7d), spatial projector and block spatiality");

    program_b::checks::run_dual_checks(report);
    program_b::checks::run_s23_checks(report);
    program_b::checks::run_s5_checks(report);

    report.print(std::cout);

    std::string error;
    if (!report.write(out_dir, error)) {
        std::cerr << "failed to write reports: " << error << "\n";
        return 3;
    }

    if (!report.all_passed()) {
        std::cout << "\nB1C-S23 GATE: FAIL (" << report.num_failed() << " of "
                  << (report.num_cases() - report.num_skipped())
                  << " cases failed)\n";
        return 1;
    }

    std::cout << "\nB1C-S23 GATE: PASS (" << report.num_passed() << "/"
              << (report.num_cases() - report.num_skipped()) << ")";
    if (report.num_skipped() > 0) {
        std::cout << ", " << report.num_skipped()
                  << " recorded as unresolved (zeta normalisation, step 5)";
    }
    std::cout << "\n";
    return 0;
}
