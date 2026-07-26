#include <filesystem>
#include <iostream>
#include <string>

#include "derivative_checks.hpp"

// Milestone A2 gate: numerical derivative operators.
// Exit 0 only if every stencil test passes AND the observed convergence order
// is demonstrated, interior and boundary separately.

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/milestone_02";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            std::cerr << "usage: milestone2_verify [--out <dir>]\n";
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
        "milestone_02",
        "Finite-difference operators: first, second and mixed derivatives at "
        "second order, interior and boundary");

    checks::run_derivative_checks(report, out_dir);
    report.print(std::cout);

    std::string error;
    if (!report.write(out_dir, error)) {
        std::cerr << "failed to write verification reports: " << error << "\n";
        return 3;
    }

    if (!report.all_passed()) {
        std::cout << "\nMILESTONE A2 GATE: FAIL (" << report.num_failed() << " of "
                  << report.num_cases() << " cases failed)\n"
                  << "Milestone A3 must not begin until the expected order is "
                     "demonstrated.\n";
        return 1;
    }

    std::cout << "\nMILESTONE A2 GATE: PASS (" << report.num_passed() << "/"
              << report.num_cases() << ")\n";
    return 0;
}
