#include <filesystem>
#include <iostream>
#include <string>

#include "temporal_checks.hpp"

// Milestone A2T gate: temporal derivative operators.
//
// Blocks the time-dependent metric and everything downstream of it. The
// rigid-translation identity verified here is the technical basis of the
// formulation choice in docs/time_dependence_decision.md.

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/milestone_02t";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            std::cerr << "usage: milestone2t_verify [--out <dir>]\n";
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
        "milestone_02t",
        "Temporal derivative operators, and the rigid-translation identity "
        "d_t g = -v_s d_x g");

    checks::run_temporal_checks(report, out_dir);
    report.print(std::cout);

    std::string error;
    if (!report.write(out_dir, error)) {
        std::cerr << "failed to write verification reports: " << error << "\n";
        return 3;
    }

    if (!report.all_passed()) {
        std::cout << "\nMILESTONE A2T GATE: FAIL (" << report.num_failed()
                  << " of " << (report.num_cases() - report.num_skipped())
                  << " cases failed)\n";
        return 1;
    }
    std::cout << "\nMILESTONE A2T GATE: PASS (" << report.num_passed() << "/"
              << (report.num_cases() - report.num_skipped()) << ")";
    if (report.num_skipped() > 0)
        std::cout << ", " << report.num_skipped() << " skipped";
    std::cout << "\n";
    return 0;
}
