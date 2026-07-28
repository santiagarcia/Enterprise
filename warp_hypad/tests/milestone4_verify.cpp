#include <filesystem>
#include <iostream>
#include <string>

#include "curvature_checks.hpp"

// Milestone A4 gate: the covariant curvature pipeline.
//
//   g_{mu nu} -> Gamma^rho_{mu nu} -> R^rho_{sigma mu nu} -> R_{mu nu} -> R
//             -> R_{mu nu} R^{mu nu}, R_{abcd} R^{abcd}
//
// Exit 0 only if every case passes. Milestone A5 (Einstein tensor and
// stress-energy) must not begin otherwise: an Einstein tensor built on an
// unverified Ricci tensor would put a sign error into every energy density the
// project ever reports, and no later check would be able to attribute it.

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/milestone_04";
    std::string group = "all";

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else if (a == "--group" && i + 1 < argc) {
            group = argv[++i];
        } else {
            std::cerr << "usage: milestone4_verify [--out <dir>] "
                         "[--group connection|analytic|temporal|alcubierre|"
                         "convergence|all]\n";
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
        "milestone_04",
        "Covariant curvature: Christoffel symbols, Riemann and Ricci tensors, "
        "curvature invariants, and the temporal-derivative policy they require");

    const bool all = (group == "all");
    if (all || group == "connection")
        checks::run_connection_checks(report, out_dir);
    if (all || group == "analytic")
        checks::run_analytic_curvature_checks(report, out_dir);
    if (all || group == "temporal")
        checks::run_curvature_temporal_checks(report, out_dir);
    if (all || group == "alcubierre")
        checks::run_alcubierre_curvature_checks(report, out_dir);
    if (all || group == "convergence")
        checks::run_curvature_convergence_checks(report, out_dir);

    if (report.num_cases() == 0) {
        std::cerr << "no cases ran: unknown group '" << group << "'\n";
        return 2;
    }

    report.print(std::cout);

    // Only the full run writes the record. A per-group invocation is a debugging
    // convenience, and letting it overwrite the reports would leave a directory
    // claiming a complete milestone on the evidence of one group.
    if (all) {
        std::string error;
        if (!report.write(out_dir, error)) {
            std::cerr << "failed to write verification reports: " << error << "\n";
            return 3;
        }
    }

    if (!report.all_passed()) {
        std::cout << "\nMILESTONE A4 GATE: FAIL (" << report.num_failed() << " of "
                  << report.num_cases() << " cases failed)\n"
                  << "Milestone A5 must not begin: an Einstein tensor built on an "
                     "unverified Ricci tensor cannot be debugged later.\n";
        return 1;
    }

    std::cout << "\nMILESTONE A4 GATE: PASS (" << report.num_passed() << "/"
              << report.num_cases() << ")\n";
    return 0;
}
