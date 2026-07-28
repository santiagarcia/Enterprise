#include <filesystem>
#include <iostream>
#include <string>

#include "observer_checks.hpp"

// Milestone A6 gate: observer frames and energy conditions.
//
//   g -> ADM/tetrad -> T_ab -> rho, q_i, S_ij -> NEC/WEC/DEC/SEC
//
// SCOPE, printed on every run and written into every record:
//
//   A6 evaluates observer-dependent energy density and classical pointwise
//   energy conditions for the stress-energy tensor inferred from the selected
//   metric. It does not demonstrate that the source can be generated, that
//   quantum energy inequalities are satisfied, that the geometry is dynamically
//   stable, or that operational FTL travel is possible.
//
// Exit 0 only if every case passes. Metric optimization must not begin
// otherwise: the objective it would minimise is defined here, and an
// unverified energy density would be optimised into a confident wrong answer.

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/milestone_06";
    std::string group = "all";

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else if (a == "--group" && i + 1 < argc) {
            group = argv[++i];
        } else {
            std::cerr << "usage: milestone6_verify [--out <dir>] [--group "
                         "frames|observers|fixtures|nec|alcubierre|integrals|"
                         "all]\n";
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
        "milestone_06",
        "Observer frames and energy conditions: ADM decomposition, orthonormal "
        "tetrads, typed observer families, the frame stress-energy "
        "decomposition, the exact NEC minimiser on the sphere, WEC/DEC/SEC "
        "certification, and the integrated optimization-ready objectives");

    const bool all = (group == "all");
    if (all || group == "frames") checks::run_frame_checks(report, out_dir);
    if (all || group == "observers") checks::run_observer_checks(report, out_dir);
    if (all || group == "fixtures") checks::run_fixture_checks(report, out_dir);
    if (all || group == "nec") checks::run_nec_checks(report, out_dir);
    if (all || group == "alcubierre")
        checks::run_alcubierre_energy_checks(report, out_dir);
    if (all || group == "integrals") checks::run_integral_checks(report, out_dir);

    if (report.num_cases() == 0) {
        std::cerr << "no cases ran: unknown group '" << group << "'\n";
        return 2;
    }

    report.print(std::cout);

    if (all) {
        std::string error;
        if (!report.write(out_dir, error)) {
            std::cerr << "failed to write verification reports: " << error << "\n";
            return 3;
        }
    }

    std::cout
        << "\nSCOPE: A6 evaluates observer-dependent energy density and\n"
           "classical pointwise energy conditions for the stress-energy tensor\n"
           "inferred from the selected metric. It does not demonstrate that the\n"
           "source can be generated, that quantum energy inequalities are\n"
           "satisfied, that the geometry is dynamically stable, or that\n"
           "operational FTL travel is possible.\n";

    if (!report.all_passed()) {
        std::cout << "\nMILESTONE A6 GATE: FAIL (" << report.num_failed()
                  << " of " << report.num_cases() << " cases failed)\n"
                  << "Metric optimization must not begin: the objective it "
                     "would minimise is defined here.\n";
        return 1;
    }

    std::cout << "\nMILESTONE A6 GATE: PASS (" << report.num_passed() << "/"
              << report.num_cases() << ")\n";
    return 0;
}
