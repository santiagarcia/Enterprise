#include <filesystem>
#include <iostream>
#include <string>

#include "einstein_checks.hpp"

// Milestone A5 gate: the Einstein tensor and the inferred stress-energy.
//
//   R_{mu nu}, R, g_{mu nu} -> G_{mu nu} -> T_{mu nu}^inferred
//
// SCOPE, printed on every run and written into every record:
//
//   A5 computes the geometric Einstein tensor and the stress-energy tensor
//   inferred through Einstein's equation for the selected metric. It does not
//   establish that this source can be generated physically, that energy
//   conditions are satisfied, or that FTL travel is realizable.
//
// Exit 0 only if every case passes. Milestone A6 must not begin otherwise: it
// contracts this tensor against observer four-velocities, and an error here
// would appear there as a wrong energy density with no way to attribute it.

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/milestone_05";
    std::string group = "all";

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else if (a == "--group" && i + 1 < argc) {
            group = argv[++i];
        } else {
            std::cerr << "usage: milestone5_verify [--out <dir>] [--group "
                         "reference|desitter|alcubierre|identities|bianchi|"
                         "convergence|units|policy|all]\n";
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
        "milestone_05",
        "Einstein tensor and inferred stress-energy: G = R_{mu nu} - g R / 2, "
        "raw and symmetric, its trace and square, the contracted Bianchi "
        "residual, cosmological-constant treatment and the units layer");

    const bool all = (group == "all");
    if (all || group == "reference")
        checks::run_einstein_reference_checks(report, out_dir);
    if (all || group == "desitter")
        checks::run_einstein_desitter_checks(report, out_dir);
    if (all || group == "alcubierre")
        checks::run_einstein_alcubierre_checks(report, out_dir);
    if (all || group == "identities")
        checks::run_einstein_identity_checks(report, out_dir);
    if (all || group == "bianchi")
        checks::run_bianchi_checks(report, out_dir);
    if (all || group == "convergence")
        checks::run_einstein_convergence_checks(report, out_dir);
    if (all || group == "units")
        checks::run_units_checks(report, out_dir);
    if (all || group == "policy")
        checks::run_einstein_policy_checks(report, out_dir);
    if (all || group == "lambda")
        checks::run_lambda_source_checks(report, out_dir);

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

    std::cout << "\nSCOPE: A5 computes the geometric Einstein tensor and the\n"
                 "stress-energy tensor inferred through Einstein's equation for\n"
                 "the selected metric. It does not establish that this source\n"
                 "can be generated physically, that energy conditions are\n"
                 "satisfied, or that FTL travel is realizable.\n";

    if (!report.all_passed()) {
        std::cout << "\nMILESTONE A5 GATE: FAIL (" << report.num_failed() << " of "
                  << report.num_cases() << " cases failed)\n"
                  << "Milestone A6 must not begin: it contracts this tensor "
                     "against observer four-velocities.\n";
        return 1;
    }

    std::cout << "\nMILESTONE A5 GATE: PASS (" << report.num_passed() << "/"
              << report.num_cases() << ")\n";
    return 0;
}
