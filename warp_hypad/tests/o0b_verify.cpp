#include <filesystem>
#include <iostream>
#include <string>

#include "o0b_checks.hpp"

// Milestone O0B gate: constrained ADM metric and analytic metric jets.
//
//   verified spline coordinates -> ADM design fields -> g_mu_nu, g^mu_nu
//                               -> analytic derivative jet through third order
//
// SCOPE, printed on every run and written into every record:
//
//   O0B maps verified spline coordinates into structurally valid ADM metrics
//   and analytic metric derivatives. It does not yet minimize energy, certify
//   the spline family for energy evaluation, demonstrate a causal shortcut or
//   establish physical realizability.
//
// Exit 0 only if every case passes. No energy is evaluated and no optimizer
// runs anywhere in this milestone.

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/o0b_metric";
    std::string group = "all";

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) out_dir = argv[++i];
        else if (a == "--group" && i + 1 < argc) group = argv[++i];
        else {
            std::cerr << "usage: o0b_verify [--out <dir>] "
                         "[--group metric|routes|all]\n";
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
        "o0b_metric",
        "Constrained ADM metric and analytic metric jets: typed coefficient "
        "spaces, structurally valid ADM assembly, the analytic inverse, "
        "independent references against Minkowski, a constant shift and the "
        "production Alcubierre metric, the third-order derivative jet with "
        "genuine temporal channels, coefficient gradients through the OTI "
        "backend, and a costed comparison of two future evaluation routes");

    const bool all = (group == "all");
    if (all || group == "metric") checks::run_o0b_metric_checks(report, out_dir);
    if (all || group == "routes") checks::run_o0b_route_checks(report, out_dir);

    if (report.num_cases() == 0) {
        std::cerr << "no cases ran: unknown group '" << group << "'\n";
        return 2;
    }

    report.print(std::cout);

    if (all) {
        std::string error;
        if (!report.write(out_dir, error)) {
            std::cerr << "failed to write verification reports: " << error
                      << "\n";
            return 3;
        }
        checks::write_o0b_frame_records(out_dir);
        checks::write_o0a_consistency_record(out_dir);
        checks::write_o0b_readiness(report, out_dir);
        // Last: the manifest digests every record except itself.
        checks::write_o0b_manifest(report, out_dir);
    }

    std::cout << "\nSCOPE: O0B maps verified spline coordinates into\n"
                 "structurally valid ADM metrics and analytic metric\n"
                 "derivatives. It does not yet minimize energy, certify the\n"
                 "spline family for energy evaluation, demonstrate a causal\n"
                 "shortcut or establish physical realizability.\n";

    if (!report.all_passed()) {
        std::cout << "\nMILESTONE O0B GATE: FAIL (" << report.num_failed()
                  << " of " << report.num_cases() << " cases failed)\n";
        return 1;
    }
    std::cout << "\nMILESTONE O0B GATE: PASS (" << report.num_passed() << "/"
              << report.num_cases() << ", " << report.num_skipped()
              << " skipped)\n";
    return 0;
}
