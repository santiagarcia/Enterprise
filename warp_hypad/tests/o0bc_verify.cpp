#include <filesystem>
#include <iostream>
#include <string>

#include "o0b_checks.hpp"

// Gate O0B-C-CORRECTED-DESIGN-SEMANTICS.
//
// A bounded corrective milestone over the verified O0B design layer. Each case
// exists because an O0B claim was wrong, weaker than stated, or measured
// something other than what its name said.
//
// SCOPE: O0B-C corrects derivative, optimizer-coordinate and
// numerical-integration semantics in the verified ADM design layer. It does not
// evaluate energy, run an optimizer, demonstrate a causal shortcut or establish
// physical realizability.

int main(int argc, char** argv) {
    using namespace warp_hypad;
    std::string out_dir = "results/verification/o0bc_corrections";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) out_dir = argv[++i];
        else { std::cerr << "usage: o0bc_verify [--out <dir>]\n"; return 2; }
    }
    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if (ec) { std::cerr << "cannot create " << out_dir << "\n"; return 3; }

    VerificationReport report(
        "o0bc_corrections",
        "O0B-C corrected design semantics: exact axis jets with no production "
        "finite difference, separated raw / constrained / whitened design "
        "gradients, directly measured constraint closure, a floating-point-safe "
        "exponential parameterisation, corrected quadrature exactness "
        "semantics, and a recorded optimizer-facing API surface");

    checks::run_o0bc_checks(report, out_dir);
    report.print(std::cout);

    std::string error;
    if (!report.write(out_dir, error)) {
        std::cerr << "failed to write reports: " << error << "\n";
        return 3;
    }

    std::cout << "\nSCOPE: O0B-C corrects derivative, optimizer-coordinate and\n"
                 "numerical-integration semantics in the verified ADM design\n"
                 "layer. It does not evaluate energy, run an optimizer,\n"
                 "demonstrate a causal shortcut or establish physical\n"
                 "realizability.\n";

    if (!report.all_passed()) {
        std::cout << "\nO0B-C GATE: FAIL (" << report.num_failed() << " of "
                  << report.num_cases() << " cases failed)\n";
        return 1;
    }
    std::cout << "\nO0B-C GATE: PASS (" << report.num_passed() << "/"
              << report.num_cases() << ", " << report.num_skipped()
              << " skipped)\n";
    return 0;
}
