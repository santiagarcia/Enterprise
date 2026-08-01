#include <filesystem>
#include <iostream>
#include <string>

#include "o0d_checks.hpp"

// Milestone O0D gate: certify the resolution of the F1W objective.
//
// O0C left exactly one flag false. Its convergence study certified the BASELINE
// objective on the elementwise route -- the F1 family on F1's lattice -- while
// the objective O1 would optimize is the F1W objective on F1W's own lattice.
// O0D runs that study, certifies the gradient's integration resolution, and
// exercises the candidate-validation machinery on deterministic nearby points.
//
// SCOPE, printed on every run and written into every record:
//
//   O0D certifies the numerical resolution of the F1W negative-energy
//   objective and of its whitened gradient on the F1W lattice. It does not
//   optimize, demonstrate a causal shortcut, establish physical realizability
//   or claim any reduction in energy.
//
// Exit 0 only if every case passes. No optimizer runs anywhere in this
// milestone; O0D decides whether O1 MAY run.

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/o0d_f1w_certification";
    std::string group = "all";

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) out_dir = argv[++i];
        else if (a == "--group" && i + 1 < argc) group = argv[++i];
        else {
            std::cerr << "usage: o0d_verify [--out <dir>] "
                         "[--group convergence|gradient|probes|all]\n";
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
        "o0d_f1w_certification",
        "F1W objective-resolution certification: the complete order, "
        "subdivision, domain-extent and independent-reference convergence "
        "study on the F1W family's own lattice at y = 0, certification that "
        "the gradient the optimizer would use resolves the same objective the "
        "validation evaluator resolves, and a deterministic exercise of the "
        "candidate-validation machinery on smooth directions and on every "
        "adversary O0C identified");

    const bool all = (group == "all");
    if (all || group == "convergence")
        checks::run_o0d_convergence_checks(report, out_dir);
    if (all || group == "gradient")
        checks::run_o0d_gradient_checks(report, out_dir);
    if (all || group == "probes")
        checks::run_o0d_probe_checks(report, out_dir);
    if (all) checks::write_o0d_readiness(report, out_dir);

    if (report.num_cases() == 0) {
        std::cerr << "no cases ran: unknown group '" << group << "'\n";
        return 2;
    }

    report.print(std::cout);

    std::string err;
    if (!report.write(out_dir, err)) {
        std::cerr << "cannot write reports: " << err << "\n";
        return 4;
    }
    if (all) checks::write_o0d_manifest(report, out_dir);

    std::cout << "\nSCOPE: O0D certifies the numerical resolution of the F1W\n"
                 "negative-energy objective and of its whitened gradient on\n"
                 "the F1W lattice. It does not optimize, demonstrate a causal\n"
                 "shortcut, establish physical realizability or claim any\n"
                 "reduction in energy.\n";

    if (all) {
        const checks::O0dState& st = checks::o0d_state();
        std::cout << "\nREADINESS: resolution_certified = "
                  << (st.resolution_certified ? "true" : "false")
                  << ", ready_for_optimization = "
                  << (st.ready_for_optimization ? "true" : "false")
                  << ", o1_may_begin = "
                  << (st.ready_for_optimization ? "true" : "false") << "\n"
                  << "  " << st.blocking_reason << "\n";
    }

    const bool ok = report.all_passed();
    std::cout << "\nMILESTONE O0D GATE: " << (ok ? "PASS (" : "FAIL (")
              << report.num_passed() << "/" << report.num_cases();
    if (!ok) std::cout << ", " << report.num_failed() << " failed";
    std::cout << ", " << report.num_skipped() << " skipped)\n";
    return ok ? 0 : 1;
}
