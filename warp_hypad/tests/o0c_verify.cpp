#include <filesystem>
#include <iostream>
#include <string>

#include "o0c_checks.hpp"

// Milestone O0C gate: converged stress-energy and energy-condition objectives.
//
//   y -> c = c0 + PW y -> g, dg, d2g, d3g -> Gamma, dGamma -> R^rho_{sigma mu nu}
//     -> R_{mu nu}, R -> G_{mu nu} -> T_{a b} -> rho, N_min
//     -> E_net, E_plus, E_neg, E_abs, J_NEC and the smooth surrogate E_neg,tau
//
// Every step is an existing A4/A5/A6/O0A/O0B kernel except the analytic
// connection derivative, which A4 has no equivalent of because its grid route
// differences Gamma. No production finite difference appears anywhere in the
// geometry route.
//
// SCOPE, printed on every run and written into every record:
//
//   O0C evaluates and verifies curvature, inferred stress-energy and converged
//   energy-condition objectives for behavior-preserving metric families. It
//   does not yet optimize those families, demonstrate a causal shortcut or
//   establish physical realizability.
//
// Exit 0 only if every case passes. No optimizer runs anywhere in this
// milestone.

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/o0c_energy";
    std::string group = "all";

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) out_dir = argv[++i];
        else if (a == "--group" && i + 1 < argc) group = argv[++i];
        else {
            std::cerr << "usage: o0c_verify [--out <dir>] "
                         "[--group geometry|energy|convergence|gradients|all]\n";
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
        "o0c_energy",
        "Converged stress-energy and energy-condition objectives: an analytic "
        "curvature bridge over the O0B metric jet with no production finite "
        "difference, the F1W wall-only family whose protected regions are "
        "enforced by support exclusion rather than by penalty, observer-frame "
        "energy against the analytic Alcubierre reference, the exact global "
        "null minimum, integrated objectives converged in quadrature order, "
        "element subdivision and domain extent against two independent "
        "references, and whitened objective gradients through deterministic "
        "OTI batches");

    const bool all = (group == "all");
    if (all || group == "geometry")
        checks::run_o0c_geometry_checks(report, out_dir);
    if (all || group == "energy")
        checks::run_o0c_energy_checks(report, out_dir);
    if (all || group == "convergence")
        checks::run_o0c_convergence_checks(report, out_dir);
    if (all || group == "gradients")
        checks::run_o0c_gradient_checks(report, out_dir);
    if (all) {
        checks::write_o0c_field_records(out_dir);
        checks::run_o0c_product_checks(report, out_dir);
    }

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
    if (all) checks::write_o0c_manifest(report, out_dir);

    std::cout << "\nSCOPE: O0C evaluates and verifies curvature, inferred\n"
                 "stress-energy and converged energy-condition objectives for\n"
                 "behavior-preserving metric families. It does not yet\n"
                 "optimize those families, demonstrate a causal shortcut or\n"
                 "establish physical realizability.\n";

    const bool ok = report.all_passed();
    std::cout << "\nMILESTONE O0C GATE: " << (ok ? "PASS (" : "FAIL (")
              << report.num_passed() << "/" << report.num_cases();
    if (!ok) std::cout << ", " << report.num_failed() << " failed";
    std::cout << ", " << report.num_skipped() << " skipped)\n";
    return ok ? 0 : 1;
}
