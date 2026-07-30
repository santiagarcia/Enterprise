#include <filesystem>
#include <iostream>
#include <string>

#include "o0a_checks.hpp"

// Milestone O0A gate: the constrained quintic spline basis.
//
//   knots -> univariate quintic kernel -> axisymmetric tensor product
//         -> structurally regular field families -> exact constraints
//         -> Sobolev whitening -> sparse stencils -> coefficient derivatives
//
// SCOPE, printed on every run and written into every record:
//
//   O0A constructs and verifies a constrained spline coordinate system for
//   future metric design. It does not yet optimize a spacetime metric, reduce
//   energy requirements, demonstrate a causal shortcut or establish physical
//   realizability.
//
// Exit 0 only if every case passes. A family may be verified here and still be
// classified as not ready for optimization: the readiness report records that
// separately, and an honestly measured resolution limitation is a RESULT of this
// milestone rather than a failure of it.

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/o0a_basis";
    std::string group = "all";

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else if (a == "--group" && i + 1 < argc) {
            group = argv[++i];
        } else {
            std::cerr << "usage: o0a_verify [--out <dir>] [--group "
                         "kernel|tensor|constraints|projection|derivatives|all]\n";
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
        "o0a_basis",
        "Constrained quintic spline basis: knot vectors and C4 continuity, the "
        "axisymmetric tensor product and its Cartesian derivatives, "
        "structurally regular field families, exact constraint elimination, "
        "Sobolev whitening, sparse stencils, and coefficient derivatives "
        "against the OTI backend");

    const bool all = (group == "all");
    if (all || group == "kernel") checks::run_o0a_kernel_checks(report, out_dir);
    if (all || group == "tensor") checks::run_o0a_tensor_checks(report, out_dir);
    if (all || group == "constraints")
        checks::run_o0a_constraint_checks(report, out_dir);
    if (all || group == "projection")
        checks::run_o0a_projection_checks(report, out_dir);
    if (all || group == "derivatives")
        checks::run_o0a_derivative_checks(report, out_dir);

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
        // Readiness first, then the manifest: the manifest digests every record
        // in the directory except itself, so it must be written last.
        checks::write_o0a_frame_records(out_dir);
        checks::write_o0a_readiness(report, out_dir);
        checks::write_o0a_manifest(report, out_dir);
    }

    std::cout
        << "\nSCOPE: O0A constructs and verifies a constrained spline\n"
           "coordinate system for future metric design. It does not yet\n"
           "optimize a spacetime metric, reduce energy requirements,\n"
           "demonstrate a causal shortcut or establish physical realizability.\n";

    if (!report.all_passed()) {
        std::cout << "\nMILESTONE O0A GATE: FAIL (" << report.num_failed()
                  << " of " << report.num_cases() << " cases failed)\n";
        return 1;
    }

    std::cout << "\nMILESTONE O0A GATE: PASS (" << report.num_passed() << "/"
              << report.num_cases() << ", " << report.num_skipped()
              << " skipped)\n";
    return 0;
}
