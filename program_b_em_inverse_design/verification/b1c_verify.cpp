#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "b1c_checks.hpp"
#include "warp_hypad/core/metadata.hpp"

// Program B Milestone B1C: correction gate.
//
// B1 concluded "no equivalent medium exists" where g_00 >= 0. That was an
// observer-dependent result stated as an invariant one. This gate establishes
// the corrected picture and must pass before B2 begins.

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/program_b1c";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) out_dir = argv[++i];
        else { std::cerr << "usage: b1c_verify [--out <dir>]\n"; return 2; }
    }

    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if (ec) { std::cerr << "cannot create " << out_dir << "\n"; return 3; }

    VerificationReport report(
        "program_b1c",
        "Correction gate: covariant Z^abcd versus observer-dependent 3+1 "
        "decomposition; the v_s f = 1 surface is an ergosurface, not an "
        "invariant realizability wall");

    program_b::checks::VelocityMetrics velocity;
    program_b::checks::run_dual_checks(report);
    program_b::checks::run_b1c_checks(report, out_dir, velocity);
    report.print(std::cout);

    std::string error;
    if (!report.write(out_dir, error)) {
        std::cerr << "failed to write reports: " << error << "\n"; return 3;
    }

    {
        JsonValue m = JsonValue::object();
        m.set("milestone", JsonValue("B1C"));
        m.set("corrects", JsonValue("B1's claim that no equivalent medium exists "
                                    "where g_00 >= 0"));
        m.set("corrected_statement",
              JsonValue("The Alcubierre metric is Lorentzian with det g = -1 "
                        "everywhere for every v_s, so Z^abcd exists everywhere. "
                        "The v_s f = 1 surface is an ERGOSURFACE: the "
                        "coordinate-time observer of g_eff ceases to be timelike "
                        "there. For the laboratory observer the PERMEABILITY "
                        "diverges -- the transverse eigenvalue of mu^-1 equals "
                        "-g_00 and reaches zero -- while eps = h and the "
                        "covariant tensor stay bounded. cond(eps) does NOT "
                        "diverge and the assertion that it must has been "
                        "withdrawn; see docs/b1c_aggregate_gate_audit.md."));
        m.set("sub_gates",
              JsonValue("B1C-S23 (Hodge duals, A7 extraction, projector, zeta "
                        "normalisation) and B1C-S456 (F->H reconstruction, "
                        "zeta_adjoint, handedness, A11) are separate binaries "
                        "and separate CTest tests. This gate depends on them "
                        "through CTest rather than repeating their cases."));
        m.set("primary_source",
              JsonValue("S. Schuster and M. Visser, Phys. Rev. D 96, 124019 "
                        "(2017), DOI 10.1103/PhysRevD.96.124019 (verified)"));
        m.set("build", build_info_json());
        m.set("timestamp_utc", JsonValue(utc_timestamp()));
        if (!json_write_file(out_dir + "/b1c_metadata.json", m, error)) return 3;
    }

    // Consolidated status report. The aggregate gate does not repeat the
    // sub-gates, so its record has to point at theirs.
    {
        std::ofstream cs(out_dir + "/consolidated_status.md", std::ios::binary);
        if (cs) {
            cs << "# B1C consolidated status\n\n"
               << "B1C is covered by three gates. This one holds only the "
                  "cross-cutting cases; the other two are separate binaries and "
                  "separate CTest tests, and this gate DEPENDS on them rather "
                  "than repeating their assertions.\n\n"
               << "| gate | binary | CTest test | record |\n"
               << "|---|---|---|---|\n"
               << "| B1C-S23 | `b1c_s23_verify` | `programb1c_s23.gate` | "
                  "`results/verification/program_b1c_s23/` |\n"
               << "| B1C-S456 | `b1c_s456_verify` | `programb1c_s456.gate` | "
                  "`results/verification/program_b1c_s456/` |\n"
               << "| B1C aggregate | `b1c_verify` | `programb1c.gate` | "
                  "`results/verification/program_b1c/` |\n\n"
               << "## This gate\n\n"
               << "- cases: " << report.num_cases() << "\n"
               << "- passed: " << report.num_passed() << "\n"
               << "- failed: " << report.num_failed() << "\n"
               << "- skipped: " << report.num_skipped() << "\n\n"
               << "## Velocity sensitivity\n\n";
            cs.setf(std::ios::scientific);
            cs.precision(6);
            cs << "- worst OTI versus analytic: " << velocity.worst_analytic
               << "\n- worst analytic versus centred finite differences: "
               << velocity.worst_fd << "\n- worst chain-rule residual: "
               << velocity.worst_chain_rule << "\n\n"
               << "## Withdrawn assertions\n\n"
               << "See `docs/b1c_aggregate_gate_audit.md` for the per-case record "
                  "of what was withdrawn, why, and what supersedes it. No active "
                  "gate includes `legacy_unverified.hpp`.\n";
        }
    }

    std::cout << "\nCORRECTED CLAIM: the v_s f = 1 surface is an ergosurface.\n"
              << "  The covariant constitutive tensor remains valid across it.\n"
              << "  eps = h and max|Z| stay bounded; the PERMEABILITY diverges.\n"
              << "  cond(eps) does NOT diverge -- that assertion is withdrawn.\n";

    std::cout << "\nVELOCITY SENSITIVITY\n";
    std::cout << std::scientific << std::setprecision(6)
              << "  worst OTI vs analytic          " << velocity.worst_analytic
              << "\n  worst analytic vs centred FD   " << velocity.worst_fd
              << "\n  worst chain-rule residual      " << velocity.worst_chain_rule
              << "\n";

    if (!report.all_passed()) {
        std::cout << "\nPROGRAM B1C GATE: FAIL (" << report.num_failed() << " of "
                  << (report.num_cases() - report.num_skipped()) << ")\n"
                  << "B2 must not begin until B1C passes.\n";
        return 1;
    }
    std::cout << "\nPROGRAM B1C GATE: PASS (" << report.num_passed() << "/"
              << (report.num_cases() - report.num_skipped()) << ")";
    if (report.num_skipped() > 0) std::cout << ", " << report.num_skipped() << " skipped";
    std::cout << "\n";
    return 0;
}
