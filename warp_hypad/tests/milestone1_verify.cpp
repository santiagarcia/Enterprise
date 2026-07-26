#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "checks.hpp"
#include "warp_hypad/core/config.hpp"
#include "warp_hypad/core/grid.hpp"
#include "warp_hypad/core/metadata.hpp"

// ---------------------------------------------------------------------------
// The Milestone 1 acceptance gate.
//
// Exit code 0 means every required verification passed. Any other exit code
// means the milestone is NOT complete, and by project policy Milestone 2 may
// not begin. That is the whole point of this executable: the gate is a process
// that returns a status, not a judgement someone makes by reading output.
//
// Usage:
//   milestone1_verify [--out <dir>] [--group <name>] [--quiet]
//
// --group runs a single verification group, which is how CTest registers each
// group as an individually reportable test. Report files are only written when
// the full suite runs, so a partial run cannot overwrite a complete record
// with a subset of it.
// ---------------------------------------------------------------------------

namespace {

void print_usage() {
    std::cout
        << "usage: milestone1_verify [--out <dir>] [--group <name>] [--quiet]\n"
           "  --out    directory for verification reports "
           "(default: results/verification/milestone_01)\n"
           "  --group  one of: grid, layout, algebra, contraction, genericity, "
           "config, all (default: all)\n"
           "  --quiet  suppress the per-case console listing\n";
}

}  // namespace

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/milestone_01";
    std::string group = "all";
    bool quiet = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else if (arg == "--group" && i + 1 < argc) {
            group = argv[++i];
        } else if (arg == "--quiet") {
            quiet = true;
        } else if (arg == "--help" || arg == "-h") {
            print_usage();
            return 0;
        } else {
            std::cerr << "unrecognised argument: " << arg << "\n";
            print_usage();
            return 2;
        }
    }

    const bool run_all = (group == "all");
    VerificationReport report(
        "milestone_01",
        "Grid, contiguous storage, fixed-size tensors and 4x4 matrix algebra");

    if (run_all || group == "grid") checks::run_grid_checks(report);
    if (run_all || group == "layout") checks::run_layout_checks(report);
    if (run_all || group == "algebra") checks::run_algebra_checks(report, out_dir);
    if (run_all || group == "contraction") checks::run_contraction_checks(report);
    if (run_all || group == "genericity") checks::run_genericity_checks(report);
    if (run_all || group == "config") checks::run_config_checks(report);

    if (report.num_cases() == 0) {
        std::cerr << "no verification group matched '" << group << "'\n";
        print_usage();
        return 2;
    }

    if (!quiet) report.print(std::cout);

    // Only a complete run may write the milestone record.
    if (run_all) {
        std::string error;
        if (!report.write(out_dir, error)) {
            std::cerr << "failed to write verification reports: " << error << "\n";
            return 3;
        }

        // The metadata-export prototype required by Milestone 1. It is written
        // from the reference configuration so that the exported layout and
        // convention description can be checked against the code by a reader
        // who has only the results directory.
        Config cfg;
        const Grid grid(cfg.grid);
        JsonValue meta = make_metadata(cfg, grid, "milestone_01_verification");
        if (!json_write_file(out_dir + "/metadata.json", meta, error)) {
            std::cerr << "failed to write metadata.json: " << error << "\n";
            return 3;
        }

        std::cout << "\nwrote:\n"
                  << "  " << out_dir << "/verification_summary.json\n"
                  << "  " << out_dir << "/verification_table.csv\n"
                  << "  " << out_dir << "/verification_report.md\n"
                  << "  " << out_dir << "/metadata.json\n"
                  << "  " << out_dir << "/inverse_residuals.csv\n";

        const BuildInfo& b = build_info();
        if (b.git_commit == "unknown") {
            std::cout << "\nWARNING: git commit could not be determined at "
                         "configure time; the verification record is not "
                         "traceable to a revision.\n";
        }
        if (b.git_dirty) {
            std::cout << "\nWARNING: the working tree had uncommitted changes "
                         "at configure time; results are not reproducible from "
                         "the recorded commit alone.\n";
        }
    }

    if (!report.all_passed()) {
        std::cout << "\nMILESTONE 1 GATE: FAIL (" << report.num_failed()
                  << " of " << report.num_cases() << " cases failed)\n"
                  << "Milestone 2 must not begin until these pass.\n";
        return 1;
    }

    std::cout << "\nMILESTONE 1 GATE: PASS (" << report.num_passed() << "/"
              << report.num_cases() << ")\n";
    return 0;
}
