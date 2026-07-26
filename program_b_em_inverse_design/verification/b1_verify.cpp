#include <filesystem>
#include <iostream>
#include <string>

#include "b1_checks.hpp"
#include "warp_hypad/core/metadata.hpp"

// Program B Milestone B1 gate: metric <-> electromagnetic constitutive mapping.

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/program_b1";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            std::cerr << "usage: b1_verify [--out <dir>]\n";
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
        "program_b1",
        "Metric to electromagnetic constitutive response: Plebanski-Tamm "
        "forward map, conformal-class inverse, and non-metric detection");

    program_b::checks::run_b1_checks(report, out_dir);
    program_b::checks::run_b1_sensitivity_checks(report, out_dir);

    report.print(std::cout);

    std::string error;
    if (!report.write(out_dir, error)) {
        std::cerr << "failed to write verification reports: " << error << "\n";
        return 3;
    }

    {
        JsonValue meta = JsonValue::object();
        meta.set("milestone", JsonValue("B1"));
        meta.set("units", JsonValue("Heaviside-Lorentz, c = 1; vacuum is eps = mu = I"));
        meta.set("metric_signature", JsonValue("-+++"));
        meta.set("forward_map", JsonValue("Plebanski-Tamm: eps^ij = mu^ij = "
                                          "-sqrt(-g) g^ij / g_00, w_i = -g_0i/g_00"));
        meta.set("inverse_map",
                 JsonValue("closed form h = det(eps) eps^-1 in the g_00 = -1 "
                           "gauge, valid for w = 0 only"));
        meta.set("conformal_ambiguity",
                 JsonValue("the forward map is invariant under g -> Omega^2 g, "
                           "so the inverse returns a conformal class "
                           "representative and the gauge is reported"));
        meta.set("validity_condition",
                 JsonValue("g_00 < 0 required; fails where v_s f > 1"));
        meta.set("build", build_info_json());
        meta.set("timestamp_utc", JsonValue(utc_timestamp()));
        if (!json_write_file(out_dir + "/b1_metadata.json", meta, error)) {
            std::cerr << "failed to write b1_metadata.json: " << error << "\n";
            return 3;
        }
    }

    if (!report.all_passed()) {
        std::cout << "\nPROGRAM B1 GATE: FAIL (" << report.num_failed() << " of "
                  << (report.num_cases() - report.num_skipped())
                  << " cases failed)\n";
        return 1;
    }
    std::cout << "\nPROGRAM B1 GATE: PASS (" << report.num_passed() << "/"
              << (report.num_cases() - report.num_skipped()) << ")";
    if (report.num_skipped() > 0)
        std::cout << ", " << report.num_skipped() << " skipped";
    std::cout << "\n";
    return 0;
}
