#include <filesystem>
#include <iostream>
#include <string>

#include "metric_checks.hpp"
#include "warp_hypad/core/metadata.hpp"
#include "warp_hypad/metrics/metric.hpp"

// Milestone A3 gate: metric generation.
//
// Exit 0 only if the Minkowski and Alcubierre metrics satisfy every structural
// identity AND (when the OTIlib backend is enabled) their parameter
// sensitivities match analytic derivatives and finite differences.

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/milestone_03";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            std::cerr << "usage: milestone3_verify [--out <dir>]\n";
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
        "milestone_03",
        "Metric generation: Minkowski and frozen-time Alcubierre, with "
        "determinant, inverse, signature and HyPAD parameter sensitivities");

    checks::run_metric_checks(report, out_dir);
    checks::run_metric_sensitivity_checks(report, out_dir);

    report.print(std::cout);

    std::string error;
    if (!report.write(out_dir, error)) {
        std::cerr << "failed to write verification reports: " << error << "\n";
        return 3;
    }

    // Record the time treatment in the results directory. A reader must not be
    // able to mistake these outputs for a full moving-bubble spacetime.
    {
        JsonValue meta = JsonValue::object();
        meta.set("milestone", JsonValue("A3"));
        meta.set("metric", JsonValue("alcubierre"));
        meta.set("time_treatment",
                 JsonValue(metrics::to_string(metrics::TimeTreatment::FrozenTimeDiagnostic)));
        meta.set("time_treatment_meaning",
                 JsonValue("The bubble centre is held at a fixed coordinate and "
                           "d_t g_{mu nu} == 0 by construction. These results "
                           "are a SPATIAL SLICE of a moving spacetime, not the "
                           "moving spacetime. Curvature computed from this "
                           "slice alone is NOT the complete physical Alcubierre "
                           "spacetime curvature."));
        meta.set("temporal_derivatives", JsonValue("set to zero, explicitly"));
        meta.set("future_path",
                 JsonValue("metric.evaluate(t, x, y, z, params) already accepts "
                           "t; TimeTreatment::MovingBubble is the unimplemented "
                           "branch that will use x_s(t)"));
        meta.set("parameters_differentiated",
                 JsonValue("velocity, radius, sigma"));
        meta.set("bubble_center_x",
                 JsonValue("fixed coordinate, not a differentiated parameter"));
        meta.set("build", build_info_json());
        meta.set("timestamp_utc", JsonValue(utc_timestamp()));
        if (!json_write_file(out_dir + "/metric_metadata.json", meta, error)) {
            std::cerr << "failed to write metric_metadata.json: " << error << "\n";
            return 3;
        }
    }

    std::cout << "\nTime treatment: FROZEN-TIME SPATIAL SLICE.\n"
              << "  d_t g_{mu nu} == 0 by construction. This is not the "
                 "complete Alcubierre spacetime.\n";

    if (!report.all_passed()) {
        std::cout << "\nMILESTONE A3 GATE: FAIL (" << report.num_failed() << " of "
                  << (report.num_cases() - report.num_skipped())
                  << " cases failed)\n"
                  << "Christoffel and curvature work must not begin until these "
                     "pass.\n";
        return 1;
    }

    std::cout << "\nMILESTONE A3 GATE: PASS (" << report.num_passed() << "/"
              << (report.num_cases() - report.num_skipped()) << ")";
    if (report.num_skipped() > 0)
        std::cout << ", " << report.num_skipped() << " skipped";
    std::cout << "\n";
    return 0;
}
