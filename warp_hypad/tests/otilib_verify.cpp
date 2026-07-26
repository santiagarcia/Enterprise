#include <filesystem>
#include <iostream>
#include <string>

#include "otilib_checks.hpp"
#include "warp_hypad/core/metadata.hpp"
#include "warp_hypad/sensitivity/otilib_traits.hpp"

// ---------------------------------------------------------------------------
// OTIlib integration gate.
//
// Exit code 0 means the OTIlib adapter is verified against analytic
// derivatives, an independent numerical method, and the matrix-inverse
// derivative identity. Anything else means OTI must not be used for the
// spacetime pipeline yet.
// ---------------------------------------------------------------------------

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/otilib";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            std::cerr << "usage: otilib_verify [--out <dir>]\n";
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
        "otilib",
        "OTIlib onumm3n2 adapter: 3 imaginary directions, truncation order 2");

    checks::run_otilib_checks(report, out_dir);
    checks::run_otilib_fd_checks(report, out_dir);

    report.print(std::cout);

    std::string error;
    if (!report.write(out_dir, error)) {
        std::cerr << "failed to write verification reports: " << error << "\n";
        return 3;
    }

    // Record what was actually integrated, so the artefacts are self-describing.
    {
        JsonValue meta = JsonValue::object();
        meta.set("upstream_repository",
                 JsonValue("https://github.com/mauriaristi/otilib"));
        meta.set("upstream_commit",
                 JsonValue("de3f7fa2b1ba667dff054ad5ba1f6f41f4ba8b54"));
        meta.set("upstream_license", JsonValue("GPL-3.0"));
        meta.set("oti_type", JsonValue("onumm3n2_t"));
        meta.set("imaginary_directions",
                 JsonValue(static_cast<double>(sensitivity::kOtiDirections)));
        meta.set("truncation_order",
                 JsonValue(static_cast<double>(sensitivity::kOtiOrder)));
        meta.set("integration_route",
                 JsonValue("thin C++ wrapper over the C API; scalar-only unity "
                           "translation unit compiled natively with MSVC"));
        JsonValue params = JsonValue::array();
        for (int i = 0; i < sensitivity::kParamCount; ++i)
            params.push_back(JsonValue(sensitivity::parameter_name(i)));
        meta.set("derivative_channel_order", params);
        meta.set("build", build_info_json());
        meta.set("timestamp_utc", JsonValue(utc_timestamp()));

        if (!json_write_file(out_dir + "/otilib_metadata.json", meta, error)) {
            std::cerr << "failed to write otilib_metadata.json: " << error << "\n";
            return 3;
        }
    }

    if (!report.all_passed()) {
        std::cout << "\nOTILIB GATE: FAIL (" << report.num_failed() << " of "
                  << report.num_cases() << " cases failed)\n"
                  << "OTI must not be used for the spacetime pipeline until "
                     "these pass.\n";
        return 1;
    }

    std::cout << "\nOTILIB GATE: PASS (" << report.num_passed() << "/"
              << report.num_cases() << ")\n";
    return 0;
}
