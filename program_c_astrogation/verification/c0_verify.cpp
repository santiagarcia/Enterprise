#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>

#include "c0_checks.hpp"
#include "program_c/catalogue.hpp"
#include "warp_hypad/core/metadata.hpp"

// Gate C0-CATALOGUE-UNIVERSE.
//
// SCOPE. This gate verifies the MACHINERY of a catalogue-grounded star database
// against SYNTHETIC fixtures whose correct answers are known in closed form. That
// is deliberate and is the right way to verify machinery; it says nothing about
// agreement with the sky.
//
// Agreement with the sky is the SEPARATE C0-REAL-CATALOGUE gate (c0_real_verify),
// which reads the committed Gaia DR3 snapshots under snapshots/ and checks them
// against independently published astronomy. Both gates are offline, deterministic
// and required. Retrieval itself happens only in tools/import_gaia.py, which no
// gate invokes.

int main(int argc, char** argv) {
    using namespace warp_hypad;

    std::string out_dir = "results/verification/program_c0";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) out_dir = argv[++i];
        else { std::cerr << "usage: c0_verify [--out <dir>]\n"; return 2; }
    }

    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if (ec) { std::cerr << "cannot create " << out_dir << "\n"; return 3; }

    VerificationReport report(
        "program_c0",
        "C0-CATALOGUE-UNIVERSE: astronomical data schema, catalogue ingestion, "
        "typed coordinate frames, epoch propagation, light time, uncertainty and "
        "a deterministic offline snapshot");

    program_c::checks::C0Metrics m;
    program_c::checks::run_c0_checks(report, out_dir, m);

    report.print(std::cout);

    std::string error;
    if (!report.write(out_dir, error)) {
        std::cerr << "failed to write reports: " << error << "\n"; return 3;
    }

    {
        JsonValue meta = JsonValue::object();
        meta.set("milestone", JsonValue("C0"));
        meta.set("data_kind", JsonValue(m.data_kind));
        meta.set("data_honesty",
                 JsonValue("Every fixture in THIS gate is SYNTHETIC or "
                           "LITERATURE_FROM_MEMORY and is tagged as such, by "
                           "design: closed-form fixtures verify the machinery "
                           "because their expected answers are exact. This gate "
                           "therefore makes no claim of agreement with Gaia DR3. "
                           "That claim belongs to the separate C0-REAL-CATALOGUE "
                           "gate, which reads the committed real snapshots under "
                           "snapshots/."));
        meta.set("real_data_gate", JsonValue("program_c0_real"));
        meta.set("pinned_catalogue", JsonValue(program_c::kPinnedCatalogue));
        meta.set("pinned_table", JsonValue(program_c::kPinnedTable));
        meta.set("reference_epoch_jyear",
                 JsonValue(program_c::kGaiaDR3ReferenceEpochJyear));
        meta.set("tap_endpoint", JsonValue(program_c::kGaiaTapEndpoint));
        meta.set("adql_query", JsonValue(m.adql_query));
        // The ingestion ledger, machine-readable, so the "0 silently dropped"
        // claim can be rechecked without rerunning the gate.
        meta.set("rows_in", JsonValue(static_cast<double>(m.rows_in)));
        meta.set("rows_accepted", JsonValue(static_cast<double>(m.accepted)));
        meta.set("rows_rejected", JsonValue(static_cast<double>(m.rejected)));
        meta.set("duplicates_collapsed",
                 JsonValue(static_cast<double>(m.duplicates_collapsed)));
        meta.set("ledger_balances",
                 JsonValue(m.rows_in ==
                           m.accepted + m.rejected + m.duplicates_collapsed));
        meta.set("snapshot_checksum_fnv1a",
                 JsonValue(static_cast<double>(m.checksum)));
        meta.set("database_format",
                 JsonValue("portable versioned text snapshot; SQLite deferred as "
                           "a dependency decision, materialised by "
                           "tools/snapshot_to_sqlite.py"));
        meta.set("propagation_model",
                 JsonValue("rectilinear barycentric, a LOCAL KINEMATIC model, not "
                           "a Galactic orbit integrator; validity 1e5 yr"));
        meta.set("aberration", JsonValue("absent: relativistic aberration is C1"));
        meta.set("build", build_info_json());
        meta.set("timestamp_utc", JsonValue(utc_timestamp()));
        if (!json_write_file(out_dir + "/c0_metadata.json", meta, error)) return 3;
    }

    std::cout << "\nC0 REQUIRED QUANTITIES\n";
    std::cout << "  data kind                          " << m.data_kind
              << "   (NOT retrieved catalogue data)\n";
    std::cout << "  pinned release                     " << m.catalogue_release << "\n";
    std::cout << "  rows in / accepted / rejected      " << m.rows_in << " / "
              << m.accepted << " / " << m.rejected << "\n";
    std::cout << "  duplicates collapsed               " << m.duplicates_collapsed
              << "  (ledger balances: " << m.rows_in << " = " << m.accepted << " + "
              << m.rejected << " + " << m.duplicates_collapsed << ")\n";
    std::cout << "  injected-duplicate case            "
              << m.injected_duplicate_rows_in << " in, "
              << m.injected_duplicate_collapsed
              << " collapsed  (a SEPARATE ledger, not the delivered snapshot)\n";
    std::cout << "  frames fully described             " << m.frames_described << "\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  worst unit round trip              " << m.worst_unit_roundtrip << "\n";
    std::cout << "  worst Cartesian round trip         " << m.worst_cartesian_roundtrip << "\n";
    std::cout << "  proper-motion convention error     " << m.proper_motion_convention_error << "\n";
    std::cout << "  stationary displacement            " << m.stationary_displacement_pc
              << " pc  (must be exactly 0)\n";
    std::cout << "  high-pm displacement, 1000 yr      " << m.high_pm_displacement_pc << " pc\n";
    std::cout << "  light-time residual                " << m.light_time_residual_years
              << " yr in " << m.light_time_iterations << " iterations\n";
    std::cout << "  position / velocity sigma          " << m.position_sigma_pc << " pc / "
              << m.velocity_sigma_pc_per_year << " pc/yr\n";
    std::cout << "  uncertainty growth at 1e4 yr       " << m.uncertainty_growth_ratio << "x\n";
    std::cout << "  snapshot round-trip error          " << m.snapshot_roundtrip_error << "\n";
    std::cout << "  literature sanity (loose)          " << m.literature_sanity_error << "\n";
    std::cout << "  within 5 pc / 50 pc / exo hosts    " << m.count_within_5pc << " / "
              << m.count_within_50pc << " / " << m.exoplanet_hosts << "\n";
    std::cout << "  snapshot checksum (FNV-1a)         " << m.checksum << "\n";

    std::cout << "\nSCOPE: C0 is a catalogue and coordinate layer. It does not "
                 "simulate a warp\n"
                 "  bubble, optimize a voyage, render a cockpit, or claim "
                 "Galactic dynamics.\n"
                 "  THIS gate's fixtures are SYNTHETIC by design, so it verifies "
                 "the machinery\n"
                 "  and claims nothing about the sky. Agreement with real Gaia DR3 "
                 "records is\n"
                 "  the separate C0-REAL-CATALOGUE gate (c0_real_verify).\n";

    if (!report.all_passed()) {
        std::cout << "\nC0 GATE: FAIL (" << report.num_failed() << " of "
                  << (report.num_cases() - report.num_skipped()) << " cases failed)\n";
        return 1;
    }
    std::cout << "\nC0 GATE: PASS (" << report.num_passed() << "/"
              << (report.num_cases() - report.num_skipped()) << ")";
    if (report.num_skipped() > 0) std::cout << ", " << report.num_skipped() << " skipped";
    std::cout << "\n";
    return 0;
}
