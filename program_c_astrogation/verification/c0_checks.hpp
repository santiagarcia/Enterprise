#pragma once

#include <cstdint>
#include <string>

#include "warp_hypad/core/verification.hpp"

namespace program_c {
namespace checks {

struct C0Metrics {
    std::string data_kind;
    std::string catalogue_release;
    std::string adql_query;

    // The ingestion ledger of the delivered snapshot. These four must satisfy
    // rows_in == accepted + rejected + duplicates_collapsed, which is what makes
    // "0 silently dropped records" a checkable statement.
    int rows_in = 0;
    int accepted = 0;
    int rejected = 0;
    int duplicates_collapsed = 0;
    // A SEPARATE ledger, from the injected-duplicate case only. It is reported
    // apart from the four above because mixing the two ledgers would print a
    // reconciliation that does not balance.
    int injected_duplicate_rows_in = 0;
    int injected_duplicate_collapsed = 0;

    double worst_unit_roundtrip = 0.0;
    double worst_cartesian_roundtrip = 0.0;
    double proper_motion_convention_error = 0.0;
    double stationary_displacement_pc = 0.0;
    double high_pm_displacement_pc = 0.0;
    double light_time_residual_years = 0.0;
    int light_time_iterations = 0;
    double position_sigma_pc = 0.0;
    double velocity_sigma_pc_per_year = 0.0;
    double uncertainty_growth_ratio = 0.0;
    double snapshot_roundtrip_error = 0.0;
    double literature_sanity_error = 0.0;

    std::uint64_t checksum = 0;
    int count_within_5pc = 0;
    int count_within_50pc = 0;
    int exoplanet_hosts = 0;
    int frames_described = 0;
};

void run_c0_checks(warp_hypad::VerificationReport& report,
                   const std::string& out_dir, C0Metrics& metrics);

}  // namespace checks
}  // namespace program_c
