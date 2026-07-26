#pragma once

#include <string>

#include "warp_hypad/core/config.hpp"
#include "warp_hypad/core/grid.hpp"
#include "warp_hypad/core/json.hpp"

// ---------------------------------------------------------------------------
// Run provenance.
//
// The requirement is that every figure be reproducible from saved output
// without rerunning the solver, and that every verification case record the
// commit, compiler and build mode it was produced by. That only works if the
// provenance is captured by the build system rather than typed in by hand, so
// the values below are injected as compile definitions from CMake and baked
// into the binary.
//
// If CMake could not determine a value (no git, or a source tree exported
// without history) the field reads "unknown" rather than being omitted. A
// missing provenance field must be visible, because a results directory that
// cannot be traced to a commit is not a scientific record.
// ---------------------------------------------------------------------------

namespace warp_hypad {

struct BuildInfo {
    std::string git_commit;       // full SHA, or "unknown"
    bool git_dirty = false;       // uncommitted changes present at configure time
    std::string compiler_id;      // e.g. "MSVC", "GNU", "Clang"
    std::string compiler_version;
    std::string build_type;       // e.g. "Release", "RelWithDebInfo"
    std::string cxx_standard;
    std::string platform;
};

const BuildInfo& build_info();

// ISO-8601 UTC, e.g. "2026-07-24T18:41:07Z". Captured at call time.
std::string utc_timestamp();

JsonValue build_info_json();
JsonValue grid_json(const Grid& grid);

// The full metadata.json payload for a run: build info, timestamp, resolved
// configuration, derived grid quantities and the conventions in force.
JsonValue make_metadata(const Config& cfg, const Grid& grid,
                        const std::string& run_name);

}  // namespace warp_hypad
