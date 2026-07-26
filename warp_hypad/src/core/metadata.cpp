#include "warp_hypad/core/metadata.hpp"

#include <ctime>

// Injected by CMake. Defaults keep the file compilable standalone, but a build
// that reaches here without them is a misconfigured build and says so.
#ifndef WARP_HYPAD_GIT_COMMIT
#define WARP_HYPAD_GIT_COMMIT "unknown"
#endif
#ifndef WARP_HYPAD_GIT_DIRTY
#define WARP_HYPAD_GIT_DIRTY 0
#endif
#ifndef WARP_HYPAD_COMPILER_ID
#define WARP_HYPAD_COMPILER_ID "unknown"
#endif
#ifndef WARP_HYPAD_COMPILER_VERSION
#define WARP_HYPAD_COMPILER_VERSION "unknown"
#endif
#ifndef WARP_HYPAD_BUILD_TYPE
#define WARP_HYPAD_BUILD_TYPE "unknown"
#endif
#ifndef WARP_HYPAD_PLATFORM
#define WARP_HYPAD_PLATFORM "unknown"
#endif

namespace warp_hypad {
namespace {

std::string cxx_standard_string() {
    // __cplusplus is unreliable on MSVC without /Zc:__cplusplus, which the
    // build sets; _MSVC_LANG is the fallback that is always correct there.
#if defined(_MSVC_LANG)
    constexpr long v = _MSVC_LANG;
#else
    constexpr long v = __cplusplus;
#endif
    // if constexpr, not if: the operand is a compile-time constant, and MSVC
    // /W4 correctly flags a runtime branch on one.
    if constexpr (v >= 202302L) return "C++23";
    else if constexpr (v >= 202002L) return "C++20";
    else if constexpr (v >= 201703L) return "C++17";
    else return "pre-C++17";
}

BuildInfo make_build_info() {
    BuildInfo b;
    b.git_commit = WARP_HYPAD_GIT_COMMIT;
    b.git_dirty = (WARP_HYPAD_GIT_DIRTY != 0);
    b.compiler_id = WARP_HYPAD_COMPILER_ID;
    b.compiler_version = WARP_HYPAD_COMPILER_VERSION;
    b.build_type = WARP_HYPAD_BUILD_TYPE;
    b.platform = WARP_HYPAD_PLATFORM;
    b.cxx_standard = cxx_standard_string();
    return b;
}

}  // namespace

const BuildInfo& build_info() {
    static const BuildInfo info = make_build_info();
    return info;
}

std::string utc_timestamp() {
    const std::time_t now = std::time(nullptr);
    std::tm tm_utc{};
#if defined(_WIN32)
    gmtime_s(&tm_utc, &now);
#else
    gmtime_r(&now, &tm_utc);
#endif
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm_utc);
    return std::string(buf);
}

JsonValue build_info_json() {
    const BuildInfo& b = build_info();
    JsonValue v = JsonValue::object();
    v.set("git_commit", JsonValue(b.git_commit));
    v.set("git_dirty", JsonValue(b.git_dirty));
    v.set("compiler_id", JsonValue(b.compiler_id));
    v.set("compiler_version", JsonValue(b.compiler_version));
    v.set("build_type", JsonValue(b.build_type));
    v.set("cxx_standard", JsonValue(b.cxx_standard));
    v.set("platform", JsonValue(b.platform));
    return v;
}

JsonValue grid_json(const Grid& grid) {
    JsonValue v = JsonValue::object();
    v.set("nx", JsonValue(grid.nx()));
    v.set("ny", JsonValue(grid.ny()));
    v.set("nz", JsonValue(grid.nz()));
    v.set("num_points", JsonValue(static_cast<double>(grid.num_points())));
    v.set("dx", JsonValue(grid.dx()));
    v.set("dy", JsonValue(grid.dy()));
    v.set("dz", JsonValue(grid.dz()));
    v.set("cell_volume", JsonValue(grid.cell_volume()));

    JsonValue bounds = JsonValue::object();
    bounds.set("x_min", JsonValue(grid.spec().x_min));
    bounds.set("x_max", JsonValue(grid.spec().x_max));
    bounds.set("y_min", JsonValue(grid.spec().y_min));
    bounds.set("y_max", JsonValue(grid.spec().y_max));
    bounds.set("z_min", JsonValue(grid.spec().z_min));
    bounds.set("z_max", JsonValue(grid.spec().z_max));
    v.set("bounds", bounds);
    return v;
}

JsonValue make_metadata(const Config& cfg, const Grid& grid,
                        const std::string& run_name) {
    JsonValue v = JsonValue::object();
    v.set("run_name", JsonValue(run_name));
    v.set("timestamp_utc", JsonValue(utc_timestamp()));
    v.set("build", build_info_json());
    v.set("config", config_to_json(cfg));
    v.set("grid_derived", grid_json(grid));

    // Layout facts a downstream reader needs in order to interpret raw field
    // dumps without reading the C++ source.
    JsonValue layout = JsonValue::object();
    layout.set("grid_index_order", JsonValue("flat = (i*ny + j)*nz + k"));
    layout.set("tensor_layout", JsonValue("point-major; component index fastest"));
    layout.set("rank2_offset", JsonValue("16*p + 4*mu + nu"));
    layout.set("rank3_offset", JsonValue("64*p + 16*a + 4*b + c"));
    layout.set("rank4_offset", JsonValue("256*p + 64*a + 16*b + 4*c + d"));
    layout.set("coordinate_index_map", JsonValue("0=t, 1=x, 2=y, 3=z"));
    v.set("layout", layout);

    v.set("milestone_status",
          JsonValue("Milestone 1 complete; Milestones 2-10 not implemented"));
    return v;
}

}  // namespace warp_hypad
