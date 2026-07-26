#pragma once

#include <string>
#include <vector>

#include "warp_hypad/core/grid.hpp"
#include "warp_hypad/core/json.hpp"

// ---------------------------------------------------------------------------
// The run configuration contract.
//
// The schema is exactly the one specified for the project and nothing more.
// Two properties are enforced deliberately:
//
//   1. Round-trip fidelity. to_json(from_json(x)) == x byte for byte, so a
//      config recovered from a results directory reproduces the run.
//
//   2. Unknown keys are reported, not ignored. A silently dropped "sigmaa"
//      would run happily with the default and produce results nobody could
//      explain. Unknown keys are returned as warnings rather than hard errors
//      so that a config written for a later milestone still loads today.
// ---------------------------------------------------------------------------

namespace warp_hypad {

enum class MetricType { Minkowski, Alcubierre };

struct MetricConfig {
    MetricType type = MetricType::Alcubierre;

    // Alcubierre-type parameters. Ignored when type == Minkowski.
    // These are the quantities that will carry HyPAD derivative seeds.
    double velocity = 0.2;  // bubble velocity v, in units of c
    double radius = 3.0;    // bubble radius R
    double sigma = 2.0;     // wall steepness
};

struct ConventionsConfig {
    // Only "-+++" is implemented. Recorded explicitly so that a results
    // directory is self-describing and cross-code comparison in Milestone 10
    // cannot silently compare mismatched signatures.
    std::string metric_signature = "-+++";
    std::string units = "geometrized";
    std::string riemann_convention = "documented_internal";
};

struct DerivativeConfig {
    int interior_order = 2;
    std::string boundary_scheme = "one_sided";
};

struct ObjectiveConfig {
    std::string type = "smooth_negative_energy";
    double epsilon = 1e-8;
};

struct Config {
    GridSpec grid;
    MetricConfig metric;
    ConventionsConfig conventions;
    DerivativeConfig derivatives;
    ObjectiveConfig objective;
};

std::string to_string(MetricType t);
bool metric_type_from_string(const std::string& s, MetricType& out);

// Parses `root` into `cfg`. Returns false and fills `error` on a hard failure
// (wrong type, unsupported value, missing required section). Unrecognised keys
// are appended to `warnings` and do not fail the load.
bool config_from_json(const JsonValue& root, Config& cfg, std::string& error,
                      std::vector<std::string>& warnings);

JsonValue config_to_json(const Config& cfg);

// Semantic validation beyond what parsing checks: positive radius and sigma,
// a grid that can actually be allocated, a usable epsilon.
bool config_validate(const Config& cfg, std::string& error);

bool config_load_file(const std::string& path, Config& cfg, std::string& error,
                      std::vector<std::string>& warnings);

}  // namespace warp_hypad
