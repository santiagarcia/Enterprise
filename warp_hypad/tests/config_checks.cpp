#include <string>
#include <vector>

#include "checks.hpp"
#include "warp_hypad/core/config.hpp"
#include "warp_hypad/core/json.hpp"

namespace warp_hypad {
namespace checks {
namespace {

// The reference configuration exactly as specified for the project. Kept here
// as a literal rather than read from disk so the check does not depend on the
// working directory, and so a change to config/example_alcubierre.json cannot
// silently weaken the test.
const char* kReferenceConfig = R"JSON({
  "grid": {
    "nx": 41,
    "ny": 41,
    "nz": 41,
    "x_min": -10.0,
    "x_max": 10.0,
    "y_min": -10.0,
    "y_max": 10.0,
    "z_min": -10.0,
    "z_max": 10.0
  },
  "metric": {
    "type": "alcubierre",
    "velocity": 0.2,
    "radius": 3.0,
    "sigma": 2.0
  },
  "conventions": {
    "metric_signature": "-+++",
    "units": "geometrized",
    "riemann_convention": "documented_internal"
  },
  "derivatives": {
    "interior_order": 2,
    "boundary_scheme": "one_sided"
  },
  "objective": {
    "type": "smooth_negative_energy",
    "epsilon": 1e-8
  }
})JSON";

}  // namespace

void run_config_checks(VerificationReport& report) {
    // --- The reference configuration parses and validates ---
    Config reference_cfg;
    {
        Stopwatch sw;
        JsonValue root;
        std::string err;
        std::vector<std::string> warnings;
        bool ok = json_parse(kReferenceConfig, root, err);
        if (ok) ok = config_from_json(root, reference_cfg, err, warnings);
        if (ok) ok = config_validate(reference_cfg, err);

        report.check_predicate(
            "config.reference_parses",
            "the specified reference configuration parses, validates and "
            "produces no warnings",
            VerificationLevel::A, ok && warnings.empty(),
            "parsed with 0 warnings",
            ok ? ("parsed with " + std::to_string(warnings.size()) + " warnings")
               : ("failed: " + err),
            sw.seconds());
    }

    // --- Values land where they belong ---
    //
    // Parsing "successfully" into the wrong fields is the failure mode a
    // round-trip test cannot see, because it round-trips consistently.
    {
        Stopwatch sw;
        int mismatches = 0;
        if (reference_cfg.grid.nx != 41) ++mismatches;
        if (reference_cfg.grid.ny != 41) ++mismatches;
        if (reference_cfg.grid.nz != 41) ++mismatches;
        if (reference_cfg.grid.x_min != -10.0) ++mismatches;
        if (reference_cfg.grid.x_max != 10.0) ++mismatches;
        if (reference_cfg.metric.type != MetricType::Alcubierre) ++mismatches;
        if (reference_cfg.metric.velocity != 0.2) ++mismatches;
        if (reference_cfg.metric.radius != 3.0) ++mismatches;
        if (reference_cfg.metric.sigma != 2.0) ++mismatches;
        if (reference_cfg.conventions.metric_signature != "-+++") ++mismatches;
        if (reference_cfg.conventions.units != "geometrized") ++mismatches;
        if (reference_cfg.derivatives.interior_order != 2) ++mismatches;
        if (reference_cfg.derivatives.boundary_scheme != "one_sided") ++mismatches;
        if (reference_cfg.objective.type != "smooth_negative_energy") ++mismatches;
        if (reference_cfg.objective.epsilon != 1e-8) ++mismatches;

        report.check_exact(
            "config.values_correct",
            "every field of the reference configuration parses to its "
            "specified value, exactly",
            VerificationLevel::A, "n/a", 0.0, static_cast<double>(mismatches),
            sw.seconds(),
            "catches fields that parse consistently into the wrong slot");
    }

    // --- Round-trip: config -> json -> text -> json -> config -> text ---
    {
        Stopwatch sw;
        const std::string first = json_dump(config_to_json(reference_cfg));

        JsonValue reparsed;
        std::string err;
        std::vector<std::string> warnings;
        Config cfg2;
        bool ok = json_parse(first, reparsed, err);
        if (ok) ok = config_from_json(reparsed, cfg2, err, warnings);
        const std::string second = ok ? json_dump(config_to_json(cfg2)) : std::string();

        report.check_predicate(
            "config.roundtrip_byte_exact",
            "dump(to_json(from_json(dump(to_json(cfg))))) reproduces the "
            "identical byte stream",
            VerificationLevel::A, ok && first == second && !first.empty(),
            "identical byte streams",
            ok ? (first == second ? "identical"
                                  : "differed; first " +
                                        std::to_string(first.size()) +
                                        " bytes, second " +
                                        std::to_string(second.size()))
               : ("failed: " + err),
            sw.seconds(),
            "a results directory must reproduce the run that made it");
    }

    // --- Numbers survive the round trip bit-exactly ---
    //
    // This is the property the shortest-round-trip number writer exists for.
    // 0.2 and 1e-8 are not representable in binary, so a naive %g writer at the
    // wrong precision loses them.
    {
        Stopwatch sw;
        Config cfg = reference_cfg;
        cfg.metric.velocity = 0.1 + 0.2;          // 0.30000000000000004
        cfg.metric.radius = 1.0 / 3.0;
        cfg.objective.epsilon = 4.9406564584124654e-300;
        cfg.grid.x_min = -1.7976931348623157e308;  // -DBL_MAX

        const std::string text = json_dump(config_to_json(cfg));
        JsonValue root;
        std::string err;
        std::vector<std::string> warnings;
        Config back;
        bool ok = json_parse(text, root, err) &&
                  config_from_json(root, back, err, warnings);

        int mismatches = 0;
        if (!ok) {
            mismatches = 99;
        } else {
            if (back.metric.velocity != cfg.metric.velocity) ++mismatches;
            if (back.metric.radius != cfg.metric.radius) ++mismatches;
            if (back.objective.epsilon != cfg.objective.epsilon) ++mismatches;
            if (back.grid.x_min != cfg.grid.x_min) ++mismatches;
        }

        report.check_exact(
            "config.number_roundtrip_bit_exact",
            "0.1+0.2, 1/3, a subnormal-scale value and -DBL_MAX all survive "
            "write-then-read bit for bit",
            VerificationLevel::A, "n/a", 0.0, static_cast<double>(mismatches),
            sw.seconds(),
            "the shortest-representation number writer exists for exactly this");
    }

    // --- Unknown keys are reported, not silently dropped ---
    {
        Stopwatch sw;
        const char* typo = R"JSON({
          "metric": { "type": "alcubierre", "sigmaa": 2.0 },
          "grid": { "nx": 11 },
          "unknown_section": { "a": 1 }
        })JSON";

        JsonValue root;
        std::string err;
        std::vector<std::string> warnings;
        Config cfg;
        const bool ok = json_parse(typo, root, err) &&
                        config_from_json(root, cfg, err, warnings);

        // Expect exactly two warnings: metric.sigmaa and the unknown section.
        const bool right = ok && warnings.size() == 2;

        std::string got;
        for (const std::string& w : warnings) {
            if (!got.empty()) got += "; ";
            got += w;
        }

        report.check_predicate(
            "config.unknown_keys_warn",
            "a misspelled key produces a warning rather than silently taking "
            "the default value",
            VerificationLevel::A, right, "2 warnings",
            ok ? (std::to_string(warnings.size()) + " warnings: " + got)
               : ("parse failed: " + err),
            sw.seconds(),
            "'sigmaa' running silently at the default sigma would produce "
            "results nobody could explain");
    }

    // --- Invalid values are rejected ---
    {
        Stopwatch sw;
        struct BadCase {
            const char* name;
            const char* json;
        };
        const BadCase bad[] = {
            {"unknown metric type", R"({"metric":{"type":"kerr"}})"},
            {"string where number expected", R"({"metric":{"radius":"three"}})"},
            {"non-integer grid size", R"({"grid":{"nx":40.5}})"},
            {"top level is an array", R"([1,2,3])"},
            {"malformed JSON", R"({"grid":{"nx":41,})"},
        };

        int accepted = 0;
        std::string accepted_names;
        for (const BadCase& b : bad) {
            JsonValue root;
            std::string err;
            std::vector<std::string> warnings;
            Config cfg;
            const bool ok = json_parse(b.json, root, err) &&
                            config_from_json(root, cfg, err, warnings);
            if (ok) {
                ++accepted;
                if (!accepted_names.empty()) accepted_names += ", ";
                accepted_names += b.name;
            }
        }

        report.check_predicate(
            "config.invalid_rejected",
            "unknown metric type, wrong value types, a non-integer grid size, "
            "a non-object root and malformed JSON are all rejected",
            VerificationLevel::A, accepted == 0, "0 of 5 accepted",
            accepted == 0 ? "0 of 5 accepted"
                          : (std::to_string(accepted) + " accepted: " + accepted_names),
            sw.seconds());
    }

    // --- Semantic validation rejects physically meaningless parameters ---
    {
        Stopwatch sw;
        int accepted = 0;
        std::string err;

        Config cfg = reference_cfg;
        cfg.metric.radius = -1.0;
        if (config_validate(cfg, err)) ++accepted;

        cfg = reference_cfg;
        cfg.metric.sigma = 0.0;
        if (config_validate(cfg, err)) ++accepted;

        cfg = reference_cfg;
        cfg.objective.epsilon = 0.0;
        if (config_validate(cfg, err)) ++accepted;

        cfg = reference_cfg;
        cfg.conventions.metric_signature = "+---";
        if (config_validate(cfg, err)) ++accepted;

        cfg = reference_cfg;
        cfg.derivatives.interior_order = 4;
        if (config_validate(cfg, err)) ++accepted;

        report.check_predicate(
            "config.validation_rejects_unsupported",
            "negative radius, zero sigma, zero epsilon, an unimplemented "
            "signature and an unimplemented derivative order are rejected",
            VerificationLevel::A, accepted == 0, "0 of 5 accepted",
            std::to_string(accepted) + " of 5 accepted", sw.seconds(),
            "an unimplemented option must fail loudly rather than be ignored");
    }

    // --- A config naming only some fields inherits documented defaults ---
    {
        Stopwatch sw;
        JsonValue root;
        std::string err;
        std::vector<std::string> warnings;
        Config cfg;
        const bool ok =
            json_parse(R"({"metric":{"velocity":0.5}})", root, err) &&
            config_from_json(root, cfg, err, warnings) &&
            config_validate(cfg, err);

        const bool right = ok && cfg.metric.velocity == 0.5 &&
                           cfg.metric.radius == 3.0 && cfg.grid.nx == 41 &&
                           warnings.empty();

        report.check_predicate(
            "config.partial_inherits_defaults",
            "a partial config overrides only what it names and inherits the "
            "reference defaults elsewhere",
            VerificationLevel::A, right, "velocity 0.5, radius 3, nx 41",
            ok ? ("velocity " + std::to_string(cfg.metric.velocity) +
                  ", radius " + std::to_string(cfg.metric.radius) + ", nx " +
                  std::to_string(cfg.grid.nx))
               : ("failed: " + err),
            sw.seconds());
    }

    // --- JSON parser edge cases ---
    {
        Stopwatch sw;
        struct Case {
            const char* text;
            bool should_parse;
        };
        const Case cases[] = {
            {R"({})", true},
            {R"({"a":[]})", true},
            {R"({"a":[1,2,3],"b":{"c":null},"d":true})", true},
            {R"({"s":"esc \" \\ \n \t done"})", true},
            {R"({"n":-1.5e-10})", true},
            {R"({"a":1,})", false},        // trailing comma
            {R"({"a":1} garbage)", false}, // trailing content
            {R"({a:1})", false},           // unquoted key
            {R"({"a":})", false},          // missing value
            {R"("unterminated)", false},
        };

        int wrong = 0;
        std::string detail;
        for (const Case& c : cases) {
            JsonValue v;
            std::string err;
            const bool ok = json_parse(c.text, v, err);
            if (ok != c.should_parse) {
                ++wrong;
                if (!detail.empty()) detail += "; ";
                detail += std::string(c.text) + " -> " + (ok ? "accepted" : "rejected");
            }
        }

        report.check_exact(
            "config.json_parser_edge_cases",
            "the parser accepts well-formed JSON and rejects trailing commas, "
            "trailing content, unquoted keys and unterminated strings",
            VerificationLevel::A, "n/a", 0.0, static_cast<double>(wrong),
            sw.seconds(), detail.empty() ? "all 10 cases behaved as specified" : detail);
    }
}

}  // namespace checks
}  // namespace warp_hypad
