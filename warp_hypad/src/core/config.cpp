#include "warp_hypad/core/config.hpp"

#include <cmath>
#include <set>
#include <sstream>

namespace warp_hypad {
namespace {

// Reads an optional number, leaving `out` untouched when absent.
//
// Every field is optional and falls back to the struct default. That is the
// right default policy here because the defaults are the documented reference
// configuration; a config naming only "metric.velocity" is a legitimate and
// readable way to say "the reference run, but faster".
bool get_number_opt(const JsonValue& obj, const char* key, double& out,
                    std::string& error) {
    const JsonValue* v = obj.find(key);
    if (!v) return true;
    if (!v->is_number()) {
        error = std::string("key '") + key + "' must be a number";
        return false;
    }
    out = v->as_number();
    return true;
}

bool get_int_opt(const JsonValue& obj, const char* key, int& out,
                 std::string& error) {
    const JsonValue* v = obj.find(key);
    if (!v) return true;
    if (!v->is_number()) {
        error = std::string("key '") + key + "' must be a number";
        return false;
    }
    const double d = v->as_number();
    if (d != std::floor(d)) {
        error = std::string("key '") + key + "' must be an integer";
        return false;
    }
    out = static_cast<int>(d);
    return true;
}

bool get_string_opt(const JsonValue& obj, const char* key, std::string& out,
                    std::string& error) {
    const JsonValue* v = obj.find(key);
    if (!v) return true;
    if (!v->is_string()) {
        error = std::string("key '") + key + "' must be a string";
        return false;
    }
    out = v->as_string();
    return true;
}

void collect_unknown(const JsonValue& obj, const char* section,
                     const std::set<std::string>& known,
                     std::vector<std::string>& warnings) {
    if (!obj.is_object()) return;
    for (const auto& kv : obj.members()) {
        if (known.find(kv.first) == known.end()) {
            std::ostringstream os;
            os << "unrecognised key '";
            if (section[0] != '\0') os << section << ".";
            os << kv.first << "' ignored";
            warnings.push_back(os.str());
        }
    }
}

}  // namespace

std::string to_string(MetricType t) {
    switch (t) {
        case MetricType::Minkowski: return "minkowski";
        case MetricType::Alcubierre: return "alcubierre";
    }
    return "unknown";
}

bool metric_type_from_string(const std::string& s, MetricType& out) {
    if (s == "minkowski") { out = MetricType::Minkowski; return true; }
    if (s == "alcubierre") { out = MetricType::Alcubierre; return true; }
    return false;
}

bool config_from_json(const JsonValue& root, Config& cfg, std::string& error,
                      std::vector<std::string>& warnings) {
    if (!root.is_object()) {
        error = "top-level JSON value must be an object";
        return false;
    }

    collect_unknown(root, "",
                    {"grid", "metric", "conventions", "derivatives", "objective"},
                    warnings);

    // --- grid ---
    if (const JsonValue* g = root.find("grid")) {
        if (!g->is_object()) { error = "'grid' must be an object"; return false; }
        double d;
        int n;
        n = cfg.grid.nx; if (!get_int_opt(*g, "nx", n, error)) return false; cfg.grid.nx = n;
        n = cfg.grid.ny; if (!get_int_opt(*g, "ny", n, error)) return false; cfg.grid.ny = n;
        n = cfg.grid.nz; if (!get_int_opt(*g, "nz", n, error)) return false; cfg.grid.nz = n;
        d = cfg.grid.x_min; if (!get_number_opt(*g, "x_min", d, error)) return false; cfg.grid.x_min = d;
        d = cfg.grid.x_max; if (!get_number_opt(*g, "x_max", d, error)) return false; cfg.grid.x_max = d;
        d = cfg.grid.y_min; if (!get_number_opt(*g, "y_min", d, error)) return false; cfg.grid.y_min = d;
        d = cfg.grid.y_max; if (!get_number_opt(*g, "y_max", d, error)) return false; cfg.grid.y_max = d;
        d = cfg.grid.z_min; if (!get_number_opt(*g, "z_min", d, error)) return false; cfg.grid.z_min = d;
        d = cfg.grid.z_max; if (!get_number_opt(*g, "z_max", d, error)) return false; cfg.grid.z_max = d;
        collect_unknown(*g, "grid",
                        {"nx", "ny", "nz", "x_min", "x_max", "y_min", "y_max",
                         "z_min", "z_max"},
                        warnings);
    }

    // --- metric ---
    if (const JsonValue* m = root.find("metric")) {
        if (!m->is_object()) { error = "'metric' must be an object"; return false; }
        std::string type_str = to_string(cfg.metric.type);
        if (!get_string_opt(*m, "type", type_str, error)) return false;
        if (!metric_type_from_string(type_str, cfg.metric.type)) {
            error = "metric.type must be one of: minkowski, alcubierre (got '" +
                    type_str + "')";
            return false;
        }
        if (!get_number_opt(*m, "velocity", cfg.metric.velocity, error)) return false;
        if (!get_number_opt(*m, "radius", cfg.metric.radius, error)) return false;
        if (!get_number_opt(*m, "sigma", cfg.metric.sigma, error)) return false;
        collect_unknown(*m, "metric", {"type", "velocity", "radius", "sigma"},
                        warnings);
    }

    // --- conventions ---
    if (const JsonValue* c = root.find("conventions")) {
        if (!c->is_object()) { error = "'conventions' must be an object"; return false; }
        if (!get_string_opt(*c, "metric_signature", cfg.conventions.metric_signature, error)) return false;
        if (!get_string_opt(*c, "units", cfg.conventions.units, error)) return false;
        if (!get_string_opt(*c, "riemann_convention", cfg.conventions.riemann_convention, error)) return false;
        collect_unknown(*c, "conventions",
                        {"metric_signature", "units", "riemann_convention"},
                        warnings);
    }

    // --- derivatives ---
    if (const JsonValue* d = root.find("derivatives")) {
        if (!d->is_object()) { error = "'derivatives' must be an object"; return false; }
        if (!get_int_opt(*d, "interior_order", cfg.derivatives.interior_order, error)) return false;
        if (!get_string_opt(*d, "boundary_scheme", cfg.derivatives.boundary_scheme, error)) return false;
        collect_unknown(*d, "derivatives", {"interior_order", "boundary_scheme"},
                        warnings);
    }

    // --- objective ---
    if (const JsonValue* o = root.find("objective")) {
        if (!o->is_object()) { error = "'objective' must be an object"; return false; }
        if (!get_string_opt(*o, "type", cfg.objective.type, error)) return false;
        if (!get_number_opt(*o, "epsilon", cfg.objective.epsilon, error)) return false;
        collect_unknown(*o, "objective", {"type", "epsilon"}, warnings);
    }

    return true;
}

JsonValue config_to_json(const Config& cfg) {
    JsonValue grid = JsonValue::object();
    grid.set("nx", JsonValue(cfg.grid.nx));
    grid.set("ny", JsonValue(cfg.grid.ny));
    grid.set("nz", JsonValue(cfg.grid.nz));
    grid.set("x_min", JsonValue(cfg.grid.x_min));
    grid.set("x_max", JsonValue(cfg.grid.x_max));
    grid.set("y_min", JsonValue(cfg.grid.y_min));
    grid.set("y_max", JsonValue(cfg.grid.y_max));
    grid.set("z_min", JsonValue(cfg.grid.z_min));
    grid.set("z_max", JsonValue(cfg.grid.z_max));

    JsonValue metric = JsonValue::object();
    metric.set("type", JsonValue(to_string(cfg.metric.type)));
    metric.set("velocity", JsonValue(cfg.metric.velocity));
    metric.set("radius", JsonValue(cfg.metric.radius));
    metric.set("sigma", JsonValue(cfg.metric.sigma));

    JsonValue conv = JsonValue::object();
    conv.set("metric_signature", JsonValue(cfg.conventions.metric_signature));
    conv.set("units", JsonValue(cfg.conventions.units));
    conv.set("riemann_convention", JsonValue(cfg.conventions.riemann_convention));

    JsonValue deriv = JsonValue::object();
    deriv.set("interior_order", JsonValue(cfg.derivatives.interior_order));
    deriv.set("boundary_scheme", JsonValue(cfg.derivatives.boundary_scheme));

    JsonValue obj = JsonValue::object();
    obj.set("type", JsonValue(cfg.objective.type));
    obj.set("epsilon", JsonValue(cfg.objective.epsilon));

    JsonValue root = JsonValue::object();
    root.set("grid", grid);
    root.set("metric", metric);
    root.set("conventions", conv);
    root.set("derivatives", deriv);
    root.set("objective", obj);
    return root;
}

bool config_validate(const Config& cfg, std::string& error) {
    if (!validate_grid_spec(cfg.grid, error)) return false;

    if (cfg.conventions.metric_signature != "-+++") {
        error = "only the '-+++' metric signature is implemented, got '" +
                cfg.conventions.metric_signature + "'";
        return false;
    }
    if (cfg.conventions.units != "geometrized") {
        error = "only 'geometrized' units are implemented, got '" +
                cfg.conventions.units + "'";
        return false;
    }
    if (cfg.derivatives.interior_order != 2) {
        error = "only interior_order 2 is implemented (Milestone 2)";
        return false;
    }
    if (cfg.derivatives.boundary_scheme != "one_sided") {
        error = "only the 'one_sided' boundary scheme is implemented, got '" +
                cfg.derivatives.boundary_scheme + "'";
        return false;
    }
    if (cfg.objective.type != "smooth_negative_energy") {
        error = "only the 'smooth_negative_energy' objective is implemented, got '" +
                cfg.objective.type + "'";
        return false;
    }
    if (!(cfg.objective.epsilon > 0.0) || !std::isfinite(cfg.objective.epsilon)) {
        error = "objective.epsilon must be finite and > 0";
        return false;
    }

    if (cfg.metric.type == MetricType::Alcubierre) {
        if (!(cfg.metric.radius > 0.0) || !std::isfinite(cfg.metric.radius)) {
            error = "metric.radius must be finite and > 0";
            return false;
        }
        if (!(cfg.metric.sigma > 0.0) || !std::isfinite(cfg.metric.sigma)) {
            error = "metric.sigma must be finite and > 0";
            return false;
        }
        if (!std::isfinite(cfg.metric.velocity)) {
            error = "metric.velocity must be finite";
            return false;
        }
    }
    return true;
}

bool config_load_file(const std::string& path, Config& cfg, std::string& error,
                      std::vector<std::string>& warnings) {
    JsonValue root;
    if (!json_read_file(path, root, error)) return false;
    if (!config_from_json(root, cfg, error, warnings)) return false;
    return config_validate(cfg, error);
}

}  // namespace warp_hypad
