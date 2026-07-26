#include "program_b/effective_metric/metric_medium_map.hpp"

namespace program_b {

const char* to_string(ConformalGauge g) {
    switch (g) {
        case ConformalGauge::UnitLapse: return "unit_lapse_g00_eq_minus_one";
    }
    return "unknown";
}

const char* to_string(InverseStatus s) {
    switch (s) {
        case InverseStatus::Ok: return "ok";
        case InverseStatus::NotMetricRepresentable: return "not_metric_representable";
        case InverseStatus::MagnetoelectricUnsupported:
            return "magnetoelectric_unsupported_in_B1";
        case InverseStatus::Degenerate: return "degenerate";
    }
    return "unknown";
}

}  // namespace program_b
