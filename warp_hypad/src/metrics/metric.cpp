#include "warp_hypad/metrics/metric.hpp"

#include <cmath>

namespace warp_hypad {
namespace metrics {

const char* to_string(TimeTreatment t) {
    switch (t) {
        case TimeTreatment::FrozenTimeDiagnostic:
            return "frozen_time_spatial_slice_DIAGNOSTIC_ONLY";
        case TimeTreatment::InertialTimeDependent:
            return "inertial_frame_time_dependent_xs_of_t";
        case TimeTreatment::CoMoving:
            return "comoving_frame_NOT_IMPLEMENTED";
    }
    return "unknown";
}

// Closed-form derivatives of
//
//   f = (A - B) / D,   A = tanh(s(r+R)), B = tanh(s(r-R)), D = 2 tanh(sR)
//
// with sech^2(u) = 1 - tanh^2(u):
//
//   dA/dR = s sech^2(s(r+R))         dA/ds = (r+R) sech^2(s(r+R))
//   dB/dR = -s sech^2(s(r-R))        dB/ds = (r-R) sech^2(s(r-R))
//   dD/dR = 2s sech^2(sR)            dD/ds = 2R sech^2(sR)
//
// then the quotient rule. Derived by hand from the definition rather than by
// differentiating the templated implementation, so that comparing this against
// the HyPAD result is a genuine independent check rather than a tautology.
const char* to_string(MetricSensitivityMode m) {
    switch (m) {
        case MetricSensitivityMode::FixedBubbleCenter:
            return "fixed_bubble_center_shape_sensitivity";
        case MetricSensitivityMode::FixedLaboratoryCoordinates:
            return "fixed_laboratory_coordinates_total_sensitivity";
    }
    return "unknown";
}

ShapeDerivatives shape_function_derivatives(double r_s, double R, double sigma) {
    const double ta = std::tanh(sigma * (r_s + R));
    const double tb = std::tanh(sigma * (r_s - R));
    const double tr = std::tanh(sigma * R);

    const double sa = 1.0 - ta * ta;  // sech^2(sigma (r+R))
    const double sb = 1.0 - tb * tb;  // sech^2(sigma (r-R))
    const double sr = 1.0 - tr * tr;  // sech^2(sigma R)

    const double num = ta - tb;
    const double den = 2.0 * tr;

    ShapeDerivatives out;
    out.f = num / den;

    const double dnum_dR = sigma * sa + sigma * sb;
    const double dden_dR = 2.0 * sigma * sr;
    out.df_dR = (dnum_dR * den - num * dden_dR) / (den * den);

    const double dnum_ds = (r_s + R) * sa - (r_s - R) * sb;
    const double dden_ds = 2.0 * R * sr;
    out.df_dsigma = (dnum_ds * den - num * dden_ds) / (den * den);

    return out;
}

}  // namespace metrics
}  // namespace warp_hypad
