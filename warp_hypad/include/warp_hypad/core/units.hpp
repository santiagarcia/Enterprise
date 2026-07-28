#pragma once

#include <limits>
#include <string>

// ---------------------------------------------------------------------------
// Units, and the one conversion that crosses between them.
//
// THE PROBLEM THIS SOLVES. The geometry code is dimensionless: grid coordinates
// are pure numbers, the metric is dimensionless, and G = c = 1. An SI energy
// density is J/m^3. Those two live in different arithmetic, and the classic way
// to produce a confidently wrong number is to multiply a dimensionless curvature
// by c^4/(8 pi G_N) and print the result as if the coordinates had been metres.
// They were not. A coordinate of 2.0 means two of WHATEVER LENGTH UNIT the model
// chose, and until that unit is named the SI value does not exist.
//
// So conversion is a deliberate, single-entry-point operation that REQUIRES the
// scale, and every converted output records the scale that produced it.
//
// THE CHAIN.
//
//   G_{mu nu}                 dimensionless, as computed on the grid
//   G^SI_{mu nu} = G/L0^2     1/m^2   -- curvature has dimensions of 1/length^2,
//                                        and the grid's length unit is L0 metres
//   T^SI = c^4/(8 pi G_N) * G^SI      J/m^3
//
// In geometrized units the middle step does not exist and T = G/(8 pi) is
// dimensionless. Both are provided; neither is silently the other.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace units {

// CODATA 2018.
inline constexpr double kSpeedOfLight = 299792458.0;              // m/s, exact
inline constexpr double kNewtonG = 6.67430e-11;                   // m^3 kg^-1 s^-2
inline constexpr double kPi = 3.14159265358979323846;

// c^4 / (8 pi G_N), the factor relating an SI curvature (1/m^2) to an SI
// stress-energy (J/m^3 = Pa). About 4.815e42.
inline double einstein_coupling_si() {
    const double c2 = kSpeedOfLight * kSpeedOfLight;
    return (c2 * c2) / (8.0 * kPi * kNewtonG);
}

// The scale that gives the dimensionless model its physical size. Without one,
// no SI quantity can be formed, which is why there is no default.
struct CoordinateScale {
    // One unit of grid coordinate, in metres.
    double length_metres = 0.0;

    // One unit of coordinate time, in seconds. With c = 1 internally the two are
    // locked together: t and x share a unit, so time_seconds = length/c. Stored
    // explicitly anyway so a record is readable without knowing that.
    double time_seconds() const { return length_metres / kSpeedOfLight; }

    std::string name;  // e.g. "bubble radius R = 2 corresponds to 100 m"

    bool valid() const { return length_metres > 0.0; }
};

// Geometrized stress-energy: T = G / (8 pi). Dimensionless in, dimensionless
// out. This is the quantity the A5 gate verifies; SI is a presentation layer.
inline double geometrized_stress_energy(double einstein_component) {
    return einstein_component / (8.0 * kPi);
}

// SI stress-energy, in J/m^3, from a DIMENSIONLESS Einstein component and the
// scale that gives the coordinates their size.
//
//   T^SI = c^4/(8 pi G_N) * G_dimensionless / L0^2
inline double si_stress_energy(double einstein_component,
                               const CoordinateScale& scale) {
    if (!scale.valid()) return std::numeric_limits<double>::quiet_NaN();  // the question is ill-posed
    return einstein_coupling_si() * einstein_component /
           (scale.length_metres * scale.length_metres);
}

// SI curvature, 1/m^2, from a dimensionless one.
inline double si_curvature(double dimensionless, const CoordinateScale& scale) {
    if (!scale.valid()) return std::numeric_limits<double>::quiet_NaN();
    return dimensionless / (scale.length_metres * scale.length_metres);
}

// What every converted output must record, so a reader can reconstruct or
// challenge the number without the source.
struct UnitsRecord {
    const char* system = "geometrized (G = c = 1)";
    const char* coordinate_length_units = "dimensionless grid units";
    const char* coordinate_time_units =
        "dimensionless grid units, same unit as length (c = 1)";
    const char* curvature_units = "1 / (coordinate length)^2";
    const char* stress_energy_units_geometrized =
        "1 / (coordinate length)^2, i.e. G/(8 pi); dimensionless";
    const char* stress_energy_units_si = "J / m^3 (= Pa)";
    const char* conversion_note =
        "SI values require an explicit CoordinateScale. A grid coordinate is a "
        "pure number; it is not metres. No SI quantity is emitted without the "
        "scale that produced it.";
};

}  // namespace units
}  // namespace warp_hypad

