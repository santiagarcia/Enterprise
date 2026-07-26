#pragma once

#include <cmath>

// ---------------------------------------------------------------------------
// Strong units.
//
// A public interface in Program C never carries a bare `double distance` or
// `double velocity`. Astrometry is the discipline where unit confusion is not a
// stylistic worry but the standard failure mode: milliarcseconds against
// arcseconds, parsecs against light-years, km/s against pc/yr, and
// mu_alpha against mu_alpha* = mu_alpha cos(dec). Each of those has produced
// published errors, so each gets its own type here and the compiler refuses the
// mixture.
//
// INTERNAL DYNAMICAL SYSTEM: parsec, Julian year, parsec per Julian year.
// Every boundary conversion is a named function with its factor written out, and
// the factors are collected in one place so a reader can check them all at once.
// ---------------------------------------------------------------------------

namespace program_c {

// --- The conversion factors, all in one place ------------------------------

// 1 AU in km, IAU 2012 definition (exact).
inline constexpr double kAuKm = 149597870.7;
// A Julian year in seconds (exact by definition: 365.25 * 86400).
inline constexpr double kJulianYearSeconds = 31557600.0;
// 1 AU/yr in km/s. This is the constant that turns mas/yr at a distance in pc
// into a transverse velocity: 149597870.7 / 31557600 = 4.740470446...
inline constexpr double kAuPerYearInKmPerSecond = kAuKm / kJulianYearSeconds;
// 1 pc in metres, IAU 2015: 648000/pi AU.
inline constexpr double kParsecMetres = 3.0856775814913673e16;
// 1 pc in light-years, derived: kParsecMetres / (c * Julian year).
inline constexpr double kSpeedOfLightKmPerSecond = 299792.458;
inline constexpr double kParsecKm = kParsecMetres / 1000.0;
// c in pc/yr, derived rather than quoted: c[km/s] * yr[s] / pc[km].
inline constexpr double kSpeedOfLightParsecPerYear =
    kSpeedOfLightKmPerSecond * kJulianYearSeconds / kParsecKm;
// 1 pc/yr in km/s, the inverse boundary.
inline constexpr double kParsecPerYearInKmPerSecond =
    kParsecKm / kJulianYearSeconds;

inline constexpr double kPi = 3.14159265358979323846;
inline constexpr double kMasPerRadian = 180.0 * 3600.0 * 1000.0 / kPi;
inline constexpr double kDegreesPerRadian = 180.0 / kPi;

// --- The types --------------------------------------------------------------
//
// Deliberately explicit: no implicit construction from double, and no implicit
// conversion between units. The only way across a boundary is a named function.

#define PROGRAM_C_UNIT(NAME)                                              \
    struct NAME {                                                         \
        double value = 0.0;                                               \
        NAME() = default;                                                 \
        explicit constexpr NAME(double v) : value(v) {}                   \
        friend constexpr NAME operator+(NAME a, NAME b) {                 \
            return NAME(a.value + b.value);                               \
        }                                                                 \
        friend constexpr NAME operator-(NAME a, NAME b) {                 \
            return NAME(a.value - b.value);                               \
        }                                                                 \
        friend constexpr NAME operator*(NAME a, double s) {               \
            return NAME(a.value * s);                                     \
        }                                                                 \
        friend constexpr NAME operator*(double s, NAME a) {               \
            return NAME(a.value * s);                                      \
        }                                                                 \
        friend constexpr NAME operator/(NAME a, double s) {               \
            return NAME(a.value / s);                                     \
        }                                                                 \
        friend constexpr NAME operator-(NAME a) { return NAME(-a.value); }\
    }

PROGRAM_C_UNIT(Radian);
PROGRAM_C_UNIT(Degree);
PROGRAM_C_UNIT(Milliarcsecond);
PROGRAM_C_UNIT(MasPerYear);
PROGRAM_C_UNIT(Parsec);
PROGRAM_C_UNIT(Metre);
PROGRAM_C_UNIT(JulianYear);
PROGRAM_C_UNIT(Second);
PROGRAM_C_UNIT(KmPerSecond);
PROGRAM_C_UNIT(MetrePerSecond);
PROGRAM_C_UNIT(ParsecPerYear);

#undef PROGRAM_C_UNIT

// --- Named boundary conversions --------------------------------------------

inline constexpr Radian to_radian(Degree d) {
    return Radian(d.value / kDegreesPerRadian);
}
inline constexpr Degree to_degree(Radian r) {
    return Degree(r.value * kDegreesPerRadian);
}
inline constexpr Radian to_radian(Milliarcsecond m) {
    return Radian(m.value / kMasPerRadian);
}
inline constexpr Milliarcsecond to_milliarcsecond(Radian r) {
    return Milliarcsecond(r.value * kMasPerRadian);
}
inline constexpr Metre to_metre(Parsec p) { return Metre(p.value * kParsecMetres); }
inline constexpr Parsec to_parsec(Metre m) { return Parsec(m.value / kParsecMetres); }
inline constexpr Second to_second(JulianYear y) {
    return Second(y.value * kJulianYearSeconds);
}
inline constexpr JulianYear to_julian_year(Second s) {
    return JulianYear(s.value / kJulianYearSeconds);
}
inline constexpr MetrePerSecond to_metre_per_second(KmPerSecond v) {
    return MetrePerSecond(v.value * 1000.0);
}
inline constexpr KmPerSecond to_km_per_second(MetrePerSecond v) {
    return KmPerSecond(v.value / 1000.0);
}
inline constexpr ParsecPerYear to_parsec_per_year(KmPerSecond v) {
    return ParsecPerYear(v.value / kParsecPerYearInKmPerSecond);
}
inline constexpr KmPerSecond to_km_per_second(ParsecPerYear v) {
    return KmPerSecond(v.value * kParsecPerYearInKmPerSecond);
}

// Inverse parallax. NOT a general distance estimator -- see catalogue.hpp, where
// the quality gate that licenses it lives. This function performs the arithmetic
// and nothing else, and it is never called except behind that gate.
inline Parsec inverse_parallax_distance(Milliarcsecond parallax) {
    return Parsec(1000.0 / parallax.value);
}
inline Milliarcsecond parallax_from_distance(Parsec d) {
    return Milliarcsecond(1000.0 / d.value);
}

// The transverse-velocity conversion, the one astrometry gets wrong most often:
//
//     v_t [km/s] = 4.740470446 * mu["/yr] * d[pc]
//
// with mu in mas/yr the factor carries an extra 1/1000.
inline KmPerSecond transverse_velocity(MasPerYear mu, Parsec d) {
    return KmPerSecond(kAuPerYearInKmPerSecond * (mu.value / 1000.0) * d.value);
}
inline MasPerYear proper_motion_from_transverse(KmPerSecond v_t, Parsec d) {
    return MasPerYear(1000.0 * v_t.value / (kAuPerYearInKmPerSecond * d.value));
}

// Light travel time over a distance, in the internal system.
inline JulianYear light_travel_time(Parsec d) {
    return JulianYear(d.value / kSpeedOfLightParsecPerYear);
}

// --- Reporting precision ----------------------------------------------------
//
// Derived coordinates must never be printed to more precision than the
// measurements support. Given a value and its uncertainty, this is the number of
// significant digits that are actually meaningful: one digit past the leading
// digit of the uncertainty, floored at 1 and capped at 17.
inline int significant_digits_for(double value, double uncertainty) {
    if (!(uncertainty > 0.0) || !std::isfinite(uncertainty)) return 17;
    if (!(std::fabs(value) > 0.0)) return 1;
    const double d = std::log10(std::fabs(value) / uncertainty) + 1.0;
    if (!(d > 1.0)) return 1;
    if (d > 17.0) return 17;
    return static_cast<int>(d);
}

}  // namespace program_c
