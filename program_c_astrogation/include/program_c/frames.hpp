#pragma once

#include <string>

#include "program_c/units.hpp"

// ---------------------------------------------------------------------------
// Typed coordinate frames.
//
// An anonymous Vec3 never crosses a frame boundary. Every frame is a distinct
// type, and every transformation is a named function between two named types, so
// handing a galactocentric position to a function expecting a barycentric one is
// a compile error rather than a plausible-looking wrong answer.
//
// Each frame type carries the SIX FACTS a coordinate triple is meaningless
// without -- frame, origin, orientation, epoch, units and velocity convention --
// as static data, so they can be written into the results metadata rather than
// living only in a comment.
// ---------------------------------------------------------------------------

namespace program_c {

struct FrameDescription {
    const char* frame;
    const char* origin;
    const char* orientation;
    const char* epoch;
    const char* units;
    const char* velocity_convention;
};

// --- Spherical astrometric observables --------------------------------------

// A direction on the sky. No distance, deliberately: a direction and a position
// are different things and conflating them is how stale-catalogue errors happen.
struct ICRSDirection {
    Radian right_ascension;
    Radian declination;

    static FrameDescription description() {
        return {"ICRS", "no origin (direction only)",
                "ICRS axes, aligned with J2000 equatorial to within the ICRS "
                "realisation",
                "epoch-independent as a direction",
                "radians",
                "none (a direction carries no velocity)"};
    }
};

// The full six-parameter astrometric observable set, as catalogued.
//
// PROPER-MOTION CONVENTION: pm_ra_cosdec is mu_alpha* = mu_alpha cos(dec), which
// is what Gaia reports and what every formula here assumes. The name carries the
// cos(dec) so a reader cannot mistake it, and the convention is asserted by its
// own gate case rather than trusted.
struct ICRSPhaseSpace {
    Radian right_ascension;
    Radian declination;
    Milliarcsecond parallax;
    MasPerYear pm_ra_cosdec;   // mu_alpha* = mu_alpha cos(dec)
    MasPerYear pm_dec;
    KmPerSecond radial_velocity;
    JulianYear epoch;

    static FrameDescription description() {
        return {"ICRS spherical astrometric",
                "solar system barycentre",
                "ICRS axes",
                "the catalogue reference epoch carried in `epoch`",
                "radians, mas, mas/yr, km/s, Julian year",
                "mu_alpha* = mu_alpha cos(dec); radial velocity positive away "
                "from the barycentre"};
    }
};

// --- Cartesian frames -------------------------------------------------------

struct Vector3 {
    double x = 0.0, y = 0.0, z = 0.0;
};

inline double dot(const Vector3& a, const Vector3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline double norm(const Vector3& a) { return std::sqrt(dot(a, a)); }
inline Vector3 add(const Vector3& a, const Vector3& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
inline Vector3 subtract(const Vector3& a, const Vector3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
inline Vector3 scale(const Vector3& a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}

// Barycentric Cartesian: the internal working frame. Position in parsec,
// velocity in parsec per Julian year.
struct BarycentricCartesian {
    Vector3 position_pc;
    Vector3 velocity_pc_per_year;
    JulianYear epoch;

    static FrameDescription description() {
        return {"Barycentric Cartesian (ICRS-aligned)",
                "solar system barycentre",
                "ICRS axes: +x toward (ra=0, dec=0), +z toward the north "
                "celestial pole",
                "the epoch carried in `epoch`",
                "parsec; parsec per Julian year",
                "coordinate velocity dx/dt at the stated epoch, rectilinear "
                "model"};
    }
};

// Galactocentric Cartesian. Declared and transformed, but C0 uses it only for
// reporting: nothing in the propagation model depends on it, because rectilinear
// barycentric motion is not a Galactic orbit.
struct GalactocentricCartesian {
    Vector3 position_pc;
    Vector3 velocity_pc_per_year;
    JulianYear epoch;

    static FrameDescription description() {
        return {"Galactocentric Cartesian",
                "Galactic centre, at the assumed solar radius below",
                "+x from the Galactic centre toward the Sun's projected "
                "position, +z toward the north Galactic pole",
                "the epoch carried in `epoch`",
                "parsec; parsec per Julian year",
                "coordinate velocity; NO differential-rotation model is applied "
                "at C0, so this frame is for reporting only"};
    }
};

// The ship frame. Present so that C1 and C2 have a declared target type; at C0
// it is only ever the identity offset from barycentric, and the gate says so.
struct LocalShipFrame {
    Vector3 position_pc;
    Vector3 velocity_pc_per_year;
    JulianYear epoch;

    static FrameDescription description() {
        return {"Local ship frame",
                "the ship",
                "at C0: identical to the barycentric ICRS axes. C1 introduces "
                "the boosted and rotated frame.",
                "the epoch carried in `epoch`",
                "parsec; parsec per Julian year",
                "at C0 a pure translation of the barycentric frame; no boost, no "
                "aberration"};
    }
};

// --- Direction and tangent basis --------------------------------------------

// The ICRS unit vector of a direction:
//     r_hat = (cos(dec) cos(ra), cos(dec) sin(ra), sin(dec))
inline Vector3 unit_vector(const ICRSDirection& d) {
    const double cd = std::cos(d.declination.value);
    return {cd * std::cos(d.right_ascension.value),
            cd * std::sin(d.right_ascension.value),
            std::sin(d.declination.value)};
}

// The local tangent basis at a direction.
//
//   p_hat  = direction of increasing right ascension  = (-sin ra, cos ra, 0)
//   q_hat  = direction of increasing declination
//
// p_hat is d(r_hat)/d(ra) divided by cos(dec) -- which is exactly why the
// catalogued proper motion is mu_alpha* = mu_alpha cos(dec): it is the physical
// angular rate along p_hat, not the coordinate rate of ra.
inline void tangent_basis(const ICRSDirection& d, Vector3& p_hat, Vector3& q_hat) {
    const double sa = std::sin(d.right_ascension.value);
    const double ca = std::cos(d.right_ascension.value);
    const double sd = std::sin(d.declination.value);
    const double cd = std::cos(d.declination.value);
    p_hat = {-sa, ca, 0.0};
    q_hat = {-sd * ca, -sd * sa, cd};
}

// --- Galactocentric transformation ------------------------------------------
//
// The ICRS-to-Galactic rotation, from the standard pole and node:
//   north Galactic pole  ra = 192.85948 deg, dec = +27.12825 deg
//   Galactic longitude of the north celestial pole  l = 122.93192 deg
//
// RECORDED FROM MEMORY and marked as such: these are the conventional values,
// they are quoted to the digits given, and they must be checked against a
// primary source before any result depending on the Galactic frame is published.
// Nothing in the C0 gate's numerical assertions depends on them -- the frame is
// reporting-only, exactly so that this uncertainty cannot contaminate a verified
// number.
inline constexpr double kGalacticPoleRaDeg = 192.85948;
inline constexpr double kGalacticPoleDecDeg = 27.12825;
inline constexpr double kGalacticLonOfCelestialPoleDeg = 122.93192;
// Assumed solar galactocentric radius, also recorded from memory.
inline constexpr double kSolarGalactocentricRadiusPc = 8178.0;

inline void icrs_to_galactic_matrix(double m[3][3]) {
    const double ap = kGalacticPoleRaDeg / kDegreesPerRadian;
    const double dp = kGalacticPoleDecDeg / kDegreesPerRadian;
    const double lp = kGalacticLonOfCelestialPoleDeg / kDegreesPerRadian;

    // R = Rz(-l_ncp) * Ry(dp - 90 deg) * Rz(ap + 90 deg), the standard chain.
    const double sa = std::sin(ap), ca = std::cos(ap);
    const double sd = std::sin(dp), cd = std::cos(dp);
    const double sl = std::sin(lp), cl = std::cos(lp);

    // Rows are the Galactic axes expressed in ICRS.
    // z_gal is the pole direction.
    const double zg[3] = {cd * ca, cd * sa, sd};
    // x_gal points to the Galactic centre; y_gal completes the right-handed set.
    // Build from the node: the ICRS pole has Galactic longitude lp.
    const double np[3] = {-sa, ca, 0.0};             // node direction
    const double mp[3] = {-sd * ca, -sd * sa, cd};   // completes at the pole
    double xg[3], yg[3];
    for (int i = 0; i < 3; ++i) {
        xg[i] = -(cl * np[i] + sl * mp[i]);
        yg[i] = -(-sl * np[i] + cl * mp[i]);
    }
    for (int i = 0; i < 3; ++i) {
        m[0][i] = xg[i];
        m[1][i] = yg[i];
        m[2][i] = zg[i];
    }
}

inline Vector3 apply_matrix(const double m[3][3], const Vector3& v) {
    return {m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z};
}

// Barycentric -> galactocentric. Reporting only at C0.
inline GalactocentricCartesian to_galactocentric(const BarycentricCartesian& b) {
    double m[3][3];
    icrs_to_galactic_matrix(m);
    GalactocentricCartesian g;
    g.position_pc = apply_matrix(m, b.position_pc);
    g.position_pc.x -= kSolarGalactocentricRadiusPc;
    g.velocity_pc_per_year = apply_matrix(m, b.velocity_pc_per_year);
    g.epoch = b.epoch;
    return g;
}

// Barycentric -> ship frame. At C0 a pure translation, and nothing more is
// claimed.
inline LocalShipFrame to_ship_frame(const BarycentricCartesian& b,
                                   const Vector3& ship_position_pc,
                                   const Vector3& ship_velocity_pc_per_year) {
    LocalShipFrame s;
    s.position_pc = subtract(b.position_pc, ship_position_pc);
    s.velocity_pc_per_year =
        subtract(b.velocity_pc_per_year, ship_velocity_pc_per_year);
    s.epoch = b.epoch;
    return s;
}

}  // namespace program_c
