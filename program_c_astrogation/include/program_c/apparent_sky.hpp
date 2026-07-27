#ifndef PROGRAM_C_APPARENT_SKY_HPP
#define PROGRAM_C_APPARENT_SKY_HPP

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "program_c/catalogue.hpp"
#include "program_c/relativity.hpp"

// ---------------------------------------------------------------------------
// C1: relativistic aberration, Doppler shift, brightness and the retarded sky.
//
// THE SIGN CONVENTION, STATED BEFORE ANY FORMULA
// ----------------------------------------------
// Two unit vectors are involved and they are opposite:
//
//   k_hat  the direction a photon PROPAGATES  (from the star toward the ship)
//   n_hat  the direction TOWARD THE SOURCE    (from the ship toward the star)
//   n_hat = -k_hat
//
// Confusing them flips every aberration result, so this header never uses a
// bare "direction": every function names which one it takes. With theta the
// angle between k_hat and the observer's velocity,
//
//     cos theta' = (cos theta - beta) / (1 - beta cos theta)
//
// and with theta the angle between n_hat and the velocity,
//
//     cos theta' = (cos theta + beta) / (1 + beta cos theta)
//
// Both are implemented. Implementing only one and testing it against itself
// would pass with the sign inverted; the test that actually catches an inversion
// is the physical one -- apparent positions must CROWD TOWARD the direction of
// motion -- and that is a gate case.
// ---------------------------------------------------------------------------

namespace program_c {

// --- Aberration, scalar forms -----------------------------------------------

// theta measured from the PHOTON PROPAGATION direction.
inline double aberrate_cos_photon_direction(double cos_theta, double beta) {
    return (cos_theta - beta) / (1.0 - beta * cos_theta);
}

// theta measured from the DIRECTION TOWARD THE SOURCE.
inline double aberrate_cos_source_direction(double cos_theta, double beta) {
    return (cos_theta + beta) / (1.0 + beta * cos_theta);
}

// --- Aberration, vector form ------------------------------------------------
//
// Done by boosting the photon four-vector rather than by applying the scalar
// formula to a decomposed angle. The scalar form is then a TEST of this, not its
// implementation -- which is the right way round, because the four-vector boost
// is what generalises to an arbitrary velocity direction.
struct AberrationResult {
    Vector3 apparent_direction_to_source{0.0, 0.0, 0.0};
    Vector3 apparent_photon_direction{0.0, 0.0, 0.0};
    double doppler_factor = 1.0;      // nu' / nu
    double cos_angle_to_velocity = 0.0;   // of the apparent source direction
};

// direction_to_source: unit vector from the observer toward the star, in the
// reference (barycentric) frame. beta_vec: observer velocity in units of c.
inline AberrationResult aberrate(const Vector3& direction_to_source,
                                 const Vector3& beta_vec) {
    AberrationResult r;
    const double n = norm(direction_to_source);
    const Vector3 n_hat = (n > 0.0 ? scale(direction_to_source, 1.0 / n)
                                   : direction_to_source);
    const Vector3 k_hat = scale(n_hat, -1.0);      // photon propagates inward

    // Unit frequency: only the ratio matters, and carrying a real frequency here
    // would invite a units question that has no bearing on the geometry.
    const FourVector k = photon_four_vector(k_hat, 1.0);
    const FourVector k_ship = boost(k, beta_vec);

    // The boosted time component IS the Doppler factor, because nu' / nu =
    // -(u_ship . k)/(u_ref . k) and in the boosted frame u_ship = (1,0,0,0).
    r.doppler_factor = k_ship.t;

    const double inv = (k_ship.t != 0.0 ? 1.0 / k_ship.t : 0.0);
    r.apparent_photon_direction = {k_ship.x * inv, k_ship.y * inv, k_ship.z * inv};
    r.apparent_direction_to_source = scale(r.apparent_photon_direction, -1.0);

    const double b = norm(beta_vec);
    if (b > 0.0)
        r.cos_angle_to_velocity =
            dot(r.apparent_direction_to_source, scale(beta_vec, 1.0 / b));
    return r;
}

// --- Doppler ----------------------------------------------------------------
//
//     nu' / nu = - (u_ship . k) / (u_reference . k)
//
// Computed from the four-vectors, so the longitudinal limits are a CONSEQUENCE
// rather than the definition:
//
//   approaching (k_hat = -v_hat):  nu'/nu = gamma(1+beta) = sqrt((1+b)/(1-b))
//   receding    (k_hat = +v_hat):  nu'/nu = gamma(1-beta) = sqrt((1-b)/(1+b))
//   transverse  (k_hat . v_hat=0): nu'/nu = gamma, pure time dilation
inline double doppler_factor(const Vector3& photon_direction,
                             const Vector3& beta_vec) {
    const FourVector k = photon_four_vector(photon_direction, 1.0);
    const FourVector u_ship = four_velocity(beta_vec);
    const FourVector u_ref{1.0, 0.0, 0.0, 0.0};
    const double num = -minkowski_dot(u_ship, k);
    const double den = -minkowski_dot(u_ref, k);
    return (den != 0.0 ? num / den : 1.0);
}

// --- Brightness -------------------------------------------------------------
//
// Kept as SEPARATE named factors rather than one fused number, because the
// exponents differ and a single fudge factor would hide which one is wrong.
//
// The exponent 4 is derived, not transcribed. I_nu / nu^3 is a Lorentz
// invariant. Integrating over frequency, the bolometric intensity I = \int I_nu
// dnu picks up one further power from dnu, giving I' = D^4 I. The individual
// factors below are the pieces that multiply to it:
//
//   photon energy       D    each photon is blueshifted
//   arrival rate        D    they arrive more often
//   solid angle         D^2  the beam is concentrated forward
//                     ----
//   bolometric          D^4
struct BrightnessFactors {
    double frequency_shift = 1.0;      // D
    double arrival_rate = 1.0;         // D
    double solid_angle = 1.0;          // D^2
    double bolometric = 1.0;           // D^4, the product

    // A label, carried into every output, so a renderer cannot present this as
    // a photorealistic prediction.
    static const char* qualification() {
        return "physical approximation: bolometric and grey. It applies the "
               "D^4 specific-intensity transformation to a frequency-integrated "
               "flux, which is exact only for a source whose spectrum is "
               "unchanged in shape by the shift. It is NOT a photorealistic "
               "prediction and display exposure is applied separately.";
    }
};

inline BrightnessFactors brightness_from_doppler(double D) {
    BrightnessFactors b;
    b.frequency_shift = D;
    b.arrival_rate = D;
    b.solid_angle = D * D;
    b.bolometric = D * D * D * D;
    return b;
}

// --- Colour model -----------------------------------------------------------
//
// The fixed vocabulary. Every figure declares one of these and no RGB mapping is
// ever presented as a physically predicted observed colour.
enum class ColorModel {
    PhysicalSpectralModel,
    ApproximateBlackbodyVisualization,
    CataloguePhotometricProxy,
    IllustrativeOnly,
};

inline const char* to_string(ColorModel c) {
    switch (c) {
        case ColorModel::PhysicalSpectralModel: return "physical_spectral_model";
        case ColorModel::ApproximateBlackbodyVisualization:
            return "approximate_blackbody_visualization";
        case ColorModel::CataloguePhotometricProxy:
            return "catalogue_photometric_proxy";
        case ColorModel::IllustrativeOnly: return "illustrative_only";
    }
    return "illustrative_only";
}

// --- The retarded apparent sky ----------------------------------------------
//
// The ship at event (t_obs, x_obs) sees a star not where it is, but where it WAS
// when the light left. Solve for the emission time t_e:
//
//     | x_obs - x_star(t_e) | = c (t_obs - t_e)
//
// With x_star rectilinear this is a quadratic, but it is solved by fixed-point
// iteration anyway -- the same scheme C0 uses for light time -- so that the
// solver keeps working when C4 replaces the straight line with something else.
// Convergence is reported, never assumed.

struct RetardedSolution {
    BarycentricCartesian state_at_emission;
    double emission_epoch_years = 0.0;
    double light_time_years = 0.0;
    double geometric_distance_pc = 0.0;
    int iterations = 0;
    double residual_years = 0.0;
    bool converged = false;
};

inline RetardedSolution solve_retarded_emission(
    const Vector3& observer_position_pc, double observation_epoch_years,
    const BarycentricCartesian& star_state_at_reference,
    double tolerance_years = 1.0e-12, int max_iterations = 64) {
    RetardedSolution r;
    const double c = kSpeedOfLightParsecPerYear;

    // Start from zero light time and iterate. Each pass moves the star back
    // along its worldline to where it was one light time ago.
    double lt = 0.0;
    for (int i = 0; i < max_iterations; ++i) {
        const double t_e = observation_epoch_years - lt;
        const double dt = t_e - star_state_at_reference.epoch.value;
        const Vector3 x_star =
            add(star_state_at_reference.position_pc,
                scale(star_state_at_reference.velocity_pc_per_year, dt));
        const double d = norm(subtract(x_star, observer_position_pc));
        const double lt_new = d / c;
        r.residual_years = std::abs(lt_new - lt);
        lt = lt_new;
        r.iterations = i + 1;
        if (r.residual_years < tolerance_years) {
            r.converged = true;
            break;
        }
    }

    r.light_time_years = lt;
    r.emission_epoch_years = observation_epoch_years - lt;
    const double dt = r.emission_epoch_years - star_state_at_reference.epoch.value;
    r.state_at_emission.position_pc =
        add(star_state_at_reference.position_pc,
            scale(star_state_at_reference.velocity_pc_per_year, dt));
    r.state_at_emission.velocity_pc_per_year =
        star_state_at_reference.velocity_pc_per_year;
    r.state_at_emission.epoch = JulianYear(r.emission_epoch_years);
    r.geometric_distance_pc =
        norm(subtract(r.state_at_emission.position_pc, observer_position_pc));
    return r;
}

// --- One row of the apparent-sky table --------------------------------------

struct ApparentSkyRow {
    std::uint64_t source_id = 0;
    std::string primary_name;
    double emission_epoch = 0.0;
    double geometric_distance = 0.0;
    double apparent_direction_x = 0.0;
    double apparent_direction_y = 0.0;
    double apparent_direction_z = 0.0;
    double doppler_factor = 1.0;
    double brightness_factor = 1.0;
    double uncertainty_angle = 0.0;      // radians
    std::uint64_t quality_flags = 0;

    // Kept alongside the aberrated direction so a figure can show the shift
    // itself rather than only its result.
    double geometric_direction_x = 0.0;
    double geometric_direction_y = 0.0;
    double geometric_direction_z = 0.0;
    double aberration_angle = 0.0;       // radians between the two
};

// Build one row. The uncertainty angle is the catalogue position uncertainty
// divided by the distance -- a small-angle statement, which is safe here because
// a 0.05 pc uncertainty at 1.3 pc is 0.04 rad and the approximation errs in the
// fourth decimal.
inline ApparentSkyRow apparent_sky_row(const StellarObject& star,
                                       const BarycentricCartesian& star_state,
                                       const Vector3& observer_position_pc,
                                       double observation_epoch_years,
                                       const Vector3& observer_beta) {
    ApparentSkyRow row;
    row.source_id = star.measured.gaia_source_id;
    row.primary_name = star.primary_name;
    row.quality_flags = star.measured.quality_flags;

    const RetardedSolution ret = solve_retarded_emission(
        observer_position_pc, observation_epoch_years, star_state);
    row.emission_epoch = ret.emission_epoch_years;
    row.geometric_distance = ret.geometric_distance_pc;

    Vector3 to_source =
        subtract(ret.state_at_emission.position_pc, observer_position_pc);
    const double d = norm(to_source);
    if (d > 0.0) to_source = scale(to_source, 1.0 / d);
    row.geometric_direction_x = to_source.x;
    row.geometric_direction_y = to_source.y;
    row.geometric_direction_z = to_source.z;

    const AberrationResult ab = aberrate(to_source, observer_beta);
    row.apparent_direction_x = ab.apparent_direction_to_source.x;
    row.apparent_direction_y = ab.apparent_direction_to_source.y;
    row.apparent_direction_z = ab.apparent_direction_to_source.z;
    row.doppler_factor = ab.doppler_factor;
    row.brightness_factor = brightness_from_doppler(ab.doppler_factor).bolometric;

    const double cosang =
        std::max(-1.0, std::min(1.0, dot(to_source, ab.apparent_direction_to_source)));
    row.aberration_angle = std::acos(cosang);

    if (d > 0.0 && star.distance_uncertainty_pc > 0.0)
        row.uncertainty_angle = star.distance_uncertainty_pc / d;
    return row;
}

}  // namespace program_c

#endif  // PROGRAM_C_APPARENT_SKY_HPP
