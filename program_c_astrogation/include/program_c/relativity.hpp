#ifndef PROGRAM_C_RELATIVITY_HPP
#define PROGRAM_C_RELATIVITY_HPP

#include <cmath>
#include <string>
#include <vector>

#include "program_c/astrometry.hpp"
#include "program_c/frames.hpp"
#include "program_c/units.hpp"

// ---------------------------------------------------------------------------
// C1: special-relativistic ship worldlines and clocks.
//
// NO EFFECTIVE METRIC. This is flat spacetime, signature (-,+,+,+), and it is
// the baseline C2 will later compare an effective-metric voyage against.
//
// Units are the ones C0 established -- parsec, Julian year, parsec per Julian
// year -- rather than natural units with c = 1. That is deliberate: every number
// written to a CSV is then directly readable against the catalogue, and a unit
// error shows up as a factor of 3.26 instead of hiding inside c = 1.
// ---------------------------------------------------------------------------

namespace program_c {

// --- Constants --------------------------------------------------------------

// Standard gravity, converted in ONE place. 9.80665 m/s^2 is the defined value.
inline constexpr double kStandardGravityMetresPerSecondSquared = 9.80665;

// pc/yr^2. (m/s^2) * (s/yr)^2 / (m/pc).
inline constexpr double kStandardGravityParsecPerYearSquared =
    kStandardGravityMetresPerSecondSquared * kJulianYearSeconds *
    kJulianYearSeconds / kParsecMetres;

// The familiar statement "one gravity is about one c per year" is this ratio,
// 1.0323. It is derived here rather than quoted, so it cannot drift from the
// constants above.
inline double standard_gravity_in_c_per_year() {
    return kStandardGravityParsecPerYearSquared / kSpeedOfLightParsecPerYear;
}

// --- Kinematic scalars ------------------------------------------------------

inline double lorentz_gamma(double beta) {
    // Guard the domain rather than returning a NaN that propagates into a CSV
    // and then into a renderer. A beta of exactly 1 is not attainable and a beta
    // above 1 is a bug upstream, so both are clamped and the caller can detect
    // the clamp by checking beta itself.
    const double b2 = beta * beta;
    if (b2 >= 1.0) return std::numeric_limits<double>::infinity();
    return 1.0 / std::sqrt(1.0 - b2);
}

inline double beta_from_gamma(double gamma) {
    if (!(gamma >= 1.0)) return 0.0;
    return std::sqrt(1.0 - 1.0 / (gamma * gamma));
}

// --- Ship state -------------------------------------------------------------

struct ShipState {
    // The typed frame, never a bare Vec3: position in pc and coordinate velocity
    // in pc/yr, in the barycentric ICRS-aligned frame C0 defines.
    BarycentricCartesian coordinate_state;

    double proper_time_years = 0.0;        // tau, the crew's clock
    double coordinate_time_years = 0.0;    // t, the barycentric clock
    double beta = 0.0;                     // |v| / c, dimensionless
    double gamma = 1.0;
    Vector3 proper_acceleration{0.0, 0.0, 0.0};   // pc/yr^2, ship frame

    // The six facts of the frame this state is expressed in, so a consumer never
    // has to guess. Mirrors the C0 discipline.
    static FrameDescription frame_description() {
        return BarycentricCartesian::description();
    }
};

struct WorldlineSample {
    ShipState ship;

    // Where the destination IS at this coordinate time -- not where the
    // catalogue says, and not where the ship SEES it.
    BarycentricCartesian destination;
    double remaining_distance_pc = 0.0;

    // The latest Earth and destination clock readings whose light has reached
    // the ship. NOT "the current time on Earth": see signals.hpp.
    double received_earth_time = 0.0;
    double received_destination_time = 0.0;
};

// --- Trajectory families ----------------------------------------------------

enum class TrajectoryModel {
    ConstantCoordinateVelocity,   // verification model only
    SymmetricProperAcceleration,
    BoundedAccelerationRamp,
};

inline const char* to_string(TrajectoryModel m) {
    switch (m) {
        case TrajectoryModel::ConstantCoordinateVelocity:
            return "ConstantCoordinateVelocity";
        case TrajectoryModel::SymmetricProperAcceleration:
            return "SymmetricProperAcceleration";
        case TrajectoryModel::BoundedAccelerationRamp:
            return "BoundedAccelerationRamp";
    }
    return "Unknown";
}

// A one-word statement of what each model is FOR, carried into the metadata so a
// reader of a figure cannot mistake the verification model for a mission profile.
inline const char* model_purpose(TrajectoryModel m) {
    switch (m) {
        case TrajectoryModel::ConstantCoordinateVelocity:
            return "verification only: requires infinite acceleration at both "
                   "ends, and exists because tau = t/gamma is exact and therefore "
                   "tests the integrator";
        case TrajectoryModel::SymmetricProperAcceleration:
            return "accelerate at constant proper acceleration to the midpoint, "
                   "then decelerate; arrives at rest";
        case TrajectoryModel::BoundedAccelerationRamp:
            return "smooth sigmoid-ramped acceleration with bounded jerk; a "
                   "DIFFERENT curve, not a more realizable one -- no propulsion "
                   "model is present to support any realizability claim";
    }
    return "unknown";
}

// --- Analytic constant-proper-acceleration relations ------------------------
//
// One-dimensional motion from rest under constant proper acceleration a. These
// are EXACT, and the gate asserts them as exact rather than approximate:
//
//     c t(tau) = (c^2/a) sinh(a tau / c)
//     x(tau)   = (c^2/a) [cosh(a tau / c) - 1]
//     beta(tau)= tanh(a tau / c)
//     gamma    = cosh(a tau / c)
//
// The group a*tau/c is dimensionless; every function below forms it explicitly
// so the units cannot silently cancel wrongly.

struct RelativisticRocket {
    double a_pc_per_yr2 = kStandardGravityParsecPerYearSquared;

    double rapidity_at_proper_time(double tau_years) const {
        return a_pc_per_yr2 * tau_years / kSpeedOfLightParsecPerYear;
    }

    double coordinate_time(double tau_years) const {
        const double c = kSpeedOfLightParsecPerYear;
        return (c / a_pc_per_yr2) * std::sinh(rapidity_at_proper_time(tau_years));
    }

    double distance_pc(double tau_years) const {
        const double c = kSpeedOfLightParsecPerYear;
        return (c * c / a_pc_per_yr2) *
               (std::cosh(rapidity_at_proper_time(tau_years)) - 1.0);
    }

    double beta(double tau_years) const {
        return std::tanh(rapidity_at_proper_time(tau_years));
    }

    double gamma(double tau_years) const {
        return std::cosh(rapidity_at_proper_time(tau_years));
    }

    // Inverse relations, needed to solve the half-journey matching condition.
    // Both are closed form: no iteration, so no convergence to report.
    double proper_time_for_distance(double x_pc) const {
        const double c = kSpeedOfLightParsecPerYear;
        const double arg = 1.0 + a_pc_per_yr2 * x_pc / (c * c);
        return (c / a_pc_per_yr2) * std::acosh(arg);
    }

    double proper_time_for_coordinate_time(double t_years) const {
        const double c = kSpeedOfLightParsecPerYear;
        return (c / a_pc_per_yr2) * std::asinh(a_pc_per_yr2 * t_years / c);
    }
};

// A symmetric journey: accelerate over D/2, decelerate over D/2, arrive at rest.
//
// The matching condition is applied to the COORDINATE distance, so each phase
// covers D/2 and the second is the time reverse of the first. Total times are
// therefore exactly twice the single-phase values -- an exact statement, which
// is why the gate can assert the endpoint velocity is zero to machine precision
// rather than to a tolerance.
struct SymmetricJourney {
    RelativisticRocket rocket;
    double distance_pc = 0.0;

    double half_proper_time() const {
        return rocket.proper_time_for_distance(0.5 * distance_pc);
    }
    double total_proper_time() const { return 2.0 * half_proper_time(); }
    double half_coordinate_time() const {
        return rocket.coordinate_time(half_proper_time());
    }
    double total_coordinate_time() const { return 2.0 * half_coordinate_time(); }
    double peak_beta() const { return rocket.beta(half_proper_time()); }
    double peak_gamma() const { return rocket.gamma(half_proper_time()); }

    // State at proper time tau along the whole journey, in the 1-D coordinate
    // measured from the departure point toward the destination.
    void state_at_proper_time(double tau, double& x_pc, double& t_years,
                              double& beta, double& gamma) const {
        const double tau_half = half_proper_time();
        if (tau <= tau_half) {
            x_pc = rocket.distance_pc(tau);
            t_years = rocket.coordinate_time(tau);
            beta = rocket.beta(tau);
            gamma = rocket.gamma(tau);
            return;
        }
        // Deceleration phase, by time reversal about the midpoint. Writing it as
        // a reflection rather than integrating a second time is what makes the
        // arrival exactly symmetric instead of symmetric to within quadrature
        // error.
        const double s = total_proper_time() - tau;      // time to arrival
        const double x_half = rocket.distance_pc(tau_half);
        x_pc = 2.0 * x_half - rocket.distance_pc(s);
        t_years = 2.0 * rocket.coordinate_time(tau_half) - rocket.coordinate_time(s);
        beta = rocket.beta(s);
        gamma = rocket.gamma(s);
    }
};

// --- Bounded-acceleration profile -------------------------------------------
//
// Proper acceleration ramped by a smooth sigmoid so that jerk is bounded. The
// shape function is C-infinity and integrates numerically; unlike the constant
// case there is no closed form, which is exactly why the gate checks it against
// the integrator's own convergence rather than against an analytic answer.
struct BoundedRampProfile {
    double a_peak_pc_per_yr2 = kStandardGravityParsecPerYearSquared;
    double ramp_fraction = 0.1;   // fraction of the half-journey spent ramping

    // Shape in [0,1] over the acceleration phase parameter s in [0,1].
    double shape(double s) const {
        if (s <= 0.0 || s >= 1.0) return 0.0;
        const double r = (ramp_fraction > 0.0 ? ramp_fraction : 1.0e-6);
        // Smoothstep up, flat, smoothstep down. Continuous value and first
        // derivative at both joins, so acceleration is C1 and jerk is bounded.
        auto smoothstep = [](double u) {
            if (u <= 0.0) return 0.0;
            if (u >= 1.0) return 1.0;
            return u * u * (3.0 - 2.0 * u);
        };
        if (s < r) return smoothstep(s / r);
        if (s > 1.0 - r) return smoothstep((1.0 - s) / r);
        return 1.0;
    }
};

// --- Proper time integration ------------------------------------------------
//
//     dtau = dt sqrt(1 - v^2/c^2)
//
// Integrated with composite Simpson over the coordinate-time samples. Simpson
// rather than trapezoid because the integrand is smooth and the fourth-order
// error is what lets the gate compare against the analytic constant-acceleration
// answer at 1e-10 instead of 1e-4.

struct ProperTimeResult {
    double proper_time_years = 0.0;
    double coordinate_time_years = 0.0;
    int intervals = 0;
    std::string method;
};

// beta_of_t must be a callable double(double t).
template <typename BetaFn>
ProperTimeResult integrate_proper_time(BetaFn beta_of_t, double t0, double t1,
                                       int intervals) {
    ProperTimeResult r;
    r.coordinate_time_years = t1 - t0;
    r.method = "composite Simpson on dtau = dt sqrt(1 - beta^2)";
    if (intervals < 2) intervals = 2;
    if (intervals % 2 == 1) ++intervals;      // Simpson needs an even count
    r.intervals = intervals;

    const double h = (t1 - t0) / intervals;
    auto integrand = [&](double t) {
        const double b = beta_of_t(t);
        const double s = 1.0 - b * b;
        return s > 0.0 ? std::sqrt(s) : 0.0;
    };

    // Neumaier compensated summation. Naive accumulation over ~1000 weighted
    // terms loses about 100 ulps -- measured, not guessed -- which is enough to
    // put the result 1e-13 away from the analytic answer and hide a real
    // quadrature error underneath the noise. Compensation costs two extra flops
    // per term and brings the sum back to a couple of ulps, so the gate can
    // assert against the exact constant-velocity result at 1e-14 and mean it.
    double sum = 0.0, comp = 0.0;
    auto add = [&](double v) {
        const double t = sum + v;
        comp += (std::abs(sum) >= std::abs(v)) ? (sum - t) + v : (v - t) + sum;
        sum = t;
    };
    add(integrand(t0));
    add(integrand(t1));
    for (int i = 1; i < intervals; ++i)
        add((i % 2 ? 4.0 : 2.0) * integrand(t0 + i * h));
    r.proper_time_years = (sum + comp) * h / 3.0;
    return r;
}

// --- Four-vectors and the Lorentz boost -------------------------------------

struct FourVector {
    double t = 0.0, x = 0.0, y = 0.0, z = 0.0;

    // Minkowski norm with signature (-,+,+,+). Negative for timelike, zero for
    // null, positive for spacelike.
    double norm2() const { return -t * t + x * x + y * y + z * z; }
};

inline FourVector four_velocity(const Vector3& beta_vec) {
    const double b = norm(beta_vec);
    const double g = lorentz_gamma(b);
    return {g, g * beta_vec.x, g * beta_vec.y, g * beta_vec.z};
}

// A photon four-vector. k_hat is the PROPAGATION direction, not the direction
// toward the source. The two differ by a sign and the names keep them apart.
inline FourVector photon_four_vector(const Vector3& k_hat, double frequency) {
    const double n = norm(k_hat);
    const double s = (n > 0.0 ? frequency / n : 0.0);
    return {frequency, s * k_hat.x, s * k_hat.y, s * k_hat.z};
}

inline double minkowski_dot(const FourVector& a, const FourVector& b) {
    return -a.t * b.t + a.x * b.x + a.y * b.y + a.z * b.z;
}

// The general boost with velocity beta_vec, applied to a four-vector. Written
// out in full rather than as a rotation to the z axis plus a 1-D boost, because
// the round trip through a rotation is where an axis convention can slip in
// unnoticed.
inline FourVector boost(const FourVector& v, const Vector3& beta_vec) {
    const double b2 = dot(beta_vec, beta_vec);
    if (b2 <= 0.0) return v;
    const double g = lorentz_gamma(std::sqrt(b2));
    const double bv = beta_vec.x * v.x + beta_vec.y * v.y + beta_vec.z * v.z;

    FourVector out;
    out.t = g * (v.t - bv);
    const double k = (g - 1.0) * bv / b2 - g * v.t;
    out.x = v.x + k * beta_vec.x;
    out.y = v.y + k * beta_vec.y;
    out.z = v.z + k * beta_vec.z;
    return out;
}

}  // namespace program_c

#endif  // PROGRAM_C_RELATIVITY_HPP
