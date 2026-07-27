#ifndef PROGRAM_C_INTERCEPTION_HPP
#define PROGRAM_C_INTERCEPTION_HPP

#include <cmath>
#include <cstdint>
#include <string>

#include "program_c/catalogue.hpp"
#include "program_c/relativity.hpp"

// ---------------------------------------------------------------------------
// C1: moving-destination interception, and communication light cones.
//
// The destination is NOT where the catalogue says by the time the ship arrives.
// Proxima Centauri moves 3.85 arcsec/yr across the sky and closes at 22 km/s;
// over a decades-long journey that is a real miss distance, and over the 10^4 yr
// journeys a slower profile implies it is enormous.
// ---------------------------------------------------------------------------

namespace program_c {

// --- Destination ------------------------------------------------------------
//
// Refers to a real snapshot object BY IDENTIFIER. A destination is never
// identified by display name: "Proxima Centauri" is a label SIMBAD attached to
// Gaia DR3 5853498713190525696, and the number is the thing that can be
// rechecked. Names are resolved to identifiers at selection time and the
// identifier is used internally thereafter.
struct Destination {
    std::uint64_t gaia_source_id = 0;
    BarycentricCartesian state_at_reference_epoch;
    CovarianceResult uncertainty;
    Provenance provenance;

    // Carried for labelling only. Nothing keys off it.
    std::string display_name;

    BarycentricCartesian state_at(double epoch_years) const {
        BarycentricCartesian s;
        const double dt = epoch_years - state_at_reference_epoch.epoch.value;
        s.position_pc = add(state_at_reference_epoch.position_pc,
                            scale(state_at_reference_epoch.velocity_pc_per_year, dt));
        s.velocity_pc_per_year = state_at_reference_epoch.velocity_pc_per_year;
        s.epoch = JulianYear(epoch_years);
        return s;
    }
};

// --- Interception -----------------------------------------------------------

struct InterceptionResult {
    // The four positions that casual language fuses into "the destination".
    Vector3 stale_catalogue_target{0.0, 0.0, 0.0};
    Vector3 destination_at_departure{0.0, 0.0, 0.0};
    Vector3 destination_at_arrival{0.0, 0.0, 0.0};
    Vector3 interception_point{0.0, 0.0, 0.0};

    double departure_epoch_years = 0.0;
    double arrival_epoch_years = 0.0;
    double flight_coordinate_time_years = 0.0;
    double interception_distance_pc = 0.0;      // the distance actually flown

    // What ignoring proper motion costs: the separation between where the ship
    // would arrive aiming at the stale catalogue position and where the star
    // actually is at arrival.
    double stale_targeting_error_pc = 0.0;
    double stale_targeting_error_au = 0.0;

    // The catalogue covariance carried forward to the arrival epoch.
    double interception_uncertainty_pc = 0.0;

    // Solver reporting: never assumed, always stated.
    int iterations = 0;
    double residual_pc = 0.0;
    bool converged = false;
    std::string method;
};

// Solve x_ship(t_arrival) = x_star(t_arrival).
//
// The ship flies a straight coordinate path from the origin toward the
// interception point, covering coordinate distance D in coordinate time T(D)
// given by the trajectory model. The residual
//
//     f(t) = |x_star(t)| - D_flown(t - t_departure)
//
// is monotone decreasing in t for any profile whose speed is bounded by c and
// any destination receding slower than light, so BISECTION cannot diverge. That
// is why it is used in preference to Newton, which would converge faster and
// could leave the bracket.
//
// time_of_flight_for_distance must be a callable double(double distance_pc).
template <typename TimeOfFlightFn>
InterceptionResult solve_interception(const Destination& dest,
                                      double departure_epoch_years,
                                      TimeOfFlightFn time_of_flight_for_distance,
                                      double tolerance_pc = 1.0e-12,
                                      int max_iterations = 200) {
    InterceptionResult r;
    r.departure_epoch_years = departure_epoch_years;
    r.method =
        "bracketed bisection on |x_star(t)| - D_flown(t), then secant "
        "refinement; bisection because the residual is monotone and bracketing "
        "cannot diverge";

    r.stale_catalogue_target = dest.state_at_reference_epoch.position_pc;
    r.destination_at_departure = dest.state_at(departure_epoch_years).position_pc;

    // f(t) > 0 means the star is still farther away than the ship has flown.
    auto residual = [&](double t) {
        const double d_star = norm(dest.state_at(t).position_pc);
        const double flown_time = t - departure_epoch_years;
        // Invert the time-of-flight relation by asking how far the ship has gone
        // in flown_time. Rather than invert analytically for each model, bisect
        // on distance too -- one nested solve, but model-agnostic.
        double lo = 0.0, hi = d_star * 2.0 + 1.0;
        for (int i = 0; i < 80; ++i) {
            const double mid = 0.5 * (lo + hi);
            if (time_of_flight_for_distance(mid) < flown_time) lo = mid;
            else hi = mid;
        }
        return d_star - 0.5 * (lo + hi);
    };

    // Bracket: at departure the ship has flown nothing, so f > 0. Expand until
    // f < 0.
    double t_lo = departure_epoch_years;
    double t_hi = departure_epoch_years + 1.0;
    double f_lo = residual(t_lo);
    for (int i = 0; i < 200 && residual(t_hi) > 0.0; ++i)
        t_hi = departure_epoch_years + (t_hi - departure_epoch_years) * 2.0;

    for (int i = 0; i < max_iterations; ++i) {
        const double t_mid = 0.5 * (t_lo + t_hi);
        const double f_mid = residual(t_mid);
        r.iterations = i + 1;
        r.residual_pc = std::abs(f_mid);
        if (r.residual_pc < tolerance_pc) { t_lo = t_hi = t_mid; r.converged = true; break; }
        if ((f_mid > 0.0) == (f_lo > 0.0)) { t_lo = t_mid; f_lo = f_mid; }
        else t_hi = t_mid;
    }

    r.arrival_epoch_years = 0.5 * (t_lo + t_hi);
    r.flight_coordinate_time_years = r.arrival_epoch_years - departure_epoch_years;
    r.destination_at_arrival = dest.state_at(r.arrival_epoch_years).position_pc;
    r.interception_point = r.destination_at_arrival;
    r.interception_distance_pc = norm(r.interception_point);

    // The cost of aiming at the stale position: fly the same distance along the
    // stale direction, and see how far that lands from the star.
    const double stale_norm = norm(r.stale_catalogue_target);
    if (stale_norm > 0.0) {
        const Vector3 stale_aim =
            scale(r.stale_catalogue_target,
                  r.interception_distance_pc / stale_norm);
        r.stale_targeting_error_pc =
            norm(subtract(stale_aim, r.destination_at_arrival));
        // Astronomical units make the number legible: a targeting error of
        // 1e-4 pc is 20 au, which is a statement about whether you arrive in the
        // system or merely near it.
        r.stale_targeting_error_au =
            r.stale_targeting_error_pc * kParsecMetres / (kAuKm * 1000.0);
    }

    // Uncertainty at arrival: the catalogue position sigma grows with the
    // velocity sigma over the flight. Same linear-growth model C0 verified.
    r.interception_uncertainty_pc = propagated_position_rms_pc(
        dest.uncertainty, r.arrival_epoch_years -
                              dest.state_at_reference_epoch.epoch.value);
    return r;
}

// --- Communication light cones ----------------------------------------------
//
// Earth's clock as SEEN from the ship is not Earth's clock. For a ship at
// (t, x), the latest Earth event whose light has arrived is
//
//     t_earth_seen = t - |x| / c
//
// because Earth sits at the barycentric origin and does not move in this frame.
// A message the ship sends now reaches Earth at t + |x|/c, and reaches the
// destination when the outgoing light cone meets the destination's worldline --
// which requires a solve, because the destination moves.
struct SignalLedger {
    double ship_coordinate_time = 0.0;
    double ship_proper_time = 0.0;

    double earth_time_seen_from_ship = 0.0;    // NOT "the time on Earth now"
    double earth_signal_age_years = 0.0;       // how stale that reading is

    double ship_signal_arrives_at_earth = 0.0;
    double ship_signal_arrives_at_destination = 0.0;

    double destination_time_seen_from_ship = 0.0;
    double round_trip_to_earth_years = 0.0;
};

inline SignalLedger signal_ledger(const Vector3& ship_position_pc,
                                  double ship_coordinate_time,
                                  double ship_proper_time,
                                  const Destination& dest) {
    SignalLedger s;
    const double c = kSpeedOfLightParsecPerYear;
    s.ship_coordinate_time = ship_coordinate_time;
    s.ship_proper_time = ship_proper_time;

    const double d_earth = norm(ship_position_pc);
    s.earth_signal_age_years = d_earth / c;
    s.earth_time_seen_from_ship = ship_coordinate_time - s.earth_signal_age_years;
    s.ship_signal_arrives_at_earth = ship_coordinate_time + d_earth / c;
    s.round_trip_to_earth_years = 2.0 * d_earth / c;

    // Outgoing signal to a MOVING destination: solve
    //   |x_dest(t_arrive) - x_ship| = c (t_arrive - t_now)
    // by the same fixed-point scheme used for the retarded solution, run
    // forward instead of backward.
    double lt = norm(subtract(dest.state_at(ship_coordinate_time).position_pc,
                              ship_position_pc)) / c;
    for (int i = 0; i < 64; ++i) {
        const double t_arrive = ship_coordinate_time + lt;
        const double d = norm(subtract(dest.state_at(t_arrive).position_pc,
                                       ship_position_pc));
        const double lt_new = d / c;
        if (std::abs(lt_new - lt) < 1.0e-12) { lt = lt_new; break; }
        lt = lt_new;
    }
    s.ship_signal_arrives_at_destination = ship_coordinate_time + lt;

    // And what the ship currently SEES of the destination: the retarded solve.
    double back = 0.0;
    for (int i = 0; i < 64; ++i) {
        const double t_e = ship_coordinate_time - back;
        const double d =
            norm(subtract(dest.state_at(t_e).position_pc, ship_position_pc));
        const double back_new = d / c;
        if (std::abs(back_new - back) < 1.0e-12) { back = back_new; break; }
        back = back_new;
    }
    s.destination_time_seen_from_ship = ship_coordinate_time - back;
    return s;
}

}  // namespace program_c

#endif  // PROGRAM_C_INTERCEPTION_HPP
